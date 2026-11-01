"""Asynchronous q4 transistor experiment; no behavioral detector/transfer.

External sources are VDD, initialization and additive input charge only.
Continuous PMOS-input detectors set CMOS SR latches; an RC reset closes each
one-shot. Weighted, gated current mirrors physically remove/add charge.
This open-loop packet candidate must pass charge accuracy before energy counts.
"""
from dataclasses import asdict
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse, json, hashlib, itertools
import numpy as np
from gamma_sky130 import Circuit, Design, execute, compress_trace, pwl, ROOT

OUT=ROOT/'results/self_triggered_q4'
SELECTED={'guarded':True,'cap_ff':100,'input_units':6,'input_a':1e-6,
          'packet_r':1e8,'reset_r':1.6e7,'n_scale':1.18,'p_scale':1,
          'center_n_scale':.90,'center_p_scale':1.69}

def nor(c,name,inputs,out,dd,w=.42):
    for j,g in enumerate(inputs):
        c.fet(name+'n'+str(j),'n',out,g,'0','0',w)
        c.fet(name+'p'+str(j),'p',out if j==len(inputs)-1 else name+'s'+str(j),g,
              dd if j==0 else name+'s'+str(j-1),dd,2*w*len(inputs))

def build(cfg):
    d=Design(corner=cfg.get('corner','tt'),vdd=cfg.get('vdd',1.8),temp_c=cfg.get('temp',25),
             cap_ff=cfg.get('cap_ff',25),cap_kind='mim',pdk_mismatch=cfg.get('mismatch',False),seed=cfg.get('seed',101),
             capacitor_sigma=cfg.get('cap_sigma',0),resistor_sigma=cfg.get('resistor_sigma',0))
    c=Circuit(d);dd={k:c.supply(k) for k in ('detector','logic','bias','transfer','reference')}
    # Ratiometric references. Their continuous dissipation is measured.
    for name,a,b,r in [('lo','lo','0',5),('cm','cm','lo',1),('hi','hi','cm',1),('top',dd['reference'],'hi',11)]:
        c.resistor('ref'+name,a,b,r*cfg.get('ref_r',1e6),True)
    c.fet('cmpbias','p','pb','pb',dd['detector'],dd['detector'],.84,1)
    c.resistor('cmpbias','pb','0',cfg.get('cmp_r',1e7),True)
    c.fet('nbias','n','nb','nb','0','0',.42,1)
    c.resistor('nbias',dd['bias'],'nb',cfg.get('packet_r',1e8),True)
    c.fet('pbias','p','pp','pp',dd['bias'],dd['bias'],.42,1)
    c.resistor('pbias','pp','0',cfg.get('packet_r',1e8),True)
    c.source('startup','startup',pwl([(0,d.vdd),(800,d.vdd),(801,0)]),'initialization')
    c.inv('startupinv','startup','startupb',dd['logic'])
    for node in ('center','left','right'):
        c.capacitor(node)
        c.tg('init'+node,'cm',node,'startup','startupb',dd['logic'],w=.42)
    for sign,vin,ref in [('pos','center','hi'),('neg','lo','center')]:
        # Five-transistor continuously active detector, followed by 2 inverters.
        c.fet(sign+'tail','p',sign+'tail','pb',dd['detector'],dd['detector'],.84,1)
        c.fet(sign+'in1','p',sign+'m',vin,sign+'tail',dd['detector'],1,.5)
        c.fet(sign+'in2','p',sign+'o',ref,sign+'tail',dd['detector'],1,.5)
        c.fet(sign+'m1','n',sign+'m',sign+'m','0','0',.42,1)
        c.fet(sign+'m2','n',sign+'o',sign+'m','0','0',.42,1)
        c.inv(sign+'d1',sign+'o',sign+'db',dd['detector'])
        c.inv(sign+'d2',sign+'db',sign+'req',dd['detector'])
        fire=sign+'fire';fb=sign+'fb';reset=sign+'reset'
        set_signal=sign+'req'
        if cfg.get('guarded',False):
            # Four-phase guard: a held overflow request cannot fire twice.
            c.inv(sign+'qi',sign+'req',sign+'reqb',dd['logic'])
            nor(c,sign+'arm1',[fire,sign+'armedb'],sign+'armed',dd['logic'])
            nor(c,sign+'arm2',[sign+'reqb',sign+'armed'],sign+'armedb',dd['logic'])
            nor(c,sign+'setgate',[sign+'reqb',sign+'armedb'],sign+'set',dd['logic'])
            set_signal=sign+'set'
        nor(c,sign+'l1',[reset,fb,'startup'],fire,dd['logic'])
        nor(c,sign+'l2',[set_signal,fire],fb,dd['logic'])
        c.resistor(sign+'delay',fire,sign+'rc',cfg.get('reset_r',1e7),True)
        c.capacitor(sign+'rc',10,'mim')
        c.inv(sign+'r1',sign+'rc',sign+'rb',dd['logic'])
        c.inv(sign+'r2',sign+'rb',reset,dd['logic'])
        c.inv(sign+'fi',fire,sign+'fireb',dd['logic'])
        for node,mult in [('center',4),('left',1),('right',1)]:
            sink=(sign=='pos')==(node=='center');name=sign+node
            c.lines.append(f'V{name}sense {node} {name}port 0')
            c.save.add(f'i(v{name}sense)')
            port=name+'port'
            mult*=cfg.get('n_scale' if sink else 'p_scale',1)
            if node=='center':mult*=cfg.get('center_n_scale' if sink else 'center_p_scale',1)
            if sink:
                c.fet(name+'sw','n',port,fire,name+'s','0',.84*mult,.15)
                c.fet(name+'i','n',name+'s','nb','0','0',.42*mult,1)
            else:
                c.fet(name+'sw','p',port,sign+'fireb',name+'s',dd['transfer'],1.68*mult,.15)
                c.fet(name+'i','p',name+'s','pp',dd['transfer'],dd['transfer'],.42*mult,1)
    # Input charge is independent of all detector/latch signals.
    mode=cfg.get('mode','events');pulse=[(0,0)];events=[]
    if mode!='idle':
        for j,sgn in enumerate(cfg.get('signs',[1])):
            units=cfg.get('input_sequence',[cfg.get('input_units',4.4)]*4)[j]
            amp=cfg.get('input_a',1e-7)
            t=2500+j*1500;duration=units*d.cap_ff*1e-15*.025/amp*1e9
            pulse += [(t,0),(t+.1,sgn*amp),(t+duration,sgn*amp),(t+duration+.1,0)]
            events.append({'start_ns':t,'input_end_ns':t+duration+.1,'end_ns':t+1200,'sign':sgn})
    c.lines += ['Iinput 0 injected '+pwl(pulse),'Vinput injected center 0'];c.save.add('i(vinput)')
    stop=2500+max(1,len(events))*1500
    c.save.update(['v(lo)','v(hi)','v(cm)','v(center)','v(left)','v(right)'])
    keep=['lo','hi','cm','center','left','right','posfire','negfire','posreq','negreq','posreset','negreset','posrc','negrc','poso','nego']
    if cfg.get('guarded',False):keep+=['posarmed','negarmed']
    c.save={f'v({n})' for n in keep}|{'i(vinput)'}|{f'v({p["node"]})' for p in c.ports.values()}|set(c.ports)|{f'i(v{s}{n}sense)' for s in ('pos','neg') for n in ('center','left','right')}
    meta={'seed':d.seed,'design':asdict(d),'cfg':cfg,'mos':c.mos,'ports':c.ports,'events':events,'stop_ns':stop,
          'scope':'Physical detectors, SR latches, RC reset, gated weighted mirrors; ideal external additive current stimulus and startup command only.'}
    return c.finish(f'.tran .5n {stop}n 0 {cfg.get("maxstep_ns",.5)}n'),meta

def integrate(t,y,a,b):
    sel=(t>a)&(t<b)
    return float(np.trapezoid(np.r_[np.interp(a,t,y),y[sel],np.interp(b,t,y)],np.r_[a,t[sel],b]))

def analyze(result,names,x):
    meta=result['meta'];cfg=meta['cfg'];ix={v:j for j,v in enumerate(names)};t=x[:,0]
    if t[-1]<meta['stop_ns']*1e-9*(1-1e-6):raise ValueError('truncated transient')
    v=lambda n:x[:,ix['v('+n+')']]
    power={}
    for port,info in meta['ports'].items():
        if info['component']=='initialization':continue
        power[info['component']]=-v(info['node'])*x[:,ix[port]]
    energies=lambda a,b:{k:integrate(t,p,a*1e-9,b*1e-9)*1e15 for k,p in power.items()}
    idle=energies(1800,2300);idle_nw=sum(idle.values())/500*1000
    obs=[]
    for ev in meta['events']:
        a=ev['start_ns'];b=ev['end_ns'];lo=(a-5)*1e-9;hi=b*1e-9
        fires=[]
        for s in ('pos','neg'):
            sig=v(s+'fire');edges=np.where((sig[:-1]<cfg.get('vdd',1.8)/2)&(sig[1:]>=cfg.get('vdd',1.8)/2))[0]+1
            ends=np.where((sig[:-1]>=cfg.get('vdd',1.8)/2)&(sig[1:]<cfg.get('vdd',1.8)/2))[0]+1
            for i in edges:
                if lo<t[i]<hi:
                    e=next((j for j in ends if j>i),len(t)-1)
                    fires.append({'polarity':s,'time_ns':float(t[i]*1e9),'width_ns':float((t[e]-t[i])*1e9),
                                  'center_at_fire_v':float(v('center')[i])})
        delta={n:float((np.interp(hi,t,v(n))-np.interp(lo,t,v(n)))/.025) for n in ('center','left','right')}
        charge_unit=meta['design']['cap_ff']*1e-15*.025
        inj=integrate(t,x[:,ix['i(vinput)']],a*1e-9,b*1e-9)/charge_unit
        baseline_ns=ev['input_end_ns']+2
        input_finished_before_fire=bool(fires) and min(f['time_ns'] for f in fires)>baseline_ns+2
        packet={n:float((np.interp(hi,t,v(n))-np.interp(baseline_ns*1e-9,t,v(n)))/.025) for n in ('center','left','right')}
        target={'center':-4*ev['sign'],'left':ev['sign'],'right':ev['sign']}
        branch_charge={n:-sum(integrate(t,x[:,ix[f'i(v{s}{n}sense)']],a*1e-9,b*1e-9) for s in ('pos','neg'))/charge_unit for n in ('center','left','right')}
        error=max(abs(packet[n]-target[n]) for n in packet)
        charge_error=max(abs(branch_charge[n]-target[n]) for n in target)
        one_fire=input_finished_before_fire and len(fires)==1 and fires[0]['polarity']==('pos' if ev['sign']>0 else 'neg')
        energy=energies(a,b);total=sum(energy.values())
        gross={k:integrate(t,np.maximum(p,0),a*1e-9,b*1e-9)*1e15 for k,p in power.items()}
        phases={};rearm_ns=None;settle_ns=None
        if len(fires)==1:
            f=fires[0];s=f['polarity'];rise=f['time_ns'];fall=rise+f['width_ns']
            ready=(t>fall*1e-9)&(v(s+'req')<.1*cfg.get('vdd',1.8))&(v(s+'reset')<.1*cfg.get('vdd',1.8))&(v(s+'rc')<.1*cfg.get('vdd',1.8))
            if cfg.get('guarded',False):ready&=v(s+'armed')>.9*cfg.get('vdd',1.8)
            candidates=np.flatnonzero(ready&(t<hi))
            if len(candidates):
                rearm_ns=float(t[candidates[0]]*1e9)
                phases={'before_firing_fj':energies(a,rise),'active_transfer_fj':energies(rise,fall),
                        'reset_and_rearm_fj':energies(fall,rearm_ns),'after_rearm_fj':energies(rearm_ns,b)}
            err=np.maximum.reduce([abs(v(n)-np.interp(hi,t,v(n)))/.025 for n in ('center','left','right')])
            outside=np.flatnonzero((t>=rise*1e-9)&(t<hi)&(err>.1))
            settle_ns=float(t[outside[-1]+1]*1e9) if len(outside) else rise
        obs.append({**ev,'firings':fires,'input_charge_units':inj,'node_delta_units':delta,'packet_units_apparent':packet,
                    'baseline_ns':baseline_ns,'baseline_center_v':float(np.interp(baseline_ns*1e-9,t,v('center'))),
                    'input_finished_before_fire':input_finished_before_fire,
                    'branch_charge_units':branch_charge,'max_charge_error_units':charge_error,
                    'charge_passed':one_fire and charge_error<.1,
                    'max_packet_error_units':error,'passed':one_fire and error<.1,
                    'rearm_time_ns':rearm_ns,'settled_time_ns':settle_ns,'phase_supply_energy_fj':phases,
                    'energy_fj':energy,'gross_supply_fj':gross,'total_supply_fj':total,
                    'external_input_work_fj':integrate(t,v('center')*x[:,ix['i(vinput)']],a*1e-9,b*1e-9)*1e15,
                    'incremental_over_local_idle_fj':total-sum(idle.values())*(b-a)/500})
    pulses=sum(len(o['firings']) for o in obs)
    return {'cfg':cfg,'idle_power_nw':idle_nw,'idle_500ns_energy_fj':idle,'observations':obs,
            'charge_passed':bool(obs) and all(o['charge_passed'] for o in obs),
            'firings':pulses,'passed':bool(obs) and all(o['passed'] for o in obs),
            'worst_packet_error_units':max([o['max_packet_error_units'] for o in obs],default=0),
            'event_window_mean_fj':float(np.mean([o['total_supply_fj'] for o in obs])) if obs else None,
            'incremental_mean_fj':float(np.mean([o['incremental_over_local_idle_fj'] for o in obs])) if obs else None,
            'startup_supply_fj':energies(0,1800),
            'startup_command_gross_fj':integrate(t,np.maximum(-v('startup')*x[:,ix['i(vstartup)']],0),0,1800e-9)*1e15,
            'mos_count':len(meta['mos']),'stop_ns':meta['stop_ns'],
            'monitored_state_min_v':min(float(np.min(v(n))) for n in ('center','left','right')),
            'monitored_state_max_v':max(float(np.max(v(n))) for n in ('center','left','right'))}

def run(cfg):
    label=cfg.get('label') or hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest()[:12]
    directory=OUT/'spice'/label
    directory.mkdir(parents=True,exist_ok=True)
    try:
        r,n,x=execute(*build(cfg),directory,timeout=150);row=analyze(r,n,x)
        (directory/'metrics.json').write_text(json.dumps(row,indent=2));compress_trace(directory/'output.raw')
        print(label,row['passed'],row['firings'],round(row['worst_packet_error_units'],3),row['event_window_mean_fj'],flush=True)
        return {**row,'run':str(directory.relative_to(ROOT))}
    except Exception as e:
        row={'cfg':cfg,'passed':False,'unresolved':str(e),'run':str(directory.relative_to(ROOT))}
        (directory/'unresolved.json').write_text(json.dumps(row,indent=2));print(label,'UNRESOLVED',str(e)[-200:],flush=True);return row

def characterize(cfg):
    """Small declared robustness grid, held fixed after nominal tuning."""
    from gamma_sky_batch import batch_execute
    cases=[];configs=[]
    pvt=[{'corner':k,'temp':25,'vdd':1.8} for k in ('tt','ss','ff','sf','fs')]
    pvt += [{'corner':'tt','temp':t,'vdd':1.8} for t in (-40,100)]
    pvt += [{'corner':'tt','temp':25,'vdd':v} for v in (1.62,1.9)]
    for i,item in enumerate(pvt):
        for s in (1,-1):configs.append({**cfg,**item,'signs':[s],'label':f'pvt_{i}_s{s}','group':'pvt'})
    for seed in range(2001,2009):
        for s in (1,-1):configs.append({**cfg,'mismatch':True,'seed':seed,'signs':[s],'label':f'mc_{seed}_s{s}','group':'mismatch'})
    for s in (1,-1):
        for units in (3.2,4.8):configs.append({**cfg,'input_units':units,'signs':[s],'label':f'hold_s{s}_u{units}','group':'threshold'})
        configs.append({**cfg,'signs':[s,s],'input_sequence':[6,4],'label':f'repeat_s{s}','group':'repeat'})
        configs.append({**cfg,'signs':[s],'maxstep_ns':.1,'label':f'fine_s{s}','group':'timestep'})
    configs.append({**cfg,'mode':'idle','label':'idle','group':'idle'})
    for item in configs:
        deck,meta=build(item);cases.append((deck,meta,OUT/'spice'/item['label']))
    def process(result,names,x,directory):
        row=analyze(result,names,x);row['run']=str(directory.relative_to(ROOT))
        (directory/'metrics.json').write_text(json.dumps(row,indent=2));compress_trace(directory/'output.raw')
        print(directory.name,'Qpass',row['charge_passed'],'Vpass',row['passed'],'fires',row['firings'],flush=True)
        return {**result,'observations':row['observations'],'summary':row}
    results=batch_execute(cases,workers=2,processor=process)
    rows=[r['summary'] for r in results]
    (OUT/'characterization.json').write_text(json.dumps(rows,indent=2))
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--sweep',action='store_true');p.add_argument('--characterize',action='store_true');p.add_argument('--initial',action='store_true');args=p.parse_args()
    if args.characterize:
        characterize(SELECTED);return
    configs=[{**SELECTED,'label':f'selected_s{s}','signs':[s]} for s in (1,-1)]
    if args.initial:configs=[{'label':'initial','reset_r':1e7}]
    if args.sweep:configs=[{'label':f'packet{r:g}_reset{s:g}','packet_r':r,'reset_r':s} for r,s in itertools.product([3e7,1e8,3e8],[3e6,1e7,3e7])]
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(run,configs))
    OUT.mkdir(parents=True,exist_ok=True);(OUT/('sweep.json' if args.sweep else 'initial.json' if args.initial else 'selected.json')).write_text(json.dumps(rows,indent=2))

if __name__=='__main__':main()
