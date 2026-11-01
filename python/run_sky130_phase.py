"""Narrow q=4 SKY130 experiment, with explicit measured/unbuilt energy scope."""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json,math,random,statistics,time
from collections import defaultdict
from dataclasses import asdict,replace
from pathlib import Path
import numpy as np
from gamma_sky130 import *
from gamma_sky_batch import batch_execute
from gamma_normalizer_statistics import confidence

OUT=BASE/'results'


def write_csv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def summarize(result,path,group):
    o=result['observations'];d=result['meta']['design']
    row={'group':group,'run':str(path.relative_to(ROOT)),**d,
        'simulation_status':result.get('simulation_status','complete'),
        'unresolved_checkpoints':len(result['meta']['checkpoints'])-len(o),
        'mode':result['meta']['mode'],'residual_input_units':result['meta']['residual_units_in_stage1'],
        'positions_tested':len(o),'wrong_next_mac':sum(not p['next_mac_correct'] for p in o),
        'wrong_current_class':sum(p['current_class']!=p['digit'] for p in o),
        'mean_energy_fj':result['mean_energy_fj'],'worst_energy_fj':max(p['measured_known_path_fj'] for p in o),
        'max_abs_error_units':result['max_abs_driver_error_units'],
        'max_storage_error_units':max(abs(p['stored_error_units']) for p in o),
        'max_retention_drift_units':max(abs(p['retention_drift_units']) for p in o),
        'max_reference_droop_units':max(p['peak_reference_droop_units'] for p in o),
        'max_settling_ns':max((p['settling_ns_to_point1digit'] or 1e9) for p in o),
        'cycles_not_settled_to_point1digit':sum(p['settling_ns_to_point1digit'] is None for p in o),
        'max_observed_settling_ns':max((p['settling_ns_to_point1digit'] for p in o if p['settling_ns_to_point1digit'] is not None),default=None),
        'soa_violations':len(result['device_operating_limit_violations'])}
    for key in sorted(set().union(*(p['energy_components_fj'] for p in o))):
        row['energy_'+key+'_fj']=statistics.mean(p['energy_components_fj'].get(key,0) for p in o)
    row['feasible_nominal']=row['wrong_next_mac']==0 and row['unresolved_checkpoints']==0 and row['soa_violations']==0 and row['max_abs_error_units']<.20
    return row


def run_boundaries(jobs,workers=3):
    cases=[(*boundary_deck(d,sequence,residual,direct),path) for group,d,path,sequence,residual,direct in jobs]
    def process(r,names,x,path):
        meta=r['meta']
        if r.get('simulation_status')=='incomplete_transient':
            # Only full cycles have a valid per-operation energy integral.
            checks=[p for p in meta['checkpoints'] if (p['start_ns']+meta['design']['period_ns'])*1e-9<=x[-1,0]]
            if not checks:raise RuntimeError('No completed cycles before nonconvergence: '+str(path))
            r['meta']={**meta,'checkpoints':checks}
        r=evaluate_boundary(r,names,x);r['meta']=meta
        if r.get('simulation_status')=='incomplete_transient':r['all_next_mac_correct']=False
        dump(path/'metrics.json',r)
        # Checkpoints, reference levels, integrated supplies and limit checks
        # survive for every case. Only representative full traces are needed.
        keep=r.get('simulation_status')=='incomplete_transient' or path.name in ('s000','s014','penergy004','pmargin004','l094','l095')
        raw=path/'output.raw'
        if keep:compress_trace(raw);r['waveform_storage']='lossless gzip'
        else:
            if raw.exists():raw.unlink()
            gz=Path(str(raw)+'.gz')
            if gz.exists():gz.unlink()
            r['waveform_storage']='metrics/checkpoints retained; standalone netlist regenerates full trace'
        return r
    results=batch_execute(cases,workers,processor=process)
    summaries=[]
    for r,(deck,meta,path),(group,*_) in zip(results,cases,jobs):
        if not r.get('observations'):
            names,x=fast_raw(path/'output.raw');r=process(r,names,x,path)
        dump(path/'metrics.json',r);summaries.append(summarize(r,path,group))
    return summaries


def sizing_jobs():
    designs=[];base=Design()
    for cap,tg,scale in itertools.product((25,50,100,200),(.42,.84,1.68),(.5,1,2)):
        designs.append(replace(base,cap_ff=cap,tg_w=tg,buffer_scale=scale))
    for iw,ratio in itertools.product((1,4,12),(1,2,4)):
        designs.append(replace(base,cap_ff=50,tg_w=.42,input_w=iw,pn_ratio=ratio))
    for bias,load,period in itertools.product((2e6,1e7,4e7),(2e5,1e6),(100,400)):
        designs.append(replace(base,cap_ff=50,tg_w=.42,bias_r=bias,load_r=load,period_ns=period,gated_buffer=True))
    for ladder,buffered,period in itertools.product((5e5,2e6,1e7),(False,True),(100,400)):
        designs.append(replace(base,cap_ff=50,tg_w=.42,bias_r=1e7,load_r=1e6,ladder_r=ladder,
            period_ns=period,gated_buffer=True,write_buffer=buffered))
    for bias,iw in itertools.product((2e5,2e6,1e7),(4,16)):
        designs.append(replace(base,cap_ff=50,tg_w=.42,buffer_kind='five_transistor',input_w=iw,bias_r=bias,load_r=1e6))
    for cap,period in itertools.product((25,50,100),(100,400)):
        designs.append(replace(base,cap_ff=cap,tg_w=.42,selector_kind='nmos',bias_r=1e7,load_r=1e6,
            ladder_r=2e6,period_ns=period,gated_buffer=True))
    # Clock-gating can lose on wake-up settling. Test continuous low-bias
    # operation too; do not reject low-power cells using only a gated driver.
    designs.extend(replace(d,gated_buffer=False) for d in list(designs) if d.gated_buffer)
    # Coupled charge-injection/RC tradeoffs: smaller C needs smaller switches
    # and the appropriate ladder/write time, not an isolated C shrink.
    for cap,ratio,period,gate,selector in itertools.product((25,50),(2,4),(100,400),(False,True),('tg','nmos')):
        designs.append(replace(base,cap_ff=cap,tg_w=.42,pn_ratio=ratio,bias_r=1e7,load_r=1e6,
            period_ns=period,ladder_r=1e6*25/cap*period/100,gated_buffer=gate,selector_kind=selector))
    # Compact-class margin / device-width axis rather than q/base exploration.
    unique={json.dumps(asdict(d),sort_keys=True):d for d in designs}
    return [('sizing',d,BASE/'regen_cell'/f's{i:03d}',DEFAULT_SEQUENCE,.35,False) for i,d in enumerate(unique.values())]


def comparator_deck(d=Design(),deltas=None):
    deltas=deltas or [-.020,-.010,-.00625,-1/240,-.002,-.001,-.00025,0,.00025,.001,.002,1/240,.00625,.010,.020]
    c=Circuit(d);c.supply('comparator');cm=d.vdd/3;uin=d.vdd/72*4/d.input_scale
    ref=cm+.5*uin;period=30;points=[(0,ref+deltas[0])];clock=[(0,0)];checks=[]
    for j,delta in enumerate(deltas):
        s=j*period
        points += [(s+.1,ref+delta),(s+period-.1,ref+delta)]
        clock += [(s+6,0),(s+6.1,d.vdd),(s+24-.1,d.vdd),(s+24,0)]
        checks.append({'delta_v':delta,'start_ns':s,'eval_start_ns':s+6.1,'time_ns':s+21})
    c.source('input','input',pwl(points),'external_analog_input')
    c.source('reference','reference',f'{ref:.12g}','external_reference')
    c.source('clock','clock',pwl(clock),'clock_input')
    c.comparator('input','reference','clock')
    meta={'kind':'comparator','design':asdict(d),'seed':d.seed,'ports':c.ports,'mos':c.mos,
        'checkpoints':checks,'stop_ns':len(deltas)*period,'unit_v':uin}
    return c.finish(f'.tran .05n {len(deltas)*period}n'),meta


def run_comparators(designs,group,workers=3):
    cases=[(*comparator_deck(d),BASE/'comparator'/f'{group}{i:03d}') for i,d in enumerate(designs)]
    results=batch_execute(cases,workers);rows=[]
    for r,(_,meta,path) in zip(results,cases):
        names,x=fast_raw(path/'output.raw');ix={n:i for i,n in enumerate(names)};t=x[:,0];d=meta['design'];obs=[]
        vp=x[:,ix['v(cmpop)']];vn=x[:,ix['v(cmpon)']]
        for p in meta['checkpoints']:
            i=np.searchsorted(t,p['time_ns']*1e-9);sign=p['delta_v']>0
            resolved=bool((vp[i]>.8*d['vdd'] and vn[i]<.2*d['vdd']) or (vn[i]>.8*d['vdd'] and vp[i]<.2*d['vdd']))
            # Refer decision delay to the 50% point of the 0.1 ns clock edge.
            clock_half_ns=p['eval_start_ns']-.05
            bit=bool(vp[i]>vn[i]);active=(t>=clock_half_ns*1e-9)&(t<p['time_ns']*1e-9)
            settled=active&((vp>.8*d['vdd'])&(vn<.2*d['vdd']) if bit else (vn>.8*d['vdd'])&(vp<.2*d['vdd']))
            delay=float((t[np.where(settled)[0][0]]-clock_half_ns*1e-9)*1e9) if np.any(settled) else None
            cycle=(t>=p['start_ns']*1e-9)&(t<=(p['start_ns']+30)*1e-9)
            energy=0.;supply=0.;pin=0.
            for name,port in meta['ports'].items():
                current=x[:,ix[name]];voltage=x[:,ix['v('+port['node']+')']]
                if port['component']=='comparator':
                    supply+=float(np.trapezoid(np.maximum(-voltage[cycle]*current[cycle],0),t[cycle]))*1e15
                elif port['component']=='clock_input':
                    pin+=d['vdd']*float(np.trapezoid(np.maximum(-current[cycle],0),t[cycle]))*1e15
            energy=supply+pin
            obs.append(dict(p,resolved=resolved,positive=bit,correct=(None if p['delta_v']==0 else bit==sign and resolved),
                decision_ns=delay,comparator_vdd_delivered_fj=supply,clock_pin_charging_fj=pin,energy_fj=energy,
                output_positive_v=float(vp[i]),output_negative_v=float(vn[i])))
        neg=[o['delta_v'] for o in obs if o['resolved'] and not o['positive']]
        pos=[o['delta_v'] for o in obs if o['resolved'] and o['positive']]
        monotonic=bool(neg and pos and max(neg)<min(pos))
        offset=(max(neg)+min(pos))/2 if monotonic else None
        atmargin=[o for o in obs if math.isclose(abs(o['delta_v']),1/240,abs_tol=1e-9)]
        row={'group':group,'run':str(path.relative_to(ROOT)),**d,'offset_bracket_mid_v':offset,
            'offset_bracket_low_v':max(neg) if monotonic else None,'offset_bracket_high_v':min(pos) if monotonic else None,
            'monotonic_resolved_transfer':monotonic,
            'legal_center_wrong':sum(not o['correct'] for o in atmargin),'legal_center_tests':len(atmargin),
            'mean_energy_fj':statistics.mean(o['energy_fj'] for o in obs),
            'max_legal_center_energy_fj':max(o['energy_fj'] for o in atmargin),
            'max_decision_ns':max((o['decision_ns'] if o['decision_ns'] is not None else 15) for o in obs if o['delta_v']!=0),
            'unresolved_nonzero':sum(not o['resolved'] for o in obs if o['delta_v']!=0)}
        r['observations']=obs;r['summary']=row;dump(path/'metrics.json',r);rows.append(row)
    return rows


def pareto(rows):
    feasible=[r for r in rows if r['feasible_nominal']]
    keys=('mean_energy_fj','max_abs_error_units','max_settling_ns')
    for r in rows:
        r['pareto']=r in feasible and not any(all(s[k]<=r[k] for k in keys) and any(s[k]<r[k] for k in keys) for s in feasible)
    return feasible


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=2);parser.add_argument('--mc-trials',type=int,default=16)
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    smokes=json.loads((OUT/'smoke_summary.json').read_text());assert len(smokes)==20 and all(r['passed'] for r in smokes)
    start=time.monotonic();jobs=sizing_jobs();sizing=run_boundaries(jobs,args.workers);valid=pareto(sizing)
    write_csv(OUT/'sizing_sweep.csv',sizing);dump(OUT/'sizing_sweep.json',sizing)
    if not valid:raise RuntimeError('No nominal design meets declared error/operating limits; inspect sizing_sweep rather than proceed silently')
    best=min(valid,key=lambda r:r['mean_energy_fj']);d0=Design(**{k:best[k] for k in asdict(Design())})
    margin=min(valid,key=lambda r:r['max_abs_error_units'])
    robust=Design(**{k:margin[k] for k in asdict(Design())})
    dump(OUT/'selected_design.json',{'selection':'minimum measured energy with all12 next-MAC outputs correct, <=.20 digit residual and valid device limits',
        'design':asdict(d0),'sizing_record':best,'margin_design':asdict(robust),'margin_record':margin})
    print('selected',best['run'],best['mean_energy_fj'],'fJ',flush=True)
    comp_designs=[replace(d0,corner=corner,comparator_w=w,comparator_l=l) for corner,w,l in itertools.product(CORNERS,(1,4,16),(.15,.5))]
    comp=run_comparators(comp_designs,'c',args.workers);write_csv(OUT/'comparator_summary.csv',comp)
    options=[r for r in comp if r['corner']=='tt' and r['legal_center_wrong']==0]
    compbest=min(options,key=lambda r:r['max_legal_center_energy_fj']) if options else min(comp,key=lambda r:r['legal_center_wrong'])
    dh=replace(d0,detector='hybrid',comparator_w=compbest['comparator_w'],comparator_l=compbest['comparator_l'])
    dr=replace(robust,detector='hybrid',comparator_w=16,comparator_l=.5)
    # Frozen choice: no post-PVT retuning or failure deletion.
    follow=[]
    for candidate,base in [('energy',dh),('margin',dr)]:
        for index,(corner,voltage,temp) in enumerate(itertools.product(CORNERS,(1.62,1.8,1.9),(-40,25,100))):
            d=replace(base,corner=corner,vdd=voltage,temp_c=temp,cap_kind='mim')
            follow.append(('pvt_'+candidate,d,BASE/'pvt'/f'p{candidate}{index:03d}',DEFAULT_SEQUENCE,.1,False))
    for kind in ('ideal','mim'):
        for residual in (0,.1,.35):
            for direct in (False,True):
                d=replace(dh,cap_kind=kind)
                follow.append(('chaining',d,BASE/'loaded_mac'/f'l{len(follow):03d}',DEFAULT_SEQUENCE,residual,direct))
    # Common-scale control removes the direct path's additional gain3/VCM
    # loading. Isolate actual writer/driver restoration with latched ideal
    # decisions; hybrid S12 controls above still retain the physical threshold.
    for residual,direct in itertools.product((0,.1,.35),(False,True)):
        d=replace(d0,cap_kind='mim',input_scale=4,detector='ideal')
        follow.append(('chaining_matched_scale',d,BASE/'loaded_mac'/f'matched{residual:g}_{int(direct)}',DEFAULT_SEQUENCE,residual,direct))
    for i in range(args.mc_trials):
        seed=81000+i
        for mode in ('pdk','architecture','margin'):
            rng=random.Random(seed)
            d=replace(dr if mode=='margin' else dh,seed=seed,cap_kind='mim' if mode!='architecture' else 'ideal',pdk_mismatch=mode!='architecture',
                resistor_sigma=.005 if mode=='architecture' else 0,capacitor_sigma=.02 if mode=='architecture' else 0,
                systematic_offset_mv=rng.gauss(0,.5) if mode=='architecture' else 0)
            follow.append(('mc_'+mode,d,BASE/'monte_carlo'/f'm{mode[0]}{i:03d}',DEFAULT_SEQUENCE,.1,False))
    # All49 transition pairs, each actual previously-written state retained.
    sequence=[]
    for a in range(-3,4):
        for b in range(-3,4):sequence.extend((a,b))
    follow.append(('transitions',dh,BASE/'regen_cell/transition_matrix',tuple(sequence),.1,False))
    followrows=run_boundaries(follow,args.workers)
    pvt=[r for r in followrows if r['group'].startswith('pvt_')];write_csv(OUT/'pvt_summary.csv',pvt)
    mc=[r for r in followrows if r['group'].startswith('mc_')];write_csv(OUT/'monte_carlo_trials.csv',mc)
    mcsummary=[]
    for mode in ('mc_pdk','mc_architecture','mc_margin'):
        group=[r for r in mc if r['group']==mode];fail=sum(r['wrong_next_mac']>0 or r['unresolved_checkpoints']>0 for r in group)
        mcsummary.append({'mode':mode,'network_trials':len(group),'networks_with_any_wrong_class':fail,
            'wrong_checkpoints':sum(r['wrong_next_mac'] for r in group),'checkpoints':sum(r['positions_tested'] for r in group),
            'unresolved_checkpoints':sum(r['unresolved_checkpoints'] for r in group),
            **confidence(fail,len(group))})
    write_csv(OUT/'monte_carlo_summary.csv',mcsummary)
    transition_path=BASE/'regen_cell/transition_matrix/metrics.json';transition=json.loads(transition_path.read_text())
    matrix=[]
    # Explicit (a,b) evaluations are the odd entries; intervening writes prepare a.
    for obs in transition['observations'][1::2]:matrix.append({'from':obs['previous_digit'],'to':obs['digit'],
        'energy_fj':obs['measured_known_path_fj'],'error_units':obs['driver_error_units'],'next_mac_correct':obs['next_mac_correct'],
        'settling_ns':obs['settling_ns_to_point1digit']})
    write_csv(OUT/'regen_transition_matrix.csv',matrix)
    compmc=run_comparators([replace(base,pdk_mismatch=True,seed=91000+i) for base in (dh,dr) for i in range(args.mc_trials)],'m',args.workers)
    write_csv(OUT/'comparator_monte_carlo.csv',compmc)
    allrows=sizing+followrows
    breakdown=[]
    for r in allrows:
        for k,v in r.items():
            if k.startswith('energy_') and k.endswith('_fj'):breakdown.append({'run':r['run'],'group':r['group'],'component':k[7:-3],'energy_fj_mean':v})
    write_csv(OUT/'energy_breakdown.csv',breakdown)
    summary={'q':4,'simulator':'ngspice47+SKY130A','sizing_cases':len(sizing),'nominal_feasible_cases':len(valid),
        'selected_design':asdict(d0),'hybrid_design':asdict(dh),'robust_hybrid_design':asdict(dr),
        'selected_measured_known_path_fj':best['mean_energy_fj'],
        'selected_comparator':compbest,'pvt_cases':len(pvt),'pvt_all_states_correct':sum(r['wrong_next_mac']==0 and r['unresolved_checkpoints']==0 for r in pvt),
        'incomplete_transients':[r['run'] for r in followrows if r['simulation_status']!='complete'],
        'pvt_valid_operating_limits':sum(r['soa_violations']==0 and r['unresolved_checkpoints']==0 for r in pvt),'monte_carlo':mcsummary,
        'energy_gate_ratios':{str(e):best['mean_energy_fj']/e for e in (1,2.74,5,10,20)},
        'sizing':sizing,'followup':followrows,'comparators':comp,'comparator_mc':compmc,
        'elapsed_seconds':time.monotonic()-start,'scope':'Measured VDD write/driver/reference/interface plus one comparator in hybrid. Five decisions/decoder/sequencer remain ideal. NOT full-boundary measured energy.',
        'temperature_c':[-40,25,100],'supply_v':[1.62,1.8,1.9],'supply_1p98_excluded':True,
        'unbuilt':['remaining region decisions','class latch/decoder logic','autonomous completion','reference distribution layout','full transistor normalizer'],
        'next_step_gate':'No crossbar until supply-priced boundary meets a matched baseline; evaluate Pareto and explicit unknown costs.'}
    dump(OUT/'sky130_phase_summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('sizing','followup','comparators','comparator_mc')},indent=2),flush=True)


if __name__=='__main__':main()
