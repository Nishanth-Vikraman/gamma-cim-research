"""Signed binary KCL / downward residual physical falsification, independent of q4.
Only existing PDK MOS primitives and literal-deck runner are reused.
"""
from __future__ import annotations
import argparse,csv,itertools,json,math
from fractions import Fraction as F
from dataclasses import asdict
from signed_binary_pdk import selected,focused
from pathlib import Path
import numpy as np
from gamma_sky130 import Circuit,Design,ROOT,pwl,dump,compress_trace
import signed_binary_backend as backend
OUT=ROOT/'results/signed_binary_kcl_20261026'
NS=(64,128,256,512,1024,2048)
KB=1.380649e-23

def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)

def workload(n,pattern,seed=20261026,bits=6):
    rng=np.random.default_rng(seed)
    if pattern=='zero_weights':
        a=np.ones((n,bits),dtype=np.int64);b=np.zeros_like(a)
    elif pattern.startswith('aligned'):
        a=np.ones((n,bits),dtype=np.int64);b=a.copy()
        if pattern=='aligned_negative':a=-a
    elif pattern in ('dense','sparse'):
        probs=None if pattern=='dense' else [.05,.9,.05]
        a=rng.choice([-1,0,1],(n,bits),p=probs);b=rng.choice([-1,0,1],(n,bits),p=probs)
    elif pattern=='cancellation':
        a=np.ones((n,bits),dtype=np.int64);a[::2]*=-1;b=np.ones_like(a)
    elif pattern in ('int8','repository_int8'):
        if pattern=='repository_int8':
            prior=json.loads((ROOT/'results/sar_normalizer_phase/exact_workload.json').read_text())
            pairs=prior['pairs'][:n];av=np.array([v[0] for v in pairs]);bv=np.array([v[1] for v in pairs]);bits=8
        else:av=rng.integers(-127,128,n);bv=rng.integers(-127,128,n);bits=7
        a=np.array([[abs(int(v))>>j&1 for j in range(bits)] for v in av])*np.sign(av)[:,None]
        b=np.array([[abs(int(v))>>j&1 for j in range(bits)] for v in bv])*np.sign(bv)[:,None]
    else:raise ValueError(pattern)
    counts=[]
    for k in range(2*bits-1):
        pairs=[(i,k-i) for i in range(bits) if 0<=k-i<bits]
        prod=np.concatenate([a[:,i]*b[:,j] for i,j in pairs])
        programmed=np.concatenate([b[:,j]!=0 for i,j in pairs])
        counts.append({'k':k,'terms':len(pairs)*n,'programmed':int(sum(programmed)),
            'positive':int(sum(prod==1)),'negative':int(sum(prod==-1)),
            'zero':int(sum((prod==0)&programmed)),'c':int(sum(prod))})
    d=1<<int((2*bits*n-1).bit_length())
    exact=sum(x['c']*(1<<x['k']) for x in counts)
    assert exact==sum(int(x)*int(y) for x,y in zip(a@2**np.arange(bits),b@2**np.arange(bits)))
    assert d>=2*max(x['terms'] for x in counts)
    return {'n':n,'pattern':pattern,'seed':seed,'bits':bits,'D':d,'Q':d.bit_length()-1,
        'counts':counts,'dot':exact,'product_planes':bits*bits,
        'physical_product_connections':2*sum(x['programmed'] for x in counts),
        'worst_case_product_connections':2*n*bits*bits}

def golden(w,precision=12,guards=4):
    top=2*w['bits']-1;positions=list(range(top,top-precision-guards,-1))
    coeff={v['k']:v['c'] for v in w['counts']};r=F(0);rows=[];value=F(0)
    for k in positions:
        z=F(coeff.get(k,0),w['D']);t=z+2*r
        digit=1 if t>F(1,2) else -1 if t<F(-1,2) else 0
        r=t-digit;value+=w['D']*digit*F(2)**k
        assert abs(r)<=F(1,2)
        rows.append({'k':k,'z':float(z),'t':float(t),'digit':digit,'r':float(r),
            'margin':float(min(abs(t-F(1,2)),abs(t+F(1,2))))})
    error=F(w['dot'])-value
    discarded=sum(F(v['c'])*F(2)**v['k'] for v in w['counts'] if v['k']<positions[-1])
    discarded_bound=sum(F(v['terms'])*F(2)**v['k'] for v in w['counts'] if v['k']<positions[-1])
    assert error==w['D']*r*F(2)**positions[-1]+discarded
    return {'stages':rows,'decoded':float(value),'error':float(error),
        'max_truncation_error':float(w['D']*F(2)**positions[-1]/2+discarded_bound),
        'exact_guards_sufficient':w['Q']+max(0,2*w['bits']-precision)}

class BinaryCircuit(Circuit):
    def __init__(self,cfg):
        d=Design(corner=cfg.get('corner','tt'),vdd=min(cfg.get('vdd',1.8),1.95),
            temp_c=cfg.get('temp',25),comparator_w=cfg.get('cmp_w',.42),
            comparator_l=cfg.get('cmp_l',.15),capacitor_sigma=cfg.get('cap_sigma',0),
            pdk_mismatch=cfg.get('pdk_mm',False),seed=cfg.get('seed',20261026),cap_kind='mim')
        super().__init__(d);self.d.vdd=cfg.get('vdd',1.8);self.cfg=cfg
        self.lines[0]='Signed binary passive KCL and physical downward residual chain'
        self.lines.append('.option klu noinit reltol=1e-3 abstol=1e-12 vntol=1e-6')
        if not cfg.get('full_library',False):self.lines[1]='.include "'+focused(self.d.corner)+'"'
        self.caps=[];self.res=[];self.stages=[];self.selected_models=set();self.instance_multiplier=1
        for comp in ('reference','row','row_clock','clock','comparator','logic','subtract'):self.supply(comp)
        self.source('one','one',str(self.d.vdd),'constant_pin');self.source('zero','zero','0','constant_pin')
    def fet(self,name,kind,drain,gate,source,body,w,l=.15):
        super().fet(name,kind,drain,gate,source,body,w,l)
        if self.cfg.get('native_models',False):return
        model,path=selected(kind,self.d.corner,w,l)
        mult=self.instance_multiplier
        if mult!=1:self.lines[-1]+=' m='+str(mult)+' mult='+str(mult)
        self.mos[-1]['multiplicity']=mult
        self.lines[-1]=self.lines[-1].replace('sky130_fd_pr__nfet_01v8' if kind=='n' else 'sky130_fd_pr__pfet_01v8',model)
        if model not in self.selected_models:
            self.lines.insert(2,'.include "'+path+'"');self.selected_models.add(model)
    def invnode(self,node,tag,comp='logic',w=None):
        self.inv(tag,node,tag+'b','dd_'+comp,w=w or self.cfg.get('logic_w',.42));return tag+'b'
    def clock(self,tag,a,b):
        edge=min(.08,(b-a)/4)
        row=tag=='roweval';comp='row_clock' if row else 'clock'
        self.source(tag+'pin',tag+'pin',pwl([(0,0),(a,0),(a+edge,self.d.vdd),(b-edge,self.d.vdd),(b,0)]),'row_clock_pin' if row else 'clock_pin')
        if row:self.instance_multiplier=max(1,math.ceil(math.sqrt(6*self.cfg.get('n',128))))
        nb=self.invnode(tag+'pin',tag+'drv1',comp)
        if row:self.instance_multiplier=6*self.cfg.get('n',128)
        out=self.invnode(nb,tag+'drv2',comp,self.cfg.get('driver_w',.84))
        self.instance_multiplier=1;return out
    def switch(self,tag,a,b,en,comp='logic',w=None):
        eb=self.invnode(en,tag+'eb',comp)
        self.tg(tag,a,b,en,eb,'dd_'+comp,w=w or self.cfg.get('switch_w',.42))
    def cap(self,tag,a,b,ff):
        actual=ff*(1+self.rng.gauss(0,self.d.capacitor_sigma))
        side=max(2,(-.76+math.sqrt(.76**2+8*actual))/4+.025)
        actual=2*(side-.025)**2+.76*(side-.025)
        self.lines.append(f'Xcap{tag} {a} {b} sky130_fd_pr__cap_mim_m3_1 w={side:.12g} l={side:.12g}')
        self.caps.append({'tag':tag,'a':a,'b':b,'ff':actual,'area_um2':side**2})
        self.save.update((f'v({a})',f'v({b})'))
    def poly(self,tag,a,b,unit_r,m=1,kind='array',sigma=0):
        if m<=0:return
        # Native PDK poly model, micrometer dimensions under global scale=1u.
        # Explicit sensitivity is NOT a foundry-calibrated poly mismatch model.
        factor=1+self.rng.gauss(0,sigma/math.sqrt(m))
        width=.35;length=unit_r*width/48.2*factor
        high=self.cfg.get('res_model')=='high' and kind in ('array','ballast')
        if high:length=max(.5,(unit_r-1179.98)*.35/1112.41)*factor
        self.lines.append(f'Vmeasure{tag} {a} {tag}port 0')
        if high:
            self.lines.append(f'Xpoly{tag} {tag}port {b} 0 sky130_fd_pr__res_high_po_0p35 w={width:.12g} l={length:.12g} m={m} mult={m}')
            # PDK voltage-dependent resistor internals are simulator B/R objects;
            # zero-volt port metering includes their physical terminal work.
        else:self.lines.append(f'R{tag} {tag}port {b} sky130_fd_pr__res_generic_po w={width:.12g} l={length:.12g} m={m}')
        # PDK poly-to-substrate capacitance, split between both terminals.
        # Bare resistor cards alone contain no geometry parasitic capacitance.
        if self.cfg.get('poly_parasitics',True) and not high:
            for suffix,node in (('a',a),('b',b)):
                self.lines.append(f'Cpoly{tag}{suffix} {node} 0 mcp1f w={width:.12g} l={length:.12g} m={m/2:.12g}')
        self.save.add(f'i(vmeasure{tag.lower()})');self.save.update((f'v({a})',f'v({b})'))
        self.res.append({'tag':tag,'a':a,'b':b,'r_design':unit_r/factor,'m':m,'kind':kind,
            'model':'high' if high else 'generic','w_um':width,'l_um':length,'area_um2':width*(length+(4.16 if high else 0))*m})
    def latch(self,data,clock,tag,reset):
        cb=self.invnode(clock,tag+'cb');h=tag+'h';q=tag+'q';qb=tag+'qb'
        self.tg(tag+'in',data,h,clock,cb,'dd_logic',w=self.cfg.get('switch_w',.42))
        self.inv(tag+'i1',h,qb,'dd_logic',w=self.cfg.get('regen_w',.42))
        self.inv(tag+'i2',qb,q,'dd_logic',w=self.cfg.get('regen_w',.42))
        self.tg(tag+'fb',q,h,cb,clock,'dd_logic',w=.42)
        self.fet(tag+'reset','n',h,reset,'0','0',.42)
        self.cap(tag+'load',q,'0',self.cfg.get('digit_load_ff',2));return q,qb
    def mux(self,tag,a,b,en,out,comp='subtract'):
        if comp=='row':
            key=next(k for k in ('positive','negative','zero') if k in tag)
            self.instance_multiplier=self.row_multiplicities[key]
        eb=self.invnode(en,tag+'eb',comp)
        self.tg(tag+'a',a,out,eb,en,'dd_'+comp,w=self.cfg.get('switch_w',.42))
        self.tg(tag+'b',b,out,en,eb,'dd_'+comp,w=self.cfg.get('switch_w',.42))
        self.instance_multiplier=1

def build(cfg):
    w=workload(cfg.get('n',128),cfg.get('pattern','dense'),cfg.get('work_seed',20261026))
    g=golden(w,cfg.get('precision',12),cfg.get('guards',4))
    if 'fixture_r' in cfg:
        rr=F(str(cfg['fixture_r']));zz=F(str(cfg.get('fixture_z',0)));vv=F(0)
        for j,entry in enumerate(g['stages']):
            z=zz if j==0 else F(0);tt=z+2*rr
            digit=1 if tt>F(1,2) else -1 if tt<F(-1,2) else 0
            rr=tt-digit;entry.update(z=float(z),t=float(tt),r=float(rr),digit=digit,
                margin=float(min(abs(tt-F(1,2)),abs(tt+F(1,2)))))
            vv+=w['D']*digit*F(2)**entry['k']
        g['decoded']=float(vv);g['error']=0;g['max_truncation_error']=0
    c=BinaryCircuit(cfg);period=cfg.get('period_ns',400);start=cfg.get('stage_start_ns',.25*period)
    dt=cfg.get('stage_dt_ns',(period-start-.05*period)/len(g['stages']));stop=period+100
    reset=c.clock('reset',.001*period,.025*period)
    awake=c.clock('refawake',.001*period,period-.005*period);ab=c.invnode(awake,'awakeb','reference')
    c.tg('refpower','dd_reference','refrail',awake,ab,'dd_reference',w=cfg.get('ref_w',8.4))
    ladder_r=cfg.get('ref_r',1e3)
    for j in range(16):
        c.poly('ref'+str(j),'0' if j==0 else 'tap'+str(j),'refrail' if j==15 else 'tap'+str(j+1),ladder_r/16,kind='reference')
    cm='tap8';lo='tap6';hi='tap10';threshold='tap9'
    for tag,node in (('refcm',cm),('refhi',hi),('reflo',lo)):c.cap(tag,node,'0',cfg.get('ref_ff',100))
    # Tied-row fixture is optimistic: it measures all its sources but excludes
    # individual row memory, decoder and distributed wires. Never an all-in win.
    rowen=c.clock('roweval',.05*period,period-.01*period)
    total=sum(v['programmed'] for v in w['counts']) or 1
    c.row_multiplicities={key:max(1,round(w['n']*w['bits']*sum(v[key] for v in w['counts'])/total)) for key in ('positive','negative','zero')}
    for sign,name in ((1,'positive'),(-1,'negative'),(0,'zero')):
        for branch in ('p','n'):
            src=cm if sign==0 else hi if (sign>0)==(branch=='p') else lo
            c.mux('row'+name+branch,cm,src,rowen,'row'+name+branch,'row')
    z_nodes={};r0=cfg.get('unit_r',1e4)
    for entry in w['counts']:
        k=entry['k'];ballast=w['D']-entry['programmed'];z_nodes[k]=('zp'+str(k),'zn'+str(k))
        for branch,node in zip(('p','n'),z_nodes[k]):
            for name in ('positive','negative','zero'):
                c.poly(f'k{k}{branch}{name}',node,'row'+name+branch,r0,entry[name],sigma=cfg.get('res_sigma',0))
            if cfg.get('ballast','bank')=='bank':
                c.poly(f'k{k}{branch}ballast',node,cm,r0,ballast,'ballast',cfg.get('res_sigma',0))
            else:c.poly(f'k{k}{branch}ballast',node,cm,r0/max(1,ballast),1,'ballast',cfg.get('res_sigma',0))
            c.cap(node,node,'0',cfg.get('node_ff',20))
    prev=(cm,cm);capff=cfg.get('cap_ff',100)
    if 'fixture_r' in cfg:
        prev=('fixture_rp','fixture_rn')
        for suffix,sg in (('p',1),('n',-1)):
            c.source('fixture_r'+suffix,'fixture_r'+suffix,str(c.d.vdd/2+sg*c.d.vdd/8*cfg['fixture_r']),'external_primitive_input')
            c.source('fixture_z'+suffix,'fixture_z'+suffix,str(c.d.vdd/2+sg*c.d.vdd/8*cfg.get('fixture_z',0)),'external_primitive_input')
    for index,gold in enumerate([] if cfg.get('array_only') else g['stages']):
        tag='s'+str(index);s=start+index*dt;k=gold['k']
        sample=c.clock(tag+'sample',s,s+.24*dt);sb=c.invnode(sample,tag+'sampleb')
        stack=c.clock(tag+'stack',s+.32*dt,period-.015*period)
        evaluate=c.clock(tag+'evaluate',s+.47*dt,s+.78*dt)
        capture=c.clock(tag+'capture',s+.69*dt,s+.76*dt)
        subtract=c.clock(tag+'subtract',s+.84*dt,period-.015*period)
        arrays=z_nodes.get(k,(cm,cm));outs=[];bottoms=[]
        if 'fixture_r' in cfg:arrays=('fixture_zp','fixture_zn') if index==0 else (cm,cm)
        for branch,pr,znode in zip(('p','n'),prev,arrays):
            prefix=tag+branch;nodes=[prefix+'a',prefix+'b',prefix+'t'];bottom=[prefix+'base',prefix+'cbot',prefix+'zbot']
            for j in range(3):
                ff=capff*(cfg.get('residue_ratio',1) if j<2 else 1)
                c.cap(prefix+str(j),nodes[j],bottom[j],ff);c.switch(prefix+'preset'+str(j),nodes[j],cm,reset)
            for j,src in enumerate((pr,pr,znode)):
                c.tg(prefix+'sample'+str(j),src,nodes[j],sample,sb,'dd_logic',w=cfg.get('switch_w',.42))
            for j in (1,2):
                c.tg(prefix+'bottomsample'+str(j),cm,bottom[j],sample,sb,'dd_logic',w=.42)
                c.switch(prefix+'stack'+str(j),nodes[j-1],bottom[j],stack)
                c.switch(prefix+'bottomreset'+str(j),cm,bottom[j],reset)
            c.cap(prefix+'load',nodes[2],'0',cfg.get('res_load_ff',2));outs.append(nodes[2]);bottoms.append(bottom[0])
        qs=[]
        reference_nodes=[]
        refsample=c.clock(tag+'refsample',s+.35*dt,s+.43*dt)
        refshift=c.clock(tag+'refshift',s+.45*dt,s+.79*dt)
        for branch,node,opposite in zip(('p','n'),outs,reversed(outs)):
            ref=threshold
            if cfg.get('comparator_mode','differential')=='differential':
                ref=tag+branch+'reference';bottom=tag+branch+'refbottom'
                c.cap(tag+branch+'reference',ref,bottom,cfg.get('offset_ff',20))
                c.switch(tag+branch+'refs',opposite,ref,refsample,'reference')
                c.switch(tag+branch+'refr',cm,ref,reset,'reference')
                c.mux(tag+branch+'refshift',cm,hi,refshift,bottom,'reference')
            reference_nodes.append(ref)
            op,on=c.comparator(node,ref,evaluate,tag+branch+'cmp')
            q,qb=c.latch(op,capture,tag+branch+'digit',reset);qs.append(q)
        for branch,base in zip(('p','n'),bottoms):
            ptarget=lo if branch=='p' else hi;ntarget=hi if branch=='p' else lo
            c.mux(tag+branch+'pos',cm,ptarget,qs[0],tag+branch+'poslevel')
            c.mux(tag+branch+'neg',tag+branch+'poslevel',ntarget,qs[1],tag+branch+'level')
            c.mux(tag+branch+'sub',cm,tag+branch+'level',subtract,base)
        c.stages.append({'index':index,'k':k,'outp':outs[0],'outn':outs[1],'qp':qs[0],'qn':qs[1],
            'reference_nodes':reference_nodes,'ideal':gold,'t_ns':s+.67*dt,'digit_ns':s+.80*dt,'r_ns':s+.98*dt,'sample_ns':s+.23*dt,'previous':prev})
        prev=tuple(outs)
    c.save.discard('v(0)')
    meta={'cfg':cfg,'kind':'signed_binary_kcl','workload':w,'golden':g,'design':asdict(c.d),'seed':c.d.seed,
        'ports':c.ports,'mos':c.mos,'resistors':c.res,'caps':c.caps,'stages':c.stages,'z_nodes':z_nodes,
        'stop_ns':stop,'period_ns':period,'start_ns':start,'normalization_end_ns':.95*period,
        'voltage_stress':c.d.vdd>1.95,'fixture_limit':'Tied rows; no per-row decoder, wire IR or sequencer.'}
    return c.finish(f'.tran {cfg.get("maxstep_ns",2):.10g}n {stop:.10g}n'),meta

def process(result,names,x,directory):
    meta=result['meta'];cfg=meta['cfg'];ix={n:i for i,n in enumerate(names)};t=x[:,0]
    def v(n):return np.zeros(len(t)) if n=='0' else x[:,ix['v('+str(n).lower()+')']]
    def at(n,ns):return float(np.interp(ns*1e-9,t,v(n)))
    active=t<=meta['period_ns']*1e-9;idle=t>=(meta['period_ns']+25)*1e-9
    energies={};net={};powers=np.zeros(len(t))
    for name,p in meta['ports'].items():
        power=-x[:,ix[name]]*v(p['node']);key=p['component'];powers+=power
        energies[key]=energies.get(key,0)+float(np.trapezoid(np.maximum(power[active],0),t[active]))*1e12
        net[key]=net.get(key,0)+float(np.trapezoid(power[active],t[active]))*1e12
    diss={}
    for r in meta['resistors']:
        power=(v(r['a'])-v(r['b']))*x[:,ix['i(vmeasure'+r['tag'].lower()+')']]
        diss[r['kind']]=diss.get(r['kind'],0)+float(np.trapezoid(power[active],t[active]))*1e12
    wd=meta['workload'];scale=cfg.get('vdd',1.8)/4;cm=cfg.get('vdd',1.8)/2
    stages=[];word=F(0);idealword=F(str(meta['golden']['decoded']));maxdevice=0
    for s in meta['stages']:
        tp=at(s['outp'],s['t_ns']);tn=at(s['outn'],s['t_ns']);tt=(tp-tn)/scale
        qp=at(s['qp'],s['digit_ns']);qn=at(s['qn'],s['digit_ns'])
        dp=int(qp>cfg.get('vdd',1.8)/2);dn=int(qn>cfg.get('vdd',1.8)/2);digit=dp-dn
        residual=(at(s['outp'],s['r_ns'])-at(s['outn'],s['r_ns']))/scale
        word+=wd['D']*digit*F(2)**s['k']
        stages.append({'index':s['index'],'k':s['k'],'t':tt,'ideal_t':s['ideal']['t'],'digit':digit,
            'ideal_digit':s['ideal']['digit'],'r':residual,'ideal_r':s['ideal']['r'],
            'ideal_margin':s['ideal']['margin'],'physical_margin':min(abs(tt-.5),abs(tt+.5)),
            'common_mode_error_v':(tp+tn)/2-cm,'both_high':bool(dp and dn),
            'r_error':residual-s['ideal']['r'],'t_error':tt-s['ideal']['t'],
            'potential_skip':abs(s['ideal']['t'])<=.5,
            'positive_threshold_error_v':at(s['reference_nodes'][0],s['t_ns'])-tn-scale/2,
            'negative_threshold_error_v':at(s['reference_nodes'][1],s['t_ns'])-tp-scale/2})
    for m in meta['mos']:
        vd,vg,vs,vb=(v(m[key]) for key in ('d','g','s','b'))
        maxdevice=max(maxdevice,*[float(np.max(abs(a-b))) for a,b in ((vg,vs),(vg,vd),(vd,vs),(vg,vb))])
    array=[]
    for entry in wd['counts']:
        p,n=meta['z_nodes'].get(str(entry['k']),meta['z_nodes'].get(entry['k']))
        z=(at(p,meta['start_ns']-1)-at(n,meta['start_ns']-1))/scale
        target=entry['c']/wd['D'];array.append({'k':entry['k'],'physical_z':z,'ideal_z':target,'error':z-target})
    area_res=sum(r['area_um2'] for r in meta['resistors']);area_cap=sum(c['area_um2'] for c in meta['caps'])
    wl=sum(m['w_um']*m['l_um']*m.get('multiplicity',1) for m in meta['mos']);allin=sum(energies.values())
    norm=sum(v for k,v in energies.items() if k in ('clock','clock_pin','comparator','logic','subtract'))
    complete=result['simulation_status']=='complete';valuecorrect=word==idealword and not any(s['both_high'] for s in stages)
    metric={'n':wd['n'],'pattern':wd['pattern'],'D':wd['D'],'Q':wd['Q'],'precision':cfg.get('precision',12),
        'guards':cfg.get('guards',4),'stages':len(stages),'corner':cfg.get('corner','tt'),'vdd':cfg.get('vdd',1.8),
        'temp':cfg.get('temp',25),'seed':cfg.get('seed',20261026),'pdk_mm':cfg.get('pdk_mm',False),
        'energy_pj':allin,'net_energy_pj':sum(net.values()),'fj_mac':1000*allin/wd['n'],'normalization_pj':norm,
        'array_dissipation_pj':diss.get('array',0),'ballast_dissipation_pj':diss.get('ballast',0),
        'reference_dissipation_pj':diss.get('reference',0),
        'array_and_ballast_dissipation_pj':diss.get('array',0)+diss.get('ballast',0),
        'idle_nw':float(np.mean(powers[idle]))*1e9,'latency_ns':meta['normalization_end_ns'],'period_ns':meta['period_ns'],
        'decoded':float(word),'ideal_quantized':float(idealword),'dot':wd['dot'],
        'word_error':float(word)-meta['golden']['decoded'],'truncation_error':meta['golden']['error'],
        'truncation_bound':meta['golden']['max_truncation_error'],'word_correct':bool(valuecorrect),
        'digit_failures':sum(s['digit']!=s['ideal_digit'] for s in stages),
        'max_residual_error':max((abs(s['r_error']) for s in stages),default=0),'max_kcl_error':max(abs(a['error']) for a in array),
        'max_common_mode_error_v':max((abs(s['common_mode_error_v']) for s in stages),default=0),
        'min_ideal_margin':min((s['ideal_margin'] for s in stages),default=0),
        'median_ideal_margin':float(np.median([s['ideal_margin'] for s in stages])) if stages else 0,
        'potential_skip_stages':sum(s['potential_skip'] for s in stages),'mos_count':sum(m.get('multiplicity',1) for m in meta['mos']),'simulated_mos_instances':len(meta['mos']),
        'mos_wl_um2':wl,'resistor_body_um2':area_res,'mim_um2':area_cap,'area_proxy_um2':area_res+area_cap+10*wl,
        'physical_product_connections':wd['physical_product_connections'],'product_planes':wd['product_planes'],
        'worst_connections':wd['worst_case_product_connections'],'peak_terminal_v':maxdevice,'simulation_complete':complete,
        'operating_limit_ok':maxdevice<=1.95 and not meta['voltage_stress'],
        'passed':bool(complete and valuecorrect and maxdevice<=1.95 and not meta['voltage_stress']),
        'run':str(Path(directory).relative_to(ROOT))}
    metric.update({'energy_'+k+'_pj':v for k,v in energies.items()})
    metric['fixture_r']=cfg.get('fixture_r','');metric['fixture_z']=cfg.get('fixture_z','')
    result['metrics']=metric;result['stage_observations']=stages;result['array_observations']=array
    dump(Path(directory)/'metrics.json',metric)
    write_csv(Path(directory)/'stages.csv',stages);write_csv(Path(directory)/'kcl.csv',array)
    compress_trace(Path(directory)/'output.raw');return result

def run_configs(configs,category,workers=2):
    backend.BASE=OUT/'backend';backend.BASE.mkdir(parents=True,exist_ok=True);cases=[]
    for i,cfg in enumerate(configs):
        deck,meta=build(cfg);cases.append((deck,meta,OUT/'spice'/category/f'{i:04d}'))
    results=backend.batch_execute(cases,workers,process)
    write_csv(OUT/(category+'.csv'),[r['metrics'] for r in results]);return results

def main():
    p=argparse.ArgumentParser();p.add_argument('--phase',choices=['smoke','nominal','sizing','pvt','mc','pareto','exact'],default='smoke')
    p.add_argument('--workers',type=int,default=2);args=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if args.phase=='exact':
        rows=[]
        for n,pattern in itertools.product(NS,('aligned','aligned_negative','dense','sparse','cancellation','int8','repository_int8')):
            w=workload(n,pattern)
            for guards in (0,4,w['Q']+max(0,2*w['bits']-12)):
                g=golden(w,12,guards)
                rows.append({**{k:v for k,v in w.items() if k!='counts'},'guards':guards,'stages':len(g['stages']),
                    'truncation_error':g['error'],'truncation_bound':g['max_truncation_error'],
                    'min_margin':min(r['margin'] for r in g['stages']),
                    'skip_stages':sum(abs(r['t'])<=.5 for r in g['stages'])})
        write_csv(OUT/'exact_workloads.csv',rows);dump(OUT/'exact_workloads.json',rows);print('Exact cases',len(rows));return
    if args.phase=='smoke':configs=[{'n':64,'pattern':k} for k in ('aligned','dense','aligned_negative','cancellation')]
    elif args.phase=='nominal':configs=[{'n':n,'pattern':k} for n,k in itertools.product(NS,('aligned','aligned_negative','dense','sparse','cancellation'))]
    elif args.phase=='sizing':
        base={'n':128,'pattern':'aligned'};configs=[base]
        sweeps={'unit_r':[1e3,1e5,1e6],'ref_r':[1,10,100,1e4,1e5],'node_ff':[5,100],'cap_ff':[20,500,2000],
            'switch_w':[.84,1.68],'cmp_w':[.84,1.68],'regen_w':[.84,1.68],'driver_w':[.42,1.68],
            'residue_ratio':[.5,2,4],'res_load_ff':[.1,10],'offset_ff':[5,100,500],'ref_w':[4.2,16.8],'ballast':['single'],'period_ns':[100,1600]}
        for key,vals in sweeps.items():configs.extend({**base,key:v} for v in vals)
        for cap,r,ref in itertools.product((20,500),(1e3,1e5),(100,1e4)):
            configs.append({**base,'cap_ff':cap,'unit_r':r,'ref_r':ref})
    elif args.phase=='pvt':configs=[{'n':128,'pattern':'aligned','corner':co,'vdd':v,'temp':t}
        for co,v,t in itertools.product(('tt','ss','ff','sf','fs'),(1.62,1.8,1.98),(-40,25,85,125))]
    elif args.phase=='mc':configs=[{'n':128,'pattern':'aligned','pdk_mm':True,'seed':260100+j,'cap_sigma':0,'res_sigma':0} for j in range(16)]
    else:configs=[{'n':n,'pattern':'aligned','period_ns':period} for n,period in itertools.product(NS,(25,100,400,1600))]
    results=run_configs(configs,args.phase,args.workers)
    print('Completed',args.phase,len(results),'passes',sum(r['metrics']['passed'] for r in results),flush=True)
    for r in results[:4]:print(json.dumps(r['metrics']),flush=True)
if __name__=='__main__':main()
