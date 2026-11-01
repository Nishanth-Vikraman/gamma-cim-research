"""Fully clocked, zero-intentional-bias q4 quotient/one-shot packet circuits.

PWL sources specify request phase timing and external analog test inputs only.
All quotient-dependent reference selection, latching and packet enables are MOS.
MIM capacitors and ordinary SKY130 devices are used without PDK modification.
"""
from dataclasses import asdict
from math import sqrt,ceil
import json,hashlib
import numpy as np
from gamma_sky130 import Circuit,Design,ROOT,pwl,fast_raw,compress_trace
from gamma_sar_exact import OUT,quotient

def nand(c,tag,inputs,out,dd,w=.42):
    for j,g in enumerate(inputs):
        c.fet(tag+'p'+str(j),'p',out,g,dd,dd,2*w)
        c.fet(tag+'n'+str(j),'n',out,g,'0' if j==len(inputs)-1 else tag+'s'+str(j),
              '0',w*len(inputs)) if j==0 else c.fet(tag+'n'+str(j),'n',tag+'s'+str(j-1),g,
              '0' if j==len(inputs)-1 else tag+'s'+str(j),'0',w*len(inputs))

class Dynamic(Circuit):
    def __init__(self,cfg):
        self.cfg=cfg
        d=Design(corner=cfg.get('corner','tt'),vdd=cfg.get('vdd',1.8),temp_c=cfg.get('temp',25),
          pdk_mismatch=cfg.get('pdk_mm',False),seed=cfg.get('seed',41001),
          capacitor_sigma=cfg.get('cap_sigma',0),resistor_sigma=cfg.get('ref_sigma',0),
          comparator_w=cfg.get('cmp_w',.42),comparator_l=cfg.get('cmp_l',.15),cap_kind='mim',cap_ff=cfg.get('store_ff',1000))
        super().__init__(d)
        self.lines.append('.option reltol=1e-4 abstol=1e-13 vntol=1e-7')
        for k in ('comparator','logic','reference','transfer','sample'):self.supply(k)
        self.tags=0;self.phase_points={};self.checks=[];self.ref_nodes=[];self.sense=[]
        self.probes=set();self.sleep_extra=[]
        self.source('zero','zero','0','constant_pin')
        self.source('one','one',str(d.vdd),'constant_pin')
    def not_(self,a,tag=None,component='logic'):
        tag=tag or self.tag('not');self.inv(tag,a,tag+'o','dd_'+component);return tag+'o'
    def tag(self,prefix):self.tags+=1;return prefix+str(self.tags)
    def and_(self,inputs,tag=None,component='logic'):
        if len(inputs)==1:return inputs[0]
        tag=tag or self.tag('and');nand(self,tag,inputs,tag+'b','dd_'+component)
        return self.not_(tag+'b',tag+'i',component)
    def or_(self,inputs,tag=None):
        return self.not_(self.and_([self.not_(v) for v in inputs]),tag)
    def mux(self,tag,a,b,sel,out,component='logic'):
        sb=self.not_(sel,tag+'sb',component);self.tg(tag+'a',a,out,sb,sel,'dd_'+component,w=.42)
        self.tg(tag+'b',b,out,sel,sb,'dd_'+component,w=.42)
    def latch(self,tag,data,clk,reset):
        cb=self.not_(clk,tag+'cb');h=tag+'h';q=tag+'q';qb=tag+'qb'
        self.tg(tag+'in',data,h,clk,cb,'dd_logic',w=.42)
        self.inv(tag+'i1',h,qb,'dd_logic');self.inv(tag+'i2',qb,q,'dd_logic')
        self.tg(tag+'fb',q,h,cb,clk,'dd_logic',w=.42)
        self.fet(tag+'rst','n',h,reset,'0','0',.42)
        self.probes.update((h,q,qb));return q,qb
    def cap_between(self,tag,a,b,ff):
        ff*=1+self.rng.gauss(0,self.d.capacitor_sigma)
        if ff<=0:raise ValueError('nonpositive MIM realization')
        side=(-.76+sqrt(.76**2+8*ff))/4+.025
        self.lines.append(f'Xcap{tag} {a} {b} sky130_fd_pr__cap_mim_m3_1 w={side:.12g} l={side:.12g}')
    def capacitor(self,node,cap_ff=None,kind=None):
        # Mismatch perturbs a fixed layout: it must not change instance count
        # when a nominal 1pF capacitor happens to vary just above 1pF.
        nominal=cap_ff or self.d.cap_ff;count=max(1,ceil(nominal/1000))
        actual=nominal*(1+self.rng.gauss(0,self.d.capacitor_sigma))
        if actual<=0:raise ValueError('nonpositive capacitor realization')
        side=(-.76+sqrt(.76**2+8*actual/count))/4+.025
        for j in range(count):self.lines.append(f'Xc{node}_{j} {node} 0 sky130_fd_pr__cap_mim_m3_1 w={side:.12g} l={side:.12g}')
        self.save.add(f'v({node})')
    def clock(self,name,intervals,component='logic'):
        points=[(0,0)]
        for a,b in intervals:points += [(a,0),(a+.1,self.d.vdd),(b-.1,self.d.vdd),(b,0)]
        self.source(name+'pin',name+'pin',pwl(points),'clock_pin')
        b=self.not_(name+'pin',name+'drv1',component);q=name+'drv2o'
        drive=4.2 if component in ('sample','transfer') else (.42 if component!='reference' else self.cfg.get('ref_driver_w',.42))
        self.inv(name+'drv2',b,q,'dd_'+component,w=drive)
        self.phase_points[name]=intervals;return q
    def dynamic_ref(self,tag,threshold,prepare,charge):
        """Two-cap divider with both bottom plates initially zero; top reset
        opens before the DD plate rises. Reference error includes comparator
        loading, reset injection and top parasitics. No voltage-source tap.
        """
        ff=self.cfg.get('ref_ff',1000);ratio=1/3+threshold/144
        ratio+=self.cfg.get('ref_offset_units',0)/144
        ratio*=1+self.rng.gauss(0,self.cfg.get('ratio_sigma',0))
        node=tag+'v';reset=self.not_(prepare,tag+'resetctl','reference')
        if self.cfg.get('ref_kind','capacitive')=='capacitive':
            self.cap_between(tag+'a',node,charge,ff*ratio)
            self.cap_between(tag+'b',node,'0',ff*(1-ratio))
            self.tg(tag+'reset','0',node,reset,prepare,'dd_reference',w=.42)
        else:
            # Gated ladder precharges a hold capacitor, then is disconnected.
            # A permanently buffered/static reference is an explicit control.
            self.capacitor(node,ff,'mim');en=prepare
            if self.cfg['ref_kind']=='shared_static':en='one'
            eb=self.not_(en,tag+'eb','reference')
            self.tg(tag+'supply','dd_reference',tag+'rail',en,eb,'dd_reference',w=self.cfg.get('ref_switch_w',.42))
            r=self.cfg.get('ref_r',2e4)
            self.resistor(tag+'rt',tag+'rail',node,r*(1-ratio));self.resistor(tag+'rb',node,'0',r*ratio)
            # Bottom leg must also be gated for the dynamic precharge case.
            if self.cfg['ref_kind']!='shared_static':
                self.lines[-1]=f'R{tag}rb {node} {tag}bottom {r*ratio:.12g}'
                self.tg(tag+'bottom',tag+'bottom','0',en,eb,'dd_reference',w=self.cfg.get('ref_switch_w',.42))
        self.ref_nodes.append({'node':node,'threshold_units':threshold,'prepare':prepare,'charge':charge})
        self.probes.add(node);return node
    def packet(self,tag,node,m,sgn,enable,pre,fire,sleep,kind):
        """One packet phase; all selected magnitude elements act in parallel."""
        cff=self.d.cap_ff;scale=m
        sink=(sgn>0)==(node=='center');mult=4 if node=='center' else 1
        if kind=='gated_current' and self.cfg.get('shared_current_bias',False):
            if not hasattr(self,'shared_bias'):
                if self.cfg.get('current_early_bias',False):
                    intervals=[(r[0]+9,w[1]) for r,w in zip(self.phase_points['reset'],self.phase_points['packetawake'])]
                    active=self.clock('sharedbiasearly',intervals,'transfer');self.sleep_extra.append(active);self.probes.add(active)
                else:active=self.not_(sleep,'sharedbiasawake','transfer')
                ab=self.not_(active,'sharedbiasab','transfer')
                self.tg('sharednbpower','dd_transfer','sharednr',active,ab,'dd_transfer',w=4.2)
                self.resistor('sharednb','sharednr','sharednb',self.cfg.get('packet_r',1.85e6),False)
                self.fet('sharednb','n','sharednb','sharednb','0','0',.42,1)
                self.tg('sharedpbpower','sharedpr','0',active,ab,'dd_transfer',w=4.2)
                self.resistor('sharedpb','sharedpb','sharedpr',self.cfg.get('packet_r',1.85e6),False)
                self.fet('sharedpb','p','sharedpb','sharedpb','dd_transfer','dd_transfer',.42,1)
                self.shared_bias=True
            sw=self.and_([fire,enable],tag+'swctl');sb=self.not_(sw,tag+'swb')
            target=tag+'port';self.lines.append(f'V{tag}sense {node} {target} 0');self.sense.append((tag,node,sgn,m))
            if sink:
                self.fet(tag+'ib','n',tag+'s','sharednb','0','0',.42*mult*scale,1)
                self.fet(tag+'sw','n',target,sw,tag+'s','0',.84*mult*scale)
            else:
                self.fet(tag+'ib','p',tag+'s','sharedpb','dd_transfer','dd_transfer',.42*mult*scale,1)
                self.fet(tag+'sw','p',target,sb,tag+'s','dd_transfer',1.68*mult*scale)
            return
        if kind=='gated_current':
            active=self.and_([pre,enable],tag+'biasenable');ab=self.not_(active,tag+'ab')
            if sink:
                self.tg(tag+'biaspower','dd_transfer',tag+'vr',active,ab,'dd_transfer',w=.42)
                self.resistor(tag+'bias',tag+'vr',tag+'b',self.cfg.get('packet_r',3e6),False)
                self.fet(tag+'diode','n',tag+'b',tag+'b','0','0',.42,1)
                self.fet(tag+'ib','n',tag+'s',tag+'b','0','0',.42*mult*scale,1)
                sw=self.and_([fire,enable],tag+'swctl')
                target=tag+'port';self.lines.append(f'V{tag}sense {node} {target} 0');self.sense.append((tag,node,sgn,m))
                self.fet(tag+'sw','n',target,sw,tag+'s','0',.84*mult*scale)
            else:
                self.tg(tag+'biaspower',tag+'low','0',active,ab,'dd_transfer',w=.42)
                self.resistor(tag+'bias',tag+'b',tag+'low',self.cfg.get('packet_r',3e6),False)
                self.fet(tag+'diode','p',tag+'b',tag+'b','dd_transfer','dd_transfer',.42,1)
                self.fet(tag+'ib','p',tag+'s',tag+'b','dd_transfer','dd_transfer',.42*mult*scale,1)
                sw=self.not_(self.and_([fire,enable],tag+'swctl'),tag+'swb')
                target=tag+'port';self.lines.append(f'V{tag}sense {node} {target} 0');self.sense.append((tag,node,sgn,m))
                self.fet(tag+'sw','p',target,sw,tag+'s','dd_transfer',1.68*mult*scale)
            return
        # Nominal charge-sharing calibration at residual 1.5, neighbors zero.
        if node=='center':ratio=4*scale/(49.5 if sgn>0 else 97.5)
        else:ratio=scale/((96-scale) if sgn>0 else (48-scale))
        if kind=='binary_sc':
            # Fixed binary-weighted unit capacitors; no nonlinear m calibration.
            ratio=scale*(4/(49.5 if sgn>0 else 97.5) if node=='center' else 1/(95 if sgn>0 else 47))
        if kind=='segmented_sc':
            ratio=scale*(4/(49.5 if sgn>0 else 97.5) if node=='center' else 1/(95 if sgn>0 else 47))
        self.capacitor(tag+'cap',cff*ratio,'mim')
        pc=self.and_([pre,enable],tag+'pc');pb=self.not_(pc,tag+'pcb')
        source='0' if sink else 'dd_transfer'
        self.tg(tag+'precharge',source,tag+'cap',pc,pb,'dd_transfer',w=.42)
        sw=self.and_([fire,enable],tag+'swctl');sb=self.not_(sw,tag+'swb')
        target=tag+'port';self.lines.append(f'V{tag}sense {node} {target} 0');self.sense.append((tag,node,sgn,m))
        self.tg(tag+'join',tag+'cap',target,sw,sb,'dd_transfer',w=.42)
        self.fet(tag+'reset','n',tag+'cap',sleep,'0','0',.42)

def build(cfg):
    if cfg.get('scheme')=='bank':return bank_build(cfg)
    c=Dynamic(cfg);limit=cfg.get('limit',15);bits=(limit//4).bit_length();M=limit//4
    states=cfg.get('states',list(range(-limit,limit+1)));period=cfg.get('period_ns',400)
    unit=c.d.vdd/144;cm=c.d.vdd/3;starts=[100+j*period for j in range(len(states))]
    reset=c.clock('reset',[(s,s+8) for s in starts])
    sample=c.clock('sample',[(s+9,s+20) for s in starts],component='sample');sampleb=c.not_(sample,'sampleb','sample')
    # Same nominal physical storage and same sample TG on all three positions.
    for node in ('center','left','right'):
        c.capacitor(node,c.d.cap_ff,'mim')
        vals=[(0,cm)]
        for s,value in zip(starts,states):
            coeff=value if node=='center' else cfg.get('neighbor_'+node,0)
            vals.extend([(s,cm+coeff*unit),(s+period-.1,cm+coeff*unit)])
        c.source(node+'stim',node+'stim',pwl(vals),'test_state_write')
        c.tg(node+'sample',node+'stim',node,sample,sampleb,'dd_sample',w=8.4)
        # Half-sized dummy switches compensate input-fixture charge injection.
        c.fet(node+'dummyn','n',node,sampleb,node,'0',4.2)
        c.fet(node+'dummyp','p',node,sample,node,'dd_sample',8.4)
    cmpclk=c.clock('evaluate',[(s+35+40*j,s+51+40*j) for s in starts for j in range(bits+1)])
    rawop,rawon=c.comparator('center','reference',cmpclk)
    op=c.not_(c.not_(rawop,'outp1'),'outp2');on=c.not_(c.not_(rawon,'outn1'),'outn2')
    sampleclks=[c.clock('capture'+str(j),[(s+48+40*j,s+50+40*j) for s in starts]) for j in range(bits+1)]
    sign,signb=c.latch('sign',op,sampleclks[0],reset)
    # A MOS mux chooses the correct comparator polarity for magnitude decisions.
    c.mux('magnitudepolarity',on,op,sign,'magdecision')
    qb={};q={}
    for j in reversed(range(bits)):
        stage=bits-j;q[j],qb[j]=c.latch('bit'+str(j),'magdecision',sampleclks[stage],reset)
    refcharges={};refprepares={}
    for stage in range(bits+1):
        prepare=c.clock('refprepare'+str(stage),[(s+23+40*stage,s+57+40*stage) for s in starts],'reference')
        charge=c.clock('refcharge'+str(stage),[(s+27+40*stage,s+55+40*stage) for s in starts],'reference')
        refcharges[stage]=charge;refprepares[stage]=prepare
        if cfg.get('scheme')=='binary_residual' and stage>0:continue
        possibilities=[(-.5,[])] if stage==0 else []
        if stage:
            j=bits-stage
            for prefix in range(2**(stage-1)):
                trial=(prefix<<(j+1))+(1<<j);threshold=4*trial-.5
                conditions=[q[k] if prefix&(1<<(k-j-1)) else qb[k] for k in range(j+1,bits)]
                possibilities.extend([(threshold,[sign,*conditions]),(-threshold,[signb,*conditions])])
        for index,(threshold,conds) in enumerate(possibilities):
            tag=f'r{stage}x{index}'
            pr=c.and_([prepare,*conds],tag+'pr',component='reference')
            ch=c.and_([charge,*conds],tag+'ch',component='reference')
            ref=c.dynamic_ref(tag,threshold,pr,ch);prb=c.not_(pr,tag+'prb','reference')
            c.tg(tag+'mux',ref,'reference',pr,prb,'dd_reference',w=.42)
    if cfg.get('scheme')=='binary_residual':
        active=c.clock('binaryactive',[(s+63,s+40*bits+57) for s in starts],'reference')
        base=c.clock('binarybase',[(s+67,s+40*bits+55) for s in starts],'reference')
        resetref=c.not_(active,'binaryreset','reference');ff=cfg.get('ref_ff',1000)
        for sg,node,baseunits in ((1,'binarypositive',47.5),(-1,'binarynegative',48.5)):
            c.cap_between(node+'base',node,base,ff*baseunits/144)
            rest=ff*(1-(baseunits+4*M)/144);c.cap_between(node+'ground',node,'0',rest)
            for j in range(bits):
                stage=bits-j;later=[refcharges[k] for k in range(stage+1,bits+1)]
                past=c.and_([q[j],c.or_(later)],f'binarypast{sg+1}b{j}') if later else 'zero'
                drive=c.or_([refcharges[stage],past],f'binarydrive{sg+1}b{j}')
                if sg<0:drive=c.not_(drive,f'binaryinvert{j}','reference')
                c.cap_between(node+'b'+str(j),node,drive,ff*4*(1<<j)/144)
            c.tg(node+'reset','0',node,resetref,active,'dd_reference',w=.42)
            en=c.and_([active,sign if sg>0 else signb],node+'enable',component='reference')
            eb=c.not_(en,node+'eb','reference');c.tg(node+'mux',node,'reference',en,eb,'dd_reference',w=.42)
            c.probes.add(node)
    transfer_start=40*(bits+1)+25
    pre=c.clock('packetpre',[(s+transfer_start,s+transfer_start+15) for s in starts],'transfer')
    fire=c.clock('packetfire',[(s+transfer_start+18,s+transfer_start+38) for s in starts],'transfer')
    # Packet reset is enabled at sleep and opens before precharge.
    awake=c.clock('packetawake',[(s+transfer_start-2,s+transfer_start+42) for s in starts],'transfer')
    sleep=c.not_(awake,'packetsleep','transfer')
    kind=cfg.get('transfer_kind','redistribution')
    if kind in ('binary_sc','gated_current'):
        enables=[(1<<j,q[j]) for j in range(bits)]
    else:
        dec=[]
        for m in range(1,M+1):dec.append((m,c.and_([q[j] if m&(1<<j) else qb[j] for j in range(bits)],'decode'+str(m))))
        if kind=='segmented_sc':enables=[(1,c.or_([en for m,en in dec if m>=j],f'thermo{j}')) for j in range(1,M+1)]
        else:enables=dec
    for j,(m,en) in enumerate(enables):
        for sgn,sg in ((1,sign),(-1,signb)):
            active=c.and_([en,sg],f'packet_enable_{j}_{sgn+1}')
            for node in ('center','left','right'):
                c.packet(f'p{j}s{sgn+1}{node}',node,m,sgn,active,pre,fire,sleep,kind)
    for s,value in zip(starts,states):
        c.checks.append({'start_ns':s,'sample_ns':s+21,'decision_ns':s+transfer_start-5,
          'before_transfer_ns':s+transfer_start+16,'after_transfer_ns':s+transfer_start+80,
          'end_ns':s+period-1,'coefficient':value,'golden_m':quotient(value),
          'neighbor_initial':{n:cfg.get('neighbor_'+n,0) for n in ('left','right')},
          'cmp_times_ns':[s+49+40*j for j in range(bits+1)],'eval_times_ns':[s+35+40*j for j in range(bits+1)]})
    stop=starts[-1]+period+300
    c.probes.update(['center','left','right','reference',rawop,rawon,op,on,cmpclk,sample,pre,fire,sleep,sign,signb,*q.values(),*qb.values()])
    c.save={f'v({n})' for n in c.probes}|set(c.ports)|{f'v({v["node"]})' for v in c.ports.values()}|{f'i(v{tag}sense)' for tag,*_ in c.sense}
    meta={'kind':'dynamic_sar','cfg':cfg,'seed':c.d.seed,'design':asdict(c.d),'ports':c.ports,'mos':c.mos,
      'checks':c.checks,'refs':c.ref_nodes,'sense':c.sense,'bits':bits,'q_nodes':q,'sign_node':sign,
      'sleep_controls':[cmpclk,sample,pre,fire,*c.sleep_extra],'sleep_high':sleep,
      'stop_ns':stop,'unit_v':unit,'nominal_unit_charge_fc':c.d.cap_ff*unit,
      'normalizer_external_controls':'Input-independent phase clock pins; MOS phase buffers, latches, reference mux, magnitude decode and packet enables are physical. Global scheduler/clock generator not transistorized.'}
    names=[m['name'].lower() for m in c.mos]
    if len(names)!=len(set(names)):raise ValueError('Duplicate physical MOS instance name')
    return c.finish(f'.tran .1n {stop}n 0 {cfg.get("maxstep_ns",.2)}n'),meta

def integral(t,y,a,b):
    m=(t>a)&(t<b);return float(np.trapezoid(np.r_[np.interp(a,t,y),y[m],np.interp(b,t,y)],np.r_[a,t[m],b]))

def bank_build(cfg):
    """Parallel coarse thresholds with physical thermometer-to-binary logic."""
    c=Dynamic(cfg);limit=cfg.get('limit',15);M=limit//4;bits=M.bit_length()
    states=cfg.get('states',list(range(-limit,limit+1)));period=cfg.get('period_ns',400)
    starts=[100+j*period for j in range(len(states))];dd=c.d.vdd;cm=dd/3;unit=dd/144
    reset=c.clock('reset',[(s,s+8) for s in starts]);sample=c.clock('sample',[(s+9,s+20) for s in starts],'sample')
    sb=c.not_(sample,'sampleb','sample')
    for node in ('center','left','right'):
        c.capacitor(node,c.d.cap_ff,'mim');vals=[(0,cm)]
        for s,vv in zip(starts,states):
            coeff=vv if node=='center' else cfg.get('neighbor_'+node,0)
            vals.extend([(s,cm+coeff*unit),(s+period-.1,cm+coeff*unit)])
        c.source(node+'stim',node+'stim',pwl(vals),'test_state_write');c.tg(node+'sample',node+'stim',node,sample,sb,'dd_sample',w=8.4)
        c.fet(node+'dummyn','n',node,sb,node,'0',4.2);c.fet(node+'dummyp','p',node,sample,node,'dd_sample',8.4)
    clock=c.clock('evaluate',[(s+35,s+51) for s in starts]);capclk=c.clock('capture',[(s+48,s+50) for s in starts])
    prepare=c.clock('refprepare',[(s+23,s+57) for s in starts],'reference');charge=c.clock('refcharge',[(s+27,s+55) for s in starts],'reference')
    flags={};negative={};sleep_nodes=[]
    for j in range(1,M+1):
        ps=[]
        for sg in (1,-1):
            tag=f'bank{j}s{sg+1}';ref=c.dynamic_ref(tag,sg*(4*j-.5),prepare,charge)
            cmpname='cmp' if j==1 and sg==1 else tag+'cmp'
            if sg==1:op,on=c.comparator('center',ref,clock,cmpname)
            else:op,on=c.comparator(ref,'center',clock,cmpname)
            # Balanced output loading prevents data-dependent dynamic offsets.
            pop=c.not_(c.not_(op,tag+'p1'),tag+'p2');nop=c.not_(c.not_(on,tag+'n1'),tag+'n2')
            q,qb=c.latch(tag,pop,capclk,reset);ps.append(q)
            c.probes.update((op,on,pop,nop,q,qb));sleep_nodes.extend((op,on))
            if sg<0:negative[j]=q
        flags[j]=c.or_(ps,f'flag{j}')
    sign=c.not_(negative[1],'sign');signb=c.not_(sign,'signb');q={};qb={}
    for bit in range(bits):
        terms=[]
        for m in range(1,M+1):
            if m&(1<<bit):terms.append(c.and_([flags[m],c.not_(flags[m+1])]) if m<M else flags[m])
        q[bit]=c.or_(terms,f'encodebit{bit}') if len(terms)>1 else terms[0];qb[bit]=c.not_(q[bit],f'bit{bit}b')
    pre=c.clock('packetpre',[(s+65,s+80) for s in starts],'transfer');fire=c.clock('packetfire',[(s+83,s+103) for s in starts],'transfer')
    awake=c.clock('packetawake',[(s+63,s+107) for s in starts],'transfer');sleep=c.not_(awake,'packetsleep','transfer')
    kind=cfg.get('transfer_kind','binary_sc')
    if kind in ('binary_sc','gated_current'):enables=[(1<<j,q[j]) for j in range(bits)]
    else:
        dec=[(m,c.and_([q[j] if m&(1<<j) else qb[j] for j in range(bits)],f'decode{m}')) for m in range(1,M+1)]
        if kind=='segmented_sc':enables=[(1,c.or_([v for m,v in dec if m>=j],f'thermo{j}')) for j in range(1,M+1)]
        else:enables=dec
    for j,(m,en) in enumerate(enables):
        for sg,sgnode in ((1,sign),(-1,signb)):
            active=c.and_([en,sgnode],f'packet_enable_{j}_{sg+1}')
            for node in ('center','left','right'):c.packet(f'p{j}s{sg+1}{node}',node,m,sg,active,pre,fire,sleep,kind)
    checks=[{'start_ns':s,'sample_ns':s+21,'decision_ns':s+60,'before_transfer_ns':s+81,'after_transfer_ns':s+145,
             'end_ns':s+period-1,'coefficient':v,'golden_m':quotient(v),'neighbor_initial':{n:cfg.get('neighbor_'+n,0) for n in ('left','right')},
             'cmp_times_ns':[s+49],'eval_times_ns':[s+35]} for s,v in zip(starts,states)]
    c.probes.update(['center','left','right',clock,sample,pre,fire,sleep,sign,*q.values()])
    c.save={f'v({n})' for n in c.probes}|set(c.ports)|{f'v({v["node"]})' for v in c.ports.values()}|{f'i(v{tag}sense)' for tag,*_ in c.sense}
    # Reference voltage export uses the +4 threshold as a representative tap.
    c.probes.add('bank1s2v');c.save.add('v(bank1s2v)')
    stop=starts[-1]+period+300
    meta={'kind':'dynamic_bank','cfg':cfg,'seed':c.d.seed,'design':asdict(c.d),'ports':c.ports,'mos':c.mos,
      'checks':checks,'refs':c.ref_nodes,'sense':c.sense,'bits':bits,'q_nodes':q,'sign_node':sign,
      'sleep_controls':[clock,sample,pre,fire,*c.sleep_extra],'sleep_high':sleep,'sleep_comparator_nodes':sleep_nodes,
      'reference_probe':'bank1s2v','stop_ns':stop,'unit_v':unit,'nominal_unit_charge_fc':c.d.cap_ff*unit,
      'normalizer_external_controls':'Input-independent phase clock pins, physical parallel threshold/latch/thermometer encoder and packet enables; global request scheduler unbuilt.'}
    names=[m['name'].lower() for m in c.mos]
    if len(names)!=len(set(names)):raise ValueError('Duplicate physical MOS instance name')
    return c.finish(f'.tran .1n {stop}n 0 {cfg.get("maxstep_ns",.5)}n'),meta

def analyze(result,names,x):
    meta=result['meta'];cfg=meta['cfg'];ix={n:j for j,n in enumerate(names)};t=x[:,0]
    if t[-1]<(meta['stop_ns']-.01)*1e-9:raise ValueError('incomplete transient')
    v=lambda n:x[:,ix['v('+n.lower()+')']]
    powers={};pin={};stimulus={}
    for name,p in meta['ports'].items():
        power=-v(p['node'])*x[:,ix[name]];component=p['component']
        if component=='clock_pin':pin[name]=power
        elif component=='test_state_write':stimulus[name]=power
        elif component=='constant_pin':pin[name]=power
        else:powers[component]=powers.get(component,0)+power
    def energy(a,b):
        out={k:integral(t,p,a*1e-9,b*1e-9)*1e15 for k,p in powers.items()}
        out['external_clock_pin_gross']=sum(integral(t,np.maximum(p,0),a*1e-9,b*1e-9)*1e15 for p in pin.values())
        return out
    idle_a=meta['checks'][-1]['end_ns']+20;idle_b=meta['stop_ns']-10
    idle=energy(idle_a,idle_b);idle_nw=sum(idle.values())/(idle_b-idle_a)*1000
    obs=[];unit=meta['unit_v'];qunit=meta['nominal_unit_charge_fc']*1e-15
    for p in meta['checks']:
        a=p['start_ns'];b=p['end_ns'];td=p['decision_ns']*1e-9
        resolved=all((np.interp(td,t,v(node))<.2*meta['design']['vdd'] or np.interp(td,t,v(node))>.8*meta['design']['vdd']) for node in meta['q_nodes'].values())
        magnitude=int(sum((1<<int(j))*(np.interp(td,t,v(node))>meta['design']['vdd']/2) for j,node in meta['q_nodes'].items()))
        sign=1 if np.interp(td,t,v(meta['sign_node']))>meta['design']['vdd']/2 else -1
        chosen=sign*magnitude;before={n:float(np.interp(p['before_transfer_ns']*1e-9,t,v(n))) for n in ('center','left','right')}
        after={n:float(np.interp(p['after_transfer_ns']*1e-9,t,v(n))) for n in before}
        charge={n:-sum(integral(t,x[:,ix[f'i(v{tag}sense)']],p['before_transfer_ns']*1e-9,p['after_transfer_ns']*1e-9) for tag,node,*_ in meta['sense'] if node==n)/qunit for n in before}
        target={'center':-4*chosen,'left':chosen,'right':chosen}
        charge_pass=all(abs(charge[n]-target[n])<=max(.1,.05*abs(target[n])) for n in charge)
        resid=(after['center']-meta['design']['vdd']/3)/unit
        decoded=int(np.floor(resid+.5));expected=p['coefficient']-4*p['golden_m']
        neighbor_classes={n:int(np.floor((after[n]-meta['design']['vdd']/3)/unit+.5)) for n in ('left','right')}
        neighbor_targets={n:p.get('neighbor_initial',{}).get(n,0)+p['golden_m'] for n in ('left','right')}
        energy_end=a+cfg.get('period_ns',400)
        supply=energy(a,energy_end)
        refs=[]
        reference_checks=[]
        for ti in p['cmp_times_ns']:
            refs.append(float(np.interp((ti-12)*1e-9,t,v(meta.get('reference_probe','reference')))))
        if meta['kind']=='dynamic_sar':
            for stage,ev in enumerate(p['eval_times_ns']):
                if stage==0:threshold=-.5
                else:
                    j=meta['bits']-stage
                    prior=sum((1<<int(k))*(np.interp((ev-.5)*1e-9,t,v(node))>meta['design']['vdd']/2)
                              for k,node in meta['q_nodes'].items() if int(k)>j)
                    threshold=(4*(prior+(1<<j))-.5)*(1 if np.interp((ev-.5)*1e-9,t,v(meta['sign_node']))>meta['design']['vdd']/2 else -1)
                measured=float(np.interp((ev-.2)*1e-9,t,v('reference')))
                reference_checks.append({'stage':stage,'threshold_units':float(threshold),'before_evaluate_v':measured,
                  'error_units':float((measured-meta['design']['vdd']/3)/unit-threshold),
                  'hold_plus_comparator_shift_units':float((np.interp((ev+14)*1e-9,t,v('reference'))-measured)/unit)})
        sleep_ok=all(abs(np.interp((b-10)*1e-9,t,v(node)))<.1*meta['design']['vdd'] for node in meta.get('sleep_comparator_nodes',('cmpop','cmpon')))
        sleep_ok=sleep_ok and all(abs(np.interp((b-10)*1e-9,t,v(node)))<.1*meta['design']['vdd'] for node in meta.get('sleep_controls',[]))
        if meta.get('sleep_high'):sleep_ok=sleep_ok and np.interp((b-10)*1e-9,t,v(meta['sleep_high']))>.9*meta['design']['vdd']
        reset_intervals=[(a,a+9)]+[(ev+16,ev+27) for ev in p['eval_times_ns']]+[(p['before_transfer_ns']+22,p['after_transfer_ns'])]
        reset_work={k:0. for k in supply}
        for lo,hi in reset_intervals:
            for k,val in energy(lo,min(hi,b)).items():reset_work[k]+=val
        nonreset={k:supply[k]-reset_work[k] for k in supply}
        quiet_start=max(p['after_transfer_ns']+10,b-95);quiet_end=b-5
        reference_drift=max([abs(float(np.interp(quiet_end*1e-9,t,v(ref['node']))-np.interp(quiet_start*1e-9,t,v(ref['node']))))/unit for ref in meta['refs']],default=0.)
        end_values={n:float(np.interp((b-5)*1e-9,t,v(n))) for n in after}
        end_classes={n:int(np.floor((val-meta['design']['vdd']/3)/unit+.5)) for n,val in end_values.items()}
        after_classes={'center':decoded,**neighbor_classes}
        retained=all(end_classes[n]==after_classes[n] for n in after_classes)
        obs.append({**p,'chosen_m':int(chosen),'quotient_correct':bool(resolved and chosen==p['golden_m']),
          'energy_end_ns':energy_end,
          'quotient_error':int(chosen-p['golden_m']),
          'bits_resolved':bool(resolved),'charge_q':charge,'target_charge_q':target,'transfer_charge_correct':charge_pass,
          'center_residual_digits':resid,'decoded_class':decoded,'golden_class':expected,
          'neighbor_classes':neighbor_classes,'golden_neighbor_class':p['golden_m'],
          'golden_neighbor_classes':neighbor_targets,
          'physical_class_correct':decoded==expected and -3<=decoded<=3 and all(vv==neighbor_targets[n] for n,vv in neighbor_classes.items()),
          'sleep_control_correct':bool(sleep_ok),'retention_class_preserved':retained,
          'sleep_correct':bool(sleep_ok and retained),
          'post_sleep_classes':end_classes,'post_sleep_drift_digits':{n:(end_values[n]-after[n])/unit for n in after},
          'state_sample_error_digits':(float(np.interp(p['sample_ns']*1e-9,t,v('center')))-meta['design']['vdd']/3)/unit-p['coefficient'],
          'voltage_transfer_digits':{n:(after[n]-before[n])/unit for n in before},
          'references_v':refs,'energy_fj':supply,'event_energy_fj':sum(supply.values()),
          'reference_checks':reference_checks,
          'quiet_reference_drift_units':reference_drift,'quiet_reference_interval_ns':quiet_end-quiet_start,
          'reset_window_energy_fj':sum(reset_work.values()),'reset_window_components_fj':reset_work,
          'nonreset_components_fj':nonreset,
          'charge_algebraic_error_pair':[charge['center']+4*charge['left'],charge['right']-charge['left']],
          'test_state_write_energy_fj':sum(integral(t,np.maximum(pw,0),a*1e-9,b*1e-9)*1e15 for pw in stimulus.values())})
    return {'analysis_schema_version':7,'cfg':cfg,'run':result.get('run'),'mos_count':len(meta['mos']),'observations':obs,
      'states':len(obs),'quotient_correct':sum(o['quotient_correct'] for o in obs),
      'transfer_correct':sum(o['transfer_charge_correct'] for o in obs),'class_correct':sum(o['physical_class_correct'] for o in obs),
      'sleep_correct':sum(o['sleep_correct'] for o in obs),'idle_power_nw':idle_nw,
      'idle_components_fj':idle,'mean_event_fj':float(np.mean([o['event_energy_fj'] for o in obs])),
      'idle_components_nw':{k:val/(idle_b-idle_a)*1000 for k,val in idle.items()},
      'per_component_mean_fj':{k:float(np.mean([o['energy_fj'][k] for o in obs])) for k in obs[0]['energy_fj']}}

def comparator_build(cfg):
    c=Dynamic(cfg);vdd=c.d.vdd;cm=cfg.get('cm',.6)*vdd/1.8
    deltas=cfg.get('deltas',[-.04,-.02,-.01,-.00625,-.003,-.001,0,.001,.003,.00625,.01,.02,.04])
    starts=[50+50*j for j in range(len(deltas))]
    vals=[(0,cm+deltas[0])]
    for s,d in zip(starts,deltas):vals.extend([(s,cm+d),(s+49.9,cm+d)])
    c.source('stim','stim',pwl(vals),'test_state_write');c.source('ref','reference_stim',str(cm),'test_reference')
    acquire=c.clock('acquire',[(s+1,s+10) for s in starts],'sample');ab=c.not_(acquire,'acqb','sample')
    c.tg('acquire','stim','input',acquire,ab,'dd_sample',w=4.2);c.capacitor('input',cfg.get('input_ff',200),'mim')
    c.tg('refacquire','reference_stim','reference',acquire,ab,'dd_sample',w=4.2);c.capacitor('reference',cfg.get('input_ff',200),'mim')
    clock=c.clock('evaluate',[(s+20,s+37) for s in starts])
    op,on=c.comparator('input','reference',clock)
    stop=starts[-1]+200
    keep={'input','reference',op,on,clock}
    c.save={f'v({n})' for n in keep}|set(c.ports)|{f'v({v["node"]})' for v in c.ports.values()}
    meta={'kind':'comparator','cfg':cfg,'design':asdict(c.d),'seed':c.d.seed,'ports':c.ports,'mos':c.mos,
      'clock_node':clock,
      'checks':[{'start_ns':s,'delta_v':d,'eval_ns':s+20,'observe_ns':s+35} for s,d in zip(starts,deltas)],'stop_ns':stop}
    return c.finish(f'.tran .05n {stop}n 0 .1n'),meta

def comparator_analyze(result,names,x):
    meta=result['meta'];ix={n:i for i,n in enumerate(names)};t=x[:,0];dd=meta['design']['vdd']
    if t[-1]<(meta['stop_ns']-.01)*1e-9 or not np.all(np.isfinite(x)):raise ValueError('incomplete/nonfinite comparator transient')
    v=lambda n:x[:,ix['v('+n.lower()+')']]
    allpower={}
    for n,p in meta['ports'].items():allpower[n]=-v(p['node'])*x[:,ix[n]]
    rows=[]
    for p in meta['checks']:
        a=p['start_ns']*1e-9;b=(p['start_ns']+50)*1e-9;at=p['observe_ns']*1e-9
        vp=np.interp(at,t,v('cmpop'));vn=np.interp(at,t,v('cmpon'))
        pos=bool(vp>vn);resolved=bool(max(vp,vn)>.8*dd and min(vp,vn)<.2*dd)
        observed_diff=float(np.interp((p['eval_ns']-.2)*1e-9,t,v('input')-v('reference')))
        period=(t>=p['eval_ns']*1e-9)&(t<=at)
        correctrail=((v('cmpop')>.8*dd)&(v('cmpon')<.2*dd)) if pos else ((v('cmpon')>.8*dd)&(v('cmpop')<.2*dd))
        settled=np.flatnonzero(period&correctrail)
        clk=v(meta.get('clock_node','evaluatedrv2o'))
        crosses=np.flatnonzero((clk[:-1]<dd/2)&(clk[1:]>=dd/2))+1
        candidates=[i for i in crosses if p['eval_ns']*1e-9<=t[i]<(p['eval_ns']+5)*1e-9]
        clock_half=float(np.interp(dd/2,clk[candidates[0]-1:candidates[0]+1],t[candidates[0]-1:candidates[0]+1])) if candidates else None
        after=(t>=p['eval_ns']*1e-9)&(t<(p['eval_ns']+17)*1e-9)
        kick=float(np.max(abs(v('input')[after]-np.interp((p['eval_ns']-.2)*1e-9,t,v('input')))))
        energy={}
        for n,port in meta['ports'].items():
            comp=port['component']
            if comp in ('test_state_write','test_reference','constant_pin'):continue
            power=allpower[n]
            if comp=='clock_pin':power=np.maximum(power,0)
            energy[comp]=energy.get(comp,0)+integral(t,power,a,b)*1e15
        rows.append({**p,'actual_differential_v':observed_diff,'positive':pos,'resolved':resolved,
          'correct':None if abs(p['delta_v'])<1e-12 else bool(resolved and pos==(p['delta_v']>0)),
          'decision_delay_ns':float((t[settled[0]]-clock_half)*1e9) if len(settled) and clock_half is not None else None,
          'kickback_v':kick,'energy_fj':energy,'total_fj':sum(energy.values())})
    idle_a=(meta['stop_ns']-130)*1e-9;idle_b=(meta['stop_ns']-10)*1e-9
    power=sum(p for n,p in allpower.items() if meta['ports'][n]['component'] not in ('test_state_write','test_reference','constant_pin'))
    idle=integral(t,power,idle_a,idle_b)/(idle_b-idle_a)*1e9
    rawcmp=sum(p for n,p in allpower.items() if meta['ports'][n]['component']=='comparator')
    rawcmp_idle=integral(t,rawcmp,idle_a,idle_b)/(idle_b-idle_a)*1e9
    negative=[r['actual_differential_v'] for r in rows if r['resolved'] and not r['positive']]
    positive=[r['actual_differential_v'] for r in rows if r['resolved'] and r['positive']]
    bracket=[max(negative),min(positive)] if negative and positive and max(negative)<min(positive) else None
    thresholds=sorted(set(abs(r['delta_v']) for r in rows if r['delta_v']))
    reliable=[margin for margin in thresholds if all(r['correct'] is True for r in rows if abs(r['delta_v'])>=margin)]
    return {'analysis_schema_version':6,'cfg':meta['cfg'],'kind':'comparator','observations':rows,'idle_power_nw':float(idle),
      'idle_comparator_rail_nw':rawcmp_idle,
      'mean_comparison_fj':float(np.mean([r['total_fj'] for r in rows])),
      'offset_bracket_actual_v':bracket,'offset_bracket_mid_v':sum(bracket)/2 if bracket else None,
      'minimum_tested_correct_abs_delta_v':min(reliable) if reliable else None,
      'correct_nonzero':sum(r['correct'] is True for r in rows),'nonzero':sum(r['correct'] is not None for r in rows),
      'max_kickback_v':max(r['kickback_v'] for r in rows),
      'max_delay_ns':max([r['decision_delay_ns'] for r in rows if r['decision_delay_ns'] is not None],default=None),
      'mos_count':len(meta['mos'])}
