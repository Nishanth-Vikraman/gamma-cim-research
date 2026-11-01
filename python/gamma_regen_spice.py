"""Seven-state q=4 capacitor restoration and genuinely connected analog chains.

No intermediate oracle controls a source. Only six finite-region comparators
(plus two overflow guards) choose one of seven references. Continuous probes
are instrumentation. Ideal buffers and references require separate pricing.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
import json, math, random, subprocess, hashlib
import numpy as np
from gamma_spice import DEFAULT_EXE, pwl, pulses, read_ascii_raw, classify
from gamma_regen_golden import chain_case, word_case
from gamma_native_model import value_pair


@dataclass
class RegenPhysical:
    scale: float = 4
    cap_f: float = 100e-15
    ron_ohm: float = 1000
    reference_sigma_units: float = 0
    reference_gain: float = 0
    supply_gain: float = 0
    comparator_sigma_units: float = 0
    threshold_bank_sigma_units: float = 0
    cap_sigma: float = 0
    clock_feedthrough_f: float = 0
    charge_injection_sigma_c: float = 0
    thermal_kelvin: float = 0
    packet_sigma: float = 0
    width_sigma: float = 0
    mac_gain_sigma: float = 0
    mac_noise_units: float = 0
    leakage_tau_ns: float = 1e12
    stage_drift_units: float = 0
    maxstep_ns: float = .04

    @property
    def unit(self): return .1/self.scale


def profile(name):
    if name=='ideal': return RegenPhysical()
    if name=='low': return RegenPhysical(reference_sigma_units=.012, comparator_sigma_units=.015,
        threshold_bank_sigma_units=.01, cap_sigma=.01, clock_feedthrough_f=.01e-15,
        charge_injection_sigma_c=1e-18, thermal_kelvin=300, packet_sigma=.005,
        width_sigma=.005, mac_gain_sigma=.002, mac_noise_units=.015, leakage_tau_ns=1e6)
    if name=='stress': return RegenPhysical(reference_sigma_units=.12, comparator_sigma_units=.12,
        threshold_bank_sigma_units=.08, cap_sigma=.06, clock_feedthrough_f=.2e-15,
        charge_injection_sigma_c=1e-16, thermal_kelvin=300, packet_sigma=.03,
        width_sigma=.03, mac_gain_sigma=.03, mac_noise_units=.12, leakage_tau_ns=2000)
    raise ValueError(name)


class Circuit:
    def __init__(self, physical, seed):
        self.p=physical; self.rng=random.Random(seed); self.real={}; self.energy=[('i(vcm)',.6)]
        nonideal=any(getattr(physical,k) for k in ('reference_sigma_units','reference_gain','supply_gain',
            'comparator_sigma_units','threshold_bank_sigma_units','cap_sigma','clock_feedthrough_f',
            'charge_injection_sigma_c','thermal_kelvin','packet_sigma','width_sigma','mac_gain_sigma','mac_noise_units'))
        nonideal=nonideal or physical.ron_ohm!=1000 or physical.leakage_tau_ns<1e12
        # Gear avoids trapezoidal ringing in tiny parasitic/control capacitors.
        # The method is chosen from model configuration, never from its answer.
        options='reltol=1e-7 abstol=1e-16 vntol=1e-10 method=gear' if nonideal else 'reltol=1e-8 abstol=1e-18 vntol=1e-11'
        self.lines=['q4 finite-region canonical regeneration',
            '* Behavioral references/buffers are ideal; port energy is not full supply energy.',
            f'.model WRITE SW(Ron={physical.ron_ohm:.12g} Roff=1e18 Vt=0.5 Vh=0.05)',
            '.model SAMPLE SW(Ron=1 Roff=1e18 Vt=0.5 Vh=0.05)',
            '.model LATCH SW(Ron=10000 Roff=1e18 Vt=0.5 Vh=0.05)',
            '.options plotwinsize=0 numdgt=15 '+options,
            'Vcm cm 0 .6','Vpacket packet 0 1.2']
        self.save={'v(cm)'}

    def normal(self, sigma): return self.rng.gauss(0,sigma)

    def clock(self,name,start,width):
        self.lines.append(f'V{name} {name} 0 {pwl(pulses([start],width))}')
        return name

    def cap(self,name,nominal,initial=.6):
        cap=nominal*(1+self.normal(self.p.cap_sigma))
        if cap<=0: raise ValueError('nonphysical capacitance')
        self.lines += [f'C{name} {name} 0 {cap:.14g} IC={initial:.14g}',
            f'R{name} {name} cm {self.p.leakage_tau_ns*1e-9/cap:.14g}']
        return cap

    def inject(self,name,node,start,charge):
        self.lines.append(f'I{name} 0 {node} {pwl(pulses([start],.2,charge/(.19e-9)))}')

    def regenerate(self, tag, signal, rawunit, start, mapping='differential', mode='regen', invalid='0'):
        """Sample finite decision at start; reset; write; isolate; probe at +5ns.
        Returns differential/centered actual capacitor output buffered for load.
        """
        p=self.p; u=p.unit; rng=self.rng
        sample=self.clock(tag+'s',start,.2)
        reset=self.clock(tag+'r',start+.5,.5)
        width=2*(1+self.normal(p.width_sigma))
        write=self.clock(tag+'w',start+1.2,width)
        common=self.normal(p.threshold_bank_sigma_units)
        offsets=[common+self.normal(p.comparator_sigma_units) for _ in range(8)]
        z=self.normal(p.mac_noise_units)
        normalized=f'(({signal})/{rawunit:.14g}+({z:.14g}))'
        code='-3' + ''.join(f'+if({normalized}>{t+offsets[i]:.14g},1,0)' for i,t in enumerate((-2.5,-1.5,-.5,.5,1.5,2.5)))
        bad=f'if({normalized}<-3.5+({offsets[6]:.14g}),1,if({normalized}>=3.5+({offsets[7]:.14g}),1,if(({invalid})>.5,1,0)))'
        self.lines += [f'B{tag}code {tag}code 0 V={code}',
            f'B{tag}bad {tag}bad 0 V={bad}',
            f'S{tag}code {tag}code {tag}h {sample} 0 LATCH',
            f'C{tag}h {tag}h 0 1f IC=0',f'R{tag}h {tag}h 0 1e18',
            f'S{tag}bad {tag}bad {tag}bh {sample} 0 LATCH',
            f'C{tag}bh {tag}bh 0 1f IC=0',f'R{tag}bh {tag}bh 0 1e18']
        self.save.update([f'v({tag}h)',f'v({tag}bh)'])
        halves=1 if mapping=='single' else 2
        nodes=[]; caps=[]
        refs={d:u*(d*(1+p.reference_gain+p.supply_gain)+self.normal(p.reference_sigma_units)) for d in range(-3,4)}
        gain=(1+self.normal(p.packet_sigma))*(1+p.reference_gain+p.supply_gain+self.normal(p.reference_sigma_units))
        self.real[tag]={'threshold_offsets_units':offsets,'references_v':refs,'packet_gain':gain,'width_ns':width}
        for j in range(halves):
            n=f'{tag}o{j}'; nodes.append(n)
            cap=self.cap(n,p.cap_f/halves); caps.append(cap)
            self.lines.append(f'S{tag}reset{j} {n} cm {reset} 0 WRITE')
            if p.clock_feedthrough_f:
                self.lines.append(f'C{tag}feed{j} {n} {write} {p.clock_feedthrough_f:.14g} IC=.6')
            charge=self.normal(p.charge_injection_sigma_c)
            charge+=self.normal(math.sqrt(1.380649e-23*p.thermal_kelvin*cap))
            if charge: self.inject(tag+'kick'+str(j),n,start+3.5,charge)
            factor=1 if halves==1 else (.5 if j==0 else -.5)
            if mode=='none':
                # Known physical scale compensation, preserving residual, not an ADC.
                # Differential headroom +/-0.3 V; a railing direct baseline
                # must not grow to physically meaningless kilovolt levels.
                rail=.3 if halves==1 else .15
                self.lines += [f'B{tag}direct{j} {tag}direct{j} 0 V=.6+limit(({signal})*{u/rawunit*factor:.14g},-{rail},{rail})',
                    f'S{tag}direct{j} {tag}direct{j} {n} {write} 0 WRITE']
            elif mapping=='charge':
                # Three thermometer-selected equal charge packets per branch.
                amp=(p.cap_f/2)*u*.5/(1.99e-9)*gain
                for m in (1,2,3):
                    expr=f'if(abs(V({tag}h))>{m-.5},if(V({tag}h)>0,1,-1),0)*V({write})*if(V({tag}bh)<.5,1,0)'
                    self.lines.append(f'B{tag}packet{j}{m} packet {n} I={amp*(1 if j==0 else -1):.14g}*({expr})')
                if ('i(vpacket)',1.2) not in self.energy:self.energy.append(('i(vpacket)',1.2))
            else:
                # Seven-entry reference selector. B source is an ideal finite
                # mux, not a general arithmetic DAC: each branch is constant.
                # The explicit Ron switch charges a real capacitor afterward.
                selected='.6'
                for d in range(-3,4):
                    selected=f'if(abs(V({tag}h)-({d}))<.25,{.6+factor*refs[d]:.14g},{selected})'
                ref=f'{tag}reference{j}';control=f'{tag}valid{j}'
                self.lines += [f'B{ref} {ref} 0 V=if(V({write})>.1,{selected},.6)',
                    f'B{control} {control} 0 V=if(V({tag}bh)<.5,V({write}),0)',
                    f'S{ref} {ref} {n} {control} 0 WRITE']
                self.energy.append((f'i(b{ref})',ref))
                self.save.add(f'v({ref})')
        out=tag+'out'
        expr=f'V({nodes[0]})-.6' if halves==1 else f'V({nodes[0]},{nodes[1]})'
        self.lines.append(f'B{tag}buffer {out} 0 V={expr}')
        self.save.add(f'v({out})')
        self.save.update(f'v({n})' for n in nodes)
        return out,tag+'bh',start+5

    def finish(self,stop):
        self.save.update(name for name,_ in self.energy)
        return '\n'.join(self.lines+['.save '+' '.join(sorted(self.save)),
            f'.tran 0 {stop:.12g}n 0 {self.p.maxstep_ns:.12g}n UIC','.end'])+'\n'


def primitive_deck(digit=1,mapping='differential',physical=None,seed=1,residual=.3):
    p=physical or RegenPhysical(); c=Circuit(p,seed)
    c.lines.append(f'Vin input 0 {(digit+residual)*p.unit:.14g}')
    out,bad,t=c.regenerate('r','V(input)',p.unit,1,mapping)
    metadata={'kind':'primitive','digit':digit,'mapping':mapping,'checkpoints':[
        {'time_ns':t,'out':out,'bad':bad,'expected':digit,'raw':None}], 'seed':seed,
        'physical':asdict(p),'realization':c.real,'energy_ports':c.energy,'unit':p.unit}
    return c.finish(t+.2),metadata


def chain_deck(digit=1,slope=2,depth=16,mode='regen',physical=None,seed=1):
    p=physical or RegenPhysical(); c=Circuit(p,seed)
    oracle=chain_case(digit,slope,depth); w1=float(slope); w2=1-w1; scale=abs(w1)+abs(w2)
    c.lines += [f'Vinput input 0 {digit*p.unit:.14g}']
    previous='input'; bad='0'; checks=[]
    for stage in range(depth):
        tag=f's{stage}'; start=stage*8+1; raw=tag+'raw'
        c.cap(raw,p.cap_f,initial=0)
        # Override the common-mode leakage used by storage caps: raw states centered at zero.
        c.lines[-1]=f'R{raw} {raw} 0 {p.leakage_tau_ns*1e-9/p.cap_f:.14g}'
        gate=c.clock(tag+'mac',start,.8)
        gm=p.cap_f/(.79e-9*scale)
        g1=1+c.normal(p.mac_gain_sigma);g2=1+c.normal(p.mac_gain_sigma)
        c.lines += [f'B{tag}learned 0 {raw} I={gm:.14g}*({w1*g1:.14g}*V({previous})+({w2*g2:.14g})*V(input))*V({gate})']
        if p.stage_drift_units: c.inject(tag+'drift',raw,start+.9,p.cap_f*p.unit/scale*p.stage_drift_units)
        c.save.add(f'v({raw})')
        out,badnode,t=c.regenerate(tag,f'V({raw})',p.unit/scale,start+1.3,mode=mode,invalid=f'V({bad})' if bad!='0' else '0')
        checks.append({'stage':stage+1,'time_ns':t,'out':out,'bad':badnode,'expected':digit,'raw':raw})
        previous=out;bad=badnode
    return c.finish(checks[-1]['time_ns']+.2), {'kind':'chain','digit':digit,'slope':str(slope),'mode':mode,
        'checkpoints':checks,'seed':seed,'physical':asdict(p),'realization':c.real,'unit':p.unit,
        'energy_ports':c.energy,'oracle':oracle,'rawunit':p.unit/scale}


def normalizer(c,tag,start,unit,nodes=range(-2,3),rounds=32):
    allnodes=range(min(nodes)-1,max(nodes)+2)
    for k in allnodes:
        n=tag+('m'+str(-k) if k<0 else 'p'+str(k));c.cap(n,c.p.cap_f,0)
        c.lines[-1]=f'R{n} {n} 0 1e18';c.save.add(f'v({n})')
    def name(k):return tag+('m'+str(-k) if k<0 else 'p'+str(k))
    c.lines += [f'V{tag}sample {tag}sample 0 {pwl(pulses([start+4*j for j in range(rounds)],.2))}',
        f'V{tag}fire {tag}fire 0 {pwl(pulses([start+4*j+.6 for j in range(rounds)],1))}']
    current=c.p.cap_f*unit/.99e-9
    for k in nodes:
        n=name(k); h=n+'h'
        c.lines += [f'B{n}cmd {n}cmd 0 V=if(V({n})>{3.5*unit:.14g},1,if(V({n})<-{3.5*unit:.14g},-1,0))',
            f'S{n}h {n}cmd {h} {tag}sample 0 SAMPLE',f'C{h} {h} 0 1f IC=0',f'R{h} {h} 0 1e18',
            f'B{n}remove {n} 0 I={4*current:.14g}*V({h})*V({tag}fire)',
            f'B{n}left 0 {name(k-1)} I={current:.14g}*V({h})*V({tag}fire)',
            f'B{n}right 0 {name(k+1)} I={current:.14g}*V({h})*V({tag}fire)']
    return name


def word_deck(a=1,b=2,bias=1,mode='regen',physical=None,seed=1):
    p=physical or RegenPhysical();c=Circuit(p,seed);u=p.unit; rawunit=u/3
    gold=word_case(a,b,bias);s1=normalizer(c,'a',5,rawunit);s2=normalizer(c,'b',150,rawunit)
    c.lines += [f'VA inputa 0 {a*u:.14g}',f'VB inputb 0 {b*u:.14g}',f'VBias inputbias 0 {bias*u:.14g}']
    g=c.clock('mac1',1,1);gm=p.cap_f/(.99e-9*3)
    c.lines.append(f'Bmac1 0 {s1(0)} I={gm:.14g}*(2*V(inputa)+V(inputb))*V({g})')
    if p.stage_drift_units:c.inject('stage1residual',s1(0),3,p.cap_f*rawunit*p.stage_drift_units)
    outputs={}; checks=[]
    for k in range(-3,4):
        tag='r'+str(k+3)
        out,bad,t=c.regenerate(tag,f'V({s1(k)})',rawunit,135,mode=mode)
        outputs[k]=out
        checks.append({'stage':1,'node':k,'time_ns':t,'out':out,'bad':bad,
            'expected':gold['stage1']['state'].get(k,0),'raw':s1(k)})
    g=c.clock('mac2',143,1)
    for k in range(-3,4):
        biasexpr='-V(inputbias)' if k==0 else '0'
        c.lines.append(f'Bmac2{k+3} 0 {s2(k)} I={gm:.14g}*(2*V({outputs[k]})+({biasexpr}))*V({g})')
    for k in range(-3,4):
        checks.append({'stage':2,'node':k,'time_ns':279,'out':s2(k),'bad':None,
            'expected':gold['stage2']['state'].get(k,0),'raw':None,'unit':rawunit})
    return c.finish(280),{'kind':'word','mode':mode,'seed':seed,'physical':asdict(p),
        'realization':c.real,'unit':u,'rawunit':rawunit,'oracle':gold,
        'energy_ports':c.energy,'checkpoints':checks}


def boundary_deck():
    """Fixed index-copy configuration selected by digital exponent sidecar.
    The low outward term reaches an explicit discarded-charge receiver.
    """
    p=RegenPhysical();c=Circuit(p,0);u=p.unit
    clk=c.clock('shift',1,2)
    c.lines += [f'Vhigh high 0 {u}',f'Vlow low 0 {u}']
    for n,src in [('kept','high'),('discard','low')]:
        c.cap(n,p.cap_f,0);c.lines[-1]=f'R{n} {n} 0 1e18'
        c.lines.append(f'S{n} {src} {n} {clk} 0 WRITE');c.save.add(f'v({n})')
    checks=[{'time_ns':4,'out':'kept','bad':None,'expected':1,'raw':None},
        {'time_ns':4,'out':'discard','bad':None,'expected':1,'raw':None}]
    return c.finish(4.2),{'kind':'boundary','seed':0,'unit':u,'physical':asdict(p),'energy_ports':c.energy,
        'checkpoints':checks,'index_action':'source k=3 to destination k=2; low source k=-2 to discard k=-3',
        'digital_exponent_delta':1,'finite_oracle':__import__('gamma_regen_golden').finite_word({3:1,-2:1})}


def fast_raw(path):
    data=Path(path).read_bytes()
    txt=data.decode('utf-16') if data[:2] in (b'\xff\xfe',b'\xfe\xff') else data.decode('utf-8-sig')
    head,body=txt.split('Values:',1)
    names=[line.split()[1].lower() for line in head.rsplit('Variables:',1)[1].strip().splitlines()]
    # First line of each record has point index + time; subsequent lines have
    # one value. Keeping the index column allows fast exact schema parsing.
    array=np.fromstring(body,sep=' ')
    if len(array)%(len(names)+1):raise ValueError('incomplete raw record')
    return names,array.reshape(-1,len(names)+1)[:,1:]


def run(deck,meta,directory,retain_raw=False,timeout=30):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    path=directory/'cell.cir';path.write_text(deck,encoding='ascii')
    (directory/'golden.json').write_text(json.dumps(meta,indent=2),encoding='ascii')
    result=subprocess.run([str(DEFAULT_EXE),'-b','-ascii',str(path.resolve())],capture_output=True,
        timeout=timeout,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    (directory/'process.log').write_bytes(result.stdout+result.stderr)
    if result.returncode or not path.with_suffix('.raw').exists():
        raise RuntimeError(f'LTspice failed: {directory}, code={result.returncode}')
    names,data=fast_raw(path.with_suffix('.raw'));lookup={n:i for i,n in enumerate(names)}
    time=data[:,0]; observations=[]
    for check in meta['checkpoints']:
        index=np.searchsorted(time,check['time_ns']*1e-9)
        v=float(data[index,lookup['v('+check['out'].lower()+')']]);unit=check.get('unit',meta['unit'])
        cls=classify(v,unit,3)
        bad=bool(check['bad'] and data[index,lookup['v('+check['bad'].lower()+')']]>.5)
        raw=float(data[index,lookup['v('+check['raw'].lower()+')']]) if check.get('raw') else None
        observations.append(dict(check,voltage_v=v,decoded_class=cls,overflow_latched=bad,
            correct=cls==check['expected'] and not bad,residual_units=v/unit-check['expected'],raw_v=raw))
    energy=0.;absorbed=0.
    for name,voltage in meta['energy_ports']:
        v=data[:,lookup['v('+voltage.lower()+')']] if isinstance(voltage,str) else voltage
        power=-v*data[:,lookup[name.lower()]]
        energy+=float(np.trapezoid(np.maximum(power,0),time))
        absorbed+=float(np.trapezoid(np.maximum(-power,0),time))
    # Settling wrt final plateau, not nominal reference; drift/reference errors remain reported separately.
    first=meta['checkpoints'][0];target=observations[0]['voltage_v']
    mask=(time>=(first['time_ns']-5)*1e-9)&(time<=first['time_ns']*1e-9)
    tv=time[mask];vv=data[mask,lookup['v('+first['out'].lower()+')']]
    outside=np.where(abs(vv-target)>.01*meta['unit'])[0]
    settle=float((tv[min(outside[-1]+1,len(tv)-1)]-tv[0])*1e9) if len(outside) else 0.
    report={'kind':meta['kind'],'seed':meta['seed'],'observations':observations,
        'reference_port_delivered_energy_fj':energy*1e15,'reference_port_absorbed_energy_fj':absorbed*1e15,
        'first_output_settling_ns_to_1pct_class_spacing':settle,
        'energy_scope':'reference ports only; ideal comparator/buffer/clock/reference generation excluded',
        'deck_sha256':hashlib.sha256(deck.encode()).hexdigest(),'physical':meta['physical']}
    if meta['kind']=='word':
        s1=[o for o in observations if o['stage']==1]
        report['stage1_raw_residual_units_max']=max(abs(o['raw_v']/meta['rawunit']-o['expected']) for o in s1)
        report['stage1_canonical_residual_units_max']=max(abs(o['residual_units']) for o in s1)
        report['stage2_residual_units_max']=max(abs(o['residual_units']) for o in observations if o['stage']==2)
        final={o['node']:o['decoded_class'] for o in observations if o['stage']==2}
        if all(isinstance(v,int) for v in final.values()):
            pair=value_pair(final,4);report['final_exact_pair']=[str(v) for v in pair]
            report['final_algebraically_correct']=report['final_exact_pair']==meta['oracle']['stage2_exact']
            report['exact_arithmetic_error_pair']=[str(x-__import__('fractions').Fraction(y)) for x,y in zip(pair,meta['oracle']['stage2_exact'])]
            report['equivalent_different_word']=report['final_algebraically_correct'] and any(not o['correct'] for o in observations if o['stage']==2)
        else:report.update(final_algebraically_correct=False,equivalent_different_word=False)
    (directory/'result.json').write_text(json.dumps(report,indent=2),encoding='ascii')
    if not retain_raw:path.with_suffix('.raw').unlink()
    return report
