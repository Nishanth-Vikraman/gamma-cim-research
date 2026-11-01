"""q=4 SKY130 primitive experiments. PDK content is never patched.

Class decisions are explicitly ideal or hybrid; all analog write/buffer paths
use foundry MOS devices. Observations and exact oracles never drive later nodes.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
from pathlib import Path
import ctypes,hashlib,json,math,os,random,subprocess,time,gzip,shutil
import numpy as np
from gamma_regen_spice import fast_raw as ascii_raw
from gamma_spice import classify

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'spice/sky130'
PDK=ROOT/'third_party/sky130-ngspice-models'
EXE=ROOT/'third_party/ngspice/Spice64/bin/ngspice_con.exe'
NF='sky130_fd_pr__nfet_01v8'
PF='sky130_fd_pr__pfet_01v8'
CORNERS=('tt','ss','ff','sf','fs')


def fast_raw(path):
    path=Path(path)
    if path.exists():data=path.read_bytes()
    else:
        with gzip.open(str(path)+'.gz','rb') as f:data=f.read()
    marker=next((s for s in (b'Binary:\n',b'Binary:\r\n') if s in data),None)
    if marker:
        header,body=data.split(marker,1);head=header.decode('utf-8-sig')
        names=[line.split()[1].lower() for line in head.rsplit('Variables:',1)[1].strip().splitlines()]
        values=np.frombuffer(body,dtype='<f8')
        if len(values)%len(names):raise ValueError('truncated binary waveform')
        return names,values.reshape(-1,len(names))
    text=data.decode('utf-16') if data[:2] in (b'\xff\xfe',b'\xfe\xff') else data.decode('utf-8-sig')
    head,body=text.split('Values:',1)
    names=[line.split()[1].lower() for line in head.rsplit('Variables:',1)[1].strip().splitlines()]
    array=np.fromstring(body,sep=' ')
    if len(array)%(len(names)+1):raise ValueError('truncated ASCII waveform')
    return names,array.reshape(-1,len(names)+1)[:,1:]


def compress_trace(path):
    path=Path(path)
    if not path.exists():return
    target=Path(str(path)+'.gz')
    with path.open('rb') as src,gzip.open(target,'wb',compresslevel=4) as dst:shutil.copyfileobj(src,dst)
    # CRC and size checked during complete decompression before removing input.
    h1=hashlib.sha256(path.read_bytes()).digest()
    with gzip.open(target,'rb') as f:h2=hashlib.sha256(f.read()).digest()
    if h1!=h2:raise RuntimeError('Lossless trace compression failed')
    path.unlink()


def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2),encoding='utf-8')


def short(path):
    path=Path(path).resolve()
    if os.name!='nt':return path.as_posix()
    buf=ctypes.create_unicode_buffer(32768)
    n=ctypes.windll.kernel32.GetShortPathNameW(str(path),buf,len(buf))
    return (buf.value if n else str(path)).replace('\\','/')


def pwl(points):
    return 'PWL('+' '.join(f'{t:.8g}n {v:.12g}' for t,v in sorted(dict(points).items()))+')'


@dataclass
class Design:
    corner:str='tt'
    vdd:float=1.8
    temp_c:float=25
    cap_ff:float=100
    cap_kind:str='ideal'
    tg_w:float=.84
    input_w:float=4
    buffer_scale:float=1
    pn_ratio:float=2
    buffer_kind:str='two_stage'
    bias_r:float=2e6
    ladder_r:float=5e5
    load_r:float=2e5
    load_cap_ff:float=20
    period_ns:float=200
    write_buffer:bool=False
    gated_buffer:bool=False
    detector:str='ideal'
    comparator_w:float=1
    pdk_mismatch:bool=False
    resistor_sigma:float=0
    reference_sigma_units:float=0
    capacitor_sigma:float=0
    systematic_offset_mv:float=0
    seed:int=101
    maxstep_ns:float=.25
    fanout:int=1
    input_scale:float=12
    selector_kind:str='tg'
    comparator_l:float=.15


class Circuit:
    def __init__(self,design):
        if not 0<design.vdd<=1.95:
            raise ValueError('SKY130 core experiment restricted to supply <= documented1.95 V model limit')
        self.d=design;self.rng=random.Random(design.seed);self.mos=[];self.refs={};self.ports={}
        self.lines=['Gamma q4 SKY130 transistor boundary',f'.lib "{short(PDK/"libs.tech/ngspice/sky130.lib.spice")}" {design.corner}',
            f'.param mc_mm_switch={int(design.pdk_mismatch)} mc_pr_switch=0',
            f'.temp {design.temp_c}', '.option reltol=1e-5 abstol=1e-14 vntol=1e-8 method=gear',
            '* Dimensions of SKY130 X instances are micrometers, not SI meters.']
        self.save=set()

    def source(self,name,node,voltage,component):
        self.lines.append(f'V{name} {node} 0 {voltage}')
        self.save.update((f'v({node})',f'i(v{name})'))
        self.ports[f'i(v{name.lower()})']={'node':node.lower(),'component':component}

    def supply(self,component):
        node='dd_'+component
        self.source(component,node,f'{self.d.vdd:.12g}',component)
        return node

    def fet(self,name,kind,drain,gate,source,body,w,l=.15):
        self.lines.append(f'X{name} {drain} {gate} {source} {body} {NF if kind=="n" else PF} w={w:.10g} l={l:.10g} ad={w*.29:.10g} as={w*.29:.10g} pd={2*(w+.29):.10g} ps={2*(w+.29):.10g}')
        self.mos.append({'name':name,'type':kind,'d':drain,'g':gate,'s':source,'b':body,'w_um':w,'l_um':l})
        self.save.update(f'v({node})' for node in (drain,gate,source,body) if node!='0')

    def resistor(self,name,a,b,value,mismatch=True):
        value*=1+self.rng.gauss(0,self.d.resistor_sigma if mismatch else 0)
        if value<=0:raise ValueError('nonphysical resistor realization')
        self.lines.append(f'R{name} {a} {b} {value:.14g}')

    def inv(self,name,vin,out,vdd,w=.42,ratio=2):
        self.fet(name+'n','n',out,vin,'0','0',w)
        self.fet(name+'p','p',out,vin,vdd,vdd,w*ratio)

    def tg(self,name,a,b,control,complement,vdd,w=None):
        w=w or self.d.tg_w
        self.fet(name+'n','n',a,control,b,'0',w)
        if self.d.selector_kind=='tg':self.fet(name+'p','p',a,complement,b,vdd,w*self.d.pn_ratio)

    def logic(self,name,source,component='interface',scale=1):
        vdd='dd_'+component
        self.inv(name+'1',source,name+'b',vdd,w=.42*scale)
        self.inv(name+'2',name+'b',name+'q',vdd,w=.42*scale)
        return name+'q',name+'b'

    def capacitor(self,node,cap_ff=None,kind=None):
        cap=(cap_ff or self.d.cap_ff)*1e-15
        cap*=1+self.rng.gauss(0,self.d.capacitor_sigma)
        if (kind or self.d.cap_kind)=='ideal':
            self.lines.append(f'C{node} {node} 0 {cap:.14g}')
        else:
            # Typical PDK: C=2 fF/um2*(side-.025)^2 + .76 fF/um*(side-.025).
            count=max(1,math.ceil(cap/1e-12))
            side=(-.76+math.sqrt(.76**2+8*cap*1e15/count))/4+.025
            for j in range(count):
                self.lines.append(f'Xc{node}_{j} {node} 0 sky130_fd_pr__cap_mim_m3_1 w={side:.12g} l={side:.12g}')
        self.save.add(f'v({node})')

    def ladder(self):
        dd=self.supply('reference');unit=self.d.ladder_r/144
        # 13 taps: seven integer states and six midpoint thresholds.
        self.resistor('bottom','tap0','0',42*unit)
        for k in range(12):self.resistor('lad'+str(k),'tap'+str(k),'tap'+str(k+1),unit)
        self.resistor('top',dd,'tap12',90*unit)
        for k in range(13):
            self.save.add('v(tap'+str(k)+')')
        for digit in range(-3,4):self.refs[digit]='tap'+str(2*(digit+3))
        # Architecture-level per-reference perturbations are resistor changes,
        # never ideal sources supplying the selected write path.
        if self.d.reference_sigma_units:
            raise ValueError('use resistor_sigma for physical ladder mismatch; no ideal offset reference sources')
        return self.refs

    def amplifier(self,tag,vin,out,component='buffer',enable=None,feedback=None):
        dd='dd_'+component;d=self.d;s=d.buffer_scale
        # Resistor-biased PMOS diode + mirror: its current is supply-accounted.
        self.fet(tag+'bias','p',tag+'pb',tag+'pb',dd,dd,2,1)
        self.resistor(tag+'bias',tag+'pb','0',d.bias_r,False)
        self.fet(tag+'tail','p',tag+'tail',tag+'pb',dd,dd,2*s,1)
        tail=tag+'tail'
        if enable:
            self.fet(tag+'enable','p',tag+'et',enable,tail,dd,max(2,8*s),.15);tail=tag+'et'
        first=out if d.buffer_kind=='five_transistor' else tag+'first'
        # Single stage has positive gain, two stage adds inversion.
        feedback=feedback or out
        g1=vin if d.buffer_kind=='five_transistor' else feedback
        g2=feedback if d.buffer_kind=='five_transistor' else vin
        self.fet(tag+'i1','p',tag+'mirror',g1,tail,dd,d.input_w,.5)
        self.fet(tag+'i2','p',first,g2,tail,dd,d.input_w,.5)
        self.fet(tag+'n1','n',tag+'mirror',tag+'mirror','0','0',max(.42,2*s),1)
        self.fet(tag+'n2','n',first,tag+'mirror','0','0',max(.42,2*s),1)
        if d.buffer_kind!='five_transistor':
            ps=dd
            if enable:
                self.fet(tag+'oenable','p',tag+'os',enable,dd,dd,max(4,16*s),.15);ps=tag+'os'
            self.fet(tag+'pload','p',out,tag+'pb',ps,dd,6*s,1)
            self.fet(tag+'nout','n',out,first,'0','0',max(.42,8*s/d.pn_ratio),.5)
            self.lines.append(f'C{tag}miller {first} {out} {max(3,8*s):.10g}f')
        self.save.update((f'v({vin})',f'v({out})'))

    def comparator(self,vin,reference,clock,tag='cmp'):
        dd='dd_comparator';d=self.d;w=d.comparator_w
        # Complementary StrongARM: PMOS input pair, outputs reset to ground.
        self.inv(tag+'clk',clock,tag+'clkb',dd)
        self.fet(tag+'tail','p',tag+'tail',tag+'clkb',dd,dd,4*w)
        self.fet(tag+'in1','p',tag+'a',reference,tag+'tail',dd,w,d.comparator_l)
        self.fet(tag+'in2','p',tag+'b',vin,tag+'tail',dd,w,d.comparator_l)
        self.fet(tag+'lp','p',tag+'op',tag+'on',tag+'a',dd,2*w)
        self.fet(tag+'ln','p',tag+'on',tag+'op',tag+'b',dd,2*w)
        self.fet(tag+'np','n',tag+'op',tag+'on','0','0',w)
        self.fet(tag+'nn','n',tag+'on',tag+'op','0','0',w)
        for node in ('op','on','a','b'):
            self.fet(tag+'reset'+node,'n',tag+node,tag+'clkb','0','0',w)
        self.lines += [f'C{tag}outp {tag}op 0 2f',f'C{tag}outn {tag}on 0 2f']
        return tag+'op',tag+'on'

    def finish(self,analysis):
        self.save.discard('v(0)')
        return ('\n'.join(self.lines+['.save '+' '.join(sorted(self.save)),analysis,'.end'])+'\n').replace('if(','ternary_fcn(')


def execute(deck,meta,directory,timeout=60,retain=True):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    harness='Sky130 reproducible seeded driver\n.control\nset noaskquit\nset num_threads=1\nset filetype=ascii\nsetseed '+str(meta.get('seed',101))+'\nsource circuit.cir\nrun\nwrite output.raw all\nquit\n.endc\n.end\n'
    fingerprint=hashlib.sha256((deck+harness+json.dumps(meta,sort_keys=True)+'parser_v1').encode()).hexdigest()
    cache=directory/'result.json'
    if cache.exists():
        old=json.loads(cache.read_text())
        if old.get('fingerprint')==fingerprint and (directory/'output.raw').exists():
            names,data=fast_raw(directory/'output.raw');return old,names,data
    (directory/'circuit.cir').write_text(deck,encoding='ascii')
    (directory/'run.cir').write_text(harness,encoding='ascii');dump(directory/'golden.json',meta)
    start=time.monotonic()
    init=EXE.parents[2]/'minimal_init';init.mkdir(exist_ok=True)
    (init/'spinit').write_text('set num_threads=1\nset filetype=ascii\n',encoding='ascii')
    env=dict(os.environ,OMP_NUM_THREADS='1',SPICE_SCRIPTS=short(init))
    r=subprocess.run([short(EXE),'-b','-o','simulator.log','run.cir'],cwd=short(directory),env=env,
        capture_output=True,timeout=timeout,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    (directory/'process.log').write_bytes(r.stdout+r.stderr)
    logfile=directory/'simulator.log';log=logfile.read_text(errors='replace') if logfile.exists() else ''
    if r.returncode or not (directory/'output.raw').exists():
        dump(directory/'failure.json',{'returncode':r.returncode,'log_tail':log[-3000:],'fingerprint':fingerprint})
        raise RuntimeError(f'ngspice failed at {directory}: {log[-1800:]}')
    names,data=fast_raw(directory/'output.raw')
    result={'fingerprint':fingerprint,'returncode':r.returncode,'seconds':time.monotonic()-start,
        'rows':len(data),'seed':meta.get('seed',101),'meta':meta}
    dump(cache,result)
    return result,names,data


def smoke_deck(kind,corner='tt'):
    d=Design(corner=corner);c=Circuit(d);vdd=c.supply('smoke')
    if kind=='nmos':
        c.source('gate','gate',0,'stimulus');c.source('drain','drain',1,'smoke_drain')
        c.fet('n','n','drain','gate','0','0',1,.15);analysis='.dc Vgate 0 1.8 .01'
    elif kind=='pmos':
        c.source('gate','gate',1.8,'stimulus');c.source('drain','drain',.8,'smoke_drain')
        c.fet('p','p','drain','gate',vdd,vdd,1,.15);analysis='.dc Vgate 1.8 0 -.01'
    elif kind=='inverter':
        c.source('input','input','PULSE(0 1.8 5n .2n .2n 10n 20n)','stimulus')
        c.inv('i','input','out',vdd);c.lines.append('Cload out 0 20f');c.save.add('v(out)')
        analysis='.tran .05n 45n'
    elif kind=='transmission_gate':
        c.source('input','input',.675,'input');c.source('control','control','PULSE(0 1.8 5n .2n .2n 15n 40n)','stimulus')
        c.inv('i','control','cb',vdd);c.tg('switch','input','hold','control','cb',vdd)
        c.capacitor('hold');c.lines.append('Rhold hold 0 1e15');analysis='.tran .05n 30n'
    else:raise ValueError(kind)
    return c.finish(analysis),{'kind':kind,'corner':corner,'seed':101,'mos':c.mos,'ports':c.ports}


def run_smokes():
    results=[]
    for corner in CORNERS:
        for kind in ('nmos','pmos','inverter','transmission_gate'):
            result,names,a=execute(*smoke_deck(kind,corner),BASE/'smoke'/f'{corner}_{kind}')
            ix={n:i for i,n in enumerate(names)}
            if kind in ('nmos','pmos'):
                signal=-a[:,ix['i(vdrain)']] if kind=='nmos' else -a[:,ix['i(vsmoke)']]
                target=1 if kind=='nmos' else .8;g=a[:,ix['v(gate)']]
                current=float(signal[np.argmin(abs(g-target))]);passed=1e-6<current<2e-3
                detail={'id_at_vgs_1v_a':current,'off_current_a':float(signal[0]),'max_current_a':float(max(signal))}
            elif kind=='inverter':
                inp=a[:,ix['v(input)']];v=a[:,ix['v(out)']];t=a[:,0]
                mask=(t>6e-9)&(t<14e-9);passed=bool(v[mask].max()<.1)
                detail={'output_min_v':float(min(v)),'output_max_v':float(max(v)),'peak_supply_a':float(max(-a[:,ix['i(vsmoke)']]))}
            else:
                voltage=float(a[np.searchsorted(a[:,0],25e-9),ix['v(hold)']]);passed=abs(voltage-.675)<.00625
                detail={'held_voltage_v':voltage,'peak_supply_a':float(max(-a[:,ix['i(vsmoke)']]))}
            row={'corner':corner,'kind':kind,'passed':passed,**detail};results.append(row)
            print('smoke',corner,kind,passed,detail,flush=True)
    dump(BASE/'results/smoke_summary.json',results)
    if not all(r['passed'] for r in results):raise RuntimeError('Smoke validation failed; gamma circuits blocked')
    return results


DEFAULT_SEQUENCE=(-3,3,-3,0,3,2,-2,-1,1,0,2,3)


def boundary_deck(d=Design(),sequence=DEFAULT_SEQUENCE,residual=.35,direct=False,threshold_cap_ff=1000):
    c=Circuit(d);refs=c.ladder()
    for component in ('interface','buffer','writer','comparator'):c.supply(component)
    t=d.period_ns;warmup=max(1000,4*d.ladder_r*d.cap_ff*1e-6)
    u=d.vdd/72;cm=d.vdd/3;input_unit=u*4/d.input_scale
    levels=[(0,cm-3*input_unit)];bias=[(0,cm+3*u)]
    write=[(0,d.vdd)];reset=[(0,0)];evals=[(0,0)];awake=[(0,0)]
    checks=[];previous=-3
    for j,digit in enumerate(sequence):
        start=warmup+j*t
        levels += [(start,cm+(previous+residual)*input_unit),(start+1,cm+(digit+residual)*input_unit)]
        bias += [(start,cm-previous*u),(start+1,cm-digit*u)]
        checks.append({'time_ns':start+.90*t,'write_end_ns':start+.65*t,'start_ns':start,
            'previous_digit':previous,'digit':digit,'exact_gamma_pair':[str(digit),'0']})
        previous=digit
    for points,on,off in ((write,.18,.65),(reset,.04,.10),(evals,.11,.98),(awake,.60,.98)):
        points.extend([(warmup-3,0)])
        for j in range(len(sequence)):
            start=warmup+j*t
            points.extend([(start+on*t,0),(start+on*t+.5,d.vdd),(start+off*t-.5,d.vdd),(start+off*t,0)])
    c.source('input','stage1',pwl(levels),'external_analog_input')
    c.source('bias','bias',pwl(bias),'external_mac_bias')
    for name,pts in [('write',write),('reset',reset),('evaluate',evals),('awake',awake)]:
        c.source(name,name,pwl(pts),'clock_input')
    rs,rsb=c.logic('resetctl','reset')
    enable,enableb=c.logic('awake_ctl','awake')
    if d.detector=='hybrid':
        if d.systematic_offset_mv:
            # Input-referred systematic offset only; no canonical-answer source.
            c.lines.append(f'Ecmpoffset cmpsignal stage1 0 0 0')
            c.lines[-1]=f'Vcmpoffset cmpsignal stage1 {d.systematic_offset_mv*1e-3:.12g}'
            vin='cmpsignal'
        else:vin='stage1'
        ratio=d.input_scale/4
        if ratio>1:
            c.resistor('thresholdbottom','threshold','tap6',1e7)
            c.resistor('thresholdtop','tap7','threshold',(ratio-1)*1e7)
            reference='threshold'
        else:reference='tap7'
        # Matching-driven input upsizing also increases reference kickback.
        # Scale decoupling with input gate area instead of handicapping the
        # larger comparator with the small comparator's reference capacitor.
        actual_threshold_cap=max(threshold_cap_ff,threshold_cap_ff*d.comparator_w*d.comparator_l/.15)
        if actual_threshold_cap:c.capacitor(reference,actual_threshold_cap,'mim')
        positive,negative=c.comparator(vin,reference,'evaluate')
        mid=f'if(V({positive})>{d.vdd*.5:.12g},1,0)'
    else:mid=f'if(V(stage1)>{cm+.5*input_unit:.12g},1,0)'
    code='-3'+''.join('+'+(mid if th==.5 else f'if(V(stage1)>{cm+th*input_unit:.12g},1,0)') for th in (-2.5,-1.5,-.5,.5,1.5,2.5))
    c.lines.append(f'Bregion region 0 V={code}')
    c.save.add('v(region)')
    chosen='chosen' if d.write_buffer else 'stored'
    for digit in range(-3,4):
        name='sel'+str(digit+3);c.lines.append(f'B{name} {name} 0 V=if(abs(V(region)-({digit}))<.25,V(write),0)')
        c.save.update((f'v({name})',f'i(b{name})'))
        c.ports[f'i(b{name})']={'node':name,'component':'ideal_decision_pin'}
        on,off=c.logic(name+'ctl',name)
        c.tg(name,refs[digit],chosen,on,off,'dd_interface')
    if d.write_buffer:
        c.lines.append('Cchosen chosen 0 2f')
        c.amplifier('writer',chosen,'writerout','writer',None)
        # Close the storage aperture before the mux opens: avoid amplifying
        # tiny input-mux charge injection into the already-held state.
        wout=[(0,d.vdd),(warmup-3,0)]
        for j in range(len(sequence)):
            s=warmup+j*t
            wout += [(s+.22*t,0),(s+.22*t+.5,d.vdd),(s+.60*t-.5,d.vdd),(s+.60*t,0)]
        c.source('writeout','writeout',pwl(wout),'clock_input')
        wo,wob=c.logic('writectl','writeout');c.tg('writeout','writerout','stored',wo,wob,'dd_interface')
    c.tg('reset','tap6','stored',rs,rsb,'dd_interface')
    c.capacitor('stored')
    if direct and d.input_scale>4:
        c.resistor('feedbacktop','driven','feedback',(d.input_scale/4-1)*1e6)
        c.resistor('feedbackbottom','feedback','tap6',1e6)
        c.amplifier('driver','stage1','driven','buffer',enableb if d.gated_buffer else None,'feedback')
    else:c.amplifier('driver','stage1' if direct else 'stored','driven','buffer',enableb if d.gated_buffer else None)
    # Actual passive weighted sum. Bias is an external operand, not an oracle.
    inputnode='driven'
    c.resistor('loadx',inputnode,'mac',d.load_r/2)
    c.resistor('loadbias','bias','mac',d.load_r)
    c.lines.append(f'Cmac mac 0 {d.load_cap_ff:.12g}f');c.save.add('v(mac)')
    stop=warmup+len(sequence)*t
    meta={'design':asdict(d),'seed':d.seed,'kind':'loaded_regenerator','sequence':list(sequence),
        'residual_units_in_stage1':residual,'mode':'direct' if direct else 'regenerated',
        'unit_v':u,'input_unit_v':input_unit,'common_mode_v':cm,'mac_unit_v':u/3,'checkpoints':checks,
        'warmup_ns':warmup,'stop_ns':stop,'ports':c.ports,'mos':c.mos,
        'scope':'actual reference/TG/MIM-or-C/OTA/inverters; ideal-or-hybrid finite classifier; ideal sequencer and external operands',
        'unit_dimension':'micrometer','mim_ideal_total_ff':d.cap_ff}
    if d.detector=='hybrid':
        meta['threshold_reference_mim_ff']=actual_threshold_cap
        meta['threshold_reference_units']=math.ceil(actual_threshold_cap/1000)
    return c.finish(f'.tran {d.maxstep_ns:.12g}n {stop:.12g}n'),meta


def evaluate_boundary(result,names,x):
    meta=result['meta'];d=meta['design'];unit=meta['unit_v'];cm=meta['common_mode_v'];ix={n:i for i,n in enumerate(names)}
    ts=x[:,0];checks=[]
    def voltage(n):return np.zeros(len(x)) if n=='0' else x[:,ix['v('+n.lower()+')']]
    currents={name:x[:,ix[name]] for name in meta['ports']}
    for p in meta['checkpoints']:
        point=np.searchsorted(ts,p['time_ns']*1e-9)
        start=p['start_ns']*1e-9;stop=(p['start_ns']+d['period_ns'])*1e-9
        mask=(ts>=start)&(ts<=stop);energy={};net_energy={};returned={}
        for name,port in meta['ports'].items():
            current=currents[name];component=port['component'];power=-voltage(port['node'])*current
            if component in ('ideal_decision_pin','clock_input'):
                # VDD*positive charging charge, not ideal ramp energy recovery.
                value=d['vdd']*float(np.trapezoid(np.maximum(-current[mask],0),ts[mask]))
                component='control_pin_charging_upper'
            else:value=float(np.trapezoid(np.maximum(power[mask],0),ts[mask]))
            energy[component]=energy.get(component,0)+value*1e15
            net_energy[port['component']]=net_energy.get(port['component'],0)+float(np.trapezoid(power[mask],ts[mask]))*1e15
            returned[port['component']]=returned.get(port['component'],0)+float(np.trapezoid(np.maximum(-power[mask],0),ts[mask]))*1e15
        driver=voltage('driven')[point];storage=voltage('stored')[point];mac=voltage('mac')[point]
        mac_class=classify(mac-cm,meta['mac_unit_v'],3)
        target=cm+p['digit']*unit
        j=np.where((ts>=start)&(ts<=p['time_ns']*1e-9))[0]
        outside=j[abs(voltage('driven')[j]-target)>.10*unit]
        settled=None if len(outside) and outside[-1]==j[-1] else (float((ts[outside[-1]+1]-start)*1e9) if len(outside) else 0)
        known=sum(value for key,value in energy.items() if key not in ('external_analog_input','external_mac_bias'))
        checks.append(dict(p,stored_v=float(storage),driven_v=float(driver),mac_v=float(mac),
            selected_class=float(voltage('region')[point]),stored_error_units=float((storage-target)/unit),
            driver_error_units=float((driver-target)/unit),mac_class=mac_class,next_mac_correct=mac_class==p['digit'],
            energy_components_fj=energy,signed_port_work_fj=net_energy,returned_port_energy_fj=returned,
            measured_known_path_fj=known,settling_ns_to_point1digit=settled,
            current_class=classify(driver-cm,unit,3),retention_drift_units=float((voltage('stored')[point]-voltage('stored')[np.searchsorted(ts,p['write_end_ns']*1e-9)])/unit)))
        checks[-1]['reference_levels_v']={str(k):float(voltage('tap'+str(2*(k+3)))[point]) for k in range(-3,4)}
        wm=(ts>start+.18*d['period_ns']*1e-9)&(ts<start+.65*d['period_ns']*1e-9)
        checks[-1]['peak_reference_droop_units']=max(float(np.max(abs(voltage('tap'+str(2*(k+3)))[wm]-(cm+k*unit))))/unit for k in range(-3,4))
    # Explicit external MOS terminals, excluding the pre-operating-point row.
    limits=[]
    mask=ts>meta['warmup_ns']*1e-9
    for m in meta['mos']:
        vd,vg,vs,vb=(voltage(m[k])[mask] for k in ('d','g','s','b'))
        maxabs=max(float(np.max(abs(a-b))) for a,b in ((vd,vg),(vs,vg),(vd,vs)))
        vbs=vb-vs
        forward=float(np.max(vbs)) if m['type']=='n' else float(-np.min(vbs))
        if maxabs>1.95+1e-6 or forward>(.3 if m['type']=='n' else .1)+1e-6:
            limits.append({'instance':m['name'],'max_abs_vgs_vgd_vds':maxabs,'forward_body_bias_v':forward})
    result.update(observations=checks,device_operating_limit_violations=limits,
        all_next_mac_correct=all(c['next_mac_correct'] for c in checks),
        mean_energy_fj=float(np.mean([p['measured_known_path_fj'] for p in checks])),
        max_abs_driver_error_units=max(abs(p['driver_error_units']) for p in checks),
        energy_scope='VDD branch work plus conservative external control-pin charging; remaining ideal detector/latch/sequencer energy NOT measured; external operand work separate')
    return result


if __name__=='__main__':run_smokes()
