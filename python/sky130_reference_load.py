"""Small passive reference-bank fanout audit, not a CIM array.

Coherently switched writer loads are lumped: TG width and MIM load scale
together. A10fF/input gate-load allowance supplements25fF storage. This is a
reference loading test; no claim of64 complete buffered regeneration cells.
"""
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace,asdict
import json,math
import numpy as np
from gamma_sky130 import Circuit,Design,execute,compress_trace,pwl
from gamma_spice import classify
from gamma_adc_baseline import ROOT,OUT,dump,write_csv


def circuit(fanout,resistance):
    d=Design(cap_ff=35*fanout,cap_kind='mim',tg_w=.42*fanout,ladder_r=resistance)
    c=Circuit(d);refs=c.ladder();c.supply('interface')
    sequence=(-3,3,-3,0);period=400;warmup=1000;write=[(0,1.8),(warmup-1,0)];reset=[(0,0)]
    selectors={digit:[(0,1.8 if digit==-3 else 0),(warmup-1,0)] for digit in range(-3,4)}
    checks=[]
    for j,digit in enumerate(sequence):
        s=warmup+j*period
        write.extend([(s+72,0),(s+72.5,1.8),(s+259.5,1.8),(s+260,0)])
        reset.extend([(s+16,0),(s+16.5,1.8),(s+39.5,1.8),(s+40,0)])
        for k,points in selectors.items():
            points.extend([(s+72,0),(s+72.5,1.8 if k==digit else 0),(s+259.5,1.8 if k==digit else 0),(s+260,0)])
        checks.append({'start_ns':s,'time_ns':s+360,'digit':digit})
    c.source('reset','reset',pwl(reset),'clock_input');rq,rb=c.logic('resetctl','reset')
    c.tg('reset','tap6','stored',rq,rb,'dd_interface')
    for digit,points in selectors.items():
        name='select'+str(digit+3);c.source(name,name,pwl(points),'ideal_decision_pin');q,b=c.logic(name+'ctl',name)
        c.tg(name,refs[digit],'stored',q,b,'dd_interface')
    c.capacitor('stored');c.save.add('v(stored)')
    meta={'seed':101,'design':asdict(d),'mos':c.mos,'ports':c.ports,'checkpoints':checks,'fanout':fanout,
        'storage_plus_gate_allowance_per_load_ff':35,'stop_ns':2600,'reference_bank_resistance':resistance}
    return c.finish('.tran .25n 2600n'),meta


def run(item):
    fanout,resistance=item;directory=OUT/'spice'/f'reference_f{fanout}_r{int(resistance)}'
    r,names,x=execute(*circuit(fanout,resistance),directory,timeout=180)
    assert x[-1,0]>=2599.99e-9
    ix={n:i for i,n in enumerate(names)};t=x[:,0];v=x[:,ix['v(stored)']];observations=[]
    for p in r['meta']['checkpoints']:
        i=np.searchsorted(t,p['time_ns']*1e-9);m=(t>=p['start_ns']*1e-9)&(t<(p['start_ns']+400)*1e-9)
        unit=.025;error=(v[i]-(.6+p['digit']*unit))/unit
        energy={}
        for name,port in r['meta']['ports'].items():
            current=x[:,ix[name]];voltage=x[:,ix['v('+port['node']+')']]
            component=port['component']
            if component in ('clock_input','ideal_decision_pin'):
                value=1.8*np.trapezoid(np.maximum(-current[m],0),t[m])*1e15;component='control_pin_charging'
            else:value=np.trapezoid(np.maximum(-voltage[m]*current[m],0),t[m])*1e15
            energy[component]=energy.get(component,0)+float(value)
        wm=(t>(p['start_ns']+72)*1e-9)&(t<(p['start_ns']+260)*1e-9)
        droop=max(float(np.max(abs(x[wm,ix['v(tap'+str(2*(k+3))+')']]-(.6+k*.025))))/.025 for k in range(-3,4))
        observations.append({**p,'stored_v':float(v[i]),'error_digit_units':float(error),'decoded_class':classify(float(v[i]-.6),.025,3),
            'correct_class':classify(float(v[i]-.6),.025,3)==p['digit'],'peak_ref_droop_digits':droop,'energy_fj':energy})
    row={'fanout_equivalent_loads':fanout,'bank_resistance_ohm':resistance,'lumped_mim_load_ff':35*fanout,
        'peak_abs_write_error_digits':max(abs(p['error_digit_units']) for p in observations),
        'correct_classes':sum(p['correct_class'] for p in observations),'checkpoints':len(observations),
        'reference_fj_per_load_mean':float(np.mean([p['energy_fj']['reference']/fanout for p in observations])),
        'interface_plus_pin_fj_per_load_mean':float(np.mean([(p['energy_fj']['interface']+p['energy_fj']['control_pin_charging'])/fanout for p in observations])),
        'max_reference_droop_digits':max(p['peak_ref_droop_digits'] for p in observations),
        'run':str(directory.relative_to(ROOT)),
        'scope':'Lumped coherently switched load, scaled TG width and MIM capacitance. No output buffers, mismatch/fanout distribution or complete64-cell implementation.'}
    dump(directory/'metrics.json',{'summary':row,'observations':observations,'seed':101});compress_trace(directory/'output.raw')
    return row


def main():
    assert (OUT/'baseline_definition.json').exists()
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(run,[(1,4e6),(8,4e6),(64,4e6),(64,4e6/64)]))
    write_csv(OUT/'reference_fanout.csv',rows);dump(OUT/'reference_fanout.json',rows)
    print(json.dumps(rows,indent=2),flush=True)


if __name__=='__main__':main()
