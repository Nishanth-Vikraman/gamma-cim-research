"""Only missing Stage-C smoke: actual SKY130 TG charge/discharge supply audit.

Runs AFTER baseline A/B files exist. Reuses the completed model environment;
does not regenerate prior SKY130 or LTspice reports/circuits.
"""
from dataclasses import asdict
import json,math
import numpy as np
from gamma_sky130 import Circuit,Design,execute,compress_trace
from gamma_adc_baseline import ROOT,OUT,dump,sha


def deck(kind):
    d=Design(cap_ff=100,cap_kind=kind,tg_w=.42)
    c=Circuit(d);charge=c.supply('charge');interface=c.supply('interface')
    c.source('enable','enable','PULSE(0 1.8 20n .5n .5n 40n 200n)','clock_input')
    c.source('discharge','discharge','PWL(0 1.8 10n 1.8 10.5n 0 80n 0 80.5n 1.8 120n 1.8 120.5n 0)','clock_input')
    c.inv('chargeinv','enable','enableb',interface);c.inv('dischargeinv','discharge','dischargeb',interface)
    c.tg('charge','dd_charge','stored','enable','enableb',interface)
    c.tg('discharge','stored','0','discharge','dischargeb',interface)
    c.capacitor('stored');c.resistor('leak','stored','0',1e15)
    c.save.add('v(stored)')
    return c.finish('.tran .05n 140n'),{'kind':'capacitor_supply_audit','seed':101,'design':asdict(d),'mos':c.mos,'ports':c.ports,'stop_ns':140}


def main():
    assert (OUT/'baseline_definition.json').exists(),'Establish ADC calibration and baseline definitions before new transistor work'
    old=json.loads((ROOT/'spice/sky130/results/smoke_summary.json').read_text());assert len(old)==20 and all(r['passed'] for r in old)
    rows=[]
    for kind in ('ideal','mim'):
        directory=OUT/'spice'/('cap_smoke_'+kind)
        r,names,x=execute(*deck(kind),directory,timeout=150)
        assert x[-1,0]>=139.99e-9
        ix={n:i for i,n in enumerate(names)};t=x[:,0];v=x[:,ix['v(stored)']]
        get=lambda ns:float(v[np.searchsorted(t,ns*1e-9)])
        mask=(t>=19e-9)&(t<=65e-9);power=-x[:,ix['i(vcharge)']]*x[:,ix['v(dd_charge)']]
        delivered=float(np.trapezoid(np.maximum(power[mask],0),t[mask]))*1e15
        expected=100*1.8*(get(59)-get(19))
        probes=lambda node:np.zeros(len(x)) if node=='0' else x[:,ix['v('+node.lower()+')']]
        peak=0.
        for m in r['meta']['mos']:
            vd,vg,vs,vb=(probes(m[k]) for k in ('d','g','s','b'))
            peak=max(peak,*[float(np.max(abs(a-b))) for a,b in ((vg,vs),(vg,vd),(vd,vs),(vg,vb))])
        result={'capacitor':kind,'corner':'tt','vdd':1.8,'temperature_c':25,'nominal_cap_ff':100,
            'charged_v':get(59),'held_v':get(75),'discharged_v':get(135),
            'charge_rail_delivered_fj':delivered,'nominal_C_VDD_deltaV_fj':expected,'supply_to_nominal_cap_ratio':delivered/expected,
            'storage_energy_after_charge_fj':.5*100*get(59)**2,'peak_monitored_terminal_differential_v':peak,
            'passed':get(59)>1.78 and abs(get(135))<.02 and .8<delivered/expected<1.3 and peak<=1.95,
            'accounting_note':'C*VDD*deltaV is supplied charge work;0.5*C*V^2 is stored energy. Difference includes switch parasitics; interface/clock rail work is separate, not claimed free.',
            'run':str(directory.relative_to(ROOT)),'circuit_sha256':sha(directory/'circuit.cir')}
        dump(directory/'metrics.json',result);compress_trace(directory/'output.raw');rows.append(result);print(result,flush=True)
    dump(OUT/'sky130_smoke_summary.json',{'retained_core_checks':old,'new_capacitor_checks':rows,
        'all_passed':all(r['passed'] for r in old+rows),'retained_source_sha256':sha(ROOT/'spice/sky130/results/smoke_summary.json'),
        'order':'Existing20 core/corner smokes already validated; only nominal capacitor supply audit was added, after Phase A/B calibration.'})
    assert all(r['passed'] for r in rows)


if __name__=='__main__':main()
