"""Independent model, batching, mismatch and switch/reference diagnostics."""
from dataclasses import replace,asdict
from pathlib import Path
import json
import numpy as np
from gamma_sky130 import *
from gamma_sky_batch import batch_execute


def pair_deck(d):
    c=Circuit(d);c.source('gate','gate',1,'external_input')
    for j in (1,2):
        c.source('drain'+str(j),'drain'+str(j),1,'drain')
        c.fet('n'+str(j),'n','drain'+str(j),'gate','0','0',1,.15)
    return c.finish('.dc Vgate .9 1.1 .01'),{'kind':'pair_diagnostic','seed':d.seed,'design':asdict(d)}


def switch_deck(d,on):
    c=Circuit(d);dd=c.supply('switch');c.source('input','input',.6,'input')
    c.source('control','control',d.vdd if on else 0,'control');c.inv('i','control','bar',dd)
    c.tg('pass','input','out','control','bar',dd);c.resistor('load','out','0',1e5)
    return c.finish('.dc Vinput .525 .675 .025'),{'kind':'switch','seed':d.seed,'design':asdict(d),'on':on}


def bank_deck(d):
    c=Circuit(d);c.ladder()
    return c.finish('.dc Vreference 1.62 1.9 .01'),{'kind':'reference','seed':d.seed,'design':asdict(d)}


def main():
    pairs=[Design(seed=17),Design(pdk_mismatch=True,seed=17),Design(pdk_mismatch=True,seed=18),Design(pdk_mismatch=True,seed=17)]
    cases=[(*pair_deck(d),BASE/'smoke'/f'mismatch_{i}') for i,d in enumerate(pairs)]
    for i,(w,on) in enumerate(( (w,on) for w in (.42,.84,1.68) for on in (False,True))):
        cases.append((*switch_deck(Design(tg_w=w),on),BASE/'regen_cell'/f'switch_probe{i}'))
    cases.append((*bank_deck(Design()),BASE/'reference_bank/dc'))
    cases.append((*pair_deck(Design(temp_c=100)),BASE/'smoke/temp_batch'))
    result=batch_execute(cases,workers=1)
    currents=[]
    for i in range(4):
        names,x=fast_raw(cases[i][2]/'output.raw');j=np.argmin(abs(x[:,names.index('v(gate)')]-1))
        currents.append([-float(x[j,names.index('i(vdrain'+str(k)+')')]) for k in (1,2)])
    report={'pair_currents_a':currents,'mismatch_disabled_pair_equal':abs(currents[0][0]-currents[0][1])<1e-13,
        'mismatch_enabled_instances_differ':abs(currents[1][0]-currents[1][1])>1e-10,
        'different_seeds_differ':not np.allclose(currents[1],currents[2],rtol=1e-8,atol=1e-14),
        'same_seed_replay_identical':np.allclose(currents[1],currents[3],rtol=1e-10,atol=1e-14)}
    switches=[]
    for _,meta,path in cases[4:10]:
        names,x=fast_raw(path/'output.raw')
        for row in x:
            vi=row[names.index('v(input)')];vo=row[names.index('v(out)')];current=-row[names.index('i(vinput)')]
            switches.append({'w_um':meta['design']['tg_w'],'on':meta['on'],'input_v':float(vi),
                'output_v':float(vo),'input_current_a':float(current),
                'ron_ohm':float((vi-vo)/current) if meta['on'] else None})
    report['selector_sweep']=switches
    names,x=fast_raw(BASE/'reference_bank/dc/output.raw')
    report['reference_bank_dc']=[{'vdd':float(row[names.index('v(dd_reference)')]),
        'supply_current_a':-float(row[names.index('i(vreference)')]),
        'levels_v':{str(d):float(row[names.index('v(tap'+str(2*(d+3))+')')]) for d in range(-3,4)}} for row in x]
    names,xb=fast_raw(BASE/'smoke/temp_batch/output.raw')
    # Standalone literal .temp, independent of batch's set-temp command.
    r,names2,xs=execute(*pair_deck(Design(temp_c=100)),BASE/'smoke/temp_standalone')
    report['batch_vs_literal_temperature_match']=bool(np.allclose(xb[:,names.index('i(vdrain1)')],xs[:,names2.index('i(vdrain1)')],rtol=1e-8,atol=1e-14))
    report['ngspice_random_scope']='Local device mismatch via PDK MC_MM_SWITCH; MC_PR_SWITCH=0. Not calibrated manufacturing yield.'
    dump(BASE/'results/model_validation.json',report)
    assert all(report[k] for k in ('mismatch_disabled_pair_equal','mismatch_enabled_instances_differ','different_seeds_differ','same_seed_replay_identical','batch_vs_literal_temperature_match'))
    print({k:v for k,v in report.items() if k not in ('selector_sweep','reference_bank_dc')},flush=True)


if __name__=='__main__':main()
