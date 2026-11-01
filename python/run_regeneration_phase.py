"""Reproducible q=4 regeneration phase; actual LTspice and labeled control models."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, statistics, time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from dataclasses import asdict, replace
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from gamma_regen_spice import *
from gamma_regen_golden import event_completion,finite_word,finite_merge
from gamma_native_model import normalize,convolve,NativeFloat,value_pair
from gamma_protocols import settle
from gamma_normalizer_statistics import confidence

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'spice'/'regen_results'
RUNS=ROOT/'spice'/'regen_runs'


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2),encoding='ascii')


def jobs(n):
    result=[]
    def add(group,fn,kwargs,**labels):result.append({'group':group,'fn':fn,'kwargs':kwargs,'labels':labels})
    for mapping in ('single','differential','charge'):
        for d in range(-3,4):
            add('primitive_ideal','primitive',dict(digit=d,mapping=mapping),mapping=mapping,digit=d)
    for mode in ('none','regen'):
        for slope in (F(2),F(1,2)):
            for d in range(-3,4):
                add('chain_ideal','chain',dict(digit=d,slope=slope,mode=mode),mode=mode,slope=str(slope),digit=d)
            add('chain_drift','chain',dict(digit=1,slope=slope,mode=mode,physical=RegenPhysical(stage_drift_units=.06)),
                mode=mode,slope=str(slope),digit=1)
        for a,b in [(1,2),(2,1),(2,2),(3,1),(3,2),(-1,-2),(-2,-1),(-3,-2)]:
            for drift in (0,.35):
                add('word','word',dict(a=a,b=b,bias=1 if a>0 else -1,mode=mode,
                    physical=RegenPhysical(stage_drift_units=drift)),mode=mode,drift=drift,a=a,b=b)
    add('boundary','boundary',{})
    for pname in ('low','stress'):
        for mapping in ('single','differential','charge'):
            for d in range(-3,4):
                for trial in range(max(4,n//4)):
                    seed=20000+trial*19+d+3
                    add('primitive_mc','primitive',dict(digit=d,mapping=mapping,physical=profile(pname),seed=seed),
                        profile=pname,mapping=mapping,digit=d,trial=trial)
        for slope in (F(2),F(1,2)):
            for trial in range(n):
                d=__import__('random').Random(70000+trial).randrange(-3,4)
                for mode in ('none','regen'):
                    add('chain_mc','chain',dict(digit=d,slope=slope,mode=mode,physical=profile(pname),seed=50000+trial),
                        profile=pname,slope=str(slope),mode=mode,digit=d,trial=trial)
    # Isolated mechanisms are sensitivity brackets, not calibrated process corners.
    mechanisms={
        'reference':('reference_sigma_units',(.02,.25)),
        'switch':('ron_ohm',(1e3,2e4)),
        'capacitance':('cap_sigma',(.02,.15)),
        'feedthrough':('clock_feedthrough_f',(.05e-15,2e-15)),
        'charge_injection':('charge_injection_sigma_c',(2e-18,5e-16)),
        'comparator':('comparator_sigma_units',(.02,.35)),
        'supply':('supply_gain',(.01,.15)),
        'thermal':('thermal_kelvin',(300,30000)),
    }
    for mechanism,(key,levels) in mechanisms.items():
        for level,value in enumerate(levels):
            for d in (-3,0,3):
                for trial in range(2):
                    add('sensitivity','primitive',dict(digit=d,mapping='differential',
                        physical=replace(RegenPhysical(),**{key:value}),seed=90000+trial*7+d+3),
                        mechanism=mechanism,level=level,value=value,digit=d,trial=trial)
    for scale in (4,8,16,32):
        for trial in range(max(4,n//4)):
            # Keep physical reference/comparator errors in volts fixed across S.
            p=profile('low');factor=scale/4
            for key in ('reference_sigma_units','comparator_sigma_units','threshold_bank_sigma_units','mac_noise_units'):
                setattr(p,key,getattr(p,key)*factor)
            p.scale=scale
            add('scale','chain',dict(digit=trial%7-3,physical=p,seed=100000+trial),scale=scale,trial=trial)
    return result


def execute(item,index):
    directory=RUNS/f'{"v" if index<86 else "w"}{index:04d}'
    fn={'primitive':primitive_deck,'chain':chain_deck,'word':word_deck,'boundary':boundary_deck}[item['fn']]
    deck,meta=fn(**item['kwargs'])
    fingerprint=hashlib.sha256((deck+json.dumps(meta,sort_keys=True)+'verifier_v1').encode()).hexdigest()
    cache=directory/'cached.json'
    if cache.exists():
        row=json.loads(cache.read_text())
        if row['fingerprint']==fingerprint:return row
    row=run(deck,meta,directory,retain_raw=(index<5 or item['group'] in ('boundary','chain_drift') or (item['group']=='word' and index%8==0)))
    row.update(group=item['group'],labels=item['labels'],run_path=str(directory.relative_to(ROOT)),fingerprint=fingerprint)
    dump(cache,row)
    return row


def completion_study():
    source=json.loads((ROOT/'spice/results/protocol_workloads.json').read_text())['records']
    source=[r for r in source if r['q']==4 and r['scheme']=='unit']
    output=[]
    for item in source:
        a,b=item['inputs'];da,_=normalize({0:a},4);db,_=normalize({0:b},4);raw=convolve(da,db)
        sync=settle(raw,4,'unit',active=(-8,8));event=event_completion(raw,range(-8,9))
        assert value_pair(event['state'],4)==value_pair(raw,4)
        output.append({'inputs':[a,b],'raw':raw,'synchronous':sync,'event':event,
            'fixed32_comparator_evaluations':32*17*2+17*8,
            'ideal_until_settled_comparator_evaluations':(sync['rounds']+1)*17*2+17*8,
            'explicit_event_comparator_evaluations':event['counts']['comparator_evaluations']+17*8,
            'fixed32_time_ns':128,'ideal_until_settled_time_ns':4*(sync['rounds']+1),
            'event_time_ns':event['counts']['completion_time_ns'],
            'fixed32_budget_exceeded':sync['rounds']>32,
            'event_equivalent_different_word':event['state']!=sync['state'],
            'event_final_exact_pair':[str(v) for v in value_pair(event['state'],4)],
            'digit_products':item['digit_products']})
    dump(OUT/'completion_traces.json',output)
    keys=['fixed32_comparator_evaluations','ideal_until_settled_comparator_evaluations',
        'explicit_event_comparator_evaluations','fixed32_time_ns','ideal_until_settled_time_ns','event_time_ns','digit_products']
    report={'cases':len(output),**{k:statistics.mean(r[k] for r in output) for k in keys},
        'unit_rewrites_mean':statistics.mean(r['event']['counts'].get('unit_rewrites',0) for r in output),
        'synchronous_rewrites_mean':statistics.mean(r['synchronous']['unit_rewrites'] for r in output),
        'fixed32_budget_exceeded':sum(r['fixed32_budget_exceeded'] for r in output),
        'event_equivalent_different_word':sum(r['event_equivalent_different_word'] for r in output),
        'control_model':'heap event queue + per-request revalidation + charged request/ack/tree; no asynchronous circuit claim'}
    counts=set().union(*(r['event']['counts'] for r in output))
    report['event_mean_counts']={k:statistics.mean(r['event']['counts'].get(k,0) for r in output) for k in counts}
    return report


def aggregate(rows):
    groups=defaultdict(list)
    for row in rows:
        excluded=('trial','a','b') if row['kind']=='primitive' else ('trial','digit','a','b')
        labels={k:v for k,v in row['labels'].items() if k not in excluded}
        key=json.dumps([row['group'],labels],sort_keys=True)
        groups[key].append(row)
    out=[]
    for key,items in groups.items():
        group,labels=json.loads(key)
        depths=(2,4,8,16) if items[0]['kind']=='chain' else (None,)
        for depth in depths:
            obs=[r['observations'][depth-1] if depth else r['observations'][0] for r in items]
            if items[0]['kind']=='word':
                wrong=sum(not r['final_algebraically_correct'] for r in items)
                strict_wrong=sum(any(not o['correct'] for o in r['observations'] if o.get('stage')==2) for r in items)
                equivalent=sum(r['equivalent_different_word'] for r in items)
                unresolved=sum(any(isinstance(o['decoded_class'],str) for o in r['observations'] if o.get('stage')==2) for r in items)
            else:
                wrong=sum(not o['correct'] for o in obs);strict_wrong=wrong;equivalent=0
                unresolved=sum(o['overflow_latched'] or isinstance(o['decoded_class'],str) for o in obs)
            residuals=[abs(o['residual_units']) for o in obs]
            c=confidence(wrong,len(items))
            out.append({'group':group,**labels,'depth':depth,'circuits':len(items),'wrong_final':wrong,
                'strict_word_wrong':strict_wrong,'equivalent_different_word':equivalent,'unresolved_overflow':unresolved,
                'p_correct':1-wrong/len(items),'error_confidence':c,
                'abs_residual_units_median':statistics.median(residuals),'abs_residual_units_max':max(residuals),
                'reference_port_delivered_fj_mean':statistics.mean(r['reference_port_delivered_energy_fj'] for r in items),
                'settling_ns_mean':statistics.mean(r['first_output_settling_ns_to_1pct_class_spacing'] for r in items)})
    return out


def energy_study(completion):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n=17; ec=completion['event_mean_counts']
    ctrl=sum(ec.get(k,0) for k in ('queue_insertions','queue_pops','requests','acknowledgments',
        'activity_tree_edge_toggles','completion_tree_initialization_edges','completion_guard_evaluations','completion_ack_root_events'))
    schedules={
        'fixed32':(completion['fixed32_comparator_evaluations'],completion['synchronous_rewrites_mean'],32*n),
        'ideal completion lower bound':(completion['ideal_until_settled_comparator_evaluations'],completion['synchronous_rewrites_mean'],0),
        'explicit event model':(completion['explicit_event_comparator_evaluations'],completion['unit_rewrites_mean'],ctrl)}
    records=[]
    for schedule,(detect,rewrite,control) in schedules.items():
        for er in (.1,1,10,100):
            for ed in (.01,.1,1,10):
                for eregen in (.01,.1,1,10,100):
                    for econtrol in (.001,.01,.1,1):
                        for conversions in (1/64,1/16,1,4,16):
                            base=detect*ed+rewrite*er*6/5+control*econtrol
                            # Eregen includes reference generation, reset/write, driver and local clocks.
                            # Routing/CIM delta and exponent metadata explicitly parameterized below.
                            for overhead in (0,100,1000):
                                records.append({'schedule':schedule,'rewrite_q3_unit_fj':er,'comparator_fj':ed,
                                'regen_per_position_fj':eregen,'control_event_fj':econtrol,'adc_conversions':conversions,
                                'gamma_extra_fj':base+n*eregen+overhead,'adc_break_even_fj':(base+n*eregen+overhead-20-5)/conversions,
                                'binary_accum_fj':20,'requantization_fj':5,'delta_cim_routing_sidecar_fj':overhead})
    dump(OUT/'energy_sweep.json',records)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    x=np.logspace(-2,2,120)
    reference={}
    for schedule,(detect,rewrite,control) in schedules.items():
        base=detect*.1+rewrite*1.2+control*.01
        axes[0].plot(x,base+17*x-25,label=schedule)
        reference[schedule]={'detect_fj':detect*.1,'rewrite_fj':rewrite*1.2,'completion_control_fj':control*.01,
            'before_regen_extra_fj':base,'adc_break_even_at_1fj_regen':base+17-25}
    axes[0].set(xscale='log',xlabel='All-in regeneration per physical position (fJ)',
        ylabel='ADC break-even energy (fJ / matched output)',title='One ADC conversion/output; binary+requant = 25 fJ')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    detect,rewrite,control=schedules['explicit event model'];base=detect*.1+rewrite*1.2+control*.01
    adc=np.logspace(0,4,120)
    for conv in (1/16,1,4):
        axes[1].plot(adc,(conv*adc+25-base)/17,label=f'{conv:g} ADC conversion(s)/output')
    axes[1].axhline(0,color='black',lw=.7)
    axes[1].set(xscale='log',xlabel='Energy per ADC conversion (fJ)',ylabel='Maximum all-in regeneration (fJ / position)',
        title='Explicit event model; negative budget means no win')
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
    fig.tight_layout();fig.savefig(OUT/'energy_break_even.png',dpi=170);plt.close(fig)
    gx,gy=np.meshgrid(np.linspace(0,15,151),np.linspace(0,250,151))
    minimum_adc=detect*.1+rewrite*1.2+17*gx+gy-25
    fig,ax=plt.subplots(figsize=(7,5))
    fill=ax.contourf(gx,gy,minimum_adc,levels=np.linspace(0,600,31),cmap='viridis',extend='both')
    contour=ax.contour(gx,gy,minimum_adc,levels=(50,100,200,500),colors='white')
    ax.clabel(contour,fmt='%g fJ ADC',fontsize=8)
    ax.plot([1],[control*.01],'rx',label='Illustrative parameter point')
    ax.set(xlabel='All-in regeneration (fJ / position)',ylabel='Completion + scheduling control (fJ / product)',
        title='ADC break-even: 17 positions, one conversion/output\nDelta CIM/routing/sidecar = 0; binary+requant = 25 fJ')
    fig.colorbar(fill,ax=ax,label='Minimum ADC energy for gamma to win (fJ)')
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'energy_regen_completion.png',dpi=170);plt.close(fig)
    return {'compiled_positions':17,'q':4,'reference_point':reference,'sweep_rows':len(records),
        'assumptions':'Costs parameterized, not measured. Comparator .1 fJ, q3-unit rewrite 1 fJ scaled6/5, control activity .01 fJ, binaryacc20, requant5. Delta CIM/routing/sidecar0 is optimistic. Eregen must include buffers, reference generation and local clocks. No supply-energy superiority established.',
        'equation':'ADC_min=(Ndetect*Edetect+Nrewrite*Erewrite+Ncontrol*Econtrol+17*Eregen+delta_CIM+routing+sidecar-Ebinary-Erequant)/Nadc'}


def plots(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10,8))
    for col,slope in enumerate(('2','1/2')):
        for pname in ('low','stress'):
            for mode in ('none','regen'):
                rows=[r for r in summary if r['group']=='chain_mc' and r['profile']==pname and r['mode']==mode and r['slope']==slope]
                rows.sort(key=lambda r:r['depth'])
                x=[r['depth'] for r in rows]
                axes[0,col].plot(x,[r['p_correct'] for r in rows],marker='o',ls='-' if mode=='regen' else '--',label=pname+' '+mode)
                axes[1,col].plot(x,[r['abs_residual_units_median'] for r in rows],marker='o',ls='-' if mode=='regen' else '--',label=pname+' '+mode)
        axes[0,col].plot([2,4,8,16],[1]*4,':k',label='exact digital gamma')
        axes[0,col].set(title=f'Fixed-point affine MAC, slope {slope}',ylabel='Correct final class probability',ylim=(-.02,1.02))
        axes[1,col].set(xlabel='Connected analog stages',ylabel='Median absolute residual (digit units)')
        axes[1,col].set_yscale('symlog',linthresh=.001)
        axes[0,col].legend(fontsize=8)
        for row in (0,1):axes[row,col].grid(alpha=.2)
    fig.suptitle('Actual LTspice chains; synthetic parameter brackets, not process-qualified')
    fig.tight_layout();fig.savefig(OUT/'chain_depth.png',dpi=170);plt.close(fig)
    with (OUT/'summary.csv').open('w',newline='',encoding='ascii') as f:
        keys=sorted(set().union(*(r.keys() for r in summary)))
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(summary)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--trials',type=int,default=16)
    parser.add_argument('--workers',type=int,default=1);args=parser.parse_args()
    if args.trials<4:parser.error('at least4 trials required')
    OUT.mkdir(parents=True,exist_ok=True);work=jobs(args.trials);rows=[None]*len(work);start=time.monotonic()
    dump(OUT/'manifest.json',{'q':4,'seed_policy':'seed stored in every golden.json','jobs':len(work),
        'trials':args.trials,'workers':args.workers,'source_hashes':{f:hashlib.sha256((ROOT/'python'/f).read_bytes()).hexdigest()
            for f in ('gamma_regen_spice.py','gamma_regen_golden.py','run_regeneration_phase.py')},
        'simulator':str(DEFAULT_EXE),'physical_profiles':{k:asdict(profile(k)) for k in ('ideal','low','stress')}})
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(execute,item,i):i for i,item in enumerate(work)}
        for count,future in enumerate(as_completed(futures),1):
            i=futures[future]
            try:rows[i]=future.result()
            except Exception as exc:
                dump(OUT/'incomplete_run.json',{'index':i,'group':work[i]['group'],'labels':work[i]['labels'],'error':str(exc)})
                for pending in futures:pending.cancel()
                pool.shutdown(wait=False,cancel_futures=True)
                raise
            if count%20==0 or count==len(work):print(f'{count}/{len(work)} actual SPICE jobs complete ({time.monotonic()-start:.1f}s)',flush=True)
    summary=aggregate(rows);completion=completion_study();energy=energy_study(completion);plots(summary)
    finite=[finite_word(w,exponent=e,policy=policy,renormalize_small=small)
        for w,e,small in [({3:1,-2:1},0,False),({-1:1},0,True),({-3:3,-4:2,0:1},0,False)]
        for policy in ('truncate','nearest')]
    finite.append(finite_merge(NativeFloat({0:2},4,q=4),NativeFloat({0:1},0,q=4)))
    dump(OUT/'finite_word_golden.json',finite)
    report={'q':4,'actual_spice_circuits':len(rows),'trial_count':args.trials,'summary':summary,
        'completion':completion,'energy':energy,'elapsed_seconds':time.monotonic()-start,
        'limitations':['Behavioral seven-state references and ideal buffers, not transistor hardware.',
            'Reference-port energy excludes comparator, clock, buffer and reference-generation supplies.',
            'Monte Carlo profiles are synthetic parameter brackets, not foundry statistics.',
            'Event-driven completion is a separately labeled explicit control model, not asynchronous SPICE.',
            'No 1e-3 or1e-6 SPICE reliability certification; intervals use independent circuits, not nodes.',
            'Depth affine workloads have exact fixed points; they do not establish neural-network accuracy.',
            'New word path tests ideal toppling plus controlled residual; prior phase covers toppling mismatch.']}
    dump(OUT/'phase_summary.json',report)
    incomplete=OUT/'incomplete_run.json'
    if incomplete.exists():
        previous=json.loads(incomplete.read_text())
        previous['status']='Earlier numerical configuration; excluded. All final configured jobs now complete.'
        dump(OUT/'resolved_numerical_diagnostic.json',previous)
        incomplete.unlink()
    print(json.dumps({'circuits':len(rows),'completion':completion,'elapsed_seconds':report['elapsed_seconds']},indent=2))


if __name__=='__main__':main()
