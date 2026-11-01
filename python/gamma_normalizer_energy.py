"""Parameterized break-even accounting, not measured transistor energy."""
from __future__ import annotations
from collections import defaultdict
import math
import statistics
from pathlib import Path
import numpy as np
from gamma_protocols import bands


def summarize(workloads):
    from gamma_native_model import convolve,normalize
    layouts={}
    for q in (3,4,5):
        words=[normalize({0:n},q)[0] for n in range(-128,128)]
        positions=sorted({k for word in words for k in word})
        envelope=convolve({k:q-1 for k in positions},{k:q-1 for k in positions})
        final,_=normalize(envelope,q)
        lo=min([*envelope,*final]);hi=max([*envelope,*final])
        layouts[q]=(lo,hi)
    groups=defaultdict(list)
    for row in workloads:
        groups[(row['q'],row['scheme'])].append(row)
    out=[]
    for (q,scheme),group in groups.items():
        keys=['digit_products','width','rounds','trigger_events','unit_rewrites',
              'peak_units','bank_observations_until_settled_plus_confirmation',
              'comparator_evaluations_until_settled_plus_confirmation','final_readout_comparators']
        summary={'q':q,'scheme':scheme,'cases':len(group),
            'persistent_states':2*q-1}
        for key in keys:
            vals=[r[key] for r in group if r[key] is not None]
            summary[key+'_mean']=statistics.mean(vals) if vals else None
            summary[key+'_median']=statistics.median(vals) if vals else None
        # Detector activity cannot shrink to the observed result's span.
        # Freeze a physical layout from full INT8 encoder positions and a
        # positive absolute-coefficient envelope, independent of observations.
        lo,hi=layouts[q];width=hi-lo+1
        summary['observed_minimum_span_mean']=summary['width_mean']
        summary['compiled_layout']=[lo,hi]
        summary['compiled_layout_width']=width
        summary['width_mean']=summary['width_median']=width
        summary['final_readout_comparators_mean']=summary['final_readout_comparators_median']=width*2*(q-1)
        for kind in ('mean','median'):
            banks=width*(summary['rounds_'+kind]+1)
            summary['bank_observations_until_settled_plus_confirmation_'+kind]=banks
            summary['comparator_evaluations_until_settled_plus_confirmation_'+kind]=banks*2*len(bands(scheme)) if scheme!='bulk' else None
        summary['input_span_median']=statistics.median(v for r in group for v in r['input_spans'])
        out.append(summary)
    return out


def break_even(summary,unit_rewrite_fj=1,detect_comparator_fj=0.1,
               control_per_trigger_fj=0.5,delta_cim_fj=0,routing_per_unit_fj=0,
               scale_sidecar_fj=0,binary_accum_fj=20,adc_conversions=1,
               schedule='until_settled',phase_budget=32):
    if summary['scheme']=='bulk':
        raise ValueError('ADC-like unrestricted quotient excluded from circuit energy claim')
    q=summary['q']
    if schedule=='until_settled':
        detect=summary['comparator_evaluations_until_settled_plus_confirmation_mean']
    elif schedule=='fixed_budget':
        from gamma_protocols import bands
        detect=phase_budget*summary['width_mean']*2*len(bands(summary['scheme']))
    else:
        raise ValueError('unknown scheduling model')
    detect += summary['final_readout_comparators_mean']
    # q-dependence is explicit: an assumed equal per-unit charge-packet cost
    # gives q center packets plus two neighbor packets. Unit_rewrite_fj is
    # referenced to q=3 (five packets). This is an assumption to sweep.
    rewrite_cost=summary['unit_rewrites_mean']*unit_rewrite_fj*(q+2)/5
    detector_cost=detect*detect_comparator_fj
    control=summary['trigger_events_mean']*control_per_trigger_fj
    routing=2*summary['unit_rewrites_mean']*routing_per_unit_fj
    extra=delta_cim_fj+rewrite_cost+detector_cost+control+routing+scale_sidecar_fj
    return {'gamma_extra_fj':extra,'rewrite_fj':rewrite_cost,'detector_fj':detector_cost,
        'control_fj':control,'routing_fj':routing,'delta_cim_fj':delta_cim_fj,
        'scale_sidecar_fj':scale_sidecar_fj,'detect_comparator_evaluations':detect,
        'adc_break_even_fj_per_conversion':max(0,(extra-binary_accum_fj)/adc_conversions),
        'adc_conversions':adc_conversions,'binary_accum_fj':binary_accum_fj,
        'schedule':schedule,'phase_budget':phase_budget if schedule=='fixed_budget' else None,
        'assumptions':'Common E_CIM cancels only when delta_cim_fj=0. Until-settled scheduling is a lower-bound scenario requiring unpriced completion logic. Fixed-budget matches the demonstrated clocked architecture.'}


def energy_study(workloads,destination):
    destination=Path(destination); destination.mkdir(parents=True,exist_ok=True)
    summaries=summarize(workloads)
    rows=[]
    for summary in summaries:
        if summary['scheme']=='bulk':
            continue
        for x in (0.1,1,10,100):
            for y in (0.01,0.1,1,10):
                for adc in (1/64,1/16,1,4,16):
                    for schedule in ('until_settled','fixed_budget'):
                        for delta in (0,100):
                            r=break_even(summary,x,y,delta_cim_fj=delta,adc_conversions=adc,schedule=schedule)
                            rows.append({'q':summary['q'],'scheme':summary['scheme'],
                                'unit_rewrite_q3_reference_fj':x,'comparator_decision_fj':y,**r})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),sharey=True)
    xs=np.logspace(-1,2,100)
    for ax,schedule in zip(axes,['until_settled','fixed_budget']):
        for q,color in [(3,'#1b6ca8'),(4,'#d47a12'),(5,'#208652')]:
            s=next(v for v in summaries if v['q']==q and v['scheme']=='unit')
            for y,style in [(0.1,'-'),(1,'--')]:
                ys=[break_even(s,x,y,schedule=schedule)['adc_break_even_fj_per_conversion'] for x in xs]
                ax.plot(xs,ys,color=color,ls=style,label=f'q={q}, detection={y:g} fJ')
        ax.set_xscale('log');ax.set_yscale('symlog',linthresh=1)
        ax.set_xlabel('Unit rewrite cost at q=3 (fJ; swept assumption)')
        ax.set_title('Completion-gated lower bound' if schedule=='until_settled' else '32-phase fixed budget')
        ax.grid(True,alpha=.25)
    axes[0].set_ylabel('ADC energy needed to break even (fJ/conversion)')
    axes[1].legend(fontsize=7,loc='upper left')
    fig.suptitle('One conversion/output; binary accumulation 20 fJ; common array cost assumed equal')
    fig.tight_layout();fig.savefig(destination/'energy_break_even.png',dpi=170);plt.close(fig)
    return {'workload_summaries':summaries,'parameter_grid':rows,
        'workload_scope':'128 paired scalar integer products per q/protocol, not an amortized neural dot. Detector counts use a fixed layout compiled from the complete input encoding positions and a positive absolute-coefficient envelope, never observed result span. ADC counts 1/64 and 1/16 model per-term baseline amortization. Gamma counters still normalize each product; other schedules require new counters.',
        'equation':'N_ADC*E_ADC + E_binary_accum > delta_E_CIM + N_comparator*E_detect + N_unit_rewrite*E_unit(q) + N_trigger*E_control + E_routing + E_scale_sidecar',
        'energy_note':'All fJ quantities are swept assumptions, not measured component energies. Behavioral ideal current sources do not measure supply energy.',
        'additional_terms':'Routing and scale/sidecar terms are exposed parameters; plotted slices set them to zero. No zero-cost physical conclusion follows.'}


def recommendation(energy,spice_records,scaling):
    from gamma_native_model import normalize
    from gamma_protocols import bands
    candidates=[]
    for q in (3,4,5):
        s=next(v for v in energy['workload_summaries'] if v['q']==q and v['scheme']=='unit')
        enc=[normalize({0:n},q)[0] for n in range(-128,128)]
        # Weight-only, complete-input-layout worst-case channel envelope.
        scales=[max(1,(q-1)*sum(map(abs,w.values()))*0.1/0.6) for w in enc]
        margin=statistics.median(0.05/v for v in scales)
        ideal=[r for r in spice_records if r['q']==q and r.get('profile')=='ideal']
        other=[r for r in spice_records if r['q']==q and r.get('profile')!='ideal']
        e=break_even(s,schedule='fixed_budget')
        candidates.append({'q':q,'span':s['input_span_median'],
            'digit_products':s['digit_products_mean'],'rewrite_energy_parameterized_fj':e['rewrite_fj'],
            'detector_energy_parameterized_fj':e['detector_fj'],'persistent_states':2*q-1,
            'transient_peak_mean_units':s['peak_units_mean'],'median_compiled_half_margin_v':margin,
            'ideal_spice_failures':sum(r['network_class_failure'] for r in ideal),
            'stress_failure_fraction':statistics.mean(r['network_class_failure'] for r in other) if other else 0,
            'fixed_budget_reference_extra_fj':e['gamma_extra_fj']})
    metrics={'span':.10,'digit_products':.10,'rewrite_energy_parameterized_fj':.15,
        'detector_energy_parameterized_fj':.15,'persistent_states':.10,
        'transient_peak_mean_units':.10,'inverse_margin':.10,'stress_failure_fraction':.10}
    for c in candidates:
        c['inverse_margin']=1/c['median_compiled_half_margin_v']
    for c in candidates:
        c['combined_score']=sum(weight*c[key]/max(max(x[key] for x in candidates),1e-30) for key,weight in metrics.items())
    winner=min(candidates,key=lambda c:c['combined_score'])
    winners={q:0 for q in (3,4,5)}
    # Energy-only sensitivity separately, never disguise the chosen weights.
    for x in (0.1,1,10,100):
        for y in (0.01,.1,1,10):
            for schedule in ('until_settled','fixed_budget'):
                scores={q:break_even(next(v for v in energy['workload_summaries'] if v['q']==q and v['scheme']=='unit'),x,y,schedule=schedule)['gamma_extra_fj'] for q in winners}
                winners[min(scores,key=scores.get)]+=1
    return {'recommended_next_prototype_q':winner['q'],'candidates':candidates,'score_weights':metrics,
        'reference_energy_parameters':'q3 unit rewrite 1 fJ, comparator decision 0.1 fJ, trigger control 0.5 fJ; fixed 32 phases; assumed, not measured',
        'energy_only_sensitivity_winner_counts':winners,
        'qualification':'Provisional behavioral-prototype choice under disclosed weights and parameter assumptions, not an optimized chip choice. Stress profiles have arbitrary severities, not manufacturing distributions. No cost-below-ADC conclusion without matched device energies and timing.'}
