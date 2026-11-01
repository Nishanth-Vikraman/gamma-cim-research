"""Run the local normalization engineering phase (LTspice required).

Default: exact protocol study, 108 ideal circuit cases, physical sweeps, a
separately labeled fast scaling model, and parameterized energy/recommendation.
Every completed circuit job is cached by input/source hash for resumability.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import statistics
import time
import numpy as np

from gamma_protocols import (SCHEMES, default_cases, export_golden, golden_case,
    transient_study, workload_study)
from gamma_spice import DEFAULT_EXE, Physical, realization, run_spice
from gamma_normalizer_statistics import charge_monte_carlo, confidence, sampled_predict
from gamma_normalizer_energy import energy_study, recommendation

ROOT=Path(__file__).resolve().parents[1]


def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def settling_metadata(record,case):
    last_injection=max(event['time_ns'] for event in case['stimuli'])
    last=record['last_active_phase_time_ns']
    delta=max(0,(last or 0)-last_injection)
    record['operational_settling_time_after_last_injection_ns']=delta
    record['last_active_round_after_last_injection']=max(0,int(round((delta-3.85)/4))+1) if delta else 0
    record['settling_definition']='Last active transfer-phase observation within the fixed run budget; not guaranteed perpetual analog settling under ongoing noise.'
    if record['headroom_violation']:
        mode='outside_declared_physical_headroom'
    elif record['exact_class_value_failures']:
        mode='arithmetic_output_error_or_unresolved_overflow'
    elif record['network_class_failure']:
        mode='different_oracle_digits_but_equivalent_value'
    elif record['phase_class_failures']:
        mode='transient_path_changed_final_value_correct'
    elif record['boundary_requires_extension']:
        mode='correct_with_passive_boundary_receivers_active_window_needs_extension'
    elif record['max_abs_residual_units']>.25:
        mode='correct_class_reduced_margin'
    else:
        mode='correct'
    record['failure_mode']=mode
    return record


def run_job(job,exe):
    case,physical,seed,directory,profile,severity=job
    directory=Path(directory)
    digest=hashlib.sha256(json.dumps([case,asdict(physical),seed],sort_keys=True).encode())
    for name in ['gamma_spice.py','gamma_protocols.py','gamma_normalizer_statistics.py']:
        digest.update((ROOT/'python'/name).read_bytes())
    fingerprint=digest.hexdigest()
    cache=directory/'job_summary.json'
    if cache.exists():
        saved=json.loads(cache.read_text())
        if saved.get('fingerprint')==fingerprint:
            return settling_metadata(saved,case)
    result=run_spice(case,directory,physical,seed,exe,
        retain_raw=profile=='ideal' and case['name']=='positive_steps' and case['scheme']=='unit')
    predictor=sampled_predict(case,physical,realization(case,physical,seed))
    agreement=None
    if predictor is not None:
        agreement=all({str(k):v for k,v in predicted.items()}=={str(k):v for k,v in actual['observed_regions'].items()}
                      for predicted,actual in zip(predictor,result['checkpoints']))
    # Full phases/waveform residuals remain in verification.json, not duplicated.
    keep=['q','scheme','case','active_nodes','seed','physical','wrong_classes','class_observations',
          'network_class_failure','exact_class_value_failures','phase_class_failures',
          'max_abs_residual_units','peak_abs_voltage_v','headroom_violation',
          'observed_trigger_events','observed_unit_packets','expected_trigger_events',
          'induced_extra_trigger_events','last_active_phase_time_ns',
          'observation_stop_time_ns','boundary_requires_extension','elapsed_wall_s']
    summary={k:result[k] for k in keep}
    summary.update(profile=profile,severity=severity,fingerprint=fingerprint,
        sampled_charge_class_agreement=agreement,
        run_directory=str(directory.relative_to(ROOT)))
    settling_metadata(summary,case)
    write_json(cache,summary)
    return summary


def jobs_for_sweeps(trials):
    jobs=[]
    for q in (3,4,5):
        for scheme in SCHEMES:
            for ci,(name,inject,nodes) in enumerate(default_cases(q)):
                case=golden_case(q,name,inject,nodes,scheme,32)
                export_golden(case,ROOT/'spice'/'golden')
                # Compile safe physical units from a declared injection L1
                # budget, not from a measured/normalized answer. Local rewrites
                # cannot increase global L1; reserve 10% voltage headroom.
                budget=sum(abs(event['delta_units']) for event in case['stimuli'])
                safe_scale=max(4,budget*0.1/(0.9*0.6))
                jobs.append((case,Physical(scale=safe_scale),q*100+ci,ROOT/'spice'/'runs'/f'i{q}{scheme}{ci}', 'ideal','ideal'))
        case=golden_case(q,'mc_simultaneous', [{-1:q,0:2*q,1:-q}],5,'unit',12)
        for scale in (4,32,128):
            for severity in ('low','high'):
                low=severity=='low'
                profiles={
                    'threshold_offset':{'threshold_offset_v':0.0002 if low else 0.003},
                    'comparator_offset':{'comparator_offset_sigma_v':0.0002 if low else 0.002},
                    'capacitor_mismatch':{'capacitor_sigma_rel':0.005 if low else 0.05},
                    'transfer_mismatch':{'transfer_sigma_rel':0.005 if low else 0.05},
                    'pulse_width_mismatch':{'pulse_width_sigma_rel':0.005 if low else 0.05},
                    'leakage':{'leakage_tau_ns':10000 if low else 200},
                    'gain_error':{'injection_gain_sigma_rel':0.005 if low else 0.05},
                    'voltage_noise':{'sense_noise_sigma_v':0.0002 if low else 0.002},
                    'current_noise':{'current_noise_sigma_a':0.2e-9 if low else 2e-9},
                    'correlated_threshold':{'correlated_threshold_sigma_v':0.0002 if low else 0.002},
                    'combined':{'sense_noise_sigma_v':0.0002 if low else 0.002,
                        'capacitor_sigma_rel':0.005 if low else 0.05,
                        'transfer_sigma_rel':0.005 if low else 0.05,
                        'pulse_width_sigma_rel':0.005 if low else 0.05,
                        'comparator_offset_sigma_v':0.0001 if low else 0.001,
                        'correlated_threshold_sigma_v':0.0001 if low else 0.001,
                        'injection_gain_sigma_rel':0.005 if low else 0.05}}
                for pi,(profile,kw) in enumerate(profiles.items()):
                    physical=Physical(scale=scale,**kw)
                    for trial in range(trials):
                        # Fixed independent device draw per trial, paired across S.
                        seed=q*100000+pi*1000+trial
                        physical_trial=replace(physical,threshold_offset_v=physical.threshold_offset_v*(-1 if trial%2 else 1))
                        directory=ROOT/'spice'/'runs'/f'm{q}_{scale}_{pi}_{severity[0]}_{trial}'
                        jobs.append((case,physical_trial,seed,directory,profile,severity))
    return jobs


def aggregate_spice(records):
    groups={}
    for r in records:
        key=(r['q'],r['physical']['scale'],r['profile'],r['severity'])
        groups.setdefault(key,[]).append(r)
    rows=[]
    for (q,scale,profile,severity),group in sorted(groups.items()):
        n=len(group);k=sum(r['network_class_failure'] for r in group)
        rows.append({'q':q,'scale':scale,'profile':profile,'severity':severity,
            **confidence(k,n),'node_class_errors':sum(r['wrong_classes'] for r in group),
            'arithmetic_failure_trials':sum(r['exact_class_value_failures']>0 for r in group),
            'node_class_observations':sum(r['class_observations'] for r in group),
            'max_residual_units':max(r['max_abs_residual_units'] for r in group),
            'headroom_violations':sum(r['headroom_violation'] for r in group),
            'extra_trigger_events_total':sum(r['induced_extra_trigger_events'] for r in group),
            'max_last_activity_ns':max((r['last_active_phase_time_ns'] or 0 for r in group)),
            'max_operational_settling_time_after_last_injection_ns':max(r['operational_settling_time_after_last_injection_ns'] for r in group),
            'failure_modes':{mode:sum(r['failure_mode']==mode for r in group) for mode in sorted({r['failure_mode'] for r in group})},
            'statistical_note':'Nonideal groups have independent device/noise seeds within each condition. Ideal-case heterogeneous counts are deterministic coverage, not a Monte Carlo reliability bound.'})
    return rows


def scaling_study(destination,rare_trials):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    rows=[]
    for q in (3,4,5):
        for combined in (False,True):
            for scale in (1,2,4,8,16,24,32,48,64,128,256):
                path=destination/f'scale_q{q}_{int(combined)}_{scale}.json'
                if path.exists():
                    row=json.loads(path.read_text())
                else:
                    row=charge_monte_carlo(q,scale,10000,seed=700+q,combined=combined)
                    write_json(path,row)
                rows.append(row)
        if rare_trials:
            # One declared conservative grid point, not an adaptive maximum.
            path=destination/f'rare_q{q}_{rare_trials}.json'
            if path.exists():row=json.loads(path.read_text())
            else:
                print(f'Rare-event sampled-model run: q={q}, S=16, n={rare_trials:,}',flush=True)
                row=charge_monte_carlo(q,16,rare_trials,seed=900+q,combined=True)
                row['rare_event_confirmation']=True;write_json(path,row)
            rows.append(row)
    certified=[]
    for q in (3,4,5):
        for target in (1e-3,1e-6):
            valid=[r['scale'] for r in rows if r['q']==q and r['profile']=='voltage_noise_plus_relative_mismatch'
                   and r['one_sided_95pct_upper']<=target and not r['headroom_failures_at_injection']]
            certified.append({'q':q,'target_network_error':target,
                'largest_tested_scale_pointwise_bound_passes':max(valid) if valid else None,
                'source':'sampled charge-domain model only; not a SPICE/manufacturing Smax',
                'qualification':'Largest certified point on tested grid; no continuous S maximum or simultaneous multiple-comparison guarantee.'})
    return {'records':rows,'certified_grid_points':certified,
            'spice_tail_certification':'Circuit Monte Carlo sample sizes are insufficient to certify either 1e-3 or 1e-6; see exact upper bounds per condition.'}


def figures(summary,destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    destination=Path(destination)
    fig,axes=plt.subplots(1,3,figsize=(12,4.5),sharey=True)
    for ax,q in zip(axes,(3,4,5)):
        for profile,label,color in [('voltage_noise_only','Sense noise only','#1b6ca8'),('voltage_noise_plus_relative_mismatch','Noise + mismatch','#ca6c14')]:
            rows=sorted([r for r in summary['scaling']['records'] if r['q']==q and r['profile']==profile and not r.get('rare_event_confirmation')],key=lambda r:r['scale'])
            ax.plot([r['scale'] for r in rows],[r['one_sided_95pct_upper'] for r in rows],'-o',ms=3,color=color,label=label+' (95% upper)')
        rare=[r for r in summary['scaling']['records'] if r['q']==q and r.get('rare_event_confirmation')]
        if rare:ax.scatter([16],[rare[0]['one_sided_95pct_upper']],marker='v',color='#222222',label='Rare-model confirmation')
        actual=[r for r in summary['spice_groups'] if r['q']==q and r['profile']=='voltage_noise' and r['severity']=='low']
        ax.scatter([r['scale'] for r in actual],[r['one_sided_95pct_upper'] for r in actual],marker='s',facecolors='none',edgecolors='#555555',label='Actual SPICE (95% upper)')
        ax.axhline(1e-3,color='#777',ls=':',lw=.8);ax.axhline(1e-6,color='#777',ls=':',lw=.8)
        ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_ylim(5e-7,1)
        ax.set_xlabel('Physical scale S');ax.set_title(f'q={q}');ax.grid(True,alpha=.2)
    axes[0].set_ylabel('Network class-error probability\n(one-sided 95% upper bound)')
    handles,labels=axes[2].get_legend_handles_labels()
    fig.legend(handles,labels,fontsize=7,ncol=4,loc='upper center',bbox_to_anchor=(.5,.92),frameon=False)
    fig.suptitle('Synthetic parameters; accelerated-model reliability is separate from SPICE evidence')
    fig.tight_layout(rect=(0,0,1,.83));fig.savefig(destination/'scaling_reliability.png',dpi=170);plt.close(fig)
    selected=next(r for r in summary['spice_records'] if r['q']==4 and r['scheme']=='unit' and r['case']=='positive_steps' and r['profile']=='ideal')
    run_dir=ROOT/selected['run_directory']
    raw=run_dir/'normalizer.raw'
    if raw.exists():
        from gamma_spice import read_ascii_raw,node_name
        names,values=read_ascii_raw(raw)
        case=json.loads((run_dir/'case.json').read_text())
        origin=max(event['time_ns'] for event in case['stimuli'])
        time_ns=values[:,0]*1e9-origin
        take=(time_ns>=-1)&(time_ns<=8)
        fig,ax=plt.subplots(figsize=(8,4))
        for k,label,color,style in [(0,'Center','#bc4328','-'),(-1,'Left neighbor','#1b6ca8','-'),(1,'Right neighbor','#208652','--')]:
            ax.plot(time_ns[take],values[take,names.index(f'v({node_name(k)})')]*1e3,label=label,color=color,ls=style,lw=2)
        ax.axvspan(2,2.2,color='#777',alpha=.14,label='Sample / hold')
        ax.axvspan(2.6,3.6,color='#d6a34a',alpha=.16,label='Fixed charge transfer')
        ax.text(.65,106,'Raw: (0, 4, 0)',fontsize=10)
        ax.text(4.5,75,'Settled: (1, 0, 1)',fontsize=10)
        ax.set(xlabel='Time since final input pulse (ns)',ylabel='Physical capacitor state (mV)',
            title='Actual LTspice trace: q=4 local charge-packet rewrite',ylim=(-5,115))
        ax.grid(True,alpha=.2);ax.legend(loc='upper right',fontsize=8)
        fig.tight_layout();fig.savefig(destination/'q4_local_rewrite.png',dpi=170);plt.close(fig)


def report(summary,path):
    ideal=[r for r in summary['spice_records'] if r['profile']=='ideal']
    nonideal=[r for r in summary['spice_records'] if r['profile']!='ideal']
    agree=[r for r in summary['spice_records'] if r['sampled_charge_class_agreement'] is not None]
    out=['# Local gamma-normalizer engineering phase','', 'Run date: 2026-10-23. Actual simulator: installed LTspice XVII; no FPGA or full crossbar.', '',
        '## Outcome','',f"Recommended next behavioral prototype: **q={summary['recommendation']['recommended_next_prototype_q']}**, under the explicit score weights and energy assumptions below. This is not a silicon winner or a demonstrated energy advantage.", '',
        f"- Ideal circuit cases: {len(ideal)}; final class failures: {sum(r['network_class_failure'] for r in ideal)}; phase-region failures: {sum(r['phase_class_failures'] for r in ideal)}.",
        f"- Nonideal actual-SPICE trials: {len(nonideal)}; trials with a wrong output class: {sum(r['network_class_failure'] for r in nonideal)}. These deliberately mixed low/high stresses are not a manufacturing failure probability.",
        f"- Of those strict state failures, {sum(r['exact_class_value_failures']>0 for r in nonideal)} are arithmetic errors/unresolved overflow and {sum(r['network_class_failure'] and not r['exact_class_value_failures'] for r in nonideal)} emit a different but algebraically equivalent digit word.",
        f"- Sampled charge-model cross-checks against identical circuit realizations: {sum(r['sampled_charge_class_agreement'] for r in agree)}/{len(agree)} exact agreement of checkpoint classes. Leakage/current noise are excluded from this model comparison.",
        f"- Captured active-window boundary transfers: {sum(r['boundary_requires_extension'] for r in ideal)} ideal cases. Passive guard receivers preserve the exact value; these cases are not successes of a sealed three-node system.", '',
        'The behavioral circuit uses capacitor states, finite comparator bands, a clocked switch/sample-hold and event-controlled current packets. Python provides only injection events and a fixed clock budget; it does not provide firing decisions. Readout uses the declared small digit/overflow regions. Continuous residuals are instrumentation. Exact output values are reconstructed only from emitted digit classes.', '',
        '![Actual q4 LTspice charge-packet trace](../spice/results/q4_local_rewrite.png)', '',
        '## Reproduce','', 'From python/:','', '```powershell', 'python -B -m pytest -q -p no:cacheprovider', 'python -B run_normalizer_phase.py --spice-trials '+str(summary['spice_trials_per_condition'])+' --rare-trials '+str(summary['rare_trials_per_q']), '```', '',
        'The runner caches completed jobs by inputs and source hashes. spice/golden contains CSV event streams and exact JSON phase/checkpoint traces. spice/runs contains each executed netlist, golden case, simulator log and region verification. Selected ideal raw waveforms are retained; discarded sweep waveforms can be regenerated from their retained seeds and decks.', '',
        '## Protocol tradeoffs','', '| q | Scheme | Median rounds | Median triggers | Median unit rewrites | Overflow comparators/node |', '| --- | --- | --- | --- | --- | --- |']
    from gamma_protocols import protocol_geometry
    for s in summary['energy']['workload_summaries']:
        g=protocol_geometry(s['q'],s['scheme'])
        out.append(f"| {s['q']} | {s['scheme']} | {s['rounds_median']:g} | {s['trigger_events_median']:g} | {s['unit_rewrites_median']:g} | {g['overflow_comparators']} |")
    out += ['', 'Workload: 128 identical seeded INT8 scalar product pairs per candidate. Bulk is an unrestricted-quotient lower-bound control, excluded from circuit and energy claims. Capped schemes can reduce control events/rounds without reducing total moved charge proportionally. Comparator evaluations include idle synchronous observations; final digit readout is additional.', '',
        'Transient-limit sweeps are in the JSON. They stream signed injections, defer normalization until the limit, and always flush to the persistent alphabet at the end. A larger temporary alphabet trades fewer normalization sessions for more range/classes; it is not a free precision gain.', '',
        '## Physical model and boundaries','',
        'Reference parameters are illustrative and unoptimized: C=100 fF, unscaled unit step=0.1 V, packet top duration near 1 ns, phase period=4 ns. Ideal S is at least 4 and is compiled from the declared injection L1 budget with 10% margin to a 0.6 V absolute state limit; it is not selected from an observed answer. This is testbench range planning, not an input-normalization operation inside a weight-only MAC. Signed voltages stand for differential storage. Neither these component values nor the ideal behavioral sources establish a competitive energy implementation. The sample/hold, switch, reference generator and mirrored current sources are idealized.', '',
        'Persistent alphabet is radius q-1 (5/7/9 states). A unit event removes q identical packets at the center and injects one at each neighbor, with all physical units divided by the same S. Positional significance is carried by node index, never exponential physical voltage. Three/five/nine refers to active normalizers; two extra passive receivers explicitly collect outward boundary transfers.', '',
        '## Nonideal sweeps and scaling','',
        'Actual SPICE sweeps include threshold offset, independent comparator offset, capacitor mismatch, center/neighbor transfer mismatch, pulse-width mismatch, leakage, injection gain error, sampled voltage noise, piecewise-linear injected current noise, correlated threshold-bank error, and a combined condition. The MC input produces simultaneous neighboring firings and secondary events. Per-condition counts, confidence bounds, residuals, headroom and extra events are in the JSON.', '',
        'Comparator offsets are input-referred per active cell and shared by that cell\'s signed threshold pair; the correlated component is shared across cells. Independent mismatch of every threshold inside a multiband detector is not covered by this unit-cell MC study. Final class discrimination uses the same declared offset/noise model. No canonical-amplitude output driver or subsequent learned-weight MAC is simulated; capacitor residuals remain real analog errors until an explicit regeneration driver is implemented.', '',
        '| q | Error target | Largest tested scale with pointwise bound passing | Evidence |', '| --- | --- | --- | --- |']
    for row in summary['scaling']['certified_grid_points']:
        out.append(f"| {row['q']} | {row['target_network_error']:g} | {row['largest_tested_scale_pointwise_bound_passes']} | sampled charge model only |")
    out += ['', 'No actual-SPICE condition has enough independent trials to certify 1e-3 or 1e-6. The accelerated-model table has pointwise one-sided exact binomial confidence bounds, not a continuous Smax, not simultaneous bounds across all scales, and not a manufacturing guarantee. Its scope excludes continuous leakage/current noise, detailed analog sampling dynamics and transistor reference errors. Its 1e-6 confirmation, when enabled, uses a predeclared conservative S=16 rather than an adaptively selected maximum.', '',
        '![Scaling confidence bounds](../spice/results/scaling_reliability.png)', '',
        '## Energy break-even','',
        summary['energy']['equation'], '',
        'Common array energy cancels only when explicitly assumed equal. The grid includes incremental array cost and baseline ADC conversions per scalar product of 1/64, 1/16, 1, 4 and 16. Subunit counts make baseline conversion amortization explicit; the gamma counters still normalize each scalar product. No inference-level dot-product energy win can be inferred from this local workload.', '',
        'Two scheduling cases are exposed: until-settled plus confirmation (requires unpriced completion logic), and fixed 32-phase observation. Unit-rewrite energy scales with the assumed packet count q+2; trigger control and detection are separate. Routing, scale/sidecar and extra array costs remain explicit parameters. No single ADC energy is assumed.', '',
        'Detector activity is charged across a fixed compiled output/normalization layout: 21 positions at q=3 and 17 at q=4/5, from full INT8 encoding windows and an absolute-coefficient envelope. It does not shrink to each observed answer\'s active span. The global resource count is an extrapolation of local cells, not a SPICE simulation of a complete 17/21-node MAC.', '',
        '![Parameterized energy break-even](../spice/results/energy_break_even.png)', '',
        '## Candidate recommendation','', '| q | Combined score (lower is preferred) | Median compiled half-margin, V | Reference fixed-budget extra cost, fJ | Stress failure fraction |', '| --- | --- | --- | --- | --- |']
    for row in summary['recommendation']['candidates']:
        out.append(f"| {row['q']} | {row['combined_score']:.4f} | {row['median_compiled_half_margin_v']:.6g} | {row['fixed_budget_reference_extra_fj']:.3g} | {row['stress_failure_fraction']:.3f} |")
    out += ['', 'Weights: '+json.dumps(summary['recommendation']['score_weights'],sort_keys=True)+'. All metrics are normalized to the largest candidate value. These weights are an engineering preference, not a physical law.', '',
        summary['recommendation']['reference_energy_parameters']+'. Energy-only sensitivity winner counts: '+json.dumps(summary['recommendation']['energy_only_sensitivity_winner_counts'])+'.', '',
        'The parameterized score includes representation, product count, charge-transfer and comparator cost, persistent regions, transient range, weight-compiled class margin, and stress stability. Change the assumptions before treating the selection as a chip decision.', '',
        '## Floating-point sidecar','',
        'The existing exact NativeFloat/block-dot model remains in place. Exponents are digital metadata; exponent addition/comparison, signed residual/tile merge, positional renormalization and finite-window truncation are covered by regression tests. Truncation drops out-of-window coefficients explicitly; no exponent-CIM hardware is introduced.', '',
        '## Answer to the physical feasibility question','',
        'Local coarse normalization is behaviorally stable on the tested ideal schedules and can be verified from region outputs without a broad ADC. Actual transistor cost below ADC plus binary accumulation is still unestablished. The present physical requirements to resolve next are matched q:1:1 charge transfer, reliable scaled references, held comparator decisions, finite-boundary management, comparator/control activity and the cost of fixed clock observations. These are measurable engineering requirements, not yet an identified impossibility or a proven accelerator advantage.', '']
    Path(path).write_text('\n'.join(out),encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spice-exe',type=Path,default=DEFAULT_EXE)
    parser.add_argument('--spice-trials',type=int,default=12)
    parser.add_argument('--rare-trials',type=int,default=3000000)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if not args.spice_exe.is_file():parser.error('LTspice executable unavailable; circuit results cannot be substituted with Python')
    if args.spice_trials<1 or args.rare_trials<0 or args.workers<1:parser.error('invalid run volumes')
    out=ROOT/'spice'/'results';out.mkdir(parents=True,exist_ok=True)
    print('Exact protocol and transient studies',flush=True)
    workloads=workload_study();transient=transient_study()
    write_json(out/'protocol_workloads.json',{'records':workloads,'transient_records':transient})
    jobs=jobs_for_sweeps(args.spice_trials)
    print(f'Actual LTspice jobs: {len(jobs)}',flush=True)
    records=[];failures=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending={pool.submit(run_job,job,args.spice_exe):job for job in jobs}
        for future in as_completed(pending):
            job=pending[future]
            try:records.append(future.result())
            except Exception as e:
                failures.append({'directory':str(job[3]),'error':str(e)})
            if (len(records)+len(failures))%50==0:
                print(f'Circuit jobs complete: {len(records)+len(failures)}/{len(jobs)}; execution failures={len(failures)}',flush=True)
    write_json(out/'spice_execution_failures.json',failures)
    if failures:
        write_json(out/'partial_spice_records.json',records)
        raise SystemExit(f'{len(failures)} circuit jobs failed to execute; retained logs and partial records')
    records.sort(key=lambda r:(r['q'],r['profile'],r['physical']['scale'],r['severity'],r['case'],r['scheme'],r['seed']))
    print('Circuit sweeps complete; accelerated scaling and energy studies',flush=True)
    scale=scaling_study(out/'scaling',args.rare_trials)
    energy=energy_study(workloads,out)
    summary={'schema_version':1,'date':'2026-10-23','spice_trials_per_condition':args.spice_trials,
        'rare_trials_per_q':args.rare_trials,'spice_records':records,'spice_groups':aggregate_spice(records),
        'scaling':scale,'energy':energy,'transient_records':transient,
        'recommendation':recommendation(energy,records,scale)}
    write_json(out/'normalizer_phase_summary.json',summary)
    figures(summary,out);report(summary,ROOT/'docs'/'normalizer_phase_results.md')
    print('Recommended next prototype q:',summary['recommendation']['recommended_next_prototype_q'],flush=True)
    print('Saved docs/normalizer_phase_results.md and spice/results/normalizer_phase_summary.json',flush=True)
    ideal=[r for r in records if r['profile']=='ideal']
    if any(r['network_class_failure'] or r['phase_class_failures'] for r in ideal):
        raise SystemExit('Ideal circuit validation failed; inspect retained cases before physical conclusion')


if __name__=='__main__':main()
