"""Report literal SPICE evidence; no retuning and no oracle-driven hardware."""
from pathlib import Path
import json,csv,hashlib,zipfile,shutil
import numpy as np
from sky130_self_topple import ROOT,OUT,analyze
from gamma_sky130 import fast_raw

def refresh(directory):
    result=json.loads((directory/'result.json').read_text())
    names,x=fast_raw(directory/'output.raw')
    row=analyze(result,names,x);row['run']=str(directory.relative_to(ROOT))
    (directory/'metrics.json').write_text(json.dumps(row,indent=2))
    return row,names,x

def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])

def main():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    selected=[];traces=[]
    for s in (1,-1):
        row,names,x=refresh(OUT/'spice'/f'selected_s{s}');selected.append(row);traces.append((names,x))
    details=json.loads((OUT/'characterization.json').read_text())
    allrows=selected+details
    pvt=[r for r in details if r['cfg'].get('group')=='pvt']
    mc=[r for r in details if r['cfg'].get('group')=='mismatch']
    repeats=[r for r in details if r['cfg'].get('group')=='repeat']
    threshold=[r for r in details if r['cfg'].get('group')=='threshold']
    fine=[r for r in details if r['cfg'].get('group')=='timestep']
    nets={seed:[r for r in mc if r['cfg']['seed']==seed] for seed in sorted({r['cfg']['seed'] for r in mc})}
    service=[];dynamic=[];transfer=[];gross=[]
    for r in selected:
        o=r['observations'][0];ph=o['phase_supply_energy_fj']
        service.append(sum(sum(ph[k].values()) for k in ('before_firing_fj','active_transfer_fj','reset_and_rearm_fj')))
        dynamic.append(o['incremental_over_local_idle_fj']);transfer.append(o['energy_fj']['transfer'])
        gross.append(sum(o['gross_supply_fj'].values()))
    meanservice=float(np.mean(service));meandynamic=float(np.mean(dynamic));meanwindow=float(np.mean([r['event_window_mean_fj'] for r in selected]))
    baseline=json.loads((ROOT/'results/sky130_phase/phase_summary.json').read_text())
    comparison=[]
    for row in baseline['primary_budget_rows']:
        n=row['accumulation_length'];firings=n*32.3828125
        comparison.append({'N':n,'firings':firings,'ceiling_fj':row['rewrite_ceiling_fj_if_other_norm_costs_zero'],
          'measured_service_fj':meanservice,'service_over_ceiling':meanservice/row['rewrite_ceiling_fj_if_other_norm_costs_zero'],
          'normalizer_projection_pj':firings*meanservice/1000,
          'boundary_plus_service_projection_pj':(firings*meanservice+row['boundary_cost_per_position_fj']*row['gamma_positions'])/1000,
          'optimistic_dynamic_projection_pj':firings*meandynamic/1000,
          'adc_pj':row['adc_energy_per_output_fj']/1000})
    summary={'date':'2026-10-24','status':'SELF_TRIGGERED_CHARGE_TRANSFER_DEMONSTRATED_CANDIDATE_FAILS_ENERGY_AND_ROBUSTNESS',
      'q':4,'mos_count':selected[0]['mos_count'],'storage_ff':100,'unit_q_fc':2.5,
      'polling_in_circuit':False,'ideal_detector_or_transfer':False,'nominal_charge_passes':sum(r['charge_passed'] for r in selected),
      'nominal_voltage_passes':sum(r['passed'] for r in selected),'nominal_cases':2,
      'signed_supply_window_mean_fj':meanwindow,'gross_supply_window_mean_fj':float(np.mean(gross)),
      'signed_supply_input_to_rearm_mean_fj':meanservice,'local_idle_subtracted_window_mean_fj':meandynamic,
      'transfer_rail_signed_window_mean_fj':float(np.mean(transfer)),
      'idle_power_nw':float(np.mean([r['idle_power_nw'] for r in selected])),
      'pvt_charge_passes':sum(r['charge_passed'] for r in pvt),'pvt_cases':len(pvt),
      'mismatch_networks':len(nets),'mismatch_networks_passing_both_signs':sum(all(r['charge_passed'] for r in group) for group in nets.values()),
      'comparison':comparison,
      'limitations':['Not a functional full normalizer: isolated neighbor nodes lack their own detector/transfer loads.',
        'Nominal branch-charge tolerance is 0.1 Q per leg; voltage-digit checks are separate and fail.',
        'No process-qualified resistor/capacitor variation, transient noise, extracted routing or layout.',
        'Existing workload counts are a retained per-product policy; no electrically replayed full dot.',
        'All-in projected totals still omit array expansion, other routing, metadata and completion.',
        'A continuously biased detector is one implementation; no universal energy lower bound is claimed.']}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    with (OUT/'comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(comparison[0]));w.writeheader();w.writerows(comparison)
    flat=[]
    for r in allrows:
        obs=r['observations'];flat.append({'label':r['cfg']['label'],'group':r['cfg'].get('group','selected'),
          'charge_passed':r['charge_passed'],'voltage_passed':r['passed'],'firings':r['firings'],'idle_nw':r['idle_power_nw'],
          'window_fj':r['event_window_mean_fj'],'worst_charge_error_q':max([o['max_charge_error_units'] for o in obs],default=0),
          'worst_voltage_error_digits':r['worst_packet_error_units'],'run':r['run']})
    with (OUT/'characterization.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    fig,axes=plt.subplots(3,2,figsize=(12,9),sharex=True)
    for j,((names,x),r) in enumerate(zip(traces,selected)):
        ix={v:k for k,v in enumerate(names)};t=x[:,0]*1e9-2500;m=(t>=-30)&(t<=900)
        for node in ('center','left','right'):axes[0,j].plot(t[m],x[m,ix['v('+node+')']],label=node)
        axes[0,j].axhline(.7,color='grey',ls=':',lw=1);axes[0,j].axhline(.5,color='grey',ls=':',lw=1)
        axes[0,j].set_title('Positive overflow' if j==0 else 'Negative overflow');axes[0,j].set_ylabel('Physical voltage (V)');axes[0,j].legend(fontsize=8)
        s='pos' if j==0 else 'neg'
        for node in (s+'req',s+'fire',s+'reset',s+'armed'):axes[1,j].plot(t[m],x[m,ix['v('+node+')']],label=node)
        axes[1,j].set_ylabel('Control voltage (V)');axes[1,j].legend(fontsize=8,ncol=2)
        for comp in ('detector','logic','transfer'):
            p=-x[:,ix['v(dd_'+comp+')']]*x[:,ix['i(v'+comp+')']]
            axes[2,j].plot(t[m],p[m]*1e6,label=comp)
        axes[2,j].set_ylabel('Signed rail power (uW)');axes[2,j].set_xlabel('Time after additive stimulus begins (ns)');axes[2,j].legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'waveforms.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4))
    ns=[r['N'] for r in comparison];ax.loglog(ns,[r['ceiling_fj'] for r in comparison],'o-',label='Available fJ per firing; all other extra costs free')
    ax.axhline(meanservice,color='#b33333',label=f'Measured candidate service: {meanservice:.0f} fJ')
    ax.axhline(float(np.mean(transfer)),color='#b67923',ls='--',label='Measured transfer rail alone; signed window mean')
    ax.set_xlabel('Useful MACs per dot');ax.set_ylabel('fJ per firing');ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'budget.png',dpi=160);plt.close(fig)
    nom=[]
    for r,svc in zip(selected,service):
        o=r['observations'][0];q=o['branch_charge_units'];v=o['packet_units_apparent'];f=o['firings'][0]
        nom.append(['+' if o['sign']>0 else '-',f"{q['center']:.3f}, {q['left']:.3f}, {q['right']:.3f}",f"{v['center']:.3f}, {v['left']:.3f}, {v['right']:.3f}",f"{o['total_supply_fj']/1000:.3f}",f"{svc/1000:.3f}",f"{o['rearm_time_ns']-o['start_ns']:.1f}"])
    budgetrows=[[r['N'],f"{r['firings']:.0f}",f"{r['ceiling_fj']:.3f}",f"{r['normalizer_projection_pj']/1000:.2f}",f"{r['adc_pj']:.2f}"] for r in comparison]
    holdrows=[]
    for r in threshold:
        o=r['observations'][0];holdrows.append([r['cfg']['signs'][0],r['cfg']['input_units'],f"{(o['baseline_center_v']-.6)/.025:.2f}",r['firings']])
    repeatrows=[[r['cfg']['signs'][0],r['firings'],str([o['charge_passed'] for o in r['observations']])] for r in repeats]
    lines=[
      '# Self-triggered q=4 normalizer: transistor experiment',
      '2026-10-24. ngspice 47, pinned SKY130A 1.8 V transistor and MIM models. This phase continues the matched-130nm report; all earlier reports are preserved.',
      '## Result',
      '**Self-triggering is demonstrated without polling. This particular circuit fails the energy and robustness gates. A robust all-in Gamma advantage is not established.**',
      f'The selected circuit has {summary["mos_count"]} MOSFETs, two continuously biased overflow detectors, guarded SR event latches, physical RC reset paths and gated weighted current mirrors. It uses 100 fF storage and Q = 2.5 fC. There is no behavioral threshold, event controller, ADC or ideal packet-current source inside the cell. The sole independent current stimulus adds test charge; it never receives an oracle answer.',
      'The two nominal isolated-cell tests each produce one firing. Charge transfers meet the declared 0.1 Q absolute tolerance on each leg. That tolerance is an engineering test, not exact arithmetic or a process-qualified error rate. Physical voltage changes fail a separate 0.1-digit check because the center and passive neighbor nodes have different parasitic loads. Neither success may be substituted for the other.',
      table(['Overflow','Measured branch charge: center, left, right (Q)','Voltage change / 25 mV: center, left, right','1.2 us window supply (pJ)','Input-to-rearm supply (pJ)','Re-arm after input begins (ns)'],nom),
      '## What was measured',
      'Every circuit rail is a separate 1.8 V source: detector, logic, packet bias, transfer and references. Signed supplied energy is -integral(Vrail * I(Vsource) dt), using ngspice current orientation. Gross delivered energy is also retained before subtracting returned work. The event window begins with additive input and lasts 1.2 us; input-to-rearm stops when request and reset are low, the RC node is below 0.1 VDD and the guard has re-armed. All RC reset and latch work inside that interval is counted.',
      f'Mean signed window energy is **{meanwindow/1000:.3f} pJ/firing**; gross delivered window energy is {np.mean(gross)/1000:.3f} pJ. Mean signed input-to-rearm energy is **{meanservice/1000:.3f} pJ/firing**. Subtracting the measured local pre-event idle baseline from the whole window gives {meandynamic/1000:.3f} pJ, a diagnostic rather than a universal dynamic-event constant: idle power depends on state.',
      f'At the quiet initial state, idle power is **{summary["idle_power_nw"]:.1f} nW/cell**. A 400 ns quiet interval therefore costs {summary["idle_power_nw"]*.4:.1f} fJ/cell. There are no periodic comparisons, but this continuously biased detector is not energy-free while idle. References are real resistor ladders; their cost is neither omitted nor divided by an assumed fanout.',
      'Per-rail, active-transfer, reset/re-arm, startup, command-pin and external additive-input work are retained in metrics.json. Input work is separate from normalizer rail work; it belongs to the omitted MAC/stimulus implementation, not an uncounted packet-transfer source.',
      '## Why it costs this much',
      f'The signed transfer-rail window costs are {transfer[0]:.2f} fJ for positive overflow and {transfer[1]:.2f} fJ for negative overflow. Detector output transitions and the slowly changing RC reset drive CMOS gates through their conducting region; those detector/logic costs dominate. This is a measured limitation of this topology, not a 100 fJ comparator proxy and not a lower bound on every possible threshold element.',
      'The initial 73-transistor one-shot could re-fire while an old request was still high. The selected 101-transistor version adds an armed latch that requires that request to clear before accepting another event. The nominal event is therefore not being credited with free re-arm logic.',
      'The earlier 25 fF experiments are retained, including missing fires, overshoot and opposite-polarity/repeated firing. Increasing storage to 100 fF and tuning fixed transistor widths makes nominal branch packets close to the target; this is not a proof of optimal capacitance or optimized topology.',
      '## Frozen-design robustness and re-arm',
      f'The declared PVT grid contains {len(pvt)} signed cases: five corners at 1.8 V/25 C, tt at -40/100 C, and tt at 1.62/1.9 V, both polarities. **{sum(r["charge_passed"] for r in pvt)}/{len(pvt)}** meet the nominal branch-charge criterion. Resistors and MIM coefficients remain nominal in this PVT grid; it is not a complete process-qualified passive-component sweep.',
      f'Eight PDK transistor-mismatch seeds were replayed with both input signs. **{sum(all(r["charge_passed"] for r in group) for group in nets.values())}/{len(nets)} independent networks** pass both signs. The two signs share a seed and count as one network, not independent yield samples. No yield or 1e-3/1e-6 reliability claim follows.',
      table(['Input sign','Input charge / nominal Q','Observed state just after input / 25 mV relative to 0.6 V','Firings'],holdrows),
      'The third column is measured, not assumed equal to supplied charge. These are hold/overflow probes, not a calibrated DC trip-point sweep. Detector latency and kickback permit excursions beyond the nominal +/-4 voltage thresholds.',
      table(['Input sign','Firings over two requests','Charge checks by request'],repeatrows),
      'Repeated requests test the guard and packet behavior locally. Passive neighbors are not independent normalizers and no full-word cascade is inferred. A 0.1 ns maximum timestep replay is retained alongside the 0.5 ns nominal run; numerical agreement is exported in summary.json.',
      '## Workload combination',
      'Self-triggering removes the historical 197.25 N scheduled observation count from this comparison. It does not remove supply bias. Retaining the earlier unit rewrite count gives N_firings = 32.3828125 N. The table multiplies that count by measured mean input-to-rearm supply work, with no new comparator-evaluation charge added.',
      table(['N','Retained firings/dot','Available fJ/firing if every other extra cost is free','Projected normalizer alone (nJ/dot)','Matched ADC baseline (pJ/dot)'],budgetrows),
      f'The most favorable N=128 allowance is {comparison[0]["ceiling_fj"]:.3f} fJ/firing; this candidate costs {comparison[0]["service_over_ceiling"]:.0f} times that allowance. Even the idle-subtracted diagnostic yields {comparison[0]["optimistic_dynamic_projection_pj"]/1000:.2f} nJ/dot at N=128. Adding the existing 4.81-6.73 pJ estimated Gamma readout cannot rescue this candidate.',
      '**These totals are cost projections, not a measured all-in dot product.** They still omit other routing/array expansion, metadata and global completion, and the cell fails physical-state and robustness checks. The measured local energy is enough to reject this implementation under the retained workload; it is not enough to certify a correct complete normalizer.',
      '## What remains open',
      'A different threshold/transfer topology could have lower cost. The user\'s self-triggering proposal remains distinct from clocked seven-state readout. This result specifically rules out treating continuously biased PMOS detectors, an RC/CMOS one-shot and rail-fed current-mirror packets as sufficiently cheap here. A future experiment needs uniform loaded nodes, a physical charge/state mapping, a near-zero-idle threshold mechanism, accurate self-limiting packets and fast low-energy re-arm. Its supply work must still fit 2.42 down to 0.098 fJ per retained firing; event counts may change only with a demonstrated bounded schedule.',
      '## Replay and evidence',
      '`python -B python/sky130_self_topple.py` reproduces both selected nominal tests. Add `--characterize` for the fixed robustness grid. selected_config.json freezes the selected dimensions; each case also has a standalone circuit.cir and seeded run.cir or replay.cir. report_self_topple.py only reads completed waveforms and writes the report. Historical broad sweeps are retained as exploratory controls and are not required to reproduce the decision.',
      'The bundle includes literal decks, model/simulator provenance, seeds, logs, all case metrics, selected and representative robustness waveforms, comparison CSVs, source files and hashes. All completed waveforms remain in the project; only representative waveforms are duplicated in the download. PDK files are unmodified. Fetch the pinned PDK and simulator from the recorded manifests for a replay outside this workspace.',
      'Device documentation: https://skywater-pdk.readthedocs.io/en/main/rules/device-details.html . Model provenance: spice/sky130/pdk_manifest.json in the project.',
      '## Plots',
      '![Selected signed overflow waveforms](waveforms.png)',
      '![Measured cost versus workload ceiling](budget.png)',
    ]
    # Timestep check is an observation on the executed traces, not a threshold re-fit.
    checks=[]
    for a,b in zip(selected,sorted(fine,key=lambda r:-r['cfg']['signs'][0])):
        oa=a['observations'][0];ob=b['observations'][0]
        checks.append({'sign':a['cfg']['signs'][0],'window_energy_relative_difference':abs(ob['total_supply_fj']/oa['total_supply_fj']-1),
          'max_branch_charge_difference_q':max(abs(ob['branch_charge_units'][n]-oa['branch_charge_units'][n]) for n in ('center','left','right'))})
    summary['timestep_checks']=checks
    summary['all_requested_cases_completed']=len(details)==43
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    report='\n\n'.join(lines)+'\n';(OUT/'engineering_report.md').write_text(report,encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
