"""New matched-node report/figures. Never regenerates historical reports."""
import json,math,statistics
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gamma_adc_baseline import ROOT,BASE,OUT,dump,sha,write_csv
from gamma_matched_130nm import boundary_cost


def table(headers,rows):
    safe=lambda v:str(v).replace('|','/').replace('\n',' ')
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(safe(x) for x in row)+' |' for row in rows])


def primary():
    budgets=pd.read_csv(OUT/'regeneration_budgets.csv')
    return budgets[(budgets.baseline_mode=='optimistic_anchor_envelope')&(budgets.physical_decisions_estimated==3)].sort_values('accumulation_length')


def figures():
    directory=OUT/'plots';directory.mkdir(exist_ok=True)
    def save(fig,name):fig.tight_layout();fig.savefig(directory/name,dpi=170);plt.close(fig)
    grid=pd.read_csv(OUT/'adc_130nm_pareto.csv');anchors=pd.read_csv(OUT/'adc_130nm_survey_anchors.csv')
    fig,ax=plt.subplots(figsize=(8,5))
    for b in (4,6,8,10):
        r=grid[(grid.enob==b)&(grid.n_adcs==1)&grid.within_model_hard_frequency_limit].sort_values('frequency_per_adc_hz')
        ax.plot(r.frequency_per_adc_hz,r.energy_per_conversion_pj,'o-',label=f'{b} ENOB, fitted 130nm')
    ax.scatter(anchors.sampling_clock_hz,anchors.survey_energy_per_nyquist_pj,c='black',marker='x',label='12 training anchors (mixed ENOB; fs vs P/fsnyq)')
    ax.set(xscale='log',yscale='log',xlabel='Per-ADC sample clock / requested conversion rate (Hz)',ylabel='pJ / conversion or survey Nyquist value',title='Published fit and 130nm anchors; clock/value distinction retained')
    ax.legend(fontsize=8);ax.grid(alpha=.2);save(fig,'adc130_energy_and_anchors.png')
    fig,ax=plt.subplots(figsize=(8,5))
    for b in (6,7,8,9,10):
        r=grid[(grid.enob==b)&(grid.aggregate_samples_s==160e6)]
        ax.plot(r.total_area_um2/1e6,r.energy_per_conversion_pj,'o-',label=f'{b} ENOB')
        for _,row in r.iterrows():ax.annotate(str(int(row.n_adcs)),(row.total_area_um2/1e6,row.energy_per_conversion_pj),xytext=(2,3),textcoords='offset points',fontsize=6)
    ax.set(xscale='log',ylabel='pJ / conversion',xlabel='Total ADC bank area (mm²)',title='160M values/s aggregate: labels show physical ADC count')
    ax.legend();ax.grid(alpha=.2);save(fig,'adc130_energy_area_frontier.png')
    p=primary();n=p.accumulation_length.to_numpy();adc=p.adc_energy_per_output_fj.to_numpy();pos=p.gamma_positions.to_numpy();boundary=p.boundary_cost_per_position_fj.to_numpy()*pos
    fig,ax=plt.subplots(figsize=(8,5));ax.plot(n,adc/n,'o-',label='Optimistic conventional ADC / useful MAC')
    ax.plot(n,boundary/n,'s--',label='Gamma readout-only estimate; normalization/control extra')
    costs=pd.read_csv(OUT/'matched_break_even.csv')
    c=costs[(costs.baseline_mode=='optimistic_anchor_envelope')&(costs.gamma_position_policy=='precision_bound')&(costs.detector_decisions_estimated==3)&(costs.reference_fanout_hypothesis==1)&(costs.binary_accumulation_plus_requantization_fj==0)&(costs.normalization_pricing=='legacy_UNCALIBRATED_costs')].sort_values('accumulation_length')
    ax.plot(c.accumulation_length,c.gamma_fj_per_useful_mac,'^:',label='Gamma with OLD UNCALIBRATED normalizer prices')
    ax.set(xscale='log',yscale='log',xlabel='Useful INT8 MACs / output dot',ylabel='fJ / useful scalar MAC',title='Matched work: full utilization, 6/7/8/9 ENOB\n2.5M outputs/s per column')
    ax.legend(fontsize=8);ax.grid(alpha=.2);save(fig,'matched_energy_per_mac.png')
    fig,ax=plt.subplots(figsize=(8,5));ax.plot(n,p.overflow_detection_ceiling_fj_if_other_norm_costs_zero,'o-',label='Ceiling: all other normalizer costs free')
    definition=json.loads((OUT/'baseline_definition.json').read_text());comp=definition['measured_boundary']['additional_comparison_fj']
    ax.axhline(comp,c='tab:red',ls='--',label=f'Measured class comparator proxy: {comp:.1f} fJ')
    ax.set(xscale='log',yscale='log',xlabel='Useful MACs / output',ylabel='fJ / overflow comparator evaluation',title='Retained per-product normalization schedule × N; not a fundamental limit')
    ax.legend(fontsize=8);ax.grid(alpha=.2);save(fig,'normalizer_cost_ceiling.png')
    surface=pd.read_csv(OUT/'normalization_break_even_surface.csv')
    fig,axes=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
    for ax,length in zip(axes.flat,n):
        r=surface[(surface.accumulation_length==length)&(surface.reference_fanout==1)&(surface.rewrite_fj==1.2)]
        pivot=r.pivot(index='normalization_activity_density',columns='overflow_detect_fj',values='baseline_minus_gamma_fj_per_dot')
        x=pivot.columns.to_numpy();y=pivot.index.to_numpy();z=pivot.to_numpy()/1000
        vmax=max(1,float(np.max(abs(z))));scale=matplotlib.colors.SymLogNorm(linthresh=1,vmin=-vmax,vmax=vmax)
        im=ax.pcolormesh(x,y,z,cmap='RdBu',norm=scale,shading='auto')
        if z.min()<0<z.max():ax.contour(x,y,z,levels=[0],colors='black',linewidths=1)
        ax.set(xscale='log',yscale='log',title=f'N={int(length)}',xlabel='Overflow decision cost (fJ)',ylabel='Normalization activity / retained schedule')
        fig.colorbar(im,ax=ax,label='Conventional − Gamma (pJ/dot)')
    fig.suptitle('Conditional region: positive favors Gamma; rewrite1.2fJ and control0.01fJ are UNCALIBRATED\nActivity reductions require a new bounded schedule; unbuilt decoder/routing cost=0',fontsize=10)
    save(fig,'normalization_break_even_surfaces.png')
    full=pd.read_csv(OUT/'full_sensitivity.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,4.5))
    for ax,length in zip(axes,(128,8192)):
        r=full[(full.accumulation_length==length)&(full.routing_fj_per_position==0)&(full.boundary_invocations_per_dot==1)&(full.write_refresh_fraction==1)&(full.unknown_logic_fj_per_position==0)]
        for factor in sorted(r.representation_array_energy_factor.unique()):
            a=r[r.representation_array_energy_factor==factor].sort_values('common_conventional_array_fj_per_mac')
            ax.plot(a.common_conventional_array_fj_per_mac,a.baseline_minus_gamma_fj_per_dot/1000,'o-',label=f'Array energy ×{factor:.2g}')
        ax.axhline(0,c='black');ax.set(title=f'N={length}; normalization FREE',xlabel='Conventional array fJ / useful MAC',ylabel='Conventional − Gamma (pJ/dot)',xscale='symlog',xlim=(0,105))
        ax.legend(fontsize=7);ax.grid(alpha=.2)
    save(fig,'array_expansion_sensitivity.png')
    measured=definition['measured_boundary'];fig,ax=plt.subplots(figsize=(8,4))
    keys=['reference','buffer','comparator','interface','control_pin_charging_upper']
    ax.barh([k.replace('_',' ') for k in keys],[measured[k] for k in keys]);ax.set(xlabel='fJ / position',title='Implemented hybrid boundary: one real threshold\nOther decisions and controller remain unbuilt')
    save(fig,'measured_boundary_components.png')


def main():
    smoke=json.loads((OUT/'sky130_smoke_summary.json').read_text());assert smoke['all_passed']
    regression=json.loads((OUT/'regression_result.json').read_text());assert regression['exit_code']==0
    s=json.loads((ROOT/'spice/sky130/results/sky130_phase_summary.json').read_text())
    definition=json.loads((OUT/'baseline_definition.json').read_text());measured=definition['measured_boundary'];p=primary()
    fanout=pd.read_csv(OUT/'reference_fanout.csv')
    figures()
    pricing=pd.read_csv(OUT/'matched_break_even.csv')
    key=(pricing.baseline_mode=='optimistic_anchor_envelope')&(pricing.gamma_position_policy=='precision_bound')&(pricing.detector_decisions_estimated==3)&(pricing.reference_fanout_hypothesis==1)&(pricing.binary_accumulation_plus_requantization_fj==0)
    case_table=table(['MACs/dot','ADC ENOB','ADC pJ/output','Gamma positions','Readout estimate pJ/output','Remaining pJ for all extra Gamma costs','Max overflow fJ if rewrite/control free'],
        [(int(r.accumulation_length),int(r.enob),f'{r.adc_energy_per_output_fj/1000:.3f}',int(r.gamma_positions),f'{r.boundary_cost_per_position_fj*r.gamma_positions/1000:.3f}',f'{r.remaining_normalizer_routing_logic_budget_fj_per_dot/1000:.3f}',f'{r.overflow_detection_ceiling_fj_if_other_norm_costs_zero:.4f}') for _,r in p.iterrows()])
    normtable=table(['MACs/dot','No-cost-normalizer bound pJ','OLD uncalibrated prices pJ','1fJ overflow hypothesis pJ','Measured class-comparator proxy pJ'],
        [(int(n),*[f"{pricing[key&(pricing.accumulation_length==n)&(pricing.normalization_pricing==mode)].iloc[0].gamma_peripheral_fj_per_dot/1000:.3f}" for mode in ('free_normalization_lower_bound','legacy_UNCALIBRATED_costs','one_fj_coarse_detector_HYPOTHESIS','measured_comparator_PROXY_rewrite_control_free')]) for n in p.accumulation_length])
    ftable=table(['Equivalent loads','Ladder Mohm','Correct classes','Worst write error/digit','Reference fJ/load'],
        [(int(r.fanout_equivalent_loads),f'{r.bank_resistance_ohm/1e6:.4g}',f'{int(r.correct_classes)}/{int(r.checkpoints)}',f'{r.peak_abs_write_error_digits:.3f}',f'{r.reference_fj_per_load_mean:.2f}') for _,r in fanout.iterrows()])
    pvt=s['followup'];energy=[r for r in pvt if r['group']=='pvt_energy'];joint=sum(r['wrong_next_mac']==0 and r['unresolved_checkpoints']==0 and r['soa_violations']==0 for r in energy)
    comp_est=boundary_cost(measured,3);skip_saving=boundary_cost(measured,3)-boundary_cost(measured,3,regeneration_fraction=0)
    mc=s['monte_carlo'];mctable=table(['Bracket','Independent networks','Networks with any wrong class','95% one-sided upper failure bound'],
        [(r['mode'],r['network_trials'],r['networks_with_any_wrong_class'],f"{100*r['one_sided_95pct_upper']:.1f}%") for r in mc])
    precisions=pd.read_csv(OUT/'gamma_precision_layout.csv');layouttable=table(['N','Target ENOB','Fixed mantissa positions','Representation-only empirical bits','Exact fullword positions'],
        [(int(r.mac_count),int(r.enob_target),int(r.guaranteed_truncated_positions),f'{r.empirical_quantization_equivalent_bits:.2f}',int(r.exact_fullword_positions)) for _,r in precisions.iterrows()])
    source=json.loads((OUT/'baseline_manifest.json').read_text())
    summary={'date':'2026-10-24','q':4,'status':'CALIBRATED_CONDITIONAL_HEADROOM_NOT_A_ROBUST_SYSTEM_WIN',
        'main_answer':'A readout-boundary energy opportunity remains against the optimistic130nm baseline. No demonstrated total-peripheral win: physical normalization/completion cost, full-dot throughput and PVT/mismatch accuracy remain unresolved.',
        'ADC':{'grid_points':source['grid_rows'],'training_130nm_anchors':12,'raw_130nm_survey_anchors':source['all_raw_130nm_anchors'],
            'published_example_reproduced':True,'four_inputs_confirmed':True,'legacy_and_successor_fit_identical':True,
            'primary_130nm_adc_pj':[float(x/1000) for x in p.adc_energy_per_output_fj],
            'old_100fj_budget_superseded':True},
        'restoration':{'matched_common_scale_residual_point35':'12/12 regenerated;1/12 direct;both12/12 at zero','all49_transitions_correct':True},
        'energy':{'known_hybrid_fj_per_position':measured['known_hybrid_fj'],'three_decision_boundary_estimate_fj':comp_est,
            'complete_measured_detector':False,'readout_word_pj':[float(x/1000) for x in p.boundary_cost_per_position_fj*p.gamma_positions],
            'remaining_extra_cost_budget_pj':[float(x/1000) for x in p.remaining_normalizer_routing_logic_budget_fj_per_dot],
            'overflow_detector_price_ceiling_fj_if_other_norm_costs_zero':[float(x) for x in p.overflow_detection_ceiling_fj_if_other_norm_costs_zero]},
        'robustness':{'energy_candidate_joint_pvt_passes':joint,'pvt_points':len(energy),'network_mismatch':mc,'full_physical_accuracy_matched':False},
        'dominant_cost':'Known hybrid: reference324fJ, buffer249fJ, comparator100fJ. Estimated3-decision detector312fJ. At retained per-product event counts, normalization decisions can dominate; actual overflow primitive energy is unmeasured.',
        'recommendation':'Do not build an array yet. Keep q4; characterize a low-energy bounded overflow/transfer primitive and fix driver/reference matching. Test sparse/shared normalization against explicit cost/activity ceilings; write skipping alone saves little static power.',
        'reference_fanout':fanout.to_dict(orient='records'),'software_regression':regression['output'].strip().splitlines()[-1],
        'reused_historical_circuits':351,'new_transistor_cases':6,'no_full_framework_or_ASIC_flow_installed':True,
        'limitations':definition['limitations'],'primary_budget_rows':p.to_dict(orient='records')}
    dump(OUT/'phase_summary.json',summary)
    text=f'''# Gamma-CIM: matched130nm engineering result

Date:2026-10-24. Lead q=4, beta=2+sqrt(3), seven persistent digits. Prior LTspice and SKY130 reports/data are preserved; this is a new calibration and comparison phase.

## Decision

**The old100fJ ADC /2.74fJ-position gate was not a technology-matched test and is superseded. The measured Gamma readout boundary is not automatically too expensive at130nm. A full peripheral-energy advantage is still unproved.**

At the declared matched work/rate, the optimistic conventional comparator costs14.86–32.64pJ per output. A5–7-position Gamma readout is estimated at4.81–6.73pJ before normalization, missing detector/control blocks and extra array/routing costs. The17-position historical product layout must not be compared directly with a6-bit final output without matching precision.

The remaining10.05–25.91pJ/output is a cost allowance, not measured savings. Under the retained per-product normalization schedule, overflow evaluation alone must cost below0.398 down to0.016fJ as N grows128→8192, even when rewrite and completion are free. The existing106fJ class comparator is a comparison point, not a mandatory overflow topology. A much cheaper overflow primitive and/or lower event density is required. Existing transistor PVT/mismatch failures also prevent a robust win claim.

## A — Calibrated ADC methodology

The pinned historical ADC plug-in and maintained successor have identical fitted coefficients. All four inputs are verified from executable code: ENOB, node, aggregate throughput and interleaved ADC count. The historical wrapper's example(8ENOB,16nm,5.12G aggregate values/s,32ADCs) reproduces0.465575720551pJ/conversion and46678.81um² total area. This16nm point is an API reproduction only; no Gamma advanced-node projection is made. The entire decision comparison uses130nm.

Numerical model output is pJ. ADCRequest has a stale Joules docstring, but the training CSV, Accelergy p-unit wrapper and successor's explicit1e-12 conversion establish the actual unit. The successor's80% action/20% assumed leakage split must be recombined, not used as a20% discount. No coefficients were modified/refitted, and no full CiMLoop/Accelergy stack was installed.

The588-point grid sweeps4–10ENOB,1M through10G aggregate values/s plus actual Gamma-related2.5M/5M/42.5M/85M/160M/320M rates, and1/2/4/8/16/32/64ADCs. A20G diagnostic retains strict-model rejections. The within-(ENOB,aggregate-rate) frontier trades total area against energy; different precision/throughput contracts are not mixed into one misleading frontier. Physical ADC count changes per-ADC frequency and total area, never divides energy/conversion again.

The code's10GHz cap is hard-coded rather than a measured130nm feasibility contour. Near-resolution130nm support is sparse and flagged in the grid. The original dataset uses sample clock fs but energy P/fsnyq; oversampled internal-clock energy must not be called independent conversion energy. All rates, OSR and both energy definitions are retained. Model warnings, including infeasible *alternative* configurations printed by upstream, are separate from validity of the chosen row.

Twelve130nm entries are in the distributed215-row training CSV. All{source['all_raw_130nm_anchors']}130nm records from the pinned survey are exported separately. Some actual ADCs beat the fitted trend: for8ENOB we also use a compatible130nm8.345ENOB SAR anchor's18.4pJ/sample, rather than its9.2pJ/Nyquist-value number. It can meet the160M aggregate bank rate with at least4ADCs running at≤50MHz under optimistic idle gating. Its measured area, not fitted area, follows that anchor. Model and empirical-envelope scenarios remain separate. Individual anchors are not claimed as optimized CIM peripherals.

Sources: [published methodology](https://arxiv.org/html/2404.06553v2), [historical estimator](https://github.com/Accelergy-Project/accelergy-adc-plug-in), [maintained successor](https://github.com/Accelergy-Project/hwcomponents-adc), [CiMLoop](https://github.com/mit-emze/cimloop), [pinned survey](https://github.com/bmurmann/ADC-survey/tree/79a4b6b1cfd0a337371804b0898cfba2e8a495ac). Exact commits/file hashes/model version are in baseline_manifest.json and ../../baseline_adc/source_manifest.json.

## B — Work, precision and amortization

Primary comparison:64 output columns,2.5M dot results/s/column,160M conventional one-shot conversions/s bank aggregate. A dot has N useful signed INT8 MACs. Conventional cost is conversions×energy/conversion divided by useful MACs, with ADC sharing included through actual per-ADC rate. Gamma needs64×P parallel boundary positions to match that *boundary* cadence. Full-dot normalization throughput is not demonstrated.

The128/512/2048/8192 and6/7/8/9 pairings are motivated by the published architecture study; they are not exact accumulator-bit requirements. The strong one-shot conventional assumption here is full multilevel INT8 multiplication followed by b-ENOB quantization over the declared full-scale range. Input DAC/programming costs are omitted. A separate exact64-bitplane reference for signed INT8 uses enough chunks to keep each unsigned partial count in0..2^b−1; invalid converter-rate requirements are retained. No lossy6-bit dot is equated with a complete exact INT8 accumulator.

Gamma precision uses the existing exact finite-word/exponent oracle. For R=16384N and H=floor(log_beta R), sign-uniform input encodings bound both beta-weighted absolute coefficient masses by R; legal local rewrites cannot increase them. Hence nonzero positions satisfy|k|≤H. Truncating below L discards at most3 beta^L/(beta−1). Integer/rational algebraic pairs compare that bound with the ideal ADC quantization RMS; floating beta does not establish correctness. Each layout has545 exact finite-word checks, with a common block exponent and an explicit loss ledger. This is a representation/error bound, not physical ENOB or application inference accuracy.

{layouttable}

B1 compares peripherals only. B2 explicitly cancels equal analog-array energy. Primary conventional binary accumulation/requantization costs are set to zero to favor the conventional baseline;100/1000fJ combined costs are labeled sensitivities, not calibrated130nm logic measurements. B3 sweeps array energy expansion, physical positions, routing, missing logic, reference sharing, refresh rate and normalization prices/activity. Array expansion cost is not inferred from a digit-product count alone.

{case_table}

The Gamma boundary estimate includes the measured hybrid plus two measured primitive-comparator energies to model3 decisions. That assumes an unbuilt adaptive reference-selection/latch architecture;6-decision and1-real-decision lower-bound scenarios are exported too. These figures are not a fully measured seven-state detector. Decoder, completion, reference multiplexing, all added clocks and routing must fit inside the remaining allowance.

The normalization comparison retains197.25 overflow-comparator evaluations and32.3828 unit rewrites per product from the prior event model, multiplied by N. Any merge normalization is additional. This per-product schedule is one architecture, not an inherent Gamma lower bound. Shared/sparse schedules can change counts but must demonstrate bounded transient states and completion. Deferred normalization of an arbitrarily wide raw analog sum is not a free or validated alternative.

{normtable}

The old0.1fJ detector/1.2fJ rewrite/0.01fJ control numbers remain **uncalibrated** reference assumptions; only N=128 has a small modeled win at those prices with all other costs omitted. A1fJ overflow hypothesis loses all four primary scenarios under this schedule. Reusing the measured class comparator as every overflow detector is far worse; this is a conditional cost proxy, not proof that a coarse self-triggering normalizer must cost106fJ.

The complete physical ENOB/BER match, total-normalizer throughput, analog input/weight implementation and equal array cost are not established. Mathematical break-even regions in the CSV/plots therefore authorize a focused primitive experiment, not a system efficiency claim. High converter counts, low utilization and the exact-bitserial contract can favor Gamma arithmetically; they are labeled rather than presented as wins against the strongest one-conversion/output baseline.

## C — Retained transistor evidence and new checks

The actual SKY130 regeneration cell already existed at entry. Its351 sizing/PVT/mismatch/comparator/transition/control runs are reused with original seeds, literal decks and logs. No old report or q sweep was regenerated. After ADC calibration and baseline definition, only two capacitor-energy and four reference-loading circuits were added. The primitive-only PDK remains unchanged, with ordinary1.8V NMOS/PMOS and MIM caps; no ASIC flow or layout was installed. Device paths, commit, bytes and exact device names are in sky130_device_manifest.json.

The20 retained core/corner smokes pass. The new ideal100fF and MIM charge/discharge tests both deliver331.08fJ on the charging rail versus324.43fJ expected from C×VDD×deltaV, a2.05% difference including switch parasitics. Stored energy is162fJ; it is not mistaken for total supplied energy. These full-swing diagnostics do not replace or inflate the measured25fF small-swing Gamma boundary cost. All new capacitor smoke checks pass.

Q1, restoration: the equal-scale actual transistor path returns12/12 correct with+0.35 digit input residual; direct returns1/12, and both return12/12 at zero. All49 reused state transitions classify correctly nominally. The older S12→S4 direct control has an extra finite-gain common-mode loading problem and is not used as clean evidence for error-reset benefit.

Q2, energy: known hybrid cost is748.585fJ/position:720.398fJ signed1.8V supply work,727.865fJ gross delivered rail energy, plus20.720fJ external control-pin charging allowance. Reference generation is324.314fJ, loaded driver248.635fJ, one physical comparator99.715fJ, and interface55.200fJ. Five decisions and controller remain ideal. A3-decision estimate is{comp_est:.3f}fJ/position. Typical storage is25fF, selected from the existing25–200fF sizing grid, with the actual transistor buffer driving500kohm/1Mohm MAC resistors and20fF load.

The measured voltage still drives the non-unity resistor MAC; Python inserts no intermediate answer. Per-position classification uses the declared gamma regions. Position significance never changes physical amplitude. Exponents remain digital. No high-resolution ADC/general DAC appears in the Gamma computational path.

Reference sharing is explicitly tested as a lumped coherently switched35fF/load writer load (25fF storage plus10fF input allowance), with TG width and MIM capacitance scaled together. This is not64 complete buffered cells or a mismatch-qualified distribution network:

{ftable}

Reducing resistance to recover write settling raises static reference cost; a free1/F power division is not established for this passive topology. The F=8/64 sensitivity columns are hypothetical alternatives, not validated energy savings. Sparse write skipping also leaves continuous references, driver bias and detection alive: skipping every write saves only{skip_saving:.2f}fJ/position({100*skip_saving/comp_est:.1f}% of the3-decision estimate) in this accounting. Sparse **normalization activity** has a different and potentially larger effect, but requires a new bounded schedule.

Q3, competitiveness: readout-only opportunity exists at matched node/work/rate without unusually high ADC resolution or deliberately poor amortization. Total peripheral advantage is unresolved because normalizer/control and physical accuracy are missing. The former100fJ absolute rejection is withdrawn for this matched comparison; it does not remove the new quantitative cost ceilings.

Q4, dominant costs: references and output bias dominate the implemented hybrid; a modeled3-decision detector becomes comparably large. With the retained per-product schedule, normalization can dominate the entire dot-product cost. Representation expansion and extra routing/array energy can erase readout headroom, as the sensitivity surfaces show.

PVT/mismatch remains a separate failure gate: the low-energy hybrid passes all tested output classes at29/45 PVT points, and both classes plus monitored device limits at{joint}/45. Four larger-device runs are numerical nonconvergences, retained as unresolved. Core supply sweep remains1.62/1.8/1.9V;1.98V exceeds the documented1.95V model limit. Temperatures are−40/25/100C. Existing drift, reference kickback and driver mismatch failures are not repaired by recalibrating ADC energy.

{mctable}

These are independent circuit seeds; twelve correlated checkpoints inside each trial are not twelve yield samples. Bounds are one-sided exact95% binomial limits under the stated uncalibrated model. No1e−3/1e−6 foundry BER, thermal-noise certification, silicon measurement or layout signoff is claimed.

The reused CSV's max_settling_ns=1e9 is a censored optimization sentinel: the0.1-digit settling target was not reached inside observation. It is not a measured one-second delay. Explicit unsettled-cycle counts and actual observed settling fields are retained.

## Q5 — Next move

**Keep q=4 and continue with one local primitive; do not build a crossbar yet.** The useful next circuit is a bounded overflow decision plus actual4:1:1 transfer, with honest detector/packet/completion supply energy and event counts. Evaluate whether its price/activity can fit the0.398→0.016fJ/evaluation ceiling for the retained schedule, or demonstrate a cheaper bounded shared schedule. Fix canonical driver/reference mismatch in parallel with that gate. Buffered/shared or charge-oriented storage/reference topologies are valid experiments; passive free sharing and write skipping alone are insufficient evidence.

Proceed to the smallest resistor-array experiment only after both a positive *all-in* matched budget and useful PVT/mismatch margins are demonstrated. The current evidence supports a focused optimization question, neither a system-win claim nor a universal kill of the datatype.

## Outputs and reproduce

New outputs are confined to results/sky130_phase/ and baseline_adc/. Required tables/JSON are present, with additional precision golden traces, exact-bitserial comparisons, reference-fanout checks and sensitivity surfaces. Reused physical case directories remain under spice/sky130/; prior LTspice artifacts remain under spice/regen_results/ and spice/regen_runs/.

From python/: run gamma_adc_baseline.py, gamma_matched_130nm.py, sky130_matched_cap_smoke.py, sky130_reference_load.py, then report_matched_130nm.py. Use python -B to avoid bytecode writes into the pinned source trees. The last regression result is:{regression['output'].strip().splitlines()[-1]}. No new transistor optimization precedes baseline calibration.

![ADC model and anchors](plots/adc130_energy_and_anchors.png)
![ADC energy/area](plots/adc130_energy_area_frontier.png)
![Matched work](plots/matched_energy_per_mac.png)
![Normalizer ceiling](plots/normalizer_cost_ceiling.png)
![Normalization sensitivity](plots/normalization_break_even_surfaces.png)
![Array expansion](plots/array_expansion_sensitivity.png)
![Measured boundary blocks](plots/measured_boundary_components.png)
'''
    (OUT/'engineering_report.md').write_text(text,encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ('status','ADC','energy','recommendation')},indent=2),flush=True)


if __name__=='__main__':main()
