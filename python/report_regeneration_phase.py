"""Supplemental physical checks and a concise evidence-linked phase report."""
from __future__ import annotations
import csv, json, hashlib, statistics, random
import shutil,sys
from collections import defaultdict
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from gamma_regen_spice import *
from run_regeneration_phase import ROOT,OUT,RUNS,dump


def cached_supplement(deck,meta,name):
    directory=RUNS/name;cache=directory/'supplement.json'
    sha=hashlib.sha256(deck.encode()).hexdigest()
    if cache.exists():
        row=json.loads(cache.read_text())
        if row['deck_sha256']==sha:return row
    retain=name.startswith(('a','t')) or name in ('c00','c01','c20','c21','c60','c61')
    row=run(deck,meta,directory,retain);row['run_path']=str(directory.relative_to(ROOT))
    dump(cache,row);return row


def supplements():
    convergence=[]
    cases=[('nominal direct',chain_deck,dict(digit=-3,mode='none')),
        ('drift direct',chain_deck,dict(physical=RegenPhysical(stage_drift_units=.06),mode='none')),
        ('drift regenerated',chain_deck,dict(physical=RegenPhysical(stage_drift_units=.06))),
        ('stress regenerated',chain_deck,dict(digit=random.Random(70005).randrange(-3,4),physical=profile('stress'),seed=50005)),
        ('stress direct',chain_deck,dict(digit=random.Random(70005).randrange(-3,4),physical=profile('stress'),seed=50005,mode='none')),
        ('word regenerated',word_deck,dict(physical=RegenPhysical(stage_drift_units=.35))),
        ('word direct',word_deck,dict(physical=RegenPhysical(stage_drift_units=.35),mode='none')),
        ('feedthrough initialization',primitive_deck,dict(digit=-3,mapping='single',physical=profile('low'),seed=20133))]
    for i,(name,fn,kwargs) in enumerate(cases):
        results=[]
        for j,step in enumerate((.04,.02)):
            kw=dict(kwargs);kw['physical']=replace(kw.get('physical',RegenPhysical()),maxstep_ns=step)
            results.append(cached_supplement(*fn(**kw),f'c{i}{j}'))
        same=all(a['decoded_class']==b['decoded_class'] and a['overflow_latched']==b['overflow_latched']
            for a,b in zip(results[0]['observations'],results[1]['observations']))
        convergence.append({'case':name,'all_checkpoint_classes_identical':same,'max_continuous_difference_units':
            max(abs(a['residual_units']-b['residual_units']) for a,b in zip(results[0]['observations'],results[1]['observations'])),
            'paths':[r['run_path'] for r in results]})
        print('convergence',name,same,flush=True)
    reuse=[]
    for m,mapping in enumerate(('single','differential','charge')):
        for d in range(-3,4):
            for initial in (-3,3):
                deck,meta=primitive_deck(d,mapping)
                for j in range(1 if mapping=='single' else 2):
                    factor=1 if mapping=='single' else (.5 if j==0 else -.5)
                    old=next(line for line in deck.splitlines() if line.startswith(f'Cro{j} '))
                    deck=deck.replace(old,old.replace('IC=0.6',f'IC={.6+factor*initial*.025:.14g}'))
                meta['initial_stored_digit']=initial
                r=cached_supplement(deck,meta,f'u{m}{d+3}{0 if initial<0 else 1}')
                reuse.append({'mapping':mapping,'digit':d,'initial_digit':initial,
                    'correct':r['observations'][0]['correct'],'residual_units':r['observations'][0]['residual_units'],
                    'reference_ports_fj':r['reference_port_delivered_energy_fj'],
                    'path':r['run_path']})
    adverse=[]
    for i,(label,p) in enumerate([('insufficient write settling',RegenPhysical(ron_ohm=20000)),
            ('15pct reference/supply gain',RegenPhysical(supply_gain=.15))]):
        r=cached_supplement(*chain_deck(digit=3,physical=p),f'a{i}')
        adverse.append({'mechanism':label,'run_path':r['run_path'],
            'first_output_class':r['observations'][0]['decoded_class'],
            'first_output_residual_units':r['observations'][0]['residual_units'],
            'first_wrong_stage':next((o['stage'] for o in r['observations'] if not o['correct']),None),
            'final_correct':r['observations'][-1]['correct']})
    diagnosis=[]
    for trial in range(4):
        digit=random.Random(70000+trial).randrange(-3,4)
        deck,meta=chain_deck(digit=digit,physical=profile('stress'),seed=50000+trial)
        r=cached_supplement(deck,meta,f't{trial}')
        directory=ROOT/r['run_path'];names,x=fast_raw(directory/'cell.raw')
        trace=[]
        for obs in r['observations']:
            i=np.searchsorted(x[:,0],obs['time_ns']*1e-9);tag=obs['bad'][:-2]
            code=float(x[i,names.index('v('+tag+'h)')])
            labels=[d for d in range(-3,4) if abs(code-d)<.05]
            label=labels[0] if len(labels)==1 else 'UNRESOLVED_CONTROL'
            rawclass=classify(obs['raw_v'],meta['rawunit'],3)
            failure=('correct' if obs['correct'] else 'raw MAC/input outside intended region' if rawclass!=digit
                else 'detector/overflow decision' if label!=digit or obs['overflow_latched']
                else 'canonical write/reference/noise outside intended region')
            trace.append({'stage':obs['stage'],'nominal_raw_class':rawclass,'latched_region_label':label,
                'output_class':obs['decoded_class'],'overflow_latched':obs['overflow_latched'],
                'output_residual_units':obs['residual_units'],'failure_category':failure})
        diagnosis.append({'seed':meta['seed'],'digit':digit,'path':r['run_path'],'trace':trace,
            'first_failure':next((v for v in trace if v['failure_category']!='correct'),None)})
    result={'convergence':convergence,'reused_storage':reuse,'adverse_chains':adverse,
        'instrumented_stress_diagnosis':diagnosis,
        'additional_actual_spice_runs':2*len(cases)+len(reuse)+len(adverse)+len(diagnosis)}
    dump(OUT/'supplemental_checks.json',result)
    return result


def examples_and_csv():
    destination=ROOT/'spice/behavioral'
    for name,fn,kw in [('single',primitive_deck,dict(mapping='single')),
        ('differential',primitive_deck,dict(mapping='differential')),
        ('charge',primitive_deck,dict(mapping='charge')),
        ('two_stage',word_deck,dict(physical=RegenPhysical(stage_drift_units=.35))),
        ('depth16',chain_deck,{}),('boundary',boundary_deck,{})]:
        deck,meta=fn(**kw)
        (destination/f'q4_regen_{name}.cir').write_text(deck,encoding='ascii')
        dump(destination/f'q4_regen_{name}.golden.json',meta)
    # Explicit exact JSON is authoritative; CSV is a small phase/checkpoint index.
    records=[]
    snapshot=OUT/'source_snapshot';snapshot.mkdir(exist_ok=True)
    sources=('gamma_regen_golden.py','gamma_regen_spice.py','run_regeneration_phase.py','report_regeneration_phase.py',
        'gamma_native_model.py','gamma_protocols.py','gamma_spice.py','gamma_normalizer_statistics.py',
        'test_gamma_regeneration.py','requirements.txt','conftest.py')
    for name in sources:shutil.copyfile(ROOT/'python'/name,snapshot/name)
    dump(OUT/'reproduction.json',{'python':sys.version,'simulator_executable':str(DEFAULT_EXE),
        'source_sha256':{name:hashlib.sha256((snapshot/name).read_bytes()).hexdigest() for name in sources}})
    manifest=json.loads((OUT/'manifest.json').read_text())
    for i in range(manifest['jobs']):
        directory=RUNS/f'{"v" if i<86 else "w"}{i:04d}'
        meta=json.loads((directory/'golden.json').read_text())
        with (directory/'checkpoints.csv').open('w',newline='',encoding='ascii') as f:
            writer=csv.writer(f);writer.writerow(['time_ns','stage','node','expected_digit','physical_output'])
            for c in meta['checkpoints']:writer.writerow([c['time_ns'],c.get('stage',1),c.get('node',0),c['expected'],c['out']])
        records.append(json.loads((directory/'cached.json').read_text()))
    return records


def word_plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(9,7))
    for index,mode in ((40,'direct analog'),(72,'regenerated')):
        directory=RUNS/f'v{index:04d}'
        names,x=fast_raw(directory/'cell.raw');t=x[:,0]*1e9
        for ax,node,u in zip(axes,('ap0','r3out','bp0'),(.025/3,.025,.025/3)):
            ax.plot(t,x[:,names.index('v('+node+')')]/u,label=mode)
            ax.axhline(1,color='black',ls=':',lw=1);ax.grid(alpha=.2)
    axes[0].set(xlim=(0,145),ylim=(-.2,5.8),ylabel='Stage 1 raw\ncoefficient units')
    axes[1].set(xlim=(132,148),ylim=(-.1,1.6),ylabel='Restored interface\ndigit units')
    axes[2].set(xlim=(140,280),ylim=(-.2,2.2),ylabel='Stage 2 raw\ncoefficient units',xlabel='Time (ns)')
    axes[2].axhline(1.5,color='red',ls='--',lw=1,label='class 1/2 boundary')
    axes[0].legend();axes[2].legend(fontsize=8)
    fig.suptitle('Connected weighted MACs [2,1] then [2,-1]; injected residual +0.35')
    fig.tight_layout();fig.savefig(OUT/'two_stage_restoration.png',dpi=170);plt.close(fig)


def report(rows,supplement):
    phase=json.loads((OUT/'phase_summary.json').read_text());summary=phase['summary'];completion=phase['completion'];energy=phase['energy']
    def lines_table(headers,lines):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in lines])
    chain=[r for r in summary if r['group']=='chain_mc']
    depth_table=lines_table(['Profile / slope','Mode','Correct at depths 2 / 4 / 8 / 16','Depth16 overflow','Median absolute residual at depth16'],
        [(f'{p} / {s}',m,' / '.join(f"{next(r for r in chain if r['profile']==p and r['slope']==s and r['mode']==m and r['depth']==d)['p_correct']:.1%}" for d in (2,4,8,16)),
          next(r for r in chain if r['profile']==p and r['slope']==s and r['mode']==m and r['depth']==16)['unresolved_overflow'],
          f"{next(r for r in chain if r['profile']==p and r['slope']==s and r['mode']==m and r['depth']==16)['abs_residual_units_median']:.4g}")
            for p in ('low','stress') for s in ('2','1/2') for m in ('none','regen')])
    primitive_table=[]
    level_measurements=[]
    for d in range(-3,4):
        record=[d]
        for mapping in ('single','differential','charge'):
            group=[r for r in rows if r['group']=='primitive_mc' and r['labels']['digit']==d and r['labels']['mapping']==mapping and r['labels']['profile']=='low']
            record.append(f"{sum(not r['observations'][0]['correct'] for r in group)}/{len(group)}; {max(abs(r['observations'][0]['residual_units']) for r in group):.3g}")
        primitive_table.append(record)
        measure=[d]
        for mapping in ('single','differential','charge'):
            r=next(r for r in rows if r['group']=='primitive_ideal' and r['labels']['digit']==d and r['labels']['mapping']==mapping)
            measure.append(f"{r['reference_port_delivered_energy_fj']:.3g}; {r['first_output_settling_ns_to_1pct_class_spacing']:.3g}")
        level_measurements.append(measure)
    sensitivities=[]
    for mechanism in ('reference','switch','capacitance','feedthrough','charge_injection','comparator','supply','thermal'):
        record=[mechanism]
        for level in (0,1):
            group=[r for r in rows if r['group']=='sensitivity' and r['labels']['mechanism']==mechanism and r['labels']['level']==level]
            record.append(f"{sum(not r['observations'][0]['correct'] for r in group)}/{len(group)}; {max(abs(r['observations'][0]['residual_units']) for r in group):.3g}")
        sensitivities.append(record)
    ideal=[r for r in rows if r['group'] in ('primitive_ideal','chain_ideal','boundary') or (r['group']=='word' and r['labels']['drift']==0)]
    ideal_wrong=sum(any(not o['correct'] for o in r['observations']) for r in ideal)
    word=[r for r in rows if r['group']=='word'];g=[r for r in word if r['labels']['mode']=='regen' and r['labels']['drift']==.35]
    direct=[r for r in word if r['labels']['mode']=='none' and r['labels']['drift']==.35]
    maxregen=max(r['stage1_canonical_residual_units_max'] for r in g)
    stressful=next(r for r in chain if r['profile']=='stress' and r['slope']=='2' and r['mode']=='regen' and r['depth']==16)
    adc100=(125-energy['reference_point']['explicit event model']['before_regen_extra_fj'])/17
    confidence_upper=1-.05**(1/phase['trial_count'])
    regression=OUT/'regression_result.json'
    regression_note=(json.loads(regression.read_text())['output'].strip().splitlines()[-1]
        if regression.exists() else 'Regression results are recorded separately after the full suite finishes.')
    scale_table=lines_table(['S','Correct final class at depth16','Circuits'],[(r['scale'],f"{r['circuits']-r['wrong_final']}/{r['circuits']}",r['circuits'])
        for r in summary if r['group']=='scale' and r['depth']==16])
    adverse_table=lines_table(['Adverse condition','First output class','First residual (digits)','First wrong stage'],
        [(r['mechanism'],r['first_output_class'],f"{r['first_output_residual_units']:.3f}",r['first_wrong_stage']) for r in supplement['adverse_chains']])
    text=f'''# q=4 canonical regeneration: engineering results

Run: 2026-10-23, installed LTspice XVII. Lead beta=2+sqrt(3), D=-3..3.

## Finding

**Equivalence sensing can act as a regenerative boundary in the tested behavioral circuits. Physical cost and process robustness are still open.** An actual seven-reference write restores a capacitor, whose buffered voltage drives the next non-unity-weight analog MAC. No intermediate expected answer is inserted by Python. No general high-resolution ADC/DAC is in the path.

All {len(ideal)} ideal main circuits pass every declared checkpoint ({ideal_wrong} failures), including all seven digits, three storage mappings, sixteen connected stages and the complete signed word MAC. With a stage-1 residual of +0.35 coefficient units, the regenerated word path passes {len(g)}/{len(g)} cases; the direct path passes {sum(r['final_algebraically_correct'] for r in direct)}/{len(direct)}. Maximum restored stage-1 residual is {maxregen:.3g} digit units. This near-zero ideal number is a behavioral/numerical result, not a reference-accuracy specification.

**The higher-stress expanding regenerated chains have {stressful['wrong_final']}/{stressful['circuits']} incorrect final results at depth16.** Restoration is therefore conditional on analog margins; it does not make this stress profile reliable. Instrumented first-failure traces distinguish raw-MAC/input errors from detector decisions and canonical-write errors.

There are {phase['actual_spice_circuits']} main circuit runs and {supplement['additional_actual_spice_runs']} supplemental runs. Seeds, decks, exact golden JSON, checkpoint CSV, simulator logs and class-decoded results are retained. Selected waveforms are retained; other waveforms are reproducible from their netlists. Previous local-normalizer/base comparisons remain unchanged.

![Connected word restoration](../spice/regen_results/two_stage_restoration.png)

## What is connected

External canonical a,b -> conductances [2,1]/3 -> five active q=4 normalizer nodes with two passive guards -> finite-region decision -> real reset/write capacitors -> ideal output buffer -> conductances [2,-1]/3 -> second local normalizer -> final finite regions. Raw spacing is 25/3 mV; restored interface spacing is 25 mV. Physical weights alone are normalized. Neighbor positions use identical packet sizes and identical local levels.

The seven-reference selector, comparator and buffer are behavioral. The class-hold scalar -3..3 is a simulator control label, not a proposed +/-3 V logic supply. Real one-hot/binary latches, reference generation, clocks and output buffers must be implemented and priced. The word circuit uses clocked 32-phase normalization. Event completion is modeled separately.

## Depth and robustness

Each Monte Carlo row below contains {phase['trial_count']} independent seeded whole circuits, with one sample at each depth. Paired variants share seeds; checkpoints within a circuit are correlated. Digits are sampled from a fixed uniform seven-digit workload. The affine MAC y=2x-d amplifies perturbations; y=(x+d)/2 is a contracting control. The external initial d is a fixed second input, not oracle feedback.

{depth_table}

![Depth robustness](../spice/regen_results/chain_depth.png)

Continuous residual is relative to the intended digit and includes clean wrong-digit errors. Direct outputs rail at +/-0.3 V; beyond overflow, reported residuals must not be read as a valid wider numerical representation. Scalar depth tests flag unresolved overflow rather than silently normalizing into an unmodeled second position. Equivalent-but-different scalar words are structurally impossible; that frequency is zero by construction. The complete multi-position word test separately checks exact A+B beta equivalence and records equivalent-word differences.

Low/stress profiles are **synthetic brackets, not foundry distributions**. Low: reference sigma0.012 digit, comparator sigma0.015 raw digit, bank sigma0.01, C sigma1%, MAC gain sigma0.2%, sampled MAC/read disturbance sigma0.015, 300 K kTC kick, charge injection sigma1 aC, clock feedthrough0.01 fF, leakage time1e6 ns. Stress: reference/comparator sigma0.12, bank0.08, C6%, MAC gain3%, MAC/read disturbance0.12, injection100 aC, feedthrough0.2 fF, leakage2000 ns. The sampled disturbance is input-referred at the detector; it is not a continuous current-noise spectrum. Packet/width mismatch also affects the charge candidate. Supply variation is modeled as reference gain, not a full transistor supply perturbation. Threshold-bank offsets are correlated within a cell; random reference mismatch is independent between cells/stages. Persistent shared-bank drift is only bracketed by the deterministic gain adverse case. Profiles and realizations are in every golden record.

Zero failures in {phase['trial_count']} independent circuits would only give a pointwise one-sided 95% error upper bound of **{confidence_upper:.1%}** for that workload/profile. No 1e-3/1e-6 reliability claim follows. Per-class primitive groups are smaller still. All actual failures, overflows and bounds are retained in phase_summary.json and summary.csv.

The supplementary common-scale sweep holds absolute reference/comparator/read-noise volts fixed rather than improving noise as spacing shrinks. It uses only four expanding regenerated chains per S; it does not certify Smax:

{scale_table}

## Canonical levels and nonidealities

The [circuit specification](../spice/canonical_regeneration_spec.md) defines every level and threshold for single-ended voltage, differential voltage and differential charge storage at equal total100 fF. Voltage-reference error is additive per level; charge-reference error perturbs unit-packet gain and grows with |d|. Thus these profiles compare mechanism brackets, not equal-noise/equal-area reference implementations. The following low-profile table gives wrong-class count / trials; maximum |nominal canonical residual| in digit units for every class:

{lines_table(['Digit','Single-ended','Differential voltage','Differential charge'],primitive_table)}

For each ideal class, the next table gives partial reference-port delivered fJ; achieved-level settling ns from decision start. Charge steering can have a near-zero port integral because ideal signed current sources cancel/recover supply current; this is **not** zero cell energy or a fair transistor-energy ranking.

{lines_table(['Digit','Single-ended port fJ; ns','Differential voltage port fJ; ns','Differential charge port fJ; ns'],level_measurements)}

Isolated sensitivity brackets use digits -3,0,+3 and two independent realizations each. Entries are wrong-class count / trials; maximum nominal residual. These are sensitivity examples, not rare-event estimates:

{lines_table(['Mechanism','Lower bracket','Higher bracket'],sensitivities)}

All seven classes retain individual settling and reference-port energy measurements in the run results. Settling is measured within1% of one spacing of the achieved final plateau, separately from nominal-level error; a wrong reference does not become accurate by settling quickly. Isolated common-mode feedthrough cancels in matched differential capacitors; combined mismatch breaks that cancellation. Artificial 30,000 K samples are an exaggerated thermal sensitivity bracket, not an operating recommendation.

There is a stronger next-stage requirement than merely remaining in the present class:

\\[
\\left|\\sum_i w_i\\epsilon_i+\\eta_{{MAC}}-\\tau_{{sense}}\\right|<\\tfrac12.
\\]

For the 2x following weight, a single input residual must be below0.25 digit units before reserving margin for other errors. A reference can therefore regenerate into the correct current class and still corrupt the following MAC. Supplemental adverse chains demonstrate this for insufficient write settling and reference/supply gain error; see supplemental_checks.json. Regeneration cannot recover a wrongly selected class. Correlated reference errors must not be treated as independent noise that averages away.

{adverse_table}

Both adverse first outputs are classified correctly, yet stage2 is already wrong. The first four instrumented stress seeds additionally show two first failures from raw-MAC/input perturbation outside the intended region and two from detector decisions while the nominal raw class was still correct. These diagnoses are not population fractions. The limiting requirement is the following weighted MAC's half-class error budget. Better references, lower S, larger capacitors or calibration may restore margin, but their energy/area cost remains unproved.

## Completion and finite boundaries

On the same128 frozen q=4 INT8 product workloads, using a compiled17-position layout:

{lines_table(['Schedule','Comparator decisions/product','Latency ns (model)','Rewrite events/product'],[
('Fixed32',f"{completion['fixed32_comparator_evaluations']:.2f}",128,f"{completion['synchronous_rewrites_mean']:.2f}"),
('Ideal completion lower bound',f"{completion['ideal_until_settled_comparator_evaluations']:.2f}",f"{completion['ideal_until_settled_time_ns']:.2f}",f"{completion['synchronous_rewrites_mean']:.2f}"),
('Explicit event queue',f"{completion['explicit_event_comparator_evaluations']:.2f}",f"{completion['event_time_ns']:.2f}",f"{completion['unit_rewrites_mean']:.2f}")])}

The event controller charges initial scan, six persistent plus two guard comparisons per final node, dirty-node visits, request/ack, pending flags, commit revalidation, queue insert/pop, activity-tree toggles and final guard/root confirmation. No fixed32 workload exceeds the budget. The queue preserves the exact value and reaches quiescence. Its assumed .2 ns comparisons,1 ns transfers,.2 ns local settling and tree delays are timing assumptions; the latency advantage is not measured asynchronous hardware speed. Queue/tree energy can outweigh the comparator savings versus ideal free completion.

The exact finite-word policy uses [-2,2], declared guards, high-side index shift plus exponent update, optional exponent decrement, explicit low-side truncation and exact nearest-quantum ties-to-even rounding. Every lost value has an exact pair ledger. The physical boundary circuit routes the high digit and explicitly receives discarded low charge. Rounding, exponent comparison/addition and aligned residual/tile merges stay digital models; no automatic exponent-CIM or hardware rounding controller is claimed. Hardware exponent width/overflow trap remains a design choice.

## Energy gate

\\[
E_\\gamma=E_{{CIM}}+E_{{detect}}+E_{{rewrite}}+E_{{regen}}+E_{{completion}}+E_{{routing}}+E_{{scale/exponent}}.
\\]

The matched baseline includes CIM+ADC+binary accumulation+requantization. The sweep includes ADC amortization1/64,1/16,1,4,16 and varies detector, packet, controller and all-in regeneration costs. The all-in regeneration parameter includes reference generation, reset/write, local clock and next-stage buffer. Differential charge's ideal signed current sources can recover/cancel supply current; its reference-port integral is particularly incomplete. Port measurements cannot rank transistor energy across the three mappings.

At the illustrative assumptions comparator0.1 fJ, rewrite(q3 reference)1 fJ scaled6/5, controller activity0.01 fJ, binary accumulation20 fJ and requantization5 fJ, the explicit-event extra cost before regeneration is {energy['reference_point']['explicit event model']['before_regen_extra_fj']:.2f} fJ/product. With1 fJ/position regeneration, ADC break-even is {energy['reference_point']['explicit event model']['adc_break_even_at_1fj_regen']:.2f} fJ/conversion for one conversion/output. At100 fJ per ADC conversion, the allowed all-in regeneration budget is only **{adc100:.2f} fJ/position**. Greater completion, buffering, routing or array overhead reduces it. Strong ADC amortization can eliminate the budget even before regeneration.

![Energy break-even](../spice/regen_results/energy_break_even.png)

![Regeneration and completion budget](../spice/regen_results/energy_regen_completion.png)

These assumed costs are parameter axes, not measured realistic cell costs. Common CIM cancels only when explicitly equal; the optimistic zero delta-array/routing/sidecar reference point is stated in energy_sweep.json. Reference-port energy is partial, not total supply energy or proof that an ADC has been beaten. Reused-storage reset/write checks cover both previous extreme digits; complete repetitive throughput and reference-bank fanout remain open.

## Verification, limitations and next decision

Half-step convergence checks cover nominal/direct, drift, stress, complete word and feedthrough cases. Every reported checkpoint class is compared; exact results and any differences are retained in supplemental_checks.json. Startup feedthrough IC consistency and delayed reference selection were corrected before final results. Tiny parasitic/control capacitors caused trapezoidal integration chatter in a diagnostic nonideal case; nonideal configurations use Gear integration with timestep checks. Method selection is predetermined by configuration, not by correctness. Numerical diagnostic attempts are excluded from the experiment manifest and are not counted as successful physical trials.

Full regression: {regression_note}. The log and machine-readable exit status are in spice/regen_results/regression.log and regression_result.json. All eight half-step comparisons preserve every checkpoint class; all42 reused-storage resets/writes remain in the intended class.

Keep **q=4** as requested. Differential reference-selected voltage is the behavioral reference implementation; single-ended and charge alternatives remain controls, not eliminated hardware candidates. The restoration principle passes this behavioral gate. Stress failures establish that regeneration is conditional on the next MAC's margin, not a cure for bad references or wrong decisions. Physical energy superiority, calibrated class-error rates and asynchronous control are unestablished.

The next useful experiment is a loaded seven-state reference/write/driver cell feeding a tiny resistor MAC, with actual supply accounting and reference/comparator mismatch. Advance to a small resistor-crossbar prototype only if that cell fits the disclosed regeneration+completion budget. No FPGA or full crossbar is warranted yet.

## Reproduction and deliverables

From python/: `python -B run_regeneration_phase.py --trials16 --workers1` (write the options with spaces: `--trials 16 --workers 1`), then `python -B report_regeneration_phase.py`, then `python -B -m pytest -q -p no:cacheprovider`.

- Exact oracles, exponents, finite words, completion: gamma_regen_golden.py.
- Physical cell and connected generators: gamma_regen_spice.py; generated q4_regen_*.cir and .golden.json in spice/behavioral/.
- Main run manifest, source hashes and profiles: spice/regen_results/manifest.json.
- Seeded netlists, logs, exact golden checkpoints/CSV and decoded results: spice/regen_runs/v0000..v0085 and w0086..w{phase['actual_spice_circuits']-1:04d}, plus named supplemental runs. Other diagnostic/earlier-attempt directories are excluded by the manifest.
- Full statistics and bounds: phase_summary.json, summary.csv; exact control traces: completion_traces.json; finite boundary oracles: finite_word_golden.json.
- Supplemental convergence, reused-storage and adverse cases: supplemental_checks.json.
- Parameterized costs: energy_sweep.json; plotted results alongside this report under spice/regen_results/.

All paths above are relative to the project except the Python module names, which live in python/. Prior reports are historical and remain available.
'''
    text=text.replace('`python -B run_regeneration_phase.py --trials16 --workers1` (write the options with spaces: `--trials 16 --workers 1`)','`python -B run_regeneration_phase.py --trials 16 --workers 1`')
    (ROOT/'docs/regeneration_phase_results.md').write_text(text,encoding='utf-8')


if __name__=='__main__':
    rows=examples_and_csv();s=supplements();word_plot();report(rows,s)
    print('Regeneration report, examples, checkpoints and supplemental checks saved.',flush=True)
