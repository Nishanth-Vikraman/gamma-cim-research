# q=4 canonical regeneration: engineering results

Run: 2026-10-23, installed LTspice XVII. Lead beta=2+sqrt(3), D=-3..3.

## Finding

**Equivalence sensing can act as a regenerative boundary in the tested behavioral circuits. Physical cost and process robustness are still open.** An actual seven-reference write restores a capacitor, whose buffered voltage drives the next non-unity-weight analog MAC. No intermediate expected answer is inserted by Python. No general high-resolution ADC/DAC is in the path.

All 66 ideal main circuits pass every declared checkpoint (0 failures), including all seven digits, three storage mappings, sixteen connected stages and the complete signed word MAC. With a stage-1 residual of +0.35 coefficient units, the regenerated word path passes 8/8 cases; the direct path passes 0/8. Maximum restored stage-1 residual is 2.21e-07 digit units. This near-zero ideal number is a behavioral/numerical result, not a reference-accuracy specification.

**The higher-stress expanding regenerated chains have 16/16 incorrect final results at depth16.** Restoration is therefore conditional on analog margins; it does not make this stress profile reliable. Instrumented first-failure traces distinguish raw-MAC/input errors from detector decisions and canonical-write errors.

There are 494 main circuit runs and 64 supplemental runs. Seeds, decks, exact golden JSON, checkpoint CSV, simulator logs and class-decoded results are retained. Selected waveforms are retained; other waveforms are reproducible from their netlists. Previous local-normalizer/base comparisons remain unchanged.

![Connected word restoration](../spice/regen_results/two_stage_restoration.png)

## What is connected

External canonical a,b -> conductances [2,1]/3 -> five active q=4 normalizer nodes with two passive guards -> finite-region decision -> real reset/write capacitors -> ideal output buffer -> conductances [2,-1]/3 -> second local normalizer -> final finite regions. Raw spacing is 25/3 mV; restored interface spacing is 25 mV. Physical weights alone are normalized. Neighbor positions use identical packet sizes and identical local levels.

The seven-reference selector, comparator and buffer are behavioral. The class-hold scalar -3..3 is a simulator control label, not a proposed +/-3 V logic supply. Real one-hot/binary latches, reference generation, clocks and output buffers must be implemented and priced. The word circuit uses clocked 32-phase normalization. Event completion is modeled separately.

## Depth and robustness

Each Monte Carlo row below contains 16 independent seeded whole circuits, with one sample at each depth. Paired variants share seeds; checkpoints within a circuit are correlated. Digits are sampled from a fixed uniform seven-digit workload. The affine MAC y=2x-d amplifies perturbations; y=(x+d)/2 is a contracting control. The external initial d is a fixed second input, not oracle feedback.

| Profile / slope | Mode | Correct at depths 2 / 4 / 8 / 16 | Depth16 overflow | Median absolute residual at depth16 |
| --- | --- | --- | --- | --- |
| low / 2 | none | 100.0% / 100.0% / 25.0% / 0.0% | 16 | 12.98 |
| low / 2 | regen | 100.0% / 100.0% / 100.0% / 100.0% | 0 | 0.01463 |
| low / 1/2 | none | 100.0% / 100.0% / 100.0% / 100.0% | 0 | 0.01895 |
| low / 1/2 | regen | 100.0% / 100.0% / 100.0% / 100.0% | 0 | 0.01463 |
| stress / 2 | none | 81.2% / 18.8% / 0.0% / 0.0% | 16 | 11.09 |
| stress / 2 | regen | 81.2% / 31.2% / 12.5% / 0.0% | 16 | 1.934 |
| stress / 1/2 | none | 93.8% / 100.0% / 100.0% / 81.2% | 3 | 0.1152 |
| stress / 1/2 | regen | 87.5% / 93.8% / 100.0% / 81.2% | 2 | 0.1358 |

![Depth robustness](../spice/regen_results/chain_depth.png)

Continuous residual is relative to the intended digit and includes clean wrong-digit errors. Direct outputs rail at +/-0.3 V; beyond overflow, reported residuals must not be read as a valid wider numerical representation. Scalar depth tests flag unresolved overflow rather than silently normalizing into an unmodeled second position. Equivalent-but-different scalar words are structurally impossible; that frequency is zero by construction. The complete multi-position word test separately checks exact A+B beta equivalence and records equivalent-word differences.

Low/stress profiles are **synthetic brackets, not foundry distributions**. Low: reference sigma0.012 digit, comparator sigma0.015 raw digit, bank sigma0.01, C sigma1%, MAC gain sigma0.2%, sampled MAC/read disturbance sigma0.015, 300 K kTC kick, charge injection sigma1 aC, clock feedthrough0.01 fF, leakage time1e6 ns. Stress: reference/comparator sigma0.12, bank0.08, C6%, MAC gain3%, MAC/read disturbance0.12, injection100 aC, feedthrough0.2 fF, leakage2000 ns. The sampled disturbance is input-referred at the detector; it is not a continuous current-noise spectrum. Packet/width mismatch also affects the charge candidate. Supply variation is modeled as reference gain, not a full transistor supply perturbation. Threshold-bank offsets are correlated within a cell; random reference mismatch is independent between cells/stages. Persistent shared-bank drift is only bracketed by the deterministic gain adverse case. Profiles and realizations are in every golden record.

Zero failures in 16 independent circuits would only give a pointwise one-sided 95% error upper bound of **17.1%** for that workload/profile. No 1e-3/1e-6 reliability claim follows. Per-class primitive groups are smaller still. All actual failures, overflows and bounds are retained in phase_summary.json and summary.csv.

The supplementary common-scale sweep holds absolute reference/comparator/read-noise volts fixed rather than improving noise as spacing shrinks. It uses only four expanding regenerated chains per S; it does not certify Smax:

| S | Correct final class at depth16 | Circuits |
| --- | --- | --- |
| 4 | 4/4 | 4 |
| 8 | 4/4 | 4 |
| 16 | 4/4 | 4 |
| 32 | 0/4 | 4 |

## Canonical levels and nonidealities

The [circuit specification](../spice/canonical_regeneration_spec.md) defines every level and threshold for single-ended voltage, differential voltage and differential charge storage at equal total100 fF. Voltage-reference error is additive per level; charge-reference error perturbs unit-packet gain and grows with |d|. Thus these profiles compare mechanism brackets, not equal-noise/equal-area reference implementations. The following low-profile table gives wrong-class count / trials; maximum |nominal canonical residual| in digit units for every class:

| Digit | Single-ended | Differential voltage | Differential charge |
| --- | --- | --- | --- |
| -3 | 0/4; 0.0279 | 0/4; 0.0406 | 0/4; 0.039 |
| -2 | 0/4; 0.0126 | 0/4; 0.00887 | 0/4; 0.044 |
| -1 | 0/4; 0.0163 | 0/4; 0.0256 | 0/4; 0.0225 |
| 0 | 0/4; 0.0257 | 0/4; 0.0396 | 0/4; 0.0285 |
| 1 | 0/4; 0.0313 | 0/4; 0.0234 | 0/4; 0.0294 |
| 2 | 0/4; 0.0354 | 0/4; 0.0355 | 0/4; 0.0305 |
| 3 | 0/4; 0.0221 | 0/4; 0.0353 | 0/4; 0.0385 |

For each ideal class, the next table gives partial reference-port delivered fJ; achieved-level settling ns from decision start. Charge steering can have a near-zero port integral because ideal signed current sources cancel/recover supply current; this is **not** zero cell energy or a fair transistor-energy ranking.

| Digit | Single-ended port fJ; ns | Differential voltage port fJ; ns | Differential charge port fJ; ns |
| --- | --- | --- | --- |
| -3 | 2.65e-07; 1.78 | 1.2; 1.5 | 3.05e-25; 3.19 |
| -2 | 1.77e-07; 1.74 | 0.781; 1.48 | 3.05e-25; 3.19 |
| -1 | 8.84e-08; 1.68 | 0.383; 1.44 | 3.05e-25; 3.19 |
| 0 | 2.12e-14; 0 | 3.85e-14; 0 | 3.05e-25; 0 |
| 1 | 1.56; 1.67 | 0.383; 1.44 | 3.05e-25; 3.19 |
| 2 | 3.25; 1.75 | 0.781; 1.48 | 3.05e-25; 3.19 |
| 3 | 5.06; 1.78 | 1.2; 1.5 | 3.05e-25; 3.19 |

Isolated sensitivity brackets use digits -3,0,+3 and two independent realizations each. Entries are wrong-class count / trials; maximum nominal residual. These are sensitivity examples, not rare-event estimates:

| Mechanism | Lower bracket | Higher bracket |
| --- | --- | --- |
| reference | 0/6; 0.0111 | 0/6; 0.139 |
| switch | 0/6; 2.21e-07 | 0/6; 0.41 |
| capacitance | 0/6; 2.24e-07 | 0/6; 2.56e-07 |
| feedthrough | 0/6; 2.21e-07 | 0/6; 2.13e-07 |
| charge_injection | 0/6; 0.00473 | 2/6; 1.18 |
| comparator | 0/6; 2.21e-07 | 2/6; 3 |
| supply | 0/6; 0.03 | 0/6; 0.45 |
| thermal | 0/6; 0.0155 | 0/6; 0.155 |

All seven classes retain individual settling and reference-port energy measurements in the run results. Settling is measured within1% of one spacing of the achieved final plateau, separately from nominal-level error; a wrong reference does not become accurate by settling quickly. Isolated common-mode feedthrough cancels in matched differential capacitors; combined mismatch breaks that cancellation. Artificial 30,000 K samples are an exaggerated thermal sensitivity bracket, not an operating recommendation.

There is a stronger next-stage requirement than merely remaining in the present class:

\[
\left|\sum_i w_i\epsilon_i+\eta_{MAC}-\tau_{sense}\right|<\tfrac12.
\]

For the 2x following weight, a single input residual must be below0.25 digit units before reserving margin for other errors. A reference can therefore regenerate into the correct current class and still corrupt the following MAC. Supplemental adverse chains demonstrate this for insufficient write settling and reference/supply gain error; see supplemental_checks.json. Regeneration cannot recover a wrongly selected class. Correlated reference errors must not be treated as independent noise that averages away.

| Adverse condition | First output class | First residual (digits) | First wrong stage |
| --- | --- | --- | --- |
| insufficient write settling | 3 | -0.410 | 2 |
| 15pct reference/supply gain | 3 | 0.450 | 2 |

Both adverse first outputs are classified correctly, yet stage2 is already wrong. The first four instrumented stress seeds additionally show two first failures from raw-MAC/input perturbation outside the intended region and two from detector decisions while the nominal raw class was still correct. These diagnoses are not population fractions. The limiting requirement is the following weighted MAC's half-class error budget. Better references, lower S, larger capacitors or calibration may restore margin, but their energy/area cost remains unproved.

## Completion and finite boundaries

On the same128 frozen q=4 INT8 product workloads, using a compiled17-position layout:

| Schedule | Comparator decisions/product | Latency ns (model) | Rewrite events/product |
| --- | --- | --- | --- |
| Fixed32 | 1224.00 | 128 | 32.38 |
| Ideal completion lower bound | 380.64 | 28.78 | 32.38 |
| Explicit event queue | 333.25 | 9.32 | 32.38 |

The event controller charges initial scan, six persistent plus two guard comparisons per final node, dirty-node visits, request/ack, pending flags, commit revalidation, queue insert/pop, activity-tree toggles and final guard/root confirmation. No fixed32 workload exceeds the budget. The queue preserves the exact value and reaches quiescence. Its assumed .2 ns comparisons,1 ns transfers,.2 ns local settling and tree delays are timing assumptions; the latency advantage is not measured asynchronous hardware speed. Queue/tree energy can outweigh the comparator savings versus ideal free completion.

The exact finite-word policy uses [-2,2], declared guards, high-side index shift plus exponent update, optional exponent decrement, explicit low-side truncation and exact nearest-quantum ties-to-even rounding. Every lost value has an exact pair ledger. The physical boundary circuit routes the high digit and explicitly receives discarded low charge. Rounding, exponent comparison/addition and aligned residual/tile merges stay digital models; no automatic exponent-CIM or hardware rounding controller is claimed. Hardware exponent width/overflow trap remains a design choice.

## Energy gate

\[
E_\gamma=E_{CIM}+E_{detect}+E_{rewrite}+E_{regen}+E_{completion}+E_{routing}+E_{scale/exponent}.
\]

The matched baseline includes CIM+ADC+binary accumulation+requantization. The sweep includes ADC amortization1/64,1/16,1,4,16 and varies detector, packet, controller and all-in regeneration costs. The all-in regeneration parameter includes reference generation, reset/write, local clock and next-stage buffer. Differential charge's ideal signed current sources can recover/cancel supply current; its reference-port integral is particularly incomplete. Port measurements cannot rank transistor energy across the three mappings.

At the illustrative assumptions comparator0.1 fJ, rewrite(q3 reference)1 fJ scaled6/5, controller activity0.01 fJ, binary accumulation20 fJ and requantization5 fJ, the explicit-event extra cost before regeneration is 78.38 fJ/product. With1 fJ/position regeneration, ADC break-even is 70.38 fJ/conversion for one conversion/output. At100 fJ per ADC conversion, the allowed all-in regeneration budget is only **2.74 fJ/position**. Greater completion, buffering, routing or array overhead reduces it. Strong ADC amortization can eliminate the budget even before regeneration.

![Energy break-even](../spice/regen_results/energy_break_even.png)

![Regeneration and completion budget](../spice/regen_results/energy_regen_completion.png)

These assumed costs are parameter axes, not measured realistic cell costs. Common CIM cancels only when explicitly equal; the optimistic zero delta-array/routing/sidecar reference point is stated in energy_sweep.json. Reference-port energy is partial, not total supply energy or proof that an ADC has been beaten. Reused-storage reset/write checks cover both previous extreme digits; complete repetitive throughput and reference-bank fanout remain open.

## Verification, limitations and next decision

Half-step convergence checks cover nominal/direct, drift, stress, complete word and feedthrough cases. Every reported checkpoint class is compared; exact results and any differences are retained in supplemental_checks.json. Startup feedthrough IC consistency and delayed reference selection were corrected before final results. Tiny parasitic/control capacitors caused trapezoidal integration chatter in a diagnostic nonideal case; nonideal configurations use Gear integration with timestep checks. Method selection is predetermined by configuration, not by correctness. Numerical diagnostic attempts are excluded from the experiment manifest and are not counted as successful physical trials.

Full regression: 87 passed in 15.43s. The log and machine-readable exit status are in spice/regen_results/regression.log and regression_result.json. All eight half-step comparisons preserve every checkpoint class; all42 reused-storage resets/writes remain in the intended class.

Keep **q=4** as requested. Differential reference-selected voltage is the behavioral reference implementation; single-ended and charge alternatives remain controls, not eliminated hardware candidates. The restoration principle passes this behavioral gate. Stress failures establish that regeneration is conditional on the next MAC's margin, not a cure for bad references or wrong decisions. Physical energy superiority, calibrated class-error rates and asynchronous control are unestablished.

The next useful experiment is a loaded seven-state reference/write/driver cell feeding a tiny resistor MAC, with actual supply accounting and reference/comparator mismatch. Advance to a small resistor-crossbar prototype only if that cell fits the disclosed regeneration+completion budget. No FPGA or full crossbar is warranted yet.

## Reproduction and deliverables

From python/: `python -B run_regeneration_phase.py --trials 16 --workers 1`, then `python -B report_regeneration_phase.py`, then `python -B -m pytest -q -p no:cacheprovider`.

- Exact oracles, exponents, finite words, completion: gamma_regen_golden.py.
- Physical cell and connected generators: gamma_regen_spice.py; generated q4_regen_*.cir and .golden.json in spice/behavioral/.
- Main run manifest, source hashes and profiles: spice/regen_results/manifest.json.
- Seeded netlists, logs, exact golden checkpoints/CSV and decoded results: spice/regen_runs/v0000..v0085 and w0086..w0493, plus named supplemental runs. Other diagnostic/earlier-attempt directories are excluded by the manifest.
- Full statistics and bounds: phase_summary.json, summary.csv; exact control traces: completion_traces.json; finite boundary oracles: finite_word_golden.json.
- Supplemental convergence, reused-storage and adverse cases: supplemental_checks.json.
- Parameterized costs: energy_sweep.json; plotted results alongside this report under spice/regen_results/.

All paths above are relative to the project except the Python module names, which live in python/. Prior reports are historical and remain available.
