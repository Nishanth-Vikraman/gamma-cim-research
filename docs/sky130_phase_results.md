# q=4 SKY130 regeneration: engineering result

Run date: 2026-10-24. Simulator: ngspice47, regular SKY130A 1.8 V NMOS/PMOS; no full array or digital ASIC flow.

## Decision

**The tested voltage-storage / OTA-driver / resistor-MAC implementation does not establish an energy advantage. Redesign this boundary before a crossbar prototype.** The transistor writer and loaded driver can work, but their measured cost alone is far above the aggressive reference budget. The complete seven-state detector is not finished: one important threshold is transistorized and five remain behavioral. Missing decoder, reference selection and completion costs can only tighten the disclosed budget under the stated no-recovery assumptions.

The structured isolation sweep contains 132 designs, with 32 meeting the declared nominal MAC/error/operating checks. Its minimum known-path point costs 643.7 fJ/position before a real region decision is included. It uses 25 fF storage, 0.42 um NMOS selector W, 400 ns periods, 4 Mohm total reference ladder, 10 Mohm bias resistor and 1 Mohm MAC resistor scale. Smaller storage is an observed tradeoff, not an assumption that100 fF was optimal.

At tt/1.8 V/25 C, the selected hybrid circuit measures **720.4 fJ of signed1.8 V supply work per position**. Gross delivered supply energy is 727.9 fJ; a conservative external control-pin charging allowance adds 20.7 fJ, giving 748.6 fJ for the implemented blocks. These are per-operation means over the declared 12-transition workload. They are not claimed as fully measured seven-state detector/control-system energy.

| Budget fJ/position | Measured VDD / budget | Known path + pin allowance / budget |
| --- | --- | --- |
| 1 | 720.4x | 748.6x |
| 2.74 | 262.9x | 273.2x |
| 5 | 144.1x | 149.7x |
| 10 | 72.0x | 74.9x |
| 20 | 36.0x | 37.4x |

This triggers the requested cost gate for these tested circuits. It is not proof that every possible Gamma-CIM physical representation is impossible. The quantitative limitation is this loaded voltage interface, together with the measured reference and driver cost; the algebraic datatype remains exact.

## Electrical path and architecture integrity

The native lead remains q=4, beta=2+sqrt(3), digits -3..3. Position significance is carried by the index, never by exponentially increasing voltage. No q sweep, high-voltage devices, general ADC/DAC, RRAM, FPGA or full CIM array was introduced.

Nominal canonical levels are VCM=0.6 V plus d*25 mV. The resistor ladder has seven integer-state taps plus six half-step taps, from the 1.8 V supply. Stage 1 stimulus is an imperfect analog state at the previous word path's raw scale S=12 (8.333 mV per coefficient); the restored interface uses S=4 (25 mV). All positions share the same scale policy. Inputs are not vector-normalized. Supply sweeps track the declared common reference scale; an additional untracked reference drift would need calibration/error budget.

The physical path is ladder -> seven MOS switches -> storage -> SKY130 closed-loop output OTA -> two resistors and20 fF next-MAC load. The passive sum is (2*Vdrive+Vbias)/3, with external operand Vbias=VCM-d*Delta. Declared scale-aware output bins have spacing Delta/3, so the recovered coefficient is2*d_x-d. The external bias is fixed from the original input workload, not a Python re-encoding of an intermediate result. The actual stored/transistor-driven node supplies the next MAC.

The direct control uses an actual OTA with gain3 feedback to compensate the known raw/interface scale while preserving the analog residual. It does not receive an expected canonical answer. The isolation sweep also compares an additional shared write-buffer topology, two-stage and five-transistor-core OTAs, continuous and gated operation, transmission gates and NMOS-only switches. The shared write-buffer closes its storage aperture before disconnecting its reference input, avoiding amplified mux charge injection.

Ideal one-hot decisions drive actual transistor inverter interfaces. A hybrid case replaces the +0.5 raw-digit threshold with a complementary StrongARM comparator; the other five decisions and sequencing remain ideal. The -3..3 debug label is not a physical +/-3 V logic rail. Every MOS connection uses the declared supply range and is monitored for external gate/drain/source and source/body limits. Invalid operating points are retained and excluded from nominal feasibility.

The input contract is a previously normalized legal digit. A complete input overflow/invalid-state controller is not transistorized. Final MAC readout explicitly flags out-of-alphabet regions and does not repair them in software. External analog-input and MAC-bias port work is listed separately; including it cannot be assumed free in a full-system comparison. The reported per-boundary cost excludes it only as an explicitly separate MAC operand cost, not as an energy-recovery credit.

## PDK, simulator and Stage0

Pinned repository: `84a539fb024ab87ad9b46c9008c0deb298e250c0`; open_pdks provenance hash: `b344c97eacc2aaf8e14ae7e43e2e9dc0871de2c0`. The extraction's upstream provenance is recorded, not independently rebuilt. Apache-2.0 license and NOTICE are retained. Model data is 52.2 MB excluding Git metadata. Only libs.tech/ngspice and libs.ref/sky130_fd_pr/spice were installed, plus the portable simulator runtime. No OpenLane, Magic, KLayout, standard-cell, I/O or SRAM library was installed.

All 20 NMOS/PMOS/inverter/transmission-gate smoke checks pass across tt,ss,ff,sf,fs. The typical100 fF transmission-gate test exposes about2.8 mV of combined hold shift: correct switching is not an ideal canonical write. DC Ron/leakage and reference-bank curves are retained in model_validation.json. The detector and output-driver energy is actual supply work; behavioral-source power is never substituted for measured transistor comparator energy.

The independent diagnostics establish: identical devices match with mismatch disabled; different instances and seeds differ when MC_MM_SWITCH is enabled; the same seed replays identically; batch temperature control agrees with a literal .temp standalone run. MC_PR_SWITCH remains0 so discrete process corners and local statistical mismatch are distinct. The PDK is unmodified and Git-clean. Seeds are set before circuit reset/model evaluation.

Core-model documentation gives1.95 V VDS/VGS limits. The proposed1.98 V point was therefore excluded; the PVT sweep uses1.62,1.8,1.9 V and-40,25,100 C. This is model-level validation, not foundry signoff. The five FET corners use the library's typical R/C corner; capacitor extremes and ideal-resistor mismatch are separately bracketed.

## PVT, mismatch and chaining

Both the minimum-energy and a larger-margin Pareto candidate are carried through PVT. Preferred MIM storage is used in PVT; ideal/MIM pairs isolate storage-model effects at nominal conditions. The MIM geometry comes from the supplied area/perimeter model, with large reference decoupling implemented as parallel units of at most1 pF each. This includes the model's series parasitics, not physical-layout extraction.

| Candidate | PVT points | All12 MAC results correct | Incomplete transients | Full runs inside monitored device limits | Passes both class and monitored limits | Worst observed residual (digits) |
| --- | --- | --- | --- | --- | --- | --- |
| energy | 45 | 29 | 0 | 36 | 23 | 0.524 |
| margin | 45 | 35 | 4 | 31 | 25 | 1.1 |

Incomplete transients are retained as numerical nonconvergence, never as passes. The failing logs name ideal B-selector branches, so these cases do not by themselves prove physical transistor oscillation. Their unobserved checkpoints remain unresolved; only complete cycles contribute energy statistics. The PVT spread plot includes completed simulations only. The larger-device candidate was chosen for lower residual in the ideal-detector sizing sweep; it is not a proven robust design after transistor comparator loading is attached.

| Network bracket | Networks failing acceptance | Observed checkpoints wrong | Unresolved checkpoints | 95% one-sided upper error bound |
| --- | --- | --- | --- | --- |
| mc_pdk | 15/16 | 119/192 | 0 | 99.7% |
| mc_architecture | 4/16 | 15/192 | 0 | 48.4% |
| mc_margin | 12/16 | 100/192 | 0 | 91.0% |

Trials are independent seeded whole circuits. Twelve checkpoints inside each circuit are correlated and are not twelve independent yield samples. PDK mismatch uses the provided local FET/MIM machinery. Architecture-only trials disable it and add ideal-resistor mismatch, explicit capacitance error and an input-referred systematic comparator offset. The two mechanisms are not silently double-counted. These are uncalibrated brackets; deterministic transient analysis does not automatically include device thermal/flicker noise. No1e-3/1e-6 reliability certification follows from these sample counts.

| Storage | Initial residual (digits) | Regenerated correct | Direct correct |
| --- | --- | --- | --- |
| ideal | 0 | 12/12 | 3/12 |
| ideal | 0.1 | 12/12 | 3/12 |
| ideal | 0.35 | 12/12 | 2/12 |
| mim | 0 | 12/12 | 3/12 |
| mim | 0.1 | 12/12 | 3/12 |
| mim | 0.35 | 12/12 | 2/12 |

The S12-to-S4 direct branch already fails at zero injected residual: its extra gain3 feedback divider loads the finite common-mode ladder and its amplifier has finite gain/offset. This table therefore cannot independently prove that regeneration alone prevents accumulated error. The following control holds input/output scale atS4, uses identical unity-gain physical drivers and MIM storage, and isolates the transistor write path with ideal latched seven-state decisions. No scale-conversion amplifier or expected-answer source is inserted.

| Common-scale initial residual (digits) | Regenerated correct | Direct correct |
| --- | --- | --- |
| 0 | 12/12 | 12/12 |
| 0.1 | 12/12 | 12/12 |
| 0.35 | 12/12 | 1/12 |

## Comparator and primitive characterization

The isolated comparator sweep spans W=1/4/16 um and L=0.15/0.5 um at all five corners. Its declared legal-center overdrive is +/-4.167 mV at nominal supply. Offset is a sampled transition bracket, not a precision offset measurement or calibrated yield claim; unresolved/nonmonotonic transfers do not receive a fictitious offset. Delay is measured from the clock's 50% crossing. At exactly zero differential input either resolved polarity is permitted; unresolved nonzero cases are retained.

| W (um) | L (um) | Corner | Wrong legal-center decisions /2 | Mean energy incl. clock allowance (fJ) | Max decision delay (ns) | Offset bracket (mV) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.15 | tt | 0 | 109.2 | 0.8 | [-0.25, 0] |
| 16 | 0.5 | tt | 0 | 1422.4 | 4.74 | [0, 0.25] |
| 1 | 0.15 | ss | 0 | 112.5 | 1.45 | [0, 0.25] |
| 16 | 0.5 | ss | 0 | 1491.5 | 6.85 | [0, 0.25] |
| 1 | 0.15 | ff | 0 | 110.4 | 0.522 | [0, 0.25] |
| 16 | 0.5 | ff | 0 | 1367.2 | 3.4 | [0, 0.25] |
| 1 | 0.15 | sf | 0 | 109.0 | 1.29 | [-0.25, 0] |
| 16 | 0.5 | sf | 0 | 1380.7 | 6.31 | [-0.25, 0] |
| 1 | 0.15 | fs | 0 | 112.8 | 0.597 | [0, 0.25] |
| 16 | 0.5 | fs | 0 | 1539.0 | 4.11 | [-0.25, 0] |

The table includes actual comparator VDD work plus the external clock-pin charging allowance. Separate fields retain the two components. A larger comparator does not guarantee a better loaded boundary: input capacitance, reference kickback and decision energy grow too. The wide L=0.5 um example uses about 53.3 pF of physical threshold decoupling versus 1 pF for the smallest input pair. Those capacitors are present in the circuit, and their prebiased cold-start energy is excluded explicitly.

| W um | L um | Independent seeds | Trials with wrong legal-center decision | Wrong decisions /2n |
| --- | --- | --- | --- | --- |
| 1 | 0.15 | 16 | 13 | 13 |
| 16 | 0.5 | 16 | 5 | 5 |

The switch diagnostics retain Ron versus level and width, leakage, and loaded reference-bank DC curves. Canonical residual, write hold drift, within-cycle peak reference droop, settling and each component's energy are in the per-case metrics; storage ideal/MIM comparisons use identical external timing and loads. Clock feedthrough and charge injection arise from the MOS models and declared estimated junction geometry. Supply variation and the declared resistor/capacitor/comparator perturbations are included; layout-dependent parasitic uncertainty and transient sample/thermal noise remain unmeasured.

For the selected nominal MIM path, storage error reaches0.208 digit, driver error0.164 digit, and hold drift0.0106 digit. It passes the next-MAC class test but does not settle inside0.1 digit on every cycle. In the CSV, max_settling_ns=1e9 is an optimization sentinel for failure to reach that target during observation, not a measured one-second settling time; cycles_not_settled_to_point1digit and max_observed_settling_ns make censoring explicit. Ideal/MIM nominal path errors are nearly identical in these unextracted typical-model comparisons, so nominal MIM parasitics are not the dominant observed limitation.

All 49 previous/next digit pairs are also exercised through actual preceding writes; 49/49 return the correct next-MAC class. Their mean known-path cost is 749.0 fJ; worst case is 751.9 fJ. The matrix includes no-op pairs that still execute reset/write; sparse skip control is not implemented or assumed free. Per-trial records retain seed, corner, every MOS W/L, reference levels, stored/driven voltages, final class and component energy.

## Why the budget is difficult

The ladder tradeoff is real: static energy scales as VDD squared times period divided by R, while write settling grows with R times storage-plus-gate capacitance. Increasing R without changing the electrical write aperture is not a free saving. The sweep includes coupled R/C/timing choices rather than shrinking storage alone. Reference-bank sharing has not been electrically demonstrated. Even subtracting the entire nominal reference-bank cost leaves 424.3 fJ/position for the implemented remaining blocks: an optimistic free-reference bound, not a sharing result.

The loaded output OTA also spends bias/drive current for gain, settling and resistor-load current. Clock gating can lose its nominal saving through wake-up recovery, so continuous low-bias operation was tested too. Returning a voltage inside its own +/-0.5-digit class is insufficient: the following weight 2 needs residual comfortably below0.25 digit before MAC and sensing errors. PVT and mismatch results retain this stricter requirement.

| Bracket | Wrong next-MAC outputs | Wrong outputs despite correct selected class | Wrong with correct selection and storage error <0.1 digit | RMS driver-minus-storage error (digits) |
| --- | --- | --- | --- | --- |
| pvt_energy | 89 | 89 | 5 | 0.0708 |
| pvt_margin | 12 | 7 | 0 | 0.0244 |
| mc_pdk | 119 | 76 | 5 | 0.479 |
| mc_architecture | 15 | 13 | 0 | 0.0276 |
| mc_margin | 100 | 87 | 87 | 0.345 |

This attribution uses actual selected, stored and driven checkpoints, not an offline correction. Errors can overlap or compensate, so the columns are diagnostics rather than disjoint causal probabilities. In the larger-device PDK bracket,87 wrong outputs occur even with correct class selection and storage error below0.1 digit. The measured driver-minus-storage residual therefore identifies loaded analog buffering/matching as a concrete unresolved requirement. In the energy candidate, small-capacitor write/feedthrough error also consumes much of the next-MAC margin. Simply improving the class detector does not fix either effect.

Input-device upsizing improves a comparator's matching bracket but increases reference kickback. The physical threshold interpolation and its MIM decoupling are priced in the hybrid path; decoupling scales with input gate area. This can require much more capacitance than the storage itself. Measurements use a prebiased steady-state reference bank. Cold-start charging is not included, which favors Gamma; it must be amortized explicitly before claiming a complete operating energy.

## System break-even

The new calculation replaces the old hypothetical persistent readout rather than counting it twice. The remaining optimistic normalization/control allowance is 64.78 fJ/product on the frozen17-position workload. The one-threshold hybrid known cost is measured;3- or6-decision scenarios add actual comparator-primitive energy, with unbuilt threshold selection/latch/completion cost swept separately. A 3-comparison scheme assumes an unbuilt adaptive shared-reference controller, not demonstrated hardware.

With zero extra array/routing/control cost, one ADC conversion per product and the old25 fJ binary-accumulation+requantization allowance, the optimistic3-decision ADC break-even is **16.38 pJ/conversion**. The 6-decision estimate is 21.80 pJ. Strong ADC amortization raises these thresholds. A mathematical crossing at a sufficiently expensive ADC is not a measured, matched-accuracy energy advantage. The existing100 fJ ADC reference point is not competitive for this cell.

## Required plots

![Energy versus storage](../spice/sky130/results/energy_vs_storage.png)
![Energy versus residual](../spice/sky130/results/energy_vs_error.png)
![Settling versus buffer](../spice/sky130/results/settling_vs_buffer.png)
![Following MAC versus residual](../spice/sky130/results/next_mac_vs_residual.png)
![Energy breakdown](../spice/sky130/results/energy_breakdown.png)
![PVT spread](../spice/sky130/results/pvt_spread.png)
![Monte Carlo residual](../spice/sky130/results/monte_carlo_residual.png)
![ADC break-even](../spice/sky130/results/adc_break_even.png)

## Reproduce and inspect

From python/: `python -B gamma_sky130.py`; `python -B run_sky130_phase.py --workers 2 --mc-trials 16`; `python -B validate_sky130.py`; `python -B report_sky130_phase.py`.

Every case retains a standalone circuit.cir, seeded replay.cir, golden metadata, simulator log, parameters and metrics. Full waveforms are kept for representative cases with lossless gzip; all checkpoint voltages, energy integrals and classes are retained for every reported case. Unreported incomplete old traces were archived losslessly to control disk use. Binary output and bounded batches reduce overhead; neither changes the PDK. Local decoded gamma states, not offline numerical reconstruction, determine the result.

Software regression: `93 passed in 19.15s`. Independent literal-deck replay of the selected loaded path is recorded in results/standalone_path_audit.json, alongside the model random-seed and temperature diagnostics.

Required results are under results/: sky130_phase_summary.json, energy_breakdown.csv, pvt_summary.csv, monte_carlo_summary.csv, sizing_sweep.csv and regen_transition_matrix.csv. The previous LTspice reference remains unchanged. Do not build an array or FPGA from these results. The next decision is a redesigned shared/passive/charge-oriented reference-and-driver primitive, with real fanout, load and control accounting; a datatype novelty claim does not remove these electrical costs.

Sources: [extracted models](https://github.com/chrisklaus83188/sky130-ngspice-models), [device voltage limits](https://raw.githubusercontent.com/google/skywater-pdk/main/docs/rules/device-details.rst), [portable simulator](https://downloads.sourceforge.net/project/ngspice/ng-spice-rework/47/ngspice-47_64.7z). Local licenses, source snapshots and hashes are recorded in pdk_manifest.json.
