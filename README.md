# Gamma-CIM

## Fully dynamic SAR normalization result (2026-10-24)

[The dynamic q4 bulk-normalizer experiment](results/sar_normalizer_phase/engineering_report.md) is complete. A real dynamic comparator, switched references and one bulk 4m:m:m packet phase remove intentional continuous detector bias. The tested +/-15 SAR cell costs 4.369 pJ per invocation, passes 31/31 nominal quotients but only 22/31 transfer-charge checks and 30/31 physical-class checks. Logical-zero idle improves to 94.70 nW from the rejected topology's 691.7 nW; hot mixed-corner leakage and mismatch still fail the robustness gate.

Exact workload recomputation tests +/-7, +/-15 and +/-31 windows, three quotient schemes and three ingress schedules without floating-beta correctness. Even the favorable +/-31 schedule needs 3,186 quotient invocations at N=128, including 981 nonzero bulk transfers: only 1.99x fewer transfers than its own unit control. Its measured-price readout-plus-normalization projection is 16.765 nJ/output against the unchanged 14.857 pJ ADC baseline. This is a negative local-circuit result, not a measured functional complete dot or a universal bound on Gamma arithmetic. Missing routing/completion/array costs cannot rescue it. All required outputs, 175 completed cases, two independent netlist audits, four unresolved exploratory attempts and 118 passing software tests are documented in the report. Prior results below remain historical.

## Self-triggered normalization result (2026-10-24)

[The asynchronous q4 transistor experiment](results/self_triggered_q4/engineering_report.md) now measures threshold, charge-transfer, guard, reset/re-arm and reference supply work without periodic polling. The tested 101-MOS implementation costs 4.219 pJ/input-to-rearm firing and 691.7 nW quiet idle. Nominal signed branch packets are approximately 4:1:1, but voltage-state checks, near-threshold coverage and PVT/mismatch robustness fail. Under the retained workload, this topology exceeds the most favorable 2.425 fJ/firing allowance by 1740x. It is rejected; this is not a universal lower bound on self-triggering or a measured all-in Gamma system. See the report for 43 completed characterization cases, literal decks, all failures and conditional workload projections. Prior matched-readout results below are preserved.

## Current matched130nm decision (2026-10-24)

The technology-matched ADC calibration supersedes the old100fJ ADC /2.74fJ-position interpretation. Read [the matched130nm report](results/sky130_phase/engineering_report.md). The primary6/7/8/9ENOB cases cost about14.9–32.6pJ per output; the5–7-position Gamma readout estimate costs4.8–6.7pJ before normalization and missing control. This is conditional headroom, not a measured system win. The actual boundary still has PVT/mismatch failures.

The retained per-product normalization schedule is now the main energy gate: overflow decisions must fit about0.398→0.016fJ each as accumulation grows128→8192, with rewrite/control assumed free for those ceilings. Passive reference-sharing tests also show that reducing ladder resistance to restore settling removes the simple per-load power saving. Next characterize a bounded low-energy q=4 overflow/transfer primitive and fix driver/reference matching; do not build a crossbar yet.

All prior LTspice and SKY130 reports/data remain historical and unchanged. New outputs are in results/sky130_phase/ and baseline_adc/. q=4 stays frozen; exponents remain digital. No ASIC flow or full CiMLoop stack was installed.

Research prototype for a native datatype consisting of sign, a digital/block exponent, and a redundant algebraic mantissa:

\[
X=\epsilon\,\beta_q^e\sum_k d_k\beta_q^k,\qquad
\beta_q=(q+\sqrt{q^2-4})/2,\quad q\ge3.
\]

The current lead is **q=4, beta=2+sqrt(3), D={-3,-2,-1,0,1,2,3}**. q=3/5 remain historical comparison candidates. Analog CIM computes mantissa MACs; local normalization bounds coefficients; finite-region sensing and canonical restoration prepare the next analog stage. Exponents remain digital metadata.


## SKY130 transistor phase — 2026-10-24

The narrow q=4 regeneration experiment is now executed with ngspice 47 and pinned SKY130A primitive models. See [SKY130 results](docs/sky130_phase_results.md). All 20 device smoke tests pass. The 132-design write/driver isolation sweep finds a minimum nominal working known-path cost of 644 fJ/position; adding one real threshold yields about 749 fJ/position with MIM storage. Five decisions and decoder/completion logic remain ideal, so this is a partial measured cost, not a complete boundary energy claim.

The energy candidate returns all 12 next-MAC outputs correctly at 29/45 PVT points; 15/16 seeded PDK-mismatch circuits have at least one wrong output. The larger-device candidate also has failures, including numerical nonconvergence retained separately. The tested voltage-storage/reference/OTA interface does not establish an advantage over the illustrative 100 fJ ADC baseline. Redesign it before a full crossbar; retain LTspice as the behavioral reference. Exponents stay digital and q=4 stays frozen.

Reproduce from python/: `python -B run_sky130_phase.py --workers 2 --mc-trials 16`, then `python -B report_sky130_phase.py`. Required CSVs, eight plots, seeds, literal netlists, golden states and simulator logs are under spice/sky130/. The minimal PDK and portable ngspice provenance is in spice/sky130/pdk_manifest.json.

## Previous behavioral phase (2026-10-23)

The current phase tests **canonical regeneration and connected analog stages**. Read [regeneration results](docs/regeneration_phase_results.md) and the [canonical interface specification](spice/canonical_regeneration_spec.md). From python/:

```powershell
python -B run_regeneration_phase.py --trials 16 --workers 1
```

This creates seven-reference capacitor write cells, a two-stage non-unity-weight word path, depth 2/4/8/16 circuits, nonideality trials, exact finite-word/exponent policies, an explicitly modeled completion controller and energy budgets. No Python oracle drives an intermediate analog stage. Behavioral reference-port energy is partial; transistor-level energy superiority remains open.

### Previous local-normalizer phase

The local primitive now has an exact golden generator, signed 3/5/9-active-node LTspice circuits, four bounded quotient protocols, physical-nonideality sweeps, scaling confidence bounds, and parameterized energy curves. All 108 ideal circuit cases pass; 2,376 stress trials retain failures explicitly. The combined, assumption-dependent recommendation is q=4 for the next prototype.

Read [engineering results](docs/normalizer_phase_results.md) and [circuit contract](spice/local_normalizer_spec.md). From python/:

```powershell
python -B run_normalizer_phase.py --spice-trials 12 --rare-trials 3000000
```

The large rare-event sample is a separately labeled charge-domain model, not millions of SPICE runs. No actual-SPICE 1e-6 reliability or transistor energy advantage is claimed.

## Scale-aware readout correction

Normalize **weights only** for physical current headroom. The existing input digits are unchanged. If physical current is c/S, digit centers are d/S and thresholds are (d+1/2)/S. Thus current 0.2 at S=5 decodes to digit 1. Scale-aware output already represents c: do not multiply it by S again.

This is exact for c in the declared integer alphabet. Raw MAC coefficients can exceed that alphabet even when physical current is bounded. Five scaled bins cannot directly return arbitrary coefficients. Larger transient sensing or an actual local normalization mechanism is required. See [architecture](docs/architecture_summary.md) and [results](docs/validation_status.md).

The old fixed-unscaled-bin 3/22 result is a historical control, not a valid test of the updated scale-aware detector. The updated direct-five-bin run also scores 3/22 on this frozen suite, for a distinct reason: raw coefficients outside the alphabet. Neither result disproves the proposed normalization network.

## Run

From python/:

~~~bash
python -B -m pytest -q -p no:cacheprovider
python -B run_experiments.py --output results/native_baseline.json
~~~

The default runner now covers scaled bins; additive noise, comparator offset, gain, per-cell conductance mismatch and correlated threshold errors; q=3..6; exact floating arithmetic; block/tensor/channel/tile encoding and accumulation; sparse unity routing; finite-layout losses; and empirical INT8/FP8/BF16/FP16 conversion precision.

Historical controls remain available:

~~~bash
python -B run_fixed_bins_experiments.py --output results/fixed_bins_control.json
python -B run_legacy_experiments.py
~~~

## Project map

- docs/architecture_summary.md: current proposal, equations, and limitations.
- docs/architecture_update_20261022.md: traceability to the supplied update.
- docs/validation_status.md: reproduced results and interpretation.
- docs/novelties.md; docs/unresolved_questions.md: hypotheses and unresolved gates.
- python/gamma_native_model.py: generalized field arithmetic, scaled sensing, native floats, block dots and routing.
- python/run_native_experiments.py: expanded deterministic studies.
- python/gamma_bounded_model.py: weight-only physical range planning; its detect_five is the historical fixed-bin control.
- python/gamma_exact_model.py: original q=3 integer oracle and boundary encoder.
- python/experiment_plan.md: domains, metrics and acceptance criteria.
- spice/: updated nominal-path plan and retained optional normalizer skeleton.
- agent_handoff.md; roadmap.md: implementation priorities.

The clocked behavioral normalization loop is now exercised in LTspice. Transistor-level normalization, measured energy advantage, and GPU advantage remain unestablished.
