# Local gamma-normalizer engineering phase

Historical local-normalizer phase. The next q=4 canonical-restoration and connected-MAC phase is in [regeneration_phase_results.md](regeneration_phase_results.md); its added reference, regeneration, driver and completion costs supersede energy budgets that omit those terms below.

Run date: 2026-10-23. Actual simulator: installed LTspice XVII; no FPGA or full crossbar.

## Outcome

Recommended next behavioral prototype: **q=4**, under the explicit score weights and energy assumptions below. This is not a silicon winner or a demonstrated energy advantage.

- Ideal circuit cases: 108; final class failures: 0; phase-region failures: 0.
- Nonideal actual-SPICE trials: 2376; trials with a wrong output class: 527. These deliberately mixed low/high stresses are not a manufacturing failure probability.
- Of those strict state failures, 516 are arithmetic errors/unresolved overflow and 11 emit a different but algebraically equivalent digit word.
- Sampled charge-model cross-checks against identical circuit realizations: 2052/2052 exact agreement of checkpoint classes. Leakage/current noise are excluded from this model comparison.
- Captured active-window boundary transfers: 36 ideal cases. Passive guard receivers preserve the exact value; these cases are not successes of a sealed three-node system.

The behavioral circuit uses capacitor states, finite comparator bands, a clocked switch/sample-hold and event-controlled current packets. Python provides only injection events and a fixed clock budget; it does not provide firing decisions. Readout uses the declared small digit/overflow regions. Continuous residuals are instrumentation. Exact output values are reconstructed only from emitted digit classes.

![Actual q4 LTspice charge-packet trace](../spice/results/q4_local_rewrite.png)

## Reproduce

From python/:

```powershell
python -B -m pytest -q -p no:cacheprovider
python -B run_normalizer_phase.py --spice-trials 12 --rare-trials 3000000
```

The runner caches completed jobs by inputs and source hashes. spice/golden contains CSV event streams and exact JSON phase/checkpoint traces. spice/runs contains each executed netlist, golden case, simulator log and region verification. Selected ideal raw waveforms are retained; discarded sweep waveforms can be regenerated from their retained seeds and decks.

## Protocol tradeoffs

| q | Scheme | Median rounds | Median triggers | Median unit rewrites | Overflow comparators/node |
| --- | --- | --- | --- | --- | --- |
| 3 | unit | 10 | 72 | 72 | 2 |
| 3 | cap2 | 7 | 54 | 72 | 4 |
| 3 | cap4 | 7 | 51 | 72 | 8 |
| 3 | pow2 | 7 | 53 | 72 | 6 |
| 3 | bulk | 7 | 51 | 72 | None |
| 4 | unit | 6 | 29.5 | 29.5 | 2 |
| 4 | cap2 | 4 | 23 | 29.5 | 4 |
| 4 | cap4 | 4 | 21 | 29.5 | 8 |
| 4 | pow2 | 4 | 22.5 | 29.5 | 6 |
| 4 | bulk | 4 | 21 | 29.5 | None |
| 5 | unit | 6 | 27 | 27 | 2 |
| 5 | cap2 | 4 | 19 | 27 | 4 |
| 5 | cap4 | 3 | 16 | 27 | 8 |
| 5 | pow2 | 3 | 17 | 27 | 6 |
| 5 | bulk | 3 | 16 | 27 | None |

Workload: 128 identical seeded INT8 scalar product pairs per candidate. Bulk is an unrestricted-quotient lower-bound control, excluded from circuit and energy claims. Capped schemes can reduce control events/rounds without reducing total moved charge proportionally. Comparator evaluations include idle synchronous observations; final digit readout is additional.

Transient-limit sweeps are in the JSON. They stream signed injections, defer normalization until the limit, and always flush to the persistent alphabet at the end. A larger temporary alphabet trades fewer normalization sessions for more range/classes; it is not a free precision gain.

## Physical model and boundaries

Reference parameters are illustrative and unoptimized: C=100 fF, unscaled unit step=0.1 V, packet top duration near 1 ns, phase period=4 ns. Ideal S is at least 4 and is compiled from the declared injection L1 budget with 10% margin to a 0.6 V absolute state limit; it is not selected from an observed answer. This is testbench range planning, not an input-normalization operation inside a weight-only MAC. Signed voltages stand for differential storage. Neither these component values nor the ideal behavioral sources establish a competitive energy implementation. The sample/hold, switch, reference generator and mirrored current sources are idealized.

Persistent alphabet is radius q-1 (5/7/9 states). A unit event removes q identical packets at the center and injects one at each neighbor, with all physical units divided by the same S. Positional significance is carried by node index, never exponential physical voltage. Three/five/nine refers to active normalizers; two extra passive receivers explicitly collect outward boundary transfers.

## Nonideal sweeps and scaling

Actual SPICE sweeps include threshold offset, independent comparator offset, capacitor mismatch, center/neighbor transfer mismatch, pulse-width mismatch, leakage, injection gain error, sampled voltage noise, piecewise-linear injected current noise, correlated threshold-bank error, and a combined condition. The MC input produces simultaneous neighboring firings and secondary events. Per-condition counts, confidence bounds, residuals, headroom and extra events are in the JSON.

Comparator offsets are input-referred per active cell and shared by that cell's signed threshold pair; the correlated component is shared across cells. Independent mismatch of every threshold inside a multiband detector is not covered by this unit-cell MC study. Final class discrimination uses the same declared offset/noise model. No canonical-amplitude output driver or subsequent learned-weight MAC is simulated; capacitor residuals remain real analog errors until an explicit regeneration driver is implemented.

| q | Error target | Largest tested scale with pointwise bound passing | Evidence |
| --- | --- | --- | --- |
| 3 | 0.001 | 32 | sampled charge model only |
| 3 | 1e-06 | 16 | sampled charge model only |
| 4 | 0.001 | 48 | sampled charge model only |
| 4 | 1e-06 | 16 | sampled charge model only |
| 5 | 0.001 | 32 | sampled charge model only |
| 5 | 1e-06 | 16 | sampled charge model only |

No actual-SPICE condition has enough independent trials to certify 1e-3 or 1e-6. The accelerated-model table has pointwise one-sided exact binomial confidence bounds, not a continuous Smax, not simultaneous bounds across all scales, and not a manufacturing guarantee. Its scope excludes continuous leakage/current noise, detailed analog sampling dynamics and transistor reference errors. Its 1e-6 confirmation, when enabled, uses a predeclared conservative S=16 rather than an adaptively selected maximum.

![Scaling confidence bounds](../spice/results/scaling_reliability.png)

## Energy break-even

N_ADC*E_ADC + E_binary_accum > delta_E_CIM + N_comparator*E_detect + N_unit_rewrite*E_unit(q) + N_trigger*E_control + E_routing + E_scale_sidecar

Common array energy cancels only when explicitly assumed equal. The grid includes incremental array cost and baseline ADC conversions per scalar product of 1/64, 1/16, 1, 4 and 16. Subunit counts make baseline conversion amortization explicit; the gamma counters still normalize each scalar product. No inference-level dot-product energy win can be inferred from this local workload.

Two scheduling cases are exposed: until-settled plus confirmation (requires unpriced completion logic), and fixed 32-phase observation. Unit-rewrite energy scales with the assumed packet count q+2; trigger control and detection are separate. Routing, scale/sidecar and extra array costs remain explicit parameters. No single ADC energy is assumed.

Detector activity is charged across a fixed compiled output/normalization layout: 21 positions at q=3 and 17 at q=4/5, from full INT8 encoding windows and an absolute-coefficient envelope. It does not shrink to each observed answer's active span. The global resource count is an extrapolation of local cells, not a SPICE simulation of a complete 17/21-node MAC.

![Parameterized energy break-even](../spice/results/energy_break_even.png)

## Candidate recommendation

| q | Combined score (lower is preferred) | Median compiled half-margin, V | Reference fixed-budget extra cost, fJ | Stress failure fraction |
| --- | --- | --- | --- | --- |
| 3 | 0.7348 | 0.015 | 253 | 0.186 |
| 4 | 0.6531 | 0.0111111 | 174 | 0.245 |
| 5 | 0.7031 | 0.0075 | 177 | 0.235 |

Weights: {"detector_energy_parameterized_fj": 0.15, "digit_products": 0.1, "inverse_margin": 0.1, "persistent_states": 0.1, "rewrite_energy_parameterized_fj": 0.15, "span": 0.1, "stress_failure_fraction": 0.1, "transient_peak_mean_units": 0.1}. All metrics are normalized to the largest candidate value. These weights are an engineering preference, not a physical law.

q3 unit rewrite 1 fJ, comparator decision 0.1 fJ, trigger control 0.5 fJ; fixed 32 phases; assumed, not measured. Energy-only sensitivity winner counts: {"3": 0, "4": 28, "5": 4}.

The parameterized score includes representation, product count, charge-transfer and comparator cost, persistent regions, transient range, weight-compiled class margin, and stress stability. Change the assumptions before treating the selection as a chip decision.

## Floating-point sidecar

The existing exact NativeFloat/block-dot model remains in place. Exponents are digital metadata; exponent addition/comparison, signed residual/tile merge, positional renormalization and finite-window truncation are covered by regression tests. Truncation drops out-of-window coefficients explicitly; no exponent-CIM hardware is introduced.

## Answer to the physical feasibility question

Local coarse normalization is behaviorally stable on the tested ideal schedules and can be verified from region outputs without a broad ADC. Actual transistor cost below ADC plus binary accumulation is still unestablished. The present physical requirements to resolve next are matched q:1:1 charge transfer, reliable scaled references, held comparator decisions, finite-boundary management, comparator/control activity and the cost of fixed clock observations. These are measurable engineering requirements, not yet an identified impossibility or a proven accelerator advantage.
