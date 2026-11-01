# Validation status: expanded native architecture

## Latest validation — 2026-10-24

See [SKY130 results](sky130_phase_results.md): 20/20 smoke tests; 132 structured sizing designs; 90 PVT cases, including four retained incomplete transients; 48 independent network mismatch trials; 30 corner/sizing comparator tests plus32 comparator mismatch trials; all49 transition pairs. Software regression is93 passed. An independent literal-deck replay matches batch checkpoint voltages within1 nV and per-position energy within0.04 fJ. Physical energy superiority and robust seven-state regeneration remain unestablished. The historical numerical baseline below is unchanged.

> Historical numerical baseline. Later evidence is in [normalizer_phase_results.md](normalizer_phase_results.md) and [regeneration_phase_results.md](regeneration_phase_results.md). The latter freezes q=4 and tests actual canonical capacitor outputs driving subsequent weighted analog stages; circuit-unimplemented statements below are historical.

Date: 2026-10-22. Seed: 0. Exact software references and synthetic diagnostics only.

## What was run

- Expanded unit suite: 45 tests (see final test output).
- Default architecture-v2 runner: results/native_baseline.json, schema 2.
- Historical fixed-unscaled-bin control rerun: results/fixed_bins_control.json.
- Original results/normalized_baseline.json retained unchanged as historical evidence.
- No SPICE/circuit or application inference benchmark was run.

Commands from python/:

~~~bash
python -B -m pytest -q -p no:cacheprovider
python -B run_experiments.py --seed 0 --noise-trials 2000 --output results/native_baseline.json
python -B run_fixed_bins_experiments.py --output results/fixed_bins_control.json
~~~

## Scale-aware exactness

| Test | Result |
| --- | --- |
| Legal digits at 106 positive rational scales | 530/530 exact |
| Direct scaled five-state detector on raw dot products | 3/22 exact |
| Expanded transient detector, weight-compiled range | 22/22 exact |
| Ideal software normalization, then five-state sensing | 1542/1542 exact |

The supplied example now works: I=0.2, S=5 emits digit 1. Scaled bins return the unscaled coefficient, so the decoder does not multiply by S again.

The old 3/22 result came from deliberately fixed unscaled thresholds. It is no longer the detector acceptance test. The updated direct-five-state run also gets 3/22, but the remaining limitation is raw coefficients outside the persistent alphabet, rather than a failure to move thresholds. Physical current headroom does not ensure a raw sum is one of five legal digit values.

The direct transient reference uses between 5 and 14,337 classes, conservatively compiled from weights/layout rather than observed answers. This is not a demonstration of cheap coarse readout. The ideal normalization benchmark uses an exact software rewrite before sensing, not an implemented physical network. Its aggregate unit-event count over 257 raw coefficients is 14,060, counted once per coefficient, not once per scale.

Arbitrary positive-scale legal-class tests establish sensing algebra, not physical-range validity for every tested scale. Fractional learned-weight coefficients require a separate lattice/encoding contract.

## Synthetic sensing errors

Each row uses 2,000 balanced legal-digit trials. Additive sigma=0.002 normalized physical current; combined adds comparator offset sigma=0.002 and relative gain/cell mismatch sigma=0.01. Random draws are paired across scales. These are illustrative units, not device specifications.

| S | Spacing / additive sigma | Additive class error | Combined class error |
| --- | --- | --- | --- |
| 1 | 500 | 0.00% | 0.00% |
| 5 | 100 | 0.00% | 0.00% |
| 20 | 25 | 0.00% | 0.00% |
| 100 | 5 | 0.80% | 6.05% |
| 500 | 1 | 48.10% | 56.80% |

Separate gain, conductance, comparator-offset and correlated-bank-threshold results are in the JSON. Relative mismatch and absolute physical noise have different scale dependence. Zero observed errors is not proof of zero probability: for 2,000 independent trials, the usual one-sided 95% zero-event bound is about 0.15%.

## Base family and alphabet tradeoff

All 256 INT8 encodings and 128 shared product cases per q are exact. This reference uses radius q-1, so the growing alphabet is explicitly charged.

| q | beta | States | Median INT8 span | Median digit products | Median unit rewrites/product |
| --- | --- | --- | --- | --- | --- |
| 3 | 2.618 | 5 | 9 | 43.5 | 72 |
| 4 | 3.732 | 7 | 7 | 25 | 32 |
| 5 | 4.791 | 9 | 5 | 20 | 24 |
| 6 | 5.828 | 11 | 5 | 15 | 20 |

These event counts are exact arithmetic-resource metrics, not physical latency/energy. Fixed-five-state codebooks are evaluated separately. At five positions, measured maximum absolute error on the unit grid worsens from q=5 to q=6 because fixed-alphabet representational gaps outweigh the larger radix.

## Floating point and unity routing

1,200 multiply/add/positional-renormalization checks have 0 failures. All block dot and differently exponented block-merge oracle checks pass. A two-stage unit test consumes actual prior output digits without re-encoding the expected answer.

The illustrative signed-eight-bit exponent range overflows in 48/400 product cases; Python itself uses unbounded integers. Exponent width/overflow policy remains a format decision.

Sparse routes at exponent gaps 0/2/8/32 use four nonzero signed unity connections on an unrestricted layout. With window [-7,0], gaps 8 and 32 discard two low-order terms. Adding two leading digits of 2 produces raw coefficient 4 and needs a local rewrite. Routing and normalization therefore remain real costs.

## Block exponent granularity

Synthetic 8x8 tensor, q=3, eight mantissa positions. All values start nonzero. Relative L2 emphasizes large entries; the zeroed-value count exposes loss of small entries.

| Group | Exponents/operand | Relative L2 conversion error | Nonzero values rounded to zero |
| --- | --- | --- | --- |
| tensor | 1 | 0.00109681 | 36/64 |
| channel | 8 | 0.000771165 | 27/64 |
| tile | 4 | 0.000774306 | 35/64 |
| small_block | 16 | 0.000603005 | 12/64 |
| scalar_reference | 64 | 0.000249087 | 0/64 |

The JSON includes 4/8/12 positions, activation/weight exponents, dot error from boundary conversion, native digit products, normalization events and unity merge connections. No analog noise is injected in this block experiment; ideal normalization and unrestricted positions are assumed. This is not a trained-model accuracy result.

## Empirical precision and format conversion

Codebooks enumerate every five-state word at 3 and 5 positions and deduplicate exact algebraic values. Minimum step, largest interval gap, and error over 1,025 unit-interval grid points are separate metrics. Redundancy is not counted as independent bits.

Boundary conversions use a greedy offline encoder. INT8 and finite FP8 E4M3FN bit patterns are exhaustive; BF16 and FP16 each use 256 signed samples including extremes/subnormals. Unsupported NaN/infinity semantics are not tested. Here are q=3 results at 12 positions:

| Target | Samples | Exact conversions | Maximum relative conversion error | Empirical relative bits |
| --- | --- | --- | --- | --- |
| INT8 | 256 | 256 | 0 | exact on tested domain |
| FP8_E4M3FN | 254 | 84 | 2.61719e-05 | 15.22 |
| BF16 | 256 | 4 | 3.23043e-05 | 14.92 |
| FP16 | 256 | 4 | 3.00786e-05 | 15.02 |

Relative bits are -log2(maximum relative conversion error) on the specified sample, not guaranteed format equivalence or physical ENOB. INT8 exactness is checked with pairs, avoiding false hundreds-of-bits claims from Decimal roundoff. Fixed D with q=6 retains representational gaps even at more positions; no monotonic radix-density claim is justified.

## What is established and what remains

Established in software: scale-aware class recovery on the declared lattice; exact generalized arithmetic; native exponent/shift/block-dot identities; finite-window loss accounting; synthetic sensing-error trends and finite conversion metrics.

Open: a physically realizable low-state local normalization/readout network, scale/reference-generation cost, finite-width policies, learned fractional-weight mapping, nonlinear layers, measured device noise/energy/area, and matched-accuracy application performance. Exponent metadata separates dynamic range from analog headroom; it does not eliminate mantissa precision limits.
