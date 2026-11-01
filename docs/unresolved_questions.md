# Unresolved questions and decision gates

## Matched130nm gate — 2026-10-24

The new calibrated comparison is in [engineering_report.md](../results/sky130_phase/engineering_report.md). Stop using100fJ ADC /2.74fJ-position as an absolute physical target. Current key unknown: can bounded local overflow/transfer/control fit the remaining per-dot budget at the actual event density? Under the retained schedule, overflow-only ceilings are0.398→0.016fJ/evaluation for N128→8192. A cheaper or shared schedule must be demonstrated, not assumed. Passive reference sharing is loading-limited; real driver/mismatch failures remain. No array prototype is authorized by readout-only break-even.

## Latest physical gate — 2026-10-24

[SKY130 transistor results](sky130_phase_results.md) identify a concrete unfavorable mapping: a passive reference ladder plus voltage-storage cell and biased OTA driver costs hundreds of fJ before complete region/control logic, and seeded mismatch frequently crosses next-MAC regions. The open question is now whether a different shared/passive/charge-oriented interface can remove these costs without shifting them into equally expensive precision or driver circuitry. Sharing must be electrically simulated, not assumed. A full transistor seven-state detector, completion and normalizer remain unresolved; do not build the full array yet.

Current lead is frozen q=4. Behavioral canonical restoration and connected analog MACs are evaluated in [regeneration_phase_results.md](regeneration_phase_results.md). Preserve q=3/5 comparison records; do not reopen broad base sweeps at this gate.

## Main physical bridge

The sampled local behavioral circuit passes the declared ideal cases. Can a transistor-level implementation reproduce this without expensive comparators, reference generation, sampling/control and charge transfer? Define coarse observables, q-unit transfers, threshold scaling, routing, event sequencing, settling and boundary behavior. Ideal software normalization proves algebra only.

Reject the proposed advantage on a tested mapping if total sensing/transfer/threshold-generation cost reaches or exceeds matched-accuracy ADC plus binary accumulation.

## Scale and class geometry

Which coefficients are guaranteed to lie on the declared digit lattice? How are fractional learned weights encoded/programmed? Can scaling preserve useful physical spacing at realistic fan-in? What scale generation, calibration and drift overhead is required? Five small persistent states do not cover arbitrary raw sums.

## Base, alphabet and precision

Optimize beta, q, D and finite layout together. Growing alphabets shrink region spacing. Fixed alphabets can have gaps at larger beta. Require empirical conversion precision, redundancy-aware codebook metrics and task accuracy rather than positions*log2(beta).

## Floating/block-floating datapath

Choose exponent width and overflow/underflow policy. Choose rounding, guard positions and finite mantissa window. Quantify loss of small block values. Different block exponents can be aligned positionally; different arbitrary physical scales need calibrated compensation. Renormalization needs routing and spare positions. Accumulation can require local digit normalization even with bounded current.

## Fabric and storage costs

Quantify learned conductance precision, signed/differential cells, product routing i+j=k, unity switching, fanout, width, driver cost, scale/exponent storage and programming. Sparse mathematical permutations are not a measured efficient crossbar.

## Neural compatibility

GELU, SiLU, softmax, RMSNorm/LayerNorm, reciprocal/division, square root and comparison/max are unresolved. Even sign/ReLU on a redundant full word needs a specified method. Audit an actual network for gamma-to-binary boundaries before claiming persistent transformer execution.

## System comparison

Include array and drive energy, extra cells, sensing, threshold calibration, normalization, routing, sidecar and boundary conversion. Compare at matched fan-in, accuracy, throughput and external interfaces. No present result establishes an energy, silicon-area or GPU advantage.

## Behavioral-phase decision

q=4 is frozen for this phase. Independent threshold offsets, seven-reference writes and actual following-MAC connections now have behavioral tests. Real comparator/reference supplies, driver bandwidth/load, repeated-use reset energy, calibrated process mismatch and asynchronous completion hardware remain gaps. A wrong class becomes a clean wrong digit after regeneration. Event-queue activity is a model, not completed controller circuitry. Review stress failures and confidence bounds before crossbar integration.
