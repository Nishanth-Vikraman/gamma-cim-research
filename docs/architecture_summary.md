# Gamma-CIM architecture: native algebraic floating point

## Calibrated physical comparison — 2026-10-24

[The matched130nm phase](../results/sky130_phase/engineering_report.md) replaces the arbitrary100fJ ADC threshold with the pinned published estimator and survey anchors. The tested readout cost is not intrinsically above the matched ADC allowance. A full peripheral advantage remains unproved: normalization event cost, complete-dot throughput and physical precision/yield are unresolved. Equal-array cancellation is explicitly an assumption. The prior cell/PVT data is retained; q=4 and digital exponents remain frozen.

## Latest physical evidence — 2026-10-24

[SKY130 primitive results](sky130_phase_results.md) now supersede the prior next-hardware gate. Actual transmission gates, MIM storage, references, output OTA and resistor MAC have been tested, with one transistorized threshold. The known hybrid cost is about749 fJ/position, with five decisions and controller costs excluded; useful PVT/mismatch robustness is not established. The tested voltage interface needs redesign before crossbar integration. This is a negative result for this circuit mapping, not an arithmetic contradiction or proof against every representation. Canonical spacing stays25 mV nominal; exponents cannot repair analog offset or reference kickback.

Updated 2026-10-22 from the user-supplied architecture proposal. This document supersedes earlier fixed-bin interpretations. Proposed circuit behavior is distinguished from exact software arithmetic.

Current phase override (2026-10-23): freeze the lead at **q=4**, beta=2+sqrt(3), D=-3..3. q=3/5 remain comparison records. The new boundary is analog MAC -> bounded 4:1:1 normalization -> seven-region decision -> fresh canonical physical level -> next analog MAC. Classification alone does not restore a noncanonical capacitor. See [canonical contract](../spice/canonical_regeneration_spec.md) and [regeneration results](regeneration_phase_results.md).

Canonical output spacing is 25 mV at S=4; a weight-normalized MAC with conductance scale 3 has 8.333 mV raw spacing and scaled thresholds. Inputs are unchanged. All positions have comparable physical levels; significance comes from their indices. The reference write, output buffer and completion controller are explicit costs. Exponents do not repair comparator or reference errors. The next hardware gate is a realistic comparator/reference/write/driver primitive, before a full crossbar or FPGA.

## 1. Datatype and base family

Use

\[
X=\epsilon\,\beta_q^e M,\quad M=\sum_k d_k\beta_q^k,
\quad \epsilon\in\{-1,+1\},\quad
\beta_q=\frac{q+\sqrt{q^2-4}}2,\quad q\in\mathbb Z,\ q\ge3.
\]

Then beta^2=q beta-1 and

\[
q\beta^k=\beta^{k-1}+\beta^{k+1}.
\]

q=3,4,5 give expanding roots approximately 2.618, 3.732, 4.791. Gamma=phi^2 is the initial candidate, not a fixed architectural requirement. Exact verification reduces words to rational pairs A+B beta with the selected q, rather than evaluating floating-point beta.

The persistent q=3 alphabet is {-2,-1,0,1,2}. The generalized exact reference normalizes to D_q={-(q-1),...,q-1}; this is an explicit alphabet choice, not a proof that every base requires it. The fixed-five-state finite-precision study holds D constant when comparing bases. Larger beta need not improve useful precision if the chosen alphabet leaves representational gaps.

Physical positions carry beta^k significance. Local current/charge carries a small digit amplitude; exponentially scaled physical voltages are not required.

## 2. Three different scales

Use epsilon for sign, e for radix exponent, and S>0 for physical conductance/current scaling. The source proposal reused s for sign and analog scale; these are different quantities.

- Physical range: permitted current/charge at a node.
- Mantissa precision: representable values across a finite digit layout and alphabet.
- Floating dynamic range: exponent range, in addition to mantissa range.

Increasing e changes numerical magnitude without requiring a corresponding physical current increase. It does not repair mantissa quantization, comparator errors, or loss of digits during alignment. General S is not necessarily an integer power of beta and cannot silently be absorbed into e.

## 3. Weight-only amplitude normalization

For encoded inputs a[n,i] and weight coefficients b[n,j], raw channel coefficient k is

\[
c_k=\sum_n\sum_{i+j=k}a_{n,i}b_{n,j}.
\]

Program b/S and keep input digits unchanged. A sufficient common scale is

\[
S\ge\max_k\sum_n\sum_{i+j=k}|b_{n,j}|,\qquad |a_{n,i}|\le B.
\]

Then physical output I_k=c_k/S obeys |I_k|<=B. Signed weights require absolute-weight sums. With nonnegative weights summing to one this is a weighted average. The bound is B times normalized weight L1, not a universal bound of one for unbounded raw x. Current reference inputs are exact INT8 gamma encodings with B=2. Encoding does not normalize the value of x.

The scale must include every row and overlapping digit-product term at a node. Bounded partial sums need explicit handling when merged. Analog overshoot, mismatch, leakage and guard margin remain circuit questions.

## 4. Scale-aware equivalence bins

Known common scaling must move all detector centers and thresholds:

\[
I_d=d/S,\qquad t_{d,d+1}=(d+1/2)/S.
\]

Example: inputs [1,0,0,0,0], original weights all one, programmed weights all 0.2, S=5. Physical current 0.2 is the class d=1, rather than d=0. Ideal scaling alone introduces no error when the unscaled coefficient is one of the declared classes. Threshold ties in the reference choose the higher digit.

After this detector emits d, it represents the unscaled coefficient. Do not multiply the decoded gamma word by S again. This differs from the historical fixed-unscaled-bin control, which binned c/S and multiplied the resulting word by S.

Spacing is 1/S in normalized physical units. Absolute additive noise/offset therefore becomes more consequential as S increases. Relative gain/conductance errors have a different scaling law and are modeled separately. A configurable threshold generator and scale distribution have cost even if the comparator count is small.

### Supported-domain condition

Five scaled classes exactly recover coefficients in {-2,-1,0,1,2}. A bounded physical current does not guarantee that c_k lies in this set. For example, five input digits of 2 with weights 1/5 produce current 2 at S=5, corresponding to raw coefficient 10. A five-state detector saturates. A decoder cannot reconstruct discarded information from that emitted state alone.

Three distinct implementations must be kept separate:

1. Direct small-alphabet sensing, exact only on the declared local classes.
2. Expanded transient coefficient sensing, followed by normalization. It can recover bounded integer sums exactly, but may require many regions.
3. A local analog algebraic normalization network before/during final low-state regeneration. Its clocked capacitor/sample-hold/current-packet behavior is now tested in LTspice; transistor realization and energy remain unestablished. See normalizer_phase_results.md.

The suite tests all three as domain checks/benchmarks. The third uses an explicitly labeled software normalizer, never a hidden readout oracle. Arbitrary fractional learned coefficients are not guaranteed to lie on the integer digit lattice; boundary encoding/programming and approximation must be specified.

## 5. Local normalization and robust states

A unit rewrite subtracts q at position k and adds one at k-1 and k+1, preserving value exactly. Exact software supports bulk rewrites and counts equivalent unit events. It does not establish circuit latency, unit-transfer accuracy or threshold implementability.

Weight scaling removes physical range overflow under its assumptions. It does not remove every arithmetic need to redistribute coefficients into a small persistent alphabet. Unity additions and exponent-aligned block merges can also create excess coefficients. The updated proposal includes a normalization network for these cases; nominal normalization activity must be measured rather than assumed zero.

At fixed physical range [-Imax,Imax], N centers have spacing roughly 2 Imax/(N-1). Increasing the alphabet or transient class count eventually demands fine analog discrimination. The preferred design hypothesis is few robust states per node plus more positions. Global precision comes from the full word; its attainable accuracy must be measured.

## 6. Floating multiplication, dot products and renormalization

Multiplication gives sign epsilon1 epsilon2, exponent e1+e2, and mantissa product M1 M2. Any subsequent mantissa renormalization adds its positional shift to this exponent.

For a block with shared activation exponent eX and weight exponent eW,

\[
\sum_i X_iW_i=\beta^{e_X+e_W}\sum_i M_iN_i.
\]

The array computes the mantissa dot; the sidecar computes the exponent sum. Arbitrary per-scalar exponents need product alignment before summation; simply adding two shared exponents does not cover that case.

Renormalizing by r shifts digit indices k to k-r and changes e to e+r. This is exact on an unrestricted layout. A finite layout needs spare positions, rounding/guard rules or a reported discarded residual. Index shifts require actual routing/storage movement; they are not free merely because no general multiplier is needed.

Finite exponent storage also needs overflow/underflow policy. The reference uses unbounded Python integers and separately counts illustrative signed-eight-bit exponent overflow.

## 7. Addition and unity routing fabric

For different exponents choose e=max(eX,eY), shift each mantissa by eX-e or eY-e, then sum. A shift is a sparse permutation matrix P_delta, so aligned addition is

\[
\mathbf z=P_{\Delta_X}\mathbf x+P_{\Delta_Y}\mathbf y.
\]

Connections are 0,+1 or -1 for signed routing. Digital logic compares exponents, computes differences and selects routing. A candidate unity-weight CIM fabric performs shifts, residual adds, block accumulation and permutations. Truncating a shift on a finite layout is not a full permutation: it loses terms. Summing two valid digit vectors may require normalization.

Radix exponent alignment is positional. Differing arbitrary physical scales require calibrated gains or separate scale-aware detection, not solely a permutation.

## 8. Two-fabric pipeline and digital sidecar

Candidate pipeline:

Gamma mantissas -> learned-weight CIM -> local normalization / finite-region sensing -> canonical regeneration -> unity routing/summing -> normalization / sensing / regeneration -> next learned-weight stage. Count restoration at every required boundary, including after analog residual merges.

A clocked comparator/sample-hold/charge-transfer ordering is implemented behaviorally for 3/5/9 active nodes. The q=4 phase now adds seven-reference capacitor restoration whose actual output drives a subsequent weighted analog MAC, plus sixteen-stage controls. This is not a completed two-fabric chip. Reference generation, comparator/latch implementation, buffer/load/supply costs and asynchronous completion remain hardware gates. The wide-coefficient software detector stays a separate reference.

The digital sidecar stores exponents, compares/adds/subtracts them, selects shift configurations, and sequences state updates. Exponent addition could be explored in unity-weight CIM, but is not implemented or assumed advantageous. Small digital arithmetic is a proposed low-cost choice; total cost is not measured.

## 9. Nonlinear operations and validation boundary

Gamma-native candidates now cover mantissa multiply/dot, exponent scaling, alignment, residual adds, differently exponented block accumulation, and positional renormalization at the arithmetic-reference level. Their finite-layout and physical sensing implementations remain open.

GELU, SiLU, softmax, RMSNorm/LayerNorm, reciprocal/division, square root, and reliable comparison/max remain unresolved. Sign/clipping/ReLU may use thresholds, but a local signed digit does not alone determine the sign of a redundant full word. Do not claim gamma-native transformer execution without implementing those operations and counting binary conversions.

The feasibility test includes extra cells, array drive, sensing/threshold generation, normalization, routing, exponent/scale handling, calibration, and boundary conversion, compared at matched useful accuracy/throughput against array plus ADC and binary accumulation. No current software metric establishes a GPU or energy advantage.

## 10. Current local-circuit evidence

See [regeneration_phase_results.md](regeneration_phase_results.md) for the current q=4 connected-circuit tests, canonical levels, depth/stress failures, exact finite-word policy, completion model and all-in energy budgets. q=4 is frozen for this gate, not proven globally optimal. The prior [normalizer phase](normalizer_phase_results.md) retains 108 ideal and 2,376 nonideal LTspice cases and q=3/5 controls. Digital exponent/block-floating references remain in the regression suite.
