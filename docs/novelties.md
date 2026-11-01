# Novelty hypotheses

## Physical evidence limit — 2026-10-24

[The SKY130 boundary experiment](sky130_phase_results.md) is unfavorable for the tested voltage-reference/OTA implementation. Exact algebra and behavioral restoration do not establish a low-energy physical advantage. Neither numerical representational precision nor the limited seven-reference bank should be described as demonstrated silicon efficiency; the tested hybrid costs about749 fJ/position before full detector/control implementation and has mismatch failures.

These are research hypotheses, not patentability or performance conclusions.

Current lead q=4 adds an explicit regenerative boundary: classification -> finite reference write -> following weighted analog MAC. See [regeneration results](regeneration_phase_results.md). A seven-entry selector is permitted; a full-precision intermediate ADC is not. Reference, driver and control costs remain part of the hypothesis.

1. Co-design a redundant algebraic floating datatype with analog CIM, scale-aware coarse sensing and local digit regeneration. The representation remains native across useful stages.
2. Explore the family beta^2-q beta+1=0 jointly with digit alphabet, sensing margins and local neighbor topology. q=3 is the first candidate.
3. Separate bounded physical mantissas from global radix exponent range, using block exponents and exact positional alignment/renormalization.
4. Pair a learned-weight compute fabric with a sparse unity-weight routing/summing fabric under exponent metadata control.
5. Obtain useful global precision from multiple robust local digit states rather than indefinitely increasing per-node analog resolution.

Weight scaling, floating exponents, threshold scaling and permutation matrices are supporting mechanisms, not standalone novelty findings. The hypothesized advantage depends on a low-cost sensing/normalization network. Exact recovery using thousands of transient classes is not evidence of a coarse-readout advantage.
