# Traceability to the supplied architecture update

Source: user-pasted “Gamma-CIM Architecture Update”, 2026-10-22. This is a structured incorporation record, not a verbatim archive. The current equations and qualifications are in architecture_summary.md.

| Proposal sections | Project treatment |
| --- | --- |
| 1, 16: scale-aware regions | Implemented scaled centers/thresholds, exact legal-class gate, raw-MAC saturation control and broad transient reference; five error mechanisms swept |
| 2: generalized base | q=3..6 exact field arithmetic and local rewrites; both growing-alphabet exact encoding and fixed-five-state codebooks |
| 3: small robust alphabets | Added alphabet/spacing/headroom cost gates; no unlimited alphabet precision claim |
| 4–5: gamma floating point | NativeFloat stores signed mantissa, exponent and q; exact multiply and exponent sum |
| 6: block floating point | Tensor/channel/tile/small-block encoding, mantissa dots, exponented block accumulation and metadata counts |
| 7–8: shifts and additions | Exact positional renormalization/alignment, signed add, finite-window dropped terms |
| 9–11: two fabrics and sidecar | Sparse signed unity edge lists, active connection and width metrics; digital metadata remains the default proposal |
| 12: gamma-native operations | Demonstrated at exact reference level; circuit and finite-precision costs remain gates |
| 13: nonlinear functions | Explicitly unresolved; transformer conversion-frequency audit remains future work |
| 14–15: physical picture and thesis | Integrated sign/exponent/algebraic mantissa datatype and sensing/normalization loop |
| 16: effective precision | Exact-pair deduplicated codebooks; sampled step/error metrics; INT8/FP8 exhaustive finite inputs, BF16/FP16 seeded samples and extremes |
| 17: feasibility inequality | Includes threshold generation, scale/exponent storage, all normalization/routing, and matched-accuracy comparison |

## Necessary qualifications from implementation

- 100% ideal recovery is required and achieved for declared integer classes at arbitrary positive scales. It cannot be required from five emitted symbols for arbitrary raw coefficients outside those five classes without an additional mechanism.
- The floating exponent separates global magnitude from physical amplitude. It does not remove mantissa errors, alignment truncation, noise sensitivity or exponent overflow.
- q-family experiments vary the alphabet explicitly; larger bases with fixed five-state alphabets can leave gaps.
- The normalizer is an exact software reference. It cannot be presented as a demonstrated analog circuit or free preprocessing before sensing.
- A general physical scale S differs from a radix exponent and may require more than positional routing when combined across branches.

See validation_status.md for the executed results rather than interpreting the original proposal's intended outcomes as established results.

## Subsequent implementation

The 2026-10-23 local engineering phase now supplies behavioral SPICE and physical-stress evidence: see [normalizer_phase_results.md](normalizer_phase_results.md). Its detailed limitations replace any earlier assumption that a software normalizer alone establishes circuit feasibility.
