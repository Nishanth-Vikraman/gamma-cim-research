"""Bounded coefficient MAC reference; no post-MAC toppling or numerical ADC.

Exact rational arithmetic is an offline verification oracle. The nearest-five
detector below is the HISTORICAL fixed-unscaled-bin control. Architecture v2
uses scale-aware sensing in gamma_native_model; retain this control for paired
comparisons, not as the intended current detector.
Integer-to-gamma boundary encoding remains in gamma_exact_model.
"""
from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Mapping, Sequence

from gamma_exact_model import gamma_pow_pair

Scalar = int | Fraction
Word = Mapping[int, Scalar]
DIGITS = (-2, -1, 0, 1, 2)
THRESHOLDS = tuple(Fraction(a + b, 2) for a, b in zip(DIGITS, DIGITS[1:]))


def rational(value: Scalar) -> Fraction:
    if not isinstance(value, (int, Fraction)):
        raise TypeError('exact reference accepts int or Fraction, not float')
    return Fraction(value)


def exact_pair(word: Word, scale: Scalar = 1) -> tuple[Fraction, Fraction]:
    """Offline exact value oracle; never used by the detector."""
    a = b = Fraction(0)
    for k, coefficient in word.items():
        pa, pb = gamma_pow_pair(k)
        coefficient = rational(coefficient)
        a += coefficient * pa
        b += coefficient * pb
    scale = rational(scale)
    return scale * a, scale * b


@dataclass
class MacPlan:
    # Programmed, amplitude-normalized weights; scale restores original units.
    weights: tuple[dict[int, Fraction], ...]
    input_positions: tuple[int, ...]
    output_positions: tuple[int, ...]
    input_bound: Fraction
    utilization: Fraction
    scale: Fraction
    channel_l1: dict[int, Fraction]


def make_plan(
    weights: Sequence[Word],
    input_positions: Sequence[int],
    *,
    weight_positions: Sequence[int] | None = None,
    input_bound: Scalar = 2,
    utilization: Scalar = 1,
) -> MacPlan:
    """Compile weight-only common scale for a frozen physical word layout.

No observed input value or raw MAC output participates in scale selection.
Utilization reserves amplitude headroom, not extra detector levels. Every
row/digit product feeding a given coefficient channel contributes to its bound.
"""
    positions = tuple(sorted(set(input_positions)))
    wpositions = positions if weight_positions is None else tuple(sorted(set(weight_positions)))
    if not positions or not wpositions:
        raise ValueError('physical layouts must be nonempty')
    bound, utilization = rational(input_bound), rational(utilization)
    if bound <= 0 or not 0 < utilization <= 1:
        raise ValueError('input_bound must be positive; utilization must be in (0,1]')
    converted = []
    for word in weights:
        if any(k not in wpositions for k in word):
            raise ValueError('weight outside frozen weight layout')
        converted.append({k: rational(v) for k, v in word.items() if v})
    outputs = tuple(sorted({i + j for i in positions for j in wpositions}))
    l1 = {k: Fraction(0) for k in outputs}
    for word in converted:
        for i in positions:
            for j, weight in word.items():
                l1[i + j] += abs(weight)
    scale = max(Fraction(1), max(l1.values(), default=Fraction(0)) / utilization)
    normalized = tuple({k: value / scale for k, value in word.items()} for word in converted)
    return MacPlan(normalized, positions, outputs, bound, utilization, scale, l1)


@dataclass
class MacResult:
    coefficients: dict[int, Fraction]
    scale: Fraction
    peak_final: Fraction
    peak_partial: Fraction
    digit_products: int


def bounded_mac(plan: MacPlan, inputs: Sequence[Word]) -> MacResult:
    """Compute scaled analog coefficient sums; no digit snapping or topplings.

The partial peak follows one deterministic contribution order. The weight-only
L1 envelope guarantees the bound for every ideal ordering, without modeling
physical pulse overshoot or mismatch.
"""
    if len(inputs) != len(plan.weights):
        raise ValueError('input and weight row counts differ')
    acc = defaultdict(Fraction)
    peak = Fraction(0)
    products = 0
    for word, weights in zip(inputs, plan.weights):
        if any(k not in plan.input_positions for k in word):
            raise ValueError('input outside frozen input layout')
        for i, value in word.items():
            value = rational(value)
            if abs(value) > plan.input_bound:
                raise ValueError('input digit exceeds declared bound')
            if not value:
                continue
            for j, weight in weights.items():
                acc[i + j] += value * weight
                peak = max(peak, abs(acc[i + j]))
                products += 1
    coefficients = {k: acc[k] for k in plan.output_positions}
    return MacResult(
        coefficients, plan.scale,
        max((abs(v) for v in coefficients.values()), default=Fraction(0)),
        peak, products,
    )


def detect_five(
    coefficients: Word,
    positions: Sequence[int],
    *,
    threshold_offset: Scalar = 0,
    observation_errors: Word | None = None,
) -> dict[int, int]:
    """Candidate four-comparator detector, fixed levels -2,-1,0,1,2.

Threshold ties select the higher digit. Offsets and observation errors use
digit-amplitude units. This rule can lose value; it is not asserted to implement
the final algebraic equivalence regions. It has no access to exact_pair, raw
inputs, the desired answer, or a post-MAC adaptive scale.
"""
    offset = rational(threshold_offset)
    thresholds = tuple(t + offset for t in THRESHOLDS)
    errors = observation_errors or {}
    return {
        k: DIGITS[bisect_right(thresholds, rational(coefficients.get(k, 0)) + rational(errors.get(k, 0)))]
        for k in positions
    }
