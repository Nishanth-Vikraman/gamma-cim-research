"""Architecture-v2 reference. Exact field arithmetic is an offline oracle.

Scale-aware sensing operates only on physical observations and known scale.
The normalizer is an ideal software reference, NOT an implemented analog loop.
Decimal boundary encoders are offline conversion studies, not readout logic.
"""
from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal, ROUND_HALF_UP, localcontext
from fractions import Fraction as F
from functools import lru_cache
import math


@lru_cache(maxsize=None)
def power_pair(q: int, k: int) -> tuple[int, int]:
    if not isinstance(q, int) or q < 3:
        raise ValueError('integer q >= 3 required')
    a, b = 1, 0
    for _ in range(abs(k)):
        a, b = (-b, a + q*b) if k >= 0 else (q*a + b, -a)
    return a, b


def value_pair(word, q=3):
    a = b = F(0)
    for k, digit in word.items():
        pa, pb = power_pair(q, k)
        a += digit*pa
        b += digit*pb
    return a, b


def pair_multiply(x, y, q=3):
    a, b = x
    c, d = y
    return a*c-b*d, a*d+b*c+q*b*d


def convolve(a, b):
    out = defaultdict(int)
    for i, x in a.items():
        for j, y in b.items():
            out[i+j] += x*y
    return {k: v for k, v in out.items() if v}


def normalize(word, q=3, max_rounds=10000):
    """Bulk exact benchmark into D=[-(q-1),q-1]; no physical-cost claim.

Each unit event removes q and adds 1 at each neighbor. L1 falls by at
least q-2 per unit event, giving termination for finite integer words.
"""
    power_pair(q, 0)
    if any(not isinstance(v, int) for v in word.values()):
        raise TypeError('normalizer requires integer coefficients')
    c = defaultdict(int, {k: v for k, v in word.items() if v})
    initial = value_pair(c, q)
    events = rounds = 0
    peak = max(map(abs, c.values()), default=0)
    while any(abs(v) >= q for v in c.values()):
        delta = defaultdict(int)
        for k, v in list(c.items()):
            t = (abs(v)//q) * (1 if v >= 0 else -1)
            if t:
                delta[k] -= q*t
                delta[k-1] += t
                delta[k+1] += t
                events += abs(t)
        for k, v in delta.items():
            c[k] += v
        c = defaultdict(int, {k: v for k, v in c.items() if v})
        peak = max(peak, max(map(abs, c.values()), default=0))
        rounds += 1
        if rounds > max_rounds:
            raise RuntimeError('normalization round limit')
    assert value_pair(c, q) == initial
    return dict(c), {'unit_events': events, 'bulk_rounds': rounds,
                     'peak_coefficient': peak}


@dataclass(frozen=True)
class RegionResult:
    digit: int
    outside_centers: bool


def sense(current, scale, radius=2, offset=0):
    """Known scale moves centers to d/scale, thresholds to (d+1/2)/scale.

No oracle or expected coefficient is available. Output saturates at +/-radius;
outside_centers is diagnostic instrumentation, not an extra readout symbol.
Ties select the higher digit. All parameters are in physical current units.
"""
    current, scale, offset = F(current), F(scale), F(offset)
    if scale <= 0 or not isinstance(radius, int) or radius < 1:
        raise ValueError('positive scale and integer radius >= 1 required')
    # Equivalent comparator count: 2*radius. Arithmetic is software emulation.
    index = max(-radius, min(radius, math.floor(scale*(current-offset)+F(1, 2))))
    return RegionResult(index, abs(current-offset) > F(radius)/scale)


def sense_threshold_vector(current, scale, radius, offsets):
    """Explicit ordered comparator bank, used for threshold mismatch studies."""
    if F(scale) <= 0 or radius < 1 or len(offsets) != 2*radius:
        raise ValueError('invalid comparator-bank geometry')
    thresholds = [F(2*d+1, 2)/F(scale)+F(o)
                  for d, o in zip(range(-radius, radius), offsets)]
    if any(a >= b for a, b in zip(thresholds, thresholds[1:])):
        raise ValueError('threshold order inversion')
    return bisect_right(thresholds, F(current))-radius


def shift(word, amount):
    return {k+amount: v for k, v in word.items() if v}


@dataclass(frozen=True)
class NativeFloat:
    digits: dict[int, int]
    exponent: int = 0
    sign: int = 1
    q: int = 3

    def __post_init__(self):
        power_pair(self.q, 0)
        if self.sign not in (-1, 1) or not isinstance(self.exponent, int):
            raise ValueError('sign +/-1 and integer exponent required')
        if any(not isinstance(v, int) or abs(v) >= self.q for v in self.digits.values()):
            raise ValueError('persistent digits outside family alphabet')

    def exact(self):
        return tuple(self.sign*v for v in value_pair(shift(self.digits, self.exponent), self.q))

    def renormalize(self, amount):
        # Value preserving on an unbounded layout; finite routing must count loss.
        return NativeFloat(shift(self.digits, -amount), self.exponent+amount, self.sign, self.q)


def multiply(x: NativeFloat, y: NativeFloat):
    if x.q != y.q:
        raise ValueError('mixed bases require conversion')
    digits, stats = normalize(convolve(x.digits, y.digits), x.q)
    return NativeFloat(digits, x.exponent+y.exponent, x.sign*y.sign, x.q), stats


def block_dot(inputs, weights, input_exponent, weight_exponent, q=3):
    """Exact shared-exponent mantissa MAC with ideal reference normalization.

No boundary re-encoding of an exact answer. Raw convolution coefficients are
the actual input to normalization. Physical sensing is tested separately.
"""
    if len(inputs) != len(weights):
        raise ValueError('block row counts differ')
    raw = defaultdict(int)
    products = 0
    for x, w in zip(inputs, weights):
        NativeFloat(x, q=q)
        NativeFloat(w, q=q)
        for k, v in convolve(x, w).items():
            raw[k] += v
        products += len(x)*len(w)
    digits, stats = normalize(raw, q)
    stats['digit_products'] = products
    return NativeFloat(digits, input_exponent+weight_exponent, q=q), stats


def route_align(words, exponents, signs=None, window=None):
    """Sparse signed unity connections: (source,source_k,destination_k,gain).

Window is an optional inclusive finite output layout. Dropped terms are
returned explicitly. No assumption that adding legal digits stays legal.
"""
    if len(words) != len(exponents) or not words:
        raise ValueError('nonempty matching words/exponents required')
    signs = [1]*len(words) if signs is None else signs
    if len(signs) != len(words) or any(s not in (-1, 1) for s in signs):
        raise ValueError('invalid routing signs')
    common = max(exponents)
    acc, dropped = defaultdict(int), defaultdict(int)
    edges = []
    for source, (word, exponent, sign) in enumerate(zip(words, exponents, signs)):
        for k, value in word.items():
            if not value:
                continue
            dest = k+exponent-common
            if window and not window[0] <= dest <= window[1]:
                dropped[dest] += sign*value
            else:
                acc[dest] += sign*value
                edges.append((source, k, dest, sign))
    return common, dict(acc), edges, dict(dropped)


def add(x: NativeFloat, y: NativeFloat):
    if x.q != y.q:
        raise ValueError('mixed bases require conversion')
    exponent, raw, edges, dropped = route_align(
        [x.digits, y.digits], [x.exponent, y.exponent], [x.sign, y.sign])
    digits, stats = normalize(raw, x.q)
    stats.update(active_connections=len(edges), dropped_terms=len(dropped),
                 routing_width=max(raw, default=0)-min(raw, default=0)+1,
                 pre_normalization_peak=max(map(abs, raw.values()), default=0))
    return NativeFloat(digits, exponent, 1, x.q), stats


@lru_cache(maxsize=None)
def beta_decimal(q):
    power_pair(q, 0)
    with localcontext() as ctx:
        ctx.prec = 100
        return (Decimal(q)+Decimal(q*q-4).sqrt())/2


def decimal_value(word, q):
    with localcontext() as ctx:
        ctx.prec = 90
        beta = beta_decimal(q)
        return sum((Decimal(v)*beta**k for k, v in word.items()), Decimal(0))


def block_encode(values, q=3, positions=8, radius=2):
    """Offline greedy finite mantissa encoder, shared exponent per call.

Mantissa positions [1-positions,0]. No claim of optimal representation;
conversion errors are measured. Input scaling here is FP boundary encoding,
separate from weight-only physical conductance normalization in a MAC.
"""
    if positions < 1 or radius < 1:
        raise ValueError('positive positions/radius required')
    with localcontext() as ctx:
        ctx.prec = 90
        vals = [Decimal(v) for v in values]
        if any(not v.is_finite() for v in vals):
            raise ValueError('finite values only')
        beta = beta_decimal(q)
        peak = max(map(abs, vals), default=Decimal(0))
        exponent = 0
        # Ensure max normalized amplitude in (1/beta,1], except all-zero block.
        while peak > beta**exponent:
            exponent += 1
        while peak and peak <= beta**(exponent-1):
            exponent -= 1
        words, decoded = [], []
        for v in vals:
            residual = v/(beta**exponent)
            word = {}
            for k in range(0, -positions, -1):
                place = beta**k
                digit = max(-radius, min(radius, int((residual/place).to_integral_value(rounding=ROUND_HALF_UP))))
                if digit:
                    word[k] = digit
                residual -= digit*place
            words.append(word)
            decoded.append(decimal_value(word, q)*beta**exponent)
        return words, exponent, decoded
