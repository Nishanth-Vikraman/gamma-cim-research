from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
import math
from typing import Dict, Iterable, Tuple

DigitMap = Dict[int, int]


def gamma_pow_pair(k: int) -> Tuple[int, int]:
    """Return exact (A,B) with gamma^k = A + B*gamma, gamma^2=3gamma-1."""
    A, B = 1, 0
    if k >= 0:
        for _ in range(k):
            A, B = -B, A + 3 * B
    else:
        # gamma^-1 = 3 - gamma
        for _ in range(-k):
            A, B = 3 * A + B, -A
    return A, B


def exact_value(digits: DigitMap) -> Tuple[int, int]:
    """Evaluate sum_k d[k]*gamma^k exactly as (A,B) = A+B*gamma."""
    A = B = 0
    for k, c in digits.items():
        a, b = gamma_pow_pair(k)
        A += c * a
        B += c * b
    return A, B


def trim(d: DigitMap) -> DigitMap:
    return {int(k): int(v) for k, v in d.items() if int(v) != 0}


@dataclass
class NormalizeStats:
    topplings: int
    rounds: int
    peak_abs_coeff: int


def normalize(digits: DigitMap, synchronous: bool = True, max_rounds: int = 100000):
    """
    Normalize using 3*gamma^k = gamma^(k-1)+gamma^(k+1).
    Output digits are in {-2,-1,0,1,2}.
    """
    c = defaultdict(int, trim(digits))
    before = exact_value(dict(c))
    peak = max((abs(v) for v in c.values()), default=0)
    topplings = 0
    rounds = 0

    if synchronous:
        while True:
            qmap = {k: math.trunc(v / 3) for k, v in list(c.items()) if abs(v) >= 3}
            qmap = {k: q for k, q in qmap.items() if q}
            if not qmap:
                break
            delta = defaultdict(int)
            for k, q in qmap.items():
                delta[k] -= 3 * q
                delta[k - 1] += q
                delta[k + 1] += q
                topplings += abs(q)
            for k, dv in delta.items():
                c[k] += dv
                peak = max(peak, abs(c[k]))
                if c[k] == 0:
                    del c[k]
            rounds += 1
            if rounds > max_rounds:
                raise RuntimeError('normalization did not converge')
    else:
        while True:
            bad = [k for k, v in c.items() if abs(v) >= 3]
            if not bad:
                break
            k = max(bad, key=lambda kk: (abs(c[kk]), -kk))
            q = math.trunc(c[k] / 3)
            c[k] -= 3 * q
            c[k - 1] += q
            c[k + 1] += q
            topplings += abs(q)
            peak = max(peak, abs(c[k]), abs(c[k - 1]), abs(c[k + 1]))
            if c[k] == 0:
                del c[k]
            rounds += 1
            if rounds > max_rounds:
                raise RuntimeError('normalization did not converge')

    out = trim(dict(c))
    assert all(abs(v) <= 2 for v in out.values())
    assert exact_value(out) == before
    return out, NormalizeStats(topplings, rounds, peak)


def encode_int(n: int, synchronous: bool = True) -> DigitMap:
    """Baseline exact encoder. Deliberately simple, not assumed optimal."""
    out, _ = normalize({0: int(n)}, synchronous=synchronous)
    assert exact_value(out) == (int(n), 0)
    return out


def convolve(a: DigitMap, b: DigitMap) -> DigitMap:
    out = defaultdict(int)
    for i, ai in a.items():
        for j, bj in b.items():
            out[i + j] += ai * bj
    return trim(dict(out))


def multiply_ints(a: int, b: int, synchronous: bool = True):
    da = encode_int(a, synchronous)
    db = encode_int(b, synchronous)
    raw = convolve(da, db)
    result, stats = normalize(raw, synchronous=synchronous)
    assert exact_value(result) == (a * b, 0)
    return result, stats, raw, da, db


def span(d: DigitMap) -> int:
    return 0 if not d else max(d) - min(d) + 1


def nnz(d: DigitMap) -> int:
    return sum(v != 0 for v in d.values())


@dataclass
class DotStats:
    topplings: int
    rounds: int
    topplings_per_mac: float
    rounds_per_mac: float
    peak_raw_coeff: int
    final_span: int
    final_nnz: int


def dot_product_ints(xs: Iterable[int], ws: Iterable[int], normalize_every: int = 1, synchronous: bool = True):
    xs, ws = list(xs), list(ws)
    if len(xs) != len(ws):
        raise ValueError('xs and ws must have equal length')
    if not xs:
        return {}, DotStats(0,0,0.0,0.0,0,0,0)
    if normalize_every < 1:
        raise ValueError('normalize_every must be >= 1')

    acc = defaultdict(int)
    total_topplings = total_rounds = peak_raw = 0

    for idx, (x, w) in enumerate(zip(xs, ws), 1):
        dx = encode_int(x, synchronous)
        dw = encode_int(w, synchronous)
        prod = convolve(dx, dw)
        for k, v in prod.items():
            acc[k] += v
        peak_raw = max(peak_raw, max((abs(v) for v in acc.values()), default=0))
        if idx % normalize_every == 0:
            norm, st = normalize(dict(acc), synchronous=synchronous)
            acc = defaultdict(int, norm)
            total_topplings += st.topplings
            total_rounds += st.rounds

    if len(xs) % normalize_every:
        norm, st = normalize(dict(acc), synchronous=synchronous)
        acc = defaultdict(int, norm)
        total_topplings += st.topplings
        total_rounds += st.rounds

    result = trim(dict(acc))
    expected = sum(x*w for x, w in zip(xs, ws))
    assert exact_value(result) == (expected, 0)
    n = len(xs)
    return result, DotStats(
        total_topplings,
        total_rounds,
        total_topplings/n,
        total_rounds/n,
        peak_raw,
        span(result),
        nnz(result),
    )
