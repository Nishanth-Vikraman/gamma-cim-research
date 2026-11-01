"""Exactness gates and adversarial limits for architecture v2."""
from decimal import Decimal
from fractions import Fraction as F
import itertools
import random

import pytest

from gamma_native_model import (NativeFloat, add, block_dot, block_encode, convolve,
    multiply, normalize, pair_multiply, power_pair, route_align, sense,
    sense_threshold_vector, value_pair)


@pytest.mark.parametrize('scale', [F(1, 7), F(1), F(5), F(29, 3), F(10000)])
def test_scale_aware_five_bins_recover_every_legal_digit(scale):
    for d in range(-2, 3):
        assert sense(F(d)/scale, scale).digit == d
        assert not sense(F(d)/scale, scale).outside_centers


def test_scale_aware_example_and_tie():
    assert sense(F(1, 5), 5).digit == 1
    assert sense(F(1, 10), 5).digit == 1
    assert sense(F(1, 10), 5, offset=F(1, 1000)).digit == 0
    for bad in [0, -1]:
        with pytest.raises(ValueError):
            sense(0, bad)


def test_current_headroom_does_not_imply_raw_coefficient_in_alphabet():
    # Five legal x digits, weights 1/5: physical sum 2 is bounded,
    # but represents raw coefficient 10, which is not a legal local digit.
    assert sense(F(2), 5).digit == 2
    assert sense(F(2), 5).outside_centers
    assert sense(F(2), 5, radius=10).digit == 10
    word, _ = normalize({0: 10})
    assert value_pair(word) == (10, 0)
    assert all(abs(v) <= 2 for v in word.values())


def test_smaller_spacing_increases_sensitivity_to_same_physical_error():
    assert sense(F(1, 5)+F(1, 100), 5).digit == 1
    assert sense(F(1, 100)+F(1, 100), 100).digit == 2


def test_explicit_comparator_bank_and_inverted_threshold_rejection():
    for s in [1, 5, 20]:
        for d in range(-2, 3):
            assert sense_threshold_vector(F(d, s), s, 2, [0]*4) == d
    with pytest.raises(ValueError):
        sense_threshold_vector(0, 1, 2, [0, 2, -2, 0])


@pytest.mark.parametrize('q', [3, 4, 5, 6])
def test_family_identities_integer_encoding_products_and_bound(q):
    for k in range(-5, 6):
        assert value_pair({k: q}, q) == value_pair({k-1: 1, k+1: 1}, q)
        assert pair_multiply(power_pair(q, k), power_pair(q, -k), q) == (1, 0)
    for n in range(-128, 128):
        word, _ = normalize({0: n}, q)
        assert value_pair(word, q) == (n, 0)
        assert all(abs(d) < q for d in word.values())
    rng = random.Random(q)
    for _ in range(40):
        a, b = rng.randint(-128, 127), rng.randint(-128, 127)
        da, _ = normalize({0: a}, q)
        db, _ = normalize({0: b}, q)
        out, stats = normalize(convolve(da, db), q)
        assert value_pair(out, q) == (a*b, 0)
        assert stats['unit_events'] >= 0


@pytest.mark.parametrize('q', [3, 4, 5])
def test_native_signed_multiply_add_and_renormalize_exactly(q):
    rng = random.Random(83)
    for _ in range(50):
        a = NativeFloat({-2: 1, 0: rng.choice([-2, -1, 1, 2])}, rng.randint(-10, 10), rng.choice([-1, 1]), q)
        b = NativeFloat({-1: -1, 0: 2}, rng.randint(-10, 10), rng.choice([-1, 1]), q)
        product, _ = multiply(a, b)
        assert product.exact() == pair_multiply(a.exact(), b.exact(), q)
        result, _ = add(a, b)
        assert result.exact() == tuple(x+y for x, y in zip(a.exact(), b.exact()))
        assert a.renormalize(7).exact() == a.exact()
        assert a.renormalize(-7).exact() == a.exact()


def test_sparse_permutation_alignment_and_finite_layout_loss():
    e, raw, edges, dropped = route_align([{0: 2}, {0: 1}], [3, -7])
    assert e == 3 and raw == {0: 2, -10: 1}
    assert edges == [(0, 0, 0, 1), (1, 0, -10, 1)] and not dropped
    _, kept, _, lost = route_align([{0: 2}, {0: 1}], [3, -7], window=(-4, 0))
    assert kept == {0: 2} and lost == {-10: 1}


def test_unity_addition_requires_normalization_and_sign_cancellation():
    a = NativeFloat({0: 2})
    out, stats = add(a, a)
    assert out.exact() == (4, 0)
    assert stats['pre_normalization_peak'] == 4 and stats['unit_events'] > 0
    zero, _ = add(a, NativeFloat({0: 2}, sign=-1))
    assert zero.exact() == (0, 0)


def test_block_boundary_encoder_precision_and_shared_exponent_tradeoff():
    values = [Decimal(1), Decimal('0.000001'), Decimal(-1)]
    _, shared_e, shared = block_encode(values, positions=8)
    _, individual_e, individual = block_encode([values[1]], positions=8)
    assert shared_e != individual_e
    assert abs(individual[0]-values[1]) < abs(shared[1]-values[1])
    for p in [4, 8, 12]:
        _, _, decoded = block_encode([Decimal('0.314159')], positions=p)
        if p > 4:
            assert abs(decoded[0]-Decimal('0.314159')) < previous
        previous = abs(decoded[0]-Decimal('0.314159'))
    assert block_encode([Decimal(0)])[2] == [0]
    with pytest.raises(ValueError):
        block_encode([Decimal('NaN')])


def test_incompatible_base_and_illegal_persistent_digit_rejected():
    with pytest.raises(ValueError):
        add(NativeFloat({0: 1}, q=3), NativeFloat({0: 1}, q=4))
    with pytest.raises(ValueError):
        NativeFloat({0: 3}, q=3)
    with pytest.raises(TypeError):
        normalize({0: F(1, 2)})


def test_small_codebook_contains_real_redundancy():
    pairs = [value_pair(dict(enumerate(digits)))
             for digits in itertools.product(range(-2, 3), repeat=3)]
    assert len(set(pairs)) < len(pairs)


def test_block_dot_and_second_stage_preserve_actual_output_digits():
    xs, ws = [{0: 2}, {-1: 1}], [{0: 1}, {0: -2}]
    result, _ = block_dot(xs, ws, 4, -2)
    expected = [F(0), F(0)]
    for x, w in zip(xs, ws):
        term = pair_multiply(NativeFloat(x, 4).exact(), NativeFloat(w, -2).exact())
        expected = [a+b for a, b in zip(expected, term)]
    assert result.exact() == tuple(expected)
    stage2, _ = block_dot([result.digits], [{0: 2}], result.exponent, 3)
    assert stage2.exact() == pair_multiply(result.exact(), NativeFloat({0: 2}, 3).exact())
    with pytest.raises(ValueError):
        block_dot(xs, ws[:1], 0, 0)


def test_finite_conversion_exactness_uses_algebraic_pairs_for_full_int8():
    # Avoid attributing 1e-88 Decimal cancellation residue to representation
    # error, or presenting its inverse as hundreds of achieved bits.
    for n in range(-128, 128):
        words, exponent, _ = block_encode([Decimal(n)], positions=12)
        assert NativeFloat(words[0], exponent).exact() == (n, 0)


def test_binary_target_samples_preserve_extreme_values_and_signed_symmetry():
    from run_native_experiments import format_values
    formats = format_values(0)
    assert len(formats['FP8_E4M3FN']) == 254
    assert max(formats['FP8_E4M3FN']) == 448
    assert min(v for v in formats['FP8_E4M3FN'] if v > 0) == Decimal(2)**-9
    for name in ['BF16', 'FP16']:
        values = formats[name]
        assert all(F(values[i]) == -F(values[i+1]) for i in range(0, len(values), 2))
    assert F(max(formats['BF16'])) == F(255, 128)*2**127
    assert min(v for v in formats['FP16'] if v > 0) == Decimal(2)**-24
