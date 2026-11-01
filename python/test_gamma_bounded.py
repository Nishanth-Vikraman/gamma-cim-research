"""Physical bound tests plus historical fixed-unscaled-bin diagnostics.

Updated scale-aware and native datatype gates are in test_gamma_native.py.
"""
from fractions import Fraction as F
import itertools
import random

import pytest

from gamma_bounded_model import bounded_mac, detect_five, exact_pair, make_plan
from gamma_exact_model import encode_int, exact_value


def test_nonnegative_weights_are_programmed_as_a_weighted_average():
    plan = make_plan([{0: 1}, {0: 2}, {0: 3}], [0])
    assert [word[0] for word in plan.weights] == [F(1, 6), F(1, 3), F(1, 2)]
    result = bounded_mac(plan, [{0: -2}, {}, {0: 2}])
    assert result.coefficients == {0: F(2, 3)}
    assert exact_pair(result.coefficients, result.scale) == (4, 0)


def test_signed_sum_unity_is_not_a_range_bound():
    # 2 + (-1) = 1, but aligned signed inputs give 2*2 + (-2)*(-1) = 6.
    plan = make_plan([{0: 2}, {0: -1}], [0])
    result = bounded_mac(plan, [{0: 2}, {0: -2}])
    assert plan.scale == 3
    assert result.coefficients == {0: F(2)}
    assert exact_pair(result.coefficients, result.scale) == (6, 0)


def test_all_small_signed_inputs_and_every_partial_order_are_bounded():
    weights = [2, -1, 2]
    plan = make_plan([{0: w} for w in weights], [0], utilization=F(3, 4))
    for xs in itertools.product(range(-2, 3), repeat=3):
        result = bounded_mac(plan, [{0: x} for x in xs])
        assert result.peak_final <= F(3, 2)
        assert result.peak_partial <= F(3, 2)
        assert exact_pair(result.coefficients, result.scale) == (sum(x*w for x, w in zip(xs, weights)), 0)
        for order in itertools.permutations(range(3)):
            partial = F(0)
            for row in order:
                partial += F(xs[row] * weights[row], plan.scale)
                assert abs(partial) <= F(3, 2)


def test_full_convolution_bound_includes_overlapping_digit_products():
    word = {-1: 2, 0: 2, 1: 2}
    plan = make_plan([word, word], [-1, 0, 1])
    result = bounded_mac(plan, [word, word])
    # Central channel receives 2 rows * 3 weight contributions, each magnitude 2.
    assert plan.channel_l1[0] == 12
    assert plan.scale == 12
    assert result.coefficients[0] == 2
    assert result.peak_partial <= 2


def test_full_int8_domain_with_frozen_weights_layout_and_scale():
    words = {n: encode_int(n) for n in range(-128, 128)}
    layout = sorted({k for word in words.values() for k in word})
    plan = make_plan([words[127]], layout)
    saved_scale, saved_positions = plan.scale, plan.output_positions
    for n, word in words.items():
        result = bounded_mac(plan, [word])
        assert result.peak_final <= 2 and result.peak_partial <= 2
        assert exact_pair(result.coefficients, result.scale) == (127*n, 0)
        assert plan.scale == saved_scale and tuple(result.coefficients) == saved_positions


def test_frozen_random_int8_dots_restore_exactly_before_sensing():
    rng = random.Random(7)
    layout = list(range(-4, 5))
    for length in [1, 16, 64, 256]:
        weights = [rng.randint(-128, 127) for _ in range(length)]
        plan = make_plan([encode_int(w) for w in weights], layout)
        for _ in range(3):
            xs = [rng.randint(-128, 127) for _ in range(length)]
            result = bounded_mac(plan, [encode_int(x) for x in xs])
            assert result.peak_final <= 2 and result.peak_partial <= 2
            assert exact_pair(result.coefficients, result.scale) == (sum(x*w for x, w in zip(xs, weights)), 0)


def test_zero_weights_keep_fixed_channels_and_never_divide_by_zero():
    plan = make_plan([{}, {}], [-1, 0, 1])
    result = bounded_mac(plan, [{0: 2}, {1: -1}])
    assert plan.scale == 1
    assert tuple(result.coefficients) == (-2, -1, 0, 1, 2)
    assert result.peak_final == result.peak_partial == 0
    assert exact_pair(result.coefficients, result.scale) == (0, 0)


def test_bad_layout_bounds_and_float_oracle_inputs_are_rejected():
    with pytest.raises(ValueError):
        make_plan([{1: 1}], [0])
    for utilization in [0, -1, F(3, 2)]:
        with pytest.raises(ValueError):
            make_plan([{0: 1}], [0], utilization=utilization)
    plan = make_plan([{0: 1}], [0])
    for inputs in [[{0: 3}], [{1: 1}], []]:
        with pytest.raises(ValueError):
            bounded_mac(plan, inputs)
    with pytest.raises(TypeError):
        bounded_mac(plan, [{0: 0.5}])


def test_detector_thresholds_ties_and_offset():
    observations = {-2: F(-3, 2), -1: F(-1, 2), 0: F(1, 2), 1: F(3, 2)}
    assert detect_five(observations, [-2, -1, 0, 1]) == {-2: -1, -1: 0, 0: 1, 1: 2}
    assert detect_five({0: F(1, 2)}, [0], threshold_offset=F(1, 100)) == {0: 0}
    assert detect_five({0: 0}, [0], observation_errors={0: F(3, 5)}) == {0: 1}


def test_bounded_does_not_mean_discrete_or_value_preserving_after_detection():
    plan = make_plan([{0: 1}] * 5, [0])
    result = bounded_mac(plan, [{0: 1}, {}, {}, {}, {}])
    assert result.coefficients == {0: F(1, 5)}
    assert result.peak_final <= 2
    assert exact_pair(result.coefficients, result.scale) == (1, 0)
    detected = detect_five(result.coefficients, plan.output_positions)
    assert exact_pair(detected, result.scale) == (0, 0)


def test_same_scale_same_detected_state_can_hide_nonequivalent_values():
    plan = make_plan([{0: 1}] * 5, [0])
    zero = bounded_mac(plan, [{}] * 5)
    one = bounded_mac(plan, [{0: 1}, {}, {}, {}, {}])
    assert zero.scale == one.scale
    assert exact_pair(zero.coefficients) != exact_pair(one.coefficients)
    assert detect_five(zero.coefficients, plan.output_positions) == detect_five(one.coefficients, plan.output_positions)
    # This falsifies exactness of this candidate detector, not all gamma readout rules.


def test_exact_pair_recognizes_equivalence_not_just_identical_vectors():
    assert exact_pair({0: 3}) == exact_pair({-1: 1, 1: 1})
    assert exact_pair({0: F(3, 7)}) == exact_pair({-1: F(1, 7), 1: F(1, 7)})


def test_partial_mac_outputs_must_be_rescaled_before_merging():
    left_plan = make_plan([{0: 1}], [0])
    right_plan = make_plan([{0: 1}] * 2, [0])
    left = bounded_mac(left_plan, [{0: 2}])
    right = bounded_mac(right_plan, [{0: 2}, {0: 2}])
    assert left.coefficients[0] + right.coefficients[0] == 4  # Both bounded separately.
    merged_scale = left.scale + right.scale
    merged = {0: (left.scale * left.coefficients[0] + right.scale * right.coefficients[0]) / merged_scale}
    assert abs(merged[0]) <= 2
    assert exact_pair(merged, merged_scale) == (6, 0)


def test_two_stage_scale_propagation_and_detection_error_are_separate():
    first_plan = make_plan([{0: 1}] * 3, [0])
    first = bounded_mac(first_plan, [{0: 2}, {}, {}])
    second_plan = make_plan([{0: 2}], [0])
    continuous = bounded_mac(second_plan, [first.coefficients])
    total_scale = first.scale * continuous.scale
    assert exact_pair(continuous.coefficients, total_scale) == (4, 0)
    detected = detect_five(first.coefficients, first_plan.output_positions)
    regenerated = bounded_mac(second_plan, [detected])
    assert exact_pair(regenerated.coefficients, first.scale * regenerated.scale) == (6, 0)


def test_int16_boundary_encoding_and_normalized_product():
    values = [-32768, -12345, 0, 12345, 32767]
    words = {n: encode_int(n) for n in values}
    layout = sorted({k for word in words.values() for k in word})
    plan = make_plan([words[32767]], layout)
    for n in values:
        assert exact_value(words[n]) == (n, 0)
        result = bounded_mac(plan, [words[n]])
        assert result.peak_partial <= 2
        assert exact_pair(result.coefficients, result.scale) == (n*32767, 0)
