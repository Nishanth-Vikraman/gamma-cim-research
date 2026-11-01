"""Historical fixed-unscaled-bin control, superseded by scale-aware sensing.

Use --output PATH for deterministic JSON. No post-MAC topplings are used.
This intentionally retains the old unscaled thresholds and scale restoration.
Noise/offset figures are synthetic, not device measurements. Energy and physical
timing are deliberately not inferred from exact Python operations.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path
import random

from gamma_bounded_model import THRESHOLDS, bounded_mac, detect_five, exact_pair, make_plan
from gamma_exact_model import encode_int


@lru_cache(maxsize=None)
def encode(n):
    return encode_int(n)


def numeric(pair):
    # Reporting only; all equality decisions are made on exact rational pairs.
    with localcontext() as context:
        context.prec = 80
        gamma = (Decimal(3) + Decimal(5).sqrt()) / 2
        a, b = pair
        return float(Decimal(a.numerator) / Decimal(a.denominator)
                     + Decimal(b.numerator) / Decimal(b.denominator) * gamma)


def frozen_workloads(seed):
    rng = random.Random(seed)
    workloads = []
    for length in [1, 16, 64, 256]:
        # Keep weights fixed across multiple activation samples and all settings.
        ws = [rng.randint(-128, 127) for _ in range(length)]
        for sample in range(4):
            xs = [rng.randint(-128, 127) for _ in range(length)]
            workloads.append((f'uniform_{length}_{sample}', xs, ws))
    workloads.extend([
        ('zero_weights', [127] * 64, [0] * 64),
        ('signed_cancellation', [127] * 64, [127, -127] * 32),
        ('positive_aligned', [127] * 256, [127] * 256),
        ('single_outlier', [127] + [1] * 63, [127] + [0] * 63),
        ('small_inputs', [rng.randint(-2, 2) for _ in range(64)],
         [rng.randint(-128, 127) for _ in range(64)]),
        ('identity', [1], [1]),
    ])
    return workloads


def collision_audit(layout):
    # Layout, fixed weights and scale are identical for all 256 input values.
    plan = make_plan([encode(127)], layout)
    groups = {}
    for x in range(-128, 128):
        result = bounded_mac(plan, [encode(x)])
        digits = detect_five(result.coefficients, plan.output_positions)
        code = tuple(digits[k] for k in plan.output_positions)
        groups.setdefault(code, []).append(x)
    collisions = [xs for xs in groups.values() if len(xs) > 1]
    return {
        'detector': 'candidate_nearest_five_only',
        'fixed_weight': 127, 'fixed_scale': str(plan.scale),
        'inputs': 256, 'distinct_detected_states': len(groups),
        'nonequivalent_collision_groups': len(collisions),
        'max_group_size': max(map(len, groups.values())),
        'example_input_collision': collisions[0][:2] if collisions else None,
        'interpretation': 'Tests this detector with this layout/scale, not all gamma equivalence-region rules.',
    }


def run(seed=0):
    layout = sorted({k for n in range(-128, 128) for k in encode(n)})
    workloads = frozen_workloads(seed)
    settings = [(F(0), F(0)), (F(1, 50), F(0)), (F(1, 10), F(1, 50))]
    records = []
    for utilization in [F(1), F(3, 4), F(1, 2)]:
        for index, (name, xs, ws) in enumerate(workloads):
            plan = make_plan([encode(w) for w in ws], layout, utilization=utilization)
            result = bounded_mac(plan, [encode(x) for x in xs])
            expected = sum(x*w for x, w in zip(xs, ws))
            exact_ok = exact_pair(result.coefficients, result.scale) == (expected, 0)
            ideal = detect_five(result.coefficients, plan.output_positions)
            for setting, (sigma, offset) in enumerate(settings):
                # Identical standardized noise draws across utilization settings.
                noise_rng = random.Random(seed + 10000 + index * 31 + setting)
                errors = {k: F.from_float(noise_rng.gauss(0, float(sigma))) for k in plan.output_positions}
                detected = detect_five(result.coefficients, plan.output_positions,
                                       threshold_offset=offset, observation_errors=errors)
                decoded = exact_pair(detected, result.scale)
                error_pair = (decoded[0] - expected, decoded[1])
                nearest_margin = min(abs(v - (threshold + offset))
                                     for v in result.coefficients.values() for threshold in THRESHOLDS)
                records.append({
                    'workload': name, 'dot_length': len(xs),
                    'utilization': str(utilization), 'scale': str(result.scale),
                    'noise_sigma_digit_units': str(sigma), 'threshold_offset_digit_units': str(offset),
                    'exact_before_sensing': exact_ok,
                    'ideal_range_violation': result.peak_final > 2 or result.peak_partial > 2,
                    'peak_final': float(result.peak_final), 'peak_partial': float(result.peak_partial),
                    'nonzero_digit_products': result.digit_products,
                    'physical_output_channels': len(plan.output_positions),
                    'candidate_thresholds_per_channel': 4,
                    'expected_dot': expected, 'exact_after_sensing': decoded == (expected, 0),
                    'error_pair': [str(error_pair[0]), str(error_pair[1])],
                    'restored_error': numeric(error_pair),
                    'minimum_threshold_margin': float(nearest_margin),
                    'state_flips_vs_ideal': sum(detected[k] != ideal[k] for k in plan.output_positions),
                    'observations_outside_nominal_range': sum(abs(result.coefficients[k] + errors[k]) > 2 for k in plan.output_positions),
                })
    summaries = []
    for utilization in ['1', '3/4', '1/2']:
        for sigma, offset in settings:
            group = [r for r in records if r['utilization'] == utilization
                     and r['noise_sigma_digit_units'] == str(sigma)
                     and r['threshold_offset_digit_units'] == str(offset)]
            error_sq = sum(r['restored_error']**2 for r in group)
            reference_sq = sum(r['expected_dot']**2 for r in group)
            summaries.append({
                'utilization': utilization, 'noise_sigma': str(sigma), 'threshold_offset': str(offset),
                'cases': len(group),
                'range_violations': sum(r['ideal_range_violation'] for r in group),
                'pre_sense_exact_failures': sum(not r['exact_before_sensing'] for r in group),
                'post_sense_exact_cases': sum(r['exact_after_sensing'] for r in group),
                'relative_l2_error': math.sqrt(error_sq / reference_sq) if reference_sq else None,
                'mean_absolute_error_original_units': sum(abs(r['restored_error']) for r in group) / len(group),
                'state_flip_fraction': sum(r['state_flips_vs_ideal'] for r in group)
                    / sum(r['physical_output_channels'] for r in group),
            })
    return {
        'schema_version': 1, 'seed': seed,
        'model': 'bounded rational MAC with fixed weight-only common scale; no MAC topplings',
        'detector': 'candidate nearest-five fixed thresholds, not a validated gamma equivalence rule',
        'input_layout': layout, 'workload_count': len(workloads),
        'noise_note': 'Synthetic digit-unit noise and common threshold offset; no circuit claim.',
        'energy_note': 'No measured energy, latency, or matched binary hardware comparison yet.',
        'summaries': summaries, 'collision_audit': collision_audit(layout), 'records': records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = run(args.seed)
    print('BOUNDED MAC: fixed workloads/layouts; tracked weight-only scales')
    print('Candidate nearest-five detector only; post-sensing exactness is not assumed.')
    for summary in report['summaries']:
        print(f"util={summary['utilization']:>3} sigma={summary['noise_sigma']:>4} "
              f"offset={summary['threshold_offset']:>4}: "
              f"range_fail={summary['range_violations']} "
              f"pre_sense_fail={summary['pre_sense_exact_failures']} "
              f"post_sense_exact={summary['post_sense_exact_cases']}/{summary['cases']} "
              f"relative_L2={summary['relative_l2_error']:.4f}")
    print('Fixed-scale INT8 collision audit:', json.dumps(report['collision_audit'], sort_keys=True))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        print('Saved deterministic report:', args.output)
    if any(s['range_violations'] or s['pre_sense_exact_failures'] for s in report['summaries']):
        raise SystemExit('Bounded-MAC acceptance gate failed')


if __name__ == '__main__':
    main()
