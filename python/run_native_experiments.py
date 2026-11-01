"""Deterministic architecture-v2 studies; no physical energy claims."""
from __future__ import annotations
import argparse
from bisect import bisect_left
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
import json
import math
from pathlib import Path
import random
import statistics

from gamma_bounded_model import bounded_mac, make_plan
from gamma_exact_model import encode_int
from gamma_native_model import (NativeFloat, add, beta_decimal, block_dot, block_encode,
    convolve, multiply, normalize, pair_multiply, route_align,
    sense, sense_threshold_vector, value_pair)
from run_fixed_bins_experiments import frozen_workloads


def scale_study(seed):
    scales = [F(1, 7), F(1), F(5), F(29, 3), F(100), F(1000)]
    rng = random.Random(seed)
    scales += [F(rng.randint(1, 100000), rng.randint(1, 1000)) for _ in range(100)]
    successes = sum(sense(F(d)/s, s).digit == d for s in scales for d in range(-2, 3))
    layout = sorted({k for n in range(-128, 128) for k in encode_int(n)})
    records = []
    for name, xs, ws in frozen_workloads(seed):
        plan = make_plan([encode_int(w) for w in ws], layout)
        result = bounded_mac(plan, [encode_int(x) for x in xs])
        # Compile radius from weights/layout, never from the observed answer.
        radius = max(2, math.ceil(2*max(plan.channel_l1.values(), default=0)))
        narrow = {k: sense(v, plan.scale).digit for k, v in result.coefficients.items()}
        wide = {k: sense(v, plan.scale, radius).digit for k, v in result.coefficients.items()}
        expected = (sum(x*w for x, w in zip(xs, ws)), 0)
        # Readout returns c, not c/S: do not multiply these digits by S again.
        records.append({'workload': name, 'scale': str(plan.scale),
            'expected_dot': expected[0], 'fixed_five_exact': value_pair(narrow) == expected,
            'raw_coefficients_outside_alphabet': sum(abs(v*plan.scale)>2 for v in result.coefficients.values()),
            'compiled_transient_radius': radius, 'transient_classes': 2*radius+1,
            'equivalent_flash_thresholds_per_channel': 2*radius,
            'wide_transient_exact': value_pair(wide) == expected,
            'wide_coefficients_exact': all(F(wide[k]) == v*plan.scale for k, v in result.coefficients.items()),
            'minimum_center_spacing_physical': float(1/plan.scale)})
    trials = failures = events = 0
    for c in range(-128, 129):
        digits, stats = normalize({0: c})
        events += stats['unit_events']
        for scale in scales[:6]:
            detected = {k: sense(F(d)/scale, scale).digit for k, d in digits.items()}
            trials += 1
            failures += value_pair(detected) != (c, 0)
    return {'legal_digit_trials': 5*len(scales), 'legal_digit_exact': successes,
        'fixed_five_dot_exact': sum(r['fixed_five_exact'] for r in records),
        'wide_transient_dot_exact': sum(r['wide_transient_exact'] for r in records),
        'workloads': len(records), 'records': records,
        'ideal_normalize_then_sense': {'trials': trials, 'failures': failures,
            'reference_unit_events_once_per_coefficient': events,
            'warning': 'Software normalization precedes sensing; physical loop unimplemented.'}}


def noise_study(seed, trials=2000):
    records = []
    for profile in ['ideal', 'additive', 'comparator_offset', 'gain', 'conductance', 'correlated_threshold', 'combined']:
        for scale in [1, 5, 20, 100, 500]:
            rng = random.Random(seed+1001)
            failures = inverted = 0
            sigma = 0.002 if profile in ('additive', 'combined') else 0
            for trial in range(trials):
                d = trial % 5-2
                zn, zo, zg, zc = [rng.gauss(0, 1) for _ in range(4)]
                zcell = [rng.gauss(0, 1) for _ in range(5)]
                zth = [rng.gauss(0, 1) for _ in range(4)]
                current = d/scale
                if profile in ('gain', 'combined'):
                    current *= 1+0.01*zg
                if profile in ('conductance', 'combined'):
                    current += sum(x*0.01*z for x, z in zip([2, 2, d, -2, -2], zcell))/scale
                current += sigma*zn
                offset = 0.002*zo if profile in ('comparator_offset', 'combined') else 0
                if profile == 'correlated_threshold':
                    # Common bank offset plus small independent threshold mismatch.
                    offsets = [F(0.002*zc+0.0002*z) for z in zth]
                    try:
                        actual = sense_threshold_vector(F(current), scale, 2, offsets)
                    except ValueError:
                        inverted += 1
                        failures += 1
                        continue
                else:
                    actual = sense(F(current), scale, offset=F(offset)).digit
                failures += actual != d
            records.append({'profile': profile, 'scale': scale, 'trials': trials,
                'wrong_class_probability': failures/trials, 'wrong_class_count': failures,
                'threshold_order_inversions': inverted, 'spacing_physical': 1/scale,
                'additive_sigma_physical': sigma,
                'spacing_over_additive_sigma': 1/(scale*sigma) if sigma else None,
                'absolute_offset_sigma_physical': 0.002 if profile in ('comparator_offset', 'combined', 'correlated_threshold') else 0,
                'relative_gain_or_cell_sigma': 0.01 if profile in ('gain', 'conductance', 'combined') else 0})
    return {'records': records,
        'uncertainty': 'Monte Carlo point estimates, paired draws across scales; zero observed errors does not prove zero probability.',
        'physical_units': 'Normalized current, not a calibrated device.',
        'conductance_setup': 'Five unit weights with signed inputs [2,2,d,-2,-2]; independent relative mismatch per cell.'}


def family_study(seed):
    rng = random.Random(seed+21)
    pairs = [(rng.randint(-128, 127), rng.randint(-128, 127)) for _ in range(128)]
    rows = []
    for q in [3, 4, 5, 6]:
        enc = {n: normalize({0: n}, q)[0] for n in range(-128, 128)}
        spans = [max(d)-min(d)+1 if d else 0 for d in enc.values()]
        sparsity, events, peaks, rounds = [], [], [], []
        failures = 0
        for a, b in pairs:
            out, st = normalize(convolve(enc[a], enc[b]), q)
            sparsity.append(len(enc[a])*len(enc[b]))
            events.append(st['unit_events'])
            peaks.append(st['peak_coefficient'])
            rounds.append(st['bulk_rounds'])
            failures += value_pair(out, q) != (a*b, 0)
        rows.append({'q': q, 'beta': float(beta_decimal(q)), 'persistent_radius': q-1,
            'persistent_states': 2*q-1, 'int8_median_span': statistics.median(spans),
            'int8_max_span': max(spans), 'int8_median_nonzero': statistics.median(map(len, enc.values())),
            'observed_span_per_input_bit': statistics.median(spans)/8,
            'product_median_nonzero_digit_products': statistics.median(sparsity),
            'product_median_unit_events': statistics.median(events),
            'product_median_bulk_rounds': statistics.median(rounds),
            'product_max_transient_coefficient': max(peaks),
            'exact_product_failures': failures, 'product_cases': len(pairs)})
    return {'records': rows, 'caveat': 'Exact integer encoder uses radius q-1: fewer positions trade against more local states. Fixed radius-2 finite precision studied separately. Bulk rounds are not circuit latency.'}


def float_and_routing_study(seed):
    rng = random.Random(seed+40)
    failures = cases = overflow = 0
    for q in [3, 4, 5, 6]:
        for _ in range(100):
            da, db = [normalize({0: rng.randint(-127, 127)}, q)[0] for _ in range(2)]
            a = NativeFloat(da, rng.randint(-100, 100), rng.choice([-1, 1]), q)
            b = NativeFloat(db, rng.randint(-100, 100), rng.choice([-1, 1]), q)
            c, _ = multiply(a, b)
            summed, _ = add(a, b)
            failures += c.exact() != pair_multiply(a.exact(), b.exact(), q)
            failures += summed.exact() != tuple(x+y for x, y in zip(a.exact(), b.exact()))
            failures += a.renormalize(7).exact() != a.exact()
            overflow += not -128 <= c.exponent <= 127
            cases += 3
    routes = []
    for gap in [0, 2, 8, 32]:
        a, b = NativeFloat({0: 2, -1: 1}), NativeFloat({0: 2, -1: -1}, -gap)
        out, stats = add(a, b)
        _, _, edges, dropped = route_align([a.digits, b.digits], [a.exponent, b.exponent], window=(-7, 0))
        routes.append({'exponent_gap': gap, **stats, 'finite_window': [-7, 0],
            'finite_active_connections': len(edges), 'finite_dropped_terms': len(dropped),
            'finite_dropped_exact_pair': [str(v) for v in value_pair(dropped)],
            'exact_unbounded_add': out.exact() == tuple(x+y for x, y in zip(a.exact(), b.exact()))})
    return {'exact_identity_checks': cases, 'failures': failures,
        'signed_8bit_product_exponent_overflows': overflow, 'multiply_cases': 400,
        'routing': routes, 'limits': 'Unbounded Python integers and indices; finite-window loss separate. Exponent overflow needs handling. Nonlinear operations unimplemented.'}


def block_study(seed):
    rng = random.Random(seed+60)
    values = [Decimal(str(rng.uniform(-1, 1)))*Decimal(2)**rng.randint(-12, 12) for _ in range(64)]
    weights = [Decimal(str(rng.uniform(-1, 1))) for _ in values]
    records = []
    for p in [4, 8, 12]:
        for name, group in [('tensor', 64), ('channel', 8), ('tile', 16), ('small_block', 4), ('scalar_reference', 1)]:
            decoded, exponents, decoded_w, weight_exponents = [], [], [], []
            native_blocks = []
            dot_failures = dot_events = block_digit_products = 0
            for start in range(0, len(values), group):
                xwords, e, vals = block_encode(values[start:start+group], positions=p)
                wwords, ew, wvals = block_encode(weights[start:start+group], positions=p)
                native, st = block_dot(xwords, wwords, e, ew)
                native_blocks.append(native)
                expected_pair = [F(0), F(0)]
                for xword, wword in zip(xwords, wwords):
                    term = pair_multiply(NativeFloat(xword, e).exact(), NativeFloat(wword, ew).exact())
                    expected_pair = [a+b for a, b in zip(expected_pair, term)]
                dot_failures += native.exact() != tuple(expected_pair)
                dot_events += st['unit_events']
                block_digit_products += st['digit_products']
                decoded.extend(vals)
                decoded_w.extend(wvals)
                exponents.append(e)
                weight_exponents.append(ew)
            error = [float(a-b) for a, b in zip(decoded, values)]
            relative = [abs(float((a-b)/b)) for a, b in zip(decoded, values) if b]
            reference_dot = sum(x*w for x, w in zip(values, weights))
            converted_dot = sum(x*w for x, w in zip(decoded, decoded_w))
            accumulator = NativeFloat({})
            routing_connections = merge_events = 0
            expected_pair = [F(0), F(0)]
            for block in native_blocks:
                expected_pair = [a+b for a, b in zip(expected_pair, block.exact())]
                accumulator, st = add(accumulator, block)
                routing_connections += st['active_connections']
                merge_events += st['unit_events']
            dot_failures += accumulator.exact() != tuple(expected_pair)
            records.append({'grouping': name, 'block_size': group, 'positions': p,
                'values': len(values), 'exponent_metadata_count_per_operand': len(exponents),
                'activation_exponents': exponents, 'weight_exponents': weight_exponents,
                'relative_l2_conversion_error': math.sqrt(sum(e*e for e in error)/sum(float(x)**2 for x in values)),
                'max_relative_conversion_error_nonzero': max(relative),
                'nonzero_values_rounded_to_zero': sum(a == 0 and b != 0 for a, b in zip(decoded, values)),
                'relative_dot_error_from_boundary_conversion': abs(float((converted_dot-reference_dot)/reference_dot)),
                'upper_bound_digit_slots_per_operand': len(values)*p,
                'exact_native_block_dot_or_merge_failures': dot_failures,
                'reference_mac_normalization_events': dot_events,
                'reference_merge_normalization_events': merge_events,
                'block_merge_active_unity_connections': routing_connections,
                'nonzero_digit_products': block_digit_products})
    return {'records': records, 'shape': [8, 8], 'scope': 'Boundary encoding plus exact native block mantissa MACs and differently exponented block accumulation. Normalization is ideal software, unbounded indices, no analog noise or physical cost. Synthetic data, not a neural accuracy benchmark.'}


def format_values(seed):
    rng = random.Random(seed+70)
    fp8 = []
    for sign in [-1, 1]:
        for bits in range(127):  # E4M3FN: magnitude 0x7f NaN; 0x7e is 448.
            e, m = bits >> 3, bits & 7
            v = math.ldexp(m/8, -6) if e == 0 else math.ldexp(1+m/8, e-7)
            fp8.append(Decimal(sign)*Decimal(v))
    outputs = {'INT8': [Decimal(n) for n in range(-128, 128)], 'FP8_E4M3FN': fp8}
    for name, ebits, mbits, bias in [('BF16', 8, 7, 127), ('FP16', 5, 10, 15)]:
        max_finite = ((1 << ebits)-1)*(1 << mbits)-1
        codes = [0, 1, (1 << mbits)-1, 1 << mbits, bias << mbits, max_finite]
        codes += [rng.randint(0, max_finite) for _ in range(122)]
        values = []
        for code in codes:
            e, m = code >> mbits, code & ((1 << mbits)-1)
            v = math.ldexp(m/(1 << mbits), 1-bias) if e == 0 else math.ldexp(1+m/(1 << mbits), e-bias)
            value = Decimal(v)
            values.extend([value, value.copy_negate()])
        outputs[name] = values
    return outputs


def precision_study(seed):
    codebooks = []
    for q in [3, 4, 5, 6]:
        for p in [3, 5]:
            pairs = {value_pair(dict(zip(range(1-p, 1), digits)), q)
                     for digits in product(range(-2, 3), repeat=p)}
            with localcontext() as ctx:
                ctx.prec = 90
                beta = beta_decimal(q)
                vals = sorted(Decimal(a.numerator)/Decimal(a.denominator)
                              +Decimal(b.numerator)/Decimal(b.denominator)*beta for a, b in pairs)
                gaps = [b-a for a, b in zip(vals, vals[1:])]
                inside = [v for v in vals if -1 <= v <= 1]
                errors, relative = [], []
                for i in range(1025):
                    target = Decimal(-1)+Decimal(i)/512
                    ix = bisect_left(vals, target)
                    nearest = min(vals[max(0, ix-1):min(len(vals), ix+1)], key=lambda v: abs(v-target))
                    errors.append(abs(nearest-target))
                    if target:
                        relative.append(abs((nearest-target)/target))
                worst = float(max(errors))
                codebooks.append({'q': q, 'positions': p, 'radius': 2, 'local_states': 5,
                    'raw_patterns': 5**p, 'distinct_exact_values': len(pairs),
                    'minimum_distinct_step': float(min(gaps)),
                    'max_adjacent_gap_in_unit_interval': float(max((b-a for a, b in zip(inside, inside[1:])), default=0)),
                    'max_abs_error_on_unit_grid': worst,
                    'max_relative_error_on_nonzero_unit_grid': float(max(relative)),
                    'equivalent_absolute_bits_unit_grid': -math.log2(worst) if worst else None,
                    'grid_points': 1025})
    conversions = []
    formats = format_values(seed)
    for q in [3, 4, 5, 6]:
        for p in [4, 8, 12]:
            for name, values in formats.items():
                encodings = [block_encode([v], q=q, positions=p, radius=2) for v in values]
                decoded = [item[2][0] for item in encodings]
                # Exactness comes from rational algebraic pairs, never Decimal
                # cancellation artifacts. In particular exact INT8 gets zero error.
                exact = [NativeFloat(item[0][0], item[1], q=q).exact() == (F(v), 0)
                         for item, v in zip(encodings, values)]
                errors = [0.0 if ok else abs(float(a-b)) for a, b, ok in zip(decoded, values, exact)]
                relative = [0.0 if ok else abs(float((a-b)/b)) for a, b, ok in zip(decoded, values, exact) if b]
                worst_rel = max(relative)
                conversions.append({'format': name, 'q': q, 'positions': p, 'radius': 2,
                    'samples': len(values), 'coverage': 'exhaustive finite' if name in ('INT8', 'FP8_E4M3FN') else 'seeded bit-pattern sample plus extremes',
                    'max_abs_conversion_error': max(errors), 'max_relative_conversion_error': worst_rel,
                    'equivalent_relative_mantissa_bits': -math.log2(worst_rel) if worst_rel else None,
                    'exact_conversion_count': sum(exact),
                    'nonzero_to_zero_count': sum(a == 0 and b != 0 for a, b in zip(decoded, values))})
    return {'codebooks': codebooks, 'format_conversions': conversions,
        'definitions': 'Absolute bits=-log2(max absolute error) at unit amplitude. Relative bits=-log2(max relative error) on specified samples; neither is physical ENOB. Codebooks exhaustive; unit-grid error sampled; greedy conversion not optimal. No positions*log2(beta) precision claim.'}


def run(seed=0, noise_trials=2000):
    return {'schema_version': 2, 'seed': seed,
        'scope': 'Scale-aware regions; base family; floating/block floating point; unity routing; empirical precision.',
        'scale_aware': scale_study(seed), 'noise': noise_study(seed, noise_trials),
        'base_family': family_study(seed), 'floating_and_routing': float_and_routing_study(seed),
        'block_floating': block_study(seed), 'precision': precision_study(seed),
        'limits': ['No SPICE/device measurements or energy/area/GPU advantage.',
                   'Exact normalization is offline, not a demonstrated analog loop.',
                   'Physical scale differs from sign and radix exponent.',
                   'General fractional learned coefficients need not lie on the integer digit lattice.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--noise-trials', type=int, default=2000)
    parser.add_argument('--output', type=Path, default=Path(__file__).parent/'results/native_baseline.json')
    args = parser.parse_args()
    if args.noise_trials < 1:
        parser.error('--noise-trials must be positive')
    report = run(args.seed, args.noise_trials)
    scale = report['scale_aware']
    print(f"Legal digit scale recovery: {scale['legal_digit_exact']}/{scale['legal_digit_trials']}")
    print(f"Raw dots: fixed five {scale['fixed_five_dot_exact']}/{scale['workloads']}; wide transient {scale['wide_transient_dot_exact']}/{scale['workloads']}")
    print('Ideal normalize-then-sense:', scale['ideal_normalize_then_sense'])
    print('Exact floating/routing:', report['floating_and_routing']['exact_identity_checks'], 'failures:', report['floating_and_routing']['failures'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n', encoding='utf-8')
    print('Saved', args.output)
    if scale['legal_digit_exact'] != scale['legal_digit_trials'] or scale['wide_transient_dot_exact'] != scale['workloads'] or report['floating_and_routing']['failures'] or scale['ideal_normalize_then_sense']['failures'] or any(r['exact_product_failures'] for r in report['base_family']['records']) or any(r['exact_native_block_dot_or_merge_failures'] for r in report['block_floating']['records']):
        raise SystemExit('Exact supported-domain gate failed')


if __name__ == '__main__':
    main()
