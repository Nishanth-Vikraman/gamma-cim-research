"""Exact local-normalization protocols and SPICE golden trace generator.

No floating beta is used for correctness. Quotient selection is discrete and
explicitly advertises the required threshold bands. Boundaries never wrap.
"""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
import csv
import json
from pathlib import Path
from gamma_native_model import convolve, normalize, value_pair


SCHEMES = {'unit': (1,), 'cap2': (1, 2), 'cap4': (1, 2, 3, 4),
           'pow2': (1, 2, 4)}


def bands(scheme):
    if scheme == 'bulk':
        return None
    if scheme not in SCHEMES:
        raise ValueError('unknown scheme')
    return SCHEMES[scheme]


def quotient(c, q, scheme='unit'):
    if q not in (3, 4, 5):
        raise ValueError('this phase supports q=3,4,5')
    m = abs(c)//q
    levels = bands(scheme)
    size = m if levels is None else max((v for v in levels if v <= m), default=0)
    return size if c >= 0 else -size


def pair_json(word, q):
    return [str(v) for v in value_pair(word, q)]


def clean(word):
    return {int(k): int(v) for k, v in word.items() if v}


def phase(state, q, scheme='unit', active=None):
    """One simultaneous phase. Outside active range, passive receivers retain
    transfers but do not fire. Returned firings do not observe incoming changes.
    """
    state = clean(state)
    firings = {k: quotient(c, q, scheme) for k, c in state.items()
               if (active is None or active[0] <= k <= active[1]) and abs(c) >= q}
    delta = defaultdict(int)
    for k, m in firings.items():
        delta[k] -= q*m
        delta[k-1] += m
        delta[k+1] += m
    result = defaultdict(int, state)
    for k, d in delta.items():
        result[k] += d
    result = clean(result)
    assert value_pair(result, q) == value_pair(state, q)
    return result, firings


def settle(state, q, scheme='unit', active=None, max_rounds=10000):
    state = clean(state)
    original = value_pair(state, q)
    trace = []
    packets = triggers = 0
    peak = max(map(abs, state.values()), default=0)
    while True:
        after, firings = phase(state, q, scheme, active)
        if not firings:
            break
        packets += sum(map(abs, firings.values()))
        triggers += len(firings)
        peak = max(peak, max(map(abs, after.values()), default=0))
        trace.append({'round': len(trace)+1, 'state': after, 'firings': firings,
                      'exact_pair': pair_json(after, q)})
        state = after
        if len(trace) >= max_rounds:
            raise RuntimeError('round budget exceeded')
    assert value_pair(state, q) == original
    boundary = {k: v for k, v in state.items() if active is not None and not active[0] <= k <= active[1]}
    return {'state': state, 'trace': trace, 'rounds': len(trace),
            'unit_rewrites': packets, 'trigger_events': triggers,
            'peak_units': peak, 'boundary_receivers': boundary,
            'all_persistent': all(abs(v) < q for v in state.values()),
            'boundary_requires_extension': bool(boundary), 'exact_pair': pair_json(state, q)}


def protocol_geometry(q, scheme):
    levels = bands(scheme)
    return {'persistent_radius': q-1, 'persistent_states': 2*q-1,
            'quotient_levels': levels, 'signed_command_states': 2*len(levels)+1 if levels else None,
            'overflow_comparators': 2*len(levels) if levels else None,
            'persistent_readout_comparators': 2*(q-1),
            'positive_thresholds_units': [q*m-0.5 for m in levels] if levels else None,
            'bulk_warning': 'unrestricted quotient is an arithmetic lower-bound control, not a cheap circuit primitive' if levels is None else None}


def default_cases(q):
    return [
        ('positive_steps', [(0, 1)]*q, 3),
        ('negative_steps', [(0, -1)]*q, 3),
        ('signed_sequence', [(0, 1)]*(q+1)+[(-1, -1), (1, 1), (0, -1)], 3),
        ('positive_burst', [(0, 4*q)], 3),
        ('negative_burst', [(0, -4*q)], 3),
        ('neighbor_simultaneous', [({-1: q, 0: q, 1: q})], 3),
        ('extended_cascade', [(0, 4*q*q)], 5),
        ('wide_cascade', [(0, 4*q*q)], 9),
        ('cancellation', [({-1: q, 0: -2*q, 1: q})], 5),
    ]


def golden_case(q, name, injections, nodes=3, scheme='unit', phase_budget=32):
    if nodes not in (3, 5, 9) or scheme == 'bulk':
        raise ValueError('3/5/9 active nodes and bounded scheme required')
    radius = nodes//2
    active = (-radius, radius)
    state = {}
    checkpoints, stimuli, phases = [], [], []
    t = 5.0
    for index, item in enumerate(injections):
        injection = dict(item) if isinstance(item, dict) else {item[0]: item[1]}
        if any(not -radius <= k <= radius or not isinstance(v, int) for k, v in injection.items()):
            raise ValueError('invalid signed integer injection')
        for k, v in injection.items():
            state[k] = state.get(k, 0)+v
            stimuli.append({'time_ns': t, 'node': k, 'delta_units': v})
        before = dict(state)
        result = settle(state, q, scheme, active)
        if result['rounds'] > phase_budget:
            raise ValueError('golden reference exceeds fixed circuit phase budget')
        # Fixed phase budget independent of the expected answer. Quiet phases
        # remain real comparator observations in the circuit/cost accounting.
        working = before
        for local_round in range(phase_budget):
            working, fire = phase(working, q, scheme, active)
            phases.append({'injection_index': index, 'phase': local_round,
                'sample_time_ns': t+2+4*local_round,
                'observe_time_ns': t+3.85+4*local_round,
                'state': working, 'firings': fire, 'exact_pair': pair_json(working, q)})
        checkpoint_time = t+4*phase_budget+1
        state = dict(result['state'])
        checkpoints.append({'time_ns': checkpoint_time, 'injection_index': index,
            'state': dict(state), 'exact_pair': result['exact_pair'],
            'rewrite_events': result['trigger_events'], 'unit_rewrites': result['unit_rewrites'],
            'synchronous_rounds': result['rounds'], 'boundary_receivers': result['boundary_receivers'],
            'boundary_requires_extension': result['boundary_requires_extension'],
            'all_persistent': result['all_persistent']})
        t = checkpoint_time+5
    return {'schema_version': 1, 'q': q, 'scheme': scheme, 'name': name,
            'active_nodes': list(range(-radius, radius+1)),
            'passive_boundary_receivers': [-radius-1, radius+1],
            'phase_budget_per_injection': phase_budget, 'stimuli': stimuli,
            'phases': phases, 'checkpoints': checkpoints, 'stop_time_ns': t,
            'geometry': protocol_geometry(q, scheme),
            'boundary_policy': 'passive receivers preserve value; any outward transfer flags active-window extension, never a sealed three-node success'}


def export_golden(case, destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    stem = f"q{case['q']}_{case['scheme']}_{case['name']}"
    json_path = destination/f'{stem}.json'
    json_path.write_text(json.dumps(case, indent=2)+'\n', encoding='utf-8')
    with (destination/f'{stem}.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['time_ns', 'node', 'delta_units'])
        writer.writeheader()
        writer.writerows(case['stimuli'])
    return json_path


def workload_study(seed=0, count=128):
    import random
    rng = random.Random(seed)
    pairs = [(rng.randint(-128, 127), rng.randint(-128, 127)) for _ in range(count)]
    rows = []
    for q in (3, 4, 5):
        enc = {n: normalize({0: n}, q)[0] for n in range(-128, 128)}
        for scheme in (*SCHEMES, 'bulk'):
            for a, b in pairs:
                raw = convolve(enc[a], enc[b])
                result = settle(raw, q, scheme)
                assert value_pair(result['state'], q) == (a*b, 0)
                lo = min([*raw, *result['state']], default=0)
                hi = max([*raw, *result['state']], default=0)
                width = hi-lo+1
                levels = bands(scheme)
                rows.append({'q': q, 'scheme': scheme, 'inputs': [a, b],
                    'digit_products': len(enc[a])*len(enc[b]),
                    'input_spans': [max(d)-min(d)+1 if d else 0 for d in (enc[a], enc[b])],
                    'width': width, 'rounds': result['rounds'],
                    'trigger_events': result['trigger_events'], 'unit_rewrites': result['unit_rewrites'],
                    'peak_units': result['peak_units'],
                    'bank_observations_until_settled_plus_confirmation': width*(result['rounds']+1),
                    'comparator_evaluations_until_settled_plus_confirmation': width*(result['rounds']+1)*2*len(levels) if levels else None,
                    'final_readout_comparators': width*2*(q-1),
                    'boundary_transfers_per_trigger': 2,
                    'exact': True})
    return rows


def transient_study(seed=0):
    """Delay streaming normalization until a transient trigger, then flush to
    persistent D at the final checkpoint. Wider temporary ranges are not free
    persistent states, and final flush work is never omitted.
    """
    import random
    rng = random.Random(seed+9)
    sequences = [[rng.choice((-1, 1)) for _ in range(64)] for _ in range(32)]
    sequences += [[1]*64, [-1]*64]
    rows = []
    for q in (3, 4, 5):
        for limit in (q-1, 2*q-1, 4*q-1, 8*q-1):
            for sequence in sequences:
                state = {}
                events = packets = rounds = peak = 0
                for value in sequence:
                    state[0] = state.get(0, 0)+value
                    peak = max(peak, max(map(abs, state.values()), default=0))
                    if max(map(abs, state.values()), default=0) > limit:
                        st = settle(state, q, 'unit')
                        state = st['state']
                        events += st['trigger_events']; packets += st['unit_rewrites']; rounds += st['rounds']
                        peak = max(peak, st['peak_units'])
                st = settle(state, q, 'unit')
                events += st['trigger_events']; packets += st['unit_rewrites']; rounds += st['rounds']
                assert value_pair(st['state'], q) == (sum(sequence), 0)
                rows.append({'q': q, 'transient_limit_units': limit, 'temporary_integer_classes': 2*limit+1,
                    'input_events': len(sequence), 'trigger_events': events, 'unit_rewrites': packets,
                    'serial_sum_of_settle_rounds': rounds, 'peak_units': max(peak, st['peak_units']),
                    'final_persistent': True})
    return rows
