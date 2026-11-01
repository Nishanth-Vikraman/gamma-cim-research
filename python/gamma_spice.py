"""LTspice behavioral charge-packet normalizer and declared-region verifier.

Capacitor states drive real behavioral comparator/sample/hold sources. Python
never chooses circuit firings and never normalizes measured voltages afterward.
"""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from bisect import bisect_left
import json
import math
from pathlib import Path
import random
import subprocess
import time
import numpy as np
from gamma_protocols import bands, pair_json

DEFAULT_EXE = Path(r'C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe')


@dataclass
class Physical:
    scale: float = 4.0
    capacitance_f: float = 100e-15
    unit_voltage_v: float = 0.1
    packet_width_ns: float = 1.0
    max_step_ns: float = 0.08
    threshold_offset_v: float = 0.0
    comparator_offset_sigma_v: float = 0.0
    capacitor_sigma_rel: float = 0.0
    transfer_sigma_rel: float = 0.0
    pulse_width_sigma_rel: float = 0.0
    leakage_tau_ns: float = 1e15
    injection_gain_sigma_rel: float = 0.0
    sense_noise_sigma_v: float = 0.0
    current_noise_sigma_a: float = 0.0
    correlated_threshold_sigma_v: float = 0.0
    headroom_v: float = 0.6


def node_name(k):
    return f'nm{-k}' if k < 0 else f'np{k}'


def pwl(points):
    merged = {}
    for t, v in points:
        merged[round(t, 9)] = v
    return 'PWL('+ ' '.join(f'{t:.9g}n {v:.14g}' for t,v in sorted(merged.items()))+')'


def pulses(starts, width, amplitude=1.0, rise=0.01):
    pts=[(0,0)]
    for start in starts:
        pts.extend([(start,0), (start+rise,amplitude),
                    (start+width-rise,amplitude), (start+width,0)])
    return pts


def realization(case, physical: Physical, seed):
    rng = random.Random(seed)
    active = case['active_nodes']
    all_nodes = sorted(active+case['passive_boundary_receivers'])
    common = rng.gauss(0, physical.correlated_threshold_sigma_v)
    return {'cap_gains': {k: 1+rng.gauss(0, physical.capacitor_sigma_rel) for k in all_nodes},
        'offsets_v': {k: physical.threshold_offset_v+common+rng.gauss(0,physical.comparator_offset_sigma_v) for k in active},
        'injection_gains': {k: 1+rng.gauss(0,physical.injection_gain_sigma_rel) for k in all_nodes},
        'transfer_gains': {k: [1+rng.gauss(0,physical.transfer_sigma_rel) for _ in range(3)] for k in active},
        'width_gains': {k: 1+rng.gauss(0,physical.pulse_width_sigma_rel) for k in active},
        'sense_noise_v': {k: [rng.gauss(0,physical.sense_noise_sigma_v) for _ in case['phases']] for k in active},
        'current_noise_a': {k: [rng.gauss(0,physical.current_noise_sigma_a) for _ in case['phases']] for k in all_nodes}}


def netlist(case, physical: Physical, sample):
    q, scale = case['q'], physical.scale
    if scale <= 0 or physical.capacitance_f <= 0 or physical.unit_voltage_v <= 0:
        raise ValueError('positive physical scale/capacitance/voltage required')
    levels = bands(case['scheme'])
    if levels is None:
        raise ValueError('unrestricted bulk quotient is not a circuit')
    active=case['active_nodes']
    all_nodes=sorted(active+case['passive_boundary_receivers'])
    if any(g<=0 for g in sample['cap_gains'].values()) or any(g<=0 for g in sample['width_gains'].values()):
        raise ValueError('invalid negative physical realization')
    unit=physical.unit_voltage_v/scale
    charge=physical.capacitance_f*unit
    lines=[f"Gamma-CIM q={q} {case['scheme']} {case['name']}",
        '* Signed zero-centered states stand for differential storage.',
        '* Ideal behavioral sources do not provide a transistor supply-energy estimate.',
        '.model SAMP SW(Ron=1 Roff=1e18 Vt=0.5 Vh=0.05)',
        '.options plotwinsize=0 numdgt=15 reltol=1e-6 abstol=1e-15 vntol=1e-9',
        f"Vclk clk 0 {pwl(pulses([p['sample_time_ns'] for p in case['phases']],0.2))}"]
    for k in all_nodes:
        n=node_name(k)
        cap=physical.capacitance_f*sample['cap_gains'][k]
        lines += [f'C_{n} {n} 0 {cap:.14g} IC=0',
                  f'Rleak_{n} {n} 0 {physical.leakage_tau_ns*1e-9/cap:.14g}']
        points=[(0,0)]
        for event in case['stimuli']:
            if event['node']==k:
                amp=charge*event['delta_units']*sample['injection_gains'][k]/(0.5e-9-0.01e-9)
                points += pulses([event['time_ns']],0.5,amp)[1:]
        lines.append(f'Iinput_{n} 0 {n} {pwl(points)}')
        # Band-limited current disturbance, held between phase sample times.
        # Independent noise realization is external stimulus, never a firing oracle.
        npts=[(0,sample['current_noise_a'][k][0])]
        for p,z in zip(case['phases'], sample['current_noise_a'][k]):
            npts.extend([(p['sample_time_ns']-0.5,z), (p['observe_time_ns'],z)])
        lines.append(f'Inoise_{n} 0 {n} {pwl(npts)}')
    for k in active:
        n=node_name(k); hold='h_'+n
        # The sum of a few comparator steps selects a discrete packet size.
        positive='0'; negative='0'; previous=0
        signal=f'(V({n})+V(z_{n})-({sample["offsets_v"][k]:.14g}))'
        for m in levels:
            threshold=(q*m-0.5)*unit
            positive += f'+if({signal}>{threshold:.14g},{m-previous},0)'
            negative += f'+if({signal}<-{threshold:.14g},{m-previous},0)'
            previous=m
        zpoints=[(0,0)]
        for p,z in zip(case['phases'],sample['sense_noise_v'][k]):
            zpoints.extend([(p['sample_time_ns']-0.5,z), (p['observe_time_ns'],z)])
        lines += [f'Vz_{n} z_{n} 0 {pwl(zpoints)}',
                  f'Bcmd_{n} cmd_{n} 0 V=({positive})-({negative})',
                  f'Shold_{n} cmd_{n} {hold} clk 0 SAMP',
                  f'Chold_{n} {hold} 0 1f IC=0',
                  f'Rh_{n} {hold} 0 1e18']
        width=physical.packet_width_ns*sample['width_gains'][k]
        if not 0.04 < width < 2.5:
            raise ValueError('pulse width outside supported nonoverlap timing')
        starts=[p['sample_time_ns']+0.6 for p in case['phases']]
        lines.append(f'Vgate_{n} g_{n} 0 {pwl(pulses(starts,width))}')
        current=charge/(physical.packet_width_ns*1e-9-0.01e-9)
        gc,gl,gr=sample['transfer_gains'][k]
        ctrl=f'V({hold})*V(g_{n})'
        lines += [f'Bremove_{n} {n} 0 I={q*current*gc:.14g}*{ctrl}',
                  f'Bleft_{n} 0 {node_name(k-1)} I={current*gl:.14g}*{ctrl}',
                  f'Bright_{n} 0 {node_name(k+1)} I={current*gr:.14g}*{ctrl}']
    save=[f'V({node_name(k)})' for k in all_nodes]+[f'V(h_{node_name(k)})' for k in active]
    lines += ['.save '+' '.join(save), f".tran 0 {case['stop_time_ns']:.9g}n 0 {physical.max_step_ns:.9g}n UIC", '.end']
    return '\n'.join(lines)+'\n'


def read_ascii_raw(path):
    data=Path(path).read_bytes()
    # XVII normally emits ASCII; some releases emit UTF-16 headers/text.
    text=data.decode('utf-16') if data[:2] in (b'\xff\xfe',b'\xfe\xff') else data.decode('utf-8-sig')
    head, body=text.split('Values:',1)
    varlines=head.rsplit('Variables:',1)[1].strip().splitlines()
    names=[line.split()[1].lower() for line in varlines]
    numbers=[]
    for line in body.splitlines():
        parts=line.split()
        if not parts:
            continue
        numbers.append(float(parts[-1]))
    if len(numbers)%len(names):
        raise ValueError('incomplete raw record')
    values=np.array(numbers).reshape(-1,len(names))
    return names, values


def classify(voltage, unit, radius):
    """Final/phase persistent equivalence region; never recover a broad c."""
    if voltage >= (radius+0.5)*unit:
        return 'POS_OVERFLOW'
    if voltage < -(radius+0.5)*unit:
        return 'NEG_OVERFLOW'
    return max(-radius,min(radius,math.floor(voltage/unit+0.5)))


def expected_class(value, radius):
    return 'POS_OVERFLOW' if value>radius else 'NEG_OVERFLOW' if value < -radius else value


def verify(case, physical, names, values, sample=None):
    lookup={n:i for i,n in enumerate(names)}
    ts=values[:,0]
    unit=physical.unit_voltage_v/physical.scale
    all_nodes=sorted(case['active_nodes']+case['passive_boundary_receivers'])
    radius=case['q']-1
    def read(tns, column):
        # Interpolate instrumentation only; classify via fixed declared regions.
        return float(np.interp(tns*1e-9,ts,values[:,lookup[column.lower()]]))
    def detected_voltage(k, tns):
        voltage=read(tns,f'v({node_name(k)})')
        if sample is not None and k in case['active_nodes']:
            points=[(0.0,0.0)]
            for p,z in zip(case['phases'],sample['sense_noise_v'][k]):
                points.extend([(p['sample_time_ns']-0.5,z),(p['observe_time_ns'],z)])
            tt,zz=zip(*sorted(dict(points).items()))
            voltage += float(np.interp(tns,tt,zz))-sample['offsets_v'][k]
        return voltage
    checks=[]; phase_failures=0; phase_states=[]
    trigger_events=unit_packets=0
    observed_last_active_ns=None
    for p in case['phases']:
        row={}; mismatch=0
        for k in all_nodes:
            region=classify(detected_voltage(k,p['observe_time_ns']),unit,radius)
            expected=expected_class(p['state'].get(k,p['state'].get(str(k),0)),radius)
            row[k]=region
            mismatch += region != expected
        phase_failures += bool(mismatch)
        phase_states.append({'time_ns':p['observe_time_ns'],'regions':row,'wrong_regions':mismatch})
        for k in case['active_nodes']:
            latch=read(p['sample_time_ns']+0.4,f'v(h_{node_name(k)})')
            # Instrument only the small declared command alphabet; no broad quotient.
            m=int(round(latch))
            if m:
                trigger_events += 1; unit_packets += abs(m)
                observed_last_active_ns=p['observe_time_ns']
    wrong_total=state_count=0
    max_residual=0.0
    exact_class_failures=0
    for cp in case['checkpoints']:
        observed={}; residuals={}; wrong=0
        target={int(k):v for k,v in cp['state'].items()}
        for k in all_nodes:
            voltage=read(cp['time_ns'],f'v({node_name(k)})')
            observed[k]=classify(detected_voltage(k,cp['time_ns']),unit,radius)
            wrong += observed[k] != expected_class(target.get(k,0),radius)
            # Residual is continuous diagnostic, not used to compute output.
            residuals[k]=voltage/unit-target.get(k,0)
        wrong_total += wrong; state_count += len(all_nodes)
        max_residual=max(max_residual,max(map(abs,residuals.values()),default=0))
        persistent=all(isinstance(v,int) for v in observed.values())
        observed_pair=pair_json(observed,case['q']) if persistent else None
        exact_class_failures += observed_pair != cp['exact_pair']
        checks.append({'time_ns':cp['time_ns'], 'observed_regions':observed,
            'expected_state':target, 'exact_pair_from_classes':observed_pair,
            'expected_exact_pair':cp['exact_pair'], 'wrong_regions':wrong,
            'residual_units_instrumentation':residuals,
            'boundary_requires_extension':cp['boundary_requires_extension'],
            'all_persistent':persistent})
    state_columns=[lookup[f'v({node_name(k)})'.lower()] for k in all_nodes]
    peak=float(np.abs(values[:,state_columns]).max())
    expected_triggers=sum(cp['rewrite_events'] for cp in case['checkpoints'])
    return {'wrong_classes':wrong_total,'class_observations':state_count,
        'network_class_failure':bool(wrong_total or exact_class_failures),
        'exact_class_value_failures':exact_class_failures,
        'phase_class_failures':phase_failures, 'phase_region_traces':phase_states,
        'max_abs_residual_units':max_residual, 'peak_abs_voltage_v':peak,
        'headroom_violation':peak>physical.headroom_v,
        'observed_trigger_events':trigger_events,'observed_unit_packets':unit_packets,
        'expected_trigger_events':expected_triggers,
        'induced_extra_trigger_events':trigger_events-expected_triggers,
        'last_active_phase_time_ns':observed_last_active_ns,
        'observation_stop_time_ns':case['stop_time_ns'],
        'boundary_requires_extension':any(cp['boundary_requires_extension'] for cp in case['checkpoints']),
        'checkpoints':checks}


def run_spice(case, directory, physical=None, seed=0, exe=DEFAULT_EXE, retain_raw=False):
    physical=physical or Physical()
    directory=Path(directory).resolve(); directory.mkdir(parents=True,exist_ok=True)
    sample=realization(case,physical,seed)
    deck=directory/'normalizer.cir'
    deck.write_text(netlist(case,physical,sample),encoding='ascii')
    (directory/'case.json').write_text(json.dumps(case,indent=2)+'\n')
    started=time.perf_counter()
    command=[str(exe),'-b','-ascii',str(deck)]
    result=subprocess.run(command,capture_output=True,timeout=45,
        creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess,'CREATE_NO_WINDOW') else 0)
    raw=deck.with_suffix('.raw'); log=deck.with_suffix('.log')
    if result.returncode or not raw.exists():
        detail=log.read_text(errors='replace') if log.exists() else result.stderr.decode(errors='replace')
        raise RuntimeError(f'LTspice failed ({result.returncode}): {detail[:1500]}')
    names,values=read_ascii_raw(raw)
    report=verify(case,physical,names,values,sample)
    report.update({'q':case['q'],'scheme':case['scheme'],'case':case['name'],
        'active_nodes':len(case['active_nodes']), 'seed':seed,'physical':asdict(physical),
        'elapsed_wall_s':time.perf_counter()-started,
        'simulator_executable':str(exe),'netlist':str(deck),
        'realization':sample})
    (directory/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    if not retain_raw:
        # Only this generated run's raw file is removed; netlist, seed, golden
        # trace, parsed verification and simulator log remain reproducible.
        assert raw.parent == directory
        raw.unlink()
    return report
