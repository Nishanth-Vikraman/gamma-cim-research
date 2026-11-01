"""Confidence accounting and sampled charge-domain scaling study.

This accelerated model is NOT SPICE or transistor evidence. It is a separately
labeled model of sampled coarse comparators and fixed packets, used for rare
events beyond the affordable SPICE run count. It excludes continuous leakage
and current-noise waveforms; those are covered only by actual SPICE sweeps.
"""
from __future__ import annotations
from dataclasses import asdict
import math
import numpy as np
from scipy.stats import beta
from gamma_protocols import bands
from gamma_spice import Physical, classify, node_name


def confidence(k,n,alpha=0.05):
    if not 0 <= k <= n or n <= 0:
        raise ValueError('valid binomial count required')
    upper = -math.expm1(math.log(alpha)/n) if k == 0 else (1.0 if k==n else float(beta.ppf(1-alpha,k+1,n-k)))
    lower = 0.0 if k==0 else (alpha**(1/n) if k==n else float(beta.ppf(alpha,k,n-k+1)))
    return {'trials':n,'failures':k,'estimate':k/n,
        'one_sided_95pct_lower':lower,'one_sided_95pct_upper':upper,
        'confidence_note':'Pointwise one-sided exact Clopper-Pearson bounds on independent network trials; not simultaneous bounds across a scale sweep.'}


def sampled_predict(case,physical,sample):
    """Same discrete packet law and sampled offsets as SPICE for a *single*
    realization. Useful cross-check only, not used to verify SPICE outputs.
    Continuous current noise/leakage intentionally not predicted here.
    """
    if physical.current_noise_sigma_a or physical.leakage_tau_ns < 1e12:
        return None
    active=case['active_nodes']
    nodes=sorted(active+case['passive_boundary_receivers'])
    state={k:0.0 for k in nodes}
    unit=physical.unit_voltage_v/physical.scale
    outputs=[]; expected_index=-1
    phase_index=0
    for cp in case['checkpoints']:
        ix=cp['injection_index']
        for event in case['stimuli']:
            # Injection groups are selected by their time between prior/future checkpoints.
            lower=case['checkpoints'][ix-1]['time_ns'] if ix else 0
            if lower < event['time_ns'] < cp['time_ns']:
                k=event['node']
                state[k] += event['delta_units']*unit*sample['injection_gains'][k]/sample['cap_gains'][k]
        for phase in [p for p in case['phases'] if p['injection_index']==ix]:
            changes={k:0.0 for k in nodes}
            for k in active:
                signal=state[k]+sample['sense_noise_v'][k][phase_index]-sample['offsets_v'][k]
                m=0; previous=0
                for level in bands(case['scheme']):
                    threshold=(case['q']*level-0.5)*unit
                    m += (level-previous)*((signal>threshold)-(signal < -threshold))
                    previous=level
                factor=(physical.packet_width_ns*sample['width_gains'][k]-0.01)/(physical.packet_width_ns-0.01)
                gc,gl,gr=sample['transfer_gains'][k]
                changes[k] -= case['q']*unit*m*gc*factor/sample['cap_gains'][k]
                changes[k-1] += unit*m*gl*factor/sample['cap_gains'][k-1]
                changes[k+1] += unit*m*gr*factor/sample['cap_gains'][k+1]
            for k,dv in changes.items():
                state[k] += dv
            phase_index += 1
        outputs.append({k:classify(state[k]+(sample['sense_noise_v'][k][phase_index-1]-sample['offsets_v'][k] if k in active else 0),unit,case['q']-1) for k in nodes})
    return outputs


def charge_monte_carlo(q,scale,trials=10000,seed=0,combined=False,batch=20000):
    """Five active cells plus passive guards; one simultaneous signed burst.

    Persistent expected state is derived once by the exact Python protocol,
    not inferred from sampled observations. Array updates are physical voltage
    integration of fixed charge packets, never post-readout normalization.
    """
    from gamma_protocols import settle
    raw={-1:q,0:2*q,1:-q}
    gold=settle(raw,q,'unit',active=(-2,2))['state']
    target=np.array([gold.get(k,0) for k in range(-3,4)])
    rng=np.random.default_rng(seed)
    failures=class_errors=0
    max_residual=0.0
    sigma_voltage=0.0002
    unit=0.1/scale
    rounds=12
    headroom_failures=0
    for start in range(0,trials,batch):
        n=min(batch,trials-start)
        caps=1+rng.normal(0,0.005,(n,7)) if combined else np.ones((n,7))
        offsets=(rng.normal(0,0.0001,(n,1))+rng.normal(0,0.0001,(n,5))) / unit if combined else np.zeros((n,5))
        gains=1+rng.normal(0,0.005,(n,7)) if combined else np.ones((n,7))
        transfers=1+rng.normal(0,0.005,(n,5,3)) if combined else np.ones((n,5,3))
        widths=(1+rng.normal(0,0.005,(n,5))-0.01)/0.99 if combined else np.ones((n,5))
        state=np.zeros((n,7))
        state[:,2]=q*gains[:,2]/caps[:,2]
        state[:,3]=2*q*gains[:,3]/caps[:,3]
        state[:,4]=-q*gains[:,4]/caps[:,4]
        headroom_failures += int(np.any(np.abs(state*unit)>0.6,axis=1).sum())
        for _ in range(rounds):
            noise=rng.normal(0,sigma_voltage/unit,(n,5))
            signal=state[:,1:6]+noise-offsets
            command=(signal>(q-0.5)).astype(np.int8)-(signal < -(q-0.5)).astype(np.int8)
            packet=command*widths
            change=np.zeros_like(state)
            change[:,1:6] -= q*packet*transfers[:,:,0]/caps[:,1:6]
            change[:,0:5] += packet*transfers[:,:,1]/caps[:,0:5]
            change[:,2:7] += packet*transfers[:,:,2]/caps[:,2:7]
            state += change
        detected=state.copy()
        detected[:,1:6] += noise-offsets
        # Finite digit alphabet plus explicit overflow labels, not a broad ADC.
        classes=np.clip(np.floor(detected+0.5),-q,q).astype(int)
        wrong=classes!=target
        failures += int(np.any(wrong,axis=1).sum())
        class_errors += int(wrong.sum())
        max_residual=max(max_residual,float(np.abs(state-target).max()))
    return {'q':q,'scale':scale,'source':'sampled charge-domain surrogate, NOT SPICE',
        'profile':'voltage_noise_plus_relative_mismatch' if combined else 'voltage_noise_only',
        **confidence(failures,trials), 'wrong_node_classes':class_errors,
        'node_class_observations':trials*7,'max_residual_units':max_residual,
        'headroom_failures_at_injection':headroom_failures,
        'sense_sigma_v':sigma_voltage,'spacing_v':unit,'spacing_over_sigma':unit/sigma_voltage,
        'relative_sigma':0.005 if combined else 0,
        'bank_and_comparator_offset_sigma_v':0.0001 if combined else 0,
        'assumptions':'12 synchronous unit phases, ideal sample/hold and packet shape, no continuous leakage/current noise or temporal circuit dynamics; pointwise confidence only.'}
