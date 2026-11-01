"""Exact arithmetic, bounded physical scheduling and no-oracle boundaries."""
from collections import Counter
import numpy as np
import pytest
from gamma_sar_exact import quotient,count_dot
from gamma_native_model import value_pair
from gamma_sar_circuit import build,comparator_build,integral

@pytest.mark.parametrize('limit',[7,15,31])
def test_signed_coarse_sar_matches_exact_quotient_and_pair(limit):
    for c in range(-limit,limit+1):
        m=0
        for j in reversed(range((limit//4).bit_length())):
            trial=m+(1<<j)
            if abs(c)>=4*trial:m=trial
        m=m if c>=0 else -m
        assert m==quotient(c)
        assert value_pair({0:c-4*m,-1:m,1:m},4)==(c,0)
        assert -3<=c-4*m<=3
    assert quotient(-15)==-3

@pytest.mark.parametrize('limit',[7,15,31])
def test_bounded_streams_include_quiet_decisions_and_preserve_exact_dots(limit):
    pairs=[(-128,-128),(127,127),(-127,127),(31,-7),(0,100)]
    for order in ('stream','position_batched','position_envelope'):
        row=count_dot(pairs,limit,order=order)
        assert row['max_transient_coefficient']<=limit
        assert row['invocations_including_zero']>row['bulk_nonzero_transfers']
        assert sum(row['hist_abs_m'].values())==row['invocations_including_zero']
        assert sum(row['hist_coefficient'].values())==row['invocations_including_zero']
        assert row['unit_rewrite_equivalents']==sum(m*n for m,n in row['hist_abs_m'].items())
        assert row['exact_dot']==sum(a*b for a,b in pairs)

@pytest.mark.parametrize('scheme',['sar','binary_residual','bank'])
def test_clock_sources_do_not_encode_expected_quotients(scheme):
    a,ma=build({'scheme':scheme,'limit':15,'states':[4,13],'transfer_kind':'binary_sc'})
    b,mb=build({'scheme':scheme,'limit':15,'states':[-13,0],'transfer_kind':'binary_sc'})
    pins=lambda deck:[s for s in deck.splitlines() if s.startswith('V') and 'pin ' in s]
    assert pins(a)==pins(b)
    for deck,meta in ((a,ma),(b,mb)):
        assert not any(s.startswith(('B','I')) for s in deck.splitlines())
        assert not any(s.startswith('R') for s in deck.splitlines())
        mosnames=[v['name'].lower() for v in meta['mos']]
        assert len(set(mosnames))==len(mosnames)

def test_comparator_both_inputs_have_finite_equal_hold_capacitance():
    deck,meta=comparator_build({'input_ff':200})
    assert 'Xcinput_0 input 0 sky130_fd_pr__cap_mim_m3_1' in deck
    assert 'Xcreference_0 reference 0 sky130_fd_pr__cap_mim_m3_1' in deck
    assert 'Vref reference_stim 0 ' in deck
    assert 'Xcmpin1 cmpa reference cmptail dd_comparator' in deck

def test_energy_partition_counts_signed_supply_work_without_idle_subtraction():
    t=np.array([0,.17,.56,.91,1.])*1e-6;p=1.8*(t/1e-6-.4)*1e-6
    a,b=.13e-6,.94e-6
    analytic=1.8e-6*((b*b-a*a)/(2e-6)-.4*(b-a))
    assert integral(t,p,a,b)==pytest.approx(analytic,rel=1e-12)
    assert integral(t,p,a,.5e-6)+integral(t,p,.5e-6,b)==pytest.approx(analytic,rel=1e-12)

def test_capacitor_mismatch_keeps_layout_instance_count_fixed():
    names=[];dimensions=[]
    for seed in range(51001,51009):
        deck,_=build({'limit':15,'states':[0],'cap_sigma':.01,'seed':seed})
        caps=[s for s in deck.splitlines() if s.startswith('Xc')]
        names.append([s.split()[0] for s in caps]);dimensions.append(caps)
    assert all(v==names[0] for v in names)
    assert any(v!=dimensions[0] for v in dimensions[1:])
