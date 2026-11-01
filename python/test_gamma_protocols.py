from fractions import Fraction
from pathlib import Path
import json
import pytest
from gamma_native_model import value_pair
from gamma_protocols import (SCHEMES, default_cases, export_golden, golden_case,
    phase, protocol_geometry, quotient, settle)
from gamma_spice import DEFAULT_EXE, Physical, classify, netlist, realization, run_spice


@pytest.mark.parametrize('q',[3,4,5])
@pytest.mark.parametrize('scheme',[*SCHEMES,'bulk'])
def test_protocols_exact_signed_parallel_and_declared_alphabets(q,scheme):
    for state in [{0:4*q}, {-1:q,0:-3*q,1:q}, {0:-40*q}]:
        result=settle(state,q,scheme)
        assert value_pair(result['state'],q)==value_pair(state,q)
        assert all(abs(v)<q for v in result['state'].values())
        for t in result['trace']:
            assert value_pair(t['state'],q)==value_pair(state,q)


def test_simultaneous_phase_uses_old_state_not_sequential_updates():
    state,fire=phase({-1:3,0:3,1:3},3)
    assert fire=={-1:1,0:1,1:1}
    assert state=={-2:1,-1:1,0:2,1:1,2:1}


def test_limited_quotient_bands_and_threshold_counts():
    assert quotient(30,3,'unit')==1
    assert quotient(30,3,'cap2')==2
    assert quotient(-30,3,'cap4')==-4
    assert quotient(9,3,'pow2')==2
    assert protocol_geometry(3,'pow2')['overflow_comparators']==6
    with pytest.raises(ValueError):
        quotient(1,6)


def test_golden_checkpoints_are_immutable_snapshots(tmp_path):
    case=golden_case(3,'steps',[(0,1)]*3,phase_budget=6)
    assert case['checkpoints'][0]['state']=={0:1}
    assert case['checkpoints'][1]['state']=={0:2}
    assert case['checkpoints'][2]['state']=={-1:1,1:1}
    p=export_golden(case,tmp_path)
    saved=json.loads(p.read_text())
    assert saved['checkpoints'][0]['state']=={'0':1}
    assert p.with_suffix('.csv').read_text().splitlines()[0]=='time_ns,node,delta_units'


def test_boundary_transfers_never_disappear():
    result=settle({-1:3},3,'unit',active=(-1,1))
    assert result['boundary_receivers']=={-2:1}
    assert result['boundary_requires_extension']
    assert value_pair(result['state'])==value_pair({-1:3})


def test_scaled_netlist_has_physical_comparators_not_oracle_firings():
    c=golden_case(3,'burst',[(0,12)],phase_budget=12)
    p=Physical(scale=5)
    deck=netlist(c,p,realization(c,p,0))
    assert 'if((V(np0)' in deck and '0.05' in deck
    assert 'Shold_' in deck and 'Chold_' in deck and 'Bremove_' in deck
    assert 'expected' not in deck
    assert classify(0.02,0.02,2)==1
    assert classify(0.2,0.02,2)=='POS_OVERFLOW'


@pytest.mark.parametrize('q',[3,4,5])
@pytest.mark.skipif(not DEFAULT_EXE.exists(),reason='LTspice XVII unavailable')
def test_real_spice_unit_cell_matches_exact_regions(q,tmp_path):
    case=golden_case(q,'signed_burst',[(0,2*q),(0,-3*q)],scheme='unit',phase_budget=12)
    report=run_spice(case,tmp_path/str(q),Physical(scale=4))
    assert not report['network_class_failure']
    assert report['phase_class_failures']==0
    assert report['observed_trigger_events']==report['expected_trigger_events']
