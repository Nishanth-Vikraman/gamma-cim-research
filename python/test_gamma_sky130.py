"""Checks on the numerical experiment boundary, not on a mirrored ideal circuit."""
import gzip,json
from pathlib import Path
import numpy as np
import pytest
from gamma_sky130 import Circuit,Design,PDK,BASE,boundary_deck,fast_raw,compress_trace
from gamma_sky_batch import parameterize,numeric


def test_binary_and_gzip_waveform_values(tmp_path):
    p=tmp_path/'test.raw'
    expected=np.array([[0,.525],[1e-9,.675]],dtype='<f8')
    p.write_bytes(b'Title: probe\nFlags: real\nNo. Variables: 2\nVariables:\n0 time time\n1 v(out) voltage\nBinary:\n'+expected.tobytes())
    names,x=fast_raw(p)
    assert names==['time','v(out)'] and np.array_equal(x,expected)
    compress_trace(p)
    assert not p.exists() and Path(str(p)+'.gz').exists()
    assert np.array_equal(fast_raw(p)[1],expected)


def test_parameterization_preserves_nodes_and_units():
    deck='test\n.lib "primitive.spice" tt\n.param mc_mm_switch=1 mc_pr_switch=0\n.temp 100\nVdd dd 0 1.8\nRload node 0 100k\nXn node gate 0 0 sky130_fd_pr__nfet_01v8 w=.84 l=.15\n.tran .1n 10n\n.end\n'
    template,p=parameterize(deck)
    assert 'Rload node 0 {' in template
    assert 'Xn node gate 0 0 sky130_fd_pr__nfet_01v8 w={' in template
    assert p['MC_ENABLE']==1
    values=[p[k] for k in p if k!='MC_ENABLE']
    assert np.allclose(values,[1.8,1e5,.84,.15,1e-10,1e-8],rtol=1e-12,atol=0)
    assert numeric('2meg')==2e6 and np.isclose(numeric('10f'),1e-14,rtol=1e-14,atol=0)


def test_real_driver_and_passive_mac_connections():
    deck,meta=boundary_deck(Design(load_r=1e6),sequence=(-3,3))
    assert deck.splitlines()[1].startswith('.lib ')
    assert 'Rloadx driven mac 500000' in deck
    assert 'Rloadbias bias mac 1000000' in deck
    assert 'Xdrivernout driven driverfirst 0 0 sky130_fd_pr__nfet_01v8' in deck
    assert 'Bdriver' not in deck
    assert meta['input_unit_v']==pytest.approx(.025/3)
    assert meta['unit_v']==pytest.approx(.025)
    assert meta['checkpoints'][1]['exact_gamma_pair']==['3','0']


def test_core_supply_guard():
    with pytest.raises(ValueError,match='1.95'):
        Circuit(Design(vdd=1.98))


def test_physical_threshold_decoupling_scales_with_input_area():
    deck,meta=boundary_deck(Design(detector='hybrid',comparator_w=16,comparator_l=.5),sequence=(0,1))
    assert meta['threshold_reference_mim_ff']==pytest.approx(1000*16*.5/.15)
    assert meta['threshold_reference_units']==54
    assert deck.count('Xcthreshold_')==54
    assert 'Xcmpin1 cmpa threshold cmptail dd_comparator sky130_fd_pr__pfet_01v8' in deck
    assert not any(line.startswith('Vthreshold ') for line in deck.splitlines())


def test_executed_model_support_evidence():
    p=BASE/'results/model_validation.json'
    if not p.exists():pytest.skip('Run validate_sky130.py first')
    result=json.loads(p.read_text())
    for key in ('mismatch_disabled_pair_equal','mismatch_enabled_instances_differ','different_seeds_differ',
                'same_seed_replay_identical','batch_vs_literal_temperature_match'):
        assert result[key]
    smoke=json.loads((BASE/'results/smoke_summary.json').read_text())
    assert len(smoke)==20 and all(r['passed'] for r in smoke)
