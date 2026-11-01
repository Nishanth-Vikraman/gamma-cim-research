"""Accounting checks: interpolated window edges and physical supply ownership."""
import numpy as np
import pytest
from sky130_self_topple import integrate, build

def test_energy_windows_conserve_signed_linear_supply_work():
    # 1.8 V rail, linearly rising load through a zero-current crossing.
    t=np.array([0,.17,.56,.91,1.0])*1e-6
    power=1.8*(t/1e-6-.4)*1e-6
    a,b=.13e-6,.94e-6
    expected=1.8e-6*((b*b-a*a)/(2e-6)-.4*(b-a))
    assert integrate(t,power,a,b)==pytest.approx(expected,rel=1e-12)
    chunks=[a,.29e-6,.77e-6,b]
    assert sum(integrate(t,power,x,y) for x,y in zip(chunks,chunks[1:]))==pytest.approx(expected,rel=1e-12)

def test_guarded_cell_has_no_ideal_detector_or_current_transfer():
    deck,meta=build({'guarded':True,'cap_ff':100,'signs':[1]})
    devices=[s for s in deck.splitlines() if s and not s.startswith(('*','.'))]
    # Only the independent additive stimulus is an ideal current source.
    assert [s.split()[0] for s in devices if s.startswith('I')]==['Iinput']
    assert not any(s.startswith('B') for s in devices)
    assert not any('PULSE' in s for s in devices)
    assert {p['component'] for p in meta['ports'].values()}=={'detector','logic','bias','transfer','reference','initialization'}
    assert len(meta['mos'])==101

def test_invalid_core_supply_remains_rejected():
    with pytest.raises(ValueError,match='1.95'):
        build({'vdd':1.98})
