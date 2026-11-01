from fractions import Fraction as F
import random
import pytest
from gamma_native_model import value_pair,NativeFloat,shift
from gamma_regen_golden import pair_sign,exact_nearest_integer,finite_word,finite_merge,event_completion,word_case
from gamma_regen_spice import primitive_deck,chain_deck,word_deck,RegenPhysical,run,profile,fast_raw


def test_exact_sign_and_rounding():
    assert pair_sign((-4,1))==-1
    assert pair_sign((-3,1))==1
    assert exact_nearest_integer((F(5,2),0),5)==2
    assert exact_nearest_integer((F(7,2),0),5)==4
    assert exact_nearest_integer((F(-5,2),0),5)==-2
    assert exact_nearest_integer((-3,1),5)==1


@pytest.mark.parametrize('policy',['truncate','nearest'])
def test_finite_ledger(policy):
    rng=random.Random(610)
    for _ in range(100):
        word={k:rng.randrange(-8,9) for k in range(-5,6)}
        result=finite_word(word,exponent=3,policy=policy)
        output=NativeFloat(result['state'],result['exponent'],q=4).exact()
        residual=tuple(map(F,result['exact_discarded_or_rounding_residual']))
        assert tuple(x+y for x,y in zip(output,residual))==value_pair(shift(word,3),4)
        assert all(-2<=k<=2 and abs(d)<=3 for k,d in result['state'].items())


def test_boundary_and_merge():
    r=finite_word({3:1,-2:1})
    assert r['state']=={2:1} and r['exponent']==1
    assert tuple(map(F,r['exact_discarded_or_rounding_residual']))==value_pair({-2:1},4)
    r=finite_word({-1:1},renormalize_small=True)
    assert r['exponent']==-3 and r['state']=={2:1}
    merged=finite_merge(NativeFloat({0:1},4,q=4),NativeFloat({0:2},0,q=4),window=(-5,2))
    assert list(map(F,merged['exact_discarded_or_rounding_residual']))==[0,0]


def test_event_queue_conservation_and_quiescence():
    rng=random.Random(13)
    for _ in range(80):
        raw={k:rng.randrange(-20,21) for k in range(-3,4)}
        result=event_completion(raw,range(-3,4))
        assert result['all_persistent']
        assert value_pair(result['state'],4)==value_pair(raw,4)
        assert result['counts']['requests']==result['counts']['acknowledgments']


def test_connected_netlist_contract():
    deck,meta=chain_deck(depth=2)
    assert '*V(s0out)' in deck
    assert 'V(s0h)' in deck and 'V(s0bh)' in deck
    assert len(meta['oracle']['stages'])==2
    deck,_=word_deck()
    assert '2*V(r3out)' in deck


@pytest.mark.parametrize('mapping',['single','differential','charge'])
def test_actual_spice_canonical_levels(tmp_path,mapping):
    report=run(*primitive_deck(-2,mapping),tmp_path/'p')
    assert report['observations'][0]['correct']
    assert abs(report['observations'][0]['residual_units'])<1e-4


def test_actual_spice_connected_word_regeneration(tmp_path):
    p=RegenPhysical(stage_drift_units=.35)
    direct=run(*word_deck(physical=p,mode='none'),tmp_path/'n')
    restored=run(*word_deck(physical=p),tmp_path/'r')
    assert not direct['final_algebraically_correct']
    assert restored['final_algebraically_correct']
    assert all(o['correct'] for o in restored['observations'])


def test_fast_raw_record_indices(tmp_path):
    path=tmp_path/'probe.raw'
    path.write_text('No. Variables: 2\nVariables:\n0 time time\n1 V(out) voltage\nValues:\n0 0\n 0.025\n1 1e-9\n 0.05\n')
    names,array=fast_raw(path)
    assert names==['time','v(out)']
    assert array.tolist()==[[0,.025],[1e-9,.05]]


def test_actual_spice_feedthrough_and_latch_stability(tmp_path):
    # Previously exposed inconsistent feed-cap IC and trapezoidal chatter.
    for i,(d,mapping,seed) in enumerate([(-3,'single',20133),(-2,'differential',20058)]):
        deck,meta=primitive_deck(d,mapping,profile('low'),seed)
        assert 'IC=.6' in deck and 'method=gear' in deck
        assert 'if(V(rw)>.1,' in deck
        report=run(deck,meta,tmp_path/str(i),timeout=10)
        assert report['observations'][0]['correct']
