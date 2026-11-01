"""Independent exact-value / implementation-contract checks for binary KCL."""
import pytest
from fractions import Fraction as F
from signed_binary_kcl import workload,golden,build,NS
@pytest.mark.parametrize('n',NS)
@pytest.mark.parametrize('pattern',['aligned','aligned_negative','dense','sparse','cancellation','int8'])
def test_exact_conservation(n,pattern):
    w=workload(n,pattern);g=golden(w,12,w['Q']+2)
    assert g['error']==0
    assert all(abs(v['r'])<=.5 for v in g['stages'])
    assert w['D']&(w['D']-1)==0
    assert all(w['D']>=2*v['terms'] for v in w['counts'])
    assert sum(v['terms'] for v in w['counts'])==w['product_planes']*n
@pytest.mark.parametrize('n',NS)
def test_finite_error_bound_and_plane_cost(n):
    w=workload(n,'aligned');g=golden(w)
    assert abs(g['error'])<=g['max_truncation_error']
    assert w['physical_product_connections']==72*n
    assert g['max_truncation_error']==w['D']/32

def test_no_behavioral_datapath():
    deck,meta=build({'n':64,'pattern':'dense'})
    lines=deck.splitlines()
    assert not any(l.startswith(('B','E','F','G','H')) for l in lines)
    assert not any('ternary' in l or 'if(' in l for l in lines)
    assert len(meta['stages'])==16
    assert all('sky130' in l for l in lines if l.startswith('X'))
    assert sum(r['m'] for r in meta['resistors'] if r['kind']=='ballast')>0
    assert all(s['previous'][0]=='tap8' if j==0 else s['previous'][0]==meta['stages'][j-1]['outp'] for j,s in enumerate(meta['stages']))

def test_zero_margin_legal_not_claimed_safe():
    g=golden(workload(64,'aligned'))
    assert min(v['margin'] for v in g['stages'])==0
