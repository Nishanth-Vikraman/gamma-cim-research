"""Exact bounded q4 quotient policies, dirty requests and streaming dot counts.

No floating beta and no free analog overflow oracle. All scheduled dirty
coefficient visits, including zero quotients, count as physical invocations.
"""
from collections import defaultdict,Counter
from fractions import Fraction
from functools import lru_cache
import csv,json,random,math
from pathlib import Path
from gamma_native_model import normalize,convolve,value_pair,power_pair

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/sar_normalizer_phase'
LIMITS=(7,15,31)

def quotient(c):return (1 if c>=0 else -1)*(abs(c)//4)

def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:path.write_text('');return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def golden():
    rows=[]
    for limit in LIMITS:
        bits=(limit//4).bit_length();states=2*(limit//4)+1
        for c in range(-limit,limit+1):
            m=quotient(c);r=c-4*m
            assert -3<=r<=3 and value_pair({0:r,-1:m,1:m},4)==(c,0)
            for scheme in ('bank','sar','binary_residual'):
                rows.append({'limit':limit,'coefficient':c,'quotient':m,'golden_quotient':m,'residual':r,
                  'quotient_states':states,'magnitude_bits':bits,'scheme':scheme,
                  'comparisons_per_invocation':2*(limit//4) if scheme=='bank' else bits+1,
                  'center_charge_q':-4*m,'neighbor_charge_q':m,'exact_pair_a':c,'exact_pair_b':0})
    write_csv(OUT/'transient_range_sweep.csv',rows);return rows

@lru_cache(None)
def encode(n):return normalize({0:n},4)[0]

def count_dot(pairs,limit,unit_control=False,order='stream'):
    """Deterministic ingress schedule. L7 splits each digit-product into <=3
    unit injections (extra input-drive work is not free). L15 admits one full
    <=9 digit-product per batch; L31 admits three. No state-dependent trigger.
    Each batch starts from persistent states. Dirty neighbors are revisited
    after a synchronous rewrite even if their resulting quotient is zero.
    """
    ingress=3 if limit==7 else 9
    batch=1 if limit in (7,15) else 3
    state=defaultdict(int);invocations=nonzero=units=rounds=completions=0
    peak=0;hist=Counter();chist=Counter();checks=0;contributions=[];dirty_peak=0;input_slices=0
    def drain(items):
        nonlocal invocations,nonzero,units,rounds,completions,peak,dirty_peak
        dirty=set()
        for k,v in items:state[k]+=v;dirty.add(k)
        peak=max(peak,max((abs(state[k]) for k in dirty),default=0))
        assert all(abs(state[k])<=limit for k in dirty)
        while dirty:
            dirty_peak=max(dirty_peak,len(dirty));delta=defaultdict(int);receivers=set()
            for k in sorted(dirty):
                invocations+=1;m=quotient(state[k])
                chist[state[k]]+=1
                if unit_control:m=(m>0)-(m<0)
                hist[abs(m)]+=1
                if m:
                    nonzero+=1;units+=abs(m);delta[k]-=4*m;delta[k-1]+=m;delta[k+1]+=m
                    receivers.update((k-1,k,k+1))
            if not delta:break
            rounds+=1;dirty=receivers
            for k,v in delta.items():
                state[k]+=v
            peak=max(peak,max((abs(state[k]) for k in dirty),default=0))
            assert all(abs(state[k])<=limit for k in dirty)
        completions+=1
    inputs=[]
    for a,b in pairs:
        for i,x in sorted(encode(a).items()):
            for j,y in sorted(encode(b).items()):
                v=x*y
                while v:
                    part=(1 if v>0 else -1)*min(abs(v),ingress);v-=part
                    inputs.append((i+j,part));input_slices+=1
    if order.startswith('position_'):inputs.sort(key=lambda v:v[0])
    for k,part in inputs:
        if order.startswith('position_') and contributions and contributions[-1][0]!=k:
            drain(contributions);contributions=[]
        if order=='position_envelope' and contributions and sum(abs(v) for _,v in contributions)+abs(part)>limit-3:
            drain(contributions);contributions=[]
        contributions.append((k,part))
        if order!='position_envelope' and len(contributions)==batch:drain(contributions);contributions=[]
    if contributions:drain(contributions)
    expected=sum(a*b for a,b in pairs)
    assert value_pair({k:v for k,v in state.items() if v},4)==(expected,0)
    assert all(abs(v)<=3 for v in state.values())
    return {'N':len(pairs),'limit':limit,'input_order':order,'policy':'bounded_unit_control' if unit_control else 'bounded_bulk',
      'invocations_including_zero':invocations,'bulk_nonzero_transfers':nonzero,'unit_rewrite_equivalents':units,
      'hist_abs_m':dict(sorted(hist.items())),'hist_coefficient':dict(sorted(chist.items())),
      'max_transient_coefficient':peak,'synchronous_rounds':rounds,
      'completion_events':completions,'input_charge_slices':input_slices,'ingress_slices_per_batch':batch if order!='position_envelope' else None,
      'known_ingress_absolute_budget':limit-3 if order=='position_envelope' else None,
      'total_absolute_moved_charge_q':6*units,'dirty_bank_peak':dirty_peak,
      'final_active_span':max(state)-min(state)+1 if state else 1,'exact_dot':expected,'exact':True}

def historical_products(pairs,limit):
    """Same historical scalar-pair suite, direct raw convolution. Oversized raw
    coefficients invalidate the declared window; never silently clamp them.
    """
    from gamma_protocols import settle
    valid=0;events=units=0;peak=0
    for a,b in pairs:
        raw=convolve(encode(a),encode(b));p=max(map(abs,raw.values()),default=0);peak=max(peak,p)
        if p>limit:continue
        r=settle(raw,4,'bulk')
        if r['peak_units']>limit:continue
        valid+=1;events+=r['trigger_events'];units+=r['unit_rewrites']
    return {'limit':limit,'cases':len(pairs),'bounded_cases':valid,'raw_peak':peak,
      'nonzero_bulk_events_over_bounded_cases':events,'unit_equivalents_over_bounded_cases':units,
      'direct_product_schedule_legal_for_all':valid==len(pairs)}

def main():
    OUT.mkdir(parents=True,exist_ok=True);golden();rng=random.Random(0)
    pairs=[(rng.randint(-128,127),rng.randint(-128,127)) for _ in range(8192)]
    rows=[];details=[]
    for limit in LIMITS:
      for order in ('stream','position_batched','position_envelope'):
        for n in (128,512,2048,8192):
            b=count_dot(pairs[:n],limit,order=order);u=count_dot(pairs[:n],limit,True,order=order)
            details.extend([b,u])
            for scheme in ('bank','sar','binary_residual'):
                comps=2*(limit//4) if scheme=='bank' else 1+(limit//4).bit_length()
                rows.append({**{k:v for k,v in b.items() if k not in ('hist_abs_m','hist_coefficient')},'scheme':scheme,
                  'comparisons':b['invocations_including_zero']*comps,
                  'unit_control_invocations':u['invocations_including_zero'],
                  'unit_control_nonzero_transfers':u['bulk_nonzero_transfers'],
                  'nonzero_transfer_reduction':u['bulk_nonzero_transfers']/max(1,b['bulk_nonzero_transfers']),
                  'invocation_reduction':u['invocations_including_zero']/max(1,b['invocations_including_zero']),
                  'hist_abs_m_json':json.dumps(b['hist_abs_m'],sort_keys=True),
                  'scope':'Direct streamed dot, deterministic bounded ingress; all dirty m=0 requests charged. Physical queue/global completion and sliced-input driver overhead not yet built.'})
            print('exact',limit,order,n,b['invocations_including_zero'],b['bulk_nonzero_transfers'],b['hist_abs_m'],flush=True)
    write_csv(OUT/'quotient_event_counts.csv',rows)
    (OUT/'exact_workload.json').write_text(json.dumps({'seed':0,'pairs':pairs,'details':details,
      'historical_product_audit':[historical_products(pairs[:128],l) for l in LIMITS]},indent=2))

if __name__=='__main__':main()
