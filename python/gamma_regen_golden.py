"""Exact q=4 regeneration oracles, finite words and explicit completion model."""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as F
import heapq
import math
from gamma_native_model import normalize,value_pair,shift,power_pair,NativeFloat
from gamma_protocols import settle

Q=4
D=tuple(range(-3,4))


def pair_sign(pair):
    # beta=2+sqrt(3): sign(C+B sqrt(3)), compared with rational squares.
    a,b=map(F,pair);c=a+2*b
    if not b:return (c>0)-(c<0)
    if not c:return (b>0)-(b<0)
    if (c>0)==(b>0):return (c>0)-(c<0)
    dominant=c*c-3*b*b
    return ((c>0)-(c<0)) if dominant>0 else ((b>0)-(b<0))


def exact_nearest_integer(pair,bound):
    lo,hi=-bound-1,bound+1
    while hi-lo>1:
        m=(lo+hi)//2
        if pair_sign((pair[0]-m,pair[1]))>=0:lo=m
        else:hi=m
    cmp=pair_sign((2*pair[0]-(2*lo+1),2*pair[1]))
    return lo if cmp<0 or (cmp==0 and lo%2==0) else lo+1


def finite_word(word,exponent=0,window=(-2,2),guard_digits=2,policy='truncate',renormalize_small=False):
    """Normalize exactly, shift exponent on high spill, and explicitly book
    discarded value. Nearest rounds to a multiple of beta**window[0], with
    exact ties-to-even; this is NOT a nearest-codebook optimum claim.
    """
    if policy not in ('truncate','nearest') or window[0]>window[1] or guard_digits<0:
        raise ValueError('invalid finite-word policy')
    lo,hi=window
    original=NativeFloat(normalize(word,4)[0],exponent,q=4).exact()
    state,stats=normalize(word,4)
    updates=[];rounding=[]
    needs_extension=bool(state and (min(state)<lo-guard_digits or max(state)>hi+guard_digits))
    for _ in range(16):
        if state and (max(state)>hi or (renormalize_small and max(state)<hi)):
            r=max(state)-hi
            state=shift(state,-r);exponent+=r
            updates.append({'index_shift':-r,'exponent_delta':r})
        discarded={k:v for k,v in state.items() if k<lo}
        kept={k:v for k,v in state.items() if lo<=k<=hi}
        if policy=='nearest' and discarded:
            units=value_pair(shift(discarded,-lo),4)
            increment=exact_nearest_integer(units,sum(map(abs,discarded.values()))+1)
            kept[lo]=kept.get(lo,0)+increment
            rounding.append({'quantum_index':lo,'increment_units':increment})
        state,_=normalize(kept,4)
        if not state or (min(state)>=lo and max(state)<=hi):break
    else:raise RuntimeError('finite rounding/renormalization did not stabilize')
    result=NativeFloat(state,exponent,q=4)
    current=result.exact()
    residual=tuple(a-b for a,b in zip(original,current))
    return {'state':state,'exponent':exponent,'window':list(window),
        'guard_digits':guard_digits,'policy':policy,'exact_original':[str(v) for v in original],
        'exact_output':[str(v) for v in current],'exact_discarded_or_rounding_residual':[str(v) for v in residual],
        'shift_operations':updates,'rounding_operations':rounding,
        'guard_network_extension_needed':needs_extension,
        'reference_normalization_unit_events':stats['unit_events']}


def finite_merge(left,right,window=(-2,2),policy='truncate'):
    if left.q!=4 or right.q!=4:raise ValueError('q=4 lead only')
    e=max(left.exponent,right.exponent)
    raw=defaultdict(int)
    for item in (left,right):
        for k,v in shift(item.digits,item.exponent-e).items():raw[k]+=item.sign*v
    return finite_word(raw,e,window,policy=policy)


def word_case(a=1,b=2,bias=1):
    raw1={0:2*a+b}
    s1=settle(raw1,4,'unit',active=(-2,2))
    raw2={k:2*v for k,v in s1['state'].items()}
    raw2[0]=raw2.get(0,0)-bias
    s2=settle(raw2,4,'unit',active=(-2,2))
    return {'q':4,'inputs':{'a':a,'b':b,'bias':bias},
        'stage1_programmed_weights':[2,1],'stage2_programmed_weights':[2,-1],
        'weight_l1_scale':3,'stage1':s1,'stage2':s2,
        'stage1_exact':[str(v) for v in value_pair(s1['state'],4)],
        'stage2_exact':[str(v) for v in value_pair(s2['state'],4)],
        'active_nodes':list(range(-2,3)),'guard_receivers':[-3,3]}


def chain_case(digit=1,slope=2,depth=16):
    if digit not in D or slope not in (F(2),F(1,2)) or depth<1:
        raise ValueError('declared digit/slope/depth required')
    # The second input is the INITIAL external value, held constant. It is not
    # a Python re-encoding of an intermediate expected/observed result.
    value=F(digit);target=[]
    for stage in range(depth):
        value=slope*value+(1-slope)*digit
        assert value.denominator==1 and int(value) in D
        target.append({'stage':stage+1,'expected_state':{0:int(value)},
            'exact_pair':[str(value),'0']})
    return {'q':4,'digit':digit,'slope':str(slope),'weights':[str(slope),str(1-slope)],
        'weight_l1_scale':str(abs(slope)+abs(1-slope)), 'depth':depth,
        'input_bias':'fixed initial external gamma value, not oracle feedback',
        'stages':target}


def event_completion(raw,nodes,compare_ns=.2,transfer_ns=1,settle_ns=.2,max_operations=100000):
    """Explicit dirty-neighborhood queue, request/ack and hierarchical activity
    model. Revalidation prevents obsolete signed requests from blindly firing.
    This is a digital timing/activity MODEL, not asynchronous SPICE hardware.
    """
    state=defaultdict(int,raw);original=value_pair(raw,4)
    nodes=tuple(sorted(nodes));active=set(nodes)
    queue=[];serial=0;scheduled=set();pending=set();trace=[]
    counts=defaultdict(int);levels=math.ceil(math.log2(max(2,len(nodes))))
    def enqueue(t,k):
        nonlocal serial
        if k in active and k not in scheduled and k not in pending:
            serial+=1;heapq.heappush(queue,(t,serial,'eval',k,0));scheduled.add(k)
            counts['queue_insertions']+=1
    # Initial full scan is charged: no knowledge of a zero answer for free.
    for k in nodes:enqueue(0,k)
    time_ns=0
    while queue:
        time_ns,_,kind,k,sign=heapq.heappop(queue)
        counts['queue_pops']+=1
        counts['comparator_bank_evaluations']+=1
        c=state[k]
        if kind=='eval':
            scheduled.remove(k)
            if abs(c)>=4:
                pending.add(k);counts['requests']+=1
                counts['activity_tree_edge_toggles']+=levels
                serial+=1;heapq.heappush(queue,(time_ns+compare_ns+transfer_ns,serial,'commit',k,1 if c>0 else -1))
                counts['queue_insertions']+=1
        else:
            counts['acknowledgments']+=1
            counts['activity_tree_edge_toggles']+=levels
            pending.remove(k)
            if abs(c)>=4 and (c>0)==(sign>0):
                state[k]-=4*sign;state[k-1]+=sign;state[k+1]+=sign
                counts['trigger_events']+=1;counts['unit_rewrites']+=1
                trace.append({'time_ns':time_ns,'node':k,'sign':sign,'exact_pair':[str(v) for v in value_pair(state,4)]})
            else:counts['cancelled_stale_requests']+=1
            for other in (k-1,k,k+1):enqueue(time_ns+settle_ns,other)
        if counts['queue_pops']>max_operations:raise RuntimeError('event control activity bound exceeded')
    assert value_pair(state,4)==original
    guards={k:v for k,v in state.items() if k not in active and v}
    counts['comparator_evaluations']=2*counts['comparator_bank_evaluations']
    counts['completion_tree_initialization_edges']=max(0,len(nodes)-1)
    counts['completion_guard_evaluations']=len(nodes)
    counts['completion_ack_root_events']=1
    counts['completion_time_ns']=time_ns+compare_ns+settle_ns+levels*.05
    return {'state':{k:v for k,v in state.items() if v},'counts':dict(counts),'trace':trace,
        'boundary_receivers':guards,'all_persistent':all(abs(state[k])<=3 for k in nodes),
        'source':'explicit event-queue/control model; NOT event-driven SPICE',
        'exact_pair':[str(v) for v in original]}
