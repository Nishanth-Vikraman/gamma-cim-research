"""Direct seeded transient boundary tests of the physical comparator.
Intrinsic transient MOS thermal noise is NOT present in deterministic ngspice
transients. The extra input Gaussian samples below are an explicitly labelled
noise-injection sensitivity; PDK local-device mismatch is independently active.
"""
import itertools,math
from pathlib import Path
import numpy as np
from signed_binary_kcl import BinaryCircuit,OUT,write_csv
from gamma_sky130 import pwl,dump,compress_trace
import signed_binary_backend as backend

def build(seed,cm,inject):
    rng=np.random.default_rng(seed);cfg={'pdk_mm':True,'seed':seed};c=BinaryCircuit(cfg)
    delta=np.linspace(-.020,.020,41);period=20;points=[(0,cm+delta[0])];clocks=[(0,0)];checks=[]
    for j,dd in enumerate(delta):
        s=10+j*period;actual=float(dd+(rng.normal(0,.0002) if inject else 0))
        points += [(s-1,cm+actual),(s+period-1.1,cm+actual)]
        clocks += [(s+3,0),(s+3.2,1.8),(s+13.8,1.8),(s+14,0)]
        checks.append({'sample_ns':s+12,'delta_v':float(dd),'actual_delta_v':actual})
    c.source('input','input',pwl(points),'external_fixture');c.source('ref','ref',str(cm),'external_fixture')
    c.source('clockpin','clockpin',pwl(clocks),'clock_pin')
    cb=c.invnode('clockpin','cb','clock');clock=c.invnode(cb,'clk','clock')
    op,on=c.comparator('input','ref',clock,'testcmp')
    c.save={f'v({op})',f'v({on})','v(input)','v(ref)'};stop=10+len(delta)*period
    return c.finish(f'.tran .5n {stop}n'),{'kind':'boundary','seed':seed,'cm':cm,'inject':inject,'checks':checks,
        'op':op,'on':on,'stop_ns':stop,'design':{'temp_c':25}}

def processor(result,names,x,directory):
    meta=result['meta'];ix={n:i for i,n in enumerate(names)};rows=[];t=x[:,0]
    for ch in meta['checks']:
        q=float(np.interp(ch['sample_ns']*1e-9,t,x[:,ix['v('+meta['op']+')']]))
        qn=float(np.interp(ch['sample_ns']*1e-9,t,x[:,ix['v('+meta['on']+')']]))
        digit=int(q>.9);expect=int(ch['delta_v']>0)
        rows.append({**ch,'seed':meta['seed'],'cm':meta['cm'],'injected_noise':meta['inject'],'output_v':q,
            'output_n_v':qn,'decision':digit,'expected':expect,'wrong':digit!=expect,'tie':ch['delta_v']==0,
            'metastable':not ((q>1.6 and qn<.2) or (qn>1.6 and q<.2))})
    write_csv(Path(directory)/'boundaries.csv',rows)
    ones=[r['delta_v'] for r in rows if r['decision']];zeros=[r['delta_v'] for r in rows if not r['decision']]
    off=(min(ones)+max(zeros))/2 if ones and zeros else None
    result['observations']=True;result['metrics']={'seed':meta['seed'],'cm':meta['cm'],'inject':meta['inject'],
        'offset_v':off,'offset_censored':off is None,'wrong_non_tie':sum(r['wrong'] for r in rows if not r['tie']),
        'non_tie_trials':sum(not r['tie'] for r in rows),'metastable_trials':sum(r['metastable'] for r in rows)}
    compress_trace(Path(directory)/'output.raw');return result

def main():
    backend.BASE=OUT/'boundary_backend';backend.BASE.mkdir(parents=True,exist_ok=True);cases=[]
    for j,(cm,seed,inject) in enumerate(itertools.product((.675,.9,1.125),range(261000,261008),(False,True))):
        cases.append((*build(seed,cm,inject),OUT/'spice/boundary'/f'{j:03d}'))
    rs=backend.batch_execute(cases,2,processor);write_csv(OUT/'boundary.csv',[r['metrics'] for r in rs])
    print('Direct boundary networks',len(rs),'comparisons',len(rs)*41,flush=True)
if __name__=='__main__':main()
