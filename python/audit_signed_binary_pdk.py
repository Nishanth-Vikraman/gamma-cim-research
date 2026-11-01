"""Replay native PDK wrappers versus selected verbatim-bin wrappers."""
import numpy as np
from pathlib import Path
from signed_binary_kcl import BinaryCircuit,OUT,write_csv
from gamma_sky130 import Circuit,dump,compress_trace
import signed_binary_backend as backend

def build(corner):
    c=BinaryCircuit({'corner':corner,'full_library':True});pairs=[]
    for j,(kind,w) in enumerate((k,w) for k in ('n','p') for w in (.42,.84,1.68,3.36,8.4,16.8)):
        src='0' if kind=='n' else 'dd_row';body=src
        for native in (False,True):
            tag=f'd{j}x{int(native)}';c.source(tag,tag,'.9','audit_bias')
            if native:Circuit.fet(c,tag,kind,tag,'gate',src,body,w)
            else:c.fet(tag,kind,tag,'gate',src,body,w)
        pairs.append((f'i(vd{j}x0)',f'i(vd{j}x1)',kind,w))
    c.source('gate','gate','.9','audit_bias')
    return c.finish('.tran .1n 1n'),{'kind':'audit','seed':20261026,'corner':corner,'pairs':pairs,'design':{'temp_c':25},'stop_ns':1}

def process(result,names,x,directory):
    ix={n:i for i,n in enumerate(names)};rows=[]
    for a,b,kind,w in result['meta']['pairs']:
        va=float(x[-1,ix[a]]);vb=float(x[-1,ix[b]])
        relative=abs(va-vb)/max(abs(vb),1e-15)
        rows.append({'corner':result['meta']['corner'],'kind':kind,'w_um':w,'selected_current_a':va,
            'native_current_a':vb,'relative_error':relative,'passed':relative<1e-7})
    result['metrics']=rows;result['observations']=True
    write_csv(Path(directory)/'comparison.csv',rows);compress_trace(Path(directory)/'output.raw');return result

def main():
    backend.BASE=OUT/'audit_backend';backend.BASE.mkdir(parents=True,exist_ok=True)
    cases=[(*build(co),OUT/'spice/audit'/co) for co in ('tt','ss','ff','sf','fs')]
    rs=backend.batch_execute(cases,2,process);rows=[r for result in rs for r in result['metrics']]
    
    for row in rows:row['passed']=row['relative_error']<1e-7
    write_csv(OUT/'model_audit.csv',rows);dump(OUT/'model_audit.json',{'passed':all(r['passed'] for r in rows),'comparisons':rows})
    print('Native model audit',sum(r['passed'] for r in rows),'/',len(rows),flush=True)
    assert all(r['passed'] for r in rows)
if __name__=='__main__':main()
