"""PDK LTI array-noise measurements; not clocked comparator noise certification."""
import sys,json,math,itertools
from pathlib import Path
import numpy as np
from scipy.special import ndtr
from signed_binary_kcl import BinaryCircuit,OUT,NS,KB,write_csv,build
from gamma_sky130 import dump,compress_trace
import signed_binary_backend as backend

def noise_build(n,mos=True):
    if mos:
        deck,meta=build({'n':n,'pattern':'aligned','array_only':True})
        lines=[]
        for line in deck.splitlines():
            if line.startswith(('.save','.tran')):continue
            if line.startswith('Vrefawakepin ') or line.startswith('Vrowevalpin '):line=' '.join(line.split()[:3])+' 1.8'
            if line.startswith('Vresetpin '):line=' '.join(line.split()[:3])+' 0'
            if line.startswith('Vreference '):line+=' AC 1'
            if line=='.end':lines.extend(['.option sparse','.save onoise_spectrum inoise_spectrum','.noise v(zp5,zn5) Vreference dec 60 1e3 1e11'])
            lines.append(line)
        return '\n'.join(lines)+'\n',{'kind':'noise','seed':20261026,'n':n,'D':meta['workload']['D'],'mos':True,'design':{'temp_c':25},'stop_ns':None}
    cfg={'n':n};c=BinaryCircuit(cfg);c.lines.append('.option sparse');d=1<<(12*n-1).bit_length();count=6*n
    c.source('common','common','.9','noise_bias')
    c.source('stim','stim','DC 1.125 AC 1','noise_bias')
    c.source('negative','negative','DC .675','noise_bias')
    for sign,src in (('p','stim'),('n','negative')):
        row=src
        if mos:
            row='row'+sign;c.tg('noise'+sign,src,row,'one','zero','dd_row',w=.42)
        c.poly('product'+sign,'z'+sign,row,1e4,count)
        c.poly('ballast'+sign,'z'+sign,'common',1e4,d-count,'ballast')
        c.cap('node'+sign,'z'+sign,'0',20)
    c.save={'onoise_spectrum','inoise_spectrum'}
    deck=c.finish('.noise v(zp,zn) Vstim dec 60 1e3 1e11')
    return deck,{'kind':'noise','seed':20261026,'n':n,'D':d,'mos':mos,'design':{'temp_c':25},'stop_ns':None}

def processor(result,names,x,directory):
    meta=result['meta'];ix={v:i for i,v in enumerate(names)}
    if 'onoise_spectrum' not in ix:raise RuntimeError('Missing noise spectrum: '+repr(names))
    f=x[:,0];asd=x[:,ix['onoise_spectrum']];aperture=.24*(.7*400/16)*1e-9
    weight=np.sinc(f*aperture)**2
    sigma=float(np.sqrt(np.trapezoid(asd**2*weight,f)))
    sigma_full=float(np.sqrt(np.trapezoid(asd**2,f)))
    row={'n':meta['n'],'D':meta['D'],'mos':meta['mos'],'sigma_array_v':sigma,
        'sigma_wideband_v':sigma_full,'aperture_ns':aperture*1e9,'sigma_array_units':sigma/.45,
        'frequency_min_hz':float(f[0]),'frequency_max_hz':float(f[-1]),
        'scope':'Full array/reference/row-driver DC noise; excludes dynamic chain' if meta['mos'] else 'Resistor-only ideal-bias slice with PDK poly parasitics'}
    write_csv(Path(directory)/'spectrum.csv',[{'frequency_hz':float(a),'asd_v_sqrt_hz':float(b),'aperture_weight':float(c)} for a,b,c in zip(f,asd,weight)])
    result['metrics']=row;result['observations']=True;dump(Path(directory)/'metrics.json',row)
    compress_trace(Path(directory)/'output.raw');return result

def main():
    backend.BASE=OUT/'noise_backend';backend.BASE.mkdir(parents=True,exist_ok=True)
    cases=[]
    for j,(mos,n) in enumerate(itertools.product((False,True),NS)):
        deck,meta=noise_build(n,mos);cases.append((deck,meta,OUT/'spice/noise'/f'{j:03d}'))
    rows=backend.batch_execute(cases,2,processor)
    write_csv(OUT/'noise.csv',[r['metrics'] for r in rows]);print('Noise completed',len(rows),flush=True)
if __name__=='__main__':main()
