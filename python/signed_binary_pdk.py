"""Selected-bin elaboration: verbatim PDK model block, unmodified coefficients.
The original complete library still supplies all corner/global/mismatch params.
Discarding unused size bins avoids O(instances * 63) memory expansion. Native
wrapper replays are required and retained. Never substitutes a behavioral MOS.
"""
import re,hashlib
from functools import lru_cache
from pathlib import Path
from gamma_sky130 import PDK,short,dump
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/signed_binary_kcl_20261026/model_elaboration'
@lru_cache(None)
def selected(kind,corner,w,l):
    stem='sky130_fd_pr__'+('nfet' if kind=='n' else 'pfet')+'_01v8'
    candidates=list((PDK/'libs.ref/sky130_fd_pr/spice').glob(stem+'__'+corner+'.*spice'))
    for path in candidates:
        original=path.read_text()
        match=re.search(r'(?im)^\.subckt\s+'+stem+r'\s+d\s+g\s+s\s+b\s*$',original)
        if match:break
    else:raise ValueError('No original MOS definition '+str((kind,corner)))
    tail=original[match.start():];end=re.search(r'(?im)^\.ends\b',tail).start();body=tail[:end]
    models=list(re.finditer(r'(?im)^\.model\s+\S+\s+[np]mos\s*$',body))
    header=body[:models[0].start()];matches=[]
    for index,m in enumerate(models):
        block=body[m.start():models[index+1].start() if index+1<len(models) else len(body)]
        bounds={key:float(value) for key,value in re.findall(r'\b(lmin|lmax|wmin|wmax)\s*=\s*([-+\d.eE]+)',block)}
        if bounds['lmin']-1e-16<=l*1e-6<bounds['lmax']-1e-16 and bounds['wmin']-1e-16<=w*1e-6<bounds['wmax']-1e-16:matches.append((block,bounds))
    if len(matches)!=1:raise ValueError('Ambiguous/missing size bin '+str((kind,corner,w,l,len(matches))))
    block,bounds=matches[0];tag=stem+'_selected_'+hashlib.sha256(repr((corner,w,l)).encode()).hexdigest()[:12]
    header=re.sub(r'(?im)^(\.subckt\s+)'+stem+r'\b',lambda m:m[1]+tag,header,count=1)
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/(tag+'.spice')
    target.write_text('* Verbatim selected PDK model block; only outer wrapper renamed.\n'+header+block+'\n.ends '+tag+'\n')
    dump(target.with_suffix('.json'),{'original_file':str(path.relative_to(ROOT)),'original_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'selected_model_sha256':hashlib.sha256(block.encode()).hexdigest(),'bounds':bounds,'w_um':w,'l_um':l,'corner':corner,
        'wrapper':tag,'note':'Header parameters and native MOS invocation preserved; only subckt name changed; all model coefficients verbatim.'})
    return tag,short(target)

@lru_cache(None)
def focused(corner):
    lines=['* Focused original-PDK parameter includes; no altered coefficients.', '.option scale=1u', '.param mc_mm_switch=0 mc_pr_switch=0']
    ref=PDK/'libs.ref/sky130_fd_pr/spice';tech=PDK/'libs.tech/ngspice'
    for stem in ('sky130_fd_pr__nfet_01v8','sky130_fd_pr__pfet_01v8'):
        for path in ref.glob(stem+'__'+corner+'.*spice'):
            text=path.read_text();first=re.search(r'(?im)^\.subckt\b',text)
            prefix=text[:first.start()] if first else text
            prefix='\n'.join(line for line in prefix.splitlines() if not line.lower().startswith('.include'))
            lines.append(prefix)
        lines.append('.include "'+short(ref/(stem+'__mismatch.corner.spice'))+'"')
    for path in (tech/'parameters/lod.spice',tech/'r+c/res_typical__cap_typical.spice',tech/'r+c/res_typical__cap_typical__lin.spice',ref/'sky130_fd_pr__cap_mim_m3_1.model.spice',ref/'sky130_fd_pr__res_high_po_0p35.model.spice'):
        lines.append('.include "'+short(path)+'"')
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/('focused_'+corner+'.spice');target.write_text('\n'.join(lines)+'\n')
    return short(target)
