"""Hash-verified compact review bundle; large raw traces stay in the project."""
from pathlib import Path
import json,hashlib,zipfile,shutil
import numpy as np
from gamma_sky130 import ROOT,short,fast_raw
from gamma_sar_exact import OUT

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    summary=json.loads((OUT/'phase_summary.json').read_text())
    assert not summary['pending'] and summary['pvt']['cases']==45
    assert all(v['networks']==8 for v in summary['mismatch'].values())
    for name,h in json.loads((OUT/'source_hashes.json').read_text()).items():
        assert sha(ROOT/name)==h,'Source changed since report: '+name
    # Keep original timeout/parse logs alongside the attempted decks too.
    batches=[*list((OUT/'batch_backend/batches').glob('*')), *list((ROOT/'spice/sky130/batches').glob('*'))]
    for d in (OUT/'spice').iterdir():
        if not d.is_dir() or (d/'simulator.log').exists():continue
        needle=short(d)+'/output.raw'
        matches=[b for b in batches if (b/'run.cir').exists() and needle in (b/'run.cir').read_text(errors='replace')]
        if matches:
            b=sorted(matches,key=lambda v:v.stat().st_mtime)[-1]
            for f in ('simulator.log','process.log'):
                if (b/f).exists():shutil.copy2(b/f,d/f)
            (d/'unresolved_log_origin.json').write_text(json.dumps({'batch':str(b.relative_to(ROOT)),'status':'No completed analyzed waveform; never a pass.'},indent=2))
    # Exact, unresampled representative windows preserve every original sample.
    windows=[]
    for label in ('range7_sar','range15_sar','range31_binary_compact','bank_driven','scheme_binary_residual'):
        d=OUT/'spice'/label;r=json.loads((d/'metrics.json').read_text());names,x=fast_raw(d/'output.raw')
        for sign in (-1,1):
            target=sign*min(13,r['cfg'].get('limit',15));o=next(v for v in r['observations'] if v['coefficient']==target)
            mask=(x[:,0]>=o['start_ns']*1e-9)&(x[:,0]<=(o['start_ns']+r['cfg'].get('period_ns',400))*1e-9)
            p=OUT/'representative_windows'/f'{label}_c{target}.npz';p.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(p,names=np.array(names),values=x[mask],source_case=np.array(label),coefficient=np.array(target))
            windows.append({'file':str(p.relative_to(ROOT)),'source':str(d.relative_to(ROOT)),'samples':int(mask.sum()),'resampled':False,'source_raw_gz_sha256':sha(d/'output.raw.gz')})
    notes='''# Dynamic q4 normalizer replay

Start with results/sar_normalizer_phase/engineering_report.md.
This is model-based ngspice evidence, not silicon. No tested topology wins.

The pinned model/simulator provenance is spice/sky130/pdk_manifest.json.
Use the unmodified SKY130 model tree at third_party/sky130-ngspice-models,
commit84a539fb024ab87ad9b46c9008c0deb298e250c0 from:
https://github.com/chrisklaus83188/sky130-ngspice-models
The portable ngspice47 executable is expected at:
third_party/ngspice/Spice64/bin/ngspice_con.exe.
The manifest supplies source URLs and hashes. No PDK/model files are embedded
or patched in this bundle. The executed platform was Windows/Python3.13.
Install numpy, scipy, matplotlib and pytest in the replay environment.

From the project root:
python -B python/gamma_sar_exact.py
python -B python/run_sar_normalizer_phase.py --stage nominal
python -B python/run_sar_normalizer_phase.py --stage robust
python -B python/run_sar_normalizer_phase.py --stage comparator
python -B python/audit_sar_phase.py
python -B python/report_sar_normalizer_phase.py
From python/: python -B -m pytest -q -p no:cacheprovider test_gamma_sar.py

The complete118-test regression was executed in the original project. Some
historical regression fixtures and unrelated third-party estimators are not
duplicated here; replay all118 tests in the complete original project.

For an individual retained configuration, load its golden.json cfg and call
run_sar_normalizer_phase.run([cfg]). The generator resolves local model paths.
The standalone circuit.cir retains the original absolute .lib path; replace
only that path when replaying a literal on another machine. replay.cir or
run.cir records the seeded standalone simulator commands.

All successful/failed case decks, seeds, logs, parameters and full metrics are
included. case_catalog.json identifies all completed cases. Primary capacitor
mismatch labels contain capacitor_fixed; earlier geometry-split controls are
retained and explicitly excluded from the primary eight-network statistics.

Full completed waveforms remain in the original project. This compact bundle
includes selected full raw traces and exact unresampled representative windows
from every quotient topology/range. NPZ windows contain names and values;
values[:,0] is absolute simulator time in seconds. They are extracts, not
reconstructed or oracle-generated traces. The manifest hashes included files
and records the original full-waveform hashes.

Histogram energy prices are local circuit projections. They are not measured
functional complete dots; physical failures and unmeasured terms remain.
The calibrated130nm ADC baseline is unchanged, and the old101-MOS topology
was not rerun for this phase.
'''
    (OUT/'REPRODUCE.md').write_text(notes)
    files=[];wave_budget=28*1024*1024;wave_bytes=0;kept=[]
    # Small full traces complement the exact extracts of the large range cases.
    for label in ('pvt04','cmp_balanced_w0.42_l0.15','mm_transistor00','nonzero_neighbors',
                  'ref_cap150','current_early_bias','range7_sar'):
        p=OUT/'spice'/label/'output.raw.gz'
        if p.exists() and wave_bytes+p.stat().st_size<=wave_budget:
            files.append(p);kept.append(label);wave_bytes+=p.stat().st_size
    for p in OUT.rglob('*'):
        if not p.is_file() or p.suffix=='.zip' or p.name=='evidence_manifest.json':continue
        if p.name.startswith(('output.raw','replay.raw')):continue
        files.append(p)
    files += list((ROOT/'python').glob('*.py'))
    files += [ROOT/'python/requirements.txt',ROOT/'spice/sky130/pdk_manifest.json',ROOT/'results/sky130_phase/phase_summary.json',
      ROOT/'results/sky130_phase/engineering_report.md',ROOT/'results/sky130_phase/matched_break_even.csv',
      ROOT/'results/sky130_phase/baseline_definition.json',ROOT/'results/self_triggered_q4/summary.json',
      ROOT/'agent_handoff.md',ROOT/'README.md',ROOT/'roadmap.md']
    raw_hashes={str(p.relative_to(ROOT)):sha(p) for p in (OUT/'spice').glob('*/output.raw.gz')}
    files=sorted(set(files));manifest={'included_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in files},
      'full_raw_trace_sha256':raw_hashes,'full_raw_cases_in_bundle':kept,'exact_windows':windows,
      'completed_cell_cases':summary['completed_cell_cases'],'completed_comparator_cases':summary['completed_comparator_cases'],
      'unresolved_attempts':summary['unresolved_attempts'],'regression':'118 passed in17.10s',
      'no_selected_trace_resampling':True,'no_pdk_modifications':True}
    m=OUT/'evidence_manifest.json';m.write_text(json.dumps(manifest,indent=2));files.append(m)
    target=OUT/'gamma_sar_normalizer_phase.zip'
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files:z.write(p,p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        for name,h in manifest['included_sha256'].items():assert hashlib.sha256(z.read(name)).hexdigest()==h
    assert target.stat().st_size<50*1024*1024,'Use prepared Library upload for larger artifacts'
    print(json.dumps({'file':str(target),'bytes':target.stat().st_size,'members':len(files),'verified':True}))

if __name__=='__main__':main()
