"""Compact review bundle. All literal decks/metrics; representative raw traces."""
from pathlib import Path
import hashlib,json,zipfile
from gamma_sky130 import ROOT
OUT=ROOT/'results/self_triggered_q4'

def main():
    keep={'selected_s1','selected_s-1','fine_s1','fine_s-1','pvt_1_s1','pvt_1_s-1',
          'pvt_2_s1','pvt_2_s-1','mc_2001_s1','mc_2001_s-1','hold_s1_u3.2','hold_s1_u4.8',
          'hold_s-1_u3.2','hold_s-1_u4.8','repeat_s1','repeat_s-1','idle'}
    files=[]
    for p in OUT.rglob('*'):
        if not p.is_file() or p.suffix=='.zip' or p.name in ('evidence_manifest.json',):continue
        if p.name.startswith('output.raw') and p.parent.name not in keep:continue
        files.append(p)
    files+=list((ROOT/'python').glob('*.py'))
    files += [ROOT/'python/requirements.txt',ROOT/'spice/sky130/pdk_manifest.json',
      ROOT/'results/sky130_phase/phase_summary.json',ROOT/'results/sky130_phase/engineering_report.md',
      ROOT/'results/sky130_phase/matched_break_even.csv',ROOT/'agent_handoff.md',ROOT/'README.md',ROOT/'roadmap.md']
    notes='''# Replay self-triggered q4 experiment

Read results/self_triggered_q4/engineering_report.md first.

This is a Windows/ngspice 47 experiment, not silicon measurement. The selected
cell is rejected on energy, physical-state accuracy and PVT/mismatch robustness.

Create the unmodified model tree at third_party/sky130-ngspice-models using:
https://github.com/chrisklaus83188/sky130-ngspice-models
commit 84a539fb024ab87ad9b46c9008c0deb298e250c0.
On Windows enable git core.longpaths for that checkout.

Install the pinned portable simulator at
third_party/ngspice/Spice64/bin/ngspice_con.exe.
Archive/source hashes and download URLs: spice/sky130/pdk_manifest.json.
Python dependencies: numpy, matplotlib, pytest (the executed interpreter was 3.13).

From this extracted project root:
python -B python/sky130_self_topple.py
python -B python/sky130_self_topple.py --characterize
python -B python/report_self_topple.py
From python/: python -B -m pytest -q -p no:cacheprovider

The Python generator resolves model paths on the current machine. Existing
literal circuit.cir files retain the original absolute .lib path; replace only
that path to independently replay a literal circuit on another machine.
run.cir or replay.cir is the per-case seeded simulator harness.

All attempted decks, logs and metrics are included. Full waveforms for selected,
fine-step and representative failed robustness/threshold/repeat cases are included.
Other completed raw traces remain in the original project and can be regenerated.
Model parsing rejections are retained as unresolved, never successful simulations.
The JSON reports include all 43 fixed characterization cases plus two selected
nominal runs. Historical tuning is not part of yield evidence.
'''
    (OUT/'REPRODUCE.md').write_text(notes);files.append(OUT/'REPRODUCE.md')
    manifest={'sha256':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(files))},
              'representative_waveform_cases':sorted(keep),'all_completed_waveforms_retained_locally':True,
              'regression':'106 passed in 17.18s'}
    (OUT/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2));files.append(OUT/'evidence_manifest.json')
    target=OUT/'gamma_self_triggered_q4.zip'
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(set(files)):z.write(p,p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        for name,h in manifest['sha256'].items():assert hashlib.sha256(z.read(name)).hexdigest()==h
    print(json.dumps({'file':str(target),'bytes':target.stat().st_size,'members':len(files),'verified':True}))

if __name__=='__main__':main()
