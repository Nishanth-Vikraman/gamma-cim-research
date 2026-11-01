"""Validate, preserve and package new matched-node evidence without old reruns."""
from pathlib import Path
import hashlib,json,zipfile,shutil
from gamma_adc_baseline import ROOT,BASE,OUT,PINNED,dump,sha

REQUIRED=('adc_130nm_pareto.csv','adc_130nm_survey_anchors.csv','adc_cim_amortization.csv','baseline_manifest.json',
    'sky130_device_manifest.json','sky130_smoke_summary.json','regen_energy.csv','regen_transition_matrix.csv',
    'loaded_chain_results.csv','pvt_results.csv','mismatch_results.csv','matched_break_even.csv','phase_summary.json','engineering_report.md')


def historical_preservation():
    archive=ROOT/'spice/sky130/results/sky130_phase_bundle.zip';checked=[];changed=[]
    # The previously published archive independently fixes these report/result
    # bytes. Exclude packaging operational logs written after its creation.
    with zipfile.ZipFile(archive) as z:
        for name in z.namelist():
            historical=(name.startswith('spice/sky130/') or name.startswith('spice/regen_results/') or name=='docs/sky130_phase_results.md')
            if not historical or name.endswith('package_sky130_phase.py.log') or 'source_snapshot/' in name:continue
            file=ROOT/name
            if not file.exists():changed.append({'path':name,'reason':'missing'});continue
            expected=hashlib.sha256(z.read(name)).hexdigest();actual=sha(file)
            if expected!=actual:changed.append({'path':name,'expected':expected,'actual':actual})
            checked.append({'path':name,'sha256':actual})
    result={'reference_archive_sha256':sha(archive),'checked_historical_files':len(checked),'changed':changed,
        'all_covered_historical_reports_results_decks_logs_unchanged':not changed,'files':checked}
    dump(OUT/'historical_preservation.json',result);assert not changed,changed
    return result


def main():
    assert all((OUT/name).is_file() for name in REQUIRED)
    preservation=historical_preservation()
    sources=[ROOT/'python'/n for n in ('gamma_adc_baseline.py','gamma_matched_130nm.py','sky130_matched_cap_smoke.py',
        'sky130_reference_load.py','report_matched_130nm.py','package_matched_130nm.py','test_matched_130nm.py')]
    snapshots=OUT/'source_snapshot';snapshots.mkdir(exist_ok=True)
    for file in sources:shutil.copyfile(file,snapshots/file.name)
    dump(OUT/'source_hashes.json',{p.relative_to(ROOT).as_posix():sha(p) for p in sources})
    instructions='''# Matched130nm phase bundle

Start at results/sky130_phase/engineering_report.md and phase_summary.json. The report answers restoration, physical cost, matched energy, dominant overhead and the next experiment separately. Historical reports/results are preserved; historical_preservation.json compares their bytes against the previously published archive.

Standalone ADC calibration: python -B python/gamma_adc_baseline.py. Matched accounting: python -B python/gamma_matched_130nm.py. New capacitor/reference diagnostics are optional physical reruns: sky130_matched_cap_smoke.py and sky130_reference_load.py. Then run report_matched_130nm.py and package_matched_130nm.py from python/. See baseline_adc/requirements.txt.

Recreate the Git checkouts at the commits in baseline_adc/source_manifest.json; sources are also bundled for inspection. The estimator is numerical Python only: do not install full CiMLoop/Accelergy/HWComponents. SKY130/ngspice runtime/models remain fetched from the pinned source/device manifests; their large binaries and PDK trees are not duplicated here. Existing legacy Python exact-reference modules and compact circuit evidence are included. Standalone netlist .lib paths are refreshed by regeneration after moving the project.

All new physical seeds, decks, logs, metrics and lossless waveforms are included. Existing351 physical runs retain their compact replay/golden/log records; large older waveforms remain in the original project. The original LTspice result directories remain unchanged. This is a conditional matched-energy study, not a measured full-detector/normalizer or physical accuracy win.
'''
    (OUT/'REPRODUCE.md').write_text(instructions,encoding='utf-8')
    files=set()
    # New phase, data and unmodified estimator source; no Git internals/runtime.
    for folder in (OUT,BASE):
        files.update(p for p in folder.rglob('*') if p.is_file() and p.suffix!='.zip' and p.name not in ('package_audit.json','phase_evidence_manifest.json'))
    for name in PINNED:
        directory=ROOT/'third_party'/name
        if name=='cimloop':
            selected=[directory/'README.md',directory/'LICENSE']+list((directory/'workspace/models/arch/1_macro/raella_isca_2023').glob('*.yaml'))
        else:
            selected=[p for p in directory.rglob('*') if p.is_file() and '.git' not in p.parts and '__pycache__' not in p.parts and p.suffix.lower() in ('.py','.yaml','.csv','.txt','.md','.toml')]
        files.update(p for p in selected if p.is_file())
    files.update(p for p in (ROOT/'python').iterdir() if p.is_file() and p.suffix in ('.py','.md','.txt'))
    # Preserve all compact original circuit evidence without regenerating it.
    files.update(p for p in (ROOT/'spice/sky130').rglob('*') if p.is_file() and p.suffix in ('.cir','.json','.csv','.log','.md','.png','.txt') and 'source_snapshot' not in p.parts)
    for name in ('phase_summary.json','completion_traces.json','finite_word_golden.json','supplemental_checks.json','energy_sweep.json','summary.csv'):
        p=ROOT/'spice/regen_results'/name
        if p.is_file():files.add(p)
    files.update(ROOT/name for name in ('README.md','agent_handoff.md','roadmap.md','project_manifest.json'))
    files.add(ROOT/'docs/matched_130nm_phase.md')
    entries=[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]
    dump(OUT/'phase_evidence_manifest.json',{'required_outputs':list(REQUIRED),'file_count':len(entries),
        'historical_files_checked':preservation['checked_historical_files'],'files':entries})
    files.add(OUT/'phase_evidence_manifest.json')
    archive=OUT/'matched_130nm_phase_bundle.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for file in sorted(files):z.write(file,file.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    result={'bundle':str(archive),'bytes':archive.stat().st_size,'sha256':sha(archive),'file_count':len(files),
        'required_outputs_present':True,'historical_files_unchanged':True,'zip_crc_verified':True}
    dump(OUT/'package_audit.json',result);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
