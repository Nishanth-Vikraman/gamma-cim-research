"""Audit and package the narrow SKY130 experiment without simulator reruns."""
import json,hashlib,shutil,zipfile
import numpy as np
from gamma_sky130 import ROOT,BASE,PDK,dump,fast_raw
from run_sky130_phase import OUT,write_csv


def main():
    s=json.loads((OUT/'sky130_phase_summary.json').read_text())
    subdocs={
      'smoke':'Four independent MOS/inverter/transmission-gate circuits at tt,ss,ff,sf,fs. All20 pass. Additional identical-device mismatch and temperature diagnostics validate the harness. See results/smoke_characterization.csv for currents and switching supply work.',
      'regen_cell':'boundary_deck in python/gamma_sky130.py is the reusable transistor writer generator. s000..s131 retain the sizing sweep; switch_probe0..5 retain Ron/leakage diagnostics. transition_matrix performs98 real writes to test all49 ordered transitions. The selected loaded-cell deck is ../pvt/penergy004/circuit.cir.',
      'reference_bank':'dc/circuit.cir is the physical seven-level/six-half-tap ladder test. Supply current and levels are in results/model_validation.json. Dynamic loading is simulated jointly in regeneration cases. The energy candidate uses4 Mohm total resistance and400 ns periods: static current0.45 uA and324 fJ/period nominally. Sharing/fanout is not demonstrated.',
      'comparator':'c000..c029 cover five corners and W=1/4/16 um, L=.15/.5 um. m000..m031 are seeded mismatch runs. The complementary StrongARM primitive tests the +0.5-digit boundary. Offset brackets are sampled; delay is clock-50%-referenced. Five other decisions and decoding remain ideal.',
      'loaded_mac':'l090..l101 retain hybrid S12-to-S4 direct/regenerated controls; extra gain3 and VCM loading make the direct branch fail even at zero injected residual. matched* removes that confound with common S4, unity-gain drivers, MIM storage and ideal latched decisions. At residual .35, direct gets1/12 correct versus12/12 regenerated; both pass at zero. standalone_audit independently replays the selected literal deck.',
      'pvt':'Two frozen candidates x5 corners x1.62/1.8/1.9 V x-40/25/100 C =90 cases. MIM storage and one real threshold are used. Four numerical nonconvergences retain partial waveforms and unresolved checkpoints. Wrong classes and terminal-limit failures remain in the tables. The suggested1.98 V exceeds the documented1.95 V core-model limit and is excluded.',
      'monte_carlo':'Three brackets x16 independent seeds. mp*: PDK local device/MIM mismatch only; ma*: explicit0.5% ideal-resistor,2% ideal-capacitor and0.5 mV systematic comparator-offset standard deviations with PDK mismatch disabled; mm*: larger-device PDK bracket. Corner tt; MC_PR_SWITCH=0. Correlated checkpoints are not independent yield samples.',
      'results':'Required tables, eight plots, sources, regression and replay audits. Costs separate measured VDD work, gross delivered energy, returned energy and external pin charging allowance. Missing decisions/latches/control are unpriced. max_settling_ns=1e9 is a censored0.1-digit target, not a measured one-second latency.'}
    for folder,text in subdocs.items():
        (BASE/folder/'README.md').write_text('# '+folder.replace('_',' ')+'\n\n'+text+'\n\nSee [engineering report](../README.md). Replay a retained batch case using replay.cir in its case directory. Regenerate decks after moving the project to refresh host-absolute model paths.\n',encoding='utf-8')
    levels='| Digit | Canonical V | Lower finite-region V | Upper V |\n| --- | --- | --- | --- |\n'+'\n'.join(f'| {d:+d} | {.6+d*.025:.4f} | {.6+(d-.5)*.025:.4f} | {.6+(d+.5)*.025:.4f} |' for d in range(-3,4))
    spec='''# q=4 canonical physical contract

Exact oracle: integer pairs in beta=2+sqrt(3), with4 beta^k=beta^(k-1)+beta^(k+1). This local experiment is position0, so expected pairs are(d,0). No floating beta establishes correctness.

Nominal VDD1.8 V, VCM0.6 V, output spacing25 mV. All positions use identical levels; indices carry significance. Supply-sweep references track VDD/3 and VDD/72 under the declared common-scale policy.

'''+levels+'''

Internal output thresholds:0.5375,0.5625,0.5875,0.6125,0.6375,0.6625 V. An upstream legal-digit contract is mandatory: six internal thresholds cannot flag arbitrary raw overflow. Final MAC instrumentation flags out-of-range classes and never repairs them offline.

The raw previous path has8.333333 mV spacing (S12); its +0.5 boundary is0.604166667 V, derived from actual tap6/tap7 through10/20 Mohm interpolation. Five other decisions remain behavioral. The writer selects one of seven actual ladder references at output S4, not a general DAC or expected-answer source. Known gain and bins scale together.

Each400 ns energy-candidate slot: reset16..40 ns, comparator evaluation44..392 ns, canonical write72..260 ns, observation360 ns. Prebiased warmup1 us. Actual stored voltage feeds a transistor OTA. R_x500 kohm, R_bias1 Mohm and C_MAC20 fF implement Vmac=(2Vdrive+Vbias)/3. The fixed original second operand is Vbias=VCM-d*25 mV. Scale-aware MAC spacing25/3 mV recovers2*d_x-d, expected d on this fixed-point workload.

Actual inverters, transmission gates, bias devices, ladder supply, storage and loaded driver are priced. Ideal one-hot decisions represent the declared initially latched class. Remaining decoder/latches/sequencer/completion are unbuilt. Exponents remain digital metadata. No q sweep, array, high-resolution ADC or positional voltage scaling occurs.

Every golden.json lists MOS W/L and estimated junction geometry. Only ordinary1.8 V NMOS/PMOS and MIM devices are used. Parallel MIM units <=1 pF realize large threshold decoupling. Ideal/MIM comparisons isolate nominal model effects. No layout extraction or thermal/sample-noise yield claim is made.
'''
    (BASE/'canonical_level_spec.md').write_text(spec,encoding='utf-8')
    (BASE/'REPRODUCE.md').write_text('''# Minimal reproduction

The installed PDK is unmodified. pdk_manifest.json records the URL, pinned commit, open_pdks provenance, SHA256 and licenses. No standard cells, I/O, SRAM or full ASIC toolchain is required.

Shallow-clone https://github.com/chrisklaus83188/sky130-ngspice-models to third_party/sky130-ngspice-models. The pinned commit is84a539fb024ab87ad9b46c9008c0deb298e250c0; fetch that commit with --depth 1 and detach-checkout if HEAD has moved. Enable Git core.longpaths on Windows before checkout if needed. Keep libs.tech/ngspice and libs.ref/sky130_fd_pr/spice.

The official Windows ngspice47 archive URL/hash is in the manifest. Put Spice64/bin/ngspice_con.exe under third_party/ngspice/. A minimal spinit selects one thread and avoids unused XSPICE plugins. The scripts use hidden console execution and native BSIM4; no Spectre patching or OSDI build is needed.

From python/: gamma_sky130.py runs20 smokes; validate_sky130.py validates mismatch/temperature; run_sky130_phase.py --workers 2 --mc-trials 16 runs the phase; audit_sky130.py independently replays the selected literal circuit; report_sky130_phase.py regenerates figures/report; package_sky130_phase.py builds this bundle. Python libraries are in requirements.txt. After moving the project, regenerate literal decks to refresh absolute model includes.

Every case retains literal netlist, seed, exact expected states, metrics and original log. Representative full traces and all numerical-failure partial traces remain locally as lossless gzip. The compact bundle contains all reported checkpoints/integrals/classes and replay inputs/logs, not bulky waveforms or downloadable PDK/runtime files. Old LTspice references are unchanged. The old energy workload summary needed for break-even is included. Numerical convergence failures never count as passes.
''',encoding='utf-8')
    smoke=[]
    for row in json.loads((OUT/'smoke_summary.json').read_text()):
        path=BASE/'smoke'/(row['corner']+'_'+row['kind']);names,x=fast_raw(path/'output.raw');ix={n:i for i,n in enumerate(names)};item=dict(row)
        if row['kind'] in ('inverter','transmission_gate'):
            t=x[:,0];i=-x[:,ix['i(vsmoke)']]
            item.update(quiescent_vdd_current_a=float(np.median(i[(t>2e-9)&(t<4e-9)])),active_static_vdd_current_a=float(np.median(i[(t>10e-9)&(t<14e-9)])),vdd_delivered_entire_test_fj=float(np.trapezoid(np.maximum(i,0),t)*1.8e15))
        smoke.append(item)
    write_csv(OUT/'smoke_characterization.csv',smoke)
    runs=sorted(set(r['run'] for r in s['sizing']+s['followup']+s['comparators']+s['comparator_mc']))
    missing=[];manifest=[]
    for run in runs:
        for name in ('circuit.cir','golden.json','replay.cir','simulator.log','metrics.json'):
            file=ROOT/run/name
            if not file.is_file():missing.append(str(file))
            else:manifest.append({'path':file.relative_to(ROOT).as_posix(),'bytes':file.stat().st_size,'sha256':hashlib.sha256(file.read_bytes()).hexdigest()})
    for name in ('sky130_phase_summary.json','energy_breakdown.csv','pvt_summary.csv','monte_carlo_summary.csv','sizing_sweep.csv','regen_transition_matrix.csv'):
        if not (OUT/name).is_file():missing.append(name)
    assert not missing,missing
    dump(OUT/'evidence_manifest.json',{'reported_physical_runs':len(runs),'cases':manifest,'missing':missing,'scope':'Reported sizing/PVT/MC/comparator/transitions/controls; smoke and literal-replay diagnostics retained separately.'})
    licenses=OUT/'licenses';licenses.mkdir(exist_ok=True)
    for name in ('LICENSE','NOTICE'):
        if (PDK/name).is_file():shutil.copyfile(PDK/name,licenses/('SKY130_'+name))
    bundle=OUT/'sky130_phase_bundle.zip';files=set()
    for folder in (ROOT/'docs',ROOT/'python'):
        files.update(f for f in folder.iterdir() if f.is_file() and f.suffix in ('.py','.md','.txt'))
    files.update(f for f in BASE.rglob('*') if f.is_file() and f.suffix.lower() in ('.cir','.json','.csv','.log','.md','.png','.txt','.py') and 'source_snapshot' not in f.parts and f.name!='package_audit.json')
    files.update(licenses.iterdir());files.update(ROOT/n for n in ('README.md','agent_handoff.md','roadmap.md','project_manifest.json'));files.add(ROOT/'spice/regen_results/phase_summary.json')
    with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for file in sorted(files):z.write(file,file.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(bundle) as z:assert z.testzip() is None
    result={'bundle':str(bundle),'bytes':bundle.stat().st_size,'sha256':hashlib.sha256(bundle.read_bytes()).hexdigest(),'reported_physical_runs':len(runs),'all_required_results_present':True,'zip_crc_verified':True}
    dump(OUT/'package_audit.json',result);print(result,flush=True)


if __name__=='__main__':main()
