"""Pinned standalone CiMLoop-lineage ADC calibration. Upstream code is unmodified.

Energy from the fitted upstream routine is pJ, despite ADCRequest's stale J
docstring. The Accelergy wrapper, training CSV and successor SI conversion
establish the units. Interleaving changes frequency, not conversions charged.
"""
from __future__ import annotations
import csv,hashlib,importlib.util,importlib.metadata,json,logging,math,platform,subprocess,sys,types
from pathlib import Path
from contextlib import contextmanager
import numpy as np
import pandas as pd
import yaml

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'baseline_adc'
OUT=ROOT/'results/sky130_phase'
LEGACY=ROOT/'third_party/accelergy-adc-plug-in'
SUCCESSOR=ROOT/'third_party/hwcomponents-adc/hwcomponents_adc'
PINNED={
    'accelergy-adc-plug-in':('https://github.com/Accelergy-Project/accelergy-adc-plug-in','ff0fb48a4e55d562c3de8362fd0c74418fd27d85'),
    'hwcomponents-adc':('https://github.com/Accelergy-Project/hwcomponents-adc','22cf87b20e6696c1e3d49118fb423bbb827e97d8'),
    'cimloop':('https://github.com/mit-emze/cimloop','481c044e96034e125c1a442e92c73b5d05c1f506')}
SURVEY_COMMIT='79a4b6b1cfd0a337371804b0898cfba2e8a495ac'
THROUGHPUTS=(1e6,2.5e6,5e6,1e7,4.25e7,8.5e7,1e8,1.6e8,3.2e8,1e9,1e10,2e10)
COUNTS=(1,2,4,8,16,32,64)


def dump(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')


def write_csv(path,rows):
    rows=list(rows);path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module;spec.loader.exec_module(module);return module


@contextmanager
def isolated_modules(names):
    old={k:sys.modules.get(k) for k in names}
    try:yield
    finally:
        for k,v in old.items():
            if v is None:sys.modules.pop(k,None)
            else:sys.modules[k]=v


class ADCModel:
    def __init__(self,successor=False):
        directory=SUCCESSOR if successor else LEGACY
        names=('hwcomponents_adc','hwcomponents_adc.headers','hwcomponents_adc.model','hwcomponents_adc.optimizer') if successor else ('headers','model','optimizer')
        with isolated_modules(names):
            if successor:
                # Import the standalone numerical modules, without importing
                # __init__/main and their optional HWComponents framework.
                package=types.ModuleType('hwcomponents_adc');package.__path__=[str(directory)]
                sys.modules['hwcomponents_adc']=package;prefix='hwcomponents_adc.'
            else:prefix=''
            self.h=load_module(prefix+'headers',directory/'headers.py')
            self.model=load_module(prefix+'model',directory/'model.py')
            self.optimizer=load_module(prefix+'optimizer',directory/'optimizer.py')
        self.fit=yaml.safe_load((directory/'adc_data/model.yaml').read_text())
        self.logger=logging.getLogger('gamma_adc_standalone');self.logger.addHandler(logging.NullHandler());self.logger.propagate=False
        if BASE.exists() and not any(isinstance(h,logging.FileHandler) for h in self.logger.handlers):
            self.logger.addHandler(logging.FileHandler(BASE/'model_warnings.log',encoding='utf-8'))
        self.model.logger=self.logger
        self.max_frequency_hz=math.exp(self.fit[self.h.FREQ][self.h.MAX])

    def estimate(self,bits,technology_nm,aggregate_samples_s,n_adcs,allow_extrapolation=False):
        if n_adcs<1 or int(n_adcs)!=n_adcs or aggregate_samples_s<=0 or technology_nm<=0:
            raise ValueError('Positive technology/throughput and positive integer ADC count required')
        frequency=aggregate_samples_s/n_adcs
        params={self.h.ENOB:bits,self.h.TECH:math.log(technology_nm),self.h.FREQ:math.log(frequency)}
        direct=self.model.get_energy(params,self.fit,allow_extrapolation)
        request=self.optimizer.ADCRequest(bits,technology_nm,aggregate_samples_s,n_adcs,logger=self.logger)
        energy=request.energy_per_op(self.fit);area=request.area(self.fit)
        assert math.isclose(direct,energy,rel_tol=1e-13)
        return {'enob':bits,'technology_nm':technology_nm,'aggregate_samples_s':aggregate_samples_s,
            'n_adcs':int(n_adcs),'frequency_per_adc_hz':frequency,'energy_per_conversion_pj':energy,
            'energy_per_conversion_fj':energy*1000,'area_per_adc_um2':area/n_adcs,'total_area_um2':area,
            'bank_power_mw':energy*aggregate_samples_s*1e-9,
            'model_max_frequency_hz':self.max_frequency_hz,'within_model_hard_frequency_limit':frequency<=self.max_frequency_hz*(1+1e-14)}


def provenance():
    sources=[]
    for name,(url,commit) in PINNED.items():
        directory=ROOT/'third_party'/name
        actual=subprocess.check_output(['git','rev-parse','HEAD'],cwd=directory,text=True).strip()
        assert actual==commit,(name,actual,commit)
        dirty=subprocess.check_output(['git','status','--porcelain'],cwd=directory,text=True).strip()
        assert not dirty,(name,dirty)
        source={'name':name,'repository_url':url,'commit':actual,'unmodified':True,'license':'MIT'}
        files=list(directory.rglob('*.py'))+list(directory.glob('*.md'))+list(directory.glob('*.txt'))
        if name!='cimloop':files+=list(directory.rglob('model.yaml'))+list(directory.rglob('adc_list.csv'))
        source['file_hashes']={f.relative_to(directory).as_posix():sha(f) for f in files if '.git' not in f.parts}
        sources.append(source)
    survey=BASE/'sources/ADCsurvey_79a4b6b.xls'
    return {'date':'2026-10-24','repositories':sources,'published_method':'https://arxiv.org/abs/2404.06553v2',
        'historical_model_version':'0.1 (setup.py)','survey_gitlink_commit':SURVEY_COMMIT,
        'survey_url':f'https://raw.githubusercontent.com/bmurmann/ADC-survey/{SURVEY_COMMIT}/xls/ADCsurvey_latest.xls',
        'survey_sha256':sha(survey),'python':platform.python_version(),
        'dependencies':{p:importlib.metadata.version(p) for p in ('numpy','pandas','PyYAML','xlrd')},
        'installation':'Standalone numerical modules only; no CiMLoop/Accelergy/HWComponents stack installed.',
        'energy_unit_evidence':'adc_data/adc_list.csv energy(pJ/op), Accelergy Estimation(...,"p"), successor main.get_energy multiplies numerical model output by1e-12.',
        'docstring_defect':'ADCRequest.energy_per_op says Joules in both lineages, but its fitted numerical output is pJ. No upstream source patched.',
        'successor_note':'Same model coefficients; successor wraps pJ in SI and splits20% assumed leakage from80% action energy. Use total100% for matched steady-throughput cost; do not take an extra20% discount.',
        'training_dataset_note':'Distributed215-row adc_list.csv is the fitted model input. The pinned survey has more records; only matched distributed rows are claimed as training anchors. No model refit.'}


def survey_anchors():
    training=pd.read_csv(LEGACY/'adc_data/adc_list.csv')
    frames=[]
    for venue in ('ISSCC','VLSI'):
        d=pd.read_excel(BASE/'sources/ADCsurvey_79a4b6b.xls',sheet_name=venue);d['venue']=venue;frames.append(d)
    raw=pd.concat(frames,ignore_index=True)
    raw=raw[pd.to_numeric(raw['TECHNOLOGY'],errors='coerce')==.13]
    def number(v):
        try:r=float(v);return r if math.isfinite(r) else None
        except (TypeError,ValueError):return None
    rows=[]
    for i,r in raw.iterrows():
        fs,nyq,power,sndr=(number(r[c]) for c in ('fs [Hz]','fsnyq [Hz]','P [W]','SNDR_plot [dB]'))
        area=number(r['AREA [mm^2]']);energy=number(r['P/fsnyq [pJ]'])
        bits=(sndr-10*math.log10(1.5))/(20*math.log10(2)) if sndr is not None else None
        matches=training[(training['tech node (nm)']==130)&np.isclose(training['frequency (Hz)'],fs or -1,rtol=1e-8)&np.isclose(training['SNDR_plot [dB]'],sndr or -1,rtol=1e-8)&np.isclose(training['energy (pJ/op)'],energy or -1,rtol=1e-8)]
        row={'anchor_id':f"{r['venue']}_{r['YEAR']}_{r['ID']}",'venue':r['venue'],'year':number(r['YEAR']),'paper_id':str(r['ID']),
            'technology_nm':130,'architecture':str(r['ARCHITECTURE']),'enob_from_sndr':bits,'sndr_db':sndr,
            'sampling_clock_hz':fs,'nyquist_rate_hz':nyq,'osr':number(r['OSR']),
            'independent_value_rate_bound_hz':min(fs,nyq) if fs and nyq else None,
            'energy_per_sample_pj':power/fs*1e12 if power and fs else None,'survey_energy_per_nyquist_pj':energy,
            'area_um2':area*1e6 if area else None,'power_w':power,'title':str(r['TITLE']),
            'doi':str(r['DOI']),'link':str(r['LINK']),'in_distributed_training_csv':not matches.empty,
            'training_csv_row':int(matches.index[0]) if not matches.empty else None,
            'survey_commit':SURVEY_COMMIT,'rate_warning':'Do not equate internal oversampling clock to independent output conversions' if fs and nyq and not math.isclose(fs,nyq) else ''}
        rows.append(row)
    write_csv(BASE/'survey_all_130nm.csv',rows)
    selected=[r for r in rows if r['in_distributed_training_csv']]
    assert len(selected)==sum(training['tech node (nm)']==130)==12
    write_csv(OUT/'adc_130nm_survey_anchors.csv',selected)
    dump(BASE/'survey_anchor_notes.json',{'training_records':len(training),'training_130nm_records':len(selected),'all_pinned_survey_130nm_records':len(rows),
        'oversampling_note':'The historical fit uses fs as frequency and P/fsnyq as energy. Retain both sample/Nyquist rates and energies; a0.1 pJ internal modulator clock is not a0.1 pJ independent CIM-result conversion.',
        'anchor_envelope_rule':'Optional optimistic empirical competitor uses energy at its published rate, ENOB>=target, available value rate>=per-ADC request, OSR<=1, and known area. It is not claimed as a CIM-optimized peripheral or guaranteed clock-gatable implementation.'})
    return training,selected,rows


def reproduction(adc):
    successor=ADCModel(True)
    example=adc.estimate(8,16,512e7,32)
    same=successor.estimate(8,16,512e7,32)
    assert adc.fit==successor.fit
    assert example==same
    # Independent hand-evaluation of the pinned coefficients, without calling
    # the upstream prediction routine. Never refit or tune the model.
    f=math.log(512e7/32);m=adc.fit;h=adc.h;coeff=m[h.FOMS]
    foms=min(coeff[h.INTERCEPT]+coeff[h.FREQ]*f,coeff[h.MAX_BY_ENOB][8])
    sndr=8*20*math.log10(2)+10*math.log10(1.5)
    manual=.5*10**((sndr-foms)/10)*math.exp(coeff[h.TECH_INTERCEPT]+coeff[h.TECH_SLOPE]*math.log(16)+coeff[h.ENOB_SLOPE]*math.log(8)+coeff[h.ENRG_RESIDUAL])
    assert math.isclose(manual,example['energy_per_conversion_pj'],rel_tol=1e-12)
    rejected=False
    try:adc.estimate(8,130,2e10,1)
    except AssertionError:rejected=True
    assert rejected
    rows=[]
    for bits in range(4,11):
        for node in (16,32,130):
            rows.append(adc.estimate(bits,node,1.6e8,1))
    write_csv(BASE/'reproduction.csv',rows)
    result={'published_example_source':'accelergywrapper.py __main__:8bits,16nm,512e7 aggregate samples/s,32ADCs',
        'example':example,'independent_closed_form_pj':manual,'successor_numerical_result_identical':True,
        'coefficients_identical_across_lineages':True,'rejects_20ghz_with_extrapolation_disabled':rejected,
        'expected_published_numeric_golden':'No fixed numeric answer is supplied by upstream; reproduction executes its exact example arguments and checks independent evaluation.',
        'all_four_inputs_confirmed_from_code':True,'energy_is_per_conversion_not_divided_by_n_adcs':True,'example_energy_j':manual*1e-12}
    dump(BASE/'reproduction.json',result);return result


def pareto_table(adc,training):
    anchors=training[training['tech node (nm)']==130];rows=[]
    for bits in range(4,11):
        local=anchors[abs(anchors['number of bits']-bits)<=.5]
        for rate in THROUGHPUTS:
            group=[]
            for n in COUNTS:
                row=adc.estimate(bits,130,rate,n,True);warnings=[]
                if not row['within_model_hard_frequency_limit']:warnings.append('REJECTED by strict upstream10GHz per-ADC cap; extrapolated value retained only diagnostically')
                if local.empty:warnings.append('no130nm distributed training point within0.5 ENOB')
                elif not local['frequency (Hz)'].min()<=rate/n<=local['frequency (Hz)'].max():warnings.append('outside130nm near-resolution training frequency interval')
                if rate/n<anchors['frequency (Hz)'].min():warnings.append('below minimum130nm training clock; low-rate leakage/gating not modeled')
                if rate/n>anchors['frequency (Hz)'].max():warnings.append('above maximum130nm training clock')
                row.update(extrapolation_warning='; '.join(warnings),near_resolution_130nm_anchor_count=len(local),
                    model_validity_note='Code limits ENOB>=4 and per-ADC rate<=10GHz. No node-specific physical feasibility or confidence interval is enforced. Mixed Nyquist/oversampling survey fit.')
                group.append(row)
            for r in group:
                r['on_energy_area_frontier']=r['within_model_hard_frequency_limit'] and not any(
                    t['within_model_hard_frequency_limit'] and t['energy_per_conversion_pj']<=r['energy_per_conversion_pj'] and t['total_area_um2']<=r['total_area_um2'] and
                    (t['energy_per_conversion_pj']<r['energy_per_conversion_pj'] or t['total_area_um2']<r['total_area_um2']) for t in group)
            rows.extend(group)
    write_csv(OUT/'adc_130nm_pareto.csv',rows);write_csv(OUT/'adc_130nm_frontier.csv',[r for r in rows if r['on_energy_area_frontier']]);return rows


def main():
    BASE.mkdir(exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    manifest=provenance();dump(BASE/'source_manifest.json',manifest)
    adc=ADCModel();repro=reproduction(adc);training,anchors,allanchors=survey_anchors();rows=pareto_table(adc,training)
    manifest.update(model_frequency_limit_hz=adc.max_frequency_hz,grid_rows=len(rows),training_130nm_anchors=len(anchors),all_raw_130nm_anchors=len(allanchors),
        reproduction_passed=True,coefficients_sha256=sha(LEGACY/'adc_data/model.yaml'))
    dump(OUT/'baseline_manifest.json',manifest)
    (BASE/'README.md').write_text('''# Pinned130nm ADC baseline

The unmodified historical ADCRequest/model modules execute standalone; no CiMLoop, Accelergy or HWComponents framework is installed. See source_manifest.json for pinned repositories, model/training hashes and the exact survey commit. No coefficients were refitted.

Reproduce from python/: python -B gamma_adc_baseline.py. It executes the example in the upstream Accelergy wrapper, independently checks the formula, verifies identical successor coefficients/results and validates the hard10GHz rejection. The four parameters are resolution/ENOB, technology in nm, aggregate conversions/s and ADC count. Per-ADC frequency is aggregate/count; interleaving never divides energy/conversion. Total area is per-ADC area times count.

The fitted numerical energy is pJ, despite a stale Joules docstring in ADCRequest. The training CSV and Accelergy p-unit wrapper agree, and successor main.py explicitly converts the value by1e-12. The successor's80% action/20% leakage split does not reduce the total steady-throughput cost.

The full grid and within-(ENOB,throughput) energy/area frontier are under results/sky130_phase/. Extrapolated/refused points remain flagged. The10GHz cap is a hard-coded model limit, not a demonstrated130nm feasibility boundary. Near-resolution130nm training support is sparse. Modeled throughput/area does not include an ADC multiplexer or idle-gating implementation.

Twelve130nm records are in the distributed215-row training CSV; all51 raw130nm survey records are exported separately. The fit uses fs as frequency and P/fsnyq as energy: raw sample clock, independent Nyquist value rate and both energy definitions are retained. Low-power oversampling clock energy must not be misreported as independent full-result conversion energy. Some measured anchors beat the fitted trend; it is an optimistic statistical model, not a strict physical lower bound for every ADC.

The original100fJ ADC and2.74fJ regeneration figures are historical sensitivity assumptions only. The new comparison is in results/sky130_phase/engineering_report.md; old reports are preserved.
''',encoding='utf-8')
    print(json.dumps({'example':repro['example'],'grid_rows':len(rows),'training_130nm_anchors':len(anchors),'baseline_ready_for_phase_B':True},indent=2),flush=True)


if __name__=='__main__':main()
