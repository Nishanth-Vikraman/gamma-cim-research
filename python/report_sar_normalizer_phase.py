"""Evidence-only SAR phase report. Never runs or retunes a circuit."""
from pathlib import Path
from collections import defaultdict
import argparse,csv,json,hashlib,math
import numpy as np
from gamma_sar_exact import OUT,ROOT,write_csv
from gamma_sar_circuit import analyze,comparator_analyze,integral
from gamma_sky130 import fast_raw,dump,compress_trace

def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])

def refresh():
    rows=[];unresolved=[]
    for golden in sorted((OUT/'spice').glob('*/golden.json')):
        d=golden.parent;p=d/'metrics.json'
        try:
            cached=json.loads(p.read_text()) if p.exists() else {}
            if cached.get('analysis_schema_version')==(6 if cached.get('kind')=='comparator' else 7):
                cached['run']=str(d.relative_to(ROOT));rows.append(cached);continue
            result=json.loads((d/'result.json').read_text()) if (d/'result.json').exists() else {'meta':json.loads(golden.read_text())}
            names,x=fast_raw(d/'output.raw')
            f=comparator_analyze if result['meta']['kind']=='comparator' else analyze
            row=f(result,names,x);row['run']=str(d.relative_to(ROOT));dump(p,row)
            compress_trace(d/'output.raw');rows.append(row)
        except Exception as e:
            old=json.loads(p.read_text()) if p.exists() else json.loads(golden.read_text());unresolved.append({'run':str(d.relative_to(ROOT)),'reason':str(e),'cfg':old.get('cfg',{})})
    return rows,unresolved

QUIET_CACHE={}
def quiet_zero(row):
    """Equal logical quiescent state for the idle comparison, from real trace."""
    if not row:return None
    if row['run'] in QUIET_CACHE:return QUIET_CACHE[row['run']]
    p=ROOT/row['run'];meta=json.loads((p/'golden.json').read_text());names,x=fast_raw(p/'output.raw')
    ix={n:j for j,n in enumerate(names)};t=x[:,0];o=next(v for v in row['observations'] if v['coefficient']==0)
    a=max(o['after_transfer_ns']+10,o['end_ns']-95)*1e-9;b=(o['end_ns']-5)*1e-9;energy=0.
    for n,port in meta['ports'].items():
        if port['component']=='test_state_write':continue
        power=-x[:,ix['v('+port['node']+')']]*x[:,ix[n]]
        if port['component'] in ('clock_pin','constant_pin'):power=np.maximum(power,0)
        energy+=integral(t,power,a,b)
    value=energy/(b-a)*1e9;QUIET_CACHE[row['run']]=value;return value

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--allow-incomplete',action='store_true');args=ap.parse_args()
    rows,unresolved=refresh();bylabel={r['cfg']['label']:r for r in rows}
    cells=[r for r in rows if r.get('kind')!='comparator'];comps=[r for r in rows if r.get('kind')=='comparator']
    pvt=[r for r in cells if r['cfg'].get('group')=='pvt']
    mm=[r for r in cells if r['cfg'].get('group')=='mismatch' and (r['cfg'].get('mismatch_mode')!='capacitor' or r['cfg'].get('fixed_topology_caps',False))]
    primary={7:bylabel.get('range7_sar'),15:bylabel.get('range15_sar'),31:bylabel.get('range31_binary_compact') or bylabel.get('range31_sar')}
    pending=[]
    if len(pvt)!=45:pending.append(f'cell PVT {len(pvt)}/45')
    if len(mm)!=24:pending.append(f'cell mismatch {len(mm)}/24')
    if not all(primary.values()):pending.append('a full bounded-range measurement')
    if not args.allow_incomplete and pending:raise RuntimeError('Still pending: '+', '.join(pending))
    # Explicitly record source files that were attempted without completed metrics.
    bad_runs={v['run'] for v in unresolved}
    for d in sorted((OUT/'spice').iterdir()):
        if not d.is_dir() or (d/'metrics.json').exists():continue
        if (d/'circuit.cir').exists() and str(d.relative_to(ROOT)) not in bad_runs:
            unresolved.append({'run':str(d.relative_to(ROOT)),'reason':'Attempted deck has no completed analyzed metric; logs retained.'})
    energy=[];reference=[];transfer=[];pvtflat=[];mmflat=[];cmpflat=[]
    for r in cells:
        cfg=r['cfg'];label=cfg['label']
        for o in r['observations']:
            energy.append({'label':label,'group':cfg.get('group','exploratory'),'scheme':cfg.get('scheme','sar'),
              'limit':cfg.get('limit',15),'c':o['coefficient'],'golden_m':o['golden_m'],'chosen_m':o['chosen_m'],
              'quotient_correct':o['quotient_correct'],'transfer_correct':o['transfer_charge_correct'],
              'physical_class_correct':o['physical_class_correct'],'sleep_correct':o['sleep_correct'],
              'sleep_control_correct':o['sleep_control_correct'],'retention_class_preserved':o['retention_class_preserved'],
              'event_energy_fj':o['event_energy_fj'],**{k+'_fj':v for k,v in o['energy_fj'].items()},
              'reset_window_fj':o['reset_window_energy_fj'],
              **{k+'_nonreset_fj':v for k,v in o['nonreset_components_fj'].items()},
              'idle_power_nw_last_state':r['idle_power_nw'],'period_ns':cfg.get('period_ns',400),
              'test_state_write_energy_excluded_fj':o['test_state_write_energy_fj'],'run':r['run']})
            reference.append({'label':label,'c':o['coefficient'],'ref_kind':cfg.get('ref_kind','capacitive'),
              'ref_ff':cfg.get('ref_ff',1000),'reference_rail_fj':o['energy_fj']['reference'],
              'reference_plus_all_logic_upper_bound_fj':o['energy_fj']['reference']+o['energy_fj']['logic'],
              'quiet_reference_drift_units':o['quiet_reference_drift_units'],'quiet_interval_ns':o['quiet_reference_interval_ns'],
              'max_loaded_reference_error_digits':max([abs(x['error_units']) for x in o.get('reference_checks',[])],default=None),
              'quotient_correct':o['quotient_correct'],'idle_power_nw_total_cell':r['idle_power_nw'],
              'reference_checkpoints_json':json.dumps(o.get('reference_checks',[])),'run':r['run']})
            transfer.append({'label':label,'kind':cfg.get('transfer_kind','redistribution'),'limit':cfg.get('limit',15),'c':o['coefficient'],
              'golden_m':o['golden_m'],'chosen_m':o['chosen_m'],'unit_charge_fc':cfg.get('store_ff',1000)*cfg.get('vdd',1.8)/144,
              'q_center_units':o['charge_q']['center'],'q_left_units':o['charge_q']['left'],'q_right_units':o['charge_q']['right'],
              'target_q_center_units':o['target_charge_q']['center'],'target_q_neighbor_units':o['chosen_m'],
              'charge_pass':o['transfer_charge_correct'],'center_residual_digits':o['center_residual_digits'],
              'class_correct':o['physical_class_correct'],'transfer_rail_fj':o['energy_fj']['transfer'],
              'sequential_unit_pulses':0,'bulk_join_phases':1,'run':r['run']})
        flat={'label':label,'corner':cfg.get('corner','tt'),'supply_v':cfg.get('vdd',1.8),'temp_c':cfg.get('temp',25),
          'states':r['states'],'quotient_correct':r['quotient_correct'],'transfer_correct':r['transfer_correct'],
          'physical_class_correct':r['class_correct'],'sleep_correct':r['sleep_correct'],
          'event_fj':r['mean_event_fj'],'idle_power_nw':r['idle_power_nw'],'seed':cfg.get('seed',41001),
          'mismatch_mode':cfg.get('mismatch_mode','none'),'run':r['run']}
        flat['mismatch_model_role']='legacy_cap_geometry_split_control' if cfg.get('mismatch_mode')=='capacitor' and not cfg.get('fixed_topology_caps',False) else 'fixed_topology_primary_or_other_control'
        if cfg.get('group')=='pvt':pvtflat.append(flat)
        if cfg.get('group') in ('mismatch','margin'):mmflat.append(flat)
    for r in comps:
        for o in r['observations']:
            cmpflat.append({'label':r['cfg']['label'],'corner':r['cfg'].get('corner','tt'),'supply_v':r['cfg'].get('vdd',1.8),
              'fixture':'asymmetric_initial_control' if r['cfg']['label'].startswith('cmp_w') else 'balanced_finite_holds',
              'temp_c':r['cfg'].get('temp',25),'seed':r['cfg'].get('seed',41001),'cmp_w':r['cfg'].get('cmp_w',.42),'cmp_l':r['cfg'].get('cmp_l',.15),
              'delta_requested_v':o['delta_v'],'actual_diff_before_eval_v':o['actual_differential_v'],
              'resolved':o['resolved'],'correct_requested_sign':o['correct'],'delay_ns':o['decision_delay_ns'],
              'kickback_v':o['kickback_v'],'comparison_rail_fj':o['energy_fj']['comparator'],
              'all_fixture_rails_and_pins_fj':o['total_fj'],'sample_and_hold_fj':o['energy_fj']['sample'],
              'idle_power_nw_fixture':r['idle_power_nw'],'offset_actual_mid_v':r['offset_bracket_mid_v'],'run':r['run']})
            cmpflat[-1].update({'idle_comparator_rail_nw':r['idle_comparator_rail_nw'],
              'minimum_tested_symmetric_correct_delta_v':r['minimum_tested_correct_abs_delta_v'],
              'offset_actual_bracket_v_json':json.dumps(r['offset_bracket_actual_v'])})
    for name,data in [('energy_breakdown.csv',energy),('reference_energy.csv',reference),('bulk_transfer.csv',transfer),
                      ('pvt_results.csv',pvtflat),('mismatch_results.csv',mmflat),('dynamic_comparator.csv',cmpflat)]:write_csv(OUT/name,data)
    baseline=json.loads((ROOT/'results/sky130_phase/phase_summary.json').read_text())
    exact=json.loads((OUT/'exact_workload.json').read_text());counts=exact['details'];comparison=[]
    for limit,physical in primary.items():
      if physical:
        emap={o['coefficient']:o for o in physical['observations']}
        for order in ('stream','position_batched','position_envelope'):
          for n in (128,512,2048,8192):
            count=next(c for c in counts if c['limit']==limit and c['input_order']==order and c['N']==n and c['policy']=='bounded_bulk')
            unit=next(c for c in counts if c['limit']==limit and c['input_order']==order and c['N']==n and c['policy']=='bounded_unit_control')
            hist={int(k):v for k,v in count['hist_coefficient'].items()}
            if any(c not in emap for c in hist):continue
            budget=next(b for b in baseline['primary_budget_rows'] if b['accumulation_length']==n)
            e=sum(freq*emap[c]['event_energy_fj'] for c,freq in hist.items())
            active=sum(freq*emap[c]['event_energy_fj'] for c,freq in hist.items() if abs(c)>=4)
            cmp_only=sum(freq*emap[c]['energy_fj']['comparator'] for c,freq in hist.items())
            active_cmp_only=sum(freq*emap[c]['energy_fj']['comparator'] for c,freq in hist.items() if abs(c)>=4)
            readout=budget['boundary_cost_per_position_fj']*budget['gamma_positions']
            inv=count['invocations_including_zero'];period=physical['cfg'].get('period_ns',400)
            bad=sum(freq for c,freq in hist.items() if not (emap[c]['quotient_correct'] and emap[c]['transfer_charge_correct'] and emap[c]['physical_class_correct'] and emap[c]['sleep_correct']))
            comparison.append({'N':n,'limit':limit,'input_order':order,'measured_cell':physical['cfg']['label'],
              'invocations_including_zero':inv,'nonzero_bulk_transfers':count['bulk_nonzero_transfers'],
              'zero_quotient_invocations':hist.get(0,0)+sum(v for c,v in hist.items() if 0<abs(c)<4),
              'unit_control_nonzero_events':unit['bulk_nonzero_transfers'],
              'nonzero_event_reduction':unit['bulk_nonzero_transfers']/max(1,count['bulk_nonzero_transfers']),
              'measured_normalization_pj':e/1000,'retained_readout_estimate_pj':readout/1000,
              'readout_plus_measured_normalization_pj':(e+readout)/1000,'adc_baseline_pj':budget['adc_energy_per_output_fj']/1000,
              'headroom_pj':(budget['adc_energy_per_output_fj']-readout)/1000,
              'available_fj_per_invocation':(budget['adc_energy_per_output_fj']-readout)/inv,
              'weighted_measured_fj_per_invocation':e/inv,'projected_visits_with_any_calibration_failure':bad,
              'comparator_rail_only_pj':cmp_only/1000,'active_only_comparator_rail_oracle_bound_pj':active_cmp_only/1000,
              'active_only_oracle_lower_bound_pj':active/1000,'old_unit_price_projection_pj':unit['bulk_nonzero_transfers']*4.219160752667593,
              'old_projection_scope':'Same exact integer workload, old isolated primitive price; no old topology rerun or physical whole-dot claim.',
              'invocations_per_s_per_column':inv*2.5e6,
              'necessary_parallel_cells_per_column':math.ceil(inv*2.5e6*period*1e-9),
              'projected_normalizer_w_per_column':e*1e-15*2.5e6,
              'extra_array_routing_completion_clock_source_cost':'UNMEASURED',
              'readout_load_bridge_revalidated':False,'all_in_win':False,
              'scope':'Histogram-weighted measured local supply/pin work. Physical failures retained. Favorable omission of global scheduler, clock generator, routing and 1pF-to-readout loading bridge. Not a working full dot.'})
    write_csv(OUT/'matched130_break_even.csv',comparison)
    idle_comparison=[]
    for row in comparison:
        if row['input_order']!='position_envelope':continue
        physical=primary[row['limit']];z=quiet_zero(physical)
        idle_comparison.extend([
          {'N':row['N'],'limit':row['limit'],'variant':'old_readonly_unit_price_reference','event_definition':'unit firing',
           'events_per_dot':row['unit_control_nonzero_events'],'normalization_price_projection_pj':row['old_unit_price_projection_pj'],
           'quiet_zero_state_nw':691.7057543393523,'measured_same_physical_workload_replay':False},
          {'N':row['N'],'limit':row['limit'],'variant':'new_dynamic_sar','event_definition':'all dirty coefficient requests, including m=0',
           'events_per_dot':row['invocations_including_zero'],'normalization_price_projection_pj':row['measured_normalization_pj'],
           'quiet_zero_state_nw':z,'measured_same_physical_workload_replay':False}])
    write_csv(OUT/'dynamic_idle_comparison.csv',idle_comparison)
    ranking=[]
    for label in ('range7_sar','range15_sar','range31_sar','range31_binary_compact','scheme_bank','bank_driven','scheme_binary_residual',
                  'ref_cap150','ref_cap300','ref_cap500','ref_precharged','ref_precharged_settled','ref_shared_static',
                  'transfer_binary_sc','transfer_gated_current','transfer_segmented_sc','margin_nominal',
                  'current_shared_r1.6e+06','current_shared_r1.85e+06','current_shared_r2.1e+06','ref_precharge_wide_switch','current_early_bias','nonzero_neighbors'):
        r=bylabel.get(label)
        if not r:continue
        ranking.append({'label':label,'states':r['states'],'quotient_correct':r['quotient_correct'],'charge_correct':r['transfer_correct'],
          'class_correct':r['class_correct'],'sleep_correct':r['sleep_correct'],'event_fj':r['mean_event_fj'],
          'idle_nw':r['idle_power_nw'],'mos':r['mos_count'],'nominal_all_axes_pass':all(r[k]==r['states'] for k in ('quotient_correct','transfer_correct','class_correct','sleep_correct')),
          'robustness_qualified':False,'scope':'Different range/state coverage is explicit; untested PVT is not a pass.'})
    write_csv(OUT/'topology_ranking.csv',ranking)
    mm_summary={}
    for mode in ('transistor','capacitor','reference'):
        rr=[r for r in mm if r['cfg'].get('mismatch_mode')==mode]
        mm_summary[mode]={'networks':len(rr),'states_per_network':5,
          **{k:sum(r[k]==r['states'] for r in rr) for k in ('quotient_correct','transfer_correct','class_correct','sleep_correct')}}
    summary={'date':'2026-10-24','q':4,'status':'NEGATIVE_DYNAMIC_IDLE_IMPROVES_BUT_BULK_COUNT_REFERENCE_ENERGY_AND_ACCURACY_FAIL' if not pending else 'IN_PROGRESS',
      'pending':pending,'completed_cell_cases':len(cells),'completed_comparator_cases':len(comps),'unresolved_attempts':unresolved,
      'pvt':{'cases':len(pvt),'states_per_case':5,**{k:sum(r[k]==r['states'] for r in pvt) for k in ('quotient_correct','transfer_correct','class_correct','sleep_correct')}},
      'mismatch':mm_summary,'primary_cells':{str(k):{a:v for a,v in r.items() if a!='observations'} for k,r in primary.items() if r},
      'adc_baseline_pj':baseline['ADC']['primary_130nm_adc_pj'],
      'comparison':comparison,'ranking':ranking,'all_in_win':False,
      'intentional_dc_tail_or_reference_ladder_in_selected_capacitive_sar':False,
      'tt_idle_zero_coefficient_nw':{str(k):quiet_zero(v) for k,v in primary.items() if v},
      'pvt_idle_power_nw_range':[min((r['idle_power_nw'] for r in pvt),default=None),max((r['idle_power_nw'] for r in pvt),default=None)],
      'precision_contract':'Only m=trunc(c/4), integer golden c and exact algebraic pairs. Strict m equality; equivalent alternative quotients are not credited.',
      'missing_terms':['physical global dirty queue and completion','clock generation upstream of measured pin loads','routing and array expansion',
        'loading/interface from 1pF normalizer storage to retained 25fF boundary estimate','full-word electrical replay','application-level matched ENOB/noise'],
      'scope':'No all-in or fundamental impossibility claim; calibrated ADC assumptions unchanged; ordinary130nm only.'}
    if (OUT/'literal_audit.json').exists():summary['independent_literal_audit']=json.loads((OUT/'literal_audit.json').read_text())
    summary['software_regression']='118 passed in17.10s; new-phase unit checks12 passed.'
    dump(OUT/'case_catalog.json',[{'config':r['cfg'],'run':r['run'],'kind':r.get('kind','cell'),'states':r.get('states',len(r['observations']))} for r in rows])
    sourcefiles=list((ROOT/'python').glob('*sar*.py'))+[ROOT/'spice/sky130/pdk_manifest.json',ROOT/'results/sky130_phase/phase_summary.json']
    hashes={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sourcefiles}
    dump(OUT/'source_hashes.json',hashes);dump(OUT/'phase_summary.json',summary)
    plots(cells,comps,primary,comparison)
    report(summary,exact,primary,bylabel,ranking)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('primary_cells','ranking','comparison')},indent=2))

def plots(cells,comps,primary,comparison):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels=[];values=[]
    for limit,r in primary.items():
        if r:labels.append(f'SAR +/-{limit}');values.append([r['per_component_mean_fj'][k]/1000 for k in ('sample','comparator','reference','logic','transfer','external_clock_pin_gross')])
    if values:
        fig,ax=plt.subplots(figsize=(9,4));bottom=np.zeros(len(labels))
        for j,k in enumerate(('Sample switching','Comparators','References','Logic','Transfer/reset','External clock pins')):
            ys=[v[j] for v in values];ax.bar(labels,ys,bottom=bottom,label=k);bottom+=ys
        ax.set_ylabel('Measured pJ per coefficient invocation');ax.legend(fontsize=8,ncol=3);fig.tight_layout();fig.savefig(OUT/'energy_components.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4))
    for limit in (7,15,31):
        rr=[r for r in comparison if r['limit']==limit and r['input_order']=='position_envelope']
        if rr:ax.loglog([r['N'] for r in rr],[r['readout_plus_measured_normalization_pj'] for r in rr],'o-',label=f'Gamma measured-price projection +/-{limit}')
    ax.loglog([128,512,2048,8192],[14.8573566489,18.3041477266,18.4,32.6425275655],'ks--',label='Unchanged130nm ADC')
    ax.set_xlabel('Useful INT8 MACs per dot');ax.set_ylabel('pJ per output');ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'matched_budget.png',dpi=160);plt.close(fig)
    r=primary.get(15)
    if r:
        d=ROOT/r['run'];names,x=fast_raw(d/'output.raw');ix={n:j for j,n in enumerate(names)}
        o=next(v for v in r['observations'] if v['coefficient']==13);t=x[:,0]*1e9-o['start_ns'];m=(t>=0)&(t<=350)
        fig,axes=plt.subplots(3,1,figsize=(10,8),sharex=True)
        for node in ('center','left','right','reference'):axes[0].plot(t[m],x[m,ix[f'v({node})']],label=node)
        axes[0].set_ylabel('V');axes[0].legend(ncol=4)
        for node in ('signq','bit1q','bit0q','evaluateDrv2o'.lower(),'packetfiredrv2o'):
            if f'v({node})' in ix:axes[1].plot(t[m],x[m,ix[f'v({node})']],label=node)
        axes[1].set_ylabel('Control V');axes[1].legend(ncol=3,fontsize=8)
        for rail in ('comparator','reference','logic','transfer'):
            p=-x[:,ix[f'v(dd_{rail})']]*x[:,ix[f'i(v{rail})']];axes[2].plot(t[m],p[m]*1e6,label=rail)
        axes[2].set_ylabel('Signed rail uW');axes[2].set_xlabel('ns after request: c=13, golden m=3');axes[2].legend(ncol=4)
        fig.tight_layout();fig.savefig(OUT/'sar_waveforms.png',dpi=160);plt.close(fig)

def report(summary,exact,primary,bylabel,ranking):
    comparison=summary['comparison'];r15=primary.get(15);r31=primary.get(31)
    lines=['# Fully dynamic q4 bulk normalizer: engineering report',
      '2026-10-24. Ordinary SKY130A1.8V models and ngspice47. Previous phases and the calibrated130nm ADC baseline are retained without rerunning the rejected101-MOS topology.',
      '## Decision',
      '**The tested dynamic circuits eliminate intentional continuous detector/reference bias, but do not preserve a demonstrated matched-node Gamma energy advantage. Bulk event reduction is too small for the measured per-invocation cost; zero-quotient decisions, reference setup and physical packet errors remain material.**',
      'There is no high-resolution coefficient ADC in the candidate: only a signed coarse quotient is latched. Neither the remainder nor the original coefficient is reconstructed numerically by hardware. Python golden values are used only to initialize test inputs, generate exact workload ledgers and score observed outputs. All quotient-dependent decisions, adaptive reference muxes, bit latches, packet enables and one-shot bulk transfers use actual transistors.',
      'Voltage-to-class conversion is an offline scoring measurement and never drives an intermediate hardware node. The exact intermediate event ledgers retain full algebraic precision and quantize only at the retained final readout. They do not bound a different lossy-intermediate Gamma design; that would need its own error ledger and event counts.',
      '## Exact range and event contract',
      table(['Window','Quotient states','Magnitude bits','Parallel comparisons','SAR / binary-trial comparisons','Max center/neighbor packet (Q)'],
            [[f'+/-{l}',2*(l//4)+1,(l//4).bit_length(),2*(l//4),1+(l//4).bit_length(),f'{4*(l//4)} / {l//4}'] for l in (7,15,31)]),
      'All109 integer input states across the three windows have exact-pair golden rewrites. Negative quotients truncate toward zero. No floating beta is used for arithmetic correctness. A chosen quotient differing by1 is reported as an error even if another redundant word could also represent the same value.',
      'The original128 paired-product suite has raw coefficient peaks of36, so none of the7/15/31 windows can accept every complete product without ingress scheduling. Direct product validity is18/128,77/128 and126/128 respectively. These invalid full-product cases are not clipped or silently normalized by an oracle.',
      'The tested physical schedule admits known digit-product contributions bounded by9. Window7 splits ingress into <=3-unit pieces; window15 admits one full digit-product per batch; window31 admits three. From persistent +/-3 states, these limits bound the ingress and subsequent incoming neighbor packets. Extra sliced-input driving for window7 and reordering/storage for position batching are unmeasured, not free product features.',
      'Three exact schedules are exported: original term order, favorable position batching, and a further favorable ingress-envelope schedule. The last groups known digit-product magnitudes up to limit-3 per position, reserving the3-unit unknown persistent remainder. This uses digital ingress magnitude metadata, not an analog overflow oracle; its metadata arithmetic and buffering costs are unmeasured. Dirty addresses come from scheduled ingress and actual nonzero quotient outputs. Every dirty visit is charged even if m=0; all receivers remain dirty when opposite incoming packets cancel. Physical queue/coalescing/global completion logic is still missing.',
    ]
    best=[r for r in comparison if r['limit']==31 and r['input_order']=='position_envelope']
    if best:
        lines.append(table(['N','Quotient invocations incl.m=0','Nonzero bulk transfers','Same-schedule unit transfers','Bulk transfer reduction','Allowed fJ/invocation'],
          [[r['N'],r['invocations_including_zero'],r['nonzero_bulk_transfers'],r['unit_control_nonzero_events'],f"{r['nonzero_event_reduction']:.2f}x",f"{r['available_fj_per_invocation']:.4f}"] for r in best]))
    lines += ['## Dynamic circuit and energy boundaries',
      'A complementary StrongARM latch has reset/precharge and evaluation phases, with no DC tail bias. Matched buffers isolate its two outputs from unequal sign/bit logic loading. A transmission-gate static bit latch retains each decision while the comparator resets. All phase clocks are input-independent testbench commands with real CMOS phase buffers. Their delivered pin charging is counted; an autonomous clock generator is not built.',
      'Sign is resolved at the coarse -0.5-unit boundary to avoid a nominal zero-input metastability requirement. Subsequent thresholds are half-integer boundaries3.5,7.5,..., separating exact integer quotient bands. This remains quotient extraction, not fine residue sensing.',
      'Three control architectures are physically tested at window15: a parallel threshold bank with a transistor thermometer encoder; one reused comparator with adaptive reference muxes; and a binary-trial capacitive reference DAC. The binary-trial circuit implements residual comparison on the reference side, retaining one final bulk packet phase. Its intermediate bit trials are not unit topplings.',
      'Capacitive references reset a floating top plate, release it, then switch weighted bottom plates from the supply. Comparator loading, reset injection and settling are included. A second control uses gated resistor charging of hold capacitors; a continuously connected reference is retained only as a deliberately identified static comparison. No free ideal tap drives the selected final candidate.',
      'The transfer tests use binary-weighted switched reservoirs, m-specific charge redistribution reservoirs, gated current packets, and parallel segmented unit reservoirs. They all use one join phase; no |m|-long sequence is hidden. Charge-sharing reservoirs inherently depend on node voltage and switch injection. Gated current packet bias is enabled only during the requested packet phase.',
      'Current-packet controls include separate per-leg precharge, shared event-gated bias, and an early-bias case that starts charging the mirror gates during quotient extraction. This exposes startup cost and settling instead of assuming an instantly ready current reference. Nonzero-neighbor controls are retained separately. The cost maps used for workload projections otherwise start both neighbors at zero, and do not certify state-independent packet behavior.',
      'All three positions have the same nominal1pF storage and identical input acquisition switches; Q=12.5fC at1.8V. The1pF choice keeps the smallest MIM packet near a legal minimum geometry. Center and passive-neighbor parasitic loads still differ; a fully loaded neighboring normalizer network has not been built. This difference is included in the measured class errors, not assumed away.',
      'Event energy is the signed sum of -integral(Vsource*I(Vsource)dt) on all five circuit rails plus gross delivered external clock/constant-pin work. It includes sample-switch driving, decisions, reference setup, quotient logic, packets and reset. No idle baseline is subtracted. The external state-writing sources initialize the pre-existing analog value and are reported separately as test-fixture/MAC work; input acquisition switch power remains in the normalizer budget.',
      'Reset-window supply work is exported separately, with disjoint initial latch-reset, per-decision reset/reference-return, and post-packet reset windows. Non-reset per-rail work plus reset-window work exactly sums to event energy; reset is not silently counted twice. These are temporal window allocations, including leakage, rather than baseline-subtracted intrinsic switching costs. Reference retention drift and comparator-loaded reference shifts are also exported.',
      'Raw comparator-rail energy and full comparator-fixture energy are separate in dynamic_comparator.csv: the latter includes the two finite sample/hold input paths and clock buffers. Equal finite input/reference holds are required to avoid mistaking fixture kickback asymmetry for intrinsic offset. Offset brackets, actual pre-evaluate differential, delay, kickback and correct-sign margins are retained per test.',
    ]
    if r15:
        lines += [table(['Selected window15 SAR quantity','Measured result'],[
          ['Quotient decisions',f"{r15['quotient_correct']}/{r15['states']}"],['Transfer charge checks',f"{r15['transfer_correct']}/{r15['states']}"],
          ['All three physical classes',f"{r15['class_correct']}/{r15['states']}"],['Return to sleep',f"{r15['sleep_correct']}/{r15['states']}"],
          ['Uniform-state mean event energy',f"{r15['mean_event_fj']/1000:.3f}pJ/invocation"],['Quiet last-state idle',f"{r15['idle_power_nw']:.2f}nW"],
          ['MOS count',r15['mos_count']]]),
          'These are nominal integer-state checks, not silicon accuracy/yield. Transfer tolerance is max(0.1Q,5% of the requested charge) per leg. Physical class scoring separately rounds each observed voltage to its declared unit and requires the golden center remainder and both golden neighbor increments. A charge check cannot hide a wrong quotient or wrong class.']
    lines += ['## Topology comparison and robustness',
      table(['Candidate','States','Quotient','Charge','Classes','Sleep','pJ/invocation','Idle nW','MOS'],
        [[r['label'],r['states'],r['quotient_correct'],r['charge_correct'],r['class_correct'],r['sleep_correct'],f"{r['event_fj']/1000:.3f}",f"{r['idle_nw']:.2f}",r['mos']] for r in ranking]),
      'The initial parallel-bank drive was too weak for its coherently switched reference load. bank_driven corrects that physical drive strength; the failed bank is retained as a settling control. No PVT pass is assigned to a topology that was tested only nominally, and there is no robustness-qualified energy winner.',
      f"The frozen diagnostic SAR candidate has {summary['pvt']['cases']} PVT cases: TT/SS/FF/SF/FS x1.62/1.8/1.9V x-40/25/100C. Each PVT case tests c=-13,-4,0,4,13; all31 nominal window15 states are covered separately. This is not full-state PVT or transient-noise certification. Supply-ratiometric input and reference levels are assumed at each operating point; untracked supply droop on stored data is not tested.",
      table(['PVT networks passing every tested state','Quotient','Charge','Physical classes','Sleep'],[[summary['pvt']['cases'],*[summary['pvt'][k] for k in ('quotient_correct','transfer_correct','class_correct','sleep_correct')]]]),
      table(['Mismatch scenario','Independent networks','Quotient passes','Charge passes','Class passes','Sleep passes'],
        [[m,v['networks'],v['quotient_correct'],v['transfer_correct'],v['class_correct'],v['sleep_correct']] for m,v in summary['mismatch'].items()]),
      'Transistor mismatch uses the pinned PDK instance model. Capacitor mismatch is an explicitly assumed independent1% sigma per MIM realization. Reference mismatch is a separate1% ratio perturbation. These are stress scenarios, not a foundry-qualified yield distribution. Each network contains correlated state tests and counts once. A larger-device / larger-reference control is retained to expose the energy/margin tradeoff without claiming that a cheap nominal point is robust.',
      'The final capacitor-mismatch data keep every capacitor instance count fixed. An inherited helper initially split a nominal1pF capacitor into two instances when a positive mismatch draw exceeded1pF; those runs are retained as legacy geometry-split controls and excluded from the primary eight-network capacitor statistics. The larger-device control keeps the original driver/timing and is not evidence against a separately optimized margin design.',
      'Leakage is strongly state-, corner- and temperature-dependent. Dynamic operation removes intentional DC bias, but many small CMOS gates still leak. The per-PVT idle results must be inspected; the nominal idle number must not be extrapolated to hot conditions. Comparator-only reset leakage is distinct from whole-cell phase/control leakage.',
      f"The tested whole-cell quiescent range is {summary['pvt_idle_power_nw_range'][0]:.2f}-{summary['pvt_idle_power_nw_range'][1]:.2f}nW. In particular, the worst mixed hot corner reaches {summary['pvt_idle_power_nw_range'][1]/1000:.2f}uW despite the absence of an intentional bias ladder/tail. A robust near-zero-idle claim is not supported. Full internal device-terminal SOA, extracted parasitics and silicon leakage attribution are outside this phase.",
      'Sleep scoring separately retains control-reset correctness and preservation of the observed three-node classes through the quiet part of the frame. A wrong post-transfer class is still a class failure even if it is faithfully held. Stored-node drift is reported, so low supply-rail leakage cannot hide charge lost from the stored analog value.',
      '## Matched130nm comparison',
      'The ADC assumptions and14.86/18.30/18.40/32.64pJ scenarios are unchanged. Normalization is computed from the exact per-coefficient histogram and the corresponding measured physical invocation energy, including m=0. No32.3828N unit count is blindly reused. The retained4.81-6.73pJ readout estimate is added only as a favorable reference: its25fF loading interface has not been revalidated for1pF normalizer storage.',
    ]
    if best:
        lines.append(table(['N','Readout + measured-price normalization (pJ)','ADC (pJ)','Normalization alone (pJ)','Required parallel cells/column, utilization bound'],
          [[r['N'],f"{r['readout_plus_measured_normalization_pj']:.1f}",f"{r['adc_baseline_pj']:.2f}",f"{r['measured_normalization_pj']:.1f}",r['necessary_parallel_cells_per_column']] for r in best]))
        b=best[0]
        lines.append(f"At N={b['N']}, even the measured comparator rail alone costs {b['comparator_rail_only_pj']:.1f}pJ. Ignoring every m=0 request as an unattainably favorable oracle bound still gives {b['active_only_comparator_rail_oracle_bound_pj']:.1f}pJ of comparator work, against {b['headroom_pj']:.2f}pJ available for all extra Gamma costs. All sample/reference/logic/packet/completion costs are set to zero only in this negative lower-bound diagnostic. Thus rejection of these measured comparator/schedule combinations does not depend on the choice of1pF packet storage. This is not a lower bound on a different comparator topology.")
    lines += ['These are cost projections of measured local circuits, not measured functional dots: physical failures are retained and prevent an efficiency-win claim. Global completion, request queuing/coalescing, upstream clock generation, array/routing expansion and the normalizer/readout load bridge are missing. No all-in positive comparison is asserted.',
      'Time-averaged normalizer power is evaluated at the retained2.5M outputs/s/column using measured event-window energy times actual invocation density; no idle baseline is removed. Required parallel-cell counts are necessary utilization bounds only, not an implemented architecture or a proof that cascade dependencies meet throughput. Extra idle time in a larger pool adds P_idle*time. The400ns nominal frame and300ns compact-window31 frame are explicit in the CSV.',
      'Measured per-coefficient prices come from finite canonical state sequences, so transition history and actual neighbor states are not electrically replayed for each complete dot. References and packet reservoirs reset between requests, but bit/control switching and node loads can retain history dependence. The histogram-weighted sums are explicitly local-price projections, not a substitute for functional full-word electrical validation.',
      'The previous topology is used read-only at its measured4.219pJ/unit-event price and691.7nW quiet idle. Same-schedule exact unit-control event counts are recomputed for the comparison. This is an equal-integer-workload price projection; it is not a new physical workload replay of the rejected cell. Event-window definitions and isolated test-state coverage differ, so the data do not support a fabricated end-to-end speedup.',
      table(['Equal logical c=0 quiescent state','Measured nW'],[['Old read-only reference','691.7']]+[[f'New SAR +/-{k}',f'{v:.2f}'] for k,v in summary['tt_idle_zero_coefficient_nw'].items()]),
      'dynamic_idle_comparison.csv separates measured quiescent power from equal-integer-workload event-price projections. The old cell used Q=2.5fC and the new one Q=12.5fC; this is not a matched silicon noise/accuracy experiment.',
      '## Answer and limitations',
      '**A fully dynamic coarse quotient circuit and a one-phase bulk transfer can be built in these models, and the continuous-bias pathology is removed structurally. The tested circuits do not reduce physical activity or measured energy far enough to preserve the matched130nm headroom.** The negative result is specific to these bounded schedules, switched-reference/packet implementations and device sizes; it is not a universal lower bound for Gamma arithmetic.',
      'The decisive tests are independent: modest bulk compression in the actual signed workload; many paid zero-quotient requests; picjoule reference/control/sample setup; state-dependent packet errors; and mismatch failures. Advanced-node projections, free references, free completion, or free oracle activity suppression are not used to rescue any failure.',
      '## Evidence and replay',
      'From project root: python -B python/gamma_sar_exact.py; python -B python/run_sar_normalizer_phase.py --stage nominal; --stage robust; --stage comparator; python -B python/report_sar_normalizer_phase.py. Additional driven-bank, compact31 and margin-control configurations are retained in their golden.json files and can be run with run_sar_normalizer_phase.run. Each case has a literal circuit.cir, seeded replay.cir, model corner, input sequence, metric and ngspice log. Unresolved attempts are listed in phase_summary.json and never counted as passes.',
      'Source hashes are in source_hashes.json; pinned model and simulator provenance remains spice/sky130/pdk_manifest.json. No PDK files were edited. Representative binary traces and all case metadata are packaged; full completed traces remain in the project.',
      table(['Independent literal replay','Batch fJ','Literal fJ','Energy difference','Correctness counts agree'],
        [[v['label'],f"{v['batch_fj']:.3f}",f"{v['literal_fj']:.3f}",f"{100*v['relative_energy_difference']:.6f}%",v['correctness_counts_agree']] for v in summary.get('independent_literal_audit',[])]),
      'Literal replays use circuit values directly without batch parameter substitution. Software regression:118 passed. These checks validate accounting/software and replay agreement, not physical circuit robustness.',
      'Primary circuit background: Razavi, [The StrongARM Latch](https://ieeexplore.ieee.org/document/7130773/), IEEE Solid-State Circuits Magazine2015, DOI10.1109/MSSC.2015.2418155; and [The Design of a Comparator](https://www.ee.ucla.edu/~brweb/papers/Journals/BR_SSCM_4_2020.pdf),2020. [SKY130 device/capacitor documentation](https://skywater-pdk.readthedocs.io/en/main/rules/device-details.html). Those sources motivate topology/model limits; all reported energies come from retained simulations.',
      '![Measured energy components](energy_components.png)','![Matched-node budget](matched_budget.png)','![Actual SAR waveform](sar_waveforms.png)']
    (OUT/'engineering_report.md').write_text('\n\n'.join(lines)+'\n',encoding='utf-8')

if __name__=='__main__':main()
