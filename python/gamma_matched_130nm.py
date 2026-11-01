"""Technology/throughput/work-matched accounting over retained SKY130 evidence.

No new array, arithmetic family, neural workload or transistor optimization.
Finite-word checks reuse the existing exact q=4 reference. Unknown local
normalizer/control costs remain explicit sensitivities, never calibrated facts.
"""
from __future__ import annotations
import csv,json,math,itertools,statistics,shutil
from decimal import Decimal,localcontext
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from gamma_adc_baseline import ROOT,BASE,OUT,ADCModel,COUNTS,write_csv,dump,sha
from gamma_native_model import normalize,value_pair,power_pair,pair_multiply,shift
from gamma_regen_golden import pair_sign,finite_word

SKY=ROOT/'spice/sky130'
LENGTH_BITS=((128,6),(512,7),(2048,8),(8192,9))


def read_csv(path):return list(csv.DictReader(Path(path).open(encoding='utf-8')))


def pair_decimal(pair):
    with localcontext() as c:
        c.prec=70
        def dec(x):x=F(x);return Decimal(x.numerator)/Decimal(x.denominator)
        return float(dec(pair[0])+dec(pair[1])*(Decimal(2)+Decimal(3).sqrt()))


def precision_layout(mac_count,bits):
    """Rigorous truncation bound plus empirical existing-oracle checks.

For sign-uniform encodings of INT8 scalars, weighted absolute coefficient
mass is <=16384*N. Local q:1:1 rewrites cannot increase either beta-weighted
absolute mass. Thus every nonzero digit has |k|<=floor(log_beta(16384*N)).
This bounds leading positions without using floating beta for correctness.
"""
    radius=16384*mac_count;high=0
    while pair_sign((power_pair(4,high+1)[0]-radius,power_pair(4,high+1)[1]))<=0:high+=1
    # Infinite discarded tail <=3*beta^low/(beta-1); inverse(beta-1)=(beta-3)/2.
    limit=F((2*radius)**2,12*(2**bits)**2)
    for positions in range(1,64):
        low=high-positions+1
        tail=pair_multiply(power_pair(4,low),(F(-9,2),F(3,2)),4)
        square=pair_multiply(tail,tail,4)
        if pair_sign((square[0]-limit,square[1]))<=0:break
    # Dense deterministic integer grid, including extremes and near zero.
    samples=sorted(set([-radius,0,radius]+[(2*radius*j//512)-radius for j in range(513)]+list(range(-16,17))))
    traces=[];sumsq=(F(0),F(0));peak=(F(0),F(0))
    for n in samples:
        word,_=normalize({0:n},4)
        assert all(-high<=k<=high for k in word)
        result=finite_word(shift(word,-high),exponent=high,window=(-positions+1,0),guard_digits=2,policy='truncate')
        error=tuple(F(v) for v in result['exact_discarded_or_rounding_residual'])
        squared=pair_multiply(error,error,4);sumsq=tuple(x+y for x,y in zip(sumsq,squared))
        absolute=tuple(v if pair_sign(error)>=0 else -v for v in error)
        if pair_sign(tuple(x-y for x,y in zip(absolute,peak)))>0:peak=absolute
        assert pair_sign(tuple(x-y for x,y in zip(absolute,tail)))<=0
        traces.append({'mac_count':mac_count,'enob_target':bits,'input_integer':n,'positions':positions,'common_block_exponent':high,
            'output_word_json':json.dumps(result['state'],sort_keys=True),
            'exact_output_A':result['exact_output'][0],'exact_output_B':result['exact_output'][1],
            'exact_error_A':str(error[0]),'exact_error_B':str(error[1])})
    mse=tuple(v/len(samples) for v in sumsq)
    assert pair_sign((mse[0]-limit,mse[1]))<=0
    actual_rmse=math.sqrt(max(0,pair_decimal(mse)))
    row={'mac_count':mac_count,'enob_target':bits,'int8_worst_case_fullscale_radius':radius,'high_position':high,
        'guaranteed_truncated_positions':positions,'exact_fullword_positions':2*high+1,
        'discard_tail_bound_in_original_units':pair_decimal(tail),'ideal_adc_quantization_rmse_units':float(limit)**.5,
        'empirical_rmse_units':actual_rmse,'empirical_max_error_units':pair_decimal(peak),'samples':len(samples),
        'empirical_quantization_equivalent_bits':math.log2(2*radius/(math.sqrt(12)*actual_rmse)) if actual_rmse else None,
        'exact_tail_bound_pass':True,'exact_finite_word_checks_pass':True,
        'precision_scope':'Representation/truncation only; NOT physical ENOB, BER or matched application accuracy. Common block exponent; no per-scalar exponent oracle.'}
    return row,traces


def conventional_choice(adc,bits,aggregate,anchors,mode='model'):
    options=[]
    for n in COUNTS:
        try:r=adc.estimate(bits,130,aggregate,n)
        except AssertionError:continue
        options.append({**r,'adc_source':'pinned_statistical_model','anchor_id':'','rate_basis':'upstream aggregate independent conversions contract'})
        if mode=='optimistic_anchor_envelope':
            for a in anchors:
                if a['enob_from_sndr'] is None or a['enob_from_sndr']<bits:continue
                if a['osr'] is None or a['osr']>1 or not a['area_um2'] or not a['energy_per_sample_pj']:continue
                if not a['independent_value_rate_bound_hz'] or a['independent_value_rate_bound_hz']<aggregate/n:continue
                e=a['energy_per_sample_pj']
                options.append({**r,'energy_per_conversion_pj':e,'energy_per_conversion_fj':e*1000,
                    'area_per_adc_um2':a['area_um2'],'total_area_um2':n*a['area_um2'],
                    'bank_power_mw':e*aggregate*1e-9,'adc_source':'compatible_published_130nm_anchor',
                    'anchor_id':a['anchor_id'],'rate_basis':'Published clock/value-rate ceiling; published per-sample energy, optimistic idle gating. Not measured at requested lower clock.'})
    if not options:return None
    return min(options,key=lambda r:(r['energy_per_conversion_pj'],r['total_area_um2']))


def build_amortization(adc,anchors):
    rows=[]
    # Strong conventional baseline: one full multilevel dot read, with declared
    # final quantization. RAELLA pairings motivate the range, not a universal
    # theorem that N rows require6/7/8/9 exact bits.
    for n,bits in LENGTH_BITS:
        for occupancy,outputs,word_rate,mode in itertools.product((1,.25),(1,64),(2.5e6,5e6),('model','optimistic_anchor_envelope')):
            useful=max(1,int(n*occupancy));aggregate=outputs*word_rate
            choice=conventional_choice(adc,bits,aggregate,anchors,mode)
            assert choice is not None
            rows.append({'scenario':'one_shot_quantized_dot','accumulation_capacity':n,'useful_scalar_macs_per_output':useful,
                'utilization':occupancy,'output_columns':outputs,'outputs_per_second_per_column':word_rate,
                'adc_conversions_per_output':1,'useful_scalar_macs_per_bank_step':useful*outputs,
                'adc_conversions_per_bank_step':outputs,'adc_energy_per_scalar_mac_fj':choice['energy_per_conversion_fj']/useful,
                'adc_cost_per_output_pj':choice['energy_per_conversion_pj'],'baseline_mode':mode,**choice,
                'dynamic_range_assumption':'INT8 full-dot multilevel analog input/weights; final b-ENOB quantization over +/-16384*capacity. DAC/programming costs excluded; strong optimistic conventional array assumption.',
                'precision_note':'Published128/512/2048/8192 with6/7/8/9 ENOB are statistical architecture settings, not exact-integer bit counts.'})
    # Exact bit-serial reference:64 unsigned bit-product planes for signed INT8;
    # signed binary place weights are accumulated digitally. Tile to avoid
    # count overflow with b conversion bits. ENOB alone does not prove bit BER.
    for n in (128,512,2048,8192):
        for bits in range(6,11):
            count_capacity=2**bits-1;chunks=math.ceil(n/count_capacity);conversions=64*chunks
            aggregate=64*2.5e6*conversions
            for mode in ('model','optimistic_anchor_envelope'):
                choice=conventional_choice(adc,bits,aggregate,anchors,mode)
                common={'scenario':'exact_int8_64_bitplanes','accumulation_capacity':n,'useful_scalar_macs_per_output':n,
                    'utilization':1,'output_columns':64,'outputs_per_second_per_column':2.5e6,
                    'adc_conversions_per_output':conversions,'count_capacity_per_read':count_capacity,'chunks_per_dot':chunks,
                    'useful_scalar_macs_per_bank_step':64*n,'adc_conversions_per_bank_step':64*conversions,'baseline_mode':mode,
                    'dynamic_range_assumption':'8x8 signed INT8 bitplanes; each unsigned partial count fits0..2^b-1; digital signed place-value accumulation. No statistical clipping.',
                    'precision_note':'Exact ideal arithmetic contract; real converter class error/BER not modeled.'}
                if choice is None:rows.append({**common,'enob':bits,'technology_nm':130,'aggregate_samples_s':aggregate,'status':'REJECTED: no<=64ADC design inside10GHz per-ADC cap'});continue
                rows.append({**common,**choice,'adc_energy_per_scalar_mac_fj':conversions*choice['energy_per_conversion_fj']/n,
                    'adc_cost_per_output_pj':conversions*choice['energy_per_conversion_pj'],'status':'model_accepted'})
    write_csv(OUT/'adc_cim_amortization.csv',rows);return rows


def retained_evidence():
    s=json.loads((SKY/'results/sky130_phase_summary.json').read_text());old=json.loads((ROOT/'spice/regen_results/phase_summary.json').read_text())
    pvt=[r for r in s['followup'] if r['group']=='pvt_energy']
    nominal=next(r for r in pvt if r['corner']=='tt' and r['vdd']==1.8 and r['temp_c']==25)
    metrics=json.loads((ROOT/nominal['run']/'metrics.json').read_text())
    energy=[{'checkpoint':i,**o} for i,o in enumerate(metrics['observations'])]
    write_csv(OUT/'regen_energy.csv',[{'checkpoint':r['checkpoint'],'previous_digit':r['previous_digit'],'digit':r['digit'],
        'energy_gross_vdd_plus_pin_fj':r['measured_known_path_fj'],
        'signed_vdd_work_fj':sum(r['signed_port_work_fj'].get(k,0) for k in ('reference','interface','writer','buffer','comparator')),
        'stored_v':r['stored_v'],'loaded_v':r['driven_v'],'loaded_error_digit_units':r['driver_error_units'],
        'next_mac_class':r['mac_class'],'correct':r['next_mac_correct'],**{'energy_'+k+'_fj':v for k,v in r['energy_components_fj'].items()}} for r in energy])
    shutil.copyfile(SKY/'results/regen_transition_matrix.csv',OUT/'regen_transition_matrix.csv')
    shutil.copyfile(SKY/'results/pvt_summary.csv',OUT/'pvt_results.csv')
    shutil.copyfile(SKY/'results/monte_carlo_trials.csv',OUT/'mismatch_results.csv')
    write_csv(OUT/'loaded_chain_results.csv',[r for r in s['followup'] if r['group'].startswith('chaining')])
    manifest=json.loads((SKY/'pdk_manifest.json').read_text())
    manifest.update(existing_model_reused=True,new_device_installation=False,exact_devices_used=['sky130_fd_pr__nfet_01v8','sky130_fd_pr__pfet_01v8','sky130_fd_pr__cap_mim_m3_1'],
        pdk_disk_bytes_including_git=sum(p.stat().st_size for p in (ROOT/'third_party/sky130-ngspice-models').rglob('*') if p.is_file()),
        evidence_source='spice/sky130/results/sky130_phase_summary.json',evidence_source_sha256=sha(SKY/'results/sky130_phase_summary.json'))
    dump(OUT/'sky130_device_manifest.json',manifest)
    comp=s['selected_comparator']['max_legal_center_energy_fj']
    measured={k:nominal.get('energy_'+k+'_fj',0) for k in ('reference','interface','writer','buffer','comparator','control_pin_charging_upper')}
    measured['known_hybrid_fj']=nominal['mean_energy_fj'];measured['additional_comparison_fj']=comp
    counts=old['completion'];counts={**counts,'control_activity_units_per_product':old['energy']['reference_point']['explicit event model']['completion_control_fj']/.01}
    return s,measured,counts


def boundary_cost(measured,decisions=3,reference_fanout=1,regeneration_fraction=1,driver_duty=1,extra_logic_fj=0):
    # Reference and continuously biased driver do NOT disappear when a write
    # is skipped. Even all current control-pin allowance remains charged.
    return (measured['reference']/reference_fanout+measured['buffer']*driver_duty+
        measured['comparator']+(decisions-1)*measured['additional_comparison_fj']+
        measured['control_pin_charging_upper']+regeneration_fraction*(measured['interface']+measured['writer'])+extra_logic_fj)


def norm_cost(n,counts,detect_fj,rewrite_fj,control_fj,density=1):
    return n*density*(counts['event_mean_counts']['comparator_evaluations']*detect_fj+
        counts['unit_rewrites_mean']*rewrite_fj+counts['control_activity_units_per_product']*control_fj)


def matched_tables(amortization,layouts,measured,counts):
    layout={r['mac_count']:r for r in layouts};rows=[];budgets=[]
    norm_modes={'free_normalization_lower_bound':(0,0,0),'legacy_UNCALIBRATED_costs':(.1,1.2,.01),
        'one_fj_coarse_detector_HYPOTHESIS':(1,1.2,.01),'measured_comparator_PROXY_rewrite_control_free':(measured['additional_comparison_fj'],0,0)}
    bases=[r for r in amortization if r['scenario']=='one_shot_quantized_dot' and r['utilization']==1 and r['output_columns']==64 and r['outputs_per_second_per_column']==2.5e6]
    for b in bases:
        n=b['accumulation_capacity'];p=layout[n]['guaranteed_truncated_positions'];adc_fj=b['adc_cost_per_output_pj']*1000
        for decisions in (1,3,6):
            cell=boundary_cost(measured,decisions)
            available=adc_fj-p*cell
            budgets.append({'accumulation_length':n,'enob':b['enob'],'baseline_mode':b['baseline_mode'],'adc_source':b['adc_source'],
                'gamma_positions':p,'physical_decisions_estimated':decisions,'adc_energy_per_output_fj':adc_fj,
                'boundary_cost_per_position_fj':cell,'per_position_budget_if_other_gamma_costs_zero_fj':adc_fj/p,
                'remaining_normalizer_routing_logic_budget_fj_per_dot':available,
                'overflow_detection_ceiling_fj_if_other_norm_costs_zero':available/(n*counts['event_mean_counts']['comparator_evaluations']),
                'rewrite_ceiling_fj_if_other_norm_costs_zero':available/(n*counts['unit_rewrites_mean']),
                'budget_scope':'Common analog array energy canceled; binary accumulation/requantization priced0 to favor conventional. Normalizer count model is retained per-product policy times N; shared policy is unbuilt.'})
        for p_mode,positions in [('precision_bound',p),('historical_17_positions',17)]:
            for decisions,fanout,mode,binary_fj in itertools.product((1,3,6),(1,8,64),norm_modes,(0,100,1000)):
                de,rw,ct=norm_modes[mode];cell=boundary_cost(measured,decisions,fanout);normalization=norm_cost(n,counts,de,rw,ct)
                gamma=positions*cell+normalization;conventional=adc_fj+binary_fj
                rows.append({'comparison_level':'B1 peripheral / B2 equal-array','accumulation_length':n,'enob':b['enob'],'useful_scalar_macs':n,
                    'baseline_mode':b['baseline_mode'],'adc_source':b['adc_source'],'adc_count':b['n_adcs'],'adc_frequency_hz':b['frequency_per_adc_hz'],
                    'adc_conversions_per_output':1,'adc_energy_per_output_fj':adc_fj,'binary_accumulation_plus_requantization_fj':binary_fj,
                    'gamma_position_policy':p_mode,'gamma_positions':positions,'detector_decisions_estimated':decisions,
                    'reference_fanout_hypothesis':fanout,'regenerations_per_dot':1,'normalization_pricing':mode,
                    'overflow_detect_fj':de,'rewrite_fj':rw,'control_activity_fj':ct,'gamma_boundary_fj_per_position':cell,
                    'gamma_normalization_fj_per_dot':normalization,'gamma_peripheral_fj_per_dot':gamma,
                    'conventional_peripheral_fj_per_dot':conventional,'baseline_minus_gamma_fj_per_dot':conventional-gamma,
                    'gamma_fj_per_useful_mac':gamma/n,'conventional_fj_per_useful_mac':conventional/n,
                    'normalization_activity_density':1,'extra_array_routing_fj':0,'physical_robustness_demonstrated':False,
                    'interpretation':'Conditional energy accounting; only one threshold transistorized. Three/six decisions model unbuilt extra reference/latch/control. Fanout>1 unvalidated; no physical-normalizer energy measured.'})
    write_csv(OUT/'matched_break_even.csv',rows);write_csv(OUT/'regeneration_budgets.csv',budgets)
    surfaces=[]
    for b in bases:
        if b['baseline_mode']!='optimistic_anchor_envelope':continue
        n=b['accumulation_capacity'];positions=layout[n]['guaranteed_truncated_positions'];adc_fj=b['adc_cost_per_output_pj']*1000
        for fanout,rewrite,density,detect in itertools.product((1,64),(0,1.2,10),np.logspace(-3,0,25),np.logspace(-3,2.3,32)):
            boundary=positions*boundary_cost(measured,3,fanout)
            gamma=boundary+norm_cost(n,counts,float(detect),rewrite,.01,float(density))
            surfaces.append({'accumulation_length':n,'enob':b['enob'],'gamma_positions':positions,'reference_fanout':fanout,
                'normalization_activity_density':float(density),'overflow_detect_fj':float(detect),'rewrite_fj':rewrite,
                'control_activity_fj':.01,'baseline_minus_gamma_fj_per_dot':adc_fj-gamma,
                'assumption':'Density<1 requires a new bounded shared/sparse normalization schedule; not established. Fanout>1 is an optimistic sharing hypothesis.'})
    write_csv(OUT/'normalization_break_even_surface.csv',surfaces)
    full=[]
    for b in bases:
        if b['baseline_mode']!='optimistic_anchor_envelope':continue
        n=b['accumulation_capacity'];positions=layout[n]['guaranteed_truncated_positions'];adc_fj=b['adc_cost_per_output_pj']*1000
        for expansion,array_fj,routing_fj,frequency,refresh,logic in itertools.product((1,4,16,counts['digit_products']),(0,1,10,100),(0,10,100,1000),sorted(set((1,math.ceil(n/128),n))),(1,.1,0),(0,20,100,1000)):
            gamma_boundary=positions*boundary_cost(measured,3,regeneration_fraction=refresh,extra_logic_fj=logic)*frequency
            delta_array=(expansion-1)*array_fj*n
            # Resource costs zero here: this is the optimistic array/routing
            # sensitivity, kept separate from the normalizer price surface.
            gamma=gamma_boundary+positions*routing_fj+delta_array
            full.append({'accumulation_length':n,'enob':b['enob'],'gamma_positions':positions,'representation_array_energy_factor':expansion,
                'common_conventional_array_fj_per_mac':array_fj,'extra_array_fj':delta_array,'routing_fj_per_position':routing_fj,
                'boundary_invocations_per_dot':frequency,'write_refresh_fraction':refresh,'reference_fanout':1,
                'gamma_minus_common_array_fj':gamma,'baseline_minus_gamma_fj_per_dot':adc_fj-gamma,
                'normalizer_energy_fj':0,'unknown_logic_fj_per_position':logic,'unknown_logic_energy_fj':positions*logic*frequency,'warning':'Optimistic sensitivity. Array energy factor is not digit-product count alone; refresh skips do not remove continuous buffer/reference cost. No measured sharing/skip controller.'})
    write_csv(OUT/'full_sensitivity.csv',full)
    # Exact bit-plane baseline is a separate precision contract. No assertion
    # that low-ENOB one-shot outputs equal a complete exact INT8 dot.
    exact=[]
    for b in amortization:
        if b['scenario']!='exact_int8_64_bitplanes' or b['baseline_mode']!='optimistic_anchor_envelope':continue
        n=b['accumulation_capacity'];positions=layout[n]['exact_fullword_positions']
        if b.get('status')!='model_accepted':exact.append({**b,'gamma_positions':positions});continue
        for name,(de,rw,ct) in norm_modes.items():
            gamma=positions*boundary_cost(measured,3)+norm_cost(n,counts,de,rw,ct)
            exact.append({**b,'gamma_positions':positions,'normalization_pricing':name,'gamma_peripheral_fj_per_dot':gamma,
                'baseline_minus_gamma_fj_per_dot':b['adc_cost_per_output_pj']*1000-gamma,
                'gamma_complete_dot_throughput_proven':False,'scope':'Equal useful signed INT8 MAC count and ideal exact arithmetic. No physical full-dot normalizer/parallelism demonstrated.'})
    write_csv(OUT/'exact_bitserial_comparison.csv',exact)
    return rows,budgets,surfaces


def main():
    manifest=json.loads((OUT/'baseline_manifest.json').read_text());assert manifest['reproduction_passed']
    from gamma_adc_baseline import survey_anchors
    _,_,anchors=survey_anchors();adc=ADCModel();amortization=build_amortization(adc,anchors)
    layouts=[];traces=[]
    for n,b in LENGTH_BITS:
        row,trace=precision_layout(n,b);layouts.append(row);traces.extend(trace)
    write_csv(OUT/'gamma_precision_layout.csv',layouts);write_csv(OUT/'finite_word_precision_golden.csv',traces)
    s,measured,counts=retained_evidence();rows,budgets,surfaces=matched_tables(amortization,layouts,measured,counts)
    definition={'ordered_phase':'ADC calibration and baseline definition completed before any new transistor run or optimization.',
        'q':4,'beta_exact':'2+sqrt(3)','digit_alphabet':list(range(-3,4)),
        'work_unit':'One scalar signed INT8 MAC; one output dot contains N useful MACs. Output columns cancel symmetrically. No per-event-vs-per-conversion comparison.',
        'primary_throughput':'64 output columns,2.5M dot results/s/column,160M conventional one-shot conversions/s bank aggregate. Gamma would need positions*64 parallel boundary cells; only boundary cadence, not full-dot normalization throughput, is demonstrated.',
        'ADC_cost_formula':'Nconversions*Econversion/NusefulMACs; ADC count alters per-ADC rate and area, never divides conversion energy again.',
        'B1':'Peripheral-only comparison including known boundary blocks; missing normalizer, extra detector, completion, routing and precision-interface costs are explicit.',
        'B2':'Equal analog array energy cancels explicitly. This is an assumption, not a circuit result.',
        'B3':'Separate normalization/energy-density and array expansion/routing/sharing/refresh sensitivities. No single winner forced.',
        'quantization_contract':'One-shot multilevel INT8 dot, b-ENOB final full-scale quantization; Gamma finite-word error bounded no worse than ideal ADC quantization RMS in representation only. Physical BER/ENOB match is NOT established.',
        'exact_alternative':'64 signed INT8 bitplanes and count-limited chunks; separate from lossy one-shot quantization. Invalid required converter rates are retained.',
        'normalization_count_source':'Retained q4 per-product event model, scaled by N; merged-dot normalization/control cost is additional and unmeasured. Alternative shared/sparse schedules require proof.',
        'old_threshold_status':'100fJ ADC /2.74fJ per-position budget is superseded for decision-making, retained unchanged in historical reports.',
        'measured_boundary':measured,'retained_counts':counts,'layouts':layouts,
        'new_transistor_optimization':False,'old_physical_phase_reused':True,
        'limitations':['five ideal threshold decisions','unbuilt class latch/decoder/completion','no transistor local-normalizer energy','no electrically validated reference fanout sharing','no demonstrated whole-dot throughput','no physical precision/yield match']}
    dump(OUT/'baseline_definition.json',definition)
    print(json.dumps({'layouts':layouts,'measured_boundary':measured,'baseline_B_ready':True,'matched_rows':len(rows),'surface_rows':len(surfaces)},indent=2),flush=True)


if __name__=='__main__':main()
