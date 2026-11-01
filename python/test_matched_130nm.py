"""Unit/rate/work-accounting guards for the new calibrated comparison."""
import csv,json,math
from fractions import Fraction as F
import pytest
from gamma_adc_baseline import ADCModel,OUT,BASE
from gamma_matched_130nm import boundary_cost,precision_layout,conventional_choice


@pytest.fixture(scope='module')
def adc():return ADCModel()


def test_published_example_energy_units(adc):
    r=adc.estimate(8,16,512e7,32)
    assert math.isclose(r['energy_per_conversion_pj'],.4655757205510538,rel_tol=1e-12)
    assert r['frequency_per_adc_hz']==160e6
    assert math.isclose(r['energy_per_conversion_fj'],465.5757205510538,rel_tol=1e-12)


def test_interleaving_does_not_double_amortize_energy(adc):
    a=adc.estimate(8,130,160e6,1);b=adc.estimate(8,130,160e6,64)
    assert b['frequency_per_adc_hz']==2.5e6
    assert math.isclose(a['energy_per_conversion_pj'],b['energy_per_conversion_pj'],rel_tol=1e-12)
    assert b['total_area_um2']==64*b['area_per_adc_um2']


def test_hard_frequency_cap_and_resolution(adc):
    with pytest.raises(AssertionError):adc.estimate(8,130,20e9,1)
    assert adc.estimate(8,130,20e9,2)['within_model_hard_frequency_limit']
    with pytest.raises(AssertionError):adc.estimate(3,130,1e6,1)
    with pytest.raises(ValueError):adc.estimate(8,130,1e6,0)


def test_same_fit_successor(adc):
    modern=ADCModel(True)
    assert adc.fit==modern.fit
    assert adc.estimate(7,130,160e6,4)==modern.estimate(7,130,160e6,4)


def test_skipping_write_cannot_remove_reference_driver_or_detection():
    measured=dict(reference=100,buffer=100,comparator=100,additional_comparison_fj=100,
        control_pin_charging_upper=10,interface=50,writer=5)
    full=boundary_cost(measured,3);skip=boundary_cost(measured,3,regeneration_fraction=0)
    assert full-skip==55
    assert skip==510


def test_precision_is_exact_tail_bound_not_positions_log_beta():
    row,traces=precision_layout(128,6)
    assert row['guaranteed_truncated_positions']==5
    assert row['exact_tail_bound_pass'] and row['exact_finite_word_checks_pass']
    assert row['empirical_rmse_units']<=row['ideal_adc_quantization_rmse_units']
    assert len(traces)>500
    # Actual exact algebraic residuals, including noninteger irrational tails.
    assert any(F(r['exact_error_B'])!=0 for r in traces)


def test_oversampling_clock_is_not_an_independent_conversion(adc):
    anchor=dict(enob_from_sndr=12,osr=100,area_um2=1,energy_per_sample_pj=.001,
        independent_value_rate_bound_hz=1e12,anchor_id='invalid_internal_clock_comparison')
    r=conventional_choice(adc,8,160e6,[anchor],'optimistic_anchor_envelope')
    assert r['adc_source']=='pinned_statistical_model'


def test_useful_mac_amortization_and_exact_bitplane_contract():
    rows=list(csv.DictReader((OUT/'adc_cim_amortization.csv').open()))
    row=next(r for r in rows if r['scenario']=='one_shot_quantized_dot' and r['accumulation_capacity']=='128' and r['utilization']=='1')
    assert math.isclose(float(row['adc_energy_per_scalar_mac_fj']),float(row['adc_cost_per_output_pj'])*1000/128)
    exact=next(r for r in rows if r['scenario']=='exact_int8_64_bitplanes' and r['accumulation_capacity']=='128' and r['enob']=='8')
    assert int(exact['adc_conversions_per_output'])==64
    assert any(r.get('status','').startswith('REJECTED') for r in rows)


def test_exported_extrapolation_is_retained_not_frontier():
    rows=list(csv.DictReader((OUT/'adc_130nm_pareto.csv').open()))
    rejected=[r for r in rows if r['within_model_hard_frequency_limit']=='False']
    assert len(rejected)==7
    assert all(r['on_energy_area_frontier']=='False' and 'REJECTED' in r['extrapolation_warning'] for r in rejected)


def test_training_anchors_keep_sample_and_nyquist_work_distinct():
    rows=list(csv.DictReader((OUT/'adc_130nm_survey_anchors.csv').open()))
    assert len(rows)==12
    row=next(r for r in rows if float(r['osr'])>1)
    assert float(row['survey_energy_per_nyquist_pj'])>float(row['energy_per_sample_pj'])
