import math
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from gamma_normalizer_statistics import confidence,charge_monte_carlo
from gamma_normalizer_energy import break_even,summarize
from gamma_spice import DEFAULT_EXE,read_ascii_raw,classify


def test_zero_failure_confidence_never_means_zero_probability():
    assert confidence(0,100)['one_sided_95pct_upper'] > 0.02
    assert confidence(0,3000000)['one_sided_95pct_upper'] < 1e-6
    assert confidence(10,100)['one_sided_95pct_lower'] < .1 < confidence(10,100)['one_sided_95pct_upper']
    with pytest.raises(ValueError):confidence(1,0)


def test_sampled_scaling_model_reports_noise_failures_without_calling_spice():
    safe=charge_monte_carlo(3,4,1000,combined=True)
    unsafe=charge_monte_carlo(3,256,1000,combined=True)
    assert safe['failures']==0
    assert unsafe['estimate']>.1
    assert 'NOT SPICE' in safe['source']
    assert unsafe['one_sided_95pct_upper']>unsafe['estimate']


def test_break_even_accounts_for_baseline_adc_amortization_and_extra_array_cost():
    s={'q':3,'scheme':'unit','width_mean':5,
       'comparator_evaluations_until_settled_plus_confirmation_mean':20,
       'final_readout_comparators_mean':20,'unit_rewrites_mean':10,'trigger_events_mean':10}
    r=break_even(s,1,1,control_per_trigger_fj=0,binary_accum_fj=0)
    assert r['gamma_extra_fj']==50
    amortized=break_even(s,1,1,control_per_trigger_fj=0,binary_accum_fj=0,adc_conversions=1/16)
    assert amortized['adc_break_even_fj_per_conversion']==16*r['adc_break_even_fj_per_conversion']
    extra=break_even(s,1,1,control_per_trigger_fj=0,binary_accum_fj=0,delta_cim_fj=100)
    assert extra['gamma_extra_fj']==150
    assert break_even(s,schedule='fixed_budget')['detect_comparator_evaluations']>r['detect_comparator_evaluations']


def test_detector_counts_use_compiled_layout_not_observed_answer_span():
    from gamma_protocols import workload_study
    rows=summarize(workload_study(seed=11,count=4))
    expected={3:21,4:17,5:17}
    for row in rows:
        width=expected[row['q']]
        assert row['compiled_layout_width']==width
        assert row['width_mean']==width
        assert row['final_readout_comparators_mean']==width*2*(row['q']-1)


@pytest.mark.parametrize('q',[3,4,5])
@pytest.mark.skipif(not DEFAULT_EXE.exists(),reason='LTspice unavailable')
def test_reusable_behavioral_subcircuit_is_executable_and_signed(q,tmp_path):
    source=Path(__file__).resolve().parents[1]/'spice'/'behavioral'
    shutil.copyfile(source/'gamma_unit_cell.sub',tmp_path/'gamma_unit_cell.sub')
    deck=tmp_path/'demo.cir'
    for sign in (1,-1):
        text=(source/f'q{q}_unit_demo.cir').read_text().replace('SIGN=1',f'SIGN={sign}')
        deck.write_text(text)
        subprocess.run([str(DEFAULT_EXE),'-b','-ascii',str(deck)],check=True,timeout=30,
            creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess,'CREATE_NO_WINDOW') else 0)
        names,data=read_ascii_raw(deck.with_suffix('.raw'))
        result={n:classify(float(data[-1,i]),.025,q-1) for i,n in enumerate(names) if n!='time'}
        assert result=={'v(sm1)':sign,'v(bl)':0,'v(s0)':0,'v(sp1)':sign,'v(br)':0}
