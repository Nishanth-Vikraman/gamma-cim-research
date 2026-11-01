"""Evidence-linked SKY130 report and eight required engineering plots."""
from pathlib import Path
import json,math,statistics,hashlib,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gamma_sky130 import ROOT,BASE,dump
from run_sky130_phase import OUT,write_csv


def table(headers,rows):
    def clean(v):return str(v).replace('|','/').replace('\n',' ')
    return '\n'.join(['| '+' | '.join(map(clean,headers))+' |','| '+' | '.join(['---']*len(headers))+' |']+
        ['| '+' | '.join(map(clean,row))+' |' for row in rows])


def physical_energy(row):
    record=json.loads((ROOT/row['run']/'metrics.json').read_text())
    rails=('reference','interface','writer','buffer','comparator')
    obs=record['observations']
    net=statistics.mean(sum(p['signed_port_work_fj'].get(k,0) for k in rails) for p in obs)
    delivered=statistics.mean(sum(p['energy_components_fj'].get(k,0) for k in rails) for p in obs)
    allowance=statistics.mean(p['energy_components_fj'].get('control_pin_charging_upper',0) for p in obs)
    return record,{'vdd_net_fj':net,'vdd_delivered_fj':delivered,'external_pin_allowance_fj':allowance,
        'known_path_conservative_fj':delivered+allowance}


def failure_attribution(s):
    out=[]
    for group in ('pvt_energy','pvt_margin','mc_pdk','mc_architecture','mc_margin'):
        rows=[r for r in s['followup'] if r['group']==group];obs=[]
        for r in rows:obs.extend(json.loads((ROOT/r['run']/'metrics.json').read_text())['observations'])
        bad=[o for o in obs if not o['next_mac_correct']]
        out.append({'group':group,'observed_checkpoints':len(obs),'wrong_mac':len(bad),
            'wrong_selected_class':sum(round(o['selected_class'])!=o['digit'] for o in obs),
            'wrong_mac_with_correct_selection':sum(round(o['selected_class'])==o['digit'] for o in bad),
            'wrong_mac_with_correct_selection_and_storage_within_point1digit':sum(round(o['selected_class'])==o['digit'] and abs(o['stored_error_units'])<.1 for o in bad),
            'median_abs_storage_residual':statistics.median(abs(o['stored_error_units']) for o in obs),
            'median_abs_driven_residual':statistics.median(abs(o['driver_error_units']) for o in obs),
            'rms_driver_minus_stored_units':statistics.mean((o['driver_error_units']-o['stored_error_units'])**2 for o in obs)**.5,
            'whole_networks_correct_and_inside_monitored_limits':sum(r['wrong_next_mac']==0 and r['unresolved_checkpoints']==0 and r['soa_violations']==0 for r in rows)})
    dump(OUT/'failure_attribution.json',out)


def plots(s):
    sizes=s['sizing'];follow=s['followup'];pvt=[r for r in follow if r['group'].startswith('pvt_')]
    good=[r for r in sizes if r['feasible_nominal']];bad=[r for r in sizes if not r['feasible_nominal']]
    def save(fig,name):fig.tight_layout();fig.savefig(OUT/name,dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5))
    for rows,label,color in ((bad,'Fails margin/class/operating limit','gray'),(good,'Meets nominal margin/class limits','tab:green')):
        ax.scatter([r['cap_ff'] for r in rows],[r['mean_energy_fj'] for r in rows],label=label,c=color,s=26,alpha=.65)
    for budget in (1,2.74,5,10,20):ax.axhline(budget,ls=':',lw=.8);ax.text(.98,budget,f'{budget:g} fJ',transform=ax.get_yaxis_transform(),ha='right',va='bottom',fontsize=8)
    ax.set(yscale='log',xlabel='Nominal storage capacitance (fF)',ylabel='Known-path cost (fJ / position)',title='Storage/sizing tradeoff; detector ideal in this isolation sweep')
    ax.legend(fontsize=8);ax.grid(alpha=.2);save(fig,'energy_vs_storage.png')
    fig,ax=plt.subplots(figsize=(8,5))
    ax.scatter([r['max_abs_error_units'] for r in sizes],[r['mean_energy_fj'] for r in sizes],c=['tab:green' if r['feasible_nominal'] else 'gray' for r in sizes],s=26)
    frontier=[r for r in sizes if r['pareto']];frontier.sort(key=lambda r:r['max_abs_error_units'])
    ax.plot([r['max_abs_error_units'] for r in frontier],[r['mean_energy_fj'] for r in frontier],'o-',color='tab:red',label='Energy/error/settling nondominated points')
    ax.axvline(.25,color='black',ls='--',label='Following gain2 limit before other errors')
    ax.set(xscale='log',yscale='log',xlabel='Worst absolute canonical error (digit units)',ylabel='Known-path cost (fJ / position)',title='Finite reference writes and actual SKY130 output drivers')
    ax.legend(fontsize=8);ax.grid(alpha=.2);save(fig,'energy_vs_error.png')
    fig,ax=plt.subplots(figsize=(8,5))
    ready=[r for r in sizes if r['max_settling_ns']<1e8]
    pts=ax.scatter([r['buffer_scale'] for r in ready],[r['max_settling_ns'] for r in ready],c=[r['cap_ff'] for r in ready],cmap='viridis',alpha=.7)
    ax.set(xlabel='Output-stage width scale',ylabel='Worst settling to 0.1 digit (ns)',title='Points that reach the declared residual target')
    fig.colorbar(pts,ax=ax,label='Storage fF');ax.grid(alpha=.2);save(fig,'settling_vs_buffer.png')
    allobs=[]
    for r in follow:
        if r['mode']=='regenerated':allobs.extend(json.loads((ROOT/r['run']/'metrics.json').read_text())['observations'])
    error=np.array([abs(r['driver_error_units']) for r in allobs]);correct=np.array([r['next_mac_correct'] for r in allobs])
    bins=np.array([0,.025,.05,.1,.15,.20,.25,.30,.4,.5,.75,1,2,4,20])
    xx=[];yy=[];nn=[]
    for lo,hi in zip(bins[:-1],bins[1:]):
        m=(error>=lo)&(error<hi)
        if np.any(m):xx.append((lo+hi)/2);yy.append(float(correct[m].mean()));nn.append(int(m.sum()))
    fig,ax=plt.subplots(figsize=(8,5));ax.plot(xx,yy,'o-');ax.axvline(.25,color='red',ls='--')
    for x,y,n in zip(xx,yy,nn):ax.annotate(str(n),(x,y),xytext=(2,6),textcoords='offset points',fontsize=7)
    ax.set(xscale='symlog',xlim=(0,max(1,max(xx))),ylim=(-.05,1.1),xlabel='Absolute regenerated residual (digits)',ylabel='Observed next-MAC correct fraction',title='Mixed PVT/mismatch conditions; labels show correlated checkpoint counts')
    ax.grid(alpha=.2);save(fig,'next_mac_vs_residual.png')
    nom=[r for r in pvt if r['corner']=='tt' and r['vdd']==1.8 and r['temp_c']==25]
    components=('reference','interface','writer','buffer','comparator','control_pin_charging_upper')
    fig,ax=plt.subplots(figsize=(8,5));bottom=np.zeros(len(nom));labels=[r['group'].replace('pvt_','') for r in nom]
    for key in components:
        values=np.array([r.get('energy_'+key+'_fj',0) for r in nom]);bars=ax.bar(labels,values,bottom=bottom,label=key.replace('_',' '))
        if key=='control_pin_charging_upper':
            for bar in bars:bar.set_hatch('//')
        bottom+=values
    ax.set(ylabel='fJ / regenerated position',title='Actual VDD delivered energy + external pin charging allowance\nOne real threshold; remaining detector/control costs excluded')
    ax.legend(fontsize=8);save(fig,'energy_breakdown.png')
    fig,axes=plt.subplots(2,1,figsize=(9,7))
    corners=('tt','ss','ff','sf','fs')
    for j,candidate in enumerate(('pvt_energy','pvt_margin')):
        rows=[r for r in pvt if r['group']==candidate and r['simulation_status']=='complete'];pos=np.arange(5)+(j-.5)*.25
        for ax,key in zip(axes,('mean_energy_fj','max_abs_error_units')):
            values=[[r[key] for r in rows if r['corner']==corner] for corner in corners]
            ax.boxplot(values,positions=pos,widths=.2,patch_artist=True,boxprops={'facecolor':('lightblue','salmon')[j]},manage_ticks=False)
            ax.plot([],[],color=('tab:blue','tab:red')[j],label=candidate.replace('pvt_',''))
    axes[0].set(yscale='log',ylabel='Known-path fJ',title='Sampled PVT spread: 3 supplies x3 temperatures/corner; not yield')
    axes[1].set(ylabel='Worst error (digits)');axes[1].axhline(.25,ls='--',c='black')
    for ax in axes:ax.set_xticks(np.arange(5),corners);ax.legend();ax.grid(alpha=.2)
    save(fig,'pvt_spread.png')
    fig,ax=plt.subplots(figsize=(8,5))
    for group in ('mc_pdk','mc_architecture','mc_margin'):
        values=[r['max_abs_error_units'] for r in follow if r['group']==group]
        ax.hist(values,bins=12,histtype='step',linewidth=2,label=f'{group}: n={len(values)} networks')
    ax.axvline(.25,color='black',ls='--');ax.set(xlabel='Maximum absolute canonical residual / network (digits)',ylabel='Independent network trials',title='Seeded mismatch brackets; transient thermal/flicker noise excluded')
    ax.legend(fontsize=8);save(fig,'monte_carlo_residual.png')
    # Replace the old hypothetical persistent readout instead of double pricing it.
    legacy=json.loads((ROOT/'spice/regen_results/phase_summary.json').read_text())
    nonboundary=legacy['energy']['reference_point']['explicit event model']['before_regen_extra_fj']-17*8*.1
    energy_nom=next(r for r in nom if r['group']=='pvt_energy')
    _,measured=physical_energy(energy_nom)
    comp=s['selected_comparator']['max_legal_center_energy_fj']
    modeled3=measured['known_path_conservative_fj']+2*comp
    modeled6=measured['known_path_conservative_fj']+5*comp
    threshold3=nonboundary+17*modeled3-25;threshold6=nonboundary+17*modeled6-25
    adc=np.logspace(1,6,150);fig,ax=plt.subplots(figsize=(9,5))
    for conversions in (1/64,1/16,1,4,16):
        ax.plot(adc,conversions*adc-threshold3,label=f'{conversions:g} ADC conversion(s)/product')
    ax.axhline(0,color='black');ax.axvline(100,ls=':',c='red',label='100 fJ illustrative ADC')
    ax.set(xscale='log',yscale='symlog',xlabel='Assumed ADC energy (fJ / conversion)',ylabel='Baseline minus gamma cost (fJ / product)',title='Optimistic 3-decision estimate; positive means energy room\nMissing controller/reference selection costs = 0; delta array/routing = 0')
    ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.2);save(fig,'adc_break_even.png')
    breakrows=[]
    for decisions,cost in ((1,measured['known_path_conservative_fj']),(3,modeled3),(6,modeled6)):
        for logic in (0,5,20,100):
            for delta in (0,100,1000):
                for conversions in (1/64,1/16,1,4,16):
                    breakrows.append({'physical_decisions_estimated':decisions,'known_boundary_plus_decisions_fj':cost,
                        'extra_latch_control_reference_selection_per_position_fj':logic,'extra_array_routing_fj':delta,
                        'adc_conversions_per_product':conversions,'adc_break_even_fj':(nonboundary+17*(cost+logic)+delta-25)/conversions})
    write_csv(OUT/'system_break_even.csv',breakrows)
    return {'nominal_hybrid_energy':measured,'nonboundary_optimistic_fj':nonboundary,'extra_decision_fj':comp,
        'three_decision_boundary_lower_estimate_fj':modeled3,'six_decision_boundary_lower_estimate_fj':modeled6,
        'adc_break_even_three_decisions_fj_one_conversion':threshold3,'adc_break_even_six_decisions_fj_one_conversion':threshold6,
        'known_path_only_adc_break_even_fj':nonboundary+17*measured['vdd_net_fj']-25}


def report(s,energy):
    sizes=s['sizing'];follow=s['followup'];pvt=[r for r in follow if r['group'].startswith('pvt_')]
    manifest=json.loads((BASE/'pdk_manifest.json').read_text());validation=json.loads((OUT/'model_validation.json').read_text())
    selected=json.loads((OUT/'selected_design.json').read_text());d=selected['design']
    nominal=next(r for r in pvt if r['group']=='pvt_energy' and r['corner']=='tt' and r['vdd']==1.8 and r['temp_c']==25)
    record,physical=physical_energy(nominal)
    gate=table(['Budget fJ/position','Measured VDD / budget','Known path + pin allowance / budget'],
        [(g,f"{physical['vdd_net_fj']/g:.1f}x",f"{physical['known_path_conservative_fj']/g:.1f}x") for g in (1,2.74,5,10,20)])
    pvttable=table(['Candidate','PVT points','All12 MAC results correct','Incomplete transients','Full runs inside monitored device limits','Passes both class and monitored limits','Worst observed residual (digits)'],
        [(label,len(rows),sum(r['wrong_next_mac']==0 and r['unresolved_checkpoints']==0 for r in rows),sum(r['unresolved_checkpoints']>0 for r in rows),sum(r['soa_violations']==0 and r['unresolved_checkpoints']==0 for r in rows),sum(r['wrong_next_mac']==0 and r['soa_violations']==0 and r['unresolved_checkpoints']==0 for r in rows),f"{max(r['max_abs_error_units'] for r in rows):.3g}")
            for label in ('energy','margin') for rows in [[r for r in pvt if r['group']=='pvt_'+label]]])
    mctable=table(['Network bracket','Networks failing acceptance','Observed checkpoints wrong','Unresolved checkpoints','95% one-sided upper error bound'],
        [(r['mode'],f"{r['networks_with_any_wrong_class']}/{r['network_trials']}",f"{r['wrong_checkpoints']}/{r['checkpoints']}",r['unresolved_checkpoints'],f"{r['one_sided_95pct_upper']:.1%}") for r in s['monte_carlo']])
    chain=[r for r in follow if r['group']=='chaining']
    chaintable=table(['Storage','Initial residual (digits)','Regenerated correct','Direct correct'],
        [(kind,res,f"{12-next(r for r in chain if r['cap_kind']==kind and r['residual_input_units']==res and r['mode']=='regenerated')['wrong_next_mac']}/12",
            f"{12-next(r for r in chain if r['cap_kind']==kind and r['residual_input_units']==res and r['mode']=='direct')['wrong_next_mac']}/12") for kind in ('ideal','mim') for res in (0,.1,.35)])
    matched=[r for r in follow if r['group']=='chaining_matched_scale']
    matchedtable=table(['Common-scale initial residual (digits)','Regenerated correct','Direct correct'],
        [(res,f"{12-next(r for r in matched if r['residual_input_units']==res and r['mode']=='regenerated')['wrong_next_mac']}/12",
          f"{12-next(r for r in matched if r['residual_input_units']==res and r['mode']=='direct')['wrong_next_mac']}/12") for res in (0,.1,.35)])
    transitions=list(__import__('csv').DictReader((OUT/'regen_transition_matrix.csv').open()))
    trE=[float(r['energy_fj']) for r in transitions]
    attribution=json.loads((OUT/'failure_attribution.json').read_text())
    failuretable=table(['Bracket','Wrong next-MAC outputs','Wrong outputs despite correct selected class','Wrong with correct selection and storage error <0.1 digit','RMS driver-minus-storage error (digits)'],
        [(r['group'],r['wrong_mac'],r['wrong_mac_with_correct_selection'],r['wrong_mac_with_correct_selection_and_storage_within_point1digit'],f"{r['rms_driver_minus_stored_units']:.3g}") for r in attribution])
    comptable=table(['W (um)','L (um)','Corner','Wrong legal-center decisions /2','Mean energy incl. clock allowance (fJ)','Max decision delay (ns)','Offset bracket (mV)'],
        [(r['comparator_w'],r['comparator_l'],r['corner'],r['legal_center_wrong'],f"{r['mean_energy_fj']:.1f}",f"{r['max_decision_ns']:.3g}",
          ('unresolved/nonmonotonic' if r.get('offset_bracket_low_v') is None else f"[{1000*r['offset_bracket_low_v']:.3g}, {1000*r['offset_bracket_high_v']:.3g}]"))
         for r in s['comparators'] if (r['comparator_w'],r['comparator_l']) in ((s['hybrid_design']['comparator_w'],s['hybrid_design']['comparator_l']),(16,.5))])
    regression=json.loads((OUT/'regression_result.json').read_text())
    compmc=[]
    for w,l in sorted(set((r['comparator_w'],r['comparator_l']) for r in s['comparator_mc'])):
        rows=[r for r in s['comparator_mc'] if r['comparator_w']==w and r['comparator_l']==l]
        compmc.append((w,l,len(rows),sum(r['legal_center_wrong']>0 for r in rows),sum(r['legal_center_wrong'] for r in rows)))
    compmctable=table(['W um','L um','Independent seeds','Trials with wrong legal-center decision','Wrong decisions /2n'],compmc)
    source_urls=['https://github.com/chrisklaus83188/sky130-ngspice-models',manifest['voltage_limit_source_url'],manifest['ngspice_download']]
    text=f'''# q=4 SKY130 regeneration: engineering result

Run date: 2026-10-24. Simulator: ngspice47, regular SKY130A 1.8 V NMOS/PMOS; no full array or digital ASIC flow.

## Decision

**The tested voltage-storage / OTA-driver / resistor-MAC implementation does not establish an energy advantage. Redesign this boundary before a crossbar prototype.** The transistor writer and loaded driver can work, but their measured cost alone is far above the aggressive reference budget. The complete seven-state detector is not finished: one important threshold is transistorized and five remain behavioral. Missing decoder, reference selection and completion costs can only tighten the disclosed budget under the stated no-recovery assumptions.

The structured isolation sweep contains {len(sizes)} designs, with {s['nominal_feasible_cases']} meeting the declared nominal MAC/error/operating checks. Its minimum known-path point costs {s['selected_measured_known_path_fj']:.1f} fJ/position before a real region decision is included. It uses {d['cap_ff']:g} fF storage, {d['tg_w']:g} um NMOS selector W, {d['period_ns']:g} ns periods, {d['ladder_r']/1e6:g} Mohm total reference ladder, {d['bias_r']/1e6:g} Mohm bias resistor and {d['load_r']/1e6:g} Mohm MAC resistor scale. Smaller storage is an observed tradeoff, not an assumption that100 fF was optimal.

At tt/1.8 V/25 C, the selected hybrid circuit measures **{physical['vdd_net_fj']:.1f} fJ of signed1.8 V supply work per position**. Gross delivered supply energy is {physical['vdd_delivered_fj']:.1f} fJ; a conservative external control-pin charging allowance adds {physical['external_pin_allowance_fj']:.1f} fJ, giving {physical['known_path_conservative_fj']:.1f} fJ for the implemented blocks. These are per-operation means over the declared12-transition workload. They are not claimed as fully measured seven-state detector/control-system energy.

{gate}

This triggers the requested cost gate for these tested circuits. It is not proof that every possible Gamma-CIM physical representation is impossible. The quantitative limitation is this loaded voltage interface, together with the measured reference and driver cost; the algebraic datatype remains exact.

## Electrical path and architecture integrity

The native lead remains q=4, beta=2+sqrt(3), digits -3..3. Position significance is carried by the index, never by exponentially increasing voltage. No q sweep, high-voltage devices, general ADC/DAC, RRAM, FPGA or full CIM array was introduced.

Nominal canonical levels are VCM=0.6 V plus d*25 mV. The resistor ladder has seven integer-state taps plus six half-step taps, from the1.8 V supply. Stage1 stimulus is an imperfect analog state at the previous word path's raw scale S=12 (8.333 mV per coefficient); the restored interface uses S=4 (25 mV). All positions share the same scale policy. Inputs are not vector-normalized. Supply sweeps track the declared common reference scale; an additional untracked reference drift would need calibration/error budget.

The physical path is ladder -> seven MOS switches -> storage -> SKY130 closed-loop output OTA -> two resistors and20 fF next-MAC load. The passive sum is (2*Vdrive+Vbias)/3, with external operand Vbias=VCM-d*Delta. Declared scale-aware output bins have spacing Delta/3, so the recovered coefficient is2*d_x-d. The external bias is fixed from the original input workload, not a Python re-encoding of an intermediate result. The actual stored/transistor-driven node supplies the next MAC.

The direct control uses an actual OTA with gain3 feedback to compensate the known raw/interface scale while preserving the analog residual. It does not receive an expected canonical answer. The isolation sweep also compares an additional shared write-buffer topology, two-stage and five-transistor-core OTAs, continuous and gated operation, transmission gates and NMOS-only switches. The shared write-buffer closes its storage aperture before disconnecting its reference input, avoiding amplified mux charge injection.

Ideal one-hot decisions drive actual transistor inverter interfaces. A hybrid case replaces the +0.5 raw-digit threshold with a complementary StrongARM comparator; the other five decisions and sequencing remain ideal. The -3..3 debug label is not a physical +/-3 V logic rail. Every MOS connection uses the declared supply range and is monitored for external gate/drain/source and source/body limits. Invalid operating points are retained and excluded from nominal feasibility.

The input contract is a previously normalized legal digit. A complete input overflow/invalid-state controller is not transistorized. Final MAC readout explicitly flags out-of-alphabet regions and does not repair them in software. External analog-input and MAC-bias port work is listed separately; including it cannot be assumed free in a full-system comparison. The reported per-boundary cost excludes it only as an explicitly separate MAC operand cost, not as an energy-recovery credit.

## PDK, simulator and Stage0

Pinned repository: `{manifest['git_commit']}`; open_pdks provenance hash: `{manifest['open_pdks_build_hash']}`. The extraction's upstream provenance is recorded, not independently rebuilt. Apache-2.0 license and NOTICE are retained. Model data is {manifest['model_file_bytes']/1e6:.1f} MB excluding Git metadata. Only libs.tech/ngspice and libs.ref/sky130_fd_pr/spice were installed, plus the portable simulator runtime. No OpenLane, Magic, KLayout, standard-cell, I/O or SRAM library was installed.

All20 NMOS/PMOS/inverter/transmission-gate smoke checks pass across tt,ss,ff,sf,fs. The typical100 fF transmission-gate test exposes about2.8 mV of combined hold shift: correct switching is not an ideal canonical write. DC Ron/leakage and reference-bank curves are retained in model_validation.json. The detector and output-driver energy is actual supply work; behavioral-source power is never substituted for measured transistor comparator energy.

The independent diagnostics establish: identical devices match with mismatch disabled; different instances and seeds differ when MC_MM_SWITCH is enabled; the same seed replays identically; batch temperature control agrees with a literal .temp standalone run. MC_PR_SWITCH remains0 so discrete process corners and local statistical mismatch are distinct. The PDK is unmodified and Git-clean. Seeds are set before circuit reset/model evaluation.

Core-model documentation gives1.95 V VDS/VGS limits. The proposed1.98 V point was therefore excluded; the PVT sweep uses1.62,1.8,1.9 V and-40,25,100 C. This is model-level validation, not foundry signoff. The five FET corners use the library's typical R/C corner; capacitor extremes and ideal-resistor mismatch are separately bracketed.

## PVT, mismatch and chaining

Both the minimum-energy and a larger-margin Pareto candidate are carried through PVT. Preferred MIM storage is used in PVT; ideal/MIM pairs isolate storage-model effects at nominal conditions. The MIM geometry comes from the supplied area/perimeter model, with large reference decoupling implemented as parallel units of at most1 pF each. This includes the model's series parasitics, not physical-layout extraction.

{pvttable}

Incomplete transients are retained as numerical nonconvergence, never as passes. The failing logs name ideal B-selector branches, so these cases do not by themselves prove physical transistor oscillation. Their unobserved checkpoints remain unresolved; only complete cycles contribute energy statistics. The PVT spread plot includes completed simulations only. The larger-device candidate was chosen for lower residual in the ideal-detector sizing sweep; it is not a proven robust design after transistor comparator loading is attached.

{mctable}

Trials are independent seeded whole circuits. Twelve checkpoints inside each circuit are correlated and are not twelve independent yield samples. PDK mismatch uses the provided local FET/MIM machinery. Architecture-only trials disable it and add ideal-resistor mismatch, explicit capacitance error and an input-referred systematic comparator offset. The two mechanisms are not silently double-counted. These are uncalibrated brackets; deterministic transient analysis does not automatically include device thermal/flicker noise. No1e-3/1e-6 reliability certification follows from these sample counts.

{chaintable}

The S12-to-S4 direct branch already fails at zero injected residual: its extra gain3 feedback divider loads the finite common-mode ladder and its amplifier has finite gain/offset. This table therefore cannot independently prove that regeneration alone prevents accumulated error. The following control holds input/output scale atS4, uses identical unity-gain physical drivers and MIM storage, and isolates the transistor write path with ideal latched seven-state decisions. No scale-conversion amplifier or expected-answer source is inserted.

{matchedtable}

## Comparator and primitive characterization

The isolated comparator sweep spans W=1/4/16 um and L=0.15/0.5 um at all five corners. Its declared legal-center overdrive is +/-4.167 mV at nominal supply. Offset is a sampled transition bracket, not a precision offset measurement or calibrated yield claim; unresolved/nonmonotonic transfers do not receive a fictitious offset. Delay is measured from the clock's 50% crossing. At exactly zero differential input either resolved polarity is permitted; unresolved nonzero cases are retained.

{comptable}

The table includes actual comparator VDD work plus the external clock-pin charging allowance. Separate fields retain the two components. A larger comparator does not guarantee a better loaded boundary: input capacitance, reference kickback and decision energy grow too. The wide L=0.5 um example uses about53.3 pF of physical threshold decoupling versus1 pF for the smallest input pair. Those capacitors are present in the circuit, and their prebiased cold-start energy is excluded explicitly.

{compmctable}

The switch diagnostics retain Ron versus level and width, leakage, and loaded reference-bank DC curves. Canonical residual, write hold drift, within-cycle peak reference droop, settling and each component's energy are in the per-case metrics; storage ideal/MIM comparisons use identical external timing and loads. Clock feedthrough and charge injection arise from the MOS models and declared estimated junction geometry. Supply variation and the declared resistor/capacitor/comparator perturbations are included; layout-dependent parasitic uncertainty and transient sample/thermal noise remain unmeasured.

For the selected nominal MIM path, storage error reaches0.208 digit, driver error0.164 digit, and hold drift0.0106 digit. It passes the next-MAC class test but does not settle inside0.1 digit on every cycle. In the CSV, max_settling_ns=1e9 is an optimization sentinel for failure to reach that target during observation, not a measured one-second settling time; cycles_not_settled_to_point1digit and max_observed_settling_ns make censoring explicit. Ideal/MIM nominal path errors are nearly identical in these unextracted typical-model comparisons, so nominal MIM parasitics are not the dominant observed limitation.

All49 previous/next digit pairs are also exercised through actual preceding writes; {sum(r['next_mac_correct']=='True' for r in transitions)}/49 return the correct next-MAC class. Their mean known-path cost is {statistics.mean(trE):.1f} fJ; worst case is {max(trE):.1f} fJ. The matrix includes no-op pairs that still execute reset/write; sparse skip control is not implemented or assumed free. Per-trial records retain seed, corner, every MOS W/L, reference levels, stored/driven voltages, final class and component energy.

## Why the budget is difficult

The ladder tradeoff is real: static energy scales as VDD squared times period divided by R, while write settling grows with R times storage-plus-gate capacitance. Increasing R without changing the electrical write aperture is not a free saving. The sweep includes coupled R/C/timing choices rather than shrinking storage alone. Reference-bank sharing has not been electrically demonstrated. Even subtracting the entire nominal reference-bank cost leaves {physical['known_path_conservative_fj']-nominal['energy_reference_fj']:.1f} fJ/position for the implemented remaining blocks: an optimistic free-reference bound, not a sharing result.

The loaded output OTA also spends bias/drive current for gain, settling and resistor-load current. Clock gating can lose its nominal saving through wake-up recovery, so continuous low-bias operation was tested too. Returning a voltage inside its own +/-0.5-digit class is insufficient: the following weight2 needs residual comfortably below0.25 digit before MAC and sensing errors. PVT and mismatch results retain this stricter requirement.

{failuretable}

This attribution uses actual selected, stored and driven checkpoints, not an offline correction. Errors can overlap or compensate, so the columns are diagnostics rather than disjoint causal probabilities. In the larger-device PDK bracket,87 wrong outputs occur even with correct class selection and storage error below0.1 digit. The measured driver-minus-storage residual therefore identifies loaded analog buffering/matching as a concrete unresolved requirement. In the energy candidate, small-capacitor write/feedthrough error also consumes much of the next-MAC margin. Simply improving the class detector does not fix either effect.

Input-device upsizing improves a comparator's matching bracket but increases reference kickback. The physical threshold interpolation and its MIM decoupling are priced in the hybrid path; decoupling scales with input gate area. This can require much more capacitance than the storage itself. Measurements use a prebiased steady-state reference bank. Cold-start charging is not included, which favors Gamma; it must be amortized explicitly before claiming a complete operating energy.

## System break-even

The new calculation replaces the old hypothetical persistent readout rather than counting it twice. The remaining optimistic normalization/control allowance is {energy['nonboundary_optimistic_fj']:.2f} fJ/product on the frozen17-position workload. The one-threshold hybrid known cost is measured;3- or6-decision scenarios add actual comparator-primitive energy, with unbuilt threshold selection/latch/completion cost swept separately. A3-comparison scheme assumes an unbuilt adaptive shared-reference controller, not demonstrated hardware.

With zero extra array/routing/control cost, one ADC conversion per product and the old25 fJ binary-accumulation+requantization allowance, the optimistic3-decision ADC break-even is **{energy['adc_break_even_three_decisions_fj_one_conversion']/1000:.2f} pJ/conversion**. The6-decision estimate is {energy['adc_break_even_six_decisions_fj_one_conversion']/1000:.2f} pJ. Strong ADC amortization raises these thresholds. A mathematical crossing at a sufficiently expensive ADC is not a measured, matched-accuracy energy advantage. The existing100 fJ ADC reference point is not competitive for this cell.

## Required plots

![Energy versus storage](results/energy_vs_storage.png)
![Energy versus residual](results/energy_vs_error.png)
![Settling versus buffer](results/settling_vs_buffer.png)
![Following MAC versus residual](results/next_mac_vs_residual.png)
![Energy breakdown](results/energy_breakdown.png)
![PVT spread](results/pvt_spread.png)
![Monte Carlo residual](results/monte_carlo_residual.png)
![ADC break-even](results/adc_break_even.png)

## Reproduce and inspect

From python/: `python -B gamma_sky130.py`; `python -B run_sky130_phase.py --workers 2 --mc-trials 16`; `python -B validate_sky130.py`; `python -B report_sky130_phase.py`.

Every case retains a standalone circuit.cir, seeded replay.cir, golden metadata, simulator log, parameters and metrics. Full waveforms are kept for representative cases with lossless gzip; all checkpoint voltages, energy integrals and classes are retained for every reported case. Unreported incomplete old traces were archived losslessly to control disk use. Binary output and bounded batches reduce overhead; neither changes the PDK. Local decoded gamma states, not offline numerical reconstruction, determine the result.

Software regression: `{regression['output'].strip().splitlines()[-1]}`. Independent literal-deck replay of the selected loaded path is recorded in results/standalone_path_audit.json, alongside the model random-seed and temperature diagnostics.

Required results are under results/: sky130_phase_summary.json, energy_breakdown.csv, pvt_summary.csv, monte_carlo_summary.csv, sizing_sweep.csv and regen_transition_matrix.csv. The previous LTspice reference remains unchanged. Do not build an array or FPGA from these results. The next decision is a redesigned shared/passive/charge-oriented reference-and-driver primitive, with real fanout, load and control accounting; a datatype novelty claim does not remove these electrical costs.

Sources: [extracted models]({source_urls[0]}), [device voltage limits]({source_urls[1]}), [portable simulator]({source_urls[2]}). Local licenses, source snapshots and hashes are recorded in pdk_manifest.json.
'''
    for old,new in {'All20':'All 20','All49':'All 49','the1.8':'the 1.8','declared12':'declared 12','Stage1':'Stage 1','stage1':'stage 1','weight2':'weight 2','about53.3':'about 53.3','versus1 pF':'versus 1 pF','The6-decision':'The 6-decision','A3-comparison':'A 3-comparison','the3-decision':'the 3-decision','with one threshold and400':'with one threshold and 400','not twelve':'not twelve'}.items():text=text.replace(old,new)
    (BASE/'README.md').write_text(text,encoding='utf-8')
    (ROOT/'docs/sky130_phase_results.md').write_text(text.replace('(results/','(../spice/sky130/results/'),encoding='utf-8')
    return {'nominal_hybrid_energy':physical,'pvt_table':pvttable,'energy_break_even':energy}


def main():
    s=json.loads((OUT/'sky130_phase_summary.json').read_text())
    failure_attribution(s)
    energy=plots(s);summary=report(s,energy);dump(OUT/'engineering_conclusion.json',summary)
    source=OUT/'source_snapshot';source.mkdir(exist_ok=True)
    files=('gamma_sky130.py','gamma_sky_batch.py','run_sky130_phase.py','report_sky130_phase.py','validate_sky130.py','audit_sky130.py','package_sky130_phase.py','test_gamma_sky130.py')
    for name in files:shutil.copyfile(ROOT/'python'/name,source/name)
    dump(OUT/'source_manifest.json',{name:hashlib.sha256((ROOT/'python'/name).read_bytes()).hexdigest() for name in files})
    print(json.dumps(energy,indent=2),flush=True)


if __name__=='__main__':main()
