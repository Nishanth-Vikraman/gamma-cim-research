"""Independent literal decks versus parameterized batch evidence, new phase only."""
import json,hashlib
from concurrent.futures import ThreadPoolExecutor
from gamma_sky130 import execute,dump,compress_trace
from gamma_sar_exact import OUT
from gamma_sar_circuit import analyze,comparator_analyze

def one(label):
    src=OUT/'spice'/label;meta=json.loads((src/'golden.json').read_text())
    target=OUT/'audit'/label
    r,n,x=execute((src/'circuit.cir').read_text(),meta,target,timeout=300)
    f=comparator_analyze if meta['kind']=='comparator' else analyze
    obs=f(r,n,x);ref=json.loads((src/'metrics.json').read_text())
    field='mean_comparison_fj' if meta['kind']=='comparator' else 'mean_event_fj'
    error=abs(obs[field]/ref[field]-1)
    if meta['kind']=='comparator':same=obs['correct_nonzero']==ref['correct_nonzero']
    else:same=all(obs[k]==ref[k] for k in ('quotient_correct','transfer_correct','class_correct','sleep_correct'))
    summary={'label':label,'literal_deck_sha256':hashlib.sha256((src/'circuit.cir').read_bytes()).hexdigest(),
      'batch_fj':ref[field],'literal_fj':obs[field],'relative_energy_difference':error,'correctness_counts_agree':same,
      'passed':error<.005 and same,'seed':meta['seed'],'run':str(target.relative_to(OUT))}
    dump(target/'metrics.json',obs);compress_trace(target/'output.raw')
    print(summary,flush=True);return summary

def main():
    with ThreadPoolExecutor(2) as p:rows=list(p.map(one,['pvt04','cmp_balanced_w0.42_l0.15']))
    dump(OUT/'literal_audit.json',rows)
    if not all(v['passed'] for v in rows):raise RuntimeError('Literal replay discrepancy: inspect retained logs and metrics')

if __name__=='__main__':main()
