"""Independent literal-deck replay of the selected loaded transistor path."""
import json,numpy as np
from gamma_sky130 import *

def main():
    source=BASE/'pvt/penergy004'; target=BASE/'loaded_mac/standalone_audit'
    meta=json.loads((source/'golden.json').read_text()); deck=(source/'circuit.cir').read_text()
    expected=json.loads((source/'metrics.json').read_text())
    r,names,x=execute(deck,meta,target,timeout=150)
    assert x[-1,0]>=(meta['stop_ns']-.01)*1e-9
    actual=evaluate_boundary(r,names,x);dump(target/'metrics.json',actual)
    errors={}
    for key in ('stored_v','driven_v','mac_v','measured_known_path_fj'):
        a=np.array([o[key] for o in actual['observations']]); b=np.array([o[key] for o in expected['observations']])
        errors[key]=float(max(abs(a-b)))
    classes_match=[o['mac_class'] for o in actual['observations']]==[o['mac_class'] for o in expected['observations']]
    passed=classes_match and max(errors[k] for k in ('stored_v','driven_v','mac_v'))<1e-6 and errors['measured_known_path_fj']<.1
    out={'passed':passed,'literal_deck_vs_batch_max_absolute_differences':errors,'gamma_classes_identical':classes_match,
         'scope':'Independent source of own literal circuit, no numeric parameterization, same frozen PDK.'}
    dump(BASE/'results/standalone_path_audit.json',out);compress_trace(target/'output.raw')
    print(out,flush=True);assert passed

if __name__=='__main__':main()
