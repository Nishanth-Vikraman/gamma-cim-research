"""Frozen characterization queue. All primary cases use real PDK devices.
Run after inspecting primary smoke. Does not rerun any historical q4 circuit.
"""
import itertools,sys,subprocess
from pathlib import Path
from signed_binary_kcl import OUT,NS,run_configs
OUT.mkdir(parents=True,exist_ok=True)
configs=[{'n':n,'pattern':pat} for n,pat in itertools.product(NS,('aligned','aligned_negative','dense','sparse','cancellation'))]
rs=run_configs(configs,'primary',6)
print('PRIMARY FINISHED',sum(r['metrics'].get('passed',False) for r in rs),len(rs),flush=True)
for phase in ('sizing','pvt','mc','pareto'):
    print('STARTING',phase,flush=True)
    subprocess.run([sys.executable,'-u','-B',str(Path(__file__).with_name('signed_binary_kcl.py')),'--phase',phase,'--workers','6'],check=True)
configs=[{'n':n,'pattern':'aligned','res_model':'high','unit_r':r} for n,r in itertools.product((64,2048),(1e4,1e5))]
run_configs(configs,'high_resistor',4)
print('ALL MAIN CHARACTERIZATION FINISHED',flush=True)
