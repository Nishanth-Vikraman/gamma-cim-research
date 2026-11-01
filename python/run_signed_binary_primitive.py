"""Physical loaded/unloaded downward-stage transfer grid; external inputs only."""
import itertools
from signed_binary_kcl import run_configs
configs=[]
for precision,cap,r,z in itertools.product((1,2),(20,100,500,2000),(-.2,0,.2),(0,.2)):
    configs.append({'n':64,'pattern':'zero_weights','precision':precision,'guards':0,
        'fixture_r':r,'fixture_z':z,'cap_ff':cap,'period_ns':400,'stage_start_ns':100,'stage_dt_ns':20})
rs=run_configs(configs,'primitive',2)
print('Primitive transfer cases',len(rs),flush=True)
