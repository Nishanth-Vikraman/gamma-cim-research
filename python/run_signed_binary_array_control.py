"""Paired physical array-only cost, for incremental full-chain accounting."""
from signed_binary_kcl import run_configs,NS
rs=run_configs([{'n':n,'pattern':'aligned','array_only':True} for n in NS],'array_only',2)
print('Array-only physical controls complete',len(rs),flush=True)
