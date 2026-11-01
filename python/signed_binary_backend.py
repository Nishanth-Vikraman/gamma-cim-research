"""Fast ngspice batching: source the UNMODIFIED PDK once, reset per case.

Every per-case circuit.cir remains a standalone literal, reproducible deck.
The normalized template replaces only instance/source numbers with parameters.
This avoids repeated filesystem library loading; it does not alter PDK models.
"""
from __future__ import annotations
import hashlib,json,os,re,subprocess,time,shutil
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from gamma_sky130 import EXE,short,dump,fast_raw,compress_trace
BASE=Path(__file__).resolve().parents[1]/'results/signed_binary_kcl_20261026/backend'


def monitored_run(argv,cwd,env,timeout,creationflags=0,**kwargs):
    import psutil
    logpath=Path(cwd)/'progress.log'
    with logpath.open('wb') as out:
        proc=subprocess.Popen(argv,cwd=cwd,env=env,stdout=out,stderr=subprocess.STDOUT,creationflags=creationflags)
        start=time.monotonic();last=start
        while proc.poll() is None:
            try:proc.wait(timeout=15)
            except subprocess.TimeoutExpired:pass
            now=time.monotonic()
            if now-last>=30:
                try:
                    pp=psutil.Process(proc.pid);stat=pp.cpu_times()
                    print('simulator progress',Path(cwd).name,'wall_s',round(now-start),'cpu_s',round(stat.user+stat.system),'rss_mb',round(pp.memory_info().rss/1e6),flush=True)
                except psutil.Error:pass
                last=now
            if now-start>timeout:
                proc.kill();proc.wait();raise subprocess.TimeoutExpired(argv,timeout)
    return subprocess.CompletedProcess(argv,proc.returncode,logpath.read_bytes(),b'')

NUMBER=re.compile(r'(?<![\w.])[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?(?:meg|[fpnumkgt])?(?![\w.])',re.I)
SCALE={'f':1e-15,'p':1e-12,'n':1e-9,'u':1e-6,'m':1e-3,'k':1e3,'meg':1e6,'g':1e9,'t':1e12}


def numeric(s):
    m=re.fullmatch(r'(.+?)(meg|[fpnumkgt])?',s,re.I)
    return float(m[1])*SCALE.get((m[2] or '').lower(),1)


def parameterize(deck):
    values={};out=[]
    def replace(match):
        key='GX'+str(len(values));values[key]=numeric(match[0]);return '{'+key+'}'
    for line in deck.splitlines():
        if line.lower().startswith('.param mc_mm_switch='):
            enable=re.search(r'mc_mm_switch=(\d)',line,re.I)[1]
            values['MC_ENABLE']=float(enable);out.append('.param mc_mm_switch={MC_ENABLE} mc_pr_switch=0');continue
        if line.lower().startswith('.temp '):out.append('.temp 25');continue
        if line.lower().startswith('.tran '):out.append('.tran '+NUMBER.sub(replace,line[6:]));continue
        if line.startswith(('*','.')) or not line.strip():out.append(line);continue
        parts=line.split(maxsplit=3)
        if line[0].upper() in ('R','C','V','B') and len(parts)==4:
            out.append(' '.join(parts[:3])+' '+NUMBER.sub(replace,parts[3]));continue
        if line[0].upper()=='X' and ' w=' in line:
            a,b=line.split(' w=',1);out.append(a+' w='+NUMBER.sub(replace,b));continue
        out.append(line)
    return '\n'.join(out)+'\n',values


def batch_execute(cases,workers=2,processor=None):
    """cases=[(deck,metadata,absolute_case_directory), ...]. Return in input order."""
    groups=defaultdict(list);results=[None]*len(cases)
    for i,(deck,meta,directory) in enumerate(cases):
        directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
        replay='SKY130 seeded standalone replay\n.control\nset noaskquit\nset num_threads=1\nset filetype=binary\nsetseed '+str(meta.get('seed',101))+'\nsource circuit.cir\nsetseed '+str(meta.get('seed',101))+'\nreset\nset temp='+str(meta.get('design',{}).get('temp_c',25))+'\nrun\nwrite replay.raw all\nquit\n.endc\n.end\n'
        (directory/'replay.cir').write_text(replay,encoding='ascii')
        fingerprint=hashlib.sha256((deck+json.dumps(meta,sort_keys=True)+'sky_batch_v1').encode()).hexdigest()
        cache=directory/'result.json'
        if cache.exists():
            try:result=json.loads(cache.read_text())
            except json.JSONDecodeError:result={}
            log_ok=False
            if result.get('fingerprint')==fingerprint:
                ownlog=directory/'simulator.log'
                own=ownlog.read_text(errors='replace') if ownlog.exists() else ''
                needle=short(directory)+'/output.raw'
                if needle in own and 'COMPLETED_CASE_' in own:log_ok=True
                else:
                    oldbatch=Path(result.get('batch_path','__absent__'))/'simulator.log'
                    if oldbatch.exists():
                        text=oldbatch.read_text(errors='replace')
                        if needle in text and 'COMPLETED_CASE_' in text:
                            ownlog.write_text(text,encoding='utf-8');log_ok=True
            if result.get('fingerprint')==fingerprint and log_ok:
                if result.get('observations') or (directory/'output.raw').exists() or (directory/'output.raw.gz').exists():
                    if processor and not result.get('observations'):
                        names,x=fast_raw(directory/'output.raw');result=processor(result,names,x,directory);dump(cache,result)
                    results[i]=result;continue
        # Recover already-completed writes from an interrupted batch; never
        # recover an incomplete transient or a different standalone netlist.
        raw=directory/'output.raw'
        trace_exists=raw.exists() or Path(str(raw)+'.gz').exists()
        try:old_meta=json.loads((directory/'golden.json').read_text())
        except (OSError,json.JSONDecodeError):old_meta=None
        if trace_exists and (directory/'circuit.cir').read_text()==deck and old_meta==meta:
            matching=[]
            needle=short(directory)+'/output.raw'
            for logfile in (BASE/'batches').glob('*/simulator.log'):
                log=logfile.read_text(errors='replace')
                if needle in log and f'COMPLETED_CASE_{i}' in log:matching.append(logfile)
            if matching:
                try:
                    names,x=fast_raw(raw)
                    complete=not meta.get('stop_ns') or x[-1,0]>=(meta['stop_ns']-.01)*1e-9
                    if len(x)>2 and (complete or processor):
                        result={'fingerprint':fingerprint,'returncode':0,'rows':len(x),'seed':meta.get('seed',101),
                            'meta':meta,'recovered_from_log':str(matching[-1]),'batch_path':str(matching[-1].parent),
                            'simulation_status':'complete' if complete else 'incomplete_transient',
                            'last_time_ns':float(x[-1,0]*1e9)}
                        if processor:result=processor(result,names,x,directory)
                        (directory/'simulator.log').write_text(matching[-1].read_text(errors='replace'),encoding='utf-8')
                        dump(cache,result);results[i]=result;continue
                except (ValueError,IndexError):pass
        template,values=parameterize(deck)
        key=hashlib.sha256(template.encode()).hexdigest()
        (directory/'circuit.cir').write_text(deck,encoding='ascii');dump(directory/'golden.json',meta)
        groups[key].append((i,template,values,meta,directory,fingerprint))
    # Bound transient disk usage and circuit rebuilds; exit each process often.
    free=shutil.disk_usage(BASE).free
    chunk=max(1,min(2,int(max(0,free-100_000_000)/(max(1,workers)*35_000_000))))
    if free<100_000_000:raise RuntimeError('Less than100 MB free; stopping before creating transient outputs')
    groups={f'{key}_{offset//chunk}':group[offset:offset+chunk] for key,group in groups.items() for offset in range(0,len(group),chunk)}
    def run_group(key,group):
        identity=hashlib.sha256(''.join(item[5] for item in group).encode()).hexdigest()[:16]
        directory=BASE/'batches'/identity;directory.mkdir(parents=True,exist_ok=True)
        # Inline constants shared by a batch. Avoid thousands of needless
        # per-MOS parameters without modifying any PDK model contents.
        constants={k:v for k,v in group[0][2].items() if all(item[2][k]==v for item in group)}
        template=group[0][1]
        for k,v in constants.items():template=template.replace('{'+k+'}',f'{v:.16g}')
        group=[(i,t,{k:v for k,v in values.items() if k not in constants},m,target,fp) for i,t,values,m,target,fp in group]
        lines=template.splitlines()
        param='.param '+' '.join(f'{key}={value:.16g}' for key,value in group[0][2].items())
        if param=='.param ':param='* no varying instance parameters'
        # Actual first-case defaults avoid invalid time ordering during source.
        lines.insert(2,param)
        (directory/'template.cir').write_text('\n'.join(lines)+'\n',encoding='ascii')
        controls=['Sky130 batch; per-case standalone decks retained','.control','set noaskquit','set num_threads=2',
            'set filetype=binary','source template.cir']
        for index,_,values,meta,target,fp in group:
            controls+=['setseed '+str(meta.get('seed',101))]
            controls += [f'alterparam {k}={v:.16g}' for k,v in values.items()]
            temp=meta.get('design',{}).get('temp_c',25)
            controls+=['reset',f'set temp={temp:.12g}',*(['optran 0 0 0 10p 300p 0'] if meta.get('kind')=='signed_binary_kcl' else []),'run',*(['setplot noise1'] if meta.get('kind')=='noise' else []),f'write {short(target)}/output.raw all',f'echo COMPLETED_CASE_{index}','destroy all']
        controls+=['quit','.endc','.end']
        (directory/'run.cir').write_text('\n'.join(controls)+'\n',encoding='ascii')
        init=BASE/'init';init.mkdir(parents=True,exist_ok=True)
        (init/'spinit').write_text('set num_threads=2\nset filetype=ascii\n',encoding='ascii')
        env=dict(os.environ,OMP_NUM_THREADS='2',SPICE_SCRIPTS=short(init))
        start=time.monotonic()
        r=monitored_run([short(EXE),'-b','-o','simulator.log','run.cir'],cwd=short(directory),capture_output=True,
            env=env,timeout=max(600,90+len(group)*100),creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        (directory/'process.log').write_bytes(r.stdout+r.stderr)
        log=(directory/'simulator.log').read_text(errors='replace')
        done=[]
        for index,_,values,meta,target,fp in group:
            raw=target/'output.raw'
            if not raw.exists() or f'COMPLETED_CASE_{index}' not in log:
                result={'fingerprint':fp,'meta':meta,'simulation_status':'failed','failure':'Missing complete transient/waveform','batch_path':str(directory),'observations':True}
                result['metrics']={'n':meta.get('workload',{}).get('n',meta.get('n')),'pattern':meta.get('cfg',{}).get('pattern'),'passed':False,'simulation_complete':False,'run':str(target)}
                dump(target/'result.json',result);(target/'simulator.log').write_text(log,encoding='utf-8');done.append((index,result));continue
            names,x=fast_raw(raw)
            if x.shape[0]<2:raise RuntimeError('No transient samples: '+str(target))
            complete=not meta.get('stop_ns') or x[-1,0]>=(meta['stop_ns']-.01)*1e-9
            if not complete and not processor:
                raise RuntimeError('Incomplete transient: '+str(target))
            # Shared original simulator log is linked by path, not copied N times.
            (target/'simulator.log').write_text(log,encoding='utf-8')
            dump(target/'batch_parameters.json',{'template':str(directory/'template.cir'),'values':values,
                'seed':meta.get('seed',101),'temperature_c':meta.get('design',{}).get('temp_c',25)})
            result={'fingerprint':fp,'returncode':r.returncode,'rows':len(x),'seed':meta.get('seed',101),
                'meta':meta,'batch_path':str(directory),'batch_seconds':time.monotonic()-start,
                'simulation_status':'complete' if complete else 'incomplete_transient',
                'last_time_ns':float(x[-1,0]*1e9)}
            if processor:result=processor(result,names,x,target)
            dump(target/'result.json',result);done.append((index,result))
        print('batch complete',key[:8],len(group),'cases',round(time.monotonic()-start,1),'s',flush=True)
        return done
    print('ngspice new batches',len(groups),'new cases',sum(map(len,groups.values())),'cached',sum(r is not None for r in results),flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(run_group,key,group) for key,group in groups.items()]
        for future in as_completed(futures):
            for index,result in future.result():results[index]=result
    return results
