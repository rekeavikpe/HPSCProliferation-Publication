"""Targeted stop and checked migration for the no-automatic-reduction correction."""
import argparse,datetime,hashlib,json,os,pathlib,shutil,signal,subprocess,tempfile,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
OLD={'scripts/run_october_cells.py':'241c24e26c0b5c87b472ebc6b6a559109b39086941f05b6cd8cdd93ea112fc54',
     'scripts/report_october20.py':'ef3d5f5876da4c9a793cecfe23049ef4693301ebea98f664b269a639e86d3f91'}
def processes():
    rows={}
    for line in subprocess.check_output(['ps','-eo','pid=,ppid=,stat=,args='],text=True).splitlines():
        p=line.strip().split(None,3)
        if len(p)==4: rows[int(p[0])]=(int(p[1]),p[2],p[3])
    return rows

def targets(out):
    rows=processes(); hits=set()
    for pid,(_,state,cmd) in rows.items():
        if pid==os.getpid() or state.startswith('Z'): continue
        # Launcher with default output, or explicit requested output. Also catch
        # orphaned workers whose arguments contain this suite output directory.
        words=cmd.split()
        launcher=next((i for i,w in enumerate(words) if w.endswith('/run_october20_revision.sh')),None)
        if launcher is not None:
            arg=words[launcher+1] if len(words)>launcher+1 else 'results/david_comments_20261020'
            if (ROOT/arg).resolve()==out: hits.add(pid)
        if any((ROOT/w).resolve()==out or str((ROOT/w).resolve()).startswith(str(out)+os.sep)
               for w in words[1:] if not w.startswith('-')):
            if any(x in cmd for x in ['julia','run_october_cells.py','tee ']): hits.add(pid)
    changed=True
    while changed:
        added={pid for pid,(parent,_,_) in rows.items() if parent in hits}-hits
        changed=bool(added);hits.update(added)
    return hits

def send(pid,sig):
    try: os.kill(pid,sig)
    except ProcessLookupError: pass

def stop(out):
    frozen=set()
    # Freeze launchers and workers so no next stage starts while stopping.
    for _ in range(10):
        new=targets(out)-frozen
        if not new: break
        for pid in new: send(pid,signal.SIGSTOP)
        frozen.update(new)
    if not frozen:
        print('No running October 20 suite found for',out);return
    for pid in frozen: send(pid,signal.SIGTERM)
    for pid in frozen: send(pid,signal.SIGCONT)
    for _ in range(30):
        live=processes()
        remaining={p for p in frozen if p in live and not live[p][1].startswith('Z')}
        if not remaining:
            print('Stopped suite. Completed checkpoints remain available.');return
        time.sleep(1)
    raise SystemExit('Processes still exiting. Do not migrate yet. Retry stop. PIDs: '+str(sorted(remaining)))

def migrate(out):
    if targets(out): raise SystemExit('Suite is running. Run the stop command first.')
    target=out/'suite_manifest.json';old=json.loads(target.read_text())
    sources=old['sources'];env=os.environ.copy()
    for key,value in old['settings'].items():
        if value is None: env.pop(key,None)
        else: env[key]=value
    # Recover the actual source directories from the locked manifest.
    def unique(matches,label):
        if len(matches)!=1: raise SystemExit('Cannot identify '+label+' uniquely from manifest')
        return matches[0]
    fit=unique([str(pathlib.Path(p).parent) for p in sources if p.endswith('/selection_optimizer_bic.csv') and not p.endswith('/separate_ablations/selection_optimizer_bic.csv')],'reference fit')
    split=unique([str(pathlib.Path(p).parent.parent) for p in sources if p.endswith('/model1/summary.csv')],'split fit')
    p75=[str(pathlib.Path(p).parent.parent) for p in sources if p.endswith('/lag_p75/summary.csv')]
    if p75: env['HPSC_P75_SOURCE']=unique(p75,'P75 source')
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(['python3','scripts/october20_manifest.py',fit,split,tmp],cwd=ROOT,env=env,check=True)
        current=json.loads((pathlib.Path(tmp)/'suite_manifest.json').read_text())
    if old==current:
        print('Manifest already current. Resume normally.');return
    if old['settings']!=current['settings']: raise SystemExit('Settings changed. Migration refused.')
    changed={p for p in sources.keys()|current['sources'].keys() if sources.get(p)!=current['sources'].get(p)}
    allowed={str(ROOT/p):h for p,h in OLD.items()}
    helper=str(pathlib.Path(__file__).resolve())
    for p in changed:
        if p==helper and p not in sources: continue
        if p not in allowed or sources.get(p)!=allowed[p] or p not in current['sources']:
            raise SystemExit('Unrelated source/input change. Migration refused: '+p)
    archive=out.parent/(out.name+'_before_no_drop_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    archive.mkdir();shutil.copy2(target,archive/'suite_manifest.json')
    # Archive outside the suite so recursive report/collection cannot mix rounds.
    items=[out/p for p in ['cells.DONE','figures.DONE','report.DONE','COMPUTATIONS_COMPLETE','index.html','review_status.json','figures','cells/status.json','cells/final_rounds.csv','cells/COMPUTATIONS_COMPLETE']]
    items += [p for p in (out/'cells').glob('*/round*') if p.name!='round1']
    for p in items:
        if p.exists():
            dest=archive/p.relative_to(out);dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(p),dest)
    record={'changed_sources':sorted(changed),'archive':str(archive),'preserved':'all upstream stages and cell round1 checkpoints','policy':'no automatic parameter removal'}
    (archive/'migration.json').write_text(json.dumps(record,indent=2))
    temporary=out/'suite_manifest.json.new';temporary.write_text(json.dumps(current,indent=2));temporary.replace(target)
    print('Migration complete. Preserved completed fits. Archive:',archive)
    print('Resume scripts/run_october20_revision.sh with the SAME original budgets and output path.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['stop','migrate']);parser.add_argument('--out',default='results/david_comments_20261020');args=parser.parse_args()
    os.chdir(ROOT);out=pathlib.Path(args.out).resolve()
    (stop if args.action=='stop' else migrate)(out)
