"""Exhaustive reductions of currently unbounded parameters, with no dropping order."""
import csv,itertools,json,os,pathlib,shlex,subprocess,sys
PARAMETERS=['mu_max','K_d','K_oG','Y_GLac','alpha']
def candidates(audit):
    if any(r['status']=='numerical_review' for r in audit):
        raise RuntimeError('Resolve numerical profile flags before parameter reduction')
    if {r['parameter'] for r in audit}!=set(PARAMETERS):
        raise RuntimeError('Expected the five requested original profiles')
    bad=[n for n in PARAMETERS if next(r['status'] for r in audit if r['parameter']==n)=='not_bounded']
    return [(list(fixed),[n for n in PARAMETERS if n not in fixed])
            for count in range(1,len(bad)+1) for fixed in itertools.combinations(bad,count)]
def read(path):
    with open(path) as f:return list(csv.DictReader(f))
def main():
    fit,source,out=map(pathlib.Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=True)
    j=shlex.split(os.environ.get('HPSC_JULIA_COMMAND','julia +1.12.4'))
    def run(project,script,*args):subprocess.run(j+['--project='+project,script,*map(str,args)],check=True)
    results=[]
    for line in ['585A1','SCTi003-A']:
        original=source/'cells'/line/'round1'
        original_audit=read(original/'profiles/audit.csv')
        designs=[([],PARAMETERS,original)]+[(fixed,free,out/line/('fixed_'+'_'.join(fixed))) for fixed,free in candidates(original_audit)]
        for fixed,free,dest in designs:
            if fixed:
                if not (dest/'REDUCTION_COMPLETE').exists():
                    print(f'START {line}: fix {fixed}, fit {free}',flush=True)
                    run('population','revision/fit_october_cells.jl',fit,line,dest,','.join(free))
                    run('profiles','profiles/run_october_cells.jl',dest,dest/'profiles')
                    run('population','revision/export_october_cells.jl',fit,dest)
                    (dest/'REDUCTION_COMPLETE').write_text('Fit and profiles computed, inspect comparison.json\n')
            scorepath=out/(line+'_'+('fixed_'+'_'.join(fixed) if fixed else 'original5')+'_score.csv')
            run('population','revision/score_cell_reduction.jl',dest,scorepath)
            audit=read(dest/'profiles/audit.csv');score=read(scorepath)[0]
            if {a['parameter'] for a in audit}!=set(free):raise RuntimeError('Profile set mismatch')
            errors={a['state']:float(a['mrpe_percent']) for a in read(dest/'errors.csv') if a['tag']=='Recalibrated' and a['condition']=='0'}
            results.append(dict(line=line,fixed=fixed,free=free,fit=str(dest),all_profiles_bounded=all(a['status']=='bounded' for a in audit),audit=audit,nll=float(score['nll']),bic=float(score['bic']),n=int(score['n']),k=int(score['k']),errors=errors))
            (out/'comparison.json').write_text(json.dumps(results,indent=2))
    (out/'COMPUTATIONS_COMPLETE').write_text('Inspect bounded profiles, likelihood changes and fit errors. No model is automatically selected. Fit BIC uses the existing cell-line weighting, not independent selection data.\n')
if __name__=='__main__':main()
