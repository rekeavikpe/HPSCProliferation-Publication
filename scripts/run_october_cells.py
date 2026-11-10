"""Fit and profile David's five parameters. Never choose a reduction automatically."""
import csv,json,os,pathlib,shlex,subprocess,sys
fit,out=map(pathlib.Path,sys.argv[1:3]);out.mkdir(parents=True,exist_ok=True)
julia=shlex.split(os.environ.get('HPSC_JULIA_COMMAND','julia +1.12.4'))
initial=['mu_max','K_d','K_oG','Y_GLac','alpha']
summary=[]
for line in ['585A1','SCTi003-A']:
    dest=out/line/'round1'
    subprocess.run(julia+['--project=population','revision/fit_october_cells.jl',str(fit),line,str(dest),','.join(initial)],check=True)
    profiles=dest/'profiles'
    subprocess.run(julia+['--project=profiles','profiles/run_october_cells.jl',str(dest),str(profiles)],check=True)
    subprocess.run(julia+['--project=population','revision/export_october_cells.jl',str(fit),str(dest)],check=True)
    with (profiles/'audit.csv').open() as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=5 or {r['parameter'] for r in rows}!=set(initial):
        raise RuntimeError(f'{line}: expected exactly the five requested profile audits')
    if any(r['status'] not in ['bounded','not_bounded','numerical_review'] for r in rows):
        raise RuntimeError(f'{line}: unexpected profile audit status')
    unbounded=[r['parameter'] for r in rows if r['status']=='not_bounded']
    numerical=[r['parameter'] for r in rows if r['status']=='numerical_review']
    status=('profile_numerical_review_required' if numerical else
            'parameter_reduction_review_required' if unbounded else 'five_parameter_set_bounded')
    summary.append(dict(line=line,status=status,final_fit=str(dest),active_parameters=initial,
                        unbounded_parameters=unbounded,numerical_review_parameters=numerical,
                        history=[],automatic_reduction=False))
    (out/'status.json').write_text(json.dumps(summary,indent=2))
with (out/'final_rounds.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['line','status','final_fit']);w.writeheader()
    w.writerows({k:r[k] for k in ['line','status','final_fit']} for r in summary)
(out/'COMPUTATIONS_COMPLETE').write_text('Five-parameter computations finished. Inspect status.json. No parameters were automatically removed.\n')
