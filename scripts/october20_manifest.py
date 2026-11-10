import hashlib,json,os,pathlib,sys
fit,split,out=map(pathlib.Path,sys.argv[1:4]);root=pathlib.Path.cwd()
required=[fit/'calibration.toml',fit/'selection_optimizer_bic.csv',split/'split.csv']
required += [split/f'model{i}'/name for i in range(7) for name in ['selected_parameters.csv','summary.csv']]
required += [fit.parent/'separate_ablations'/n for n in ['manifest.toml','ranking.csv','selection_optimizer_bic.csv']]
required += [fit.parent/'separate_ablations'/s/'selected_parameters.csv' for s in ['baseline','full','shear_cell','shear_aggregate','shear_both','edr_cell','edr_aggregate','edr_both']]
required += [fit/f'model{i}/selected_parameters.csv' for i in [0,5]]
for f in required:
    if not f.is_file(): raise SystemExit(f'Missing required saved result: {f}')
p75=pathlib.Path(os.environ.get('HPSC_P75_SOURCE','results/david_revision_20261014/exposure_selection'))
if (p75/'lag_p75/selected_parameters.csv').exists():
    required += [p75/n for n in ['manifest.toml','split.csv','lag_p75/selected_parameters.csv','lag_p75/summary.csv']]
code=[]
for folder in ['src','revision','profiles','population','ablations','figures','scripts','config','data/processed','data/reference/neto2022']:
    code += [p for p in (root/folder).glob('*') if p.is_file() and p.suffix in ['.jl','.py','.sh','.toml','.csv','.txt','.m']]
manifest={'sources':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(required+code))},'settings':{k:os.environ.get(k) for k in ['HPSC_REVIEW_STARTS','HPSC_REVIEW_EVALUATIONS','HPSC_REVIEW_ITERATIONS','HPSC_PROFILE_RESOLUTION','HPSC_PROFILE_MAXITERS']}}
target=out/'suite_manifest.json'
if target.exists() and json.loads(target.read_text())!=manifest: raise SystemExit('Code, inputs or budgets changed. Use a new output folder to avoid mixing analyses.')
target.write_text(json.dumps(manifest,indent=2))
print('Verified required saved fits and locked inputs. Existing completed stages can be resumed.')
