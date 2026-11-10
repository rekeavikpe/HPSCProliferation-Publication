import hashlib,json,os,pathlib,sys
source,out=map(pathlib.Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=True)
files=[]
for folder in ['src','revision','profiles','population','ablations','data/processed']:
    files.extend(p for p in pathlib.Path(folder).rglob('*') if p.is_file() and p.suffix in ['.jl','.csv','.toml'])
files.extend(pathlib.Path(p) for p in ['scripts/run_cell_reductions.py','scripts/run_october20_followup.sh','scripts/lock_october20_followup.py','figures/render_day1_saved.jl','figures/render_cell_reductions.jl','figures/ManuscriptStyle.jl'])
for root in [source,pathlib.Path('results/training60_final/reference_fit'),pathlib.Path('results/split6040_refined_hybrid_lac100_kg100')]:
    files.extend(p for p in root.rglob('*') if p.is_file() and p.name in ['selected_parameters.csv','reported_parameters.csv','audit.csv','calibration.toml','split.csv','initial_states.csv'])
record={'hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(files))},'settings':{n:os.environ.get(n,'') for n in ['HPSC_REVIEW_STARTS','HPSC_REVIEW_EVALUATIONS','HPSC_REVIEW_ITERATIONS','HPSC_PROFILE_RESOLUTION','HPSC_PROFILE_MAXITERS']}}
p=out/'manifest.json'
if p.exists() and json.loads(p.read_text())!=record:raise SystemExit('Inputs or settings changed. Keep existing outputs and use a new output directory.')
p.write_text(json.dumps(record,indent=2));print('Verified follow-up inputs and budgets')
