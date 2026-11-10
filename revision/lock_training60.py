"""Reject mixed-source resumes while preserving all completed output."""
import hashlib, json, os, sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
split, exposure, out=map(lambda x:Path(x).resolve(),sys.argv[1:])
required=[split/'calibration.toml',split/'split.csv',split/'protocol.toml',
          split/'calibration_observations.csv',split/'selection_observations.csv',
          exposure/'COMPLETE',exposure/'ranking.csv',exposure/'split.csv']
required += [split/f'model{i}'/f for i in range(7) for f in ('selected_parameters.csv','summary.csv')]
required += list(exposure.glob('*/selected_parameters.csv'))
if len(list(exposure.glob('*/selected_parameters.csv')))!=11:
    raise SystemExit('Expected all eleven saved 60% exposure fits')
for folder in ('src','revision','profiles','population','ablations','scripts','figures'):
    required += [p for p in (root/folder).rglob('*') if p.is_file() and p.suffix in ('.jl','.py','.sh','.toml')]
required += list((root/'data/processed').glob('*.csv'))
manifest={'split':str(split),'exposure':str(exposure),
          'hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(required))},
          'budgets':{k:os.environ.get(k,v) for k,v in {
              'HPSC_ABLATION_STARTS':'100','HPSC_ABLATION_EVALUATIONS':'20000',
              'HPSC_ABLATION_ITERATIONS':'3000','HPSC_REVIEW_STARTS':'100',
              'HPSC_REVIEW_EVALUATIONS':'20000','HPSC_REVIEW_ITERATIONS':'3000',
              'HPSC_PROFILE_RESOLUTION':'100'}.items()}}
out.mkdir(parents=True,exist_ok=True)
path=out/'workflow_manifest.json'
if path.exists() and json.loads(path.read_text())!=manifest:
    raise SystemExit('Inputs, code or budgets changed. Preserve this output and choose a new output folder.')
path.write_text(json.dumps(manifest,indent=2)+'\n')
print('Verified workflow sources and resume settings')
