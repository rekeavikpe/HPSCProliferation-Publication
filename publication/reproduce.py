"""Verify the publication release and reproduce numerical summaries without refitting."""
from pathlib import Path
import argparse,subprocess,os,sys,hashlib,csv,math,tempfile,shutil,json
P=Path(__file__).resolve().parent;ROOT=P.parent
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--integrity-only',action='store_true');args=parser.parse_args()
manifest=ROOT/'SHA256SUMS.txt'
if manifest.exists():
 count=0
 for line in manifest.read_text().splitlines():
  digest,name=line.split('  ',1);path=ROOT/name
  if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:raise SystemExit('Integrity failure: '+name)
  count+=1
 print('PASS: verified',count,'release files',flush=True)
else:raise SystemExit('Missing SHA256SUMS.txt; use the complete release archive.')
if args.integrity_only:sys.exit(0)
julia=os.environ.get('HPSC_JULIA','julia');env=os.environ.copy();env.setdefault('OPENBLAS_NUM_THREADS','1')
(P/'verification').mkdir(exist_ok=True);(P/'generated').mkdir(exist_ok=True)
def run(project,script,*rest):
 subprocess.run([julia,'--project='+str(ROOT/project),str(ROOT/script),*map(str,rest)],check=True,cwd=ROOT,env=env)
for project,script in [('population','verify_primary.jl'),('population','verify_additional.jl'),('population','verify_cross_validation.jl'),('profiles','verify_nine_profiles.jl'),('profiles','verify_two_diagnostics.jl')]:
 print('Checking',script,flush=True);run(project,'publication/'+script)
with tempfile.TemporaryDirectory(prefix='hpsc_publication_')as td:
 for name,script in [('media','predict_media_exchange.jl'),('cell_transfer','predict_cell_lines.jl')]:
  out=Path(td)/name;run('.','scripts/'+script,P/'cache/training60/reference_fit',out)
  # Only generated output is replaced; original caches are never modified.
  dest=P/'generated'/name
  if dest.exists():shutil.rmtree(dest)
  shutil.copytree(out,dest)
run('revision/statistics','revision/october_statistics.jl',P/'generated/statistics')
subprocess.run([sys.executable,str(P/'fetch_reference.py')],cwd=ROOT,check=True)
run('population','scripts/oxygen_sensitivity.jl',P/'generated/oxygen_reference')
subprocess.run([sys.executable,str(ROOT/'scripts/october_oxygen.py'),str(P/'generated/oxygen_reference'),str(P/'generated/oxygen')],cwd=ROOT,check=True)
subprocess.run([sys.executable,str(P/'report_tables.py')],cwd=ROOT,check=True)
print('All publication numerical checks passed. Results: publication/verification and publication/generated/tables.')
