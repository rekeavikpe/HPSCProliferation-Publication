"""Render the final paper's plots from saved numerical results; no fitting."""
from pathlib import Path
import argparse,subprocess,sys,os,shutil
ROOT=Path(__file__).resolve().parent.parent
p=argparse.ArgumentParser();p.add_argument('--profiles',action='store_true',help='Also render native Julia profile plots S5/S6');p.add_argument('--ipsc-python',default=sys.executable,help='Optional Python with Matplotlib 3.10.8 for original S7/S8 typography');args=p.parse_args()
for name in ['figure2','figure3','figure4','exposure_correlation','exposure_sweeps','ipsc_fits','ipsc_profiles']:
 print('Rendering',name,flush=True)
 subprocess.run([args.ipsc_python if name.startswith('ipsc_') else sys.executable,str(ROOT/'publication/plotters'/f'{name}.py')],check=True,cwd=ROOT)
if args.profiles:
 julia=os.environ.get('HPSC_JULIA','julia')
 for name,cache in [('nine_profiles','nine_profiles'),('two_diagnostics','two_diagnostics')]:
  subprocess.run([julia,'--project='+str(ROOT/'profile_figures'),str(ROOT/'publication/plotters'/f'{name}.jl'),str(ROOT/'publication/cache'/cache),str(ROOT/'publication/generated/figures')],check=True,cwd=ROOT)
print('Saved publication/generated/figures. Figure 1 and S4 are supplied coauthor artwork in paper_figures; their original drawing/CFD projects are not included.')
