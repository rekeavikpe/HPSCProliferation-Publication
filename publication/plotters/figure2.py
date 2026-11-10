from pathlib import Path
import csv,shutil,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
CACHE=ROOT/'publication/cache'
OUT=ROOT/'publication/generated/figures'
OUT.mkdir(parents=True,exist_ok=True)
W=ROOT;O=OUT;F=OUT;M=F;S=F;Q=F
for p in [M,S,Q]:p.mkdir(parents=True,exist_ok=True)
R=CACHE/'training60';D=CACHE/'october20';B=ROOT
def rows(p):return list(csv.DictReader(p.open()))
def style(fs=14):plt.rcParams.update({'font.family':'Arial','font.size':fs,'axes.labelsize':fs,'axes.titlesize':fs,'xtick.labelsize':fs-1,'ytick.labelsize':fs-1,'legend.fontsize':fs-1,'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold','pdf.fonttype':42})
def save(f,folder,name):
 for ext in ['png','pdf']:f.savefig(folder/(name+'.'+ext),dpi=350,bbox_inches='tight')
 plt.close(f)
def copy(src,dest):
 for ext in ['png','pdf']:shutil.copy2(W/(src+'.'+ext),dest.with_suffix('.'+ext))
style(17)
a=rows(D/'reuse/candidate_selection_optimizer_bic.csv');f,ax=plt.subplots(figsize=(12.5,6));models=[0,1,3,4,5,6,7];colors=['#777777','#dd3434','#ed9b00','#c8a517','#90248e','#17851c','#99504b'];names=['No hydro.\nterms','Michaelis–\nMenten','Exponential','Growth-only\nshear','Hill\n(n = 2)','Growth-coupled\nEDR','Dual\ninteraction']
for i,(m,col) in enumerate(zip(models,colors)):
 vals=[float(v['bic']) for v in a if int(v['model'])==m];assert len(vals)==100
 bx=ax.boxplot([vals],positions=[i],widths=.5,patch_artist=True,showfliers=False,medianprops={'color':'black','lw':1.5});bx['boxes'][0].set(facecolor=col,alpha=.3);ax.scatter(i+np.random.default_rng(m).uniform(-.12,.12,len(vals)),vals,s=18,color=col,alpha=.55)
ax.set(ylabel='BIC',xticks=range(7),xticklabels=[f'Model {m}\n{n}' for m,n in zip(models,names)]);ax.grid(axis='y',alpha=.15);ax.margins(y=.08);f.tight_layout();save(f,M,'Figure2_selection_BIC')
