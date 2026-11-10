from pathlib import Path
import csv,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pypdf import PdfReader
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
CACHE=ROOT/'publication/cache'
OUT=ROOT/'publication/generated/figures'
OUT.mkdir(parents=True,exist_ok=True)
W=ROOT;O=OUT;R=CACHE/'ipsc'
plt.rcParams.update({'font.family':'Arial','font.size':16,'axes.labelsize':18,'xtick.labelsize':14,'ytick.labelsize':14,'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold','mathtext.sf':'Arial','mathtext.fallback':None,'pdf.fonttype':42})
def rows(p):return list(csv.DictReader(p.open()))
def fmt(v):
 if v<.001:return f'{v:.2e}'
 if v<1:return f'{v:.4f}'.rstrip('0')
 return f'{v:.1f}'
f,axs=plt.subplots(2,3,figsize=(15,8.4));f.subplots_adjust(left=.085,right=.985,top=.88,bottom=.09,hspace=.65,wspace=.36)
params=['mu_max','K_oG','alpha'];labels=[r'$\mu_{max}$ (h$^{-1}$)',r'$K_{oG}$ (mL/cell)',r'$\alpha$ (µm)'];threshold=-1.920729410347062
for i,line in enumerate(['585A1','SCTi003-A']):
 p=R/line/'fixed_K_d_Y_GLac'/'profiles';audit={r['parameter']:r for r in rows(p/'audit.csv')}
 for j,(param,label) in enumerate(zip(params,labels)):
  a=axs[i,j];rr=rows(p/param/'profile.csv');x=np.array([float(r['value']) for r in rr]);y=np.array([float(r['relative_loglik']) for r in rr]);info=audit[param]
  a.plot(x,y,color='#0064ba',lw=2.3);a.fill_between(x,threshold,y,where=y>=threshold,color='#7770ff',alpha=.6,interpolate=True)
  a.axhline(threshold,color='red',ls='--',lw=1.2);a.axvline(float(info['mle']),color='red',ls='--',lw=1.2)
  a.set_xscale('log');a.set_ylim(-2.9,.15);a.set_xlabel(label);a.set_ylabel([r'$\ell_p^*(\mu_{max})-\ell^*$',r'$\ell_p^*(K_{oG})-\ell^*$',r'$\ell_p^*(\alpha)-\ell^*$'][j]);a.grid(alpha=.2)
  a.set_title(f'({chr(97+i*3+j)}) 95% CI: ({fmt(float(info["lower"]))}, {fmt(float(info["upper"]))})',loc='left',fontsize=14,pad=9)
 pos=axs[i,0].get_position();f.text(.085,pos.y1+.065,line,fontsize=20,weight='bold')
for ext in ['png','pdf']:f.savefig(O/('FigureS8_iPSC_profiles.'+ext),dpi=350,bbox_inches='tight')
plt.close(f)
