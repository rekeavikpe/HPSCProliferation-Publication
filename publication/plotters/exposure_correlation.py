from pathlib import Path
import csv,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
CACHE=ROOT/'publication/cache'
OUT=ROOT/'publication/generated/figures'
OUT.mkdir(parents=True,exist_ok=True)
S=OUT;a=list(csv.DictReader((CACHE/'october20/reuse/exposure_distributions.csv').open()));conditions=[(100,30),(100,40),(100,60),(500,30)]
def val(s,v,r,field):return float(next(x[field]for x in a if x['scenario']==s and int(x['volume_ml'])==v and int(x['rpm'])==r and(s.startswith('eul')or int(float(x['size_um']))==200)))
plt.rcParams.update({'font.family':'Arial','font.size':17,'axes.labelsize':17,'xtick.labelsize':15,'ytick.labelsize':15,'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','pdf.fonttype':42})
f,axs=plt.subplots(1,2,figsize=(12.8,5.5));f.subplots_adjust(left=.09,right=.985,top=.91,bottom=.25,wspace=.37)
for j,(ax,field,title,unit)in enumerate(zip(axs,['shear_pa','edr_m2_s3'],['Shear stress','EDR'],['Pa',r'm$^2$ s$^{-3}$'])):
 x=[val('eul_va',v,r,field)for v,r in conditions];y=[val('lag_p50',v,r,field)for v,r in conditions];high=max(x+y)*1.1
 ax.plot([0,high],[0,high],'--',color='#777777',lw=1.6)
 offsets=([(-36,10),(12,2),(-50,12),(-65,-24)]if j==0 else[(-58,-25),(14,-12),(-52,12),(-62,12)])
 for i,(v,r)in enumerate(conditions):
  ax.scatter(x[i],y[i],color='#0072B2'if v==100 else'#D55E00',marker='o'if v==100 else's',s=65)
  ax.annotate(f'{r} RPM',(x[i],y[i]),xytext=offsets[i],textcoords='offset points',fontsize=13)
 ax.set(xlabel=f'Eulerian volume average ({unit})',ylabel=f'Lagrangian P50 ({unit})',title=f'({chr(97+j)}) {title}',xlim=(0,high),ylim=(0,high));ax.ticklabel_format(axis='both',style='sci',scilimits=(-3,3),useMathText=True);ax.grid(alpha=.15)
f.legend([Line2D([],[],marker='o',ls='none',color='#0072B2',ms=8),Line2D([],[],marker='s',ls='none',color='#D55E00',ms=8),Line2D([],[],ls='--',color='#777777')],['100 mL','500 mL','1:1'],loc='lower center',ncol=3,bbox_to_anchor=(.5,.015),fontsize=15)
for ext in ['png','pdf']:f.savefig(S/('FigureS1_exposure_correlation.'+ext),dpi=350,bbox_inches='tight')
