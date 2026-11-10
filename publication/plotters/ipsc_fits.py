from pathlib import Path
import csv,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pypdf import PdfReader
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
CACHE=ROOT/'publication/cache'
OUT=ROOT/'publication/generated/figures'
OUT.mkdir(parents=True,exist_ok=True)
W=ROOT;O=OUT;R=CACHE/'ipsc'
O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':16,'axes.labelsize':17,'xtick.labelsize':14,'ytick.labelsize':14,'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold','mathtext.sf':'Arial','mathtext.cal':'Arial','mathtext.fallback':None,'pdf.fonttype':42})
def rows(p):return list(csv.DictReader(p.open()))
f,axs=plt.subplots(2,4,figsize=(17.5,10.6));f.subplots_adjust(left=.075,right=.99,bottom=.22,top=.87,hspace=.64,wspace=.4)
colors=['#2C4B80','#80392D'];states=['cells','glucose','lactate','aggregate'];labels=[r'Viable cells ($\times10^6$ cells/mL)','Glucose (mM)','Lactate (mM)','Aggregate size (µm)'];lims=[(0,1.6),(0,25),(0,30),(0,350)]
for i,(line,rpm,color) in enumerate(zip(['585A1','SCTi003-A'],[50,40],colors)):
 p=R/line/'fixed_K_d_Y_GLac';cr=rows(p/'curves.csv');ob=rows(p/'observations.csv')
 for j,state in enumerate(states):
  a=axs[i,j];sc=1e6 if j==0 else 1
  for strategy,co,ls in [('Unchanged ESI-017','#777777','--'),('Recalibrated',color,'-')]:
   rr=[v for v in cr if v['state']==state and v['strategy']==strategy];a.plot([float(v['day']) for v in rr],[float(v['value'])/sc for v in rr],color=co,ls=ls,lw=2.4)
  rr=[v for v in ob if v['state']==state];a.errorbar([float(v['day']) for v in rr],[float(v['mean'])/sc for v in rr],yerr=[float(v['sd'])/sc for v in rr],fmt='o',color=color,ms=5.5,capsize=4,elinewidth=1.8)
  end=max(float(v['day']) for v in rr);a.set_xlim(-.12,end+.3);a.set_xticks(range(int(end)+1));a.set_ylim(*lims[j]);a.set_xlabel('Time (days)');a.set_ylabel(labels[j]);a.set_title(f'({chr(97+4*i+j)})',loc='left',fontweight='bold',fontsize=17);a.grid(alpha=.14)
 y=axs[i,0].get_position().y1+.045
 f.add_artist(Rectangle((.075,y),.915,.038,transform=f.transFigure,color=color))
 f.text(.5325,y+.019,f'{line}  |  {rpm} RPM',ha='center',va='center',color='white',weight='bold',fontsize=18)
# Separate legends retain each cell line's curve and experimental colours.
for x,line,col in zip([.32,.72],['585A1','SCTi003-A'],colors):
 handles=[Line2D([],[],color='#777777',ls='--',lw=2.4),Line2D([],[],color=col,lw=2.4),Line2D([],[],color=col,marker='o',ls='',markersize=6)]
 f.legend(handles,['Prediction without recalibration','Recalibration 2','Experimental data ± SD'],loc='lower center',bbox_to_anchor=(x,.015),title=line,title_fontsize=17,fontsize=15,frameon=True)
for ext in ['png','pdf']:f.savefig(O/('FigureS7_iPSC_recalibration.'+ext),dpi=350,bbox_inches='tight')
plt.close(f)
fonts={str(v.get_object()['/BaseFont']) for p in PdfReader(O/'FigureS7_iPSC_recalibration.pdf').pages for v in p['/Resources']['/Font'].values()};assert all('Arial' in x for x in fonts)
