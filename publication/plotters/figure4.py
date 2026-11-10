from pathlib import Path
import csv,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
CACHE=ROOT/'publication/cache'
OUT=ROOT/'publication/generated/figures'
OUT.mkdir(parents=True,exist_ok=True)
R=CACHE/'training60'
O=OUT
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.labelsize':13,'axes.titlesize':14,'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':11,'pdf.fonttype':42,'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold','mathtext.sf':'Arial','mathtext.fallback':None})
def rows(p):return list(csv.DictReader(p.open()))
r=rows(R/'media/trajectories.csv');p=rows(R/'media/predictions.csv')
states=['cells','glucose','lactate','aggregate'];labels=[r'Viable cells ($\times 10^6$ cells/mL)','Glucose (mM)','Lactate (mM)','Aggregate size (µm)']
fig,axes=plt.subplots(2,2,figsize=(9,6.7))
for j,(a,st,label) in enumerate(zip(axes.flat,states,labels)):
 scale=1e6 if j==0 else 1
 for protocol,color in [('20','#005A70'),('50','#A94E32')]:
  x=[v for v in r if v['protocol_percent']==protocol and v['state']==st and v['model']=='5'];a.plot([float(v['day']) for v in x],[float(v['value'])/scale for v in x],color=color,lw=2,label=protocol+'% exchange')
  x=[v for v in p if v['protocol_percent']==protocol and v['state']==st and v['model']=='5'];a.errorbar([float(v['day']) for v in x],[float(v['observed'])/scale for v in x],yerr=[float(v['sd'])/scale for v in x],fmt='o',ms=5,capsize=4,elinewidth=1.7,color=color)
 a.set(xlabel='Time (days)',ylabel=label,title=f'({chr(97+j)})');a.grid(alpha=.13)
h,l=axes.flat[0].get_legend_handles_labels()
proxy=axes.flat[0].errorbar([],[],yerr=[],fmt='o',ms=5,capsize=4,elinewidth=1.7,color='#333333')
fig.legend(h+[proxy],l+['Experimental data ± SD'],loc='lower center',ncol=3)
fig.tight_layout(rect=[0,.06,1,1])
for ext in ['png','pdf']:fig.savefig(O/('Figure4_media_exchange.'+ext),dpi=350,bbox_inches='tight')
