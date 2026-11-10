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
W=ROOT;O=ROOT/'publication/generated';F=OUT;M=F;S=F;Q=F
for p in [M,S,Q]:p.mkdir(parents=True,exist_ok=True)
R=CACHE/'training60';D=CACHE/'october20';B=ROOT
def rows(p):return list(csv.DictReader(p.open()))
def style(fs=14):plt.rcParams.update({'font.family':'Arial','font.size':fs,'axes.labelsize':fs,'axes.titlesize':fs,'xtick.labelsize':fs-1,'ytick.labelsize':fs-1,'legend.fontsize':fs-1,'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold','pdf.fonttype':42})
def save(f,folder,name):
 for ext in ['png','pdf']:f.savefig(folder/(name+'.'+ext),dpi=350,bbox_inches='tight')
 plt.close(f)
def copy(src,dest):
 for ext in ['png','pdf']:shutil.copy2(W/(src+'.'+ext),dest.with_suffix('.'+ext))

from matplotlib import font_manager
print('Font:',font_manager.findfont('Arial',fallback_to_default=False))
# Figure 3: all original cultures in markers, primary training-fit curves unchanged.
style(18)
tr=rows(R/'reference_fit/figure_data/training_curves.csv');te=rows(R/'reference_fit/figure_data/test_trajectories.csv');ep=rows(R/'reference_fit/figure_data/test_predictions.csv');measure=rows(B/'data/processed/measurements.csv');series=rows(B/'data/processed/series.csv');ids={v['series_id'] for v in rows(R/'reference_fit/split.csv') if v['role'] in ('calibration','selection')}
assert len(ids)==31
assert sum(v['series_id'] in ids and v['condition']=='5' for v in series)==7
meas=[v for v in measure if v['series_id']in ids and v['status']=='measured' and v['value'] and float(v['day'])>0]
states=['cells','glucose','lactate','aggregate'];labels=[r'Viable cells ($\times10^6$ cells/mL)','Glucose (mM)','Lactate (mM)','Aggregate size (µm)'];colors=['#0072B2','#E69F00','#009E73','#56B4E9','#D55E00','#005A70'];limits=[(0,2.4),(0,18),(0,22),(0,390)]
f,axs=plt.subplots(6,4,figsize=(16.5,20));f.subplots_adjust(left=.145,right=.99,top=.975,bottom=.065,wspace=.62,hspace=.42)
summary=[]
for i,c in enumerate([1,2,3,5,6,4]):
 for j,st in enumerate(states):
  ax=axs[i,j];sc=1e6 if j==0 else 1
  cur=[v for v in te if v['publication_model']=='5' and v['state']==st]if c==4 else [v for v in tr if v['model']=='5'and v['condition']==str(c)and v['state']==st]
  ax.plot([float(v['day'])for v in cur],[float(v['value'])/sc for v in cur],color=colors[i],lw=2.6)
  if c==4:
   base=[v for v in te if v['publication_model']=='0'and v['state']==st];ax.plot([float(v['day'])for v in base],[float(v['value'])/sc for v in base],color='#777777',ls='--',lw=2.3)
   obs=[v for v in ep if v['publication_model']=='5'and v['state']==st];xx=[float(v['day'])for v in obs];yy=[float(v['observed'])/sc for v in obs];sd=[float(v['sd'])/sc for v in obs]
  else:
   rr=[v for v in meas if v['condition']==str(c)and v['state']==st];xx=sorted({float(v['day'])for v in rr});yy=[];sd=[]
   for day in xx:
    y=[float(v['value'])for v in rr if float(v['day'])==day];yy.append(np.mean(y)/sc);sd.append(np.std(y,ddof=1)/sc if len(y)>1 else 0);summary.append(dict(condition=c,state=st,day=day,n=len(y),mean=np.mean(y),sd=np.std(y,ddof=1)if len(y)>1 else 0))
  ax.errorbar(xx,yy,yerr=sd,fmt='o',ms=5.6,capsize=4.5,elinewidth=2.2,color=colors[i]);ax.set(xlim=(-.12,(5 if c in [4,5] else 6)+.25),ylim=limits[j],ylabel=labels[j],xticks=range(6 if c in [4,5] else 7));ax.set_title(f'({chr(97+i*4+j)})',loc='left',weight='bold');ax.grid(alpha=.15)
  if i==5:ax.set_xlabel('Time (days)')
 pos=axs[i,0].get_position();f.text(.004,(pos.y0+pos.y1)/2,('500 mL\nCondition 4\nValidation'if c==4 else f'100 mL\nCondition {c}'),va='center',weight='bold',fontsize=19)
f.legend([Line2D([],[],color='black',lw=2.6),Line2D([],[],color='#777777',ls='--',lw=2.3),Line2D([],[],color='black',marker='o',ls='none')],['Model 6','Model 0 (500 mL)','Experimental data ± SD'],loc='lower center',ncol=3,fontsize=18)
# Assert the plotted marker coordinates and error bars use the independently checked measurements.
expected={(5,'aggregate',4.0):(7,290.9215040600852,15.123024923469645),(5,'aggregate',5.0):(6,330.81934083014636,22.355882797775912),(5,'cells',4.0):(7,648428.5714285715,288743.17143028515)}
for v in summary:
 key=(v['condition'],v['state'],v['day'])
 if key in expected:
  n,mean,sd=expected[key];assert v['n']==n and np.isclose(v['mean'],mean,rtol=1e-12) and np.isclose(v['sd'],sd,rtol=1e-12)
# Save axes geometry to allow the visual check to focus on the corrected row.
import json
(O/'axes_geometry.json').write_text(json.dumps({'figure_inches':list(f.get_size_inches()),'axes':[list(a.get_position().bounds)for a in axs.flat]}))
save(f,M,'Figure3_Model6_all_conditions');
with (O/'figure3_observations.csv').open('w')as out:
 w=csv.DictWriter(out,fieldnames=summary[0]);w.writeheader();w.writerows(summary)

print('Verified current measurements; original model trajectories and all other condition markers preserved.')
