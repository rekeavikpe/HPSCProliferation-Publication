from pathlib import Path
import csv,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
CACHE=ROOT/'publication/cache'
OUT=ROOT/'publication/generated/figures'
OUT.mkdir(parents=True,exist_ok=True)
W=ROOT;O=ROOT/'publication/generated';S=OUT
# Same workbook summaries as the final figure, supplied as curated CSV.
a=list(csv.DictReader((ROOT/'data/processed/exposure_summary_inputs.csv').open()))
e=list(csv.DictReader((CACHE/'october20/reuse/exposure_rpm_sweep.csv').open()));e=[r for r in e if r['scenario']=='eul_va'];rpms=[20,30,40,60,80]
plt.rcParams.update({'font.family':'Arial','font.size':16,'axes.labelsize':17,'axes.titlesize':17,'xtick.labelsize':15,'ytick.labelsize':15,'legend.fontsize':14,'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','pdf.fonttype':42})
values=[]
for field,col,name,label in [('shear_pa','Shear','FigureS2_shear_rpm','Shear stress (Pa)'),('edr_m2_s3','EDR','FigureS3_edr_rpm',r'EDR (m$^2$ s$^{-3}$)')]:
 f,axs=plt.subplots(1,2,figsize=(13,5.8));f.subplots_adjust(left=.09,right=.98,top=.91,bottom=.30,wspace=.32)
 for ax,stat,letter in zip(axs,['Mean','75th Percentile'],['a','b']):
  y=[float(next(r[field]for r in e if int(r['rpm'])==rpm))for rpm in rpms];ax.plot(rpms,y,'o-',color='black',ms=7,lw=2.1,label='Eulerian volume average')
  for size,color,marker in [(100,'#0072B2','s'),(200,'#D55E00','o'),(300,'#009E73','^')]:
   yy=[]
   for rpm in rpms:
    scenario='lag_mean' if stat=='Mean' else 'lag_p75'
    y200=next(float(r[field])for r in a if int(r['volume_ml'])==100 and int(r['rpm'])==rpm and float(r['size_um'])==200 and r['scenario']==scenario);y400=next(float(r[field])for r in a if int(r['volume_ml'])==100 and int(r['rpm'])==rpm and float(r['size_um'])==400 and r['scenario']==scenario)
    yv=y200+(size-200)/200*(y400-y200);assert yv>=0;yy.append(yv);values.append(dict(field=field,summary=stat,rpm=rpm,size_um=size,value=yv))
   ax.plot(rpms,yy,marker=marker,ms=7,lw=2,color=color,label=f'{size} µm'+(' (extrapolated)'if size==100 else' (interpolated)'if size==300 else''))
  ax.set(xlabel='Agitation (RPM)',ylabel=label,title=f'({letter}) Lagrangian '+('mean'if stat=='Mean'else'P75'),xticks=rpms,xlim=(17,83));ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True);ax.grid(alpha=.15);ax.margins(y=.12)
 f.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=2,bbox_to_anchor=(.53,.005))
 for ext in ['png','pdf']:f.savefig(S/(name+'.'+ext),dpi=350,bbox_inches='tight')
 plt.close(f)
with(O/'sweep_values.csv').open('w')as f:
 w=csv.DictWriter(f,fieldnames=values[0]);w.writeheader();w.writerows(values)
print('RPM figures regenerated from rebuilt particle summaries.')
