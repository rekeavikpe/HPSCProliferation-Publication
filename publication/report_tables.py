"""Recalculate paper table values and compare their displayed rounding to the final tables."""
from pathlib import Path
import csv,json,math,re,statistics as st
from decimal import Decimal
P=Path(__file__).resolve().parent;R=P.parent;C=P/'cache';V=P/'verification';G=P/'generated/tables';G.mkdir(parents=True,exist_ok=True)
def read(p):
 with p.open(newline='')as f:return list(csv.DictReader(f))
paper=json.loads((P/'paper_tables.json').read_text())
def table(n):return paper['Table_'+str(n)+'_v9_7']['tables'][0]
checks=[];exports={}
def number(s):
 s=s.translate(str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻','0123456789-')).replace('−','-').strip()
 if '×'in s:
  a,b=s.split('×'); exponent=int(b.strip().removeprefix('10').replace('^',''));mant=a.strip()
 else:mant=s;exponent=0
 value=float(mant)*10**exponent
 tol=.500001*10**(Decimal(mant).as_tuple().exponent+exponent)
 return value,tol
def check(item,actual,display):
 val,tol=number(display);assert math.isfinite(actual) and abs(actual-val)<=tol+1e-9*max(abs(val),1), (item,actual,display,tol)
 checks.append({'item':item,'calculated':actual,'paper_display':display,'rounding_tolerance':tol})
def export(name,rows):
 rows=list(rows);exports[name]=rows
 if rows:
  with (G/(name+'.csv')).open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
primary=read(V/'primary_predictions.csv');scores=read(V/'primary_scores.csv')
for row in table(1)[1:]:
 m=int(row[0]);r=next(x for x in scores if int(x['manuscript_model'])==m and x['cohort']=='selection')
 check('Table1 Model'+str(m),float(r['bic']),row[-1])
export('Table1_selection_scores',[x for x in scores if x['cohort']=='selection'and x['state']=='all'])
params={r['parameter']:float(r['value'])for r in read(C/'training60/reference_fit/model5/selected_parameters.csv')}
names=['mu_max','K_d','K_oG','Y_GLac','alpha','K_Agg','K_tau','K_eps','K_G','K_Lac','Y_GX'];intervals={r['parameter']:r for r in read(C/'nine_profiles/audit.csv')};out=[]
for name,row in zip(names,table(2)[1:]):
 check('Table2 '+name,params[name],row[2]);z={'parameter':name,'estimate':params[name],'lower':'','upper':''}
 if name in intervals:
  z.update(lower=float(intervals[name]['lower']),upper=float(intervals[name]['upper']))
  lo,hi=row[3].split(' to ');check('Table2 '+name+' lower',z['lower'],lo);check('Table2 '+name+' upper',z['upper'],hi)
  assert z['lower']<=z['estimate']<=z['upper']
 out.append(z)
export('Table2_parameters',out)
def mrpe(rr,obs='observed'):
 return st.mean(100*abs(float(r['predicted'])-float(r[obs]))/abs(float(r[obs]))for r in rr if float(r[obs])!=0)
def variability(rr):return st.mean(100*float(r['sd'])/abs(float(r['observed']))for r in rr if float(r['observed'])!=0 and int(r['n_vessels'])>=2)
scale=read(C/'training60/reference_fit/figure_data/test_predictions.csv');media=read(C/'training60/media/predictions.csv')
sets=[('calibration',[r for r in primary if r['manuscript_model']=='6'and r['cohort']=='training']),('selection',[r for r in primary if r['manuscript_model']=='6'and r['cohort']=='selection']),('Model0_500mL',[r for r in scale if r['publication_model']=='0']),('Model6_500mL',[r for r in scale if r['publication_model']=='5']),('media20',[r for r in media if r['model']=='5'and r['protocol_percent']=='20']),('media50',[r for r in media if r['model']=='5'and r['protocol_percent']=='50'])]
state_names=['cells','glucose','lactate','aggregate','all'];out=[]
for j,(name,rr) in enumerate(sets):
 for i,state in enumerate(state_names):
  sub=rr if state=='all'else[r for r in rr if r['state']==state];v=mrpe(sub)
  check('Table3 '+name+' '+state,v,table(3)[i+1][j+1]);out.append({'dataset':name,'state':state,'n':len(sub),'mrpe_percent':v,'experimental_variability_percent':variability(sub)})
 if j>=2:check('Table3 variability '+name,variability(rr),table(3)[6][j+1])
export('Table3_prediction_errors',out)
series=read(R/'data/processed/series.csv');split=read(C/'training60/reference_fit/split.csv');out=[]
for c in range(1,7):
 rr=[r for r in series if int(r['condition'])==c];n=len(rr);check('TableS1 condition'+str(c),n,table('S1')[c][5]);out.append({'condition':c,'runs':n,'calibration':sum(r['role']=='calibration'and r['series_id']in{v['series_id']for v in rr}for r in split),'selection':sum(r['role']=='selection'and r['series_id']in{v['series_id']for v in rr}for r in split)})
export('TableS1_cohorts',out)
expo=read(R/'data/processed/exposure_summary_inputs.csv');conditions=[(100,30),(100,40),(100,60),(500,30)]
for sn,field in [('S2','shear_pa'),('S3','edr_m2_s3')]:
 out=[]
 for i,(vol,rpm)in enumerate(conditions):
  for j,size in enumerate([None,200,400]):
   row=table(sn)[1+i*3+j];z={'volume_ml':vol,'rpm':rpm,'size_um':size or 'Eulerian'}
   for k,stat in enumerate(['mean','p50','p75','p90','max']):
    scenario=('eul_va'if stat=='mean'else'eul_'+stat)if size is None else'lag_'+stat
    x=next(x for x in expo if x['scenario']==scenario and int(x['volume_ml'])==vol and int(x['rpm'])==rpm and (size is None or float(x['size_um'])==size))
    v=float(x[field]);check(sn+str((vol,rpm,size,stat)),v,row[2+k]);z[stat]=v
   out.append(z)
 export('Table'+sn,out)
meas=[r for r in read(R/'data/processed/measurements.csv')if r['status']=='measured'and r['value']]
out=[]
for row in table('S4')[1:]:
 c=int(row[0]);a=[float(r['value'])for r in meas if int(r['condition'])==c and r['state']=='aggregate'and float(r['day'])==5]
 mean,sd=row[2].split(' ± ');check('S4 mean '+str(c),st.mean(a),mean);check('S4 SD '+str(c),st.stdev(a),sd);check('S4 n '+str(c),len(a),row[3]);out.append({'condition':c,'n':len(a),'mean_um':st.mean(a),'sd_um':st.stdev(a)})
export('TableS4_day5',out)
oxygen=read(C/'october20/oxygen/condition_summary.csv');out=[]
for row in table('S5')[1:]:
 c=int(row[0]);r=next(x for x in oxygen if int(x['condition'])==c)
 groups={d:[float(x['value'])for x in meas if int(x['condition'])==c and x['state']=='cells'and float(x['day'])==d] for d in range(1,7)}
 peak=max(st.mean(a)for a in groups.values()if a);assert math.isclose(peak,float(r['peak_cells_ml']),rel_tol=1e-12)
 kla=float(r['kla_h']);da=17.7e-18*peak*1e6*3600/(.21*kla)
 check('S5 peak '+str(c),peak/1e6,row[2]);check('S5 kLa '+str(c),kla,row[3].split(' (')[0]);check('S5 Da '+str(c),da,row[4].split(' (')[0]);out.append({'condition':c,'peak_cells_ml':peak,'kla_h':kla,'Da':da})
 if c==4:
  low=17.7e-18*peak*1e6*3600/(.21*9);high=17.7e-18*peak*1e6*3600/(.21*.44)
  check('S5 Da lower',low,'0.03');check('S5 Da upper',high,'0.71')
export('TableS5_oxygen',out)
for row in table('S6')[1:]:
 m=row[0];rr=[r for r in primary if r['manuscript_model']==m and r['cohort']=='selection'];b=next(r for r in scores if r['manuscript_model']==m and r['cohort']=='selection');check('S6 BIC '+m,float(b['bic']),row[2])
 for j,state in enumerate(state_names):check('S6 '+m+' '+state,mrpe(rr if state=='all'else[r for r in rr if r['state']==state]),row[j+3])
p75=read(C/'october20/p75/selection_scores.csv')
for model,display in zip(table('S7')[0][1:],table('S7')[1][1:]):check('S7 Model'+model,float(next(r['bic']for r in p75 if r['tag']=='Model '+model)),display)
export('TableS7_p75',p75)
rr=[r for r in primary if r['manuscript_model']=='6'and r['cohort']=='selection'];out=[]
for c,display in zip([1,2,3,5,6],table('S8')[1][1:]):
 v=mrpe([r for r in rr if int(r['condition'])==c]);check('S8 C'+str(c),v,display);out.append({'condition':c,'mrpe_percent':v})
export('TableS8_selection_by_condition',out)
cv=read(V/'loco_metrics.csv')
for i,c in enumerate([1,2,3,5,6]):
 for j,state in enumerate(state_names):check('S9 '+str(c)+' '+state,float(next(r['mrpe_percent']for r in cv if int(r['manuscript_condition'])==c and r['state']==state)),table('S9')[i+1][j+1])
for j,state in enumerate(state_names):check('S9 mean '+state,st.mean(float(r['mrpe_percent'])for r in cv if r['state']==state),table('S9')[6][j+1])
export('TableS9_loco',cv)
transfer=read(C/'training60/cell_transfer/predictions.csv');refit=read(V/'ipsc_refit_scores.csv')
for j,line in enumerate(['585A1','SCTi003-A']):
 rr=[r for r in transfer if r['cell_line']==line and r['model']=='5'];check('S10 direct '+line,mrpe(rr),table('S10')[1][j+1])
 for i,label in enumerate(['Recalibration 1','Recalibration 2']):check('S10 '+line+' '+label,float(next(r['mrpe_percent']for r in refit if r['line']==line and r['strategy']==label)),table('S10')[i+2][j+1])
 for i,name in enumerate(['mu_max','K_oG','alpha']):
  ar=next(r for r in read(C/'ipsc'/line/'fixed_K_d_Y_GLac/profiles/audit.csv')if r['parameter']==name)
  text=table('S11')[i+1][j+2];mle,ci=text.split(' [');lo,hi=ci.rstrip(']').split('–')
  for tag,key,display in [('estimate','mle',mle),('lower','lower',lo),('upper','upper',hi)]:check('S11 '+line+' '+name+' '+tag,float(ar[key]),display)
export('TableS10_ipsc_refits',refit)
export('TableS6_selection_errors',[r for r in scores if r['cohort']=='selection'])
export('TableS11_ipsc_intervals',[dict(cell_line=line,**r) for line in ['585A1','SCTi003-A'] for r in read(C/'ipsc'/line/'fixed_K_d_Y_GLac/profiles/audit.csv')])
# Frozen transfer outputs must match the saved predictions after fresh ODE evaluation.
for name,key in [('media',('model','protocol_percent','state','day')),('cell_transfer',('model','cell_line','state','day'))]:
 saved=read(C/'training60'/name/'predictions.csv');new=read(P/'generated'/name/'predictions.csv');indexed={tuple(r[k]for k in key):r for r in saved};assert len(new)==len(saved)
 for r in new:
  old=indexed[tuple(r[k]for k in key)]
  for f in ['observed','sd','predicted']:assert math.isclose(float(r[f]),float(old[f]),rel_tol=1e-7,abs_tol=1e-6),(name,r,f,old[f])
(V/'table_rounding_checks.json').write_text(json.dumps({'status':'PASS','numeric_comparisons':len(checks),'checks':checks},indent=2))
print('PASS:',len(checks),'numeric comparisons against submitted tables; fresh frozen-transfer predictions match their saved values.')

# Reproduce the accepted descriptive statistics and final analytical oxygen scaling.
matched=[]
for rel in ['statistics/summary.csv','statistics/tests.csv','statistics/pairwise_rank_tests.csv','statistics/exposure_correlations.csv','oxygen/condition_summary.csv','oxygen/scaling_scenarios.csv']:
 old=read(C/'october20'/rel);new=read(P/'generated'/rel)
 assert len(old)==len(new),(rel,'row count')
 for a,b in zip(old,new):
  assert a.keys()==b.keys(),(rel,'columns')
  for k in a:
   if a[k]==b[k]:continue
   try:x,y=float(a[k]),float(b[k])
   except ValueError:raise AssertionError((rel,k,a[k],b[k]))
   assert (math.isnan(x) and math.isnan(y)) or math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-10),(rel,k,x,y)
 matched.append(rel)
(V/'statistics_oxygen_checks.json').write_text(json.dumps({'status':'PASS','matched_outputs':matched},indent=2))
print('PASS: reproduced statistical summaries, tests and final oxygen scaling.')
