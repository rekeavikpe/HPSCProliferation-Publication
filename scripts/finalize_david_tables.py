#!/usr/bin/env python3
"""Read saved results only. No fitting, parameter changes or model reselection."""
import csv, html, json, sys
from pathlib import Path
root=Path(sys.argv[1] if len(sys.argv)>1 else 'results/david_revision_20261014').resolve()
out=root/'writing_tables';out.mkdir(exist_ok=True)
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def table(rows,cols,bold=None):
 return '<table><tr>'+''.join('<th>'+html.escape(label)+'</th>' for key,label in cols)+'</tr>'+''.join('<tr>'+''.join('<td>'+('<b>' if bold and bold(r,k) else '')+html.escape(str(r.get(k,'')))+('</b>' if bold and bold(r,k) else '')+'</td>' for k,label in cols)+'</tr>' for r in rows)+'</table>'
sections=[];status=[]
for stage in ('selection','exposure_selection'):
 folder=root/stage
 if not (folder/'COMPLETE').exists():status.append({'stage':stage,'status':'pending'});continue
 ranks=read(folder/'ranking.csv');metrics=read(folder/'state_errors.csv');assembled=[]
 for r in ranks:
  row={'model':('Model '+r['publication_model']) if stage=='selection' else r['scenario']}
  for state in ('cells','glucose','lactate','aggregate','all'):
   matches=[m for m in metrics if m['scenario']==r['scenario'] and m['cohort']=='selection40' and m['state']==state]
   if len(matches)!=1:raise ValueError(f'Missing or duplicate metric {stage}/{r["scenario"]}/{state}')
   row[state]=f'{float(matches[0]["mrpe_percent"]):.1f}'
   row[state+'_raw']=float(matches[0]['mrpe_percent'])
  row['k']=r['k'];assembled.append(row)
 # Rows are not labelled a new model selection ranking.
 if stage=='selection':assembled.sort(key=lambda r:int(r['model'].split()[-1]))
 thresholds={s:sorted(r[s+'_raw'] for r in assembled)[min(2,len(assembled)-1)] for s in ('cells','glucose','lactate','aggregate','all')}
 cols=[('model','Model' if stage=='selection' else 'Exposure summary'),('k','Fitted parameters'),('cells','Cells'),('glucose','Glucose'),('lactate','Lactate'),('aggregate','Aggregate size'),('all','Overall')]
 with (out/(stage+'_heldout_errors.csv')).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=[k for k,l in cols],extrasaction='ignore');w.writeheader();w.writerows(assembled)
 sections.append('<h2>'+('Model comparison' if stage=='selection' else 'Exposure-summary comparison')+'</h2><p>MRPE (%) on the withheld 40% of 100 mL replicates. Parameters are estimated using only the calibration 60%. Bold identifies the three lowest errors in each column, including exact ties. These are prediction errors, not calibration BIC.</p>'+table(assembled,cols,lambda r,k:k in thresholds and r[k+'_raw']<=thresholds[k]))
 status.append({'stage':stage,'status':'exported'})
cv=root/'pooled_validation'/'replicate_error_summary.csv'
if cv.exists() and (root/'pooled_validation'/'COMPLETE').exists():
 rows=[r for r in read(cv) if r['group']=='all_groups']
 sections.append('<h2>Pooled validation and exclusion sensitivity</h2><p>Mean and variance across individual replicate MRPEs. Variance units are percentage points squared. These are distinct from errors calculated against condition means. Each withheld replicate appears once within each pooled design.</p>'+table(rows,[(k,k.replace('_',' ')) for k in ('experiment','model','state','n_replicates','mean_mrpe_percent','variance_mrpe_percentage_points_squared')]))
 status.append({'stage':'pooled_validation','status':'exported'})
else:status.append({'stage':'pooled_validation','status':'pending'})
rpm=root/'rpm_statistics'/'welch_holm.csv'
if rpm.exists():
 rows=read(rpm)
 for matched in (True,False):
  rs=[r for r in rows if (r['protocol_matched'].lower()=='true')==matched]
  sections.append('<h2>RPM comparisons: '+('matched recorded protocols' if matched else 'different recorded handling')+'</h2><p>'+('Experimental comparisons with matched volume, initial density and recorded protocols.' if matched else 'These comparisons do not isolate RPM because recorded handling differs.')+'</p>'+table(rs,[(k,k.replace('_',' ')) for k in ('state','endpoint','rpm_a','rpm_b','n_a','n_b','mean_difference','p_value','p_holm','note')]))
 status.append({'stage':'rpm_statistics','status':'exported'})
(out/'status.json').write_text(json.dumps(status,indent=2))
(out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Tables for David response</title><style>body{font:16px Arial;line-height:1.5;margin:30px}table{border-collapse:collapse;margin-bottom:30px}th,td{border-bottom:1px solid #bbb;padding:8px;text-align:right}th:first-child,td:first-child{text-align:left}h2{margin-top:40px}b{font-weight:800}</style><h1>Results for response drafting</h1><p>Saved outputs only. No optimization or model reselection performed. Original ranking files remain unchanged. The held-out error tables address the request to compare errors on data excluded from calibration. Calibration BIC must not be labelled validation BIC.</p>'+''.join(sections))
print(out/'index.html')
print(json.dumps(status,indent=2))
