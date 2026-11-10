#!/usr/bin/env python3
"""Reporting only. Read completed CSVs without changing any fit or raw output."""
import csv, html, math, sys
from pathlib import Path
root=Path(sys.argv[1] if len(sys.argv)>1 else 'results/david_revision_20261014').resolve()
out=root/'reporting_corrected';out.mkdir(exist_ok=True)
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def save(name,rows):
 if not rows:return
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def table(rows):
 if not rows:return '<p>No saved results available.</p>'
 def fmt(v):return f'{v:.3g}' if isinstance(v,float) else str(v)
 return '<table><tr>'+''.join('<th>'+html.escape(k.replace('_',' '))+'</th>' for k in rows[0])+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(fmt(v))+'</td>' for v in r.values())+'</tr>' for r in rows)+'</table>'
parts=[]
for stage in ('selection','exposure_selection'):
 ranks=read(root/stage/'ranking.csv');metrics=read(root/stage/'state_errors.csv');result=[]
 for r in ranks:
  m={x['state']:x for x in metrics if x['scenario']==r['scenario'] and x['cohort']=='selection40'}
  n=int(m['all']['n']);k=int(r['k']);nll=float(r['selection_nll'])
  assert n>0
  result.append(dict(model=('Model '+r['publication_model']) if stage=='selection' else r['scenario'],n_selection=n,k=k,selection_nll=nll,selection_BIC_form_score=2*nll+k*math.log(n),**{s+'_MRPE_percent':float(m[s]['mrpe_percent']) for s in ('cells','glucose','lactate','aggregate','all')}))
 result.sort(key=lambda r:r['selection_BIC_form_score']);save(stage+'_scores.csv',result)
 parts.append('<h2>'+stage.replace('_',' ').title()+'</h2><p>The BIC-form score is 2 × selection negative log-likelihood + k × log(n selection). Parameters remain those fitted on the calibration data. This is a penalized held-out score, not standard BIC evaluated at an estimate fitted to those same observations. n counts the condition–day–state means scored by the likelihood, not experimental replicates. No model is refitted or automatically reselected.</p>'+table(result))
profiles=[]
for summary in sorted((root/'separate_ablations/profiles').glob('*/*/profile_confidence_intervals.csv')):
 for r in read(summary):
  points=read(summary.parent/r['parameter']/'profile.csv');mle=float(r['mle']);row=dict(mode=summary.parent.parent.name,formulation=summary.parent.name,parameter=r['parameter'],MLE=mle)
  for side in ('lower','upper'):
   pts=[x for x in points if (float(x['parameter_value'])<mle if side=='lower' else float(x['parameter_value'])>mle)]
   crossed=any(float(x['profile_loglik'])<=-1.920729410347062 for x in pts)
   endpoint=float(r['package_'+side]);ordered=endpoint<mle if side=='lower' else endpoint>mle
   valid=crossed and ordered
   row[side+'_limit']=endpoint if valid else ''
   row[side+'_status']='threshold crossed' if valid else ('endpoint requires review' if crossed else 'unresolved within tested range')
   row[side+'_tested_extent']=(min if side=='lower' else max)(float(x['parameter_value']) for x in points)
  row['better_minimum']=r['better_minimum'];profiles.append(row)
save('profile_intervals_for_reporting.csv',profiles)
parts.append('<h2>Profile interval reporting</h2><p>Only sides with a threshold crossing in the saved profile and an endpoint on the correct side of the MLE receive a numerical limit. An unresolved limit is not replaced by an imposed bound or the MLE. The original ProfileLikelihood.jl files and figures remain unchanged. Fixed-parameter profiles retain their fixed-parameter interpretation.</p>'+table(profiles))
cv=[r for r in read(root/'pooled_validation/replicate_error_summary.csv') if r['group']=='all_groups']
cv=[{k:(float(v) if k in ('mean_mrpe_percent','variance_mrpe_percentage_points_squared') else v) for k,v in r.items() if k in ('experiment','model','state','n_replicates','mean_mrpe_percent','variance_mrpe_percentage_points_squared')} for r in cv]
save('pooled_validation_mean_variance.csv',cv)
parts.append('<h2>Pooled cross-validation and exclusion sensitivity</h2><p>Mean and sample variance of errors across held-out experimental replicates. Variance is in percentage points squared. These are not errors against condition means. The exclusion analysis ranks training replicates by final viable-cell density in eligible experimental groups. Held-out replicates are unchanged.</p>'+table(cv))
cell=read(root/'cell_lines/state_errors.csv')
for r in cell:
 r['assessment']='prediction without refitting' if r['strategy'].startswith('transfer_') else 'cell-line calibration'
 r['mrpe_percent']=float(r['mrpe_percent']);r['experimental_noise_percent']=float(r.pop('experimental_noise_percent'))
save('cell_line_errors_and_noise.csv',cell)
parts.append('<h2>Cell-line comparisons</h2><p>Refitting results are calibration errors. Unchanged ESI-017 parameter predictions are reported separately. Experimental noise is the mean experimental SD divided by its corresponding mean, expressed as a percentage.</p>'+table(cell))
ab=read(root/'separate_ablations/state_errors.csv')
for r in ab:r['cohort']=r['cohort'].replace('retrospective ','')
save('ablation_errors_and_noise.csv',ab)
parts.append('<h2>Ablations</h2><p>Independent response parameters in the both-equation variants. Calibration uses all 100 mL data. The 500 mL errors use predictions without refitting.</p>'+table(ab))
parts.append('<h2>Remaining response work</h2><p>Pouyan’s oxygen CFD is not supplied in this bundle. Final manuscript parameter, units and figure consistency must be checked against the manuscript. The held-out BIC-form score must retain its definition. It must not be described as standard fitted-data BIC.</p>')
(out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Corrected reporting</title><style>body{font:16px Arial;margin:30px;line-height:1.5}table{border-collapse:collapse;margin-bottom:30px}td,th{border:1px solid #ccc;padding:7px;text-align:left}h2{margin-top:40px}</style><h1>Reporting from completed David revision analyses</h1><p><a href="figures_report.html">Regenerated figures</a></p><p>Saved results only. No calibration, optimization or parameter changes.</p>'+''.join(parts))
print(out/'index.html')
print(f'Profile parameters checked: {len(profiles)}')
