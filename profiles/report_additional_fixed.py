"""Build the comparison page from verified profile and fit tables."""
import csv
import html
import sys
from pathlib import Path
root=Path(sys.argv[1])
def read(name):
    with (root/name).open() as f:
        return list(csv.DictReader(f))
checks=read('profile_checks.csv'); fits=read('fit_comparison.csv'); errors=read('prediction_errors.csv')
labels={'inactive_yield_profiles':'Nine free parameters','inactive_fixed_KAgg':'Fix K_Agg at its parent MLE','inactive_fixed_KG':'Fix K_G at its parent MLE'}
content=[]
for fit in fits:
    case=fit['case']; rows=[r for r in checks if r['case']==case]
    closed=sum(r['lower_crossed']=='true' and r['upper_crossed']=='true' for r in rows)
    qualified=all(r['unconfirmed_points']=='0' and r['better_minimum']=='false' for r in rows)
    bounds=sorted(set(r['supported_nuisance_boundaries'] for r in rows if r['supported_nuisance_boundaries']))
    content.append(f"<tr><td>{html.escape(labels.get(case,case))}</td><td>{closed}/{len(rows)}</td><td>{float(fit['nll']):.6f}</td><td>{'yes' if qualified else 'review flags'}</td><td>{html.escape(', '.join(bounds) or 'none detected')}</td></tr>")
figures=''.join(f"<h2>{html.escape(labels.get(r['case'],r['case']))}</h2><a href='../{r['case']}/figures/model5_profiles.pdf'>PDF</a><img src='../{r['case']}/figures/model5_profiles.png' alt='Native package profiles'>" for r in fits)
error_table=''.join('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in (labels.get(r['case'],r['case']),r['cohort'],r['state'],f"{float(r['mrpe']):.4f}"))+'</tr>' for r in errors)
(root/'index.html').write_text('''<!doctype html><meta charset="utf-8"><title>Additional fixed-parameter comparison</title><style>body{max-width:1250px;margin:40px auto;padding:0 24px;font:17px/1.5 system-ui;color:#162c3a}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:9px;border-bottom:1px solid #ddd;text-align:left}img{width:100%;height:auto}h2{margin-top:38px}.note{padding:18px;background:#eaf2f7}a{color:#076d9e}</style><h1>Eight-parameter conditional profiles</h1><p>Inactive lactate inhibition, Y_GX held at the saved 100-mL MLE. Each additional fixation uses the nine-parameter inactive model’s MLE. All other free parameters are refitted and profiled using ProfileLikelihood.jl. No bounds are tightened.</p><div class="note"><strong>Interpretation:</strong> Closing conditional profiles after fixing K_G addresses its compensation with μ_max. It does not establish practical identifiability of the original nine-parameter model or make K_G known independently. Report the fixed values and their provenance, retain the nine-parameter diagnostics in the supplement, and retain uncertainty from that parent model when assessing sensitivity to K_G.</div><h2>Verified comparison</h2><table><tr><th>Analysis</th><th>Two-sided threshold crossings</th><th>NLL</th><th>All saved points converged; no better optimum found</th><th>Nuisance boundaries within supported sampled profiles</th></tr>'''+''.join(content)+'''</table><p>Bound checks cover the sampled profile paths within the nominal 95% likelihood-ratio reference. They are not proof of a global optimum. <a href="profile_checks.csv">Full checks</a> · <a href="fit_comparison.csv">Fits and fixed values</a></p>'''+figures+'<h2>Prediction errors (%)</h2><table><tr><th>Analysis</th><th>Cohort</th><th>State</th><th>MRPE</th></tr>'+error_table+'</table>')
print(root/'index.html')
