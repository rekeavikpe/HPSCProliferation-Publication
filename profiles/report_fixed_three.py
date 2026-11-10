"""Present the verified original-formulation fixed-three profile analysis."""
import csv
import html
import sys
from pathlib import Path
root=Path(sys.argv[1])
def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))
fixed=read(root/'fixed_parameters.csv')
checks=[r for r in read(root/'review/profile_checks.csv') if r['case']==root.name]
fits=[r for r in read(root/'review/fit_comparison.csv') if r['case']==root.name]
if len(checks)!=8 or len(fits)!=1 or {r['parameter'] for r in fixed}!={'K_G','K_Lac','Y_GX'}:
    raise ValueError('Expected the fixed-three original-formulation analysis')
closed=sum(r['lower_crossed']=='true' and r['upper_crossed']=='true' for r in checks)
converged=all(r['unconfirmed_points']=='0' for r in checks)
fixed_rows=''.join(f"<tr><td>{html.escape(r['parameter'])}</td><td>{float(r['value']):.9g}</td><td>{html.escape(r['source'])}</td></tr>" for r in fixed)
profile_rows=''.join('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in (r['parameter'],f"{float(r['mle']):.6g}",f"{float(r['package_lower']):.6g}" if r['lower_crossed']=='true' else 'unresolved',f"{float(r['package_upper']):.6g}" if r['upper_crossed']=='true' else 'unresolved',r['unconfirmed_points'],r['supported_nuisance_boundaries'] or 'none'))+'</tr>' for r in checks)
(root/'index.html').write_text('''<!doctype html><meta charset="utf-8"><title>Original formulation: three fixed parameters</title><style>body{max-width:1300px;margin:40px auto;padding:0 24px;font:17px/1.5 system-ui;color:#182c37}table{border-collapse:collapse;width:100%;font-size:15px}td,th{padding:9px;border-bottom:1px solid #ddd;text-align:left}img{width:100%;height:auto}.note{padding:18px;background:#edf3f8}a{color:#076d9e}</style><h1>Original formulation: three fixed parameters</h1><p>The finite lactate-inhibition factor K<sub>Lac</sub>/(L + K<sub>Lac</sub>) remains in the equations. K<sub>G</sub>, K<sub>Lac</sub> and Y<sub>GX</sub> are fixed at their exact saved 100-mL calibration estimates. The other eight parameters are refitted and profiled using ProfileLikelihood.jl. No bounds were tightened.</p>'''+f'<p><strong>{closed}/8 profiles cross the nominal 95% reference on both sides.</strong> All saved profile optimizations report success: {converged}. Recomputed NLL: {float(fits[0]["nll"]):.9f}.</p>'+'''<div class="note">These intervals are conditional on the three fixed estimates, all obtained from the same data. K<sub>Lac</sub> is a constrained estimate at its upper bound, not a separately identified biological constant. Preserve that qualification and the earlier open profiles in the supplementary analysis. The inactive-inhibition branch remains a separate sensitivity analysis.</div><h2>Fixed values</h2><table><tr><th>Parameter</th><th>Value</th><th>Source</th></tr>'''+fixed_rows+'''</table><h2>Native package profiles</h2><a href="figures/model5_profiles.pdf">PDF</a> · <a href="figures/model5_profiles.png">PNG</a><img src="figures/model5_profiles.png" alt="Eight conditional profile likelihoods"><h2>Verification</h2><table><tr><th>Parameter</th><th>MLE</th><th>Lower limit</th><th>Upper limit</th><th>Unconfirmed points</th><th>Supported nuisance boundaries</th></tr>'''+profile_rows+'''</table><p><a href="review/prediction_errors.csv">MRPE by state variable</a> · <a href="review/profile_checks.csv">Profile checks</a> · <a href="review/fit_comparison.csv">Fit comparison</a></p>''')
print(root/'index.html')
