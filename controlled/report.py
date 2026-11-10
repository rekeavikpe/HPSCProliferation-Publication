"""Summarize the controlled runs without treating search endpoints as confidence limits."""
import csv
import html
import sys
from pathlib import Path

root = Path(sys.argv[1])
def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))

def escape(value):
    return html.escape(str(value))

summary = []
sections = []
for case in rows(root / 'design.csv'):
    folder = root / case['name']
    if not (folder / 'COMPLETE').exists():
        raise RuntimeError(f"Incomplete case: {case['name']}")
    estimates = {r['parameter']: float(r['value']) for r in rows(folder / 'mle_parameters.csv')}
    diagnostics = rows(folder / 'profile_confidence_intervals.csv')
    fit = rows(folder / 'mle_summary.csv')[0]
    for parameter in ('K_G', 'mu_max'):
        points = rows(folder / parameter / 'profile.csv')
        mle = estimates[parameter]
        lower = any(float(p['parameter_value']) < mle and float(p['profile_loglik']) <= -1.920729410347062 for p in points)
        upper = any(float(p['parameter_value']) > mle and float(p['profile_loglik']) <= -1.920729410347062 for p in points)
        diagnostic = next(r for r in diagnostics if r['parameter'] == parameter)
        summary.append(dict(case=case['name'], parent=case['parent'], change=case['change'], parameter=parameter, mle=mle,
                            lower_threshold_crossed=lower, upper_threshold_crossed=upper,
                            package_lower=diagnostic['package_lower'], package_upper=diagnostic['package_upper'],
                            unconfirmed_points=diagnostic['unconfirmed_points'], better_minimum=diagnostic['better_minimum'],
                            mle_stopping_met=fit['stopping_met']))
    warning = '<p><strong>K_G has no upper threshold crossing within the tested range. The native title’s upper endpoint is not a verified confidence limit.</strong></p>' if any(r['case']==case['name'] and r['parameter']=='K_G' and not r['upper_threshold_crossed'] for r in summary) else ''
    sections.append(f"<section><h2>{escape(case['name'])}</h2><p>{escape(case['change'])}. Parent: {escape(case['parent'] or 'original manuscript setup')}.</p>{warning}<a href='{case['name']}/figures/model5_profiles.pdf'>PDF</a><img src='{case['name']}/figures/model5_profiles.png' alt='Native package profiles for {escape(case['name'])}'></section>")
with (root / 'comparison.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
    writer.writeheader()
    writer.writerows(summary)
body = ''.join('<tr>' + ''.join(f'<td>{escape(value)}</td>' for value in (r['case'], r['parameter'], f"{r['mle']:.6g}", 'yes' if r['lower_threshold_crossed'] else 'no', 'yes' if r['upper_threshold_crossed'] else 'no', r['unconfirmed_points'])) + '</tr>' for r in summary)
(root / 'index.html').write_text('''<!doctype html><meta charset="utf-8"><title>Controlled identifiability comparison</title><style>body{max-width:1100px;margin:40px auto;font:17px/1.5 system-ui;color:#18232b;padding:0 24px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:8px;border-bottom:1px solid #ddd;text-align:left}img{width:100%;height:auto}section{margin-top:48px}a{color:#086c9c}.note{background:#edf4f7;padding:18px}</style><h1>Why did the profiles change?</h1>
<p>Original Model 6 (current Model 5), original physical-parameter BOBYQA multistart calibration, and native ProfileLikelihood.jl profiling and plots. Only K<sub>G</sub> and μ<sub>max</sub> are profiled; every remaining free parameter is reoptimized.</p>
<div class="note"><strong>Interpretation:</strong> A_original reproduces the published fixed-three setup. A_fullrange changes only the plotting/profile search window. B_100mL removes the 500-mL condition from the original dataset. All remaining cases branch separately from B_100mL. C changes all three fixed values; H–J separate those changes. These historical fixed values are diagnostic controls, not a leakage-free 100-mL calibration.</div>
<h2>Interpretation</h2><p>Compare each case with its recorded parent. H–J isolate changes to the fixed lactate, yield and aggregate parameters. A closed conditional profile does not establish identifiability when the fixed values are also uncertain. Historical combined-data fixed values must not be selected simply to produce narrower intervals. The μ<sub>max</sub> profiles in cases with K<sub>G</sub> ≤ 40 remain conditional on that bound.</p><h2>Profile checks</h2><p>“Crossed” means a computed point on that side reaches the nominal 95% likelihood-ratio reference (−1.92073). It does not prove a global optimum, account for uncertainty in fixed values, or validate boundary asymptotics. Native package titles are retained, but a returned endpoint without a threshold crossing is not a confidence limit.</p><table><tr><th>Case</th><th>Parameter</th><th>MLE</th><th>Lower crossed</th><th>Upper crossed</th><th>Unconfirmed points</th></tr>''' + body + '</table><p><a href="comparison.csv">Detailed CSV</a> · <a href="design.csv">Comparison design</a></p>' + ''.join(sections) + '<p>The original-data controls retain original feeding and initial-glucose assumptions. This comparison does not yet isolate curated-data changes, interactions between changes, or replacement of the optimizer. Raw likelihoods across datasets are not comparable.</p>')
print(root / 'index.html')
