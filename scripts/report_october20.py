import csv,html,json,pathlib,sys
root=pathlib.Path(sys.argv[1]);sections=[]
for p in sorted(root.rglob('*.csv')):
    if any(x.startswith('run') and x[3:].isdigit() for x in p.parts) or 'profiles' in p.parts: continue
    if p.name not in ['selection_scores.csv','candidate_selection_scores.csv','scores.csv','ablation_scores.csv','model6_errors.csv','tests.csv','summary.csv','condition_summary.csv','aggregate_core_sensitivity.csv','final_rounds.csv','parameter_policy.csv','errors.csv','fold_predictions.csv','pairwise_rank_tests.csv','exposure_correlations.csv','shear_inhibition_summary.csv','scaling_scenarios.csv']: continue
    rows=list(csv.reader(p.open()))
    sections.append('<h2>'+html.escape(str(p.relative_to(root)))+'</h2><table>'+''.join('<tr>'+''.join('<'+('th' if i==0 else 'td')+'>'+html.escape({"baseline":"Model 0","full":"Model 6"}.get(c,c))+'</'+('th' if i==0 else 'td')+'>' for c in row)+'</tr>' for i,row in enumerate(rows))+'</table>')
status=json.loads((root/'cells/status.json').read_text())
review=[r for r in status if r['status']=='profile_numerical_review_required']
reductions=[r for r in status if r.get('unbounded_parameters')]
notes='Five-parameter fits and profiles are reported. No parameters were automatically removed. Numerical or parameter-reduction review remains required.' if review or reductions else 'All five requested parameter profiles have bounded intervals under the stated error model.'
(root/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>David October 20 analyses</title><style>body{font:16px Arial;max-width:1250px;margin:40px auto}table{border-collapse:collapse}td,th{padding:7px;border:1px solid #ccc}img{max-width:100%}</style><h1>David’s October 20 requests</h1><p>'+notes+'</p><p>The original 19/12 split is retained. Day 1 scores exclude the measured initial states. Cell-line errors are fit errors. The 500 mL oxygen coefficient is a conditional scaling estimate. Computation completion does not establish the draft’s proposed conclusions.</p>'+''.join('<h2>'+p.stem+'</h2><img src="'+str(p.relative_to(root))+'">' for p in sorted((root/'figures').glob('*.png')))+''.join(sections))
(root/'review_status.json').write_text(json.dumps({'computations':'complete','cell_profile_numerical_flags':review,'cell_parameter_reduction_review':reductions,'oxygen_geometry_and_surface_CFD':'not independently verified, scaling assumptions documented','manuscript':'wording and placeholders require reconciliation against these outputs'},indent=2))
print(root/'index.html')
