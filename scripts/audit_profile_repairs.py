"""Check repaired intervals against saved evaluated curves. Never certify from COMPLETE alone."""
import csv, math, sys
from pathlib import Path
root=Path(sys.argv[1] if len(sys.argv)>1 else 'results/training60_final/profile_repairs')
rows=[]
for table in sorted(root.glob('*/*/profile_confidence_intervals.csv')):
    for v in csv.DictReader(table.open()):
        points=list(csv.DictReader((table.parent/v['parameter']/'profile.csv').open()))
        mle=float(v['mle']);issues=[]
        for side in ['lower','upper']:
            branch=[p for p in points if (float(p['parameter_value'])<mle if side=='lower' else float(p['parameter_value'])>mle)]
            crosses=any(float(p['profile_loglik'])<=-1.920729410347062 for p in branch)
            bound=any((float(p['normalized_parameter'])<=1e-7 if side=='lower' else float(p['normalized_parameter'])>=1-1e-7) for p in branch)
            if not crosses and not bound: issues.append(side+' incomplete')
            if crosses and math.isclose(float(v['package_'+side]),mle,rel_tol=1e-8): issues.append(side+' invalid interval endpoint')
            if not crosses and bound: issues.append(side+' bound-limited')
        inside=[float(p['profile_loglik'])>=-1.920729410347062 for p in points]
        components=sum(value and (j==0 or not inside[j-1]) for j,value in enumerate(inside))
        if components>1: issues.append(f'{components} separate confidence-region segments: inspect curve')
        dips=0
        for j in range(1,len(points)-1):
            a,b,c=[float(points[k]['profile_loglik']) for k in (j-1,j,j+1)]
            if min(a,c)-b>0.2 and max(a,c)>-3: dips+=1
        if dips: issues.append(f'{dips} sharp dips near confidence region: inspect curve')
        if max(float(p['profile_loglik']) for p in points)>1e-5: issues.append('better likelihood than reference: investigate MLE')
        if any(p['stopping_met']!='true' for p in points): issues.append('nuisance optimization requires review')
        rows.append(dict(mode=table.parent.parent.name,scenario=table.parent.name,parameter=v['parameter'],status=' | '.join(issues) or 'thresholds crossed on both sides'))
if not rows: raise SystemExit('No completed repair tables found')
with (root/'audit.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(f'Audited {len(rows)} repaired profiles. Review {root / "audit.csv"} before publication.')
