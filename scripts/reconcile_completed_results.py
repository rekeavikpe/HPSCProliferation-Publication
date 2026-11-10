"""Audit an exported server bundle without refitting or changing saved results."""
import argparse
import csv
import hashlib
import math
from pathlib import Path


def read(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bundle', type=Path)
    parser.add_argument('out', type=Path)
    args = parser.parse_args()
    bundle, out = args.bundle.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    fit = bundle / 'calibration100_hybrid_lac100_kg100'
    profiles = bundle / 'active_fixed_three_hybrid_lac100_kg100'
    checks = []
    sources = {}

    def load(path):
        sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return read(path)

    def check(name, passed, detail):
        checks.append(dict(check=name, passed=bool(passed), detail=detail))

    mapping = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 5}
    for rank in load(fit / 'ranking.csv'):
        rows = load(fit / f'model{rank["model"]}' / 'residuals.csv')
        nll = sum(math.log(float(r['sigma'])) + math.log(2*math.pi)/2
                  + ((float(r['predicted'])-float(r['observed']))/float(r['sigma']))**2/2
                  for r in rows)
        bic = 2*nll + int(rank['k'])*math.log(len(rows))
        check(f'Model {mapping[int(rank["model"])]} saved likelihood',
              abs(nll-float(rank['nll'])) < 1e-6 and abs(bic-float(rank['bic'])) < 2e-6,
              f'n={len(rows)}; NLL={nll:.12g}; BIC={bic:.12g}')

    p = {r['parameter']: float(r['value']) for r in load(fit/'model5/selected_parameters.csv')}
    for row in load(profiles/'fixed_parameters.csv'):
        check('Fixed '+row['parameter'], math.isclose(float(row['value']), p[row['parameter']], rel_tol=1e-12), row['value'])
    intervals = load(profiles/'profile_confidence_intervals.csv')
    for r in intervals:
        check('Conditional profile '+r['parameter'],
              float(r['package_lower']) < float(r['mle']) < float(r['package_upper'])
              and all(r[k]=='false' for k in ('lower_hits_bound','upper_hits_bound','better_minimum'))
              and int(r['unconfirmed_points'])==0,
              'Saved package interval and optimization flags; conditional on the three fixed values')
    for folder in ('media_exchange_hybrid_lac100_kg100', 'cell_lines_frozen_hybrid_lac100_kg100',
                   'cell_lines_recalibrated_hybrid_lac100_kg100'):
        hashes = load(bundle/folder/'source_hashes.csv')
        for model in (0, 5):
            suffix = f'model{model}/selected_parameters.csv'
            matches = [r for r in hashes if r['path'].endswith(suffix)]
            expected = hashlib.sha256((fit/suffix).read_bytes()).hexdigest()
            check(folder+' source '+suffix, len(matches)==1 and matches[0]['sha256']==expected,
                  'Exact parameter-file SHA256 match')

    # This deterministic contrast uses the final all-100-mL parameters in both vessels.
    large = load(fit/'scale_transfer/test_predictions.csv')
    contrasts = []
    for model in (0, 5):
        small = load(fit/f'model{model}/residuals.csv')
        for r in small:
            if r['condition']!='5' or r['state']!='aggregate':
                continue
            s = next(x for x in large if x['publication_model']==str(model)
                     and x['day']==r['day'] and x['state']=='aggregate')
            contrasts.append(dict(publication_model=mapping[model], day=r['day'],
                observed_100_minus_500_um=float(r['observed'])-float(s['observed']),
                predicted_100_minus_500_um=float(r['predicted'])-float(s['predicted']),
                source_100mL_condition=5, publication_100mL_condition=4,
                source_500mL_condition=4))
    write(out/'final_vessel_contrast.csv', contrasts)
    write(out/'verification.csv', checks)
    write(out/'conditional_intervals.csv', intervals)
    write(out/'source_manifest.csv', [dict(path=p, sha256=h) for p,h in sources.items()])
    failures = [r for r in checks if not r['passed']]
    print(f'{len(checks)-len(failures)}/{len(checks)} checks passed; outputs: {out}')
    if failures:
        raise RuntimeError(failures)


if __name__ == '__main__':
    main()
