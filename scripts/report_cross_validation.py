#!/usr/bin/env python3
"""Tabulate completed cross-validation results without fitting or solving models."""
import argparse
import csv
import hashlib
import html
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

STATES = ('cells', 'glucose', 'lactate', 'aggregate', 'all')
LABELS = ('Cells', 'Glucose', 'Lactate', 'Aggregate size', 'Overall')

def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = args.results.resolve()
    out = args.output or root / 'tables'
    metrics = read(root / 'state_metrics.csv')
    fits = {(r['fold'], r['publication_model']): r for r in read(root / 'fit_summary.csv')}
    folds = {r['fold']: r for r in read(root / 'folds.csv')}
    expected = {'reference100', *(f'loco_condition{i}' for i in range(1, 6)),
                *(f'partial500_{i:03d}' for i in range(1, 10))}
    if set(folds) != expected:
        raise ValueError('Expected the completed balanced-nine design: reference, five LOCO folds, nine partial-500 splits.')
    groups = defaultdict(dict)
    for r in metrics:
        key = tuple(r[k] for k in ('fold', 'experiment', 'variant', 'publication_model', 'cohort'))
        if r['state'] in groups[key]:
            raise ValueError(f'Duplicate metric: {key}, {r["state"]}')
        groups[key][r['state']] = r
    expected_groups = set()
    for fold in expected:
        experiment = folds[fold]['experiment']
        variants = ('partial500_fit', '100mL_only_reference') if experiment == 'partial500' else ('fold_fit',)
        for model in ('0', '6'):
            for variant in variants:
                for cohort in ('calibration', 'prediction'):
                    if variant == '100mL_only_reference' and cohort == 'calibration':
                        continue
                    expected_groups.add((fold, experiment, variant, model, cohort))
    if set(groups) != expected_groups:
        raise ValueError('Incomplete or unexpected model/fold/evaluation combinations.')
    # Independently reproduce the saved MRPEs from prediction rows.
    raw = defaultdict(list)
    for r in read(root / 'predictions.csv'):
        if float(r['day']) <= 0 or float(r['mean']) == 0:
            continue
        key = tuple(r[k] for k in ('fold', 'experiment', 'variant', 'publication_model', 'cohort'))
        error = 100 * abs((float(r['predicted']) - float(r['mean'])) / float(r['mean']))
        raw[key + (r['state'],)].append(error)
        raw[key + ('all',)].append(error)
    rows = []
    for key, states in groups.items():
        if set(states) != set(STATES):
            raise ValueError(f'Missing states: {key}')
        fold, experiment, variant, model, cohort = key
        source = 'reference100' if variant == '100mL_only_reference' else fold
        fit = fits[source, model]
        row = dict(zip(('fold', 'experiment', 'variant', 'publication_model', 'cohort'), key))
        row.update(training_cultures=folds[source]['training_cultures'],
                   evaluated_cultures=folds[fold]['training_cultures' if cohort == 'calibration' else 'prediction_cultures'],
                   stopping_met=fit['training_stopping_met'], training_bic=fit['training_bic'])
        for state in STATES:
            r = states[state]
            values = raw[key + (state,)]
            value = float(r['mrpe_percent'])
            if len(values) != int(r['n']) or not math.isclose(mean(values), value, rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f'MRPE mismatch: {key}, {state}')
            row[state + '_mrpe_percent'] = value
            row[state + '_experimental_noise_percent'] = float(r['experimental_noise_percent'])
            row[state + '_n'] = len(values)
        rows.append(row)
    rows.sort(key=lambda r: (r['experiment'] != 'reference', r['experiment'], r['fold'], -int(r['publication_model']), r['variant'], r['cohort']))
    summaries = []
    by_design = defaultdict(list)
    for row in rows:
        if row['experiment'] != 'reference':
            by_design[tuple(row[k] for k in ('experiment', 'variant', 'publication_model', 'cohort'))].append(row)
    for key, records in by_design.items():
        row = dict(zip(('experiment', 'variant', 'publication_model', 'cohort'), key))
        row['fold_count'] = len(records)
        for state in STATES:
            row[state + '_mean_fold_mrpe_percent'] = mean(r[state + '_mrpe_percent'] for r in records)
            row[state + '_pooled_point_mrpe_percent'] = sum(r[state + '_mrpe_percent'] * r[state + '_n'] for r in records) / sum(r[state + '_n'] for r in records)
        summaries.append(row)
    out.mkdir(parents=True, exist_ok=True)
    write(out / 'errors_by_run.csv', rows)
    write(out / 'average_errors_by_design.csv', summaries)
    metadata = {'source': str(root), 'verified_metric_rows': len(metrics), 'definition': 'MRPE = mean absolute relative prediction error, times 100; day > 0 and nonzero observed condition-day mean. Overall pools eligible state-day observations.',
                'averaging': 'Mean-fold gives equal weight to each fold. Pooled-point weights folds by eligible observation count. Overlapping partial-500 splits are not independent experiments.',
                'files': {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in ('state_metrics.csv','fit_summary.csv','folds.csv','predictions.csv')}}
    (out / 'provenance.json').write_text(json.dumps(metadata, indent=2))
    parts = ['<!doctype html><meta charset="utf-8"><title>Cross-validation error tables</title><style>body{font:16px Arial;margin:36px;color:#17212b}table{border-collapse:collapse;margin:20px 0;width:100%}th,td{padding:9px;border-bottom:1px solid #ddd;text-align:right}th:first-child,td:first-child{text-align:left}th{background:#eaf0f5}h2{margin-top:32px}p{max-width:1100px;line-height:1.5}</style><h1>Cross-validation: calibration and prediction errors</h1>',
             '<p>All entries are MRPE (%), shown to one decimal place. Each row represents a selected fit and evaluation cohort, not an individual optimizer start. Calibration errors and held-out prediction errors are kept separate. LOCO predicts an omitted 100 mL condition; the reference predicts 500 mL, so these assess different prediction tasks.</p>']
    for model in ('6', '0'):
        for experiment, title in (('reference','100 mL calibration → 500 mL prediction'),('loco','Leave one 100 mL condition out'),('partial500','Four 500 mL cultures added to calibration; five held out')):
            parts.append(f'<h2>Model {model}: {title}</h2><table><tr><th>Run / approach</th><th>Evaluation</th>'+''.join(f'<th>{s}</th>' for s in LABELS)+'<th>Stopping met</th></tr>')
            for r in rows:
                if r['publication_model'] != model or r['experiment'] != experiment: continue
                label=r['fold'] + (' — 100 mL-only reference' if r['variant']=='100mL_only_reference' else '')
                parts.append('<tr><td>'+html.escape(label)+'</td><td>'+r['cohort']+'</td>'+''.join(f'<td>{r[s+"_mrpe_percent"]:.1f}</td>' for s in STATES)+f'<td>{r["stopping_met"]}</td></tr>')
            for r in summaries:
                if r['publication_model'] != model or r['experiment'] != experiment: continue
                parts.append('<tr><th>Mean across folds — '+html.escape(r['variant'])+'</th><th>'+r['cohort']+'</th>'+''.join(f'<th>{r[s+"_mean_fold_mrpe_percent"]:.1f}</th>' for s in STATES)+'<td>See individual fits</td></tr>')
            parts.append('</table>')
    parts.append('<p>'+html.escape(metadata['definition'])+' '+html.escape(metadata['averaging'])+' Experimental noise and observation counts are retained in errors_by_run.csv. A failed stopping flag is not concealed by averaging.</p>')
    (out / 'index.html').write_text('\n'.join(parts))
    print(f'Verified {len(metrics)} saved MRPEs against predictions. Tables: {out}')

if __name__ == '__main__':
    main()
