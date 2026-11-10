"""Report prediction error and replicate variability on the saved evaluation points."""
import argparse
import csv
import html
import math
import statistics as stats
from pathlib import Path

STATES = ('cells', 'glucose', 'lactate', 'aggregate')

def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def number(value):
    return float(value) if value not in ('', None) else math.nan

def summarize(rows, cohort, model, state, scope):
    eligible = [r for r in rows if math.isfinite(r['observed']) and math.isfinite(r['predicted'])]
    relative = [r for r in eligible if r['observed'] != 0]
    variability = [r for r in relative if r['n'] >= 2 and math.isfinite(r['sd']) and r['sd'] >= 0]
    average = lambda xs: stats.mean(xs) if xs else ''
    return dict(cohort=cohort, model=model, scope=scope, state=state,
                n_points=len(eligible), n_mrpe=len(relative), n_noise=len(variability),
                mrpe_percent=average([100*abs(r['predicted']-r['observed'])/abs(r['observed']) for r in relative]),
                experimental_noise_cv_percent=average([100*r['sd']/abs(r['observed']) for r in variability]),
                mrpe_on_noise_available_points_percent=average([100*abs(r['predicted']-r['observed'])/abs(r['observed']) for r in variability]),
                rmse=math.sqrt(stats.mean([(r['predicted']-r['observed'])**2 for r in eligible])) if eligible and state!='all' else '',
                bias=average([r['predicted']-r['observed'] for r in eligible]) if state!='all' else '')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('calibration',type=Path)
    parser.add_argument('--media',type=Path)
    parser.add_argument('--cell-lines',type=Path,help='Directory containing predictions.csv with cohort, model, day, state, observed, predicted, sd, n_vessels; frozen and recalibrated cohorts must be distinct')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    groups={}; sources=[]; pending=[]
    def verify_evaluation(directory):
        import hashlib
        hashes=read(directory/'source_hashes.csv')
        for model in ('0','5'):
            suffix=f'model{model}/selected_parameters.csv'
            expected=hashlib.sha256((args.calibration/suffix).read_bytes()).hexdigest()
            matching=[r for r in hashes if r['path'].endswith(suffix)]
            if len(matching)!=1 or matching[0]['sha256']!=expected:
                raise ValueError(f'{directory} does not use the supplied calibration parameters for model {model}')

    def add(cohort, model, row, scope='cohort'):
        if number(row.get('day',1))<=0:
            return
        entry=dict(state=row['state'],observed=number(row['observed']),predicted=number(row['predicted']),sd=number(row['sd']),n=number(row['n_vessels']))
        groups.setdefault((cohort,str(model),scope),[]).append(entry)
    calibration_label='100 mL training (60%)' if (args.calibration/'source_manifest.toml').exists() else '100mL calibration'
    obs_path=args.calibration/'observations.csv'
    obs={(r['condition'],r['day'],r['state']):r for r in read(obs_path)}
    sources.append(obs_path)
    for rank in read(args.calibration/'ranking.csv'):
        model=rank['model']; path=args.calibration/f'model{model}'/'residuals.csv'
        if not path.exists():
            raise FileNotFoundError(path)
        sources.append(path)
        for row in read(path):
            key=(row['condition'],row['day'],row['state']); observation=obs[key]
            assert math.isclose(number(row['observed']),number(observation['mean']),rel_tol=1e-10,abs_tol=1e-10),key
            row.update(sd=observation['sd'],n_vessels=observation['n_vessels'])
            add(calibration_label,model,row)
            add(calibration_label,model,row,'original condition '+row['condition'])
    scale=args.calibration/'scale_transfer/test_predictions.csv'
    if scale.exists():
        sources.append(scale)
        for row in read(scale): add('500 mL prediction',row['publication_model'],row)
    else: pending.append('500mL prediction results not supplied')
    if args.media:
        verify_evaluation(args.media)
        path=args.media/'predictions.csv'; sources.append(path)
        for row in read(path): add(row['protocol_percent']+'% media exchange',row['model'],row)
    else: pending.append('20% and 50% media exchange results not supplied')
    if args.cell_lines:
        verify_evaluation(args.cell_lines)
        path=args.cell_lines/'predictions.csv'; sources.append(path)
        for row in read(path): add(row['cohort'],row['model'],row)
    else: pending.append('585A1 and SCTi003-A: evaluations pending; no values inferred or copied from the manuscript')
    results=[]
    for (cohort,model,scope),rows in groups.items():
        for state in STATES+('all',):
            subset=rows if state=='all' else [r for r in rows if r['state']==state]
            if subset: results.append(summarize(subset,cohort,model,state,scope))
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/'prediction_metrics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=results[0]);writer.writeheader();writer.writerows(results)
    import hashlib
    with (args.out/'metric_sources.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['path','sha256'])
        for p in sources+[Path(__file__)]: writer.writerow([str(p.resolve()),hashlib.sha256(p.read_bytes()).hexdigest()])
    methods=('MRPE (%) = 100 × mean(|prediction − observed mean| / |observed mean|). '
             'Experimental noise (%) = 100 × mean(sample SD / |observed mean|), using at least two measured replicates. '
             'Day-zero values are excluded. Observed means equal to zero are excluded from relative metrics; valid zero replicate measurements remain in the source means and SDs. '
             'No SD floors, SEM, assumed noise, or model residuals replace experimental SD. '
             'The all-state row weights each observed condition–day–state mean equally, following the original pooled-point summary. '
             'RMSE and bias are only calculated within a state, in its physical units. '
             'Counts identify missing variability estimates; the matched-point MRPE column allows comparison on the same subset. '
             'CV is descriptive between-replicate variability, not a validation acceptance threshold or a prediction interval. '
             'These metrics do not account for repeated-day correlation.')
    (args.out/'metric_definitions.txt').write_text(methods+'\n')
    selected=read(args.calibration/'selection_lock.csv')[0]['model']
    visible=[r for r in results if r['scope']=='cohort' and r['model']==selected]
    def fmt(v): return f'{v:.2f}' if isinstance(v,(float,int)) else html.escape(str(v))
    table=''.join('<tr>'+''.join('<td>'+fmt(r[k])+'</td>' for k in ('cohort','state','mrpe_percent','experimental_noise_cv_percent','n_mrpe','n_noise'))+'</tr>' for r in visible)
    if args.cell_lines:
        pending.append('585A1 cultures received dextran sulfate experimentally. The model does not include an explicit dextran-sulfate term. Recalibration rows are fits to the respective cell-line data.')
    notice=''.join('<li>'+html.escape(p)+'</li>' for p in pending)
    (args.out/'index.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>Prediction error and experimental variability</title><style>body{{font:16px system-ui;max-width:1150px;margin:36px auto;line-height:1.5}}table{{border-collapse:collapse;width:100%}}td,th{{padding:8px;border-bottom:1px solid #ddd;text-align:left}}tr:nth-child(even){{background:#f3f6f8}}</style><h1>Model { {'5':'6','6':'5'}.get(selected,selected) }: prediction error and experimental variability</h1><p>Calibration source: {html.escape(args.calibration.name)}</p><p>{html.escape(methods)}</p><table><tr><th>Evaluation</th><th>State</th><th>MRPE (%)</th><th>Experimental Noise (%)</th><th>MRPE points</th><th>Noise points</th></tr>{table}</table><ul>{notice}</ul><p><a href="prediction_metrics.csv">All models and calibration conditions</a> · <a href="metric_sources.csv">Input provenance</a></p>''')
    print(f'Saved {len(results)} summaries to {args.out}; '+ '; '.join(pending))

if __name__=='__main__': main()
