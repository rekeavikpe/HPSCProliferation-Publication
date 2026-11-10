"""Verify saved split and all-data results and assemble publication reporting tables."""
import argparse
import csv
import hashlib
import html
import math
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

MAP = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 5}


def read(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f))


def write(p, rows):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('split', type=Path)
    parser.add_argument('results_root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    split, root, out = (p.resolve() for p in (args.split,args.results_root,args.output))
    fit = root/'calibration100_hybrid_lac100_kg100'
    profile = root/'active_fixed_three_hybrid_lac100_kg100'
    media = root/'media_exchange_hybrid_lac100_kg100'
    frozen = root/'cell_lines_frozen_hybrid_lac100_kg100'
    adapted = root/'cell_lines_recalibrated_hybrid_lac100_kg100'
    for source in (split,fit,profile,media,frozen,adapted):
        if out == source or out in source.parents or source in out.parents:
            raise ValueError('Assembly output must be separate from source result folders')
    if not (split/'COMPLETE').is_file():
        raise ValueError('Split results incomplete')
    ranks=read(split/'selection_ranking.csv'); lock=read(split/'selection_lock.csv')[0]
    if {int(r['publication_model']) for r in ranks} != set(range(7)):
        raise ValueError('Expected seven models')
    best=min(float(r['selection_nll']) for r in ranks)
    selected=min((r for r in ranks if float(r['selection_nll'])<=best+1e-6),key=lambda r:(int(r['k']),float(r['training_bic']),int(r['publication_model'])))
    if selected['storage_model'] != lock['storage_model'] or int(lock['storage_model'])!=5:
        raise ValueError('Selection differs from Model 6; supply matching refit and evaluations')
    splitp=split/f"model{lock['storage_model']}/selected_parameters.csv"
    if sha(splitp)!=lock['parameter_sha256']:
        raise ValueError('Split selection lock parameter hash differs')
    assignments=read(split/'split.csv')
    if len({r['series_id'] for r in assignments})!=len(assignments):
        raise ValueError('Culture appears more than once in split')
    if sum(r['role']=='calibration' for r in assignments)!=19 or sum(r['role']=='selection' for r in assignments)!=12:
        raise ValueError('Unexpected partition counts')
    points=read(split/'predictions.csv'); metrics=[]
    for rank in ranks:
        for cohort,key in (('calibration','training_nll'),('selection','selection_nll')):
            group=[r for r in points if r['publication_model']==rank['publication_model'] and r['cohort']==cohort]
            if len(group)!=int(rank['training_n' if cohort=='calibration' else 'selection_n']):
                raise ValueError('Missing split prediction points')
            score=sum(math.log(float(r['sigma']))+math.log(2*math.pi)/2+((float(r['predicted'])-float(r['observed']))/float(r['sigma']))**2/2 for r in group)
            if abs(score-float(rank[key]))>1e-6:
                raise ValueError(f'Split likelihood mismatch: Model {rank["publication_model"]}, {cohort}')
            if cohort=='calibration' and abs(2*score+int(rank['k'])*math.log(len(group))-float(rank['training_bic']))>2e-6:
                raise ValueError('Training BIC mismatch')
            for state in ('cells','glucose','lactate','aggregate','all'):
                sub=[r for r in group if (state=='all' or r['state']==state) and float(r['observed'])!=0]
                noise=[r for r in sub if int(r['n_vessels'])>=2]
                metrics.append(dict(publication_model=rank['publication_model'],cohort=cohort,state=state,
                    mrpe_percent=statistics.mean(100*abs(float(r['predicted'])-float(r['observed']))/abs(float(r['observed'])) for r in sub),
                    experimental_noise_percent=statistics.mean(100*float(r['sd'])/abs(float(r['observed'])) for r in noise),n=len(sub)))
    # Verify every evaluation references the same saved all-data coefficients.
    for source,ids in ((profile,(5,)),(media,(0,5)),(frozen,(0,5)),(adapted,(0,5))):
        hashes=read(source/'source_hashes.csv')
        for i in ids:
            suffix=f'model{i}/selected_parameters.csv'
            matches=[r for r in hashes if r['path'].endswith(suffix)]
            if len(matches)!=1 or matches[0]['sha256']!=sha(fit/suffix):
                raise ValueError(f'Parameter provenance mismatch: {source.name}, {suffix}')
    p={r['parameter']:float(r['value']) for r in read(fit/'model5/selected_parameters.csv')}
    fixed=read(profile/'fixed_parameters.csv')
    if {r['parameter'] for r in fixed}!={'K_G','K_Lac','Y_GX'}:
        raise ValueError('Unexpected conditional parameter policy')
    if not all(math.isclose(float(r['value']),p[r['parameter']],rel_tol=1e-12) for r in fixed):
        raise ValueError('Conditional fixed values differ from all-data MLE')
    intervals=read(profile/'profile_confidence_intervals.csv')
    if len(intervals)!=8 or any(r['lower_hits_bound']!='false' or r['upper_hits_bound']!='false' or r['better_minimum']!='false' or int(r['unconfirmed_points'])!=0 for r in intervals):
        raise ValueError('Conditional profile qualifications changed; inspect source suite')
    out.mkdir(parents=True,exist_ok=True)
    for stale in ('COMPLETE','TABLES_VERIFIED'):
        (out/stale).unlink(missing_ok=True)
    write(out/'tables/split_selection.csv',ranks)
    write(out/'tables/split_metrics.csv',metrics)
    write(out/'tables/fixed_parameters.csv',fixed)
    write(out/'tables/conditional_intervals.csv',read(profile/'profile_confidence_intervals.csv'))
    params=[]
    for rank in read(fit/'ranking.csv'):
        for row in read(fit/f"model{rank['model']}/selected_parameters.csv"):
            params.append(dict(publication_model=MAP[int(rank['model'])],parameter=row['parameter'],value=row['value']))
    write(out/'tables/all100_parameters.csv',params)
    script=Path(__file__).with_name('report_prediction_metrics.py')
    subprocess.run([sys.executable,str(script),str(fit),'--media',str(media),'--cell-lines',str(adapted),'--out',str(out/'metrics_source')],check=True)
    report=[]
    for r in read(out/'metrics_source/prediction_metrics.csv'):
        if r['scope']!='cohort': continue
        report.append(dict(publication_model=MAP[int(r['model'])],evaluation=r['cohort'].replace('frozen transfer','parameter transfer'),
            state=r['state'],mrpe_percent=r['mrpe_percent'],experimental_noise_percent=r['experimental_noise_cv_percent'],n=r['n_points']))
    write(out/'tables/evaluation_metrics.csv',report)
    # Copy evaluation inputs; plotting cannot change their source counterparts.
    for source,name in ((media,'media'),(frozen,'cell_transfer'),(adapted,'cell_recalibration')):
        dest=out/name;dest.mkdir(exist_ok=True)
        for f in source.iterdir():
            if f.is_file(): shutil.copy2(f,dest/f.name)
        # Verify required prediction inputs rather than requiring markers omitted by the export archive.
        if not (dest/'predictions.csv').exists() or not (dest/'trajectories.csv').exists():
            raise ValueError('Missing evaluation plotting inputs')
        (dest/'COMPLETE').write_text('Existing prediction tables and source parameter hashes verified; no fitting performed.\n')
    files=[f for d in (split,fit,profile,media,frozen,adapted) for f in d.rglob('*') if f.is_file() and f.suffix in ('.csv','.toml')]
    write(out/'source_manifest.csv',[dict(path=str(f),sha256=sha(f)) for f in files])
    def table(headers,rows):
        return '<table><tr>'+''.join('<th>'+html.escape(h)+'</th>' for h in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in row)+'</tr>' for row in rows)+'</table>'
    summary=table(['Evaluation','Overall MRPE (%)','Experimental Noise (%)'],[(r['evaluation'],f"{float(r['mrpe_percent']):.2f}",f"{float(r['experimental_noise_percent']):.2f}") for r in report if r['publication_model']==6 and r['state']=='all'])
    (out/'index.html').write_text('''<!doctype html><meta charset="utf-8"><title>Final results after split selection</title><style>body{font:16px Arial;max-width:1150px;margin:40px auto;line-height:1.5}td,th{padding:10px;border-bottom:1px solid #ccc;text-align:left}img{max-width:100%}</style><h1>Final results after 60/40 selection</h1><p>Left comparison: models fitted to 19 cultures, selected on 12 other cultures. Profiles and external predictions: matching all-100-mL fit after selection. Refined split likelihoods, BIC and MRPEs have been checked from their saved predictions. Split nonconvergence flags remain in the tables.</p><h2>Model selection and scale transfer</h2><img src="combined/bic_and_500mL.svg"><p><a href="combined_cfd/cfd_bic_and_500mL.pdf">Separate combined CFD, BIC and 500-mL figure</a></p><h2>Final Model 6 errors</h2>'''+summary+'''<p>Experimental Noise is descriptive replicate SD/mean, not an acceptance threshold. Cell-line recalibration is in-sample. Conditional profiles hold three data-derived coefficients fixed and do not propagate their uncertainty.</p><ul><li><a href="tables/split_selection.csv">Final refined selection and convergence</a></li><li><a href="tables/split_metrics.csv">Final refined split MRPEs and noise</a></li><li><a href="tables/evaluation_metrics.csv">All evaluation metrics</a></li><li><a href="tables/all100_parameters.csv">All-data parameters, publication numbering</a></li><li><a href="tables/conditional_intervals.csv">Conditional intervals</a></li><li><a href="calibration/training_with_hydrodynamics.pdf">Five-condition calibration</a></li><li><a href="profiles/model6_profiles.pdf">Native profile plots</a></li><li><a href="media/index.html">Media exchange</a></li><li><a href="cell_recalibration/index.html">Both cell-line comparisons</a></li><li><a href="source_manifest.csv">Source file hashes</a></li></ul><p>Open profiles and data-derived fixing remain documented in the source suites. The 585A1 no-dextran-sulfate source cohort must be reconciled with manuscript Methods.</p>''')
    (out/'TABLES_VERIFIED').write_text('Likelihood, BIC, split membership, MRPE and parameter provenance verified. Figures are produced by the wrapper.\n')
    print('Verified final assembly tables:',out)

if __name__=='__main__':main()
