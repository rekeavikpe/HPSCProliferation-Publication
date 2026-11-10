"""Regression checks for no-drop reporting and narrowly scoped resume migration."""
import csv,importlib.util,json,pathlib,runpy,sys,tempfile,unittest
from unittest.mock import patch
ROOT=pathlib.Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('resume',ROOT/'scripts/october20_resume.py');resume=importlib.util.module_from_spec(spec);spec.loader.exec_module(resume)
class NoDrop(unittest.TestCase):
    def test_unbounded_and_numerical_keep_five(self):
        with tempfile.TemporaryDirectory() as t:
            out=pathlib.Path(t)/'cells';calls=[]
            def run(cmd,**kw):
                calls.append(cmd)
                if 'profiles/run_october_cells.jl' in cmd:
                    p=pathlib.Path(cmd[-1]);p.mkdir(parents=True)
                    with (p/'audit.csv').open('w') as f:
                        w=csv.writer(f);w.writerow(['parameter','status'])
                        for n in ['mu_max','K_d','K_oG','Y_GLac','alpha']:
                            w.writerow([n,('not_bounded' if '585A1' in str(p) else 'numerical_review') if n=='alpha' else 'bounded'])
            with patch.object(sys,'argv',['driver','fit',str(out)]),patch('subprocess.run',side_effect=run):
                runpy.run_path(str(ROOT/'scripts/run_october_cells.py'),run_name='__main__')
            rows=json.loads((out/'status.json').read_text())
            self.assertEqual(len(calls),6)
            self.assertEqual([r['status'] for r in rows],['parameter_reduction_review_required','profile_numerical_review_required'])
            self.assertTrue(all(len(r['active_parameters'])==5 and not r['automatic_reduction'] for r in rows))
            self.assertFalse(list(out.glob('*/round2')))
    def migration(self,unrelated=False):
        with tempfile.TemporaryDirectory() as t:
            root=pathlib.Path(t);out=root/'results/suite';out.mkdir(parents=True)
            sources={str(root/p):h for p,h in resume.OLD.items()}
            sources.update({str(root/'fit/selection_optimizer_bic.csv'):'fit',str(root/'split/model1/summary.csv'):'split',str(root/'src/a.jl'):'unchanged'})
            old={'sources':sources,'settings':{}}
            current=json.loads(json.dumps(old))
            for p in resume.OLD: current['sources'][str(root/p)]='corrected'
            if unrelated: current['sources'][str(root/'src/a.jl')]='changed'
            (out/'suite_manifest.json').write_text(json.dumps(old))
            for name in ['p75.DONE','day1.DONE','cells.DONE','figures.DONE','report.DONE','COMPUTATIONS_COMPLETE','cells/585A1/round1/checkpoint','cells/585A1/round2/checkpoint']:
                p=out/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('saved')
            def manifest(cmd,**kw): (pathlib.Path(cmd[-1])/'suite_manifest.json').write_text(json.dumps(current))
            with patch.object(resume,'ROOT',root),patch.object(resume,'targets',return_value=set()),patch.object(resume.subprocess,'run',side_effect=manifest):
                if unrelated:
                    with self.assertRaises(SystemExit): resume.migrate(out)
                    self.assertTrue((out/'cells.DONE').exists())
                else:
                    resume.migrate(out)
                    self.assertTrue((out/'p75.DONE').exists());self.assertTrue((out/'day1.DONE').exists())
                    self.assertTrue((out/'cells/585A1/round1/checkpoint').exists())
                    self.assertFalse((out/'cells/585A1/round2').exists());self.assertFalse((out/'cells.DONE').exists())
                    self.assertEqual(json.loads((out/'suite_manifest.json').read_text()),current)
                    resume.migrate(out) # idempotent
    def test_migration_preserves_fits(self): self.migration()
    def test_migration_rejects_other_changes(self): self.migration(True)
    def test_process_scope(self):
        out=ROOT/'results/david_comments_20261020'
        rows={101:(1,'S','bash scripts/run_october20_revision.sh'),102:(101,'R','julia worker.jl'),103:(1,'R','julia unrelated.jl'),104:(1,'S','bash scripts/run_october20_revision.sh results/another_suite')}
        with patch.object(resume,'processes',return_value=rows): self.assertEqual(resume.targets(out),{101,102})
if __name__=='__main__': unittest.main()
