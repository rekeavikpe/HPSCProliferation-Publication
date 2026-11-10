import importlib.util
from pathlib import Path
import unittest
path=Path(__file__).resolve().parents[1]/'scripts/report_prediction_metrics.py'
spec=importlib.util.spec_from_file_location('metrics',path)
metrics=importlib.util.module_from_spec(spec);spec.loader.exec_module(metrics)
class MetricsTests(unittest.TestCase):
    def test_relative_errors_and_actual_sd(self):
        rows=[dict(observed=10.,predicted=12.,sd=1.,n=3),dict(observed=20.,predicted=16.,sd=4.,n=3)]
        r=metrics.summarize(rows,'x','5','glucose','cohort')
        self.assertEqual(r['mrpe_percent'],20.)
        self.assertEqual(r['experimental_noise_cv_percent'],15.)
        self.assertEqual(r['bias'],-1.)
        self.assertAlmostEqual(r['rmse'],10**.5)
    def test_zero_means_and_singletons_not_invented(self):
        rows=[dict(observed=0.,predicted=2.,sd=0.,n=3),dict(observed=10.,predicted=12.,sd=0.,n=1)]
        r=metrics.summarize(rows,'x','5','glucose','cohort')
        self.assertEqual((r['n_points'],r['n_mrpe'],r['n_noise']),(2,1,0))
        self.assertEqual(r['experimental_noise_cv_percent'],'')
        self.assertEqual(r['mrpe_percent'],20.)
    def test_pooled_states_do_not_mix_physical_rmse_units(self):
        r=metrics.summarize([dict(observed=10.,predicted=12.,sd=1.,n=3)],'x','5','all','cohort')
        self.assertEqual((r['rmse'],r['bias']),('',''))
if __name__=='__main__': unittest.main()
