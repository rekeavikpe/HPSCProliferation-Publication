import importlib.util,unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location('r',Path(__file__).resolve().parents[1]/'scripts/run_cell_reductions.py');r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
class Designs(unittest.TestCase):
 def test_exhaustive(self):
  for bad,n in [(['K_d','Y_GLac','alpha'],7),(['K_d','Y_GLac'],3)]:
   audit=[dict(parameter=p,status='not_bounded' if p in bad else 'bounded') for p in r.PARAMETERS]
   plans=r.candidates(audit);self.assertEqual(len(plans),n)
   self.assertEqual(len({tuple(x) for x,y in plans}),n)
   for fixed,free in plans:
    self.assertFalse(set(fixed)&set(free));self.assertEqual(set(fixed)|set(free),set(r.PARAMETERS));self.assertTrue(set(fixed)<=set(bad))
 def test_numerical_blocks(self):
  with self.assertRaises(RuntimeError):r.candidates([dict(parameter=p,status='numerical_review') for p in r.PARAMETERS])
if __name__=='__main__':unittest.main()
