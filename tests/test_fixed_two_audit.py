import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location("fixed_two_audit", Path(__file__).resolve().parents[1] / "scripts/audit_fixed_two_profiles.py")
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)

def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def fixture(root):
    write_csv(root / "fixed_parameters.csv", [dict(parameter="K_Lac", value=100), dict(parameter="Y_GX", value=3.96e-6)])
    write_csv(root / "mle_summary.csv", [dict(free_parameters=9, stopping_met="true", smoke_test="false", saved_nll=10, nll=10)])
    intervals=[]
    for name in sorted(a.EXPECTED):
        intervals.append(dict(parameter=name, mle=1, package_lower=.55, package_upper=1.45))
        points=[dict(parameter_value=x, normalized_parameter=z, profile_loglik=y, stopping_met="true") for x,z,y in ((.5,0,-3),(1,.5,0),(1.5,1,-3))]
        write_csv(root/name/"profile.csv",points)
        write_csv(root/name/"nuisance_parameters.csv",[dict(point=1,parameter=n,value=1) for n in sorted(a.EXPECTED)])
        (root/name/"COMPLETE").write_text("complete")
    write_csv(root/"profile_confidence_intervals.csv",intervals)

class AuditTests(unittest.TestCase):
    def test_closed_profiles_and_bound_limited_status(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); fixture(root)
            self.assertTrue(all(r["status"]=="thresholds crossed on both sides" for r in a.audit(root)))
            rows=a.read_csv(root/"K_G/profile.csv"); rows[-1]["profile_loglik"]=-1
            write_csv(root/"K_G/profile.csv",rows)
            status=next(r["status"] for r in a.audit(root) if r["parameter"]=="K_G")
            self.assertIn("upper bound-limited",status)
    def test_missing_kg_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); fixture(root)
            rows=[r for r in a.read_csv(root/"profile_confidence_intervals.csv") if r["parameter"]!="K_G"]
            write_csv(root/"profile_confidence_intervals.csv",rows)
            with self.assertRaises(ValueError): a.audit(root)
    def test_unchecked_recovered_profile_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); fixture(root)
            for marker in ("RECOVERED_NEEDS_CHECK", "RECOMPUTE_REQUIRED"):
                (root/"K_G"/marker).write_text("check")
                with self.assertRaisesRegex(ValueError,"verification is pending"): a.audit(root)
                (root/"K_G"/marker).unlink()
    def test_better_optimum_and_optimizer_failure_are_flagged(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); fixture(root)
            rows=a.read_csv(root/"alpha/profile.csv"); rows[1]["profile_loglik"]=.2; rows[1]["stopping_met"]="false"
            write_csv(root/"alpha/profile.csv",rows)
            status=next(r["status"] for r in a.audit(root) if r["parameter"]=="alpha")
            self.assertIn("better likelihood than reference",status)
            self.assertIn("nuisance optimization requires review",status)

if __name__=="__main__": unittest.main()
