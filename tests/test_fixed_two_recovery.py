import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import csv

path=Path(__file__).resolve().parents[1]/"scripts/recover_fixed_two_intervals.py"
spec=importlib.util.spec_from_file_location("recovery",path)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name).resolve()
        self.repo=self.base/"repo";self.source=self.base/"failed";self.out=self.base/"recovered"
        (self.repo/"profiles").mkdir(parents=True);self.source.mkdir()
        (self.repo/module.RUNNER).write_text("corrected runner")
        (self.repo/module.HELPER).write_text("isolated seeds")
        self.input=self.repo/"observations.csv";self.input.write_text("unchanged experimental observations")
        for f in module.ROOT_FILES:
            (self.source/f).write_text('full_range_repair = false\nprofile_scope = "fixed-two"\n' if f=="settings.toml" else "saved numerical result")
        (self.source/"mle_solution.bin").write_bytes(b"old closure must not be copied")
        for name in ("K_G","K_tau"):
            (self.source/name).mkdir()
            for f in module.PROFILE_FILES:
                (self.source/name/f).write_text(name+f)
        (self.source/"K_d").mkdir();(self.source/"K_d"/"partial.csv").write_text("not complete")
        with (self.source/"source_hashes.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=["path","sha256"]);w.writeheader()
            w.writerows([dict(path=str(self.repo/module.RUNNER),sha256=module.LEGACY_RUNNER_SHA),dict(path=str(self.input),sha256=module.digest(self.input))])

    def snapshot(self):
        return {str(p.relative_to(self.source)):module.digest(p) for p in self.source.rglob("*") if p.is_file()}

    def test_exact_copy_source_preservation_and_recheck_markers(self):
        before=self.snapshot()
        module.recover(self.source,self.out,self.repo)
        self.assertEqual(before,self.snapshot())
        for n in ("K_G","K_tau"):
            for f in module.PROFILE_FILES:
                self.assertEqual(module.digest(self.source/n/f),module.digest(self.out/n/f))
            self.assertTrue((self.out/n/"RECOVERED_NEEDS_CHECK").is_file())
        self.assertFalse((self.out/"K_d").exists())
        self.assertFalse((self.out/"mle_solution.bin").exists())
        manifest=module.read_csv(self.out/"source_hashes.csv")
        self.assertEqual(len(manifest),3)
        self.assertEqual(manifest[1]["path"],str(self.repo/module.HELPER))
        self.assertEqual(json.loads((self.out/"recovery.json").read_text())["source_snapshot"]["K_G/profile.csv"],before["K_G/profile.csv"])
        with self.assertRaisesRegex(ValueError,"already exists"):
            module.recover(self.source,self.out,self.repo)

    def test_changed_inputs_are_rejected(self):
        self.input.write_text("changed data")
        with self.assertRaisesRegex(ValueError,"Protected input changed"):
            module.recover(self.source,self.out,self.repo)
        self.assertFalse(self.out.exists())

    def test_unknown_runner_and_overlapping_output_are_rejected(self):
        f=self.source/"source_hashes.csv"
        f.write_text(f.read_text().replace(module.LEGACY_RUNNER_SHA,"unrecognized"))
        with self.assertRaisesRegex(ValueError,"known failed runner"):
            module.recover(self.source,self.out,self.repo)
        with self.assertRaisesRegex(ValueError,"separate"):
            module.recover(self.source,self.source/"new",self.repo)

    def test_full_range_run_cannot_be_used_as_failed_interval_source(self):
        if module.tomllib is None:
            self.skipTest("Julia validates settings when Python lacks tomllib")
        (self.source/"settings.toml").write_text('full_range_repair = true\nprofile_scope = "fixed-two"\n')
        with self.assertRaisesRegex(ValueError,"confidence-interval"):
            module.recover(self.source,self.out,self.repo)


if __name__=="__main__":
    unittest.main()
