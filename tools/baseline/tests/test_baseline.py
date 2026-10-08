#!/usr/bin/env python3
"""Offline fixtures: no CERN access, no actual Condor submission."""
import array
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

HELPER = Path(__file__).resolve().parents[1]
try:
    import ROOT
    ROOT.gROOT.SetBatch(True)
except ImportError:
    ROOT = None


def embedded(script, marker):
    return script.split("<<'" + marker + "'\n", 1)[1].split("\n" + marker + "\n", 1)[0]


class BaselineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_shell_syntax(self):
        for script in HELPER.glob("*.sh"):
            subprocess.run(["bash", "-n", str(script)], check=True)

    def test_condor_transfer_and_failure_logs(self):
        text = (HELPER / "baseline.sub").read_text()
        for setting in ("when_to_transfer_output = ON_SUCCESS",
                        "transfer_output_files = $(pico_name)",
                        "transfer_executable = False", "/cmssw/el9:x86_64",
                        "on_exit_hold = (ExitBySignal == true) || (ExitCode != 0)"):
            self.assertIn(setting, text)
        self.assertNotIn("root://", next(line for line in text.splitlines() if line.startswith("error =")))
        self.assertNotIn("tar.gz", text)
        self.assertNotIn("stream_output =", text)
        self.assertNotIn("stream_error =", text)
        self.assertNotIn("export XrdSecPROTOCOL=gsi", (HELPER / "setup_runtime.sh").read_text())

    def test_json_directory_key_adaptation(self):
        code = embedded((HELPER / "run_baseline.sh").read_text(), "NORM")
        source, output = self.base / "global.json", self.base / "worker.json"
        weights = {"genEventSumw": 13, **{f"LHEScaleSumw{i}": 3 for i in range(9)}}
        source.write_text(json.dumps({"dataset": "Higgs", "weights": weights}))
        argv = ["-", str(source), "/scratch/NanoAODv12/nano/2022/mc",
                "Higgs__Run3Summer22NanoAODv12__abc.root", str(output)]
        with patch.object(sys, "argv", argv):
            exec(compile(code, "<worker normalization>", "exec"), {})
        self.assertEqual(json.loads(output.read_text()), {argv[2]: {"Higgs": weights}})
        argv[3] = "WrongDataset__abc.root"
        with patch.object(sys, "argv", argv), self.assertRaises(SystemExit):
            exec(compile(code, "<worker normalization>", "exec"), {})

    def prepare_fixture(self, depth=2, limit=1, existing=False, bad=None, bad_index=0):
        inputs, output, run, helper = [self.base / name for name in ("nano", "pico", "run", "helper")]
        for path in (inputs, output, run, helper):
            path.mkdir(exist_ok=True)
        (helper / "baseline.sub").write_text((HELPER / "baseline.sub").read_text())
        files = [inputs / "a.root", inputs / "child/b.root", inputs / "child/deeper/c.root"]
        records = {}
        for index, path in enumerate(files):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
            values = [5] if index == 0 else [-2, 10]
            entries = [SimpleNamespace(genEventSumw=value, LHEScaleSumw=[1.] * 9) for value in values]
            if bad == "scale" and index == bad_index:
                entries[0].LHEScaleSumw = [1.] * 8
            if bad == "nan" and index == bad_index:
                entries[0].genEventSumw = float("nan")
            missing = bad if index == bad_index else None
            class Runs(list):
                def GetEntries(self): return len(self)
                def GetBranch(self, name, missing=missing): return name != missing
            runs = Runs(entries)
            records["root://fixture/" + str(path)] = SimpleNamespace(
                Get=lambda name, runs=runs: runs if name == "Runs" else SimpleNamespace(GetEntries=lambda: 7),
                IsZombie=lambda: False, Close=lambda: None)
        self.opened = []
        fake = SimpleNamespace(gROOT=SimpleNamespace(SetBatch=lambda flag: None),
                               TFile=SimpleNamespace(Open=lambda url, mode: self.opened.append(url) or records[url]))
        if existing:
            (output / "old.root").touch()
        argv = ["-", "2022", "ggH", str(inputs), "root://fixture", str(depth),
                "GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8",
                str(output), str(run), str(helper), str(limit)]
        code = embedded((HELPER / "submit_baseline.sh").read_text(), "PREPARE")
        with patch.dict(sys.modules, ROOT=fake), patch.object(sys, "argv", argv), redirect_stdout(io.StringIO()):
            exec(compile(code, "<sample preparation>", "exec"), {})
        return run

    def test_normalize_selected_files_only(self):
        run = self.prepare_fixture()
        norm = json.loads((run / "normalization.json").read_text())
        self.assertEqual(norm["weights"]["genEventSumw"], 5)
        self.assertEqual([norm["weights"][f"LHEScaleSumw{i}"] for i in range(9)], [1] * 9)
        inputs = (run / "inputs.txt").read_text().splitlines()
        self.assertEqual(len(inputs), 1)
        self.assertEqual(self.opened, inputs)
        rows = (run / "jobs.tsv").read_text().splitlines()
        self.assertEqual(len(rows), 1)
        self.assertEqual([row.split("\t")[0] for row in rows], inputs)
        self.assertIn("Run3Summer22NanoAODv12", rows[0])
        self.assertIn("pico_baseline_ggH2022_", rows[0])

    def test_normalize_all_files_without_limit(self):
        run = self.prepare_fixture(limit=0)
        norm = json.loads((run / "normalization.json").read_text())
        # Includes the second file's two Runs entries, including a negative sum.
        self.assertEqual(norm["weights"]["genEventSumw"], 13)
        self.assertEqual([norm["weights"][f"LHEScaleSumw{i}"] for i in range(9)], [3] * 9)
        inputs = (run / "inputs.txt").read_text().splitlines()
        rows = (run / "jobs.tsv").read_text().splitlines()
        self.assertEqual(len(inputs), 2)
        self.assertEqual(self.opened, inputs)
        self.assertEqual([row.split("\t")[0] for row in rows], inputs)

    def test_limit_larger_than_available_inputs(self):
        run = self.prepare_fixture(limit=100)
        self.assertEqual(len((run / "jobs.tsv").read_text().splitlines()), 2)
        self.assertEqual(len(self.opened), 2)
        self.assertEqual(json.loads((run / "normalization.json").read_text())["weights"]["genEventSumw"], 13)

    def test_invalid_metadata_in_unselected_file_is_not_read(self):
        run = self.prepare_fixture(bad="scale", bad_index=1)
        self.assertEqual(len(self.opened), 1)
        self.assertEqual(json.loads((run / "normalization.json").read_text())["weights"]["genEventSumw"], 5)

    def test_private_depth_one_excludes_subdirectories(self):
        run = self.prepare_fixture(depth=1)
        self.assertEqual(len((run / "inputs.txt").read_text().splitlines()), 1)
        self.assertEqual(json.loads((run / "normalization.json").read_text())["weights"]["genEventSumw"], 5)

    def test_existing_outputs_are_not_overwritten(self):
        with self.assertRaisesRegex(SystemExit, "already contains"):
            self.prepare_fixture(existing=True)

    def test_invalid_scale_metadata_is_not_assumed(self):
        with self.assertRaisesRegex(SystemExit, "Invalid Runs"):
            self.prepare_fixture(bad="scale")

    def test_nonfinite_metadata_rejected(self):
        with self.assertRaisesRegex(SystemExit, "Invalid Runs"):
            self.prepare_fixture(bad="nan")

    def test_missing_normalization_branch_rejected(self):
        with self.assertRaisesRegex(SystemExit, "Missing genEventSumw"):
            self.prepare_fixture(bad="genEventSumw")

    def test_setup_clean_shell_without_home(self):
        source = self.base / "cmsset.sh"
        source.write_text("export PATH=/usr/bin:/bin\n")
        cmssw = self.base / "CMSSW_15_0_17/src"
        cmssw.mkdir(parents=True)
        script = (HELPER / "setup_runtime.sh").read_text()
        script = script.replace("/cvmfs/cms.cern.ch/cmsset_default.sh", str(source))
        script = script.replace("/afs/cern.ch/user/j/junhyuk/CMSSW_15_0_17/src", str(cmssw))
        setup = self.base / "setup_runtime.sh"
        setup.write_text(script)
        command = ('set -eu; scramv1() { printf "%s\\n" "export CMSSW_VERSION=CMSSW_15_0_17"; }; '
                   'before="$PWD"; source ' + shlex.quote(str(setup)) +
                   '; test "$PWD" = "$before"; test "$SCRAM_ARCH" = el9_amd64_gcc12; '
                   'test "$CMSSW_VERSION" = CMSSW_15_0_17; test "$SET_ENV_PATH" = ' + shlex.quote(str(setup)))
        subprocess.run(["bash", "-c", command], env={"PATH": "/usr/bin:/bin"}, check=True, cwd=self.base)

    @unittest.skipIf(ROOT is None, "PyROOT unavailable")
    def test_bitmap_and_all_weight_branches_preserved(self):
        raw, selected = self.base / "raw.root", self.base / "selected.root"
        output = ROOT.TFile.Open(str(raw), "RECREATE")
        tree = ROOT.TTree("tree", "tree")
        use, bits, weight, lumi, extra = (array.array("i", [1]), array.array("i", [0]),
                                         array.array("d", [2.]), array.array("d", [.5]), array.array("i", [42]))
        for name, value, leaf in (("use_event", use, "I"), ("zg_cutBitMap", bits, "I"),
                                  ("weight", weight, "D"), ("w_lumi", lumi, "D"), ("extra", extra, "I")):
            tree.Branch(name, value, name + "/" + leaf)
        for flag in (1, 0):
            use[0] = flag
            for bit in range(4096):
                bits[0] = bit
                tree.Fill()
        tree.Write()
        output.Close()
        code = embedded((HELPER / "run_baseline.sh").read_text(), "SKIM")
        with patch.object(sys, "argv", ["-", str(raw), str(selected)]), redirect_stdout(io.StringIO()):
            exec(compile(code, "<baseline selection>", "exec"), {})
        result = ROOT.TFile.Open(str(selected))
        chosen = result.Get("tree")
        self.assertEqual(chosen.GetEntries(), 4)
        self.assertEqual(sorted(int(row.zg_cutBitMap) for row in chosen), [3070, 3071, 3582, 3583])
        self.assertTrue(chosen.GetBranch("extra"))
        self.assertEqual([float(row.weight) for row in chosen], [2.] * 4)
        self.assertEqual([float(row.w_lumi) for row in chosen], [.5] * 4)
        result.Close()

    def worker_fixture(self, converter_exit=0, zero=False):
        helper = self.base / "helper"
        converter = self.base / "nano2pico/run"
        sandbox = self.base / "sandbox"
        for directory in (helper, converter, sandbox):
            directory.mkdir(parents=True)
        (helper / "run_baseline.sh").write_text((HELPER / "run_baseline.sh").read_text())
        bin_dir = self.base / "bin"
        bin_dir.mkdir()
        python = bin_dir / "python3"
        python.write_text("#!/bin/bash\nexec " + shlex.quote(sys.executable) + ' "$@"\n')
        python.chmod(0o755)
        (helper / "setup_runtime.sh").write_text(
            "export NANO2PICO_DIR=" + shlex.quote(str(converter.parent)) +
            "\nexport PATH=" + shlex.quote(str(bin_dir)) + ':$PATH\n')
        program = converter / "process_nano.exe"
        if converter_exit:
            program.write_text(f"#!/bin/bash\necho fixture-converter-failure >&2\nexit {converter_exit}\n")
        else:
            program.write_text("""#!/usr/bin/env python3
import array, json, sys
from pathlib import Path
import ROOT
args = dict(zip(sys.argv[1::2], sys.argv[2::2]))
assert args["--nent"] == "-1"
assert "NanoAODv12" in args["--in_dir"]
assert Path.cwd() == Path(__file__).resolve().parents[1]
norm = json.loads(Path(args["--norm"]).read_text())
assert norm[args["--in_dir"]]["Higgs"]["genEventSumw"] == 13
raw = Path(args["--out_dir"]) / "raw_pico" / ("raw_pico_" + args["--in_file"])
file = ROOT.TFile.Open(str(raw), "RECREATE")
tree = ROOT.TTree("tree", "tree")
use, bits, weight, lumi = array.array("i", [1]), array.array("i", [3582]), array.array("d", [-2.]), array.array("d", [.5])
for name, value, kind in [("use_event", use, "I"), ("zg_cutBitMap", bits, "I"), ("weight", weight, "D"), ("w_lumi", lumi, "D")]:
    tree.Branch(name, value, name+"/"+kind)
for flag in FLAGS:
    use[0] = flag
    tree.Fill()
tree.Write()
file.Close()
""".replace("FLAGS", "[0, 0]" if zero else "[1, 0, 1]"))
        program.chmod(0o755)
        (sandbox / "uuid.root").write_bytes(b"NanoAOD fixture")
        (sandbox / "normalization.json").write_text(json.dumps({
            "dataset": "Higgs", "weights": {"genEventSumw": 13, **{f"LHEScaleSumw{i}": 3 for i in range(9)}}}))
        return helper, sandbox

    @unittest.skipIf(ROOT is None, "PyROOT unavailable")
    def test_complete_worker_direct_conversion_and_filter(self):
        helper, sandbox = self.worker_fixture()
        result = subprocess.run(["bash", str(helper / "run_baseline.sh"), "2022", "uuid.root",
                                 "Higgs__Run3Summer22NanoAODv12__abc.root", "pico.root", "normalization.json"],
                                cwd=sandbox, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Baseline events: 2", result.stdout)
        self.assertFalse(list(sandbox.glob("baseline_work.*")))
        pico = ROOT.TFile.Open(str(sandbox / "pico.root"))
        self.assertEqual(pico.Get("tree").GetEntries(), 2)
        pico.Close()

    @unittest.skipIf(ROOT is None, "PyROOT unavailable")
    def test_zero_baseline_events_is_valid_output(self):
        helper, sandbox = self.worker_fixture(zero=True)
        result = subprocess.run(["bash", str(helper / "run_baseline.sh"), "2022", "uuid.root",
                                 "Higgs__Run3Summer22NanoAODv12__abc.root", "pico.root", "normalization.json"],
                                cwd=sandbox, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        pico = ROOT.TFile.Open(str(sandbox / "pico.root"))
        self.assertEqual(pico.Get("tree").GetEntries(), 0)
        pico.Close()

    def test_worker_preserves_converter_failure_status(self):
        helper, sandbox = self.worker_fixture(converter_exit=134)
        result = subprocess.run(["bash", str(helper / "run_baseline.sh"), "2022", "uuid.root",
                                 "Higgs__Run3Summer22NanoAODv12__abc.root", "pico.root", "normalization.json"],
                                cwd=sandbox, capture_output=True, text=True)
        self.assertEqual(result.returncode, 134, result.stderr)
        self.assertIn("exit=134", result.stderr)
        self.assertFalse((sandbox / "pico.root").exists())
        self.assertFalse(list(sandbox.glob("baseline_work.*")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
