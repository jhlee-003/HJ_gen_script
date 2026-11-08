#!/usr/bin/env python3
"""Read-only configuration checks; no CMSSW, batch submission or gridpack run.

Run with: python3 2022-PrvGrid/tests/check_private_gridpack.py
Only temporary preflight fixtures are written. Python's standard library and
bash are sufficient; WMCore/FWCore are stubbed for configuration inspection.
"""

import ast
import os
import re
import runpy
import shutil
import string
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch


CAMPAIGN = Path(__file__).resolve().parents[1]
REFERENCE = CAMPAIGN.parent / "2022"
ENV_NAME = "HJ_M-125_TuneCP5_2022_FullSim.env"
FRAGMENT_NAME = "HJ_M-125_TuneCP5_2022_fragment.py"
PRODUCER_NAME = "ProduceHJ_M-125_TuneCP5_2022.sh"
GRIDPACK_NAME = "HJ_el8_amd64_gcc11_CMSSW_13_2_5_HJ_MiNNLO_M125_Private.tgz"
EOS_DESTINATION = "root://eosuser.cern.ch//eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022-PrvGrid/"


def read_assignments(path, expand_shell_names=False):
    """Parse the repository's simple assignments, without executing shell."""
    values = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        value = value.strip().strip('"')
        if expand_shell_names:
            value = string.Template(value).substitute(values)
        values[key.strip()] = value
    return values


def fragment_tree(path):
    return ast.parse(path.read_text(), filename=str(path))


def assignment(tree, name):
    return next(node for node in tree.body if isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == name
                        for target in node.targets))


def fragment_preamble(path, gridpack_path=None, readable=True):
    """Load only path checks and ExternalLHEProducer, with a tiny cms stub."""
    tree = fragment_tree(path)
    producer = assignment(tree, "externalLHEProducer")
    preamble = ast.Module(body=tree.body[:tree.body.index(producer) + 1], type_ignores=[])
    cms = types.ModuleType("FWCore.ParameterSet.Config")
    cms.EDProducer = lambda plugin, **kwargs: types.SimpleNamespace(plugin=plugin, **kwargs)
    cms.vstring = lambda *values: list(values)
    cms.uint32 = int
    cms.string = str
    cms.FileInPath = str
    cms.untracked = types.SimpleNamespace(uint32=int, bool=bool)
    fwcore = types.ModuleType("FWCore")
    parameters = types.ModuleType("FWCore.ParameterSet")
    fwcore.ParameterSet = parameters
    parameters.Config = cms
    modules = {"FWCore": fwcore, "FWCore.ParameterSet": parameters,
               "FWCore.ParameterSet.Config": cms}
    namespace = {"__file__": str(path), "__name__": "gridpack_fragment_test"}
    with patch.dict(sys.modules, modules), patch.dict(os.environ):
        os.environ.pop("GRIDPACK_PATH", None)
        if gridpack_path is not None:
            os.environ["GRIDPACK_PATH"] = str(gridpack_path)
        with patch("os.access", return_value=readable):
            exec(compile(preamble, str(path), "exec"), namespace)
    return namespace["externalLHEProducer"]


class PrivateGridpackChecks(unittest.TestCase):
    def setUp(self):
        self.names = read_assignments(CAMPAIGN / "config" / ENV_NAME, True)
        self.fragment = CAMPAIGN / "config" / FRAGMENT_NAME
        self.producer = CAMPAIGN / PRODUCER_NAME

    def test_syntax_and_unchanged_simulation(self):
        for name in (PRODUCER_NAME, "condor_convert_wrapper.sh", "crab_convert_wrapper.sh"):
            result = subprocess.run(["bash", "-n", str(CAMPAIGN / name)],
                                    text=True, capture_output=True, check=False)
            self.assertEqual(result.returncode, 0, result.stderr)
        for path in (self.fragment, CAMPAIGN / "config" / "HJ_M-125_TuneCP5_2022_crabConfig.py"):
            compile(path.read_text(), str(path), "exec")
        current = self.producer.read_text()
        original = (REFERENCE / PRODUCER_NAME).read_text()
        marker = "cat <<EndOfTestFile"
        self.assertEqual(current.split(marker, 1)[1], original.split(marker, 1)[1],
                         "Generated CMSSW workflow changed")
        for pattern in (r"^\s*--step\s+(.+)$", r"^export SCRAM_ARCH=(.+)$",
                        r"^\s*--conditions\s+(.+)$"):
            self.assertEqual(re.findall(pattern, current, re.M),
                             re.findall(pattern, original, re.M))

    def test_shared_names_and_condor_inputs(self):
        self.assertEqual(self.names["GRIDPACK_FILENAME"], GRIDPACK_NAME)
        self.assertTrue((CAMPAIGN / GRIDPACK_NAME).is_file())
        self.assertGreater((CAMPAIGN / GRIDPACK_NAME).stat().st_size, 0)
        original = read_assignments(REFERENCE / "config" / ENV_NAME, True)
        self.assertEqual(self.names["BASE_TAG"], original["BASE_TAG"] + "__PrvGrid")
        submit = read_assignments(CAMPAIGN / "config" / "HJ_M-125_TuneCP5_2022_condor.sub")
        self.assertEqual(submit["gridpack_file"], GRIDPACK_NAME)
        expanded = submit["transfer_input_files"]
        for _ in range(8):
            new = re.sub(r"\$\(([A-Za-z_][A-Za-z_0-9]*)\)",
                         lambda m: submit.get(m.group(1), m.group(0)), expanded)
            if new == expanded:
                break
            expanded = new
        inputs = [item.strip() for item in expanded.split(",")]
        self.assertIn(GRIDPACK_NAME, inputs)
        for name in inputs:
            self.assertTrue((CAMPAIGN / name).is_file(), "Missing transfer input: " + name)
        expected = self.names["NANOAOD_NAME"] + "__condor-$(Cluster)__job-$(job_number).root"
        self.assertEqual(submit["nano_output_file"], expected)
        self.assertEqual(submit["output_destination"], EOS_DESTINATION)
        self.assertEqual(submit["transfer_output_files"], "$(nano_output_file)")

    def test_crab_gridpack_input_and_campaign(self):
        class Configuration:
            def section_(self, name):
                setattr(self, name, types.SimpleNamespace())

        wmcore = types.ModuleType("WMCore")
        configuration = types.ModuleType("WMCore.Configuration")
        configuration.Configuration = Configuration
        wmcore.Configuration = configuration
        with patch.dict(sys.modules, {"WMCore": wmcore, "WMCore.Configuration": configuration}):
            namespace = runpy.run_path(str(CAMPAIGN / "config" / "HJ_M-125_TuneCP5_2022_crabConfig.py"))
        config = namespace["config"]
        self.assertIn(str(CAMPAIGN / GRIDPACK_NAME), config.JobType.inputFiles)
        self.assertIn("PrvGrid", config.General.requestName)
        self.assertIn("PrvGrid", config.Data.outputDatasetTag)
        self.assertEqual(config.JobType.outputFiles, [self.names["NANOAOD_NAME"] + "__job.root"])
        # CRAB is intentionally still its separate KNU stageout workflow.
        self.assertEqual(config.Site.storageSite, "T3_KR_KNU")

    def test_fragment_absolute_path_and_rejection(self):
        with tempfile.TemporaryDirectory(prefix="hj private gridpack ") as temporary:
            archive = Path(temporary) / GRIDPACK_NAME
            archive.write_bytes(b"fixture only - never executed")
            producer = fragment_preamble(self.fragment, archive)
            self.assertEqual(producer.args, [str(archive.resolve())])
            self.assertTrue(Path(producer.args[0]).is_absolute())
            for value in (None, "", GRIDPACK_NAME, str(Path(temporary) / "absent.tgz")):
                with self.subTest(path=value), self.assertRaises((RuntimeError, ValueError, OSError)):
                    fragment_preamble(self.fragment, value)
            with self.assertRaises((RuntimeError, ValueError, OSError)):
                fragment_preamble(self.fragment, archive, readable=False)

    def test_shower_and_decay_settings_unchanged(self):
        current = fragment_tree(self.fragment)
        original = fragment_tree(REFERENCE / "config" / FRAGMENT_NAME)
        marker = "from Configuration.Generator.Pythia8CommonSettings_cfi"
        self.assertEqual(self.fragment.read_text().split(marker, 1)[1].rstrip("\n"),
                         (REFERENCE / "config" / FRAGMENT_NAME).read_text().split(marker, 1)[1].rstrip("\n"),
                         "Shower/decay source changed")
        for name in ("generator", "ProductionFilterSequence"):
            self.assertEqual(ast.dump(assignment(current, name)), ast.dump(assignment(original, name)))
        # Only the archive argument changes in ExternalLHEProducer.
        def properties(tree):
            return {item.arg: ast.dump(item.value)
                    for item in assignment(tree, "externalLHEProducer").value.keywords if item.arg != "args"}
        self.assertEqual(properties(current), properties(original))

    def test_producer_preflight_without_cmssw(self):
        source = self.producer.read_text()
        marker = "cat <<EndOfTestFile"
        self.assertIn(marker, source)
        # Execute only the prefix, never the generated CMSSW job. This remains
        # safe even if a missing-gridpack check regresses.
        preflight = source.split(marker, 1)[0]
        preflight += '\nprintf "TEST_GRIDPACK_PATH=%s\\n" "${GRIDPACK_PATH:-}"\n'
        with tempfile.TemporaryDirectory(prefix="hj worker path ") as temporary:
            worker = Path(temporary)
            (worker / "config").mkdir()
            shutil.copy2(self.fragment, worker / "config" / FRAGMENT_NAME)
            env_path = worker / ENV_NAME
            shutil.copy2(CAMPAIGN / "config" / ENV_NAME, env_path)
            environment = dict(os.environ, GRIDPACK_PATH="/stale/gridpack.tgz", NTHREADS="2")

            def run():
                return subprocess.run(["bash", "-c", preflight, PRODUCER_NAME, "0", "10", str(env_path)],
                                      cwd=worker, env=environment, text=True,
                                      capture_output=True, check=False)

            missing = run()
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn("gridpack", (missing.stdout + missing.stderr).lower())
            self.assertNotIn("TEST_GRIDPACK_PATH=", missing.stdout)
            archive = worker / GRIDPACK_NAME
            archive.touch()
            empty = run()
            self.assertNotEqual(empty.returncode, 0, "Empty gridpack was accepted")
            archive.write_bytes(b"fixture only - never executed")
            valid = run()
            self.assertEqual(valid.returncode, 0, valid.stdout + valid.stderr)
            self.assertIn("TEST_GRIDPACK_PATH=" + str(archive.resolve()), valid.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
