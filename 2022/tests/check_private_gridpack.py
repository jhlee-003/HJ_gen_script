#!/usr/bin/env python3
"""Check all 2022 mass configurations without CMSSW or job submission.

Run: python3 2022/tests/check_private_gridpack.py
Add --check-gridpacks to require all three real archives. Otherwise absent
payloads are permitted while scripts are being prepared. Present archives are
checked read-only. Mocked wrappers use temporary fixtures, never real CMS tools.
"""

import ast
import os
import re
import runpy
import shutil
import string
import subprocess
import sys
import tarfile
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch


CAMPAIGN = Path(__file__).resolve().parents[1]
MASSES = (120, 125, 130)
REQUIRE_GRIDPACKS = "--check-gridpacks" in sys.argv
if REQUIRE_GRIDPACKS:
    sys.argv.remove("--check-gridpacks")


def producer_name(mass):
    return "ProduceHJ_M-{}_TuneCP5_2022.sh".format(mass)


def config_name(mass, suffix):
    return "HJ_M-{}_TuneCP5_2022_{}".format(mass, suffix)


def gridpack_name(mass):
    return "HJ_el8_amd64_gcc11_CMSSW_13_2_5_HJ_MiNNLO_M{}_Private.tgz".format(mass)


def read_assignments(path, expand=False):
    values = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        value = value.strip().strip('"')
        values[key.strip()] = string.Template(value).substitute(values) if expand else value
    return values


def names_for(mass):
    return read_assignments(CAMPAIGN / "config" / config_name(mass, "FullSim.env"), True)


def assignment(tree, name):
    return next(node for node in tree.body if isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == name for target in node.targets))


def fragment_preamble(path, archive=None, readable=True):
    tree = ast.parse(path.read_text())
    producer = assignment(tree, "externalLHEProducer")
    preamble = ast.Module(body=tree.body[:tree.body.index(producer) + 1], type_ignores=[])
    cms = types.ModuleType("FWCore.ParameterSet.Config")
    cms.EDProducer = lambda plugin, **kw: types.SimpleNamespace(plugin=plugin, **kw)
    cms.vstring = lambda *args: list(args)
    cms.uint32, cms.string, cms.FileInPath = int, str, str
    cms.untracked = types.SimpleNamespace(uint32=int, bool=bool)
    fwcore, parameters = types.ModuleType("FWCore"), types.ModuleType("FWCore.ParameterSet")
    fwcore.ParameterSet, parameters.Config = parameters, cms
    modules = {"FWCore": fwcore, "FWCore.ParameterSet": parameters, "FWCore.ParameterSet.Config": cms}
    namespace = {"__file__": str(path)}
    with patch.dict(sys.modules, modules), patch.dict(os.environ):
        os.environ.pop("GRIDPACK_PATH", None)
        if archive is not None:
            os.environ["GRIDPACK_PATH"] = str(archive)
        with patch("os.access", return_value=readable):
            exec(compile(preamble, str(path), "exec"), namespace)
    return namespace["externalLHEProducer"]


def crab_config(mass):
    class Configuration:
        def section_(self, name):
            setattr(self, name, types.SimpleNamespace())
    wmcore, config_module = types.ModuleType("WMCore"), types.ModuleType("WMCore.Configuration")
    config_module.Configuration, wmcore.Configuration = Configuration, config_module
    with patch.dict(sys.modules, {"WMCore": wmcore, "WMCore.Configuration": config_module}):
        return runpy.run_path(str(CAMPAIGN / "config" / config_name(mass, "crabConfig.py")))["config"]


def fixture(worker, mass, archive=True):
    """Create a flattened worker sandbox using synthetic, non-executable data."""
    names = names_for(mass)
    env_name = config_name(mass, "FullSim.env")
    shutil.copy2(CAMPAIGN / "config" / env_name, worker / env_name)
    (worker / "config").mkdir()
    shutil.copy2(CAMPAIGN / "config" / names["Fragment_filename"], worker / "config" / names["Fragment_filename"])
    if archive:
        (worker / names["GRIDPACK_FILENAME"]).write_bytes(b"synthetic archive - never executed")
    return names, env_name


class MassConfigurationChecks(unittest.TestCase):
    def test_syntax_and_release_matrix(self):
        sources = {}
        expected_steps = ["LHE,GEN,SIM", "DIGI,DATAMIX,L1,DIGI2RAW,HLT:2022v12",
                          "RAW2DIGI,L1Reco,RECO,RECOSIM", "PAT", "NANO"]
        for mass in MASSES:
            with self.subTest(mass=mass):
                path = CAMPAIGN / producer_name(mass)
                sources[mass] = source = path.read_text()
                subprocess.run(["bash", "-n", str(path)], check=True, capture_output=True)
                self.assertTrue(os.access(path, os.X_OK))
                gen_release = "CMSSW_12_4_11_patch3" if mass == 125 else "CMSSW_12_4_14_patch2"
                releases = re.findall(r"^if \[ -r (CMSSW[^/]+)/src", source, re.M)
                self.assertEqual(releases, [gen_release, "CMSSW_12_4_11_patch3", "CMSSW_13_0_13", "CMSSW_13_0_13"])
                self.assertEqual(re.findall(r"^export SCRAM_ARCH=(\S+)", source, re.M),
                                 ["el8_amd64_gcc10"] * 2 + ["el8_amd64_gcc11"] * 2)
                self.assertEqual(re.findall(r"--step ([^\s]+)", source), expected_steps)
                self.assertEqual(re.findall(r"--conditions ([^\s]+)", source),
                                 ["124X_mcRun3_2022_realistic_v12"] * 3 + ["130X_mcRun3_2022_realistic_v5"] * 2)
                self.assertIn("--eventcontent NANOAODSIM", source)
                self.assertNotIn("NANOEDMAODSIM", source)
                for suffix in ("fragment.py", "crabConfig.py"):
                    path = CAMPAIGN / "config" / config_name(mass, suffix)
                    compile(path.read_text(), str(path), "exec")
        reference = sources[125].split("cat <<EndOfTestFile", 1)[1]
        for mass, source in sources.items():
            self.assertEqual(source.split("cat <<EndOfTestFile", 1)[1].replace(
                "CMSSW_12_4_14_patch2", "CMSSW_12_4_11_patch3"), reference)
        for name in ("condor_convert_wrapper.sh", "crab_convert_wrapper.sh"):
            subprocess.run(["bash", "-n", str(CAMPAIGN / name)], check=True, capture_output=True)

    def test_mass_names_transfers_and_stageout(self):
        outputs = set()
        for mass in MASSES:
            with self.subTest(mass=mass):
                names = names_for(mass)
                self.assertEqual(names["HIGGS_MASS"], str(mass))
                self.assertEqual(names["PRODUCER_SCRIPT"], producer_name(mass))
                self.assertEqual(names["Fragment_filename"], config_name(mass, "fragment.py"))
                self.assertEqual(names["GRIDPACK_FILENAME"], gridpack_name(mass))
                self.assertIn("_M-{}_".format(mass), names["BASE_TAG"])
                self.assertTrue(names["BASE_TAG"].endswith("__Run3Summer22__PrvGrid"))
                submit = read_assignments(CAMPAIGN / "config" / config_name(mass, "condor.sub"))
                self.assertEqual(submit["gridpack_file"], gridpack_name(mass))
                inputs = submit["transfer_input_files"].replace("$(gridpack_file)", submit["gridpack_file"])
                inputs = [name.strip() for name in inputs.split(",")]
                expected = [producer_name(mass), "config/" + config_name(mass, "FullSim.env"),
                            "config/" + config_name(mass, "fragment.py"), gridpack_name(mass)]
                self.assertEqual(inputs, expected)
                for name in inputs[:-1]:
                    self.assertTrue((CAMPAIGN / name).is_file())
                self.assertEqual(submit["nano_output_file"], names["NANOAOD_NAME"] + "__condor-$(Cluster)__job-$(job_number).root")
                self.assertEqual(submit["output_destination"],
                                 "root://eosuser.cern.ch//eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022/M{}/".format(mass))
                self.assertEqual(submit["transfer_output_files"], "$(nano_output_file)")
                self.assertEqual(submit["when_to_transfer_output"], "ON_SUCCESS")
                self.assertEqual(submit["events_per_job"], "10")
                self.assertEqual(submit["number_of_jobs"], "1")
                self.assertEqual(submit["number_of_threads"], "2")
                self.assertEqual(submit["+JobBatchName"], "HJ_M{}_2022".format(mass))
                for stream in ("log", "output", "error"):
                    self.assertTrue(submit[stream].startswith("logs/M{}/".format(mass)))
                self.assertIn(config_name(mass, "FullSim.env"), submit["arguments"])
                outputs.add(names["NANOAOD_NAME"])
        self.assertEqual(len(outputs), 3)

    def test_crab_mass_dispatch_and_safe_pilots(self):
        for mass in MASSES:
            with self.subTest(mass=mass):
                cfg, names = crab_config(mass), names_for(mass)
                expected = [CAMPAIGN / "config" / config_name(mass, "FullSim.env"),
                            CAMPAIGN / "config" / config_name(mass, "fragment.py"),
                            CAMPAIGN / producer_name(mass), CAMPAIGN / gridpack_name(mass)]
                self.assertEqual(cfg.JobType.inputFiles, list(map(str, expected)))
                self.assertIn("script=" + producer_name(mass), cfg.JobType.scriptArgs)
                self.assertIn("events=10", cfg.JobType.scriptArgs)
                self.assertEqual(cfg.Data.unitsPerJob, 10)
                self.assertEqual(cfg.Data.totalUnits, 10)
                self.assertEqual(cfg.JobType.outputFiles, [names["NANOAOD_NAME"] + "__job.root"])
                self.assertIn("M-{}_".format(mass), cfg.General.requestName)
                self.assertIn("M{}_".format(mass), cfg.Data.outputDatasetTag)
                self.assertEqual(cfg.Site.storageSite, "T3_KR_KNU")
                self.assertEqual(cfg.JobType.psetName, str(CAMPAIGN / "PSet.py"))

    def test_fragment_paths_and_minnlo_physics(self):
        reference = None
        for mass in MASSES:
            with self.subTest(mass=mass), tempfile.TemporaryDirectory(prefix="hj fragment path ") as temporary:
                path = CAMPAIGN / "config" / config_name(mass, "fragment.py")
                tree = ast.parse(path.read_text())
                generator = assignment(tree, "generator")
                canonical = ast.dump(generator).replace("25:m0 = {}.0".format(mass), "25:m0 = MASS.0")
                reference = canonical if reference is None else reference
                self.assertEqual(canonical, reference)
                constants = {node.value for node in ast.walk(generator) if isinstance(node, ast.Constant) and isinstance(node.value, str)}
                for value in ("25:m0 = {}.0".format(mass), "25:onIfMatch = 23 22", "23:mMin = 50.0",
                              "23:onIfAny = 11 13 15", "SpaceShower:pTmaxMatch = 1", "TimeShower:pTmaxMatch = 1"):
                    self.assertIn(value, constants)
                self.assertNotIn("Pythia8PowhegEmissionVetoSettings", path.read_text())
                archive = Path(temporary) / gridpack_name(mass)
                archive.write_bytes(b"fixture only")
                producer = fragment_preamble(path, archive)
                self.assertEqual(producer.args, [str(archive.resolve())])
                self.assertEqual(producer.nEvents, 5000)
                self.assertTrue(producer.generateConcurrently)
                for value in (None, "", gridpack_name(mass), str(Path(temporary) / "absent" / gridpack_name(mass))):
                    with self.assertRaises(RuntimeError):
                        fragment_preamble(path, value)
                other = Path(temporary) / gridpack_name(130 if mass == 120 else 120)
                other.write_bytes(b"wrong mass fixture")
                with self.assertRaisesRegex(RuntimeError, "matching"):
                    fragment_preamble(path, other)
                with self.assertRaises(RuntimeError):
                    fragment_preamble(path, archive, readable=False)

    def test_producer_preflight_and_generated_script(self):
        for mass in MASSES:
            with self.subTest(mass=mass), tempfile.TemporaryDirectory(prefix="hj worker path ") as temporary:
                worker = Path(temporary)
                names, env_name = fixture(worker, mass, archive=False)
                source = (CAMPAIGN / producer_name(mass)).read_text()
                preflight = source.split("cat <<EndOfTestFile", 1)[0]
                preflight += '\nprintf "EXPORTED_GRIDPACK=%s\\n" "$(printenv GRIDPACK_PATH)"\n'
                environment = dict(os.environ, GRIDPACK_PATH="/stale/gridpack.tgz", NTHREADS="3")
                environment.pop("BASH_ENV", None)

                def run(script, env_file=env_name):
                    return subprocess.run(["bash", "-c", script, producer_name(mass), "0", "10", env_file],
                                          cwd=worker, env=environment, text=True, capture_output=True)

                missing = run(preflight)
                self.assertNotEqual(missing.returncode, 0)
                self.assertIn("gridpack", missing.stderr.lower())
                archive = worker / names["GRIDPACK_FILENAME"]
                archive.touch()
                self.assertNotEqual(run(preflight).returncode, 0)
                archive.write_bytes(b"fixture only - never extracted")
                valid = run(preflight)
                self.assertEqual(valid.returncode, 0, valid.stderr)
                self.assertIn("EXPORTED_GRIDPACK=" + str(archive.resolve()), valid.stdout)
                other_mass = 130 if mass == 120 else 120
                other_env = config_name(other_mass, "FullSim.env")
                shutil.copy2(CAMPAIGN / "config" / other_env, worker / other_env)
                mixed = run(preflight, other_env)
                self.assertEqual(mixed.returncode, 3)
                self.assertIn("matching", mixed.stderr)
                # Include the heredoc but STOP before chmod/execution. No CMS
                # command is executed: only the generated file is inspected.
                generation = source.split('\necho "Made ${TAG}_cmd.sh"', 1)[0]
                self.assertTrue(generation.rstrip().endswith("EndOfTestFile"))
                result = run(generation)
                self.assertEqual(result.returncode, 0, result.stderr)
                generated = worker / "job_scripts" / (names["BASE_TAG"] + "__job-1_cmd.sh")
                self.assertTrue(generated.is_file())
                subprocess.run(["bash", "-n", str(generated)], check=True, capture_output=True)
                text = generated.read_text()
                self.assertEqual(text.count("--nThreads 3"), 5)
                self.assertIn("-n 10", text)
                self.assertIn("firstRun=cms.untracked.uint32(1)", text)
                self.assertEqual(len(re.findall(r"^cmsRun -j FrameworkJobReport.xml", text, re.M)), 5)
                self.assertFalse(any(path.name.startswith("CMSSW_") for path in worker.iterdir()))

    def test_wrappers_dispatch_each_mass_with_mocked_tools(self):
        fake_producer = '''#!/bin/bash
set -euo pipefail
source "$3"
printf '%s,%s,%s,%s\\n' "$HIGGS_MASS" "$1" "$2" "$NTHREADS" > dispatch.txt
out="${NANOAOD_NAME}__job-$(($1 + 1)).root"
printf 'synthetic output\\n' > "$out"
printf '<File>%s</File>\\n' "$out" > FrameworkJobReport.xml
'''
        for mass in MASSES:
            for scheduler in ("condor", "crab"):
                with self.subTest(mass=mass, scheduler=scheduler), tempfile.TemporaryDirectory(prefix="hj wrapper ") as temporary:
                    worker = Path(temporary)
                    names, env_name = fixture(worker, mass)
                    (worker / names["PRODUCER_SCRIPT"]).write_text(fake_producer)
                    tools_dir = worker / "bin"
                    tools_dir.mkdir()
                    for tool, body in {
                        "root": "#!/bin/bash\nprintf '10\\n'\n",
                        "cmsRun": "#!/bin/bash\nprintf '<Bootstrap/>\\n' > FrameworkJobReport.xml\n",
                    }.items():
                        path = tools_dir / tool
                        path.write_text(body)
                        path.chmod(0o755)
                    proxy = worker / "mock-proxy"
                    proxy.write_text("test fixture, not a credential")
                    environment = dict(os.environ, X509_USER_PROXY=str(proxy), PATH=str(tools_dir) + os.pathsep + os.environ["PATH"])
                    environment.pop("BASH_ENV", None)
                    wrapper = CAMPAIGN / (scheduler + "_convert_wrapper.sh")
                    if scheduler == "condor":
                        output = names["NANOAOD_NAME"] + "__condor-99__job-13.root"
                        args = ["2", "10", env_name, "3", "10", output]
                        expected_index = 12
                    else:
                        output = names["NANOAOD_NAME"] + "__job.root"
                        args = ["12", "script=" + names["PRODUCER_SCRIPT"], "events=10", "names=" + env_name, "threads=3"]
                        expected_index = 12
                    result = subprocess.run(["bash", str(wrapper)] + args, cwd=worker,
                                            env=environment, text=True, capture_output=True)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertEqual((worker / "dispatch.txt").read_text().strip(),
                                     "{},{},10,3".format(mass, expected_index))
                    self.assertTrue((worker / output).is_file())
                    self.assertFalse((worker / (names["NANOAOD_NAME"] + "__job-13.root")).exists())
                    if scheduler == "crab":
                        self.assertIn(output, (worker / "FrameworkJobReport.xml").read_text())

    def test_crab_rejects_mismatched_script_before_bootstrap(self):
        with tempfile.TemporaryDirectory(prefix="hj mismatch ") as temporary:
            worker = Path(temporary)
            names, env_name = fixture(worker, 125)
            wrong = producer_name(120)
            (worker / wrong).write_text("#!/bin/bash\nexit 99\n")
            result = subprocess.run(["bash", str(CAMPAIGN / "crab_convert_wrapper.sh"), "0",
                                     "script=" + wrong, "events=10", "names=" + env_name, "threads=2"],
                                    cwd=worker, text=True, capture_output=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn("does not match", result.stderr)
            self.assertFalse((worker / "FrameworkJobReport.xml").exists())

    def test_gridpack_payloads_when_present(self):
        for mass in MASSES:
            with self.subTest(mass=mass):
                path = CAMPAIGN / gridpack_name(mass)
                if not path.exists() and not REQUIRE_GRIDPACKS:
                    continue
                self.assertTrue(path.is_file(), "Required gridpack is missing: " + str(path))
                self.assertGreater(path.stat().st_size, 0)
                with tarfile.open(str(path), "r:gz") as archive:
                    self.assertTrue(archive.getmember("pwhg_main").mode & 0o111)
                    self.assertTrue(archive.getmember("runcmsgrid.sh").mode & 0o111)
                    card = archive.extractfile("powheg.input").read().decode(errors="replace")
                    masses = []
                    for line in card.splitlines():
                        fields = line.split("!", 1)[0].split("#", 1)[0].split()
                        if fields and fields[0].lower() == "hmass":
                            masses.append(float(fields[1].lower().replace("d", "e")))
                    self.assertTrue(masses, "No hmass setting in gridpack card")
                    self.assertEqual(set(masses), {float(mass)})


if __name__ == "__main__":
    absent = ["M{}".format(mass) for mass in MASSES if not (CAMPAIGN / gridpack_name(mass)).exists()]
    if absent and not REQUIRE_GRIDPACKS:
        print("Payloads not installed yet: " + ", ".join(absent) + "; configuration checks use temporary fixtures.", flush=True)
    unittest.main(verbosity=2)
