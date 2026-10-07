# NanoAOD -> nano2pico -> H-to-Zgamma baseline (lxplus)

This workflow uses the existing checkout at
`/afs/cern.ch/user/j/junhyuk/nano2pico_sequoia_v1`. It prepares a normal CERN
HTCondor DAG; it does not use UCSB's `auto_submit_jobs.py`, modify the physics
sources, submit jobs automatically, or modify the NanoAOD generation/plotters.
All workflow sources and generated campaign files live under `tools/baseline`.

## Inputs and outputs

| Sample | Input discovery | Final baseline pico directory |
| --- | --- | --- |
| HJ | ALL top-level `.root` files in `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022` | `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico` |
| ggH (central) | ALL `.root` files at maxdepth 2 under `/eos/cms/store/mc/Run3Summer22NanoAODv12/GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2` | `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/ggH2022pico` |

No Condor-cluster filtering is applied. `maxdepth 2` includes the dataset root
and its immediate child directories, not deeper directories or `logs/` in HJ.
Neither input sample is copied to AFS. Each conversion job stages ONE NanoAOD
into worker scratch, runs the converter, filters the raw pico, and transfers
only the baseline pico to EOS. Raw picos stay in worker scratch and disappear
with the job sandbox. No merge or branch slimming is performed.

## Physics definition

The [upstream README](https://github.com/richstu/nano2pico/blob/htozgamma_sequoia_v1/README.md)
describes conversion, normalization, skimming, and optional slimming. Its `llg`
skim only requires `nll >= 1 && nphoton > 0`, NOT the full baseline.

We convert without `--skim` or `--nent`, then select the baseline encoded by
[`src/zgamma_producer.cpp`](https://github.com/richstu/nano2pico/blob/1403ed07a4220457a6dbd7a218b98a10c081e6e5/src/zgamma_producer.cpp):

* Reconstructed OS same-flavor dilepton and selected photon, using nano2pico's
  object definitions, corrections, and candidate reconstruction.
* Relevant electron/muon triggers AND their lepton-pT turn-on requirements.
* `80 <= m(ll) <= 100 GeV`.
* `pT(gamma) / m(llgamma) >= 15/110`.
* `m(ll) + m(llgamma) > 185 GeV`.
* `100 <= m(llgamma) <= 180 GeV`.
* Event filters (`pass`).

The shared implementation is:

```text
nll >= 1 && nphoton >= 1 &&
(zg_cutBitMap & 2558) == 2558 && (zg_cutBitMap & 1536) != 0
```

Bits 11 and 8..1 must be set; either channel bit 10 (ee) or 9 (mumu) must
be set. Bit 0 is the data-blinding flag and is deliberately NOT required for
these MC samples: the 120--130 GeV Higgs region remains. No ggF category,
BDT, or other category-specific cuts are applied; ggF production is not the
same thing as the ggF reconstructed category. All pico branches, including
weights and systematic variations, are retained; acceptance is nominal only.

The source SHA256 is checked against the reviewed branch implementation.
If your checkout differs, preparation stops; do not bypass this by changing
the hash without reviewing the baseline definition. The manifest records the
actual checkout revision, source and binary hashes, and dirty status.

Normalization uses ALL discovered files in each sample independently, even
for a one-file pilot or `--skip-existing`. It reproduces the schema and sums
from upstream `scripts/find_normalization.py`: `genEventSumw` and nine
`LHEScaleSumw` sums over every Runs entry. Missing, nonfinite, zero, or
non-nine LHE scale metadata is an error, not an assumed default. Never
normalize using only baseline-passing events. This preserves upstream weight
behavior; it is not a new prescription for LHE uncertainty normalization.
Before using `weight` for yield/systematics studies, review upstream's ggF
NNLO reweighting and scale-weight conventions for the private MiNNLO sample.
The baseline selection itself does not depend on these event weights.

## 1. Build once with a matching runtime

The upstream `set_env.sh` and `SConstruct` assume UCSB and SL7 dependencies.
Our helper bypasses `set_env.sh`, but does not rewrite `SConstruct`: it uses
the matching `slc7_amd64_gcc12` / `CMSSW_14_2_2` runtime and correctionlib
2.6.4. This is the **analysis runtime**, not the generation runtime.
Condor workers use the corresponding CC7 container. Do not build with a
different ROOT/CMSSW runtime and then mix it with this worker configuration.

On lxplus, enter CMS's CC7 container:

```bash
/cvmfs/cms.cern.ch/common/cmssw-cc7
```

Inside it:

```bash
cd ~/HJ_gen_script/tools/baseline
source setup_runtime.sh

# Only if scons is not installed in this Python environment:
python3 -m pip install --user scons

bash build_nano2pico.sh
exit
```

The helper builds the existing checkout; it does not clone, pull, or switch
its branch. Ensure the checkout contains the intended
`htozgamma_sequoia_v1` physics sources before building. If it was already built
in this exact runtime and sources have not changed, skip the build. Preparation
packages the real ELF binary (not the kernel-dependent `run/` shell launcher)
and `data/` directory. If multiple kernel builds exist, use `--binary` to
select `.../kernel/<build-kernel>/run/process_nano.exe` explicitly.

## 2. Proxy on AFS

Run on the normal lxplus host. Do not use an EosSubmit-only schedd: workflow
inputs and DAG logs are on AFS, while successful pico outputs use the XRootD
plugin via `output_destination`.

```bash
mkdir -p "$HOME/tmp"
chmod 700 "$HOME/tmp"
export X509_USER_PROXY="$HOME/tmp/x509up"
voms-proxy-init --voms cms --valid 168:00 --out "$X509_USER_PROXY"
voms-proxy-info --file "$X509_USER_PROXY" --timeleft
```

Do not use `/tmp/x509up_*`: the remote schedd cannot read the login host's
local `/tmp`. Preparation embeds the resolved AFS proxy path in submit files;
renew the proxy at that same path before it expires. Do not commit a proxy.

## 3. Prepare a one-file pilot for EACH sample

```bash
cd ~/HJ_gen_script/tools/baseline
python3 prepare_baseline.py 2022 --sample both --campaign pilot_2022 --limit 1
cd runs/pilot_2022
```

Preparation lists the EOS inputs, creates the requested output directories
and their `logs/` subdirectories, snapshots the converter/helpers, and writes
`manifest.json`, two normalization submit files, two conversion submit files,
and `workflow.dag`. It does NOT submit. Each sample gets one normalization
parent job; conversion+baseline child jobs start only after its normalization
JSON has returned successfully to the campaign directory.

Basic checks before submitting:

```bash
condor_submit -dry-run normalization_HJ.classad normalize_HJ.sub
condor_submit -dry-run normalization_ggH.classad normalize_ggH.sub
condor_submit_dag -no_submit workflow.dag
```

`-no_submit` checks DAG generation, not EOS I/O or physics execution. Conversion
submit files depend on normalization JSONs that do not exist until the parent
jobs complete; do not create dummy normalization files just to dry-run them.

Then submit the pilot:

```bash
condor_submit_dag -maxjobs 100 workflow.dag
condor_q
```

The pilot has four execution jobs (two normalizations + one conversion per
sample), plus DAGMan. Even the pilot normalizes all inputs. The output names
are deterministic, `pico_baseline_HJ2022_<input-hash>.root` and
`pico_baseline_ggH2022_<input-hash>.root`, not new timestamps/cluster IDs.
Original UUID/source filenames are recorded in `baseline_provenance`.

## 4. Full production after the pilot succeeds

```bash
cd ~/HJ_gen_script/tools/baseline
python3 prepare_baseline.py 2022 --sample both --campaign production_2022 --skip-existing
cd runs/production_2022
condor_submit_dag -maxjobs 100 workflow.dag
```

Omit `--limit` for all files. `--skip-existing` avoids regenerating pilot
outputs, but normalization still includes them. Preparation otherwise refuses
existing final filenames. Do not submit duplicate/concurrent campaigns for
the same inputs: they have the same EOS output names. The input set is frozen
in the manifest; adding more NanoAODs changes global normalization, so outputs
from different input sets should not be mixed without renormalization.
`--sample HJ` or `--sample ggH` processes just one sample. Only 2022 is configured.

## Monitoring and validation

Scheduler `.log` files and successful normalization stdout/stderr are under
`runs/<campaign>/logs` on AFS. Successful conversion stdout/stderr go into
`logs/` under the respective EOS pico directory. Worker stdout/stderr are
combined in the `.out` log. On a worker failure the wrapper explicitly uploads
`logs/baseline_failure_<sample>_<mode>_<task>_<cluster.proc>.log` to that sample's
EOS pico directory. This is best-effort: proxy/network or runtime-setup failures
can prevent uploading diagnostics. Failed jobs are held on their nonzero exit
code, without attempting to transfer a nonexistent pico:

```bash
condor_q -hold -af ClusterId ProcId HoldReasonCode HoldReason ExitCode
# While a worker is still running, stream its combined stdout:
condor_tail CLUSTER.PROCESS
```

For an exited/held worker, read the failure log in the appropriate EOS `logs/`
directory. No automatic retries are configured. Inspect the actual error; do not release
a completed/transfer-failed job blindly and accidentally regenerate its output.
DAGMan rescue files support resuming incomplete nodes once failures are resolved.

A valid zero-entry output is allowed if no events pass. Workers check that
raw pico entries match NanoAOD entries, that saved entries match the baseline
count, and that every saved event passes the same baseline. On lxplus, inspect
an actual final file with a compatible ROOT environment:

```bash
root -l '/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico/ACTUAL_OUTPUT.root'
```

At the ROOT prompt:

```cpp
auto t = (TTree*)_file0->Get("tree");
t->GetEntries();
t->GetEntries("nll >= 1 && nphoton >= 1 && (zg_cutBitMap & 2558) == 2558 && (zg_cutBitMap & 1536) != 0");
((TNamed*)_file0->Get("baseline_provenance"))->GetTitle();
```

The two entry counts must agree. This confirms baseline membership, not a
complete physics validation or a cross-section/weighting validation.

Local tests (standard Python; PyROOT tests run additionally when available):

```bash
python3 tests/test_baseline.py
```
