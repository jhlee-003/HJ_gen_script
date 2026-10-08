# NanoAOD to baseline pico on lxplus

The workflow follows the working DY structure, without truth matching:

| File | Role |
| --- | --- |
| `setup_runtime.sh` | Load CMSSW_15_0_17 / el9_amd64_gcc12 and the existing nano2pico paths. |
| `submit_baseline.sh` | Choose inputs, calculate common normalization for those selected files once, prepare one job per file, optionally submit. |
| `baseline.sub` | Condor resources, input transfer, AFS logs and successful pico transfer to EOS. |
| `run_baseline.sh` | Stage one file, run process_nano directly, select baseline events, save the pico. |

All generated lists, normalization JSONs, submit files and logs are under
`tools/baseline/runs/CAMPAIGN/SAMPLE/`. No DAG, tarball packaging, checkout
hash guards, provenance records, timestamp labels or extra Python worker.

## Inputs and outputs

Only 2022 is configured. There is no Condor-cluster filtering.

| Sample | Inputs | Output |
| --- | --- | --- |
| HJ | All top-level ROOT files in `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022` (maxdepth 1) | `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico` |
| ggH | All ROOT files at maxdepth 2 under `/eos/cms/store/mc/Run3Summer22NanoAODv12/GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2` | `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/ggH2022pico` |

Discovery requires the EOS mounts on lxplus. Condor stages input files over
XRootD into worker scratch. UUID filenames get a local parser-compatible
dataset/year alias; the original files on EOS are not renamed.

Only the final baseline pico is transferred to EOS. NanoAODs and raw picos
are temporary scratch files. stdout/stderr return to AFS at job completion,
not by streaming (which CERN batch does not permit). Modern HTCondor also
returns these diagnostic logs on failure with `ON_SUCCESS`; only a successful
baseline pico is sent to EOS.

## Baseline and normalization

The [nano2pico README](https://github.com/richstu/nano2pico/blob/htozgamma_sequoia_v1/README.md#higgs-to-z-gamma-related-variables)
and working DY workflow use:

```text
use_event && (zg_cutBitMap==3582 || zg_cutBitMap==3583 ||
              zg_cutBitMap==3070 || zg_cutBitMap==3071)
```

Both states of the data-blinding bit are accepted, retaining the Higgs mass
region for MC. No truth matching, category/BDT selection, merging or branch
slimming is added. All existing pico branches and weights are retained.
Conversion uses `--nent -1`, without `--skim`.

Before submission, files are sorted and `FILE_LIMIT`, if supplied, selects
the first N inputs. Normalization reads the Runs metadata of exactly those
selected files, separately for HJ and ggH. It sums the signed `genEventSumw`
and the nine `LHEScaleSumw` entries in the same schema as upstream
`find_normalization.py`. These sums cover all generated events in the selected
NanoAOD files BEFORE baseline cuts, not only the surviving pico events.
This preserves the baseline selection efficiency. With no limit, all
discovered inputs are selected and normalized together.

`inputs.txt` and `jobs.tsv` describe the same selected input set. Unselected
files are not opened for normalization. Missing, nonfinite, zero or non-nine
metadata stops preparation. Workers use these common selected-input sums via
`--norm`; they only adapt the directory key required by process_nano.

A limited submission is now an independently normalized subset, not a partial
production carrying full-directory normalization. Do not combine picos from
separately normalized subsets as one yield sample; select the complete desired
input set in one campaign so all its jobs share the same normalization. Old
prepared campaigns and picos retain their old normalization and must be
regenerated to use this behavior.

For yield/systematics analysis, review upstream's Higgs NNLO reweighting and
LHE weight conventions for MiNNLO; this workflow preserves upstream behavior.

## Setup once

Run on normal lxplus, not an EosSubmit-only schedd. The existing checkout is
`/afs/cern.ch/user/j/junhyuk/nano2pico_sequoia_v1`; use its intended
`htozgamma_sequoia_v1` sources. The setup uses your existing CMSSW project
and SCons installation.

```bash
cd ~/HJ_gen_script/tools/baseline
source setup_runtime.sh
cd "$NANO2PICO_DIR"
scons -j2
cd ~/HJ_gen_script/tools/baseline
```

Skip the build if the current sources are already built in this runtime.
Like DY, the worker calls the existing `run/process_nano.exe` directly and
runs from the checkout to find correction data. Keep this checkout and the
baseline helper files unchanged while jobs are queued/running. If the
launcher cannot find a build on a worker, its stderr will report that error;
this workflow does not package or rebuild the converter per job.

## Proxy

Use an AFS path, never a login node's `/tmp`:

```bash
mkdir -p "$HOME/tmp"
chmod 700 "$HOME/tmp"
export X509_USER_PROXY="$HOME/tmp/x509up"
voms-proxy-init --voms cms --valid 168:00 --out "$X509_USER_PROXY"
```

Export this same path in every submitting shell, and renew it there before
expiry. Preparation requires at least 30 minutes remaining; for production,
ensure the lifetime covers queue time and execution.

Do not force `XrdSecPROTOCOL=gsi`: EOS redirects may negotiate other
authentication protocols. If an old setup exported it, run
`unset XrdSecPROTOCOL` before preparation/submission.

## Prepare, dry-run, submit

Syntax: `bash submit_baseline.sh YEAR HJ|ggH|both CAMPAIGN [FILE_LIMIT] [--dry-run|--submit]`.
Without a flag, it only prepares. Each invocation needs a new campaign name.

```bash
cd ~/HJ_gen_script/tools/baseline

# One-file pilot for each sample; normalize only that selected input file.
bash submit_baseline.sh 2022 both pilot_2022 1 --dry-run

# Submit those already-prepared jobs, without repeating normalization.
condor_submit runs/pilot_2022/HJ/baseline.sub
condor_submit runs/pilot_2022/ggH/baseline.sub
condor_q
```

For an HJ sample consisting of 100 files, with an empty HJ pico output directory:

```bash
bash submit_baseline.sh 2022 HJ selected100_2022 100 --dry-run
condor_submit runs/selected100_2022/HJ/baseline.sub
```

This uses the same 100 selected inputs for normalization and conversion.

A dry-run validates submit syntax, not the converter, EOS I/O or physics.
Check the two real pilot jobs before proceeding.

If CERN rejects `stream_output` / `stream_error` in an old prepared submit
file, remove those lines from both that file and the template. Reuse the
existing normalization JSON and job list; do not repeat preparation.

```bash
# Full samples: one conversion+baseline job for every discovered file.
bash submit_baseline.sh 2022 both production_2022 --submit
```

Preparation refuses an output directory that already contains ROOT files.
Archive/remove pilot outputs before full production, or edit the output
directory constants to fresh paths. Do not submit overlapping campaigns
for the same sample/output path. No existing pico is skipped or overwritten.

For a failed job, inspect `runs/CAMPAIGN/SAMPLE/logs/baseline_CLUSTER.PROC.err`
and `condor_q -hold -af ClusterId ProcId HoldReasonCode HoldReason`.
After fixing a transient issue, release the held job using its exact ID;
do not rerun preparation into an existing campaign. To retry under changed
code, remove the old jobs and prepare a fresh campaign/output set.

## Offline checks

```bash
for script in tools/baseline/*.sh; do bash -n "$script" || exit; done
python3 tools/baseline/tests/test_baseline.py
```

ROOT-dependent tests skip if PyROOT is unavailable. These checks do not submit
jobs or access CERN services; an lxplus pilot is still required.
