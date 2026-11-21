# 2022 private-gridpack HJ production: M120, M125 and M130

This is the production directory for all three 2022 Higgs mass points. Each
mass has its own sequential, DY-style producer and configuration set; Condor
and CRAB wrappers are shared. Existing EOS data are not moved or regenerated
by these scripts. Gridpack contents, including the approved `fakevirt 1`
setting, are used as supplied and are not rewritten by the submission setup.

## Directory layout

The `{120,125,130}` notation below represents three separate files:

    2022/
      ProduceHJ_M-{120,125,130}_TuneCP5_2022.sh
      HJ_el8_amd64_gcc11_CMSSW_13_2_5_HJ_MiNNLO_M{120,125,130}_Private.tgz
      condor_convert_wrapper.sh
      crab_convert_wrapper.sh
      PSet.py
      config/
        HJ_M-{120,125,130}_TuneCP5_2022_FullSim.env
        HJ_M-{120,125,130}_TuneCP5_2022_fragment.py
        HJ_M-{120,125,130}_TuneCP5_2022_condor.sub
        HJ_M-{120,125,130}_TuneCP5_2022_crabConfig.py
      tests/check_private_gridpack.py
      logs/M120/
      logs/M125/
      logs/M130/

Place the three actual gridpack archives beside the producers, using exactly
these names. Do not rename the M125 archive to stand in for another mass.
The archive's POWHEG `hmass` and the fragment's `25:m0` must match the mass.

## Mass-specific settings and shared wrappers

Each `FullSim.env` declares `HIGGS_MASS`, `PRODUCER_SCRIPT`, `GRIDPACK_FILENAME`,
`Fragment_filename` and the output names. `BASE_TAG` contains both the mass and
`__PrvGrid`; AOD, MiniAOD and NanoAOD names continue to derive from it.

The Condor wrapper reads the producer name from the selected environment file.
CRAB passes that same name explicitly, and its wrapper verifies the match.
Each input sandbox contains only the matching producer, environment file,
fragment and gridpack. The producers reject mixed-mass configurations.

The producer resolves and exports an absolute worker-local `GRIDPACK_PATH`
before CMSSW setup, checks the archive, and logs its SHA-256. The fragment
checks its mass-specific basename. This works after Condor/CRAB flatten the
input sandbox and when concurrent LHE generation enters per-thread directories.
No developer home path is embedded in the fragment.

## CMSSW pipeline

| Stage | M120 | M125 | M130 | SCRAM_ARCH |
| --- | --- | --- | --- | --- |
| LHE/GEN/SIM | CMSSW_12_4_14_patch2 | CMSSW_12_4_11_patch3 | CMSSW_12_4_14_patch2 | el8_amd64_gcc10 |
| DIGIPREMIX and RECO | CMSSW_12_4_11_patch3 | CMSSW_12_4_11_patch3 | CMSSW_12_4_11_patch3 | el8_amd64_gcc10 |
| MiniAODv4 and NanoAODv12 | CMSSW_13_0_13 | CMSSW_13_0_13 | CMSSW_13_0_13 | el8_amd64_gcc11 |

The LHE/GEN/SIM release difference follows the central mass-point setups.
The gridpacks' own CMSSW_13_2_5/el8_amd64_gcc11 runtime remains separate.
Global tags remain `124X_mcRun3_2022_realistic_v12` through RECO and
`130X_mcRun3_2022_realistic_v5` for MiniAOD/NanoAOD. The beamspot, premix
dataset, HLT:2022v12, era/modifiers and step ordering are unchanged.

The fragments preserve the HJ MiNNLO Pythia matching, CP5 tune and PS weights,
with H -> Z gamma and Z -> ee, mumu or tautau (`mZ > 50 GeV`). Only the Higgs
mass changes. NanoAOD event content and datatier are both `NANOAODSIM`.
FrameworkJobReport handling, inherited proxies, thread arguments and final
output validation remain in place. There are no timestamp-based names.

## Condor and EOS stageout

Each Condor file starts with a safe pilot: **10 events, one job, two threads**.
The existing requests are 5000 MB memory, 8000000 KB disk and the `tomorrow`
job flavour. Tune production sizes after checking a pilot's measured usage.

| Mass | Final NanoAOD destination | Local event/stdout/stderr log paths |
| --- | --- | --- |
| M120 | /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022/M120/ | logs/M120/ |
| M125 | /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022/M125/ | logs/M125/ |
| M130 | /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022/M130/ | logs/M130/ |

The NanoAOD is transferred directly from worker scratch to EOS, with no extra
persistent local copy. Filenames contain the mass, Condor cluster and job number.
Log directories are ignored by Git and must be created in every new checkout.
CERN submit-side transforms may relocate stdout/stderr to EOS; query the actual
`Out`, `Err` and `Iwd` job attributes rather than assuming a returned log path.

On LXPLUS, with the usual CMS tools and a valid proxy:

    cd ~/HJ_gen_script/2022
    mkdir -p logs/M120 logs/M125 logs/M130
    export X509_USER_PROXY="$HOME/tmp/x509up"
    voms-proxy-info --file "$X509_USER_PROXY" --timeleft
    condor_submit -dry-run /tmp/HJ_M120_2022.classad config/HJ_M-120_TuneCP5_2022_condor.sub
    condor_submit -dry-run /tmp/HJ_M130_2022.classad config/HJ_M-130_TuneCP5_2022_condor.sub

After checking the dry runs, submit the two new pilots explicitly:

    condor_submit config/HJ_M-120_TuneCP5_2022_condor.sub
    condor_submit config/HJ_M-130_TuneCP5_2022_condor.sub

M125 has its own equivalent submission file; the existing successful M125
dataset does not need to be regenerated. For additional production, use the
matching mass file and adjust `events_per_job`, `number_of_jobs` and `job_offset`
per mass. Command-line `-append` overrides may also be used. Human-facing job
numbers are `Process + job_offset + 1`. Do not resubmit already completed jobs
as an extension of a campaign. No other mass's gridpack needs to be transferred.

The old central-gridpack sample is now in EOS `2022-test/`, outside the three
production mass directories. No archival data or repository `2022-test/` files
are modified by this setup. Analysis code should select the intended mass
directory explicitly rather than combine all three.

## CRAB compatibility

There is one CRAB configuration per mass, also defaulting to **10 events and
one job**. Input files, request names, dataset tags and producer arguments are
mass-specific; the bootstrap PSet and wrapper remain shared.

CRAB retains the existing `T3_KR_KNU` storage-site stageout. It does **not** use
the CERNBox destinations above. Use Condor for the requested CERNBox workflow.

## Validation without CMSSW

    python3 tests/check_private_gridpack.py

This checks all three configuration sets, release boundaries, generated-job
syntax, physics settings, per-mass output/input names, and mocked worker
dispatch. Temporary fixtures let configuration checks run before all archives
are added. Real archives that are present are opened read-only to check their
POWHEG mass and expected runtime files; no gridpack is executed or extracted.

Once all three payloads are installed, require them explicitly:

    python3 tests/check_private_gridpack.py --check-gridpacks

These checks do not run actual CMSSW, authenticate to EOS, or submit jobs.
Real dry runs and complete CMSSW pilots are still required for execution and
stageout validation. The existing 12_4_11 pixel-reconstruction assertion is not
patched or claimed fixed by this mass-point configuration change.
