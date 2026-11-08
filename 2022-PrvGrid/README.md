# 2022 production with the private HJ gridpack

This directory is independent of the original `2022/` campaign. The CMSSW
stages, releases, global tags, Pythia settings, event numbering and proxy
handling are unchanged. The archive is used as supplied, including `fakevirt 1`.

## Gridpack input

Keep this archive in the same directory as the producer and wrappers:

    HJ_el8_amd64_gcc11_CMSSW_13_2_5_HJ_MiNNLO_M125_Private.tgz

`config/HJ_M-125_TuneCP5_2022_FullSim.env` declares `GRIDPACK_FILENAME`.
The producer checks it, logs its SHA-256, and exports its absolute worker-local
`GRIDPACK_PATH`. The fragment uses that path instead of the central CVMFS
gridpack. Do not hard-code `/home/jhlee/...` or an LXPLUS home path: concurrent
LHE generation runs in per-thread directories on remote workers.

Both Condor and CRAB input lists include the archive. Transfer it alongside
the scripts when moving this directory to LXPLUS; a checkout of the old
`2022/` campaign alone is not sufficient. No gridpack contents are modified
or extracted by the submission setup itself.

## Condor output

Submit from this directory. The existing pilot remains **10 events, one job,
two threads**. Output names contain `__PrvGrid` to distinguish the campaign.
Condor transfers the final NanoAOD directly to:

    root://eosuser.cern.ch//eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022-PrvGrid/

No extra persistent local NanoAOD copy is configured; the worker scratch
output is the transfer source. Event logs retain their existing configuration.
There are no timestamp-based names.

After setting up the usual CMS proxy and Condor environment on LXPLUS:

    cd ~/HJ_gen_script/2022-PrvGrid
    mkdir -p logs
    export X509_USER_PROXY="$HOME/tmp/x509up"
    voms-proxy-info --file "$X509_USER_PROXY" --timeleft
    condor_submit -dry-run /tmp/HJ_2022_PrvGrid.classad config/HJ_M-125_TuneCP5_2022_condor.sub

After a successful dry run, submit explicitly:

    condor_submit config/HJ_M-125_TuneCP5_2022_condor.sub

## CRAB compatibility

The CRAB input sandbox and dataset tag also identify the private gridpack.
CRAB **still uses its existing `T3_KR_KNU` storage-site stageout**, not the
CERNBox destination above. Use the Condor route for the requested CERNBox
output. The existing CRAB event/job counts were not changed to pilot settings.

## Read-only checks

Without CMSSW, run:

    python3 tests/check_private_gridpack.py

These checks do not execute the gridpack, run simulation or submit jobs.
A real CMSSW/Condor pilot is still required to validate execution and stageout.
