# Central versus private NanoAOD kinematic plots

Use a Python environment with PyROOT (for example, CMSSW on LXPLUS). Run from
the repository root:

    python3 tools/plot_nanoaod_kinematics.py /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022
    python3 tools/plot_nanoaod_kinematics_no_selection.py /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022

Each command takes **only the private EOS directory**. An EOS directory URL,
such as `root://eosuser.cern.ch//eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022`,
also works. All private top-level `.root` files are included, without a Condor
cluster filter. This includes any pilot files or duplicate productions.

Both scripts compare against this fixed **2022** central directory:

    /eos/cms/store/mc/Run3Summer22NanoAODv12/GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2

Central file discovery has the same depth as `find DIRECTORY -maxdepth 2`:
files in the directory and its immediate subdirectories, but no deeper files.
Discovery uses the mounted EOS filesystem if available, otherwise bounded
`xrdfs ls`/directory `stat` calls. ROOT reads the files through XRootD; no
merging or copying of the NanoAOD files is needed. Use a valid proxy, for example
`export X509_USER_PROXY="$HOME/tmp/x509up"`, in the CMS/XRootD environment.

Every input must have a NanoAOD `Events` tree and the required branches. A broken
file causes an explicit error instead of silently changing the comparison.
Two ROOT worker threads fill all ten histograms in one event loop per sample.

## Output and image formatting

The commands produce one five-column, two-row PNG each:

- `plots/HJ_2022_central_vs_private_selected.png`
- `plots/HJ_2022_central_vs_private_no_selection.png`

Output paths are relative to this repository, not the working directory.
Rerunning replaces the corresponding PNG; generated PNGs are ignored by Git.
The no-selection entry point imports the shared implementation in the selected
plotter, so their file discovery, weighting, and rendering cannot drift apart.

The layout follows the supplied comparison script: red Central and blue Private
line histograms, weighted error bars, legends with unweighted entry counts, a
CMS Private Work label and 13.6 TeV in each panel, and a Private/Central ratio
pad with a dashed line at one. There is **no global title, event-summary banner,
or bottom selection/explanation text**. Selection details belong in this README
and processing counts are printed to the terminal. No DY truth-matching
requirements or data definitions from the reference script were adopted.

## Ten variables

| Variable | Definition | Selected range | No-selection range |
| --- | --- | --- | --- |
| Higgs transverse recoil | `pT(llgamma)/m(llgamma)` | 0–2 | same |
| Photon pseudorapidity | `eta(gamma)` | −2.5–2.5 | same |
| Nearest photon–lepton separation | `min(DeltaR(gamma,l1), DeltaR(gamma,l2))` | 0–6 | same |
| Farthest photon–lepton separation | `max(DeltaR(gamma,l1), DeltaR(gamma,l2))` | 0–6 | same |
| Jet multiplicity | Number of jets as defined below | 0–8, ≥9 | same |
| Leading-jet transverse momentum | Maximum jet `pT`; requires at least one jet | 0–300 GeV | same |
| Dilepton–photon invariant mass | `m(llgamma)` | 100–180 GeV | 0–300 GeV |
| Dilepton invariant mass | `m(ll)` | 50–120 GeV | 0–300 GeV |
| Photon transverse momentum | `pT(gamma)` | 0–200 GeV | same |
| Dilepton–photon transverse momentum | `pT(llgamma)` | 0–300 GeV | same |

## Candidate definitions

Both samples receive exactly the same treatment within each version. Both
versions use opposite-sign, same-flavor electron or muon pairs, as requested.
Among eligible pair/photon combinations, choose the dilepton mass closest to
91.1876 GeV, breaking equal-mass ties with the highest-pT photon. Each event
contributes at most one candidate. Neither version is generator truth-matched;
tau daughters and other reconstructed objects can contribute. The central
POWHEG-Pythia8 sample and private POWHEG-MiNNLO-Pythia8 sample can have genuinely
different production shapes.

### With the existing selection

- Both leptons have `pT >= 10 GeV`; the leading lepton has `pT >= 20 GeV`.
  Electrons require `abs(eta) < 2.5`; muons require `abs(eta) < 2.4`.
- Require `50 <= m(ll) <= 120 GeV` and `100 <= m(llgamma) <= 180 GeV`.
- The photon has `pT >= 15 GeV`, `abs(eta) < 2.5`, excludes
  `1.4442 < abs(eta) < 1.566`, and has `DeltaR(gamma,l) > 0.4` for both leptons.
  The gap cut uses `Photon_eta`, not a corrected supercluster coordinate.
- Count tight-ID jets with `pT > 30 GeV`, `abs(eta) < 4.7`, and
  `(Jet_jetId & 2) != 0`, separated by `DeltaR > 0.4` from both leptons and photon.
- No lepton/photon ID, isolation, trigger cuts, or efficiency/energy-scale
  corrections are applied, matching the previous plotter.

### Without object selection

No lepton/photon pT or eta thresholds, ECAL-gap veto, mass windows, photon–lepton
separation cut, jet pT/eta thresholds, jet ID, or jet overlap cleaning are applied.
Jets are counted directly from the NanoAOD jet collection; array ordering is not
assumed. This does not undo object thresholds already used to construct NanoAOD.

OS ee/mumu pairing and the candidate-choice rule above are retained. The
llgamma distributions necessarily require two such leptons and one photon;
the leading-jet distribution additionally requires a jet. Objects must have
finite, physical four-vectors and the llgamma mass must be positive so that
the observables are mathematically defined. These are validity checks, not
analysis object cuts. Thus “no selection” does not mean every NanoAOD event
can supply every plotted observable.

## Bin weighting and ratios

The full signed `genWeight` is retained, including negative weights. Each
sample's histogram is normalized independently to unit summed weight, with
sum-of-squared-weight statistical errors. These compare **shapes**, not yields
scaled to cross section or luminosity. The leading-jet histogram has its own
normalization from the events with jets. Empty jet histograms remain empty.
A non-empty histogram with non-positive total weight raises an error.

Underflow/overflow contents and their squared errors are folded into the first
and last bins before normalization, as in the previous plotter. Edge-bin spikes
can therefore include events outside the displayed ranges, especially without
selection; the range is not a cut. The last jet bin means at least nine jets.

The ratio is the normalized Private shape divided by the normalized Central
shape, with independent weighted-bin error propagation. Zero Central bins
are omitted, not drawn as zero ratios. Signed bins are retained. The reference
ratio range is 0–2 and expands if needed to show outlying/signed ratio points.
Errors are diagnostic per-bin statistical errors; they do not include
normalization-induced bin correlations or systematic uncertainties.

## Local validation

With PyROOT available, run:

    python3 tools/tests/check_plotter.py

This creates temporary synthetic NanoAOD files and checks both modes, multi-file
reading, maxdepth-2 discovery, negative weights/errors, candidate ranking,
jet cleaning, overflow handling, ratios, malformed inputs, and two PNG previews.
It does not access production EOS data. Previews are saved separately as
`plots/HJ_synthetic_test_central_vs_private_selected.png` and
`plots/HJ_synthetic_test_central_vs_private_no_selection.png`.
