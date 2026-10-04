# NanoAOD kinematic plots

From an environment where `python3 -c 'import ROOT'` works (such as a CMSSW
environment on LXPLUS), run from the repository root:

    python3 tools/plot_nanoaod_kinematics.py /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022

The EOS directory is the only argument. A directory URL is also accepted:

    python3 tools/plot_nanoaod_kinematics.py root://eosuser.cern.ch//eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022

All top-level `.root` files are read through XRootD with ROOT RDataFrame; there
is no Condor-cluster filter or recursive search. Pilot files and duplicate
productions in the same directory are included too. Every file must contain a
NanoAOD `Events` tree and the required branches; invalid files cause an error
rather than silently disappearing from the sample. ROOT uses two worker threads
and fills all six histograms in one event loop without merging the inputs.

The output is one three-column, two-row PNG:
`plots/HJ_2022_kinematics.png` for a directory named `2022`. Its location is
relative to this repository, regardless of the shell's working directory.
Running again for the same directory name replaces that PNG. Generated PNGs are
ignored by Git.

## Variables

These six variables come from the supplied H-to-Z-gamma variable table:

| Panel | Definition and ggF motivation |
| --- | --- |
| Higgs transverse recoil | `pT(llgamma)/m(llgamma)`, sensitive to the recoil of the Higgs against QCD radiation. |
| Photon pseudorapidity | `eta(gamma)`, describing the decay photon's angular distribution and detector acceptance. |
| Nearest photon-lepton separation | `min(DeltaR(gamma,l1), DeltaR(gamma,l2))`, describing the decay geometry after the separation cut. |
| Farthest photon-lepton separation | `max(DeltaR(gamma,l1), DeltaR(gamma,l2))`, complementing the nearest separation. |
| Jet multiplicity | Number of cleaned tight-ID jets with `pT > 30 GeV` and `abs(eta) < 4.7`, describing accompanying QCD radiation. |
| Leading-jet transverse momentum | Highest `pT` among those jets, describing the hard recoil; filled only for events with at least one such jet. |

No jet is required for the first five distributions. A sample made with an HJ
MiNNLO gridpack can still have zero reconstructed jets above the selected
threshold. Photon ID MVA and energy resolution were not selected because they
measure reconstruction quality rather than production kinematics.

## Reconstructed candidate and weights

This is a diagnostic reconstructed-level plotter with a loose kinematic
selection, not a complete H-to-Z-gamma analysis selection or a truth-matched
measurement. It does not require lepton/photon IDs, isolation, triggers, or apply
efficiency/energy-scale corrections. Tau decays can enter if their reconstructed
electron/muon daughters satisfy the selection; there is no generator ancestry
requirement.

- Use opposite-sign, same-flavor electron or muon pairs. Both leptons need
  `pT >= 10 GeV`, with the leading lepton `pT >= 20 GeV`; electrons require
  `abs(eta) < 2.5` and muons `abs(eta) < 2.4`.
- Require `50 <= m(ll) <= 120 GeV` and `100 <= m(llgamma) <= 180 GeV`.
- Require a photon with `pT >= 15 GeV`, `abs(eta) < 2.5`, outside
  `1.4442 < abs(eta) < 1.566`, and `DeltaR(gamma,l) > 0.4` for each lepton.
  The gap cut uses `Photon_eta`, not a corrected supercluster coordinate.
- Among eligible pair/photon combinations, choose the dilepton mass closest to
  91.1876 GeV. For the same pair/mass, choose the highest-pT eligible photon.
  Each event contributes at most one candidate.
- Jets require the tight-ID bit `(Jet_jetId & 2) != 0` and separation
  `DeltaR(j,l1) > 0.4`, `DeltaR(j,l2) > 0.4`, `DeltaR(j,gamma) > 0.4`.
  Their order in the NanoAOD arrays is not assumed.

Histograms use the full signed `genWeight`, including negative weights. Each
panel is normalized to unit summed weight with weighted statistical error bars;
these are shape plots, not yields scaled to a luminosity. Underflow and overflow
are folded into the first/last bins before normalization (the last jet-count bin
means at least nine jets). The leading-jet panel has its own normalization using
the events with at least one jet. Empty jet samples leave that panel empty.

## Local validation

With PyROOT available, run `python3 tools/tests/check_plotter.py`. It builds two
small synthetic NanoAOD files with different cluster names and checks candidate
selection, negative weights, jet cleaning, overflow handling, and PNG rendering.
Its preview is labelled synthetic and saved separately as
`plots/HJ_synthetic_test_kinematics.png`; it does not read EOS production data.
