# Baseline-pico ggF BDT-input comparison

With CMSSW/PyROOT active, run from the repository root:

    python3 tools/plotter/plot_pico_kinematics.py 2022

Only the year is required. Currently 2022 uses all top-level `.root` files in:

- HJ: `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico`
- ggH: `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/ggH2022pico`

Files are read through XRootD without merging, copying, or filtering by Condor
cluster. Keep duplicate productions and unrelated ROOT files out of these
sample directories. Add future paths to `SAMPLE_DIRECTORIES`.

## Three output images

The 19 variables follow Figure 48 in the supplied ggF BDT-input reference PDF,
read left-to-right and then top-to-bottom:

1. `plots/HJ_2022_pico_central_vs_private_1.png`: eight panels, pT(llgamma)/m(llgamma) through phi, arranged 4 x 2.
2. `plots/HJ_2022_pico_central_vs_private_2.png`: eight panels, eta(l1) through pT(j1), arranged 4 x 2.
3. `plots/HJ_2022_pico_central_vs_private_3.png`: three panels, DeltaPhi(Zgamma,j1), system balance, and photon Zeppenfeld, arranged 3 x 1.

Each panel has an approximately 3:2 width/height ratio. Existing styling is
retained: CMS Private Work / 13.6 TeV, red ggH_qme (central), blue HJ (private),
top-aligned y titles, weighted errors, unit-integral shapes, folded edge bins,
and 0-2 Private/Central ratios. There is no legend heading, global explanation,
or bottom note. The event-count legend is shifted slightly left. The old
m(llgamma), pT(llgamma), and jet-count panels are not part of the PDF's 19-panel
set, so there is no mass panel or 125 GeV marker in these images.
Older PNGs are not deleted automatically.

## Reference axes and binning

Ranges are display limits, not new event cuts. Values outside them are folded
into the first/last bins with their weights and errors retained.

| Image | Variable | Range | Bins | Bin width |
| --- | --- | --- | --- | --- |
| 1 | pT(llgamma)/m(llgamma) | 0 to 2.5 | 40 | 0.0625 |
| 1 | Photon IDMVA | 0 to 1 | 45 | 0.0222222 |
| 1 | Photon sigmaE/E | 0.01 to 0.25 | 40 | 0.006 |
| 1 | min DeltaR(gamma,l) | 0 to 3.5 | 35 | 0.1 |
| 1 | max DeltaR(gamma,l) | 0 to 6 | 60 | 0.1 |
| 1 | cosTheta | -1 to 1 | 40 | 0.05 |
| 1 | costheta | -1 to 1 | 40 | 0.05 |
| 1 | phi | -3.2 to 3.2 | 40 | 0.16 |
| 2 | eta(l1) | -2.6 to 2.6 | 40 | 0.13 |
| 2 | eta(l2) | -2.6 to 2.6 | 40 | 0.13 |
| 2 | eta(gamma) | -2.6 to 2.6 | 40 | 0.13 |
| 2 | DeltaPhi(HmissT,gamma) | 0 to 3.15 | 40 | 0.07875 |
| 2 | DeltaR(gamma,j1) | 0.4 to 6 | 40 | 0.14 |
| 2 | eta(j1) | -5 to 5 | 50 | 0.2 |
| 2 | m(j1) [GeV] | 0 to 40 | 50 | 0.8 GeV |
| 2 | pT(j1) [GeV] | 30 to 150 | 40 | 3 GeV |
| 3 | DeltaPhi(Zgamma,j1) | 0 to 3.15 | 40 | 0.07875 |
| 3 | System balance | 0 to 1 | 40 | 0.025 |
| 3 | Photon Zeppenfeld | 0 to 6 | 40 | 0.15 |

## Pico definitions

Use `tree`, not NanoAOD `Events`. No object ID, isolation, pT/eta, trigger,
mass-window, truth matching, `use_event`, bitmap, or ggF-category cuts are
reapplied. Candidate zero is nominal: photon index `llphoton_iph[0]`, pair
index `llphoton_ill[0]`, then `ll_i1[ill]`, `ll_i2[ill]`, and `ll_lepid[ill]`.
Every index is checked. Undefined/sentinel/non-finite observables are omitted
only from the affected panel; no-jet events remain in the other panels.

- Recoil: `llphoton_pt[0]/llphoton_m[0]`.
- Photon IDMVA: `photon_idmva[iph]`.
- Photon resolution: `photon_energyErr[iph]/(photon_pt[iph]*cosh(photon_eta[iph]))`, i.e. sigmaE/E for a massless photon as labeled in the PDF. This is not sigmaE/pT. The current upstream `ggF_input_plots` code uses pT in the denominator despite the E label; this distinction matters for literal reproduction of a particular upstream curve.
- Lepton separations: wrapped DeltaR to the two indexed candidate leptons, not other leptons in the event.
- Angles: stored `llphoton_cosTheta[0]`, `llphoton_costheta[0]`, and `llphoton_psi[0]` (the plotted phi). No sign flips or reconstruction.
- Lepton eta: sort the candidate's two leptons by stored `el_pt` or `mu_pt`; l1 is higher-pT, l2 lower-pT. Eta remains floating-point.
- Photon eta and missing-HT DeltaPhi: `photon_eta[iph]`, `photon_mht_dphi[iph]`.
- Leading jet: highest-pT stored `jet_isgood` jet, requiring a defined jet with `njet > 0`. Eta, phi, mass, and pT refer to that same index. Stored good-jet flags are the producer's definition, not new cuts.
- Photon-jet DeltaR and Zgamma-jet DeltaPhi: calculated from that jet and the nominal candidate photon/`llphoton_phi[0]`, with wrapped phi differences.
- System balance: stored `mht/ht`, the magnitude of the vector pT sum divided by the scalar pT sum over nano2pico's signal leptons/photons and good jets. Do not substitute `llphoton_dijet_balance`, which groups different objects.
- Photon Zeppenfeld in the PDF's single-jet definition: `abs(photon_eta[iph]-jet_eta[j1])`. Do not use the two-jet mean-eta definition in `photon_zeppenfeld` for multi-jet events.

Definitions are checked against the sequoia
[pico schema](https://github.com/richstu/nano2pico/blob/htozgamma_sequoia_v1/txt/variables/pico)
and [zgamma producer](https://github.com/richstu/nano2pico/blob/htozgamma_sequoia_v1/src/zgamma_producer.cpp).
Reference bins are also consistent with the applicable
[draw_pico histogram definitions](https://github.com/richstu/draw_pico/blob/master/src/zgamma/categorization_utilities.cpp).

## Weights and interpretation

Use the signed stored pico `weight`, which already contains normalization and
corrections. Do not multiply by `genWeight`, `w_lumi`, or correction factors
again. Each histogram is independently normalized to unit summed weight;
negative weights and sum-of-squared-weight errors are retained. Ratios compare
private/central shapes, not signal/background as in the PDF.

`N=` is the unweighted number of events with a defined value in that panel,
including events folded into edge bins. It is not a file/job count or an
expected yield. The same counts and undefined-value counts are printed.

Matching the variables and axes does not guarantee literal reproduction of
the reference red curve: the PDF combines Run 2/Run 3 and describes ggF signal
stacked on VBF, with signal normalized to background. This plotter compares
only your baseline-selected 2022 HJ and central ggH samples, without adding
those samples, additional selections, or background normalization.
For yield studies, subset normalization must also be correct; a common
normalization constant cancels only in these shape plots.

Every file's required branches are checked. Zero-entry files are supported
inside a nonempty sample. Broken files, entirely empty samples, or non-finite
weights cause explicit errors. All 19 histograms and counts are booked before
execution: one event loop per sample, not one per output image.

Synthetic validation (no EOS access):

    python3 tools/plotter/tests/check_pico_plotter.py

Tests cover candidate/flavor indices, float eta and pT ordering, photon energy
resolution, stored angles, leading good jets, balance, signed weights/errors,
PDF ranges/binning, empty/malformed inputs, year routing, and all three rendered
layouts. Previews are `/tmp/HJ_pico_synthetic_check_1.png` through `_3.png`.
