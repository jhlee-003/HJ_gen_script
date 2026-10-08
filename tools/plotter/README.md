# Baseline-pico kinematic comparison

## Run the plotter

For the already baseline-selected HJ and central ggF picos, run:

    python3 tools/plotter/plot_pico_kinematics.py 2022

This also takes **only the year**. Currently only 2022 is configured, with:

- HJ: `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico`
- ggH: `/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/ggH2022pico`

It reads every top-level `.root` file in each directory through XRootD,
irrespective of Condor cluster or campaign, without merging/copying inputs.
Keep duplicate productions and unrelated ROOT files out of these directories.
Add future year paths to `SAMPLE_DIRECTORIES` in the pico plotter.

Output: `plots/HJ_2022_pico_central_vs_private.png`. The standalone
five-column/two-row PNG keeps the same ten variables and appearance. Pico-specific
display ranges are **100-150 GeV for m(llgamma)** and **0-200 GeV for
pT(llgamma)**; bin counts and the other eight ranges are unchanged.
These are display ranges, not additional event cuts; overflow folding and
unit-integral normalization are retained. It keeps the red/blue labels,
top-aligned y titles, process legends in the upper-left corner, right-aligned
`N=` event counts in the upper-right corner, weighted error bars, overflow folding,
and 0-2 Private/Central ratios. Both legends are inset from the frame, and the
m(llgamma) overlay has a black dashed vertical reference line at 125 GeV.
The sample-heading legend is removed. There is no global title
or bottom explanation. No draw_pico installation/build or other plotting script is required.

The input is a pico **`tree`**, not a NanoAOD `Events` tree. No object ID,
isolation, pT/eta, trigger, mass, truth-matching, `use_event`, or bitmap cuts
are reapplied. The supplied `plot_vars_central_vs_private.py` was used only
as a reference for pico array/index handling, not its selections, unweighted
histograms, variables, file discovery, or formatting.

### Pico variable mapping

Use nominal llgamma candidate index 0, as used by the baseline producer.
`iph = llphoton_iph[0]`, `ill = llphoton_ill[0]`; the lepton indices are
`ll_i1[ill]` and `ll_i2[ill]`, with flavor from `ll_lepid[ill]` (11 or 13).
All indices are checked before dereferencing.

| Existing panel | Pico source |
| --- | --- |
| Higgs transverse recoil | `llphoton_pt[0] / llphoton_m[0]` |
| Photon eta | `photon_eta[iph]`, not necessarily photon index 0 |
| Nearest photon-lepton separation | Minimum DeltaR to the two indexed candidate leptons |
| Farthest photon-lepton separation | Maximum DeltaR to those same two leptons |
| Jet multiplicity | Stored scalar `njet` |
| Leading-jet pT | Maximum `jet_pt` among stored `jet_isgood` flags, only when `njet > 0` |
| Dilepton-photon mass | `llphoton_m[0]` (nominal, not the refitted mass) |
| Z production angle | Stored `llphoton_cosTheta[0]` |
| Lepton polar angle | Stored `llphoton_costheta[0]` (lowercase theta) |
| Dilepton-photon pT | `llphoton_pt[0]` |

DeltaR uses wrapped delta-phi and the indexed `el_eta/el_phi` or
`mu_eta/mu_phi` arrays. Do not substitute `photon_drmin`: that branch
can include leptons outside the chosen pair. Likewise, the pico jet
collection also contains non-good jets; using the stored flags is not
a new object selection but follows the producer's nominal jet definition.
These definitions are checked against the sequoia branch's
[pico schema](https://github.com/richstu/nano2pico/blob/htozgamma_sequoia_v1/txt/variables/pico),
[candidate producer](https://github.com/richstu/nano2pico/blob/htozgamma_sequoia_v1/src/zgamma_producer.cpp),
and [jet producer](https://github.com/richstu/nano2pico/blob/htozgamma_sequoia_v1/src/jetmet_producer.cpp).

**Angle conventions differ from the old NanoAOD helper.** Here both samples
use nano2pico's stored `CalculateAngles` results unchanged: cosTheta uses
the producer's beam-derived incoming-parton axis; costheta uses the
**positive** lepton relative to the photon in the Z rest frame. The old
NanoAOD plotter uses the laboratory H flight direction and the negative
lepton. No sign flips or angle reconstruction are added in the pico plotter.
An angle outside its physical domain (except tiny float round-off),
a sentinel, a missing index, or a non-finite observable is omitted only
from its own panel. No-jet events still enter the other nine panels.
Counts of undefined/no-object values are printed for every panel.

### Pico weights and checks

Use the signed, stored pico **`weight`** branch, which already includes
luminosity/generator normalization and the correction factors produced
by nano2pico. Do not multiply by `genWeight`, `w_lumi`, or those correction
factors again. Each histogram is independently normalized to unit summed
weight, not divided by its unweighted event count. The upper-right `N=`
counts the unweighted events contributing a defined observable to that panel,
including events folded into the edge bins. For example, the leading-jet
panel excludes events without a good jet, while other panels can retain them.
It is not a file/job count or an expected yield. These counts are also printed
per panel in the terminal.
Negative weights and sum-of-squared-weight errors are retained.

This follows the weighted overlay/shape, overflow, and ratio concepts in
[draw_pico](https://github.com/richstu/draw_pico#plot-options-explanation),
while keeping the existing PyROOT renderer and unit-integral convention.
These are **shape comparisons**, so a common overall normalization cancels;
the plots do not fix subset yield normalization. The baseline submitter now
normalizes exactly the selected input files. Older pilot outputs normalized
to the full input directory must be regenerated before use for subset yields.

Every file's required branches are checked. Zero-entry files are supported
inside a non-empty sample. Broken files, entirely empty samples, and
non-finite stored weights cause explicit errors rather than silent skips.
Two ROOT threads fill all ten histograms in one event loop per sample.

Synthetic validation (no EOS access):

    python3 tools/plotter/tests/check_pico_plotter.py

It tests candidate/flavor indices, the existing good-jet flags, signed
weights/errors, empty or malformed inputs, angle sentinels, year routing,
the one-loop requirement, and actual rendered canvas formatting.
The synthetic preview is written to `/tmp/HJ_pico_synthetic_check.png`.
