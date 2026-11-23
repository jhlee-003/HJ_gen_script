# HJ analysis tools

- [Baseline workflow](baseline/README.md): convert NanoAOD to baseline-selected picos.
- Pico plotting: the three standalone scripts in `tools/plotter` compare BDT inputs, kinematics and NNLO/NLO ratios.

## BDT-input comparison

With CMSSW/PyROOT active, run the single consolidated plotter:

```bash
python3 tools/plotter/plot_ggF_BDT_var.py \
  --sample1-name "HJ (central gridpack)" \
  --sample1-location /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico \
  --sample1-type HJ \
  --sample2-name "HJ (private gridpack)" \
  --sample2-location /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022PrvGridpico \
  --sample2-type HJ \
  --output-prefix HJ_gridpack_2022 \
  --max-depth 1
```

Samples 1 and 2 require all six sample arguments. To overlay an optional third
sample, supply `--sample3-name`, `--sample3-location`, and `--sample3-type`
together. For example, append:

```bash
  --sample3-name "ggH (central)" \
  --sample3-location /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/ggH2022pico \
  --sample3-type ggH
```

Add a continuation backslash to the preceding command line when appending
these flags. Omitting all three keeps the existing two-sample behavior.
Each `ggH` sample uses `weight`; each
`HJ` sample uses `weight / w_nnlo`, regardless of its position. Sample 1 is red,
sample 2 blue, and optional sample 3 green. Sample 1 is always the reference:
ratios are sample 2 / sample 1 and sample 3 / sample 1. With three samples,
both ratio curves use their matching overlay colors and legend/count rows
expand to three. The same 19 variables, ranges,
binning and three PNG layouts (8 + 8 + 3) are retained. No new selections.

The two optional arguments are `--output-prefix` (default derived from sample
names, including sample 3 when present) and `--max-depth` (default 1;
use 2 to include one subdirectory level, applied to every sample).
The example outputs `plots/HJ_gridpack_2022_{1,2,3}.png`. Existing PNGs are
protected from overwrite. Each histogram is normalized independently to unit
summed plotting weight; negative weights are retained and overflow is folded
into the last visible bin. Run the scripts with `--help` for their arguments.

The old separate BDT variants and their year-only CLI are retired.

## Kinematics comparison

Use the same sample arguments with the standalone kinematics plotter:

```bash
python3 tools/plotter/plot_ggF_kinematics.py \
  --sample1-name "HJ (central gridpack)" \
  --sample1-location /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico \
  --sample1-type HJ \
  --sample2-name "HJ (private gridpack)" \
  --sample2-location /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022PrvGridpico \
  --sample2-type HJ \
  --output-prefix HJ_kinematics_2022 \
  --max-depth 1
```

One output: `plots/HJ_kinematics_2022.png`. The existing six variables, binning,
ranges, 3 x 2 layout and 91.2/125 GeV mass reference lines are retained.
Each `ggH` sample uses `weight`; each `HJ` sample uses `weight / w_nnlo`.
Sample names become legends. The same optional `--sample3-name`,
`--sample3-location`, `--sample3-type` group is supported; all three must be
supplied together. Ratios compare each additional sample to sample 1.
The third sample is green and uses `weight / w_nnlo` if declared HJ, or `weight`
if declared ggH. Defaults are
depth 1 and a filename prefix derived from both names. Existing output PNGs
are protected from overwrite. The two old kinematics variants and year-only
CLI are retired. The new script does not depend on any other plotter.

## NNLO/NLO ratio comparison

```bash
python3 tools/plotter/plot_ggF_NNLO_to_NLO_ratio.py \
  --sample1-name "Central gridpack" \
  --sample1-nnlo-location /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico \
  --sample1-nlo-location /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/ggH2022pico \
  --sample2-name "Private gridpack" \
  --sample2-nnlo-location /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022PrvGridpico \
  --sample2-nlo-location /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/ggH2022pico
```

Sample 1's name/NNLO(HJ)/NLO(ggH) locations are required. The three `--sample2-*`
arguments and three `--sample3-*` arguments are optional all-or-none groups.
One PNG overlays up to three named ratios against the embedded Run 3 lookup,
with ratio/lookup below. It uses generator-Higgs pT and the existing ratio
prescription: BOTH HJ and ggH use `weight / w_nnlo`, separately normalized.
This is distinct from BDT/kinematics weighting. The script is standalone and
needs neither of the deleted Higgs-pT/w_nnlo diagnostics.

Output: `plots/ggF_NNLO_to_NLO_ratio_<comparison-names>.png`.
The embedded Run 3 lookup covers 16 generator-Higgs-pT intervals from 0 to
500 GeV. Events without a saved Higgs or with pT >= 500 GeV are reported and
omitted. Use `--help` for the full argument list.
