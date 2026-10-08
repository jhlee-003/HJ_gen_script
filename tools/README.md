# HJ analysis tools

- [Baseline workflow](baseline/README.md): convert NanoAOD to pico and save baseline-selected events.
- [Pico plotter](plotter/README.md): compare 19 weighted ggF BDT-input shapes between HJ and central ggF picos.

With CMSSW/PyROOT active, run from the repository root:

    python3 tools/plotter/plot_ggF_BDT_var.py 2022

Outputs: `plots/ggF_BDT_var_2022_1.png`, `_2.png`, and `_3.png` (8 + 8 + 3 panels).
The plotter is standalone; the retired NanoAOD plotting scripts are no longer required.

For a six-variable kinematics comparison (photon, dilepton and llgamma pT;
dilepton and llgamma mass; primary-vertex count), run:

    python3 tools/plotter/plot_ggF_kinematics.py 2022

Output: `plots/ggF_kinematics_2022.png`, one 3 x 2 image using the older, taller
panel proportions. Both scripts use the same baseline pico directories.
