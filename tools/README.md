# HJ analysis tools

- [Baseline workflow](baseline/README.md): convert NanoAOD to pico and save baseline-selected events.
- [Pico plotter](plotter/README.md): compare ten weighted kinematic shapes between HJ and central ggF picos.

With CMSSW/PyROOT active, run from the repository root:

    python3 tools/plotter/plot_pico_kinematics.py 2022

Output: `plots/HJ_2022_pico_central_vs_private.png`.
The plotter is standalone; the retired NanoAOD plotting scripts are no longer required.
