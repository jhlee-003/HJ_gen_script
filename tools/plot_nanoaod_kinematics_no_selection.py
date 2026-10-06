#!/usr/bin/env python3
"""Central/private comparison without object cuts; retain OS ee/mumu pairing.

Usage: python3 tools/plot_nanoaod_kinematics_no_selection.py PRIVATE_EOS_DIRECTORY
Output: plots/HJ_<directory-name>_central_vs_private_no_selection.png.
The shared implementation keeps binning/weighting/formatting consistent.
"""

from plot_nanoaod_kinematics import main


if __name__ == "__main__":
    main(apply_selection=False)
