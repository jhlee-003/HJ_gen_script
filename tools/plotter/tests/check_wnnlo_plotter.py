#!/usr/bin/env python3
"""Synthetic unit-weight correction-factor comparison; no EOS access."""

from array import array
from contextlib import redirect_stdout
import io
import math
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_ggF_wnnlo as plotter

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT(2)
ROOT.TH1.SetDefaultSumw2(True)


def write_file(path, values, branch_name="w_nnlo", tree_name="tree", vector=False,
               double=False, include_weight=True):
    source = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree(tree_name, "Synthetic correction factors")
    buffer = ROOT.std.vector("float")() if vector else array("d" if double else "f", [0.])
    if vector:
        tree.Branch(branch_name, buffer)
    else:
        tree.Branch(branch_name, buffer, branch_name + ("/D" if double else "/F"))
    # An invalid stored weight proves it is NOT used to weight the histogram.
    weight = array("f", [float("nan")])
    if include_weight:
        tree.Branch("weight", weight, "weight/F")
    for value in values:
        if vector:
            buffer.clear()
            buffer.push_back(value)
        else:
            buffer[0] = value
        tree.Fill()
    tree.Write()
    source.Close()


def expect_error(phrase, function, *args):
    try:
        function(*args)
    except ValueError as error:
        assert phrase in str(error), str(error)
    else:
        raise AssertionError("Expected " + phrase)


def main():
    assert plotter.parse_arguments(["2022"])[1].year == "2022"
    run3_groups = plotter.protected_groups((*plotter.lookup_weights("2022"), 1.))
    run3_grid = plotter.optimal_binning([{"min": 1., "max": 1.}], "2022")
    assert run3_grid[0] == 84, run3_grid
    assert plotter.separates_values(run3_grid, run3_groups,
                                   run3_groups[0][0], run3_groups[-1][1])
    assert plotter.lookup_weights("2016APV") == plotter.LOOKUP_WEIGHTS["run2"]
    run2_groups = plotter.protected_groups((*plotter.lookup_weights("2018"), 1.))
    run2_grid = plotter.optimal_binning([{"min": 1., "max": 1.}], "2018")
    assert plotter.separates_values(run2_grid, run2_groups,
                                   run2_groups[0][0], run2_groups[-1][1])
    # An independent lower-bound example: 0,.01,.02 require two boundaries
    # inside (0,.02], so spacing < .02. With an endpoint at 2, <=101 bins
    # cannot cover everything; 102 is feasible. A sampled-width heuristic
    # or anchoring the lower edge would miss this optimum.
    hard_groups = plotter.protected_groups((0., .01, .02, 1., 2.))
    assert plotter.minimum_uniform_grid(hard_groups, 0., 2.)[0] == 102
    simple_groups = plotter.protected_groups((0., 1., 2.))
    assert plotter.minimum_uniform_grid(simple_groups, 0., 2.)[0] == 3
    with tempfile.TemporaryDirectory(prefix="ggf_wnnlo_") as tmp:
        directory = Path(tmp)
        central, private, empty = [directory / name for name in ("central.root", "private.root", "empty.root")]
        write_file(central, [0.8, 1., 1.2])
        write_file(private, [0.8, 0.8, 1.2, 1.5])
        write_file(empty, [])
        files = {"Central": [str(central), str(empty)], "Private": [str(private)]}
        frames, summaries = {}, {}
        for sample, paths in files.items():
            frames[sample], summaries[sample] = plotter.inspect_sample(ROOT, paths, sample)
            assert plotter.validate_files(ROOT, paths) == summaries[sample]["total"]
            assert frames[sample].GetNRuns() == 1
        binning = plotter.optimal_binning(list(summaries.values()), "2022")
        bins, low, high = binning
        assert low < 0.8 and high > 1.5
        histograms = {sample: plotter.book_histogram(frame, sample, binning)
                      for sample, frame in frames.items()}
        for sample, hist in histograms.items():
            assert hist.Integral() == hist.GetEntries() == summaries[sample]["total"]
            assert hist.GetBinContent(0) == hist.GetBinContent(hist.GetNbinsX() + 1) == 0.
            for b in range(1, hist.GetNbinsX() + 1):
                assert math.isclose(hist.GetBinError(b) ** 2, hist.GetBinContent(b), abs_tol=1e-12)
            assert frames[sample].GetNRuns() == 2
        assert histograms["Private"].GetBinContent(histograms["Private"].FindBin(0.8)) == 2.
        # Use ROOT's actual bin finder to verify JSON doubles and stored
        # floats: all 17 protected corrections, including unity, are distinct.
        peaks = (*plotter.lookup_weights("2022"), 1.)
        float_file, double_file = directory / "float.root", directory / "double.root"
        write_file(float_file, peaks, include_weight=False)
        write_file(double_file, peaks, double=True)
        # TChain requires one branch type within a sample; check Float_t
        # and Double_t samples separately against the same histogram grid.
        assert plotter.validate_files(ROOT, [str(float_file)]) == 17
        assert plotter.validate_files(ROOT, [str(double_file)]) == 17
        peak_frame, peak_summary = plotter.inspect_sample(ROOT, [str(float_file)], "peaks")
        peak_grid = plotter.optimal_binning([peak_summary], "2022")
        peak_hist = plotter.book_histogram(peak_frame, "peaks", peak_grid)
        double_frame, _ = plotter.inspect_sample(ROOT, [str(double_file)], "doubles")
        double_hist = plotter.book_histogram(double_frame, "doubles", peak_grid)
        peak_hist.Add(double_hist)
        occupied = [peak_hist.FindBin(value) for value in peaks]
        assert len(set(occupied)) == len(peaks)
        assert peak_hist.Integral() == peak_hist.GetEntries() == 34
        for b in occupied:
            assert peak_hist.GetBinContent(b) == 2.
            assert math.isclose(peak_hist.GetBinError(b) ** 2, 2.)
        for b in range(1, peak_hist.GetNbinsX() + 1):
            assert math.isclose(peak_hist.GetBinWidth(b),
                                (peak_grid[2] - peak_grid[1]) / peak_grid[0])
        assert sum(peak_hist.GetBinContent(b) > 0
                   for b in range(1, peak_hist.GetNbinsX() + 1)) == 17
        for value in (float("nan"), float("inf")):
            bad = directory / "bad.root"
            write_file(bad, [value])
            expect_error("non-finite w_nnlo", plotter.inspect_sample, ROOT, [str(bad)], "Private")
        signed = directory / "signed.root"
        write_file(signed, [-1., 0., 1.])
        _, summary = plotter.inspect_sample(ROOT, [str(signed)], "Private")
        assert summary["total"] == 3 and summary["min"] == -1.
        expect_error("no events", plotter.validate_files, ROOT, [str(empty)])
        expect_error("no events", plotter.inspect_sample, ROOT, [str(empty)], "Private")
        missing, vector, wrong = [directory / name for name in ("missing.root", "vector.root", "wrong.root")]
        write_file(missing, [1.], branch_name="other")
        write_file(vector, [1.], vector=True)
        write_file(wrong, [1.], tree_name="Events")
        for path in (missing, vector):
            expect_error("w_nnlo", plotter.validate_files, ROOT, [str(path)])
        expect_error("Missing pico tree", plotter.validate_files, ROOT, [str(wrong)])
        with patch.object(sys, "argv", ["plot_ggF_wnnlo.py", "2022"]), \
             patch.object(plotter, "list_pico_files", side_effect=[(files["Central"], "central"),
                                                                   (files["Private"], "private")]), \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(io.StringIO()):
            plotter.main()
        assert draw.call_count == 1
        call = draw.call_args
        assert call.args[4].name == "ggF_wnnlo_2022.png"
        assert call.kwargs["event_counts"] == {"Central": (3,), "Private": (4,)}
        assert call.kwargs["reference_lines"] == {"w_nnlo": 1.}
        assert call.kwargs["panels"][0][3:] == binning
        # Render realistic lookup peaks as well as the out-of-table inputs
        # above, so visual QA exercises the optimized 84-bin grid.
        float_hist = peak_hist.Clone("central_lookup_peaks")
        private_hist = peak_hist.Clone("private_lookup_peaks")
        private_hist.SetBinContent(private_hist.FindBin(peaks[0]), 4.)
        private_hist.SetBinError(private_hist.FindBin(peaks[0]), 2.)
        private_hist.SetEntries(36.)
        output = Path("/tmp/ggF_wnnlo_synthetic_check.png")
        plotter.draw_plots(*call.args[:4], output, **call.kwargs)
        assert output.is_file() and output.stat().st_size > 0
        for hist in call.args[1] + call.args[2]:
            assert math.isclose(hist.Integral(), 1.)
        peak_output = Path("/tmp/ggF_wnnlo_lookup_bins_check.png")
        plotter.draw_plots(ROOT, [float_hist], [private_hist], "2022", peak_output,
                           event_counts={"Central": (34,), "Private": (36,)},
                           panels=(("w_nnlo", "NNLO correction factor", "w_{NNLO}", *peak_grid),),
                           columns=1, panel_size=(900, 800), legend_text_size=0.040,
                           reference_lines={"w_nnlo": 1.}, y_title_x=0.11)
        assert peak_output.is_file() and peak_output.stat().st_size > 0
        print("Unit-weight w_nnlo comparison checks passed; previews: " +
              str(output) + ", " + str(peak_output))


if __name__ == "__main__":
    main()
