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


def write_file(path, values, branch_name="w_nnlo", tree_name="tree", vector=False):
    source = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree(tree_name, "Synthetic correction factors")
    buffer = ROOT.std.vector("float")() if vector else array("f", [0.])
    if vector:
        tree.Branch(branch_name, buffer)
    else:
        tree.Branch(branch_name, buffer, branch_name + "/F")
    # An invalid stored weight proves it is NOT used to weight the histogram.
    weight = array("f", [float("nan")])
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
        low, high = plotter.common_range(list(summaries.values()))
        assert low < 0.8 and high > 1.5
        histograms = {sample: plotter.book_histogram(frame, sample, low, high)
                      for sample, frame in frames.items()}
        for sample, hist in histograms.items():
            assert hist.Integral() == hist.GetEntries() == summaries[sample]["total"]
            assert hist.GetBinContent(0) == hist.GetBinContent(hist.GetNbinsX() + 1) == 0.
            for b in range(1, hist.GetNbinsX() + 1):
                assert math.isclose(hist.GetBinError(b) ** 2, hist.GetBinContent(b), abs_tol=1e-12)
            assert frames[sample].GetNRuns() == 2
        assert histograms["Private"].GetBinContent(histograms["Private"].FindBin(0.8)) == 2.
        assert plotter.common_range([{"min": 1., "max": 1.}]) == (0.95, 1.05)
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
        output = Path("/tmp/ggF_wnnlo_synthetic_check.png")
        plotter.draw_plots(*call.args[:4], output, **call.kwargs)
        assert output.is_file() and output.stat().st_size > 0
        for hist in call.args[1] + call.args[2]:
            assert math.isclose(hist.Integral(), 1.)
        print("Unit-weight w_nnlo comparison checks passed; preview: " + str(output))


if __name__ == "__main__":
    main()
