#!/usr/bin/env python3
"""Synthetic checks for unweighted numbered NNLO correction categories."""

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
    # A non-finite stored weight must not influence this diagnostic.
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
    values = plotter.category_values("2022")
    assert len(values) == 17
    assert values == tuple(sorted((*plotter.LOOKUP_WEIGHTS["run3"], 1.)))
    assert values[0] == 0.9251468016571035
    assert values[7] == 1.
    assert values[-1] == 1.6961873040033373
    assert plotter.lookup_weights("2016APV") == plotter.LOOKUP_WEIGHTS["run2"]

    with tempfile.TemporaryDirectory(prefix="ggf_wnnlo_") as tmp:
        directory = Path(tmp)
        central, private, empty = [directory / name for name in
                                   ("central.root", "private.root", "empty.root")]
        write_file(central, [values[0], 1., values[-1]])
        write_file(private, [values[0], values[0], 1., values[-1]], include_weight=False)
        write_file(empty, [])
        files = {"Central": [str(central), str(empty)], "Private": [str(private)]}
        histograms, summaries = {}, {}
        for sample, paths in files.items():
            frame, summary = plotter.inspect_sample(ROOT, paths, sample)
            assert plotter.validate_files(ROOT, paths) == summary["total"]
            assert frame.GetNRuns() == 1
            hist = plotter.book_histogram(frame, sample, "2022")
            assert frame.GetNRuns() == 2
            assert hist.Integral() == hist.GetEntries() == summary["total"]
            assert hist.GetNbinsX() == 17
            assert hist.GetBinContent(0) == hist.GetBinContent(18) == 0.
            for b in range(1, 18):
                assert hist.GetBinWidth(b) == 1.
                assert hist.GetXaxis().GetBinLabel(b) == str(b)
                assert math.isclose(hist.GetBinError(b) ** 2, hist.GetBinContent(b),
                                    abs_tol=1e-12)
            assert hist.GetBinContent(8) == 1.
            histograms[sample], summaries[sample] = hist, summary
        assert histograms["Private"].GetBinContent(1) == 2.
        assert histograms["Central"].GetBinContent(17) == 1.

        # Every JSON double and Float_t representation maps to exactly one
        # category; ROOT chains are kept homogeneous in branch type.
        for year in ("2022", "2018"):
            expected = plotter.category_values(year)
            peak_histograms = []
            for double in (False, True):
                path = directory / ("peaks_double.root" if double else "peaks_float.root")
                write_file(path, expected, double=double, include_weight=False)
                assert plotter.validate_files(ROOT, [str(path)]) == 17
                frame, _ = plotter.inspect_sample(ROOT, [str(path)], "peaks")
                hist = plotter.book_histogram(frame, "peaks", year)
                assert hist.Integral() == hist.GetEntries() == 17
                for b in range(1, 18):
                    assert hist.GetBinContent(b) == hist.GetBinError(b) == 1.
                peak_histograms.append(hist)
            for b in range(19):
                assert peak_histograms[0].GetBinContent(b) == peak_histograms[1].GetBinContent(b)

        bad = directory / "bad.root"
        for value in (float("nan"), float("inf")):
            write_file(bad, [value])
            expect_error("non-finite w_nnlo", plotter.inspect_sample, ROOT, [str(bad)], "bad")
        # Do not silently discard unexpected/rescaled values or assign them
        # to the nearest known correction.
        for value in (-1., 0., 0.8, 1.2, values[0] + 1.e-4):
            write_file(bad, [value])
            frame, _ = plotter.inspect_sample(ROOT, [str(bad)], "unknown")
            expect_error("unknown w_nnlo", plotter.book_histogram, frame, "unknown", "2022")
        expect_error("no events", plotter.validate_files, ROOT, [str(empty)])
        expect_error("no events", plotter.inspect_sample, ROOT, [str(empty)], "empty")
        missing, vector, wrong = [directory / name for name in
                                  ("missing.root", "vector.root", "wrong.root")]
        write_file(missing, [1.], branch_name="other")
        write_file(vector, [1.], vector=True)
        write_file(wrong, [1.], tree_name="Events")
        for path in (missing, vector):
            expect_error("w_nnlo", plotter.validate_files, ROOT, [str(path)])
        expect_error("Missing pico tree", plotter.validate_files, ROOT, [str(wrong)])

        # Presentation-style preview with varied populations in all 17 bins.
        central_counts = (8, 9, 12, 10, 9, 8, 7, 1, 4, 8, 10, 3, 2, 1, 1, 1, 1)
        private_counts = (10, 11, 11, 9, 8, 7, 6, 1, 4, 8, 10, 3, 2, 1, 1, 1, 1)
        write_file(central, [value for value, count in zip(values, central_counts)
                             for _ in range(count)])
        write_file(private, [value for value, count in zip(values, private_counts)
                             for _ in range(count)], include_weight=False)
        output_text = io.StringIO()
        with patch.object(sys, "argv", ["plot_ggF_wnnlo.py", "2022"]), \
             patch.object(plotter, "list_pico_files",
                          side_effect=[(files["Central"], "central"),
                                       (files["Private"], "private")]), \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(output_text):
            plotter.main()
        assert draw.call_count == 1
        assert "0.9251468016571035" in output_text.getvalue()
        assert "unity reference" in output_text.getvalue()
        call = draw.call_args
        assert call.args[4].name == "ggF_wnnlo_2022.png"
        assert call.kwargs["event_counts"] == {
            "Central": (sum(central_counts),), "Private": (sum(private_counts),)}
        assert call.kwargs["reference_lines"] == {"w_nnlo": 8.}
        assert call.kwargs["panels"][0][3:] == (17, 0.5, 17.5)
        output = Path("/tmp/ggF_wnnlo_synthetic_check.png")
        plotter.draw_plots(*call.args[:4], output, **call.kwargs)
        assert output.is_file() and output.stat().st_size > 0
        for hist in call.args[1] + call.args[2]:
            assert math.isclose(hist.Integral(), 1.)
            assert hist.GetNbinsX() == 17  # No automatic label-axis expansion.
        print("Numbered unweighted w_nnlo checks passed; preview: " + str(output))


if __name__ == "__main__":
    main()
