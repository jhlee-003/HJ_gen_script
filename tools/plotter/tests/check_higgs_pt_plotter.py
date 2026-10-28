#!/usr/bin/env python3
"""Synthetic checks for the last-record generator Higgs pT comparison."""

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
import plot_ggF_Higgs_pT as plotter
from plot_ggF_BDT_var import normalize_histogram

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT(2)


def write_file(path, rows, omit=(), tree_name="tree"):
    source = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree(tree_name, "Synthetic generator Higgs records")
    buffers = {"weight": array("f", [0.]),
               "mc_id": ROOT.std.vector("int")(),
               "mc_pt": ROOT.std.vector("float")()}
    for name, buffer in buffers.items():
        if name not in omit:
            if name == "weight":
                tree.Branch(name, buffer, "weight/F")
            else:
                tree.Branch(name, buffer)
    for weight, ids, pts in rows:
        buffers["weight"][0] = weight
        for name, values in (("mc_id", ids), ("mc_pt", pts)):
            buffers[name].clear()
            for value in values:
                buffers[name].push_back(value)
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
    rows = [
        (2., [25, 22, 25, 23], [210., 10., 35., 20.]),
        (-0.5, [25, 25], [50., 75.]),
        (1., [25], [0.]),
        (3., [21, 23], [100., 30.]),
        (1., [25, 25], [float("nan"), 650.]),
        (0.5, [-25, 25, 22], [40., 120., 999.]),
        (0.5, [], []),
    ]
    with tempfile.TemporaryDirectory(prefix="ggf_higgs_pt_") as tmp:
        directory = Path(tmp)
        central, private, empty = [directory / name for name in
                                   ("central.root", "private.root", "empty.root")]
        write_file(central, rows[:3])
        write_file(private, rows[3:])
        write_file(empty, [])
        files = [str(central), str(empty), str(private)]
        assert plotter.validate_files(ROOT, files) == 7
        histograms, summary = plotter.book_histograms(ROOT, files, "test")
        assert summary == {
            "total": 7, "no_higgs": 2, "one_higgs": 2, "multiple_higgs": 3,
            "counts": (5,), "event_loops": 1,
        }, summary
        hist = histograms[0]
        assert hist.GetEntries() == 5
        # Select last matching index, not first, highest pT, or last particle.
        for pt in (35., 75., 0., 120.):
            assert hist.GetBinContent(hist.FindBin(pt)) == 1.
            assert math.isclose(hist.GetBinError(hist.FindBin(pt)) ** 2, 1.)
        assert hist.Integral(0, hist.GetNbinsX() + 1) == hist.GetEntries() == 5.
        assert hist.GetBinContent(hist.FindBin(210.)) == 0.
        assert hist.GetBinContent(hist.GetNbinsX() + 1) == 1.
        folded = hist.Clone("folded_higgs_pt")
        normalize_histogram(folded)
        assert math.isclose(folded.Integral(), 1.)
        assert math.isclose(folded.GetBinContent(folded.GetNbinsX()), 0.2)
        assert math.isclose(folded.GetBinContent(folded.FindBin(75.)), 0.2)

        # A bad final Higgs cannot be replaced with an earlier good copy.
        bad = directory / "bad.root"
        for pt in (float("nan"), float("inf"), -999.):
            write_file(bad, [(1., [25, 25], [50., pt])])
            expect_error("invalid last Higgs pT", plotter.book_histograms,
                         ROOT, [str(bad)], "bad")
        write_file(bad, [(1., [25, 25], [50.])])
        expect_error("vector lengths", plotter.book_histograms, ROOT, [str(bad)], "bad")
        # Pico weights are never read: even non-finite values are irrelevant.
        for stored_weight in (float("nan"), float("inf"), -10., 0., 1234.):
            write_file(bad, [(stored_weight, [25], [50.])])
            unit_hist, unit_summary = plotter.book_histograms(ROOT, [str(bad)], "unit")
            assert unit_summary["counts"] == (1,)
            assert unit_hist[0].GetBinContent(unit_hist[0].FindBin(50.)) == 1.
        # A missing weight branch is supported, not just ignored at fill time.
        write_file(bad, rows, omit=("weight",))
        assert plotter.validate_files(ROOT, [str(bad)]) == 7
        no_weight_hist, no_weight_summary = plotter.book_histograms(ROOT, [str(bad)], "no_weight")
        assert no_weight_summary == summary
        for b in range(hist.GetNbinsX() + 2):
            assert no_weight_hist[0].GetBinContent(b) == hist.GetBinContent(b)
        write_file(bad, [(1., [22], [50.])])
        expect_error("no saved Higgs", plotter.book_histograms, ROOT, [str(bad)], "bad")
        expect_error("no events", plotter.validate_files, ROOT, [str(empty)])
        for branch in plotter.REQUIRED_BRANCHES:
            write_file(bad, rows[:1], omit=(branch,))
            expect_error(branch, plotter.validate_files, ROOT, [str(bad)])
        write_file(bad, rows[:1], tree_name="Events")
        expect_error("Missing pico tree", plotter.validate_files, ROOT, [str(bad)])

        # Exercise full CLI routing without EOS and render with the shared style.
        write_file(private, rows)
        real_draw = plotter.draw_plots
        with patch.object(sys, "argv", ["plot_ggF_Higgs_pT.py", "2022"]), \
             patch.object(plotter, "list_pico_files",
                          side_effect=[([str(central)], "central"),
                                       ([str(private), str(empty)], "private")]), \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(io.StringIO()):
            plotter.main()
        assert draw.call_count == 1
        call = draw.call_args
        assert call.args[4].name == "ggF_Higgs_pT_2022.png"
        assert call.kwargs["event_counts"] == {"Central": (3,), "Private": (5,)}
        assert call.kwargs["panels"] == plotter.PANELS
        output = Path("/tmp/ggF_Higgs_pT_synthetic_check.png")
        real_draw(*call.args[:4], output, **call.kwargs)
        assert output.is_file() and output.stat().st_size > 0
        for plotted_hist in call.args[1] + call.args[2]:
            assert math.isclose(plotted_hist.Integral(), 1.)
        print("Generator Higgs pT checks passed; preview: " + str(output))


if __name__ == "__main__":
    main()
