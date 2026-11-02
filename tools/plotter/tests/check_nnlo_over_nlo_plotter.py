#!/usr/bin/env python3
"""Check NNLO/NLO shape ratios and lookup overlays without EOS access."""

from contextlib import redirect_stdout
import io
import math
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_ggF_NNLO_over_NLO as plotter

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT(2)


def expect_error(phrase, function, *args):
    try:
        function(*args)
    except ValueError as error:
        assert phrase in str(error), str(error)
    else:
        raise AssertionError("Expected " + phrase)


def write_file(path, counts):
    source = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree("tree", "Synthetic last-Higgs records")
    ids, pts = ROOT.std.vector("int")(), ROOT.std.vector("float")()
    tree.Branch("mc_id", ids)
    tree.Branch("mc_pt", pts)
    # There is deliberately no weight or w_nnlo branch.
    for index, count in enumerate(counts):
        midpoint = (plotter.PT_BIN_EDGES[index] + plotter.PT_BIN_EDGES[index + 1]) / 2.
        for _ in range(count):
            ids.clear()
            pts.clear()
            for value in (25, 22, 25):
                ids.push_back(value)
            for value in (499., 7., midpoint):
                pts.push_back(value)
            tree.Fill()
    # These events must not enter the shape normalization.
    ids.clear()
    pts.clear()
    tree.Fill()
    ids.push_back(25)
    pts.push_back(500.)
    tree.Fill()
    tree.Write()
    source.Close()


def main():
    assert plotter.parse_arguments(["2022"])[1].year == "2022"
    zeros = [0] * 14
    weights = [2., 1.] + [1.] * 14
    central, private = [100, 100] + zeros, [200, 100] + zeros
    rows, summary = plotter.build_comparison(central, private, weights)
    assert summary["central_total"] == 200 and summary["private_total"] == 300
    assert summary["mean_lookup"] == 1.5
    assert math.isclose(rows[0]["ratio"], 4. / 3.)
    assert math.isclose(rows[1]["ratio"], 2. / 3.)
    for row in rows[:2]:
        assert math.isclose(row["ratio"], row["shape_lookup"])
        assert math.isclose(row["closure"], 1.)
    fc, fh = 0.5, 2. / 3.
    variance = fh * (1. - fh) / 300 / fc ** 2 + fh ** 2 * fc * (1. - fc) / 200 / fc ** 4
    assert math.isclose(rows[0]["error"] ** 2, variance)
    assert summary["undefined_bins"] == 14
    observed, closure = plotter.make_graphs(ROOT, rows)
    assert observed.GetN() == closure.GetN() == 2
    scaled, _ = plotter.build_comparison([n * 2 for n in central],
                                         [n * 2 for n in private], weights)
    assert math.isclose(scaled[0]["ratio"], rows[0]["ratio"])
    assert math.isclose(scaled[0]["error"], rows[0]["error"] / math.sqrt(2.))
    expect_error("intervals", plotter.build_comparison, central[:-1], private, weights)
    expect_error("no events", plotter.build_comparison, [0] * 16, private, weights)
    expect_error("integer", plotter.build_comparison, [-1, 101] + zeros, private, weights)
    expect_error("integer", plotter.build_comparison, [1.5, 100] + zeros, private, weights)
    expect_error("positive", plotter.build_comparison, central, private, [0.] * 16)

    with tempfile.TemporaryDirectory(prefix="ggf_nnlo_over_nlo_") as tmp:
        directory = Path(tmp)
        central_path, private_path = directory / "central.root", directory / "private.root"
        nc = tuple(20 + i for i in range(16))
        nh = tuple(22 + i + i % 3 for i in range(16))
        write_file(central_path, nc)
        write_file(private_path, nh)
        text = io.StringIO()
        with patch.object(sys, "argv", ["plot_ggF_NNLO_over_NLO.py", "2022"]), \
             patch.object(plotter, "list_pico_files", side_effect=[
                 ([str(central_path)], "central"), ([str(private_path)], "private")]), \
             patch.object(plotter, "draw_comparison") as draw, redirect_stdout(text):
            plotter.main()
        assert draw.call_count == 1
        args = draw.call_args.args
        assert args[3].name == "ggF_NNLO_over_NLO_2022.png"
        measured, summary = args[1:3]
        assert summary["central_total"] == sum(nc)
        assert summary["private_total"] == sum(nh)
        assert summary["undefined_bins"] == 0
        assert [row["central"] for row in measured] == list(nc)
        assert [row["private"] for row in measured] == list(nh)
        assert [row["lookup"] for row in measured] == list(plotter.lookup_weights("2022"))
        assert "no Higgs=1" in text.getvalue()
        assert "pT >= 500 GeV omitted=1" in text.getvalue()
        output = Path("/tmp/ggF_NNLO_over_NLO_synthetic_check.png")
        plotter.draw_comparison(ROOT, measured, summary, output)
        assert output.is_file() and output.stat().st_size > 0
        print("NNLO/NLO comparison checks passed; preview: " + str(output))


if __name__ == "__main__":
    main()
