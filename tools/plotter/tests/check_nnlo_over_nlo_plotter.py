#!/usr/bin/env python3
"""Check signed weight/w_nnlo in both samples and the raw-lookup comparison."""

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


def float32(value):
    return array("f", [value])[0]


def event(pt, weight, nnlo):
    return {"mc_id": [25, 22, 25], "mc_pt": [499., 7., pt],
            "weight": weight, "w_nnlo": nnlo}


def sample_events(private=False):
    rows, sums, squared = [], [], []
    for i, nnlo in enumerate(plotter.lookup_weights("2022")):
        pt = (plotter.PT_BIN_EDGES[i] + plotter.PT_BIN_EDGES[i + 1]) / 2.
        targets = (1.5 + i * .07, -.1) if private else (2. + i * .05, -.25)
        actual = [float32(w * nnlo) / float32(nnlo) for w in targets]
        sums.append(sum(actual))
        squared.append(sum(w * w for w in actual))
        rows.extend(event(pt, w * nnlo, nnlo) for w in targets)
    # Large valid weights outside the plotted selection cannot change its normalization.
    rows.append({"mc_id": [], "mc_pt": [], "weight": 10000., "w_nnlo": 1.})
    rows.append(event(500., 20000., 1.))
    return rows, sums, squared


def write_file(path, rows, omit=()):
    source = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree("tree", "Synthetic weighted last-Higgs records")
    buffers = {"mc_id": ROOT.std.vector("int")(), "mc_pt": ROOT.std.vector("float")(),
               "weight": array("f", [0.]), "w_nnlo": array("f", [0.])}
    for name, buffer in buffers.items():
        if name not in omit:
            if name in ("weight", "w_nnlo"):
                tree.Branch(name, buffer, name + "/F")
            else:
                tree.Branch(name, buffer)
    for row in rows:
        for name, buffer in buffers.items():
            if name in ("weight", "w_nnlo"):
                buffer[0] = row[name]
            else:
                buffer.clear()
                for value in row[name]:
                    buffer.push_back(value)
        tree.Fill()
    tree.Write()
    source.Close()


def main():
    assert plotter.parse_arguments(["2022"])[1].year == "2022"
    zeros = [0.] * 14
    weights = [2., 1.] + [1.] * 14
    central, private = [100., 100.] + zeros, [200., 100.] + zeros
    rows, summary = plotter.build_comparison(central, private, central, private, weights)
    assert summary["central_sumw"] == 200 and summary["private_sumw"] == 300
    assert math.isclose(rows[0]["ratio"], 4. / 3.)
    assert math.isclose(rows[1]["ratio"], 2. / 3.)
    # The lower panel now uses RAW lookup, with no shape-normalization factor.
    assert math.isclose(rows[0]["closure"], 2. / 3.)
    assert math.isclose(rows[1]["closure"], 2. / 3.)
    assert all("shape_lookup" not in row for row in rows)
    assert "mean_lookup" not in summary
    fc, fh = 0.5, 2. / 3.
    variance = fh * (1. - fh) / 300 / fc ** 2 + fh ** 2 * fc * (1. - fc) / 200 / fc ** 4
    assert math.isclose(rows[0]["error"] ** 2, variance)
    observed, closure = plotter.make_graphs(ROOT, rows)
    assert observed.GetN() == closure.GetN() == 2
    # A constant rescaling of each sample's event weights cancels in shapes AND errors.
    scaled, _ = plotter.build_comparison(
        [s * 5 for s in central], [s * 7 for s in private],
        [q * 25 for q in central], [q * 49 for q in private], weights)
    assert math.isclose(scaled[0]["ratio"], rows[0]["ratio"])
    assert math.isclose(scaled[0]["error"], rows[0]["error"])
    # Negative bin sums remain signed; variance includes the shared weighted total.
    fractions, variances, total = plotter.normalized_fractions(
        [2., -.5] + zeros, [4., .25] + zeros, "signed")
    assert total == 1.5 and fractions[1] < 0
    expected = ((1. - fractions[0]) ** 2 * 4. + fractions[0] ** 2 * .25) / 1.5 ** 2
    assert math.isclose(variances[0], expected)
    expect_error("intervals", plotter.build_comparison, central[:-1], private,
                 central, private, weights)
    expect_error("non-positive", plotter.normalized_fractions, [1., -1.], [1., 1.], "zero")
    expect_error("squared", plotter.normalized_fractions, [1.], [-1.], "bad")
    expect_error("positive", plotter.build_comparison, central, private,
                 central, private, [0.] * 16)

    with tempfile.TemporaryDirectory(prefix="ggf_nnlo_over_nlo_") as tmp:
        directory = Path(tmp)
        central_path, private_path = directory / "central.root", directory / "private.root"
        central_rows, sc, qc = sample_events()
        private_rows, sh, qh = sample_events(private=True)
        write_file(central_path, central_rows)
        write_file(private_path, private_rows)
        assert plotter.validate_files(ROOT, [str(central_path)]) == 34
        histograms, stats = plotter.book_histograms(ROOT, [str(central_path)], "direct")
        assert stats["event_loops"] == 1
        assert stats["counts"] == (32,) and stats["no_higgs"] == stats["above_range"] == 1
        for b in range(1, 17):
            assert math.isclose(histograms[0].GetBinContent(b), sc[b - 1], rel_tol=1.e-12)
            assert math.isclose(histograms[0].GetBinError(b) ** 2, qc[b - 1], rel_tol=1.e-12)
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
        assert math.isclose(summary["central_sumw"], sum(sc))
        assert math.isclose(summary["private_sumw"], sum(sh))
        for i, row in enumerate(measured):
            assert math.isclose(row["central_sumw"], sc[i])
            assert math.isclose(row["private_sumw"], sh[i])
            assert math.isclose(row["central_sumw2"], qc[i])
            assert math.isclose(row["private_sumw2"], qh[i])
            assert math.isclose(row["closure"], row["ratio"] / row["lookup"])
        assert "no Higgs=1" in text.getvalue()
        assert "pT >= 500 GeV omitted=1" in text.getvalue()
        assert "shape lookup" not in text.getvalue().lower()

        bad = directory / "bad.root"
        for omit in ("weight", "w_nnlo"):
            write_file(bad, central_rows, omit=(omit,))
            expect_error(omit, plotter.validate_files, ROOT, [str(bad)])
        for weight, nnlo in ((1., 0.), (1., float("nan")), (float("inf"), 1.)):
            write_file(bad, [event(10., weight, nnlo)])
            expect_error("invalid weight/w_nnlo", plotter.book_histograms,
                         ROOT, [str(bad)], "bad")
        malformed = event(10., 1., 1.)
        malformed["mc_pt"] = [10.]
        write_file(bad, [malformed])
        expect_error("vector lengths", plotter.book_histograms, ROOT, [str(bad)], "bad")
        write_file(bad, [event(-999., 1., 1.)])
        expect_error("invalid last Higgs", plotter.book_histograms, ROOT, [str(bad)], "bad")

        output = Path("/tmp/ggF_NNLO_over_NLO_synthetic_check.png")
        # ROOTFacade dynamically resolves C++ names; use a forwarding proxy
        # to count graphics constructors without patching the facade itself.
        class GraphicsProxy:
            def __init__(self):
                self.legend_count = 0

            def TLegend(self, *args):
                self.legend_count += 1
                return ROOT.TLegend(*args)

            def __getattr__(self, name):
                return getattr(ROOT, name)

        graphics = GraphicsProxy()
        plotter.draw_comparison(graphics, measured, summary, output)
        assert graphics.legend_count == 1  # No event-size legend.
        assert output.is_file() and output.stat().st_size > 0
        print("Weighted NNLO/NLO comparison checks passed; preview: " + str(output))


if __name__ == "__main__":
    main()
