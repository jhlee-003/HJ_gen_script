#!/usr/bin/env python3
"""Synthetic tests of both-HJ NNLO removal and the six-panel layout; no EOS."""

from contextlib import redirect_stdout
import io
import math
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_HJ_kinematics_gridpack as plotter
import plot_HJ_BDT_var_gridpack as bdt
import plot_ggF_kinematics_wnnlo_removed as reference
from check_kinematics_plotter import event, write_file, expect_error
from check_wnnlo_removed_plotter import add_nnlo


def main():
    assert plotter.PANELS == reference.PANELS
    assert plotter.REQUIRED_BRANCHES == reference.REQUIRED_BRANCHES
    assert plotter.CPP_HELPERS == reference.CPP_HELPERS
    assert plotter.SAMPLE_DIRECTORIES is bdt.SAMPLE_DIRECTORIES
    assert plotter.parse_arguments(["2022"])[1].year == "2022"
    assert plotter.SAMPLE_DIRECTORIES["2022"] == {
        "previous": "/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico",
        "new": "/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022PrvGridpico"}
    rows = []
    for weight, pt in ((2., 40.), (6., 80.), (-0.5, 20.)):
        row = event()
        row.update(weight=weight, photon_pt=[10., pt])
        rows.append(row)
    with tempfile.TemporaryDirectory(prefix="hj_kinematics_gridpack_") as tmp:
        directory = Path(tmp)
        previous, new, missing = [directory / name for name in ("previous.root", "new.root", "missing.root")]
        for filename in (previous, new, missing):
            write_file(filename, rows)
        add_nnlo(previous, [2., 3., 0.5])
        add_nnlo(new, [1., 2., 0.25])
        for sample, filename, expected_sum, weights in (
                ("Previous", previous, 2., (1., 2., -1.)),
                ("New", new, 3., (2., 3., -2.))):
            assert plotter.validate_files(ROOT, [str(filename)]) == 3
            histograms, summary = plotter.book_histograms(ROOT, [str(filename)], sample)
            assert summary["sumw"] == expected_sum
            assert summary["total"] == 3 and summary["counts"] == (3,) * 6
            assert summary["event_loops"] == 1
            for pt, weight in zip((40., 80., 20.), weights):
                b = histograms[0].FindBin(pt)
                assert histograms[0].GetBinContent(b) == weight
                assert math.isclose(histograms[0].GetBinError(b) ** 2, weight ** 2)
            expected_hists, expected_summary = reference.book_histograms(ROOT, [str(filename)], sample, True)
            assert summary == expected_summary
            for hist, expected_hist in zip(histograms, expected_hists):
                for b in range(hist.GetNbinsX() + 2):
                    assert hist.GetBinContent(b) == expected_hist.GetBinContent(b)
                    assert hist.GetBinError(b) == expected_hist.GetBinError(b)
            expect_error(ValueError, "w_nnlo", plotter.validate_files, ROOT, [str(missing)])
            for divisor in (0., float("nan"), float("inf")):
                invalid = directory / (sample + "_invalid.root")
                write_file(invalid, [event()])
                add_nnlo(invalid, [divisor])
                expect_error(ValueError, "invalid w_nnlo", plotter.book_histograms,
                             ROOT, [str(invalid)], sample)
        with patch.object(sys, "argv", ["plot_HJ_kinematics_gridpack.py", "2022"]), \
             patch.object(plotter, "list_pico_files", side_effect=[([str(previous)], "previous"),
                                                                   ([str(new)], "new")]) as discover, \
             patch.object(plotter, "validate_files", wraps=plotter.validate_files) as validate, \
             patch.object(plotter, "book_histograms", wraps=plotter.book_histograms) as book, \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(io.StringIO()):
            plotter.main()
        assert [call.args[0] for call in discover.call_args_list] == [
            plotter.SAMPLE_DIRECTORIES["2022"][key] for key in ("previous", "new")]
        for function in (validate, book):
            assert [call.kwargs["remove_nnlo"] for call in function.call_args_list] == [True, True]
        assert draw.call_count == 1
        call = draw.call_args
        assert call.args[4].name == "HJ_kinematics_gridpack_2022.png"
        assert call.kwargs["panels"] == reference.PANELS
        assert call.kwargs["columns"] == 3 and call.kwargs["panel_size"] == (600, 700)
        assert call.kwargs["legend_text_size"] == 0.04
        assert call.kwargs["reference_lines"] == {"m_ll": 91.2, "m_llgamma": 125.}
        assert call.kwargs["event_counts"] == {"Previous": (3,) * 6, "New": (3,) * 6}
        output = Path("/tmp/HJ_kinematics_gridpack_synthetic_check.png")
        original_close = ROOT.TCanvas.Close
        def inspect_canvas(canvas):
            for index in range(1, 7):
                top = canvas.GetPad(index).GetPrimitive("hj_top_" + str(index))
                legends = [p for p in top.GetListOfPrimitives() if p.InheritsFrom("TLegend")]
                assert [str(entry.GetLabel()) for entry in legends[0].GetListOfPrimitives()] == [
                    " HJ (central gridpack)", " HJ (private gridpack)"]
                assert math.isclose(legends[0].GetTextSize(), 0.04, rel_tol=1e-6)
                lines = [p for p in top.GetListOfPrimitives() if p.InheritsFrom("TLine")]
                assert len(lines) == (1 if index in (4, 5) else 0)
                if lines:
                    assert lines[0].GetX1() == lines[0].GetX2() == {4: 91.2, 5: 125.}[index]
                    assert lines[0].GetLineColor() == ROOT.kBlack and lines[0].GetLineStyle() == 2
                bottom = canvas.GetPad(index).GetPrimitive("hj_ratio_" + str(index))
                assert "New / Previous" in [p.GetTitle() for p in bottom.GetListOfPrimitives()
                                            if p.InheritsFrom("TLatex")]
                axes = bottom.GetPrimitive("hj_ratio_axes_" + str(index))
                assert axes.GetMinimum() == 0. and axes.GetMaximum() == 2.
            original_close(canvas)
        with patch.object(ROOT.TCanvas, "Close", inspect_canvas):
            plotter.draw_plots(*call.args[:4], output, **call.kwargs)
        assert output.is_file() and output.stat().st_size > 0
        previous_hist, new_hist = call.args[1][0], call.args[2][0]
        assert math.isclose(previous_hist.Integral(), 1.) and math.isclose(new_hist.Integral(), 1.)
        assert previous_hist.GetBinContent(previous_hist.FindBin(20.)) == -0.5
        assert math.isclose(new_hist.GetBinContent(new_hist.FindBin(20.)), -2. / 3.)
        ratio = bdt.make_ratio(ROOT, new_hist, previous_hist)
        assert math.isclose(ratio.GetPointY(0), 4. / 3.)
        print("Both-HJ kinematics checks passed; preview: " + str(output))


if __name__ == "__main__":
    main()
