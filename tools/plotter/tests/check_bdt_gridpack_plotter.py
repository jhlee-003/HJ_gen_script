#!/usr/bin/env python3
"""Synthetic two-HJ gridpack comparison checks; never accesses EOS."""

from contextlib import redirect_stdout
import io
import math
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_ggF_BDT_var_wnnlo_removed as reference
import plot_ggF_BDT_var_gridpack_wnnlo_removed as plotter
from check_pico_plotter import event, write_file, expect_error
from check_wnnlo_removed_plotter import add_nnlo


def main():
    assert plotter.PANELS == reference.PANELS
    assert plotter.PLOT_GROUPS == reference.PLOT_GROUPS
    assert plotter.CPP_HELPERS == reference.CPP_HELPERS
    assert plotter.parse_arguments(["2022"])[1].year == "2022"
    paths = plotter.SAMPLE_DIRECTORIES["2022"]
    assert set(paths) == {"previous", "new"}
    assert paths["previous"] == "/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico"
    assert paths["new"] == "/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022PrvGridpico"
    rows = []
    for weight, pt in ((2., 30.), (6., 80.), (-0.5, 150.)):
        row = event()
        row.update(weight=weight, llphoton_pt=[pt])
        rows.append(row)
    with tempfile.TemporaryDirectory(prefix="ggf_bdt_gridpack_") as tmp:
        directory = Path(tmp)
        previous, new, missing = [directory / name for name in ("previous.root", "new.root", "missing.root")]
        for filename in (previous, new, missing):
            write_file(filename, rows)
        add_nnlo(previous, [2., 3., 0.5])
        add_nnlo(new, [1., 2., 0.25])
        samples, summaries = {}, {}
        for sample, filename, expected_sum, expected_weights in (
                ("Previous", previous, 2., (1., 2., -1.)),
                ("New", new, 3., (2., 3., -2.))):
            assert plotter.validate_files(ROOT, [str(filename)]) == 3
            histograms, summary = plotter.book_histograms(ROOT, [str(filename)], sample)
            assert summary["sumw"] == expected_sum
            assert summary["total"] == 3 and summary["counts"] == (3,) * 19
            assert summary["event_loops"] == 1
            for pt, weight in zip((30., 80., 150.), expected_weights):
                b = histograms[0].FindBin(pt / 125.)
                assert histograms[0].GetBinContent(b) == weight
                assert math.isclose(histograms[0].GetBinError(b) ** 2, weight ** 2)
            expected_hists, expected_summary = reference.book_histograms(ROOT, [str(filename)], sample, True)
            assert summary == expected_summary
            for hist, expected_hist in zip(histograms, expected_hists):
                for b in range(hist.GetNbinsX() + 2):
                    assert hist.GetBinContent(b) == expected_hist.GetBinContent(b)
                    assert hist.GetBinError(b) == expected_hist.GetBinError(b)
            samples[sample], summaries[sample] = histograms, summary
            expect_error(ValueError, "w_nnlo", plotter.validate_files, ROOT, [str(missing)])
            for divisor in (0., float("nan"), float("inf")):
                invalid = directory / (sample + "_invalid.root")
                write_file(invalid, [event()])
                add_nnlo(invalid, [divisor])
                expect_error(ValueError, "invalid w_nnlo", plotter.book_histograms,
                             ROOT, [str(invalid)], sample)
        with patch.object(sys, "argv", ["plot_ggF_BDT_var_gridpack_wnnlo_removed.py", "2022"]), \
             patch.object(plotter, "list_pico_files", side_effect=[([str(previous)], "previous"),
                                                                   ([str(new)], "new")]) as discover, \
             patch.object(plotter, "validate_files", wraps=plotter.validate_files) as validate, \
             patch.object(plotter, "book_histograms", wraps=plotter.book_histograms) as book, \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(io.StringIO()):
            plotter.main()
        assert [call.args[0] for call in discover.call_args_list] == [paths["previous"], paths["new"]]
        for function in (validate, book):
            assert [call.kwargs["remove_nnlo"] for call in function.call_args_list] == [True, True]
        assert draw.call_count == 3
        original_close = ROOT.TCanvas.Close
        def inspect_canvas(canvas):
            for index in range(1, len(canvas.GetListOfPrimitives()) + 1):
                top = canvas.GetPad(index).GetPrimitive("hj_top_" + str(index))
                legends = [p for p in top.GetListOfPrimitives() if p.InheritsFrom("TLegend")]
                labels = [str(entry.GetLabel()) for entry in legends[0].GetListOfPrimitives()]
                assert labels == [" HJ (previous)", " HJ (private gridpack)"]
                ratio = canvas.GetPad(index).GetPrimitive("hj_ratio_" + str(index))
                titles = [p.GetTitle() for p in ratio.GetListOfPrimitives() if p.InheritsFrom("TLatex")]
                assert "New / Previous" in titles
            original_close(canvas)
        for page, call in enumerate(draw.call_args_list, 1):
            assert call.args[4].name == "ggF_BDT_var_gridpack_wnnlo_removed_2022_{}.png".format(page)
            panels = call.kwargs["panels"]
            assert panels == reference.PLOT_GROUPS[page - 1]
            assert call.kwargs["y_title_x"] == 0.11
            assert call.kwargs["event_counts"] == {"Previous": (3,) * len(panels), "New": (3,) * len(panels)}
            output = Path("/tmp/ggF_BDT_gridpack_synthetic_check_{}.png".format(page))
            with patch.object(ROOT.TCanvas, "Close", inspect_canvas):
                plotter.draw_plots(*call.args[:4], output, **call.kwargs)
            assert output.is_file() and output.stat().st_size > 0
        # The rendered histograms are independently unit-normalized; the ratio
        # numerator is NEW and the denominator is PREVIOUS, never the reverse.
        previous_hist = draw.call_args_list[0].args[1][0]
        new_hist = draw.call_args_list[0].args[2][0]
        ratio = plotter.make_ratio(ROOT, new_hist, previous_hist)
        assert ratio.GetN() == 3
        assert math.isclose(ratio.GetPointY(0), 4. / 3.)
        assert math.isclose(previous_hist.Integral(), 1.)
        assert math.isclose(new_hist.Integral(), 1.)
        print("Both-HJ NNLO removal, input routing, ratio, and three rendered layouts passed.")


if __name__ == "__main__":
    main()
