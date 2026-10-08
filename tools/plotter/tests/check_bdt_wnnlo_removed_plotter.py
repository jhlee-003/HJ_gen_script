#!/usr/bin/env python3
"""Verify HJ-only w_nnlo removal for all 19 BDT panels, without EOS."""

from contextlib import redirect_stdout
import io
import math
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_ggF_BDT_var as original
import plot_ggF_BDT_var_wnnlo_removed as plotter
from check_pico_plotter import event, write_file, expect_error
from check_wnnlo_removed_plotter import add_nnlo


def main():
    assert plotter.PANELS == original.PANELS
    assert plotter.PLOT_GROUPS == original.PLOT_GROUPS
    assert plotter.CPP_HELPERS == original.CPP_HELPERS
    assert plotter.REQUIRED_BRANCHES == original.REQUIRED_BRANCHES
    assert plotter.parse_arguments(["2022"])[1].year == "2022"
    rows = []
    for weight, pt in ((2., 30.), (6., 80.), (-0.5, 150.)):
        row = event()
        row.update(weight=weight, llphoton_pt=[pt])
        rows.append(row)
    with tempfile.TemporaryDirectory(prefix="ggf_bdt_wnnlo_removed_") as tmp:
        directory = Path(tmp)
        central, private = directory / "central.root", directory / "private.root"
        write_file(central, rows)
        write_file(private, rows)
        add_nnlo(private, [2., 3., 0.5])
        assert plotter.validate_files(ROOT, [str(central)]) == 3
        assert plotter.validate_files(ROOT, [str(private)], True) == 3
        expect_error(ValueError, "w_nnlo", plotter.validate_files, ROOT, [str(central)], True)
        central_hists, central_summary = plotter.book_histograms(ROOT, [str(central)], "central")
        private_hists, private_summary = plotter.book_histograms(ROOT, [str(private)], "private", True)
        assert central_summary["sumw"] == 7.5 and private_summary["sumw"] == 2.
        for summary in (central_summary, private_summary):
            assert summary["total"] == 3 and summary["counts"] == (3,) * 19
            assert summary["event_loops"] == 1
        for pt, cw, pw in ((30., 2., 1.), (80., 6., 2.), (150., -0.5, -1.)):
            for hist, weight in ((central_hists[0], cw), (private_hists[0], pw)):
                b = hist.FindBin(pt / 125.)
                assert hist.GetBinContent(b) == weight
                assert math.isclose(hist.GetBinError(b) ** 2, weight ** 2)
        # Every central histogram must match the original implementation exactly.
        references, reference_summary = original.book_histograms(ROOT, [str(central)], "reference")
        assert central_summary == reference_summary
        for hist, reference in zip(central_hists, references):
            for b in range(hist.GetNbinsX() + 2):
                assert hist.GetBinContent(b) == reference.GetBinContent(b)
                assert hist.GetBinError(b) == reference.GetBinError(b)
        for value in (0., float("nan"), float("inf")):
            bad = directory / "bad.root"
            write_file(bad, [event()])
            add_nnlo(bad, [value])
            expect_error(ValueError, "invalid w_nnlo", plotter.book_histograms,
                         ROOT, [str(bad)], "private", True)
            _, summary = plotter.book_histograms(ROOT, [str(bad)], "central")
            assert summary["sumw"] == 2. and summary["counts"] == (1,) * 19
        negative = directory / "negative.root"
        write_file(negative, [event()])
        add_nnlo(negative, [-2.])
        _, summary = plotter.book_histograms(ROOT, [str(negative)], "private", True)
        assert summary["sumw"] == -1.
        with patch.object(sys, "argv", ["plot_ggF_BDT_var_wnnlo_removed.py", "2022"]), \
             patch.object(plotter, "list_pico_files", side_effect=[([str(central)], "central"),
                                                                   ([str(private)], "private")]), \
             patch.object(plotter, "book_histograms", wraps=plotter.book_histograms) as book, \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(io.StringIO()):
            plotter.main()
        assert [call.kwargs["remove_nnlo"] for call in book.call_args_list] == [False, True]
        assert draw.call_count == 3
        for page, call in enumerate(draw.call_args_list, 1):
            assert call.args[4].name == "ggF_BDT_var_wnnlo_removed_2022_{}.png".format(page)
            assert call.kwargs["panels"] == original.PLOT_GROUPS[page - 1]
            assert call.kwargs["y_title_x"] == 0.11
            assert call.kwargs["event_counts"] == {
                "Central": (3,) * len(call.kwargs["panels"]),
                "Private": (3,) * len(call.kwargs["panels"])}
        start = 0
        for page, panels in enumerate(plotter.PLOT_GROUPS, 1):
            stop = start + len(panels)
            output = Path("/tmp/ggF_BDT_var_wnnlo_removed_synthetic_check_{}.png".format(page))
            plotter.draw_plots(ROOT, central_hists[start:stop], private_hists[start:stop],
                               "2022", output, panels=panels, y_title_x=0.11,
                               event_counts={"Central": (3,) * len(panels),
                                             "Private": (3,) * len(panels)})
            assert output.is_file() and output.stat().st_size > 0
            start = stop
        assert private_hists[0].GetBinContent(private_hists[0].FindBin(150. / 125.)) == -0.5
        for hist in central_hists + private_hists:
            assert math.isclose(hist.Integral(), 1.)
        print("BDT HJ-only w_nnlo removal checks passed; all three previews rendered.")


if __name__ == "__main__":
    main()
