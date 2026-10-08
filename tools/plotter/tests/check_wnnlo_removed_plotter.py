#!/usr/bin/env python3
"""Check the HJ-only removal of w_nnlo using synthetic picos, without EOS."""

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
import plot_ggF_kinematics as original
import plot_ggF_kinematics_wnnlo_removed as plotter
from check_kinematics_plotter import event, write_file, expect_error


def add_nnlo(path, values):
    source = ROOT.TFile(str(path), "UPDATE")
    tree = source.Get("tree")
    assert tree.GetEntries() == len(values)
    buffer = array("f", [0.])
    branch = tree.Branch("w_nnlo", buffer, "w_nnlo/F")
    for value in values:
        buffer[0] = value
        branch.Fill()
    tree.Write("", ROOT.TObject.kOverwrite)
    source.Close()


def main():
    assert plotter.PANELS == original.PANELS
    assert plotter.CPP_HELPERS == original.CPP_HELPERS
    assert plotter.REQUIRED_BRANCHES == original.REQUIRED_BRANCHES
    assert plotter.parse_arguments(["2022"])[1].year == "2022"
    rows = []
    for weight, pt in ((2., 40.), (6., 80.), (-0.5, 20.)):
        row = event()
        row.update(weight=weight, photon_pt=[10., pt])
        rows.append(row)
    with tempfile.TemporaryDirectory(prefix="ggf_wnnlo_removed_") as tmp:
        directory = Path(tmp)
        central, private = directory / "central.root", directory / "private.root"
        write_file(central, rows)
        write_file(private, rows)
        add_nnlo(private, [2., 3., 0.5])
        assert plotter.validate_files(ROOT, [str(central)]) == 3
        assert plotter.validate_files(ROOT, [str(private)], remove_nnlo=True) == 3
        expect_error(ValueError, "w_nnlo", plotter.validate_files, ROOT, [str(central)], True)
        central_hists, central_summary = plotter.book_histograms(ROOT, [str(central)], "central")
        private_hists, private_summary = plotter.book_histograms(ROOT, [str(private)], "private", remove_nnlo=True)
        assert central_summary["sumw"] == 7.5 and private_summary["sumw"] == 2.
        for summary in (central_summary, private_summary):
            assert summary["total"] == 3 and summary["counts"] == (3,) * 6
            assert summary["event_loops"] == 1
        for pt, central_weight, private_weight in ((40., 2., 1.), (80., 6., 2.), (20., -0.5, -1.)):
            for hist, weight in ((central_hists[0], central_weight), (private_hists[0], private_weight)):
                bin_index = hist.FindBin(pt)
                assert hist.GetBinContent(bin_index) == weight
                assert math.isclose(hist.GetBinError(bin_index) ** 2, weight ** 2)
        # Compare central output against the unchanged original implementation.
        original_hists, original_summary = original.book_histograms(ROOT, [str(central)], "reference")
        assert central_summary == original_summary
        for central_hist, reference in zip(central_hists, original_hists):
            for b in range(central_hist.GetNbinsX() + 2):
                assert central_hist.GetBinContent(b) == reference.GetBinContent(b)
                assert central_hist.GetBinError(b) == reference.GetBinError(b)
        for value in (0., float("nan"), float("inf")):
            bad = directory / "bad.root"
            write_file(bad, [event()])
            add_nnlo(bad, [value])
            expect_error(ValueError, "invalid w_nnlo", plotter.book_histograms,
                         ROOT, [str(bad)], "private", True)
            # Central ignores w_nnlo even if that branch is present and invalid.
            _, summary = plotter.book_histograms(ROOT, [str(bad)], "central")
            assert summary["sumw"] == 2. and summary["counts"] == (1,) * 6
        with patch.object(sys, "argv", ["plot_ggF_kinematics_wnnlo_removed.py", "2022"]), \
             patch.object(plotter, "list_pico_files", side_effect=[([str(central)], "central"),
                                                                   ([str(private)], "private")]), \
             patch.object(plotter, "book_histograms", wraps=plotter.book_histograms) as book, \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(io.StringIO()):
            plotter.main()
        assert [call.kwargs["remove_nnlo"] for call in book.call_args_list] == [False, True]
        assert draw.call_count == 1
        assert draw.call_args.args[4].name == "ggF_kinematics_wnnlo_removed_2022.png"
        assert draw.call_args.kwargs["panels"] == original.PANELS
        assert draw.call_args.kwargs["legend_text_size"] == 0.04
        assert draw.call_args.kwargs["reference_lines"] == {"m_ll": 91.2, "m_llgamma": 125.}
        assert draw.call_args.kwargs["event_counts"] == {"Central": (3,) * 6, "Private": (3,) * 6}
        output = Path("/tmp/ggF_kinematics_wnnlo_removed_synthetic_check.png")
        original_close = ROOT.TCanvas.Close
        def inspect_canvas(canvas):
            for index in range(1, 7):
                top = canvas.GetPad(index).GetPrimitive("hj_top_" + str(index))
                lines = [p for p in top.GetListOfPrimitives() if p.InheritsFrom("TLine")]
                assert len(lines) == (1 if index in (4, 5) else 0)
                if lines:
                    assert lines[0].GetX1() == lines[0].GetX2() == {4: 91.2, 5: 125.}[index]
                    assert lines[0].GetLineColor() == ROOT.kBlack and lines[0].GetLineStyle() == 2
            original_close(canvas)
        with patch.object(ROOT.TCanvas, "Close", inspect_canvas):
            plotter.draw_plots(ROOT, central_hists, private_hists, "2022", output,
                               event_counts={"Central": (3,) * 6, "Private": (3,) * 6},
                               panels=plotter.PANELS, columns=3, panel_size=(600, 700),
                               reference_lines={"m_ll": 91.2, "m_llgamma": 125.}, legend_text_size=0.04)
        assert output.is_file() and output.stat().st_size > 0
        # Negative weights survive normalization; the removed factor changes shape.
        assert private_hists[0].GetBinContent(private_hists[0].FindBin(20.)) == -0.5
        assert math.isclose(central_hists[0].Integral(), 1.)
        assert math.isclose(private_hists[0].Integral(), 1.)
        print("HJ-only w_nnlo removal checks passed; preview: " + str(output))


if __name__ == "__main__":
    main()
