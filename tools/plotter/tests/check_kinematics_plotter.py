#!/usr/bin/env python3
"""Synthetic checks for the six-variable pico plotter; no EOS access."""

from array import array
from contextlib import redirect_stdout, redirect_stderr
import io
import math
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_ggF_kinematics as plotter
from plot_ggF_BDT_var import normalize_histogram, make_ratio

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT(2)


def event():
    return {"weight": 2., "photon_pt": [10., 40.], "ll_pt": [900., 35.],
            "ll_m": [200., 91.], "llphoton_pt": [90., 1.], "llphoton_m": [125., 1.],
            "llphoton_iph": [1], "llphoton_ill": [1], "npv": 30}


def write_file(path, rows, omit=(), tree_name="tree"):
    source = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree(tree_name, "Synthetic kinematics")
    buffers = {}
    for name in plotter.REQUIRED_BRANCHES:
        if name in omit:
            continue
        if name == "weight":
            buffers[name] = array("f", [0.])
            tree.Branch(name, buffers[name], "weight/F")
        elif name == "npv":
            buffers[name] = array("i", [0])
            tree.Branch(name, buffers[name], "npv/I")
        else:
            buffers[name] = ROOT.std.vector("int" if name in ("llphoton_iph", "llphoton_ill") else "float")()
            tree.Branch(name, buffers[name])
    for row in rows:
        for name, buffer in buffers.items():
            if name in ("weight", "npv"):
                buffer[0] = row[name]
            else:
                buffer.clear()
                for value in row[name]:
                    buffer.push_back(value)
        tree.Fill()
    tree.Write()
    source.Close()


def expect_error(error_type, phrase, function, *args):
    try:
        function(*args)
    except error_type as error:
        assert phrase in str(error), str(error)
    else:
        raise AssertionError("Expected " + phrase)


def main():
    assert plotter.parse_arguments(["2022"])[1].year == "2022"
    with redirect_stderr(io.StringIO()):
        expect_error(SystemExit, "2", plotter.parse_arguments, ["2023"])
        expect_error(SystemExit, "2", plotter.parse_arguments, ["2022", "extra"])
    keys = [panel[0] for panel in plotter.PANELS]
    assert keys == ["pt_gamma", "pt_ll", "pt_llgamma", "m_ll", "m_llgamma", "npv"]
    assert [panel[3:] for panel in plotter.PANELS] == [
        (50, 0., 200.), (50, 0., 200.), (50, 0., 200.),
        (35, 50., 120.), (50, 100., 150.), (100, -0.5, 99.5)]
    negative = event()
    negative.update(weight=-0.5, photon_pt=[10., 70.], ll_pt=[900., 60.],
                    llphoton_pt=[150.], ll_m=[200., 100.], llphoton_m=[135.], npv=80)
    zero_pt = event()
    zero_pt.update(weight=1., photon_pt=[10., 0.], ll_pt=[900., 0.], llphoton_pt=[0.], npv=0)
    bad_photon = event()
    bad_photon.update(weight=0.25, llphoton_iph=[99])
    bad_pair = event()
    bad_pair.update(weight=0.5, llphoton_ill=[99])
    empty = {name: [] for name in plotter.REQUIRED_BRANCHES}
    empty.update(weight=1., npv=-1)

    with tempfile.TemporaryDirectory(prefix="ggf_kinematics_") as tmp:
        directory = Path(tmp)
        source, zero = directory / "sample.root", directory / "zero.root"
        write_file(source, [event(), negative, zero_pt, bad_photon, bad_pair, empty])
        write_file(zero, [])
        files = [str(source), str(zero)]
        assert plotter.validate_files(ROOT, files) == 6
        histograms, summary = plotter.book_histograms(ROOT, files, "test")
        hist = dict(zip(keys, histograms))
        assert summary == {"total": 6, "sumw": 4.25, "counts": (4, 4, 5, 4, 5, 5), "event_loops": 1}
        assert hist["pt_gamma"].GetBinContent(hist["pt_gamma"].FindBin(40.)) == 2.5
        assert hist["pt_ll"].GetBinContent(hist["pt_ll"].FindBin(35.)) == 2.25
        assert hist["pt_ll"].GetBinContent(hist["pt_ll"].GetNbinsX() + 1) == 0.
        mass_bin = hist["m_llgamma"].FindBin(125.)
        assert hist["m_llgamma"].GetBinContent(mass_bin) == 3.75
        assert math.isclose(hist["m_llgamma"].GetBinError(mass_bin) ** 2, 5.3125)
        assert hist["npv"].GetBinContent(hist["npv"].FindBin(30.)) == 2.75
        signed = hist["pt_gamma"].Clone("signed_check")
        normalize_histogram(signed)
        assert signed.GetBinContent(signed.FindBin(70.)) < 0.
        assert math.isclose(signed.Integral(), 1.)
        overflow = hist["npv"].Clone("overflow_check")
        overflow.Fill(200., 2.)
        entries = overflow.GetEntries()
        normalize_histogram(overflow)
        assert overflow.GetBinContent(overflow.GetNbinsX() + 1) == 0.
        assert overflow.GetEntries() == entries and math.isclose(overflow.Integral(), 1.)

        missing = directory / "missing.root"
        write_file(missing, [event()], omit=("npv",))
        expect_error(ValueError, "npv", plotter.validate_files, ROOT, [str(missing)])
        wrong = directory / "wrong.root"
        write_file(wrong, [event()], tree_name="Events")
        expect_error(ValueError, "pico tree", plotter.validate_files, ROOT, [str(wrong)])
        expect_error(ValueError, "no events", plotter.validate_files, ROOT, [str(zero)])
        nonfinite = directory / "nonfinite.root"
        row = event()
        row["weight"] = float("nan")
        write_file(nonfinite, [row])
        expect_error(ValueError, "non-finite stored weights", plotter.book_histograms,
                     ROOT, [str(nonfinite)], "bad")

        with patch.object(sys, "argv", ["plot_ggF_kinematics.py", "2022"]), \
             patch.object(plotter, "list_pico_files", return_value=(files, "unused")) as discovery, \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(io.StringIO()):
            plotter.main()
        assert discovery.call_count == 2 and draw.call_count == 1
        assert draw.call_args.args[4].name == "ggF_kinematics_2022.png"
        assert draw.call_args.args[4].parent == Path(__file__).resolve().parents[3] / "plots"
        assert draw.call_args.kwargs["columns"] == 3
        assert draw.call_args.kwargs["panel_size"] == (600, 700)
        assert draw.call_args.kwargs["event_counts"] == {"Central": summary["counts"], "Private": summary["counts"]}

        central = [h.Clone("central_" + key) for key, h in zip(keys, histograms)]
        private = [h.Clone("private_" + key) for key, h in zip(keys, histograms)]
        output = Path("/tmp/ggF_kinematics_synthetic_check.png")
        original_close = ROOT.TCanvas.Close
        def inspect_canvas(canvas):
            assert len(canvas.GetListOfPrimitives()) == 6
            assert math.isclose(canvas.GetWw() / canvas.GetWh(), 1800. / 1400., rel_tol=0.03)
            for index in range(1, 7):
                cell = canvas.GetPad(index)
                top = cell.GetPrimitive("hj_top_" + str(index))
                bottom = cell.GetPrimitive("hj_ratio_" + str(index))
                axes = bottom.GetPrimitive("hj_ratio_axes_" + str(index))
                assert axes.GetMinimum() == 0. and axes.GetMaximum() == 2.
                legends = [p for p in top.GetListOfPrimitives() if p.InheritsFrom("TLegend")]
                assert len(legends) == 2
                assert all(math.isclose(legend.GetTextSize(), 0.06, rel_tol=1.e-6) for legend in legends)
                lines = [p for p in top.GetListOfPrimitives() if p.InheritsFrom("TLine")]
                assert len(lines) == (1 if index == 5 else 0)
                if lines:
                    assert lines[0].GetX1() == lines[0].GetX2() == 125.
                    assert lines[0].GetLineColor() == ROOT.kBlack and lines[0].GetLineStyle() == 2
            original_close(canvas)
        with patch.object(ROOT.TCanvas, "Close", inspect_canvas):
            plotter.draw_plots(ROOT, central, private, "2022", output,
                               event_counts={"Central": summary["counts"], "Private": summary["counts"]},
                               panels=plotter.PANELS, columns=3, panel_size=(600, 700),
                               reference_lines={"m_llgamma": 125.})
        assert output.is_file() and output.stat().st_size > 0
        ratio = make_ratio(ROOT, private[0], central[0])
        assert all(math.isclose(ratio.GetPointY(i), 1.) for i in range(ratio.GetN()))
        print("Kinematic checks passed; preview: " + str(output))


if __name__ == "__main__":
    main()
