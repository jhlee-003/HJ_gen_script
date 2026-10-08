#!/usr/bin/env python3
"""Synthetic baseline-pico checks; requires PyROOT, never accesses EOS."""

from array import array
from contextlib import redirect_stderr, redirect_stdout
import copy
import io
import math
from pathlib import Path
import subprocess
import sys
import tempfile
from unittest.mock import patch

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_pico_kinematics as plotter
from plot_pico_kinematics import normalize_histogram, make_ratio

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.EnableImplicitMT(2)


def event():
    # Candidate references deliberately do not point to photon/pair zero.
    # No use_event/bitmap/ID/charge branches: no selection is reapplied.
    return {
        "weight": 2., "llphoton_pt": [30., 90.], "llphoton_m": [125., 200.],
        "llphoton_cosTheta": [-0.3, 0.9], "llphoton_costheta": [0.45, -0.9],
        "llphoton_iph": [1, 0], "llphoton_ill": [1, 0],
        "photon_eta": [9., 0.4], "photon_phi": [0., 0.7],
        "ll_i1": [0, 1], "ll_i2": [0, 2], "ll_lepid": [13, 11],
        "el_eta": [8., -0.7, 0.8], "el_phi": [0., 0.2, 2.],
        "mu_eta": [9.], "mu_phi": [0.], "njet": 2,
        "jet_pt": [500., 45., 80.], "jet_isgood": [False, True, True],
    }


def write_file(path, rows, omit=(), tree_name="tree"):
    source = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree(tree_name, "Synthetic baseline pico")
    buffers = {}
    for name in plotter.REQUIRED_BRANCHES:
        if name in omit:
            continue
        if name == "weight":
            buffers[name] = array("f", [0.])
            tree.Branch(name, buffers[name], "weight/F")
        elif name == "njet":
            buffers[name] = array("i", [0])
            tree.Branch(name, buffers[name], "njet/I")
        else:
            kind = "bool" if name == "jet_isgood" else (
                "int" if name in ("llphoton_iph", "llphoton_ill", "ll_i1", "ll_i2", "ll_lepid") else "float")
            buffers[name] = ROOT.std.vector(kind)()
            tree.Branch(name, buffers[name])
    for row in rows:
        for name, buffer in buffers.items():
            if name in ("weight", "njet"):
                buffer[0] = row[name]
            else:
                buffer.clear()
                for item in row[name]:
                    buffer.push_back(item)
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
    assert "plot_nanoaod_kinematics" not in sys.modules
    assert plotter.parse_arguments(["2022"])[1].year == "2022"
    for arguments in (["2023"], ["2022", "extra"]):
        with redirect_stderr(io.StringIO()):
            expect_error(SystemExit, "2", plotter.parse_arguments, arguments)
    keys = [panel[0] for panel in plotter.PANELS]
    assert keys == ["pt_over_mass", "photon_eta", "dr_min", "dr_max", "n_jets",
                    "leading_jet_pt", "m_llgamma", "cos_Theta", "cos_theta", "pt_llgamma"]
    ranges = {panel[0]: panel[4:] for panel in plotter.PANELS}
    assert ranges["m_llgamma"] == (100., 150.)
    assert ranges["pt_llgamma"] == (0., 200.)
    assert [panel[3] for panel in plotter.PANELS] == [50, 50, 50, 50, 10, 50, 60, 50, 50, 50]

    electron = event()
    muon = event()
    muon.update(weight=-0.5, llphoton_iph=[0], llphoton_ill=[0],
                photon_eta=[0.7], photon_phi=[-3.05], ll_i1=[0], ll_i2=[1],
                ll_lepid=[13], mu_eta=[-0.5, 1.2], mu_phi=[3.05, 0.9])
    nojet = event()
    nojet.update(weight=3., njet=0, jet_isgood=[False] * 3, llphoton_pt=[0.])
    bad_photon = event()
    bad_photon.update(weight=0.25, llphoton_iph=[99])
    bad_angles = event()
    bad_angles.update(weight=1., llphoton_cosTheta=[float("nan")], llphoton_costheta=[-999.])
    empty = {name: [] for name in plotter.REQUIRED_BRANCHES}
    empty.update(weight=1., njet=0)
    rows = [electron, muon, nojet, bad_photon, bad_angles, empty]

    with tempfile.TemporaryDirectory(prefix="hj_pico_plotter_") as tmp:
        directory = Path(tmp)
        first, second, zero = [directory / name for name in ("first.root", "second.root", "zero.root")]
        write_file(first, rows[:3])
        write_file(second, rows[3:])
        write_file(zero, [])
        files = [str(path) for path in (first, second, zero)]
        nested = directory / "nested"
        nested.mkdir()
        write_file(nested / "excluded.root", [electron])
        with patch.object(plotter, "eos_location", return_value=("root://eosuser.cern.ch", "/eos/user/j/junhyuk/pico", directory)):
            discovered, label = plotter.list_pico_files("/eos/user/j/junhyuk/pico")
        assert [name.rsplit("/", 1)[1] for name in discovered] == ["first.root", "second.root", "zero.root"]
        assert label == "pico"
        remote_listing = subprocess.CompletedProcess([], 0, "b.root\n/eos/user/j/junhyuk/pico/a.root\nlogs\n", "")
        with patch.object(plotter.shutil, "which", return_value="xrdfs"), \
             patch.object(plotter.subprocess, "run", return_value=remote_listing) as xrdfs:
            remote_files, _ = plotter.list_pico_files("root://eosuser.cern.ch//eos/user/j/junhyuk/pico")
        assert remote_files == ["root://eosuser.cern.ch//eos/user/j/junhyuk/pico/a.root",
                                "root://eosuser.cern.ch//eos/user/j/junhyuk/pico/b.root"]
        assert xrdfs.call_count == 1
        assert plotter.validate_files(ROOT, files) == 6
        histograms, summary = plotter.book_histograms(ROOT, files, "synthetic")
        hist = dict(zip(keys, histograms))
        assert summary["total"] == 6 and summary["event_loops"] == 1
        assert math.isclose(summary["sumw"], 6.75)
        assert summary["counts"] == (5, 4, 4, 4, 6, 4, 5, 4, 4, 5), summary
        mass_bin = hist["m_llgamma"].FindBin(125.)
        assert math.isclose(hist["m_llgamma"].GetBinContent(mass_bin), 5.75)
        assert math.isclose(hist["m_llgamma"].GetBinError(mass_bin) ** 2, 14.3125)
        # Good jets, not the 500 GeV bad jet; stored candidate, not candidate 1.
        assert hist["leading_jet_pt"].GetBinContent(hist["leading_jet_pt"].FindBin(80.)) == 2.75
        assert hist["leading_jet_pt"].GetBinContent(hist["leading_jet_pt"].GetNbinsX() + 1) == 0.
        assert hist["photon_eta"].GetBinContent(hist["photon_eta"].GetNbinsX() + 1) == 0.
        expected_dr = min(math.hypot(0.4 + 0.7, 0.7 - 0.2), math.hypot(0.4 - 0.8, 0.7 - 2.))
        # The wrapped-phi muon example lands in the same bin with weight -0.5.
        assert hist["dr_min"].GetBinContent(hist["dr_min"].FindBin(expected_dr)) == 5.5
        # Stored negative-lepton/axis conventions are not recomputed or sign-flipped.
        assert hist["cos_theta"].GetBinContent(hist["cos_theta"].FindBin(0.45)) == 4.75

        missing = directory / "missing.root"
        write_file(missing, [electron], omit=("llphoton_iph",))
        expect_error(ValueError, "llphoton_iph", plotter.validate_files, ROOT, [str(missing)])
        wrong = directory / "wrong.root"
        write_file(wrong, [electron], tree_name="Events")
        expect_error(ValueError, "pico tree", plotter.validate_files, ROOT, [str(wrong)])
        expect_error(ValueError, "no events", plotter.validate_files, ROOT, [str(zero)])
        nonfinite = directory / "nonfinite.root"
        row = copy.deepcopy(electron)
        row["weight"] = float("nan")
        write_file(nonfinite, [row])
        expect_error(ValueError, "non-finite stored weights", plotter.book_histograms, ROOT, [str(nonfinite)], "bad")

        signed = directory / "signed.root"
        positive, negative = event(), event()
        positive.update(weight=2., llphoton_pt=[100.])
        negative.update(weight=-0.5, llphoton_pt=[10.])
        write_file(signed, [positive, negative])
        signed_hists, _ = plotter.book_histograms(ROOT, [str(signed)], "signed")
        signed_pt = signed_hists[-1]
        normalize_histogram(signed_pt)
        assert signed_pt.GetBinContent(signed_pt.FindBin(10.)) < 0.
        assert math.isclose(signed_pt.Integral(), 1.)
        # Check overflow merging, preserved counts and propagated errors.
        jet_hist = hist["n_jets"].Clone("overflow_test")
        entries = jet_hist.GetEntries()
        jet_hist.Fill(12., 2.)
        normalize_histogram(jet_hist)
        assert jet_hist.GetBinContent(11) == 0. and jet_hist.GetEntries() == entries + 1
        assert math.isclose(jet_hist.Integral(), 1.)

        # Validate main's automatic routes, no extra CLI arguments or selection.
        with patch.object(sys, "argv", ["plot_pico_kinematics.py", "2022"]), \
             patch.object(plotter, "list_pico_files", return_value=(files, "unused")) as discovery, \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(io.StringIO()):
            plotter.main()
        assert [item.args[0] for item in discovery.call_args_list] == [
            plotter.SAMPLE_DIRECTORIES["2022"]["central"], plotter.SAMPLE_DIRECTORIES["2022"]["private"]]
        assert all(not item.kwargs for item in discovery.call_args_list)
        assert draw.call_args.args[4].name == "HJ_2022_pico_central_vs_private.png"
        assert draw.call_args.args[4].parent == Path(__file__).resolve().parents[3] / "plots"
        assert draw.call_args.kwargs["job_counts"] == {"Central": 3, "Private": 3}

        central = [h.Clone("central_" + key) for key, h in zip(keys, histograms)]
        private = [h.Clone("private_" + key) for key, h in zip(keys, histograms)]
        output = Path("/tmp/HJ_pico_synthetic_check.png")
        original_close = ROOT.TCanvas.Close
        def inspect_canvas(canvas):
            for index in range(1, 11):
                cell = canvas.GetPad(index)
                top = cell.GetPrimitive("hj_top_" + str(index))
                bottom = cell.GetPrimitive("hj_ratio_" + str(index))
                axes = bottom.GetPrimitive("hj_ratio_axes_" + str(index))
                assert axes.GetMinimum() == 0. and axes.GetMaximum() == 2.
                assert axes.GetXaxis().GetXmin() == plotter.PANELS[index - 1][4]
                assert axes.GetXaxis().GetXmax() == plotter.PANELS[index - 1][5]
                titles = [pad.GetPrimitive("hj_y_title_" + kind + "_" + str(index))
                          for pad, kind in ((top, "top"), (bottom, "ratio"))]
                assert titles[0].GetX() == titles[1].GetX() == 0.065
                legends = [p for p in top.GetListOfPrimitives() if p.InheritsFrom("TLegend")]
                assert len(legends) == 2
                assert legends[0].GetX1NDC() == 0.22
                assert legends[1].GetX2NDC() == 0.945
                assert legends[0].GetY2NDC() == legends[1].GetY2NDC() == 0.91
                assert legends[0].GetX2NDC() < legends[1].GetX1NDC()
                assert legends[1].GetTextAlign() == 32
                labels = [entry.GetLabel() for legend in legends for entry in legend.GetListOfPrimitives()]
                assert "2022 ggF signal sample" not in labels
                assert " ggH_qme (central)" in labels and " HJ (private)" in labels
                assert "Jobs=19" in labels and "Jobs=100" in labels
            original_close(canvas)
        with patch.object(ROOT.TCanvas, "Close", inspect_canvas):
            plotter.draw_plots(ROOT, central, private, "2022", output,
                               job_counts={"Central": 19, "Private": 100}, panels=plotter.PANELS)
        assert output.is_file() and output.stat().st_size > 0
        ratio = make_ratio(ROOT, private[0], central[0])
        assert all(math.isclose(ratio.GetPointY(i), 1.) for i in range(ratio.GetN()))
        print("Synthetic pico checks passed; preview: " + str(output))


if __name__ == "__main__":
    main()
