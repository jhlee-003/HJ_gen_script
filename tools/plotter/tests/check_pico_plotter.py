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
import plot_ggF_BDT_var as plotter
from plot_ggF_BDT_var import normalize_histogram, make_ratio

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
        "jet_eta": [9., -1.5, 1.2], "jet_phi": [0., -0.8, -2.4],
        "jet_m": [100., 8., 12.], "photon_pt": [1., 40.],
        "photon_idmva": [-0.9, 0.8], "photon_energyErr": [100., 2.],
        "llphoton_psi": [-1.2, 2.], "llphoton_phi": [1.4, 0.],
        "el_pt": [1., 20., 40.], "mu_pt": [1.],
        "photon_mht_dphi": [0.2, 2.4], "ht": 200., "mht": 60.,
    }


def write_file(path, rows, omit=(), tree_name="tree"):
    source = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree(tree_name, "Synthetic baseline pico")
    buffers = {}
    for name in plotter.REQUIRED_BRANCHES:
        if name in omit:
            continue
        if name in ("weight", "ht", "mht"):
            buffers[name] = array("f", [0.])
            tree.Branch(name, buffers[name], name + "/F")
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
            if name in ("weight", "njet", "ht", "mht"):
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
    assert keys == ["pt_over_mass", "photon_idmva", "photon_rel_energy_err", "dr_min", "dr_max",
                    "cos_Theta", "cos_theta", "phi", "lepton1_eta", "lepton2_eta", "photon_eta",
                    "photon_mht_dphi", "photon_jet_dr", "leading_jet_eta", "leading_jet_mass",
                    "leading_jet_pt", "llgamma_jet_dphi", "system_balance", "photon_zeppenfeld"]
    assert [panel[3:] for panel in plotter.PANELS] == [
        (40, 0., 2.5), (45, 0., 1.), (40, 0.01, 0.25), (35, 0., 3.5), (60, 0., 6.),
        (40, -1., 1.), (40, -1., 1.), (40, -3.2, 3.2), (40, -2.6, 2.6), (40, -2.6, 2.6),
        (40, -2.6, 2.6), (40, 0., 3.15), (40, 0.4, 6.), (50, -5., 5.), (50, 0., 40.),
        (40, 30., 150.), (40, 0., 3.15), (40, 0., 1.), (40, 0., 6.)]
    assert [len(group) for group in plotter.PLOT_GROUPS] == [8, 8, 3]

    electron = event()
    muon = event()
    muon.update(weight=-0.5, llphoton_iph=[0], llphoton_ill=[0],
                photon_eta=[0.7], photon_phi=[-3.05], ll_i1=[0], ll_i2=[1],
                ll_lepid=[13], mu_eta=[-0.5, 1.2], mu_phi=[3.05, 0.9], mu_pt=[50., 25.],
                photon_pt=[40.], photon_idmva=[0.8], photon_energyErr=[2.], photon_mht_dphi=[2.4])
    nojet = event()
    nojet.update(weight=3., njet=0, jet_isgood=[False] * 3, llphoton_pt=[0.])
    bad_photon = event()
    bad_photon.update(weight=0.25, llphoton_iph=[99])
    bad_angles = event()
    bad_angles.update(weight=1., llphoton_cosTheta=[float("nan")], llphoton_costheta=[-999.])
    empty = {name: [] for name in plotter.REQUIRED_BRANCHES}
    empty.update(weight=1., njet=0, ht=0., mht=0.)
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
        assert summary["counts"] == (5, 4, 4, 4, 4, 4, 4, 5, 5, 5, 4, 4, 3, 4, 4, 4, 4, 5, 3), summary
        recoil_bin = hist["pt_over_mass"].FindBin(30. / 125.)
        assert math.isclose(hist["pt_over_mass"].GetBinContent(recoil_bin), 2.75)
        assert math.isclose(hist["pt_over_mass"].GetBinError(recoil_bin) ** 2, 5.3125)
        # Good jets, not the 500 GeV bad jet; stored candidate, not candidate 1.
        assert hist["leading_jet_pt"].GetBinContent(hist["leading_jet_pt"].FindBin(80.)) == 2.75
        assert hist["leading_jet_pt"].GetBinContent(hist["leading_jet_pt"].GetNbinsX() + 1) == 0.
        assert hist["photon_eta"].GetBinContent(hist["photon_eta"].GetNbinsX() + 1) == 0.
        expected_dr = min(math.hypot(0.4 + 0.7, 0.7 - 0.2), math.hypot(0.4 - 0.8, 0.7 - 2.))
        # The wrapped-phi muon example lands in the same bin with weight -0.5.
        assert hist["dr_min"].GetBinContent(hist["dr_min"].FindBin(expected_dr)) == 5.5
        # Stored negative-lepton/axis conventions are not recomputed or sign-flipped.
        assert hist["cos_theta"].GetBinContent(hist["cos_theta"].FindBin(0.45)) == 4.75
        # New observables use the candidate's indices and the highest-pT good jet.
        expr = (
            "hjpicoplot::build(llphoton_pt, llphoton_m, llphoton_iph, llphoton_ill, "
            "photon_eta, photon_phi, ll_i1, ll_i2, ll_lepid, el_eta, el_phi, mu_eta, mu_phi, "
            "njet, jet_pt, jet_isgood, llphoton_cosTheta, llphoton_costheta, "
            "photon_pt, photon_idmva, photon_energyErr, llphoton_psi, el_pt, mu_pt, "
            "jet_eta, jet_phi, jet_m, llphoton_phi, photon_mht_dphi, ht, mht)")
        # Separate one-row file avoids Range's incompatibility with implicit MT.
        one = directory / "one.root"
        write_file(one, [electron])
        values = ROOT.RDataFrame("tree", str(one)).Define("values", expr)
        expected = {"photon_idmva": 0.8, "photon_rel_energy_err": 2. / (40. * math.cosh(0.4)),
                    "phi": -1.2, "lepton1_eta": 0.8, "lepton2_eta": -0.7,
                    "photon_mht_dphi": 2.4, "leading_jet_eta": 1.2, "leading_jet_mass": 12.,
                    "leading_jet_pt": 80., "llgamma_jet_dphi": 2.4831853071795864,
                    "system_balance": 0.3, "photon_zeppenfeld": 0.8,
                    "photon_jet_dr": math.hypot(0.4 - 1.2, 0.7 + 2.4)}
        for name, result in expected.items():
            actual = float(values.Define("check_" + name, "values." + name).Mean("check_" + name).GetValue())
            assert math.isclose(actual, result, rel_tol=1.e-6, abs_tol=1.e-6), (name, actual, result)

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
        signed_pt = signed_hists[0]
        normalize_histogram(signed_pt)
        assert signed_pt.GetBinContent(signed_pt.FindBin(10. / 125.)) < 0.
        assert math.isclose(signed_pt.Integral(), 1.)
        # Check overflow merging, preserved counts and propagated errors.
        jet_hist = hist["leading_jet_pt"].Clone("overflow_test")
        entries = jet_hist.GetEntries()
        jet_hist.Fill(250., 2.)
        normalize_histogram(jet_hist)
        assert jet_hist.GetBinContent(jet_hist.GetNbinsX() + 1) == 0. and jet_hist.GetEntries() == entries + 1
        assert math.isclose(jet_hist.Integral(), 1.)

        # Validate main's automatic routes, no extra CLI arguments or selection.
        with patch.object(sys, "argv", ["plot_ggF_BDT_var.py", "2022"]), \
             patch.object(plotter, "list_pico_files", return_value=(files, "unused")) as discovery, \
             patch.object(plotter, "draw_plots") as draw, redirect_stdout(io.StringIO()):
            plotter.main()
        assert [item.args[0] for item in discovery.call_args_list] == [
            plotter.SAMPLE_DIRECTORIES["2022"]["central"], plotter.SAMPLE_DIRECTORIES["2022"]["private"]]
        assert all(not item.kwargs for item in discovery.call_args_list)
        assert draw.call_count == 3
        for page, call in enumerate(draw.call_args_list, start=1):
            assert call.args[4].name == "ggF_BDT_var_2022_{}.png".format(page)
            assert call.args[4].parent == Path(__file__).resolve().parents[3] / "plots"
            start, stop = ((0, 8), (8, 16), (16, 19))[page - 1]
            assert call.kwargs["panels"] == plotter.PANELS[start:stop]
            assert call.kwargs["event_counts"] == {
                "Central": summary["counts"][start:stop], "Private": summary["counts"][start:stop]}

        central = [h.Clone("central_" + key) for key, h in zip(keys, histograms)]
        private = [h.Clone("private_" + key) for key, h in zip(keys, histograms)]
        original_close = ROOT.TCanvas.Close
        def inspect_canvas(canvas):
            columns = min(4, len(group))
            rows = math.ceil(len(group) / columns)
            assert math.isclose(canvas.GetWw() / canvas.GetWh(), 1.5 * columns / rows, rel_tol=0.06)
            assert len(canvas.GetListOfPrimitives()) == len(group)
            for index in range(1, len(group) + 1):
                cell = canvas.GetPad(index)
                top = cell.GetPrimitive("hj_top_" + str(index))
                bottom = cell.GetPrimitive("hj_ratio_" + str(index))
                axes = bottom.GetPrimitive("hj_ratio_axes_" + str(index))
                assert axes.GetMinimum() == 0. and axes.GetMaximum() == 2.
                assert axes.GetXaxis().GetXmin() == group[index - 1][4]
                assert axes.GetXaxis().GetXmax() == group[index - 1][5]
                titles = [pad.GetPrimitive("hj_y_title_" + kind + "_" + str(index))
                          for pad, kind in ((top, "top"), (bottom, "ratio"))]
                assert titles[0].GetX() == titles[1].GetX() == 0.065
                legends = [p for p in top.GetListOfPrimitives() if p.InheritsFrom("TLegend")]
                assert len(legends) == 2
                assert legends[0].GetX1NDC() == 0.24
                assert legends[1].GetX2NDC() == 0.895
                assert legends[0].GetY2NDC() == legends[1].GetY2NDC() == 0.88
                assert legends[0].GetY1NDC() == legends[1].GetY1NDC() == 0.72
                assert all(math.isclose(legend.GetTextSize(), 0.060, rel_tol=1.e-6)
                           for legend in legends)
                assert legends[0].GetX2NDC() <= legends[1].GetX1NDC()
                assert legends[1].GetTextAlign() == 32
                labels = [entry.GetLabel() for legend in legends for entry in legend.GetListOfPrimitives()]
                assert "2022 ggF signal sample" not in labels
                assert " [Central] ggH_qme" in labels and " [Priavate] HJ (MiNNLO)" in labels
                assert labels[-2:] == ["N={:,}".format(summary["counts"][start + index - 1])] * 2
                assert not any(label.startswith("Jobs=") for label in labels)
                vertical_lines = [p for p in top.GetListOfPrimitives()
                                  if p.InheritsFrom("TLine") and p.GetX1() == p.GetX2()]
                assert not vertical_lines
            original_close(canvas)
        start = 0
        for page, group in enumerate(plotter.PLOT_GROUPS, start=1):
            stop = start + len(group)
            output = Path("/tmp/ggF_BDT_var_synthetic_check_{}.png".format(page))
            with patch.object(ROOT.TCanvas, "Close", inspect_canvas):
                plotter.draw_plots(ROOT, central[start:stop], private[start:stop], "2022", output,
                                   event_counts={"Central": summary["counts"][start:stop],
                                                 "Private": summary["counts"][start:stop]}, panels=group)
            assert output.is_file() and output.stat().st_size > 0
            start = stop
        ratio = make_ratio(ROOT, private[0], central[0])
        assert all(math.isclose(ratio.GetPointY(i), 1.) for i in range(ratio.GetN()))
        print("Synthetic pico checks passed; previews: /tmp/ggF_BDT_var_synthetic_check_{1,2,3}.png")


if __name__ == "__main__":
    main()
