#!/usr/bin/env python3
"""End-to-end synthetic NanoAOD check; requires PyROOT, no EOS access."""

import importlib.util
import math
import tempfile
import subprocess
from array import array
from pathlib import Path
from unittest.mock import patch

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.EnableImplicitMT(2)
repo = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("plotter", repo / "tools/plot_nanoaod_kinematics.py")
plotter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plotter)


def event(flavor="Electron", charge=(1, -1), photon_eta=(0.3,), photon_phi=(1.5,), jets=(), weight=1.):
    row = {name: [] for name in plotter.REQUIRED_BRANCHES if name != "genWeight"}
    row.update({flavor + "_pt": [45., 45.], flavor + "_eta": [0., 0.],
                flavor + "_phi": [0., math.pi], flavor + "_mass": [0.000511 if flavor == "Electron" else 0.10566] * 2,
                flavor + "_charge": list(charge), "Photon_pt": [30. + 20. * i for i in range(len(photon_eta))],
                "Photon_eta": list(photon_eta), "Photon_phi": list(photon_phi), "genWeight": weight})
    for name, index in (("Jet_pt", 0), ("Jet_eta", 1), ("Jet_phi", 2), ("Jet_jetId", 3)):
        row[name] = [jet[index] for jet in jets]
    return row


def write_file(path, rows):
    output = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree("Events", "Synthetic NanoAOD for testing")
    buffers = {}
    for name in plotter.REQUIRED_BRANCHES:
        if name == "genWeight":
            buffers[name] = array("f", [0.])
            tree.Branch(name, buffers[name], "genWeight/F")
        else:
            buffers[name] = ROOT.std.vector("int" if name.endswith(("_charge", "_jetId")) else "float")()
            tree.Branch(name, buffers[name])
    for row in rows:
        for name, buffer in buffers.items():
            if name == "genWeight":
                buffer[0] = row[name]
            else:
                buffer.clear()
                for value in row[name]:
                    buffer.push_back(value)
        tree.Fill()
    tree.Write()
    output.Close()


original_canvas_close = ROOT.TCanvas.Close
checked_canvases = []


def check_canvas_format(canvas):
    """Inspect actual rendered objects before ROOT closes the canvas."""
    for index in range(1, 11):
        cell = canvas.GetPad(index)
        top = cell.GetPrimitive("hj_top_" + str(index))
        bottom = cell.GetPrimitive("hj_ratio_" + str(index))
        axes = bottom.GetPrimitive("hj_ratio_axes_" + str(index))
        assert axes.GetMinimum() == 0. and axes.GetMaximum() == 2.
        upper_title = top.GetPrimitive("hj_y_title_top_" + str(index))
        ratio_title = bottom.GetPrimitive("hj_y_title_ratio_" + str(index))
        assert upper_title.GetX() == ratio_title.GetX() == 0.065
        assert upper_title.GetTextSize() == ratio_title.GetTextSize()
        assert ratio_title.GetTextSize() == 20.
        assert abs(upper_title.GetY() - (1. - top.GetTopMargin())) < 1.e-7
        assert abs(ratio_title.GetY() - (1. - bottom.GetTopMargin())) < 1.e-7
        assert upper_title.GetTextAlign() == ratio_title.GetTextAlign() == 32
        assert upper_title.GetTextFont() == ratio_title.GetTextFont() == 43
        legends = [p for p in top.GetListOfPrimitives() if p.InheritsFrom("TLegend")]
        assert len(legends) == 3
        labels = [entry.GetLabel() for legend in legends for entry in legend.GetListOfPrimitives()]
        assert "Synthetic ggF signal sample" in labels
        assert " ggH_qme (central)" in labels and " HJ (private)" in labels
        assert legends[0].GetX1NDC() == legends[1].GetX1NDC()
        assert legends[1].GetX2NDC() < legends[2].GetX1NDC()
        assert "N=10,234,567" in labels and "N=1,000,100" in labels
    checked_canvases.append(canvas.GetName())
    original_canvas_close(canvas)


def cartesian(vector):
    return ROOT.Math.PxPyPzEVector(vector.Px(), vector.Py(), vector.Pz(), vector.E())


def check_angle_geometry():
    """Known rest-frame angles, checked independently with TLorentzVector."""
    mh, mz, polar_cosine = 125., 91.1876, -0.35
    momentum = (mh * mh - mz * mz) / (2. * mh)
    z_axis = ROOT.TVector3(0.8, 0., 0.6)
    perpendicular = ROOT.TVector3(0.6, 0., -0.8)
    z = ROOT.TLorentzVector(momentum * z_axis.X(), 0., momentum * z_axis.Z(),
                            math.sqrt(momentum * momentum + mz * mz))
    photon = ROOT.TLorentzVector(-z.Px(), 0., -z.Pz(), momentum)
    # In Z rest the photon is along -z_axis. Place l- at a known angle to it.
    direction = -polar_cosine * z_axis + math.sqrt(1. - polar_cosine ** 2) * perpendicular
    negative = ROOT.TLorentzVector(direction.X() * mz / 2., 0., direction.Z() * mz / 2., mz / 2.)
    positive = ROOT.TLorentzVector(-negative.Px(), 0., -negative.Pz(), negative.E())
    negative.Boost(z.BoostVector())
    positive.Boost(z.BoostVector())
    rest_angles = ROOT.hjplot.decay_angles(cartesian(z), cartesian(photon), cartesian(negative))
    assert math.isnan(rest_angles.cos_Theta)  # Zero H lab momentum: no flight axis.
    assert abs(rest_angles.cos_theta - polar_cosine) < 1.e-12
    # Non-collinear boosts exercise full 3D boosts, not eta/phi approximations.
    for boost in (ROOT.TVector3(0., 0., 0.4), ROOT.TVector3(0.25, -0.12, 0.32)):
        z_lab, photon_lab = ROOT.TLorentzVector(z), ROOT.TLorentzVector(photon)
        negative_lab, positive_lab = ROOT.TLorentzVector(negative), ROOT.TLorentzVector(positive)
        for vector in (z_lab, photon_lab, negative_lab, positive_lab):
            vector.Boost(boost)
        angles = ROOT.hjplot.decay_angles(cartesian(z_lab), cartesian(photon_lab), cartesian(negative_lab))
        assert abs(angles.cos_Theta - z_axis.Dot(boost.Unit())) < 1.e-12
        assert abs(angles.cos_theta - polar_cosine) < 1.e-12
        swapped = ROOT.hjplot.decay_angles(cartesian(z_lab), cartesian(photon_lab), cartesian(positive_lab))
        assert abs(swapped.cos_theta + polar_cosine) < 1.e-12
        assert abs(swapped.cos_Theta - angles.cos_Theta) < 1.e-12
    # A lightlike pair has no Z rest frame, rather than an angle of zero.
    lightlike = ROOT.Math.PxPyPzEVector(10., 0., 0., 10.)
    opposite = ROOT.Math.PxPyPzEVector(-5., 0., 0., 5.)
    undefined = ROOT.hjplot.decay_angles(lightlike, opposite, lightlike)
    assert math.isnan(undefined.cos_theta) and math.isfinite(undefined.cos_Theta)


def reference_lepton_angle(row):
    """Independent boost calculation using the NanoAOD float32 inputs."""
    leptons = []
    for index in range(2):
        vector = ROOT.TLorentzVector()
        vector.SetPtEtaPhiM(*(array("f", row["Electron_" + field])[index]
                             for field in ("pt", "eta", "phi", "mass")))
        leptons.append(vector)
    z = leptons[0] + leptons[1]
    negative = leptons[0] if row["Electron_charge"][0] < 0 else leptons[1]
    photon = ROOT.TLorentzVector()
    photon.SetPtEtaPhiM(*(array("f", row["Photon_" + field])[-1]
                         for field in ("pt", "eta", "phi")), 0.)
    negative.Boost(-z.BoostVector())
    photon.Boost(-z.BoostVector())
    return negative.Vect().Dot(photon.Vect()) / (negative.P() * photon.P())


def main():
    rows = [
        event(weight=2.),
        event(flavor="Muon", jets=((60., 2., -1., 2),)),
        event(jets=((60., 2., -1., 2),), weight=-0.25),
        event(charge=(1, 1)),
        event(photon_eta=(), photon_phi=()),
        event(photon_eta=(1.5,)),
        event(photon_eta=(0.,), photon_phi=(0.,)),
        event(jets=((25., 2., -1., 2), (60., 2., -1., 2), (80., 0., 0.01, 2),
                    (100., 2., -1., 0), (350., 1., -2., 2)), weight=2.),
        event(photon_eta=(0.3, 0.3), photon_phi=(1.5, 1.5)),
    ]
    with tempfile.TemporaryDirectory(prefix="hj_plotter_check_", dir=repo) as directory:
        directory = Path(directory)
        files = [directory / "pilot__condor-1.root", directory / "production__condor-2.root"]
        write_file(files[0], rows[:4])
        write_file(files[1], rows[4:])
        # Discovery includes both clusters, ignores non-ROOT files, and uses no glob over URLs.
        listing = "\n".join(["/eos/user/j/junhyuk/2022/" + p.name for p in files] + ["/eos/user/j/junhyuk/2022/logs"])
        with patch.object(plotter.shutil, "which", return_value="xrdfs"), patch.object(plotter.subprocess, "run") as run:
            run.return_value.stdout = listing
            discovered, label = plotter.list_nanoaod_files("root://eosuser.cern.ch//eos/user/j/junhyuk/2022")
            assert len(discovered) == 2 and label == "2022"
            assert all(url.startswith("root://eosuser.cern.ch//eos/user/") for url in discovered)
        # Mounted discovery includes depth-2 files, but not depth-3 files.
        nested = directory / "0000"
        nested.mkdir()
        write_file(nested / "central.root", rows)
        deeper = nested / "too_deep"
        deeper.mkdir()
        write_file(deeper / "excluded.root", rows)
        with patch.object(plotter, "eos_location", return_value=("root://eoscms.cern.ch", "/eos/cms/test", directory)):
            assert len(plotter.list_nanoaod_files("/eos/cms/test")[0]) == 2
            depth_two = plotter.list_nanoaod_files("/eos/cms/test", maxdepth=2)[0]
            assert len(depth_two) == 3 and not any("excluded.root" in p for p in depth_two)
        # Remote discovery likewise stops after immediate subdirectories.
        remote = "/eos/cms/test"
        responses = {
            ("ls", remote): remote + "/top.root\n" + remote + "/0000\n" + remote + "/README.txt\n",
            ("stat", remote + "/0000"): "Flags: 19 (XBitSet|IsDir|IsReadable)\n",
            ("stat", remote + "/README.txt"): "Flags: 0 ()\n",
            ("ls", remote + "/0000"): remote + "/0000/central.root\n" + remote + "/0000/too_deep\n",
        }
        def run_remote(command, **kwargs):
            return subprocess.CompletedProcess(command, 0, responses[tuple(command[2:])], "")
        with patch.object(plotter.shutil, "which", return_value="xrdfs"), \
                patch.object(plotter.subprocess, "run", side_effect=run_remote) as run:
            discovered, _ = plotter.list_nanoaod_files("root://eoscms.cern.ch/" + remote, maxdepth=2)
            assert len(discovered) == 2
            assert not any("too_deep" in " ".join(call.args[0]) for call in run.call_args_list)
        assert plotter.validate_files(ROOT, [str(p) for p in files]) == 9
        histograms, counts = plotter.book_histograms(ROOT, [str(p) for p in files], sample_name="private_selected")
        assert len(histograms) == 10
        assert tuple(counts[:3]) == (9, 5, 3), counts
        assert abs(counts[3] - 5.75) < 1.e-10, counts
        jets = histograms[4]
        assert abs(jets.GetBinContent(jets.FindBin(0)) - 3.) < 1.e-10
        assert abs(jets.GetBinContent(jets.FindBin(1)) - 0.75) < 1.e-10
        assert abs(jets.GetBinContent(jets.FindBin(2)) - 2.) < 1.e-10
        leading = histograms[5]
        assert abs(leading.GetBinContent(leading.GetNbinsX() + 1) - 2.) < 1.e-10
        assert abs(ROOT.hjplot.delta_r(0., math.pi - 0.01, 0., -math.pi + 0.01) - 0.02) < 1.e-10
        assert {panel[0] for panel in plotter.PANELS}.isdisjoint({"m_ll", "photon_pt"})
        for hist in histograms[7:9]:
            assert hist.GetXaxis().GetXmin() == -1. and hist.GetXaxis().GetXmax() == 1.
            assert hist.GetEntries() == 5.
            assert abs(hist.Integral(0, hist.GetNbinsX() + 1) - 5.75) < 1.e-10
        check_angle_geometry()
        # Verify charge-based l- assignment inside candidate building, and
        # preserve the highest-pT photon tie-break despite removing its panel.
        charge_rows = [event(), event(charge=(-1, 1)), rows[-1]]
        charge_file = directory / "angle_charge_order.root"
        write_file(charge_file, charge_rows)
        test_frame = ROOT.RDataFrame("Events", str(charge_file)).Define(
            "candidate", "hjplot::build(" + ", ".join(plotter.REQUIRED_BRANCHES[:-1]) + ", false)")
        test_frame = test_frame.Define("theta", "candidate.cos_theta").Define("chosen_photon_pt", "candidate.photon_pt")
        theta_action = test_frame.Take["double"]("theta")
        photon_action = test_frame.Take["double"]("chosen_photon_pt")
        theta_values = list(theta_action.GetValue())
        assert all(abs(actual - reference_lepton_angle(row)) < 1.e-10
                   for actual, row in zip(theta_values, charge_rows))
        assert abs(theta_values[0] + theta_values[1]) < 1.e-10
        assert list(photon_action.GetValue()) == [30., 30., 50.]
        # Undefined angles are filtered only for their own histogram, with no
        # effect on the candidate count or the remaining plotted variables.
        degenerate_file = directory / "undefined_z_rest_frame.root"
        write_file(degenerate_file, [dict(event(), Electron_phi=[0., 0.], Electron_mass=[0., 0.])])
        degenerate_hist, degenerate_counts = plotter.book_histograms(
            ROOT, [str(degenerate_file)], False, "undefined_angle")
        assert tuple(degenerate_counts[:2]) == (1, 1)
        assert degenerate_hist[8].GetEntries() == 0.
        assert all(h.GetEntries() == 1. for i, h in enumerate(degenerate_hist) if i not in (5, 8))
        # A second independent sample exercises repeated C++ declarations and
        # action names; scaling every weight leaves normalized shapes unchanged.
        central_rows = [dict(row, genWeight=3. * row["genWeight"]) for row in rows]
        central_file = directory / "central_comparison.root"
        write_file(central_file, central_rows)
        central_histograms, central_counts = plotter.book_histograms(
            ROOT, [str(central_file)], sample_name="central_selected")
        assert abs(central_counts[3] - 3. * counts[3]) < 1.e-10
        output = repo / "plots/HJ_synthetic_test_central_vs_private_selected.png"
        # Stress the count column with production-scale text, without changing
        # the small synthetic histograms' bin contents or sumw2 errors.
        for hist in central_histograms:
            hist.SetEntries(10_234_567)
        for hist in histograms:
            hist.SetEntries(1_000_100)
        with patch.object(ROOT.TCanvas, "Close", check_canvas_format):
            plotter.draw_plots(ROOT, central_histograms, histograms, "Synthetic", output)
        assert output.stat().st_size > 10_000
        assert all(abs(h.Integral() - 1.) < 1.e-10 for h in histograms)
        assert leading.GetBinContent(leading.GetNbinsX()) > 0.
        for private, central in zip(histograms, central_histograms):
            ratio = plotter.make_ratio(ROOT, private, central)
            assert all(abs(ratio.GetPointY(p) - 1.) < 1.e-10 for p in range(ratio.GetN()))
        raw_histograms, raw_counts = plotter.book_histograms(
            ROOT, [str(p) for p in files], False, "private_raw")
        raw_central, _ = plotter.book_histograms(ROOT, [str(central_file)], False, "central_raw")
        assert tuple(raw_counts[:3]) == (9, 7, 3), raw_counts
        assert abs(raw_counts[3] - 7.75) < 1.e-10
        assert raw_histograms[4].GetBinContent(raw_histograms[4].FindBin(5)) == 2.
        raw_output = repo / "plots/HJ_synthetic_test_central_vs_private_no_selection.png"
        for hist in raw_central:
            hist.SetEntries(10_234_567)
        for hist in raw_histograms:
            hist.SetEntries(1_000_100)
        with patch.object(ROOT.TCanvas, "Close", check_canvas_format):
            plotter.draw_plots(ROOT, raw_central, raw_histograms, "Synthetic", raw_output, False)
        assert len(checked_canvases) == 2
        assert raw_output.stat().st_size > 10_000
        assert all(abs(h.Integral() - 1.) < 1.e-10 for h in raw_histograms)
        # Each cut is disabled in the raw version; OS same-flavor and existence
        # requirements remain. Large/low masses should be folded, not rejected.
        changes = [
            {"Electron_pt": [8., 6.]}, {"Electron_eta": [3., 3.]},
            {"Photon_pt": [5.]}, {"Photon_eta": [3.]},
            {"Electron_pt": [100., 100.]}, {"Electron_phi": [0., 0.3]},
            {"Photon_pt": [600.]},
        ]
        outside_rows = [dict(event(), **change) for change in changes]
        absent_pair = event()
        for field in ("pt", "eta", "phi", "mass", "charge"):
            absent_pair["Muon_" + field] = [absent_pair["Electron_" + field].pop()]
        outside_rows.extend([absent_pair, event(weight=float("nan"))])
        outside_file = directory / "outside_cuts.root"
        write_file(outside_file, outside_rows)
        _, outside_selected = plotter.book_histograms(ROOT, [str(outside_file)], True, "outside_selected")
        outside_hist, outside_raw = plotter.book_histograms(ROOT, [str(outside_file)], False, "outside_raw")
        assert outside_selected[1] == 0, outside_selected
        assert tuple(outside_raw[:2]) == (8, 7), outside_raw
        assert outside_hist[6].GetBinContent(outside_hist[6].GetNbinsX() + 1) > 0.
        # Signed bin contents, squared-weight errors, and negative ratios survive.
        signed = ROOT.TH1D("signed_check", "", 3, 0., 3.)
        denominator = ROOT.TH1D("denominator_check", "", 3, 0., 3.)
        signed.Fill(0.5, -1.)
        signed.Fill(1.5, 3.)
        denominator.Fill(0.5, 2.)
        denominator.Fill(1.5, 2.)
        plotter.normalize_histogram(signed)
        plotter.normalize_histogram(denominator)
        ratio = plotter.make_ratio(ROOT, signed, denominator)
        assert ratio.GetN() == 2  # Third denominator bin is empty.
        assert ratio.GetPointY(0) == -1.
        assert abs(ratio.GetErrorY(0) - math.sqrt(2.)) < 1.e-10
        assert signed.GetBinContent(1) == -0.5 and signed.GetBinError(1) == 0.5
        empty = ROOT.TH1D("empty_check", "", 10, 0., 100.)
        plotter.normalize_histogram(empty)
        assert empty.Integral() == 0. and empty.GetEntries() == 0.
        # An input with a wrong schema must cause an explicit failure.
        bad = directory / "wrong_schema.root"
        source = ROOT.TFile(str(bad), "RECREATE")
        ROOT.TTree("Events", "wrong schema").Write()
        source.Close()
        try:
            plotter.validate_files(ROOT, [str(bad)])
        except ValueError as error:
            assert "Missing NanoAOD branches" in str(error)
        else:
            raise AssertionError("Wrong-schema file was not rejected")
    print("PASS: both modes, ten variables, multi-file/maxdepth-2 discovery, OS ee/mumu pairing,")
    print("      rest-frame angles/charge signs, signed weights/errors/ratios, cuts,")
    print("      cleaned/unsorted jets, overflow, top-aligned titles, legends and PNGs.")
    print("Preview: " + str(output))
    print("Preview: " + str(raw_output))


if __name__ == "__main__":
    main()
