#!/usr/bin/env python3
"""Compare ten reconstructed H -> Z gamma variables in central/private NanoAOD.

Usage: python3 tools/plot_nanoaod_kinematics.py /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022
Requires PyROOT (for example, the Python environment supplied by CMSSW).
The argument is the private EOS directory; the 2022 central directory is fixed.
Every private top-level ROOT file is included, irrespective of cluster ID.
Central discovery uses maxdepth 2 (the directory and its immediate children).
Output: plots/HJ_<directory-name>_central_vs_private_selected.png.
See tools/README.md for the candidate selection and histogram definitions.
"""

import argparse
import math
import re
import shutil
import subprocess
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit


CPP_HELPERS = r"""
#ifndef HJ_NANOAOD_PLOT_HELPERS
#define HJ_NANOAOD_PLOT_HELPERS
#include <ROOT/RVec.hxx>
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace hjplot {
using ROOT::VecOps::RVec;
using P4 = ROOT::Math::PtEtaPhiMVector;
using CartesianP4 = ROOT::Math::PxPyPzEVector;
struct Angles {
    double cos_Theta = std::numeric_limits<double>::quiet_NaN();
    double cos_theta = std::numeric_limits<double>::quiet_NaN();
};
struct Candidate {
    bool valid = false;
    double pt_over_mass = 0., photon_eta = 0., dr_min = 0., dr_max = 0.;
    double m_llgamma = 0., m_ll = 0., photon_pt = 0., pt_llgamma = 0.;
    double cos_Theta = std::numeric_limits<double>::quiet_NaN();
    double cos_theta = std::numeric_limits<double>::quiet_NaN();
    int n_jets = 0;
    double leading_jet_pt = -1.;
};
double delta_r(double eta1, double phi1, double eta2, double phi2) {
    const double dphi = std::remainder(phi1 - phi2, 2. * std::acos(-1.));
    return std::hypot(eta1 - eta2, dphi);
}
bool finite_p4(double pt, double eta, double phi, double mass) {
    return std::isfinite(pt) && std::isfinite(eta) && std::isfinite(phi)
        && std::isfinite(mass) && pt > 0. && mass >= 0.;
}
double direction_cosine(const ROOT::Math::XYZVector& first,
                        const ROOT::Math::XYZVector& second) {
    const double norm = std::sqrt(first.Mag2()) * std::sqrt(second.Mag2());
    if (!std::isfinite(norm) || norm <= 0.) return std::numeric_limits<double>::quiet_NaN();
    const double value = first.Dot(second) / norm;
    if (!std::isfinite(value)) return std::numeric_limits<double>::quiet_NaN();
    return std::clamp(value, -1., 1.);
}
Angles decay_angles(const CartesianP4& z, const CartesianP4& photon,
                    const CartesianP4& negative_lepton) {
    Angles result;
    const auto higgs = z + photon;
    const auto has_rest_frame = [](const CartesianP4& p) {
        return std::isfinite(p.E()) && p.E() > 0. && std::isfinite(p.M2())
            && p.M2() > 0. && p.BoostToCM().Mag2() < 1.;
    };
    if (has_rest_frame(higgs)) {
        const ROOT::Math::Boost to_higgs(higgs.BoostToCM());
        // Helicity-axis prescription of arXiv:1112.1405: Z in the H rest
        // frame relative to the lab Zgamma flight direction, not lab Z eta.
        // Exactly zero lab momentum leaves this axis undefined (NaN).
        result.cos_Theta = direction_cosine(to_higgs(z).Vect(), higgs.Vect());
    }
    if (has_rest_frame(z)) {
        const ROOT::Math::Boost to_z(z.BoostToCM());
        // CMS HIG-25-010: negatively charged lepton versus photon in Z rest.
        result.cos_theta = direction_cosine(to_z(negative_lepton).Vect(), to_z(photon).Vect());
    }
    return result;
}
// Prefer the eligible pair closest to mZ, then its highest-pT photon.
template <typename Charge>
void consider_pairs(const RVec<float>& pt, const RVec<float>& eta,
                    const RVec<float>& phi, const RVec<float>& mass,
                    const RVec<Charge>& charge, double eta_limit,
                    const RVec<float>& photon_pt, const RVec<float>& photon_eta,
                    const RVec<float>& photon_phi, double& best_distance,
                    double& best_photon_pt, Candidate& result,
                    P4& best_l1, P4& best_l2, P4& best_photon, P4& best_lminus,
                    bool apply_selection) {
    for (std::size_t i = 0; i < pt.size(); ++i) {
        if (!finite_p4(pt[i], eta[i], phi[i], mass[i])) continue;
        if (apply_selection && (pt[i] < 10. || std::abs(eta[i]) >= eta_limit)) continue;
        for (std::size_t j = i + 1; j < pt.size(); ++j) {
            if (!finite_p4(pt[j], eta[j], phi[j], mass[j]) || charge[i] * charge[j] >= 0) continue;
            if (apply_selection && (pt[j] < 10. || std::abs(eta[j]) >= eta_limit
                || std::max(pt[i], pt[j]) < 20.)) continue;
            const P4 l1(pt[i], eta[i], phi[i], mass[i]);
            const P4 l2(pt[j], eta[j], phi[j], mass[j]);
            const auto z = l1 + l2;
            if (apply_selection && (z.M() < 50. || z.M() > 120.)) continue;
            const double distance = std::abs(z.M() - 91.1876);
            for (std::size_t k = 0; k < photon_pt.size(); ++k) {
                const double aeta = std::abs(photon_eta[k]);
                if (!finite_p4(photon_pt[k], photon_eta[k], photon_phi[k], 0.)) continue;
                if (apply_selection && (photon_pt[k] < 15. || aeta >= 2.5
                    || (aeta > 1.4442 && aeta < 1.566))) continue;
                const double dr1 = delta_r(eta[i], phi[i], photon_eta[k], photon_phi[k]);
                const double dr2 = delta_r(eta[j], phi[j], photon_eta[k], photon_phi[k]);
                if (apply_selection && std::min(dr1, dr2) <= 0.4) continue;
                const P4 photon(photon_pt[k], photon_eta[k], photon_phi[k], 0.);
                const auto higgs = z + photon;
                if (!std::isfinite(higgs.M()) || higgs.M() <= 0.) continue;
                if (apply_selection && (higgs.M() < 100. || higgs.M() > 180.)) continue;
                if (distance > best_distance
                    || (distance == best_distance && photon_pt[k] <= best_photon_pt)) continue;
                best_distance = distance;
                best_photon_pt = photon_pt[k];
                best_l1 = l1;
                best_l2 = l2;
                best_photon = photon;
                best_lminus = charge[i] < 0 ? l1 : l2;
                result.valid = true;
                result.pt_over_mass = higgs.Pt() / higgs.M();
                result.photon_eta = photon_eta[k];
                result.dr_min = std::min(dr1, dr2);
                result.dr_max = std::max(dr1, dr2);
                result.m_llgamma = higgs.M();
                result.m_ll = z.M();
                result.photon_pt = photon_pt[k];
                result.pt_llgamma = higgs.Pt();
            }
        }
    }
}
template <typename ElectronCharge, typename MuonCharge, typename JetId>
Candidate build(const RVec<float>& electron_pt, const RVec<float>& electron_eta,
                const RVec<float>& electron_phi, const RVec<float>& electron_mass,
                const RVec<ElectronCharge>& electron_charge,
                const RVec<float>& muon_pt, const RVec<float>& muon_eta,
                const RVec<float>& muon_phi, const RVec<float>& muon_mass,
                const RVec<MuonCharge>& muon_charge,
                const RVec<float>& photon_pt, const RVec<float>& photon_eta,
                const RVec<float>& photon_phi,
                const RVec<float>& jet_pt, const RVec<float>& jet_eta,
                const RVec<float>& jet_phi, const RVec<JetId>& jet_id,
                bool apply_selection) {
    Candidate result;
    double best_distance = std::numeric_limits<double>::infinity();
    double best_photon_pt = -1.;
    P4 l1, l2, photon, lminus;
    consider_pairs(electron_pt, electron_eta, electron_phi, electron_mass,
                   electron_charge, 2.5, photon_pt, photon_eta, photon_phi,
                   best_distance, best_photon_pt, result, l1, l2, photon, lminus, apply_selection);
    consider_pairs(muon_pt, muon_eta, muon_phi, muon_mass,
                   muon_charge, 2.4, photon_pt, photon_eta, photon_phi,
                   best_distance, best_photon_pt, result, l1, l2, photon, lminus, apply_selection);
    if (!result.valid) return result;
    const auto angles = decay_angles(CartesianP4(l1 + l2), CartesianP4(photon), CartesianP4(lminus));
    result.cos_Theta = angles.cos_Theta;
    result.cos_theta = angles.cos_theta;
    for (std::size_t j = 0; j < jet_pt.size(); ++j) {
        if (!finite_p4(jet_pt[j], jet_eta[j], jet_phi[j], 0.)) continue;
        if (apply_selection && (jet_pt[j] <= 30. || std::abs(jet_eta[j]) >= 4.7
            || (static_cast<int>(jet_id[j]) & 2) == 0)) continue;
        if (apply_selection && (delta_r(jet_eta[j], jet_phi[j], l1.Eta(), l1.Phi()) <= 0.4
            || delta_r(jet_eta[j], jet_phi[j], l2.Eta(), l2.Phi()) <= 0.4
            || delta_r(jet_eta[j], jet_phi[j], photon.Eta(), photon.Phi()) <= 0.4)) continue;
        ++result.n_jets;
        result.leading_jet_pt = std::max(result.leading_jet_pt, double(jet_pt[j]));
    }
    return result;
}
} // namespace hjplot
#endif
"""

CENTRAL_DIRECTORY = (
    "/eos/cms/store/mc/Run3Summer22NanoAODv12/"
    "GluGluHtoZG_Zto2L_M-125_TuneCP5_13p6TeV_powheg-pythia8/"
    "NANOAODSIM/130X_mcRun3_2022_realistic_v5-v2"
)

REQUIRED_BRANCHES = (
    "Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass", "Electron_charge",
    "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge",
    "Photon_pt", "Photon_eta", "Photon_phi",
    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "genWeight",
)
# column, panel title, x-axis label, bin count, lower edge, upper edge
PANELS = (
    ("pt_over_mass", "Higgs transverse recoil", "p_{T}(#it{l}#it{l}#gamma)/m_{#it{l}#it{l}#gamma}", 50, 0., 2.),
    ("photon_eta", "Photon pseudorapidity", "#eta(#gamma)", 50, -2.5, 2.5),
    ("dr_min", "Nearest photon-lepton separation", "#DeltaR_{min}(#gamma,#it{l})", 50, 0., 6.),
    ("dr_max", "Farthest photon-lepton separation", "#DeltaR_{max}(#gamma,#it{l})", 50, 0., 6.),
    ("n_jets", "Jet multiplicity", "N_{jet}", 10, -0.5, 9.5),
    ("leading_jet_pt", "Leading-jet transverse momentum", "p_{T}(j_{1}) [GeV]", 50, 0., 300.),
    ("m_llgamma", "Dilepton-photon invariant mass", "m_{#it{l}#it{l}#gamma} [GeV]", 50, 100., 180.),
    ("cos_Theta", "Z boson production angle", "cos#Theta", 50, -1., 1.),
    ("cos_theta", "Lepton production polar angle", "cos#theta", 50, -1., 1.),
    ("pt_llgamma", "Dilepton-photon transverse momentum", "p_{T}(#it{l}#it{l}#gamma) [GeV]", 50, 0., 300.),
)


def panels_for_selection(apply_selection):
    # Without mass windows, show a broader range instead of folding most events
    # into the edges of the selected mass windows. All other binning is unchanged.
    if apply_selection:
        return PANELS
    return tuple(
        (*panel[:3], 60, 0., 300.) if panel[0] == "m_llgamma" else panel
        for panel in PANELS
    )


def eos_location(directory):
    """Return EOS endpoint, absolute directory, and optional mounted path."""
    if directory.startswith("root://"):
        parsed = urlsplit(directory)
        if not parsed.netloc or parsed.query or parsed.fragment or parsed.username:
            raise ValueError("Supply a directory URL without credentials, query, or fragment.")
        endpoint = "root://" + parsed.netloc
        eos_path = "/" + parsed.path.lstrip("/")
    else:
        eos_path = str(PurePosixPath(directory))
        if eos_path.startswith("/eos/user/"):
            endpoint = "root://eosuser.cern.ch"
        elif eos_path.startswith("/eos/cms/"):
            endpoint = "root://eoscms.cern.ch"
        else:
            raise ValueError("Use an absolute /eos/user/... or /eos/cms/... directory, or a root:// EOS URL.")
    eos_path = eos_path.rstrip("/")
    if not eos_path.startswith("/eos/") or ".." in PurePosixPath(eos_path).parts:
        raise ValueError("The input must name a directory under /eos/.")
    mounted = None if directory.startswith("root://") else Path(eos_path)
    return endpoint, eos_path, mounted


def list_nanoaod_files(directory, maxdepth=1):
    """Discover files like find -maxdepth N; never recurse beyond N=2."""
    if maxdepth not in (1, 2):
        raise ValueError("File discovery supports maxdepth 1 or 2 only.")
    endpoint, eos_path, mounted = eos_location(directory)
    if mounted is not None and mounted.is_dir():
        folders = [mounted]
        if maxdepth == 2:
            folders.extend(p for p in mounted.iterdir() if p.is_dir() and not p.is_symlink())
        paths = [str(p) for folder in folders for p in folder.iterdir()
                 if p.is_file() and p.suffix == ".root"]
    else:
        if not shutil.which("xrdfs"):
            raise RuntimeError("EOS is not mounted here and xrdfs is unavailable. Run in a CMS/XRootD environment.")
        def listing(folder):
            result = subprocess.run(["xrdfs", endpoint, "ls", folder],
                                    check=True, text=True, capture_output=True)
            return [line if line.startswith("/") else folder + "/" + line
                    for line in map(str.strip, result.stdout.splitlines()) if line]

        top_entries = listing(eos_path)
        paths = [p for p in top_entries if p.endswith(".root")]
        if maxdepth == 2:
            for entry in top_entries:
                if entry.endswith(".root"):
                    continue
                result = subprocess.run(["xrdfs", endpoint, "stat", entry],
                                        check=True, text=True, capture_output=True)
                if "IsDir" in result.stdout:
                    paths.extend(p for p in listing(entry) if p.endswith(".root"))
    paths = sorted(set(paths))
    if not paths:
        raise ValueError("No .root files found within maxdepth {} in {}".format(maxdepth, eos_path))
    return [endpoint + "//" + p.lstrip("/") for p in paths], PurePosixPath(eos_path).name


def validate_files(root, files):
    """Check every file before plotting; never silently discard a broken file."""
    entries = 0
    for index, filename in enumerate(files, start=1):
        source = root.TFile.Open(filename)
        if not source or source.IsZombie():
            if source:
                source.Close()
            raise RuntimeError("Cannot open " + filename)
        try:
            tree = source.Get("Events")
            if not tree or not tree.InheritsFrom("TTree"):
                raise ValueError("Not a NanoAOD Events tree: " + filename)
            missing = [branch for branch in REQUIRED_BRANCHES if not tree.GetBranch(branch)]
            if missing:
                raise ValueError("Missing NanoAOD branches in {}: {}".format(filename, ", ".join(missing)))
            entries += int(tree.GetEntries())
        finally:
            source.Close()
        if index % 100 == 0 or index == len(files):
            print("Validated {}/{} files".format(index, len(files)), flush=True)
    if entries == 0:
        raise ValueError("The input NanoAOD files contain no events.")
    return entries


def book_histograms(root, files, apply_selection=True, sample_name="sample"):
    if not root.gInterpreter.Declare(CPP_HELPERS):
        raise RuntimeError("ROOT could not compile the candidate builder.")
    inputs = root.std.vector("string")()
    for filename in files:
        inputs.push_back(filename)
    frame = root.RDataFrame("Events", inputs)
    frame = frame.Define("plot_weight", "static_cast<double>(genWeight)")
    frame = frame.Filter("std::isfinite(plot_weight)", "Finite generator weight")
    expression = "hjplot::build(" + ", ".join(REQUIRED_BRANCHES[:-1]) + ", " + str(apply_selection).lower() + ")"
    selected = frame.Define("hj_candidate", expression).Filter("hj_candidate.valid", "Reconstructed llgamma candidate")
    panels = panels_for_selection(apply_selection)
    for column, *_ in panels:
        selected = selected.Define(column, "hj_candidate." + column)
    with_jet = selected.Filter("n_jets > 0", "At least one jet")
    actions = []
    for column, title, x_label, bins, low, high in panels:
        source = with_jet if column == "leading_jet_pt" else selected
        if column in ("cos_Theta", "cos_theta"):
            # An undefined rest-frame axis must not remove candidates from the
            # other observables, or be folded into a physical endpoint bin.
            source = source.Filter("std::isfinite(" + column + ")", "Defined " + column)
        actions.append(source.Histo1D((sample_name + "_" + column, "", bins, low, high), column, "plot_weight"))
    counts = (frame.Count(), selected.Count(), with_jet.Count(), selected.Sum("plot_weight"))
    # All actions are booked before this access: one loop fills all panels.
    counts[0].GetValue()
    histograms = []
    for action in actions:
        hist = action.GetValue().Clone()
        hist.SetDirectory(0)
        histograms.append(hist)
    return histograms, tuple(action.GetValue() for action in counts)


def normalize_histogram(hist):
    """Fold under/overflow, preserve sumw2 errors, and retain negative weights."""
    bins, entries = hist.GetNbinsX(), hist.GetEntries()
    for outside, inside in ((0, 1), (bins + 1, bins)):
        error = math.hypot(hist.GetBinError(outside), hist.GetBinError(inside))
        hist.SetBinContent(inside, hist.GetBinContent(inside) + hist.GetBinContent(outside))
        hist.SetBinError(inside, error)
        hist.SetBinContent(outside, 0.)
        hist.SetBinError(outside, 0.)
    hist.SetEntries(entries)
    integral = hist.Integral(1, bins)
    if entries > 0 and (not math.isfinite(integral) or integral <= 0.):
        raise ValueError("Non-positive sum of generator weights for " + hist.GetName())
    if integral > 0.:
        hist.Scale(1. / integral)


def make_ratio(root, private_hist, central_hist):
    """Private/central with independent weighted errors; omit zero denominators."""
    graph = root.TGraphErrors()
    graph.SetName(private_hist.GetName() + "_over_central")
    graph.SetLineColor(root.kBlack)
    graph.SetMarkerColor(root.kBlack)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(0.45)
    for index in range(1, central_hist.GetNbinsX() + 1):
        central, private = central_hist.GetBinContent(index), private_hist.GetBinContent(index)
        if central == 0.:
            continue
        value = private / central
        error = math.hypot(private_hist.GetBinError(index) / central,
                           private * central_hist.GetBinError(index) / (central * central))
        if not math.isfinite(value) or not math.isfinite(error):
            continue
        point = graph.GetN()
        graph.SetPoint(point, central_hist.GetBinCenter(index), value)
        graph.SetPointError(point, 0., error)
    return graph


def draw_y_title(root, pad, text, name):
    """Top-align both vertical titles at one x anchor and pixel font size."""
    pad.cd()
    title = root.TLatex(0.065, 1. - pad.GetTopMargin(), text)
    title.SetName(name)
    title.SetNDC(True)
    title.SetTextFont(43)
    title.SetTextSize(20.)
    title.SetTextAngle(90.)
    # Right alignment before the 90-degree rotation anchors the text's top end.
    title.SetTextAlign(32)
    title.Draw()
    return title


def draw_plots(root, central_histograms, private_histograms, label, output, apply_selection=True):
    """Reference-style overlays and ratio pads, without global title/cut notes."""
    root.gStyle.SetOptStat(0)
    root.gStyle.SetTextFont(42)
    root.gStyle.SetLegendFont(42)
    root.gStyle.SetCanvasColor(0)
    root.gStyle.SetPadColor(0)
    root.gStyle.SetEndErrorSize(0)
    canvas = root.TCanvas("hj_comparison_canvas", "", 3000, 1400)
    canvas.Divide(5, 2, 0.001, 0.001)
    keep = []
    colors = (root.TColor.GetColor("#D62728"), root.TColor.GetColor("#5B88CF"))
    for index, (central, private, panel) in enumerate(
            zip(central_histograms, private_histograms, panels_for_selection(apply_selection)), start=1):
        canvas.cd(index)
        top = root.TPad("hj_top_" + str(index), "", 0., 0.30, 1., 1.)
        bottom = root.TPad("hj_ratio_" + str(index), "", 0., 0., 1., 0.30)
        for pad in (top, bottom):
            pad.SetLeftMargin(0.20)
            pad.SetRightMargin(0.05)
            pad.SetTicks(1, 1)
        top.SetTopMargin(0.08)
        top.SetBottomMargin(0.025)
        bottom.SetTopMargin(0.025)
        bottom.SetBottomMargin(0.36)
        top.Draw()
        bottom.Draw()
        keep.extend((top, bottom))
        top.cd()
        for hist, color in zip((central, private), colors):
            normalize_histogram(hist)
            hist.SetTitle("")
            hist.SetLineColor(color)
            hist.SetMarkerColor(color)
            hist.SetLineWidth(3)
            hist.SetLineStyle(1)
            hist.SetFillStyle(0)
            for axis in (hist.GetXaxis(), hist.GetYaxis()):
                axis.SetTitleFont(42)
                axis.SetLabelFont(42)
        values = [(h.GetBinContent(b), h.GetBinError(b))
                  for h in (central, private) for b in range(1, h.GetNbinsX() + 1)]
        central.SetMaximum(max(0.05, max(v + e for v, e in values) * 1.35))
        central.SetMinimum(min(0., min(v - e for v, e in values) * 1.15))
        width = (panel[5] - panel[4]) / panel[3]
        unit = " GeV" if panel[2].endswith("[GeV]") else ""
        central.GetYaxis().SetTitle("")
        central.GetYaxis().SetLabelSize(0.049)
        central.GetYaxis().SetNdivisions(508)
        central.GetXaxis().SetLabelSize(0.)
        central.GetXaxis().SetTitleSize(0.)
        central.Draw("HIST")
        private.Draw("HIST SAME")
        central.Draw("E SAME")
        private.Draw("E SAME")
        keep.append(draw_y_title(root, top, "A.U. / {:g}{}".format(width, unit),
                                 "hj_y_title_top_" + str(index)))
        legend_left = 0.32
        header = root.TLegend(legend_left, 0.825, 0.95, 0.89)
        # Dedicated non-overlapping columns leave room for the longer labels
        # and multi-million-event counts; text no longer crosses column borders.
        names = root.TLegend(legend_left, 0.695, 0.71, 0.825)
        counts = root.TLegend(0.74, 0.695, 0.95, 0.825)
        for legend in (header, names, counts):
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextFont(42)
            legend.SetTextSize(0.040)
            legend.SetMargin(0.)
        header.SetTextSize(0.048)
        header.AddEntry(root.nullptr, "{} ggF signal sample".format(label), "")
        names.SetMargin(0.11)
        names.AddEntry(central, " ggH_qme (central)", "l")
        names.AddEntry(private, " HJ (private)", "l")
        counts.AddEntry(root.nullptr, "N={:,}".format(int(central.GetEntries())), "")
        counts.AddEntry(root.nullptr, "N={:,}".format(int(private.GetEntries())), "")
        for legend in (header, names, counts):
            legend.Draw()
        cms = root.TLatex()
        cms.SetNDC(True)
        cms.SetTextFont(42)
        cms.SetTextSize(0.060)
        cms.DrawLatex(0.195, 0.935, "#bf{CMS}#kern[-0.4]{ }#scale[0.9]{#it{Private Work}}")
        energy = root.TLatex()
        energy.SetNDC(True)
        energy.SetTextFont(42)
        energy.SetTextAlign(31)
        energy.SetTextSize(0.052)
        energy.DrawLatex(0.95, 0.935, "13.6 TeV")
        keep.extend((header, names, counts, cms, energy))
        top.RedrawAxis()
        bottom.cd()
        ratio = make_ratio(root, private, central)
        axis_hist = central.Clone("hj_ratio_axes_" + str(index))
        axis_hist.SetDirectory(0)
        axis_hist.Reset()
        # Display limits only: retain every weighted ratio point in the graph.
        axis_hist.SetMinimum(0.)
        axis_hist.SetMaximum(2.)
        axis_hist.GetXaxis().SetTitle(panel[2])
        axis_hist.GetYaxis().SetTitle("")
        for axis in (axis_hist.GetXaxis(), axis_hist.GetYaxis()):
            axis.SetTitleSize(0.125)
            axis.SetLabelSize(0.108)
        axis_hist.GetXaxis().SetTitleOffset(0.92)
        axis_hist.GetXaxis().SetLabelOffset(0.015)
        axis_hist.GetXaxis().SetNdivisions(510)
        axis_hist.GetYaxis().SetNdivisions(505)
        if panel[0] == "n_jets":
            for b in range(1, axis_hist.GetNbinsX()):
                axis_hist.GetXaxis().SetBinLabel(b, str(b - 1))
            axis_hist.GetXaxis().SetBinLabel(axis_hist.GetNbinsX(), "#geq9")
        axis_hist.Draw("AXIS")
        line = root.TLine(panel[4], 1., panel[5], 1.)
        line.SetLineColor(root.kBlack)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        ratio.Draw("P SAME")
        keep.append(draw_y_title(root, bottom, "Private / Central",
                                 "hj_y_title_ratio_" + str(index)))
        keep.extend((axis_hist, line, ratio))
        bottom.RedrawAxis()
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(output))
    canvas.Close()
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError("ROOT did not create the PNG: " + str(output))


def main(apply_selection=True):
    description = __doc__
    if not apply_selection:
        description = description.replace("_selected.png", "_no_selection.png")
        description = description.replace("Usage: python3 tools/plot_nanoaod_kinematics.py",
                                          "Usage: python3 tools/plot_nanoaod_kinematics_no_selection.py")
        description += "\nThis version keeps OS ee/mumu pairing but removes all object cuts."
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("eos_directory", help="Private NanoAOD EOS directory (all top-level ROOT files)")
    args = parser.parse_args()
    try:
        import ROOT
    except ImportError:
        parser.exit(1, "PyROOT is unavailable. Activate a CMSSW environment, then run with its python3.\n")
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.SetDefaultSumw2(True)
    ROOT.EnableImplicitMT(2)
    try:
        private_files, directory_label = list_nanoaod_files(args.eos_directory)
        central_files, _ = list_nanoaod_files(CENTRAL_DIRECTORY, maxdepth=2)
        samples = {}
        for name, files, directory in (
                ("Central", central_files, CENTRAL_DIRECTORY),
                ("Private", private_files, args.eos_directory)):
            print("{}: ALL {} ROOT files in {}".format(name, len(files), directory), flush=True)
            total_events = validate_files(ROOT, files)
            print("{}: processing {:,} events in one ROOT event loop...".format(name, total_events), flush=True)
            histograms, counts = book_histograms(ROOT, files, apply_selection, name.lower())
            if counts[1] == 0:
                raise ValueError("No reconstructed OS ee/mumu + photon candidates in " + name)
            samples[name] = histograms
            print("{}: {:,} candidates, {:,} with jets; sum(genWeight)={:.8g}".format(
                name, int(counts[1]), int(counts[2]), counts[3]), flush=True)
            if counts[0] != total_events:
                print("{}: excluded {:,} non-finite genWeight events.".format(
                    name, total_events - int(counts[0])), flush=True)
        label = re.sub(r"[^A-Za-z0-9_-]", "_", directory_label) or "NanoAOD"
        suffix = "selected" if apply_selection else "no_selection"
        output = Path(__file__).resolve().parents[1] / "plots" / (
            "HJ_" + label + "_central_vs_private_" + suffix + ".png")
        draw_plots(ROOT, samples["Central"], samples["Private"], label, output, apply_selection)
        print("Saved " + str(output))
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "Plotting failed: {}\n".format(error))


if __name__ == "__main__":
    main()
