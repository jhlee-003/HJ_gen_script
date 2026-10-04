#!/usr/bin/env python3
"""Plot six reconstructed H -> Z gamma variables from one EOS directory.

Usage: python3 tools/plot_nanoaod_kinematics.py /eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/2022
Requires PyROOT (for example, the Python environment supplied by CMSSW).
Every top-level NanoAOD ROOT file is included, irrespective of cluster ID.
Output: plots/HJ_<directory-name>_kinematics.png in this repository.
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
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace hjplot {
using ROOT::VecOps::RVec;
using P4 = ROOT::Math::PtEtaPhiMVector;
struct Candidate {
    bool valid = false;
    double pt_over_mass = 0., photon_eta = 0., dr_min = 0., dr_max = 0.;
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
// Prefer the eligible pair closest to mZ, then its highest-pT photon.
template <typename Charge>
void consider_pairs(const RVec<float>& pt, const RVec<float>& eta,
                    const RVec<float>& phi, const RVec<float>& mass,
                    const RVec<Charge>& charge, double eta_limit,
                    const RVec<float>& photon_pt, const RVec<float>& photon_eta,
                    const RVec<float>& photon_phi, double& best_distance,
                    double& best_photon_pt, Candidate& result,
                    P4& best_l1, P4& best_l2, P4& best_photon) {
    for (std::size_t i = 0; i < pt.size(); ++i) {
        if (!finite_p4(pt[i], eta[i], phi[i], mass[i]) || pt[i] < 10.
            || std::abs(eta[i]) >= eta_limit) continue;
        for (std::size_t j = i + 1; j < pt.size(); ++j) {
            if (!finite_p4(pt[j], eta[j], phi[j], mass[j]) || pt[j] < 10.
                || std::abs(eta[j]) >= eta_limit || charge[i] * charge[j] >= 0
                || std::max(pt[i], pt[j]) < 20.) continue;
            const P4 l1(pt[i], eta[i], phi[i], mass[i]);
            const P4 l2(pt[j], eta[j], phi[j], mass[j]);
            const auto z = l1 + l2;
            if (z.M() < 50. || z.M() > 120.) continue;
            const double distance = std::abs(z.M() - 91.1876);
            for (std::size_t k = 0; k < photon_pt.size(); ++k) {
                const double aeta = std::abs(photon_eta[k]);
                if (!finite_p4(photon_pt[k], photon_eta[k], photon_phi[k], 0.)
                    || photon_pt[k] < 15. || aeta >= 2.5
                    || (aeta > 1.4442 && aeta < 1.566)) continue;
                const double dr1 = delta_r(eta[i], phi[i], photon_eta[k], photon_phi[k]);
                const double dr2 = delta_r(eta[j], phi[j], photon_eta[k], photon_phi[k]);
                if (std::min(dr1, dr2) <= 0.4) continue;
                const P4 photon(photon_pt[k], photon_eta[k], photon_phi[k], 0.);
                const auto higgs = z + photon;
                if (higgs.M() < 100. || higgs.M() > 180.) continue;
                if (distance > best_distance
                    || (distance == best_distance && photon_pt[k] <= best_photon_pt)) continue;
                best_distance = distance;
                best_photon_pt = photon_pt[k];
                best_l1 = l1;
                best_l2 = l2;
                best_photon = photon;
                result.valid = true;
                result.pt_over_mass = higgs.Pt() / higgs.M();
                result.photon_eta = photon_eta[k];
                result.dr_min = std::min(dr1, dr2);
                result.dr_max = std::max(dr1, dr2);
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
                const RVec<float>& jet_phi, const RVec<JetId>& jet_id) {
    Candidate result;
    double best_distance = std::numeric_limits<double>::infinity();
    double best_photon_pt = -1.;
    P4 l1, l2, photon;
    consider_pairs(electron_pt, electron_eta, electron_phi, electron_mass,
                   electron_charge, 2.5, photon_pt, photon_eta, photon_phi,
                   best_distance, best_photon_pt, result, l1, l2, photon);
    consider_pairs(muon_pt, muon_eta, muon_phi, muon_mass,
                   muon_charge, 2.4, photon_pt, photon_eta, photon_phi,
                   best_distance, best_photon_pt, result, l1, l2, photon);
    if (!result.valid) return result;
    for (std::size_t j = 0; j < jet_pt.size(); ++j) {
        if (!finite_p4(jet_pt[j], jet_eta[j], jet_phi[j], 0.) || jet_pt[j] <= 30.
            || std::abs(jet_eta[j]) >= 4.7 || (static_cast<int>(jet_id[j]) & 2) == 0) continue;
        if (delta_r(jet_eta[j], jet_phi[j], l1.Eta(), l1.Phi()) <= 0.4
            || delta_r(jet_eta[j], jet_phi[j], l2.Eta(), l2.Phi()) <= 0.4
            || delta_r(jet_eta[j], jet_phi[j], photon.Eta(), photon.Phi()) <= 0.4) continue;
        ++result.n_jets;
        result.leading_jet_pt = std::max(result.leading_jet_pt, double(jet_pt[j]));
    }
    return result;
}
} // namespace hjplot
"""

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
    ("n_jets", "Jet multiplicity", "N_{jet} (p_{T}>30 GeV, |#eta|<4.7)", 10, -0.5, 9.5),
    ("leading_jet_pt", "Leading-jet transverse momentum", "p_{T}(j_{1}) [GeV]", 50, 0., 300.),
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


def list_nanoaod_files(directory):
    endpoint, eos_path, mounted = eos_location(directory)
    if mounted is not None and mounted.is_dir():
        paths = [str(p) for p in mounted.iterdir() if p.is_file() and p.suffix == ".root"]
    else:
        if not shutil.which("xrdfs"):
            raise RuntimeError("EOS is not mounted here and xrdfs is unavailable. Run in a CMS/XRootD environment.")
        result = subprocess.run(
            ["xrdfs", endpoint, "ls", eos_path], check=True, text=True, capture_output=True
        )
        paths = [line.strip() for line in result.stdout.splitlines() if line.strip().endswith(".root")]
        paths = [p if p.startswith("/") else eos_path + "/" + p for p in paths]
    paths = sorted(set(paths))
    if not paths:
        raise ValueError("No top-level .root files found in " + eos_path)
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


def book_histograms(root, files):
    if not root.gInterpreter.Declare(CPP_HELPERS):
        raise RuntimeError("ROOT could not compile the candidate builder.")
    inputs = root.std.vector("string")()
    for filename in files:
        inputs.push_back(filename)
    frame = root.RDataFrame("Events", inputs)
    frame = frame.Define("plot_weight", "static_cast<double>(genWeight)")
    frame = frame.Filter("std::isfinite(plot_weight)", "Finite generator weight")
    expression = "hjplot::build(" + ", ".join(REQUIRED_BRANCHES[:-1]) + ")"
    selected = frame.Define("hj_candidate", expression).Filter("hj_candidate.valid", "Reconstructed llgamma candidate")
    for column, *_ in PANELS:
        selected = selected.Define(column, "hj_candidate." + column)
    with_jet = selected.Filter("n_jets > 0", "At least one cleaned jet")
    actions = []
    for column, title, x_label, bins, low, high in PANELS:
        source = with_jet if column == "leading_jet_pt" else selected
        actions.append(source.Histo1D(("h_" + column, "", bins, low, high), column, "plot_weight"))
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


def draw_plots(root, histograms, counts, n_files, total_events, label, output):
    root.gStyle.SetOptStat(0)
    root.gStyle.SetTextFont(42)
    root.gStyle.SetCanvasColor(0)
    root.gStyle.SetPadColor(0)
    canvas = root.TCanvas("hj_canvas", "HJ NanoAOD kinematics", 1800, 1200)
    area = root.TPad("hj_panels", "", 0., 0.105, 1., 0.895)
    area.Draw()
    area.Divide(3, 2, 0.007, 0.013)
    annotations = []
    for index, (hist, panel) in enumerate(zip(histograms, PANELS), start=1):
        pad = area.cd(index)
        pad.SetLeftMargin(0.17)
        pad.SetRightMargin(0.035)
        pad.SetTopMargin(0.13)
        pad.SetBottomMargin(0.16)
        pad.SetTicks(1, 1)
        normalize_histogram(hist)
        hist.SetLineColor(root.kAzure + 2)
        hist.SetLineWidth(2)
        hist.SetFillColorAlpha(root.kAzure - 9, 0.55)
        hist.SetMarkerColor(root.kAzure + 2)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(0.35)
        hist.GetXaxis().SetTitle(panel[2])
        hist.GetYaxis().SetTitle("Normalized weighted events / bin")
        for axis in (hist.GetXaxis(), hist.GetYaxis()):
            axis.SetTitleFont(42)
            axis.SetLabelFont(42)
            axis.SetTitleSize(0.047)
            axis.SetLabelSize(0.04)
        hist.GetYaxis().SetTitleOffset(1.65)
        hist.GetYaxis().SetNdivisions(505)
        if panel[0] == "n_jets":
            for i in range(1, hist.GetNbinsX()):
                hist.GetXaxis().SetBinLabel(i, str(i - 1))
            hist.GetXaxis().SetBinLabel(hist.GetNbinsX(), "#geq9")
        high = max(hist.GetBinContent(i) + hist.GetBinError(i) for i in range(1, hist.GetNbinsX() + 1))
        low = min(hist.GetBinContent(i) - hist.GetBinError(i) for i in range(1, hist.GetNbinsX() + 1))
        hist.SetMaximum(max(0.05, high * 1.25))
        hist.SetMinimum(min(0., low * 1.15))
        hist.Draw("HIST")
        hist.Draw("E1 SAME")
        title = root.TLatex()
        title.SetNDC(True)
        title.SetTextFont(62)
        title.SetTextSize(0.044)
        title.DrawLatex(0.17, 0.925, panel[1])
        annotations.append(title)
        if panel[0] == "leading_jet_pt":
            note = root.TLatex()
            note.SetNDC(True)
            note.SetTextSize(0.036)
            note.DrawLatex(0.55, 0.80, "N_{jet} #geq 1")
            annotations.append(note)
    canvas.cd()
    text = root.TLatex()
    text.SetNDC(True)
    text.SetTextFont(62)
    text.SetTextSize(0.030)
    text.DrawLatex(0.035, 0.960, "ggF H #rightarrow Z#gamma #rightarrow #it{l}^{+}#it{l}^{-}#gamma  |  13.6 TeV")
    text.SetTextFont(42)
    text.SetTextSize(0.022)
    text.DrawLatex(0.035, 0.920, "{}  |  {:,} files  |  {:,} events  |  {:,} selected".format(label, n_files, total_events, int(counts[1])))
    text.SetTextSize(0.018)
    text.DrawLatex(0.035, 0.073, "Reco OS ee/#mu#mu; p_{T}(#it{l}_{1,2}) #geq 20,10 GeV; 50 #leq m_{#it{l}#it{l}} #leq 120 GeV; 100 #leq m_{#it{l}#it{l}#gamma} #leq 180 GeV.")
    text.DrawLatex(0.035, 0.045, "p_{T}(#gamma) #geq 15 GeV; |#eta(#gamma)|<2.5, ECAL gap excluded; #DeltaR(#gamma,#it{l})>0.4. No lepton/photon ID, isolation or trigger cuts.")
    text.DrawLatex(0.035, 0.017, "genWeight-weighted shapes; edge bins include under/overflow. Cleaned tight-ID jets: p_{T}>30 GeV, |#eta|<4.7, #DeltaR(j,#it{l}/#gamma)>0.4.")
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(output))
    canvas.Close()
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError("ROOT did not create the PNG: " + str(output))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("eos_directory", help="EOS directory containing the NanoAOD .root files (no cluster filtering)")
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
        files, directory_label = list_nanoaod_files(args.eos_directory)
        print("Using ALL {} ROOT files in {}".format(len(files), args.eos_directory), flush=True)
        total_events = validate_files(ROOT, files)
        print("Processing {:,} events in one ROOT event loop...".format(total_events), flush=True)
        histograms, counts = book_histograms(ROOT, files)
        if counts[1] == 0:
            raise ValueError("No reconstructed ee/mumu + photon candidates pass the documented baseline selection.")
        label = re.sub(r"[^A-Za-z0-9_-]", "_", directory_label) or "NanoAOD"
        output = Path(__file__).resolve().parents[1] / "plots" / ("HJ_" + label + "_kinematics.png")
        draw_plots(ROOT, histograms, counts, len(files), total_events, label, output)
        print("Selected {:,} llgamma events, {:,} with at least one cleaned jet.".format(int(counts[1]), int(counts[2])))
        print("Selected sum(genWeight): {:.8g}".format(counts[3]))
        if counts[0] != total_events:
            print("Excluded {:,} events with non-finite genWeight.".format(total_events - int(counts[0])))
        print("Saved " + str(output))
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "Plotting failed: {}\n".format(error))


if __name__ == "__main__":
    main()
