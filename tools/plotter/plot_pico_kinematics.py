#!/usr/bin/env python3
"""Compare the 19 ggF BDT-input shapes in baseline-selected HJ/central picos.

Usage: python3 tools/plotter/plot_pico_kinematics.py 2022
Read ALL top-level ROOT files in the configured HJ/ggH pico directories.
Use tree/weight, not Events/genWeight; do not apply object or baseline cuts.
Standalone PyROOT plotting; three PNGs in the reference PDF's order.
"""

import argparse
import math
from pathlib import Path
from pathlib import PurePosixPath
import shutil
import subprocess
from urllib.parse import urlsplit

# column, panel title, x-axis label, bin count, lower edge, upper edge
PANELS = (
    ("pt_over_mass", "Higgs transverse recoil", "p_{T}(#it{l}#it{l}#gamma)/m_{#it{l}#it{l}#gamma}", 40, 0., 2.5),
    ("photon_idmva", "Photon identification MVA", "#gamma IDMVA", 45, 0., 1.),
    ("photon_rel_energy_err", "Photon relative energy uncertainty", "#sigma_{E}(#gamma)/E(#gamma)", 40, 0.01, 0.25),
    ("dr_min", "Nearest photon-lepton separation", "#DeltaR_{min}(#gamma,#it{l})", 35, 0., 3.5),
    ("dr_max", "Farthest photon-lepton separation", "#DeltaR_{max}(#gamma,#it{l})", 60, 0., 6.),
    ("cos_Theta", "Z boson production angle", "cos#Theta", 40, -1., 1.),
    ("cos_theta", "Lepton production polar angle", "cos#theta", 40, -1., 1.),
    ("phi", "Lepton production azimuthal angle", "#phi", 40, -3.2, 3.2),
    ("lepton1_eta", "Leading-lepton pseudorapidity", "#eta(#it{l}_{1})", 40, -2.6, 2.6),
    ("lepton2_eta", "Subleading-lepton pseudorapidity", "#eta(#it{l}_{2})", 40, -2.6, 2.6),
    ("photon_eta", "Photon pseudorapidity", "#eta(#gamma)", 40, -2.6, 2.6),
    ("photon_mht_dphi", "Photon-missing-HT azimuthal separation", "#Delta#phi(H_{T}^{miss},#gamma)", 40, 0., 3.15),
    ("photon_jet_dr", "Photon-leading-jet separation", "#DeltaR(#gamma,j_{1})", 40, 0.4, 6.),
    ("leading_jet_eta", "Leading-jet pseudorapidity", "#eta(j_{1})", 50, -5., 5.),
    ("leading_jet_mass", "Leading-jet mass", "m(j_{1}) [GeV]", 50, 0., 40.),
    ("leading_jet_pt", "Leading-jet transverse momentum", "p_{T}(j_{1}) [GeV]", 40, 30., 150.),
    ("llgamma_jet_dphi", "Dilepton-photon-leading-jet azimuthal separation", "#Delta#phi(Z#gamma,j_{1})", 40, 0., 3.15),
    ("system_balance", "System transverse-momentum balance", "System balance", 40, 0., 1.),
    ("photon_zeppenfeld", "Photon-leading-jet eta separation", "|#eta(#gamma)-#eta(j_{1})|", 40, 0., 6.),
)
PLOT_GROUPS = (PANELS[:8], PANELS[8:16], PANELS[16:])


SAMPLE_DIRECTORIES = {
    "2022": {
        "private": "/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/HJ2022pico",
        "central": "/eos/user/j/junhyuk/ggF_MiNNLO_NanoAOD/pico/ggH2022pico",
    },
}

REQUIRED_BRANCHES = (
    "weight", "llphoton_pt", "llphoton_m", "llphoton_iph", "llphoton_ill",
    "llphoton_cosTheta", "llphoton_costheta", "photon_eta", "photon_phi",
    "ll_i1", "ll_i2", "ll_lepid", "el_eta", "el_phi", "mu_eta", "mu_phi",
    "njet", "jet_pt", "jet_isgood", "jet_eta", "jet_phi", "jet_m",
    "photon_pt", "photon_idmva", "photon_energyErr", "llphoton_psi",
    "el_pt", "mu_pt", "llphoton_phi", "photon_mht_dphi", "ht", "mht",
)

CPP_HELPERS = r"""
#ifndef HJ_PICO_PLOT_HELPERS
#define HJ_PICO_PLOT_HELPERS
#include <ROOT/RVec.hxx>
#include <algorithm>
#include <cmath>
#include <limits>
namespace hjpicoplot {
using ROOT::VecOps::RVec;
double missing() { return std::numeric_limits<double>::quiet_NaN(); }
struct Values {
    double pt_over_mass = missing(), photon_eta = missing();
    double dr_min = missing(), dr_max = missing(), n_jets = missing();
    double leading_jet_pt = missing(), m_llgamma = missing();
    double cos_Theta = missing(), cos_theta = missing(), pt_llgamma = missing();
    double photon_idmva = missing(), photon_rel_energy_err = missing(), phi = missing();
    double lepton1_eta = missing(), lepton2_eta = missing(), photon_mht_dphi = missing();
    double photon_jet_dr = missing(), leading_jet_eta = missing(), leading_jet_mass = missing();
    double llgamma_jet_dphi = missing(), system_balance = missing(), photon_zeppenfeld = missing();
};
template <typename T> bool has(const RVec<T>& values, int index) {
    return index >= 0 && static_cast<std::size_t>(index) < values.size();
}
template <typename T> double value(const RVec<T>& values, int index) {
    if (!has(values, index)) return missing();
    const double result = values[index];
    return std::isfinite(result) ? result : missing();
}
double angle(const RVec<float>& values) {
    const double result = value(values, 0);
    // Reject undefined/sentinel angles; allow float round-off at endpoints.
    if (!std::isfinite(result) || std::abs(result) > 1. + 1.e-6) return missing();
    return std::clamp(result, -1., 1.);
}
double delta_r(double eta1, double phi1, double eta2, double phi2) {
    return std::hypot(eta1 - eta2,
                      std::remainder(phi1 - phi2, 2. * std::acos(-1.)));
}
Values build(const RVec<float>& pt, const RVec<float>& mass,
             const RVec<int>& iph, const RVec<int>& ill,
             const RVec<float>& photon_eta, const RVec<float>& photon_phi,
             const RVec<int>& ll_i1, const RVec<int>& ll_i2,
             const RVec<int>& ll_lepid,
             const RVec<float>& el_eta, const RVec<float>& el_phi,
             const RVec<float>& mu_eta, const RVec<float>& mu_phi,
             int njet, const RVec<float>& jet_pt, const RVec<bool>& jet_isgood,
             const RVec<float>& cosTheta, const RVec<float>& costheta,
             const RVec<float>& photon_pt, const RVec<float>& photon_idmva,
             const RVec<float>& photon_energyErr, const RVec<float>& psi,
             const RVec<float>& el_pt, const RVec<float>& mu_pt,
             const RVec<float>& jet_eta, const RVec<float>& jet_phi,
             const RVec<float>& jet_m, const RVec<float>& llphoton_phi,
             const RVec<float>& photon_mht_dphi, double ht, double mht) {
    Values out;
    // Candidate zero is the nominal candidate used by the baseline producer.
    out.pt_llgamma = value(pt, 0);
    out.m_llgamma = value(mass, 0);
    if (out.pt_llgamma < 0.) out.pt_llgamma = missing();
    if (out.m_llgamma <= 0.) out.m_llgamma = missing();
    out.pt_over_mass = out.pt_llgamma / out.m_llgamma;
    out.cos_Theta = angle(cosTheta);
    out.cos_theta = angle(costheta);
    const double psi0 = value(psi, 0);
    if (std::abs(psi0) <= std::acos(-1.) + 1.e-6) out.phi = psi0;
    if (std::isfinite(ht) && std::isfinite(mht) && ht > 0. && mht >= 0.)
        out.system_balance = mht / ht;
    if (njet >= 0) out.n_jets = njet;
    int leading_index = -1;
    // Use the producer's stored good-jet flags, without new object cuts.
    if (njet > 0 && jet_pt.size() == jet_isgood.size()) {
        double leading = missing();
        bool invalid = false;
        for (std::size_t j = 0; j < jet_pt.size(); ++j) {
            if (!jet_isgood[j]) continue;
            if (!std::isfinite(jet_pt[j]) || jet_pt[j] <= 0.) { invalid = true; break; }
            if (!std::isfinite(leading) || jet_pt[j] > leading) {
                leading = jet_pt[j];
                leading_index = j;
            }
        }
        if (!invalid) out.leading_jet_pt = leading;
        else leading_index = -1;
    }
    if (leading_index >= 0) {
        out.leading_jet_eta = value(jet_eta, leading_index);
        out.leading_jet_mass = value(jet_m, leading_index);
        if (out.leading_jet_mass < 0.) out.leading_jet_mass = missing();
        out.llgamma_jet_dphi = std::abs(std::remainder(
            value(llphoton_phi, 0) - value(jet_phi, leading_index), 2. * std::acos(-1.)));
    }
    if (iph.empty()) return out;
    const int photon_index = iph[0];
    out.photon_eta = value(photon_eta, photon_index);
    const double gphi = value(photon_phi, photon_index);
    out.photon_idmva = value(photon_idmva, photon_index);
    if (std::abs(out.photon_idmva) > 1. + 1.e-6) out.photon_idmva = missing();
    const double energy = value(photon_pt, photon_index) * std::cosh(out.photon_eta);
    const double energy_error = value(photon_energyErr, photon_index);
    if (std::isfinite(energy) && energy > 0. && energy_error >= 0.)
        out.photon_rel_energy_err = energy_error / energy;
    const double mht_dphi = value(photon_mht_dphi, photon_index);
    if (mht_dphi >= 0. && mht_dphi <= std::acos(-1.) + 1.e-6)
        out.photon_mht_dphi = mht_dphi;
    if (leading_index >= 0) {
        out.photon_jet_dr = delta_r(out.photon_eta, gphi,
                                  out.leading_jet_eta, value(jet_phi, leading_index));
        out.photon_zeppenfeld = std::abs(out.photon_eta - out.leading_jet_eta);
    }
    if (ill.empty()) return out;
    const int pair_index = ill[0];
    if (!has(ll_i1, pair_index) || !has(ll_i2, pair_index)
        || !has(ll_lepid, pair_index)) return out;
    const int first = ll_i1[pair_index], second = ll_i2[pair_index];
    const RVec<float>* eta = nullptr;
    const RVec<float>* phi = nullptr;
    const RVec<float>* lepton_pt = nullptr;
    if (ll_lepid[pair_index] == 11) { eta = &el_eta; phi = &el_phi; lepton_pt = &el_pt; }
    else if (ll_lepid[pair_index] == 13) { eta = &mu_eta; phi = &mu_phi; lepton_pt = &mu_pt; }
    else return out;
    const double pt1 = value(*lepton_pt, first), pt2 = value(*lepton_pt, second);
    if (std::isfinite(pt1) && std::isfinite(pt2) && pt1 > 0. && pt2 > 0.) {
        out.lepton1_eta = value(*eta, pt1 >= pt2 ? first : second);
        out.lepton2_eta = value(*eta, pt1 >= pt2 ? second : first);
    }
    const double dr1 = delta_r(out.photon_eta, gphi, value(*eta, first), value(*phi, first));
    const double dr2 = delta_r(out.photon_eta, gphi, value(*eta, second), value(*phi, second));
    if (std::isfinite(dr1) && std::isfinite(dr2)) {
        out.dr_min = std::min(dr1, dr2);
        out.dr_max = std::max(dr1, dr2);
    }
    return out;
}
} // namespace hjpicoplot
#endif
"""


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
        raise ValueError("Non-positive sum of pico weights for " + hist.GetName())
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


def draw_plots(root, central_histograms, private_histograms, label, output, event_counts, panels=None):
    """Reference-style overlays and ratio pads, without global title/cut notes."""
    if panels is None:
        panels = PANELS
    if not panels or len(central_histograms) != len(panels) or len(private_histograms) != len(panels):
        raise ValueError("Each panel needs exactly one central and one private histogram.")
    if any(len(event_counts[sample]) != len(panels) for sample in ("Central", "Private")):
        raise ValueError("Each panel needs a central and private event count.")
    root.gStyle.SetOptStat(0)
    root.gStyle.SetTextFont(42)
    root.gStyle.SetLegendFont(42)
    root.gStyle.SetCanvasColor(0)
    root.gStyle.SetPadColor(0)
    root.gStyle.SetEndErrorSize(0)
    columns = min(4, len(panels))
    rows = math.ceil(len(panels) / columns)
    # The PDF's individual canvases are wider than tall (approximately 3:2).
    canvas = root.TCanvas("hj_comparison_canvas", "", 900 * columns, 600 * rows)
    canvas.Divide(columns, rows, 0.001, 0.001)
    keep = []
    colors = (root.TColor.GetColor("#D62728"), root.TColor.GetColor("#5B88CF"))
    for index, (central, private, panel) in enumerate(
            zip(central_histograms, private_histograms, panels), start=1):
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
        # Two aligned rows at opposite top corners, without a sample heading.
        names = root.TLegend(0.24, 0.76, 0.65, 0.88)
        counts = root.TLegend(0.65, 0.76, 0.895, 0.88)
        for legend in (names, counts):
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextFont(42)
            legend.SetTextSize(0.040)
            legend.SetMargin(0.)
        names.SetMargin(0.11)
        names.AddEntry(central, " ggH_qme (central)", "l")
        names.AddEntry(private, " HJ (private)", "l")
        counts.SetTextAlign(32)
        counts.AddEntry(root.nullptr, "N={:,}".format(event_counts["Central"][index - 1]), "")
        counts.AddEntry(root.nullptr, "N={:,}".format(event_counts["Private"][index - 1]), "")
        for legend in (names, counts):
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
        keep.extend((names, counts, cms, energy))
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


def eos_location(directory):
    """Return the EOS endpoint, absolute path and optional mounted directory."""
    if directory.startswith("root://"):
        parsed = urlsplit(directory)
        if not parsed.netloc or parsed.query or parsed.fragment or parsed.username:
            raise ValueError("Supply an EOS directory URL without credentials, query or fragment.")
        endpoint = "root://" + parsed.netloc
        eos_path = "/" + parsed.path.lstrip("/")
        mounted = None
    else:
        eos_path = str(PurePosixPath(directory))
        if eos_path.startswith("/eos/user/"):
            endpoint = "root://eosuser.cern.ch"
        elif eos_path.startswith("/eos/cms/"):
            endpoint = "root://eoscms.cern.ch"
        else:
            raise ValueError("Use an absolute /eos/user/... or /eos/cms/... directory, or a root:// EOS URL.")
        mounted = Path(eos_path)
    eos_path = eos_path.rstrip("/")
    if not eos_path.startswith("/eos/") or ".." in PurePosixPath(eos_path).parts:
        raise ValueError("The input must name a directory under /eos/.")
    return endpoint, eos_path, mounted


def list_pico_files(directory):
    """Include every top-level ROOT file, but never recurse into subdirectories."""
    endpoint, eos_path, mounted = eos_location(directory)
    if mounted is not None and mounted.is_dir():
        paths = [str(path) for path in mounted.iterdir()
                 if path.is_file() and path.suffix == ".root"]
    else:
        if not shutil.which("xrdfs"):
            raise RuntimeError("EOS is not mounted and xrdfs is unavailable. Run in a CMS/XRootD environment.")
        result = subprocess.run(["xrdfs", endpoint, "ls", eos_path],
                                check=True, text=True, capture_output=True)
        paths = [line if line.startswith("/") else eos_path + "/" + line
                 for line in map(str.strip, result.stdout.splitlines()) if line.endswith(".root")]
    paths = sorted(set(paths))
    if not paths:
        raise ValueError("No top-level .root files found in " + eos_path)
    return [endpoint + "//" + path.lstrip("/") for path in paths], PurePosixPath(eos_path).name


def configured_year(value):
    if value not in SAMPLE_DIRECTORIES:
        raise argparse.ArgumentTypeError(
            "Pico locations for year {!r} are not configured. Configured: {}.".format(
                value, ", ".join(sorted(SAMPLE_DIRECTORIES))))
    return value


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("year", type=configured_year, metavar="YEAR")
    return parser, parser.parse_args(argv)


def validate_files(root, files):
    """Validate every file, including zero-entry baseline picos; never skip files."""
    entries = 0
    for index, filename in enumerate(files, 1):
        source = root.TFile.Open(filename, "READ")
        if not source or source.IsZombie():
            if source:
                source.Close()
            raise RuntimeError("Cannot open " + filename)
        try:
            tree = source.Get("tree")
            if not tree or not tree.InheritsFrom("TTree"):
                raise ValueError("Missing pico tree 'tree': " + filename)
            missing = [name for name in REQUIRED_BRANCHES if not tree.GetBranch(name)]
            if missing:
                raise ValueError("Missing pico branches in {}: {}".format(filename, ", ".join(missing)))
            entries += int(tree.GetEntries())
        finally:
            source.Close()
        if index % 100 == 0 or index == len(files):
            print("Validated {}/{} pico files".format(index, len(files)), flush=True)
    if entries == 0:
        raise ValueError("The pico sample contains no events.")
    return entries


def book_histograms(root, files, sample_name):
    if not root.gInterpreter.Declare(CPP_HELPERS):
        raise RuntimeError("ROOT could not compile the pico variable helpers.")
    root.TH1.SetDefaultSumw2(True)
    inputs = root.std.vector("string")()
    for filename in files:
        inputs.push_back(filename)
    frame = root.RDataFrame("tree", inputs)
    total = frame.Count()
    frame = frame.Define("hj_pico_weight", "static_cast<double>(weight)")
    finite_weight = frame.Filter("std::isfinite(hj_pico_weight)", "Finite stored pico weight")
    usable = finite_weight.Count()
    sumw = finite_weight.Sum("hj_pico_weight")
    expression = (
        "hjpicoplot::build(llphoton_pt, llphoton_m, llphoton_iph, llphoton_ill, "
        "photon_eta, photon_phi, ll_i1, ll_i2, ll_lepid, el_eta, el_phi, mu_eta, mu_phi, "
        "njet, jet_pt, jet_isgood, llphoton_cosTheta, llphoton_costheta, "
        "photon_pt, photon_idmva, photon_energyErr, llphoton_psi, el_pt, mu_pt, "
        "jet_eta, jet_phi, jet_m, llphoton_phi, photon_mht_dphi, ht, mht)"
    )
    values = finite_weight.Define("hj_pico_values", expression)
    actions, counts = [], []
    for column, _, _, bins, low, high in PANELS:
        scalar = "hj_pico_" + column
        defined = values.Define(scalar, "hj_pico_values." + column)
        # Bounds/finite-value checks only, not physics or baseline selection.
        defined = defined.Filter("std::isfinite(" + scalar + ")", "Defined " + column)
        counts.append(defined.Count())
        actions.append(defined.Histo1D(
            (sample_name + "_pico_" + column, "", bins, low, high), scalar, "hj_pico_weight"))
    # Book all actions before execution: one event loop per sample.
    total_events = int(total.GetValue())
    if int(usable.GetValue()) != total_events:
        raise ValueError("{} contains non-finite stored weights; refusing to silently drop events.".format(sample_name))
    histograms = []
    for action in actions:
        hist = action.GetValue().Clone()
        hist.SetDirectory(0)
        histograms.append(hist)
    return histograms, {
        "total": total_events, "sumw": float(sumw.GetValue()),
        "counts": tuple(int(count.GetValue()) for count in counts),
        "event_loops": int(frame.GetNRuns()),
    }


def main():
    parser, args = parse_arguments()
    try:
        import ROOT
    except ImportError:
        parser.exit(1, "PyROOT is unavailable. Activate CMSSW, then run with its python3.\n")
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.EnableImplicitMT(2)
    try:
        samples, event_counts = {}, {}
        for sample, key in (("Central", "central"), ("Private", "private")):
            directory = SAMPLE_DIRECTORIES[args.year][key]
            files, _ = list_pico_files(directory)
            print("{}: ALL {} pico ROOT files in {}".format(sample, len(files), directory), flush=True)
            entries = validate_files(ROOT, files)
            histograms, summary = book_histograms(ROOT, files, sample.lower())
            if summary["total"] != entries:
                raise RuntimeError("Input event count changed while reading " + sample)
            print("{}: {:,} baseline events; sum(weight)={:.8g}".format(
                sample, entries, summary["sumw"]), flush=True)
            for panel, count in zip(PANELS, summary["counts"]):
                print("  {}: N={:,}, undefined/no-object={:,}".format(panel[0], count, entries - count))
            samples[sample] = histograms
            event_counts[sample] = summary["counts"]
        start = 0
        for page, panels in enumerate(PLOT_GROUPS, start=1):
            stop = start + len(panels)
            output = Path(__file__).resolve().parents[2] / "plots" / (
                "HJ_" + args.year + "_pico_central_vs_private_" + str(page) + ".png")
            draw_plots(ROOT, samples["Central"][start:stop], samples["Private"][start:stop],
                       args.year, output,
                       event_counts={sample: counts[start:stop] for sample, counts in event_counts.items()},
                       panels=panels)
            print("Saved " + str(output))
            start = stop
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "Pico plotting failed: {}\n".format(error))


if __name__ == "__main__":
    main()
