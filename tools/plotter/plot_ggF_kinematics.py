#!/usr/bin/env python3
"""Compare six kinematic shapes between two or three baseline pico samples.

Required: --sample1-name, --sample1-location, --sample1-type {ggH,HJ},
          --sample2-name, --sample2-location, --sample2-type {ggH,HJ}.
Optional: --sample3-name, --sample3-location, --sample3-type (all three together),
          --output-prefix NAME, --max-depth N (default: 1).
ggH uses weight; HJ uses weight/w_nnlo, independently for each sample.
No extra selections. One unit-area 3x2 PNG; ratios = samples 2 and 3 / sample 1.
"""

import argparse
import math
import re
from pathlib import Path
from pathlib import PurePosixPath
import shutil
import subprocess
from urllib.parse import urlsplit

# column, panel title, x-axis label, bin count, lower edge, upper edge
PANELS = (
    ("pt_gamma", "Photon transverse momentum", "p_{T}(#gamma) [GeV]", 50, 15., 100.),
    ("pt_ll", "Dilepton transverse momentum", "p_{T}(#it{l}#it{l}) [GeV]", 50, 0., 200.),
    ("pt_llgamma", "Dilepton-photon transverse momentum", "p_{T}(#it{l}#it{l}#gamma) [GeV]", 50, 0., 200.),
    ("m_ll", "Dilepton invariant mass", "m_{#it{l}#it{l}} [GeV]", 35, 80., 100.),
    ("m_llgamma", "Dilepton-photon invariant mass", "m_{#it{l}#it{l}#gamma} [GeV]", 50, 100., 150.),
    ("npv", "Primary-vertex multiplicity", "N_{PV}", 70, 0., 70.),
)
REQUIRED_BRANCHES = (
    "weight", "photon_pt", "ll_pt", "ll_m", "llphoton_pt", "llphoton_m",
    "llphoton_iph", "llphoton_ill", "npv",
)

CPP_HELPERS = r"""
#ifndef HJ_PICO_KINEMATICS_HELPERS
#define HJ_PICO_KINEMATICS_HELPERS
#include <ROOT/RVec.hxx>
#include <cmath>
#include <limits>
namespace hjpicokin {
using ROOT::VecOps::RVec;
double missing() { return std::numeric_limits<double>::quiet_NaN(); }
struct Values {
    double pt_gamma = missing(), pt_ll = missing(), pt_llgamma = missing();
    double m_ll = missing(), m_llgamma = missing(), npv = missing();
};
double positive_value(const RVec<float>& values, int index, bool allow_zero) {
    if (index < 0 || static_cast<std::size_t>(index) >= values.size()) return missing();
    const double result = values[index];
    if (!std::isfinite(result) || result < 0. || (!allow_zero && result == 0.)) return missing();
    return result;
}
Values build(const RVec<float>& photon_pt, const RVec<float>& ll_pt,
             const RVec<float>& ll_m, const RVec<float>& llphoton_pt,
             const RVec<float>& llphoton_m, const RVec<int>& iph,
             const RVec<int>& ill, int npv) {
    Values out;
    out.pt_llgamma = positive_value(llphoton_pt, 0, true);
    out.m_llgamma = positive_value(llphoton_m, 0, false);
    if (!iph.empty()) out.pt_gamma = positive_value(photon_pt, iph[0], true);
    if (!ill.empty()) {
        out.pt_ll = positive_value(ll_pt, ill[0], true);
        out.m_ll = positive_value(ll_m, ill[0], false);
    }
    if (npv >= 0) out.npv = npv;
    return out;
}
} // namespace hjpicokin
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


def draw_y_title(root, pad, text, name, x=0.065):
    """Top-align both vertical titles at one x anchor and pixel font size."""
    pad.cd()
    title = root.TLatex(x, 1. - pad.GetTopMargin(), text)
    title.SetName(name)
    title.SetNDC(True)
    title.SetTextFont(43)
    title.SetTextSize(20.)
    title.SetTextAngle(90.)
    # Right alignment before the 90-degree rotation anchors the text's top end.
    title.SetTextAlign(32)
    title.Draw()
    return title


def draw_plots(root, central_histograms, private_histograms, label, output, event_counts,
               panels=None, columns=None, panel_size=(900, 600), reference_lines=None,
               legend_text_size=0.060, y_title_x=0.065,
               sample_names=("[Central] ggH_qme", "[Priavate] HJ (MiNNLO)"),
               ratio_title="Private / Central", third_histograms=None):
    """Reference-style overlays and ratio pads, without global title/cut notes."""
    if panels is None:
        panels = PANELS
    if not panels or len(central_histograms) != len(panels) or len(private_histograms) != len(panels):
        raise ValueError("Each panel needs exactly one central and one private histogram.")
    groups = [central_histograms, private_histograms]
    slots = ["Central", "Private"]
    if third_histograms is not None:
        if len(third_histograms) != len(panels):
            raise ValueError("Sample 3 needs one histogram per panel.")
        groups.append(third_histograms)
        slots.append("Third")
    if len(sample_names) != len(groups):
        raise ValueError("Each plotted sample needs its own legend name.")
    if any(sample not in event_counts or len(event_counts[sample]) != len(panels) for sample in slots):
        raise ValueError("Each panel needs event counts for every plotted sample.")
    root.gStyle.SetOptStat(0)
    root.gStyle.SetTextFont(42)
    root.gStyle.SetLegendFont(42)
    root.gStyle.SetCanvasColor(0)
    root.gStyle.SetPadColor(0)
    root.gStyle.SetEndErrorSize(0)
    if columns is None:
        columns = min(4, len(panels))
    if not isinstance(columns, int) or not 1 <= columns <= len(panels):
        raise ValueError("Invalid plot column count.")
    rows = math.ceil(len(panels) / columns)
    # BDT plots default to the PDF's wider 3:2 panels; other callers can opt
    # into the original taller kinematic panels without changing the styling.
    canvas = root.TCanvas("hj_comparison_canvas", "", panel_size[0] * columns, panel_size[1] * rows)
    canvas.Divide(columns, rows, 0.001, 0.001)
    keep = []
    colors = tuple(root.TColor.GetColor(color) for color in ("#D62728", "#5B88CF", "#2CA02C"))
    for index, panel in enumerate(panels, start=1):
        histograms = [group[index - 1] for group in groups]
        central = histograms[0]
        canvas.cd(index)
        categorical = bool(central.GetXaxis().GetLabels())
        vertical_labels = categorical and any(
            len(central.GetXaxis().GetBinLabel(b)) > 3
            for b in range(1, central.GetNbinsX() + 1))
        top = root.TPad("hj_top_" + str(index), "", 0., 0.30, 1., 1.)
        bottom = root.TPad("hj_ratio_" + str(index), "", 0., 0., 1., 0.30)
        for pad in (top, bottom):
            pad.SetLeftMargin(0.20)
            pad.SetRightMargin(0.05)
            pad.SetTicks(1, 1)
        top.SetTopMargin(0.08)
        top.SetBottomMargin(0.025)
        bottom.SetTopMargin(0.025)
        bottom.SetBottomMargin(0.55 if vertical_labels else 0.36)
        top.Draw()
        bottom.Draw()
        keep.extend((top, bottom))
        top.cd()
        for hist, color in zip(histograms, colors):
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
                  for h in histograms for b in range(1, h.GetNbinsX() + 1)]
        central.SetMaximum(max(0.05, max(v + e for v, e in values) * (1.50 if len(groups) == 3 else 1.35)))
        central.SetMinimum(min(0., min(v - e for v, e in values) * 1.15))
        width = (panel[5] - panel[4]) / panel[3]
        unit = " GeV" if panel[2].endswith("[GeV]") else ""
        central.GetYaxis().SetTitle("")
        central.GetYaxis().SetLabelSize(0.049)
        central.GetYaxis().SetNdivisions(508)
        central.GetXaxis().SetLabelSize(0.)
        central.GetXaxis().SetTitleSize(0.)
        central.Draw("HIST")
        for hist in histograms[1:]:
            hist.Draw("HIST SAME")
        if reference_lines and panel[0] in reference_lines:
            x = reference_lines[panel[0]]
            reference = root.TLine(x, central.GetMinimum(), x, central.GetMaximum())
            reference.SetLineColor(root.kBlack)
            reference.SetLineStyle(2)
            reference.SetLineWidth(2)
            reference.Draw()
            keep.append(reference)
        for hist in histograms:
            hist.Draw("E SAME")
        # Variable-width histograms show normalized bin populations. A single
        # fixed GeV width would mislabel them; fixed-width plotters keep theirs.
        y_title = ("A.U. / bin" if categorical or central.GetXaxis().GetXbins().GetSize()
                   else "A.U. / {:g}{}".format(width, unit))
        keep.append(draw_y_title(root, top, y_title,
                                 "hj_y_title_top_" + str(index), x=y_title_x))
        # Two aligned rows at opposite top corners, without a sample heading.
        legend_bottom = 0.64 if len(groups) == 3 else 0.72
        names = root.TLegend(0.24, legend_bottom, 0.65, 0.88)
        counts = root.TLegend(0.65, legend_bottom, 0.895, 0.88)
        for legend in (names, counts):
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextFont(42)
            legend.SetTextSize(legend_text_size)
            legend.SetMargin(0.)
        names.SetMargin(0.11)
        for hist, name in zip(histograms, sample_names):
            names.AddEntry(hist, " " + name, "l")
        counts.SetTextAlign(32)
        for slot in slots:
            counts.AddEntry(root.nullptr, "N={:,}".format(event_counts[slot][index - 1]), "")
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
        ratios = []
        for sample_index, hist in enumerate(histograms[1:], start=1):
            ratio = make_ratio(root, hist, central)
            if len(groups) == 3:
                ratio.SetLineColor(colors[sample_index])
                ratio.SetMarkerColor(colors[sample_index])
                ratio.SetMarkerStyle(19 + sample_index)
            ratios.append(ratio)
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
        axis_hist.GetXaxis().SetTitleOffset(2.4 if vertical_labels else 0.92)
        if categorical:
            axis_hist.GetXaxis().LabelsOption("v" if vertical_labels else "h")
            axis_hist.GetXaxis().SetLabelSize(0.085 if vertical_labels else 0.108)
        axis_hist.GetXaxis().SetLabelOffset(0.015)
        axis_hist.GetXaxis().SetNdivisions(510)
        axis_hist.GetYaxis().SetNdivisions(505)
        axis_hist.Draw("AXIS")
        line = root.TLine(panel[4], 1., panel[5], 1.)
        line.SetLineColor(root.kBlack)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        for ratio in ratios:
            ratio.Draw("P SAME")
        keep.append(draw_y_title(root, bottom, ratio_title,
                                 "hj_y_title_ratio_" + str(index), x=y_title_x))
        keep.extend((axis_hist, line, *ratios))
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


def list_pico_files(directory, max_depth=1):
    """Discover .root files at bounded depth (top-level files have depth 1).

    Local directories are read locally. EOS directories are read through
    XRootD, using a mounted directory for discovery when it is available.
    Never traverse symlinks or unbounded remote trees.
    """
    if not isinstance(max_depth, int) or max_depth < 1:
        raise ValueError("max_depth must be a positive integer.")
    is_eos = directory.startswith(("root://", "/eos/"))
    if is_eos:
        endpoint, base, mounted = eos_location(directory)
    else:
        mounted = Path(directory).expanduser().resolve()
        if not mounted.is_dir():
            raise ValueError("Input location is not a directory: " + directory)
        endpoint, base = None, str(mounted)
    paths = set()
    if mounted is not None and mounted.is_dir():
        pending = [(mounted, 1)]
        while pending:
            folder, depth = pending.pop()
            for path in folder.iterdir():
                if path.is_symlink():
                    continue
                if path.is_file() and path.suffix == ".root":
                    paths.add(str(path))
                elif depth < max_depth and path.is_dir():
                    pending.append((path, depth + 1))
    else:
        if not shutil.which("xrdfs"):
            raise RuntimeError("EOS is not mounted and xrdfs is unavailable. Activate CMS/XRootD.")
        pending, visited = [(base, 1)], set()
        while pending:
            folder, depth = pending.pop()
            if folder in visited:
                continue
            visited.add(folder)
            result = subprocess.run(["xrdfs", endpoint, "ls", "-l", folder],
                                    check=True, text=True, capture_output=True)
            for line in filter(None, map(str.strip, result.stdout.splitlines())):
                # xrdfs ls -l: permissions, date, time, size, full path.
                fields = line.split(None, 4)
                if len(fields) != 5 or fields[0][0] not in ("d", "-", "l"):
                    raise ValueError("Cannot parse xrdfs listing: " + line)
                kind, path = fields[0][0], fields[4].strip()
                if not path.startswith("/"):
                    path = folder + "/" + path
                candidate = PurePosixPath(path)
                if ".." in candidate.parts or str(candidate.parent) != folder:
                    raise ValueError("Unexpected path in xrdfs listing: " + path)
                if kind == "-" and candidate.suffix == ".root":
                    paths.add(str(candidate))
                elif kind == "d" and depth < max_depth:
                    pending.append((str(candidate), depth + 1))
    if not paths:
        raise ValueError("No .root files found within max-depth {} of {}".format(max_depth, base))
    files = sorted(paths)
    if endpoint:
        files = [endpoint + "//" + path.lstrip("/") for path in files]
    return files, PurePosixPath(base).name


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    for index in (1, 2, 3):
        required = index < 3
        parser.add_argument("--sample{}-name".format(index), required=required, type=sample_label,
                            help="Sample {} legend label (quote spaces).".format(index))
        parser.add_argument("--sample{}-location".format(index), required=required, type=sample_label,
                            help="Local/EOS pico directory, or root:// EOS directory URL.")
        parser.add_argument("--sample{}-type".format(index), required=required, choices=("ggH", "HJ"),
                            help="ggH: weight; HJ: weight/w_nnlo.")
    parser.add_argument("--output-prefix", type=output_prefix,
                        help="PNG filename prefix inside plots/; default uses all sample names.")
    parser.add_argument("--max-depth", type=positive_depth, default=1,
                        help="Input search depth for all samples: 1=direct files, 2=one subdirectory level.")
    args = parser.parse_args(argv)
    third = (args.sample3_name, args.sample3_location, args.sample3_type)
    if any(value is not None for value in third) and any(value is None for value in third):
        parser.error("Sample 3 requires --sample3-name, --sample3-location and --sample3-type together.")
    return parser, args


def sample_label(value):
    if not value.strip():
        raise argparse.ArgumentTypeError("Sample names and locations must not be empty.")
    return value.strip()


def positive_depth(value):
    try:
        result = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("max-depth must be a positive integer.")
    if result < 1:
        raise argparse.ArgumentTypeError("max-depth must be at least 1.")
    return result


def output_prefix(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", value):
        raise argparse.ArgumentTypeError("output-prefix must be a filename prefix using letters, "
                                         "digits, dots, hyphens or underscores; no directory paths.")
    return value


def output_paths(args):
    def slug(name):
        return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_") or "sample"
    names = [args.sample1_name, args.sample2_name]
    if args.sample3_name is not None:
        names.append(args.sample3_name)
    prefix = args.output_prefix or "ggF_kinematics_" + "_vs_".join(slug(name) for name in names)
    directory = Path(__file__).resolve().parents[2] / "plots"
    return [directory / (prefix + ".png")]




def book_histograms(root, files, sample_name, sample_type="ggH"):
    if sample_type not in ("ggH", "HJ"):
        raise ValueError("Sample type must be ggH or HJ.")
    if not root.gInterpreter.Declare(CPP_HELPERS):
        raise RuntimeError("ROOT could not compile the kinematic variable helpers.")
    root.TH1.SetDefaultSumw2(True)
    inputs = root.std.vector("string")()
    for filename in files:
        inputs.push_back(filename)
    frame = root.RDataFrame("tree", inputs)
    total = frame.Count()
    weight_expression = "static_cast<double>(weight)"
    if sample_type == "HJ":
        weight_expression = (
            "std::isfinite(static_cast<double>(w_nnlo)) && w_nnlo != 0. "
            "? static_cast<double>(weight) / static_cast<double>(w_nnlo) "
            ": std::numeric_limits<double>::quiet_NaN()")
    frame = frame.Define("hj_kin_weight", weight_expression)
    frame = frame.Define("hj_kin_weight2", "hj_kin_weight * hj_kin_weight")
    finite_weight = frame.Filter("std::isfinite(hj_kin_weight) && std::isfinite(hj_kin_weight2)",
                                 "Finite plotting weight and squared weight")
    usable = finite_weight.Count()
    sumw = finite_weight.Sum("hj_kin_weight")
    sumw2 = finite_weight.Sum("hj_kin_weight2")
    negative = finite_weight.Filter("hj_kin_weight < 0.").Count()
    expression = "hjpicokin::build(photon_pt, ll_pt, ll_m, llphoton_pt, llphoton_m, llphoton_iph, llphoton_ill, npv)"
    values = finite_weight.Define("hj_kin_values", expression)
    actions, counts = [], []
    for column, _, _, bins, low, high in PANELS:
        scalar = "hj_kin_" + column
        defined = values.Define(scalar, "hj_kin_values." + column)
        # Bounds/finite-value checks only, not physics or baseline selection.
        defined = defined.Filter("std::isfinite(" + scalar + ")", "Defined " + column)
        counts.append(defined.Count())
        actions.append(defined.Histo1D(
            (sample_name + "_kin_" + column, "", bins, low, high), scalar, "hj_kin_weight"))
    # Book all actions before execution: one event loop per sample.
    total_events = int(total.GetValue())
    if total_events == 0:
        raise ValueError(sample_name + ": no events.")
    if int(usable.GetValue()) != total_events:
        raise ValueError("{} contains non-finite stored weights, invalid HJ w_nnlo, or overflowing "
                         "plotting weights; refusing to silently drop events.".format(sample_name))
    histograms = []
    for action in actions:
        hist = action.GetValue().Clone()
        hist.SetDirectory(0)
        histograms.append(hist)
    return histograms, {
        "total": total_events, "sumw": float(sumw.GetValue()),
        "sumw2": float(sumw2.GetValue()), "negative_weights": int(negative.GetValue()),
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
        outputs = output_paths(args)
        for path in outputs:
            if path.exists():
                raise ValueError("Output already exists: {}. Use another --output-prefix or "
                                 "remove the old output explicitly.".format(path))
        samples, event_counts = {}, {}
        # Internal renderer slots retain their historical names; names, input
        # paths and ggH/HJ weighting are controlled solely by the CLI, not slots.
        sample_names = (args.sample1_name, args.sample2_name)
        if args.sample3_name is not None:
            sample_names += (args.sample3_name,)
        for index, slot in enumerate(("Central", "Private", "Third")[:len(sample_names)], 1):
            name = getattr(args, "sample{}_name".format(index))
            directory = getattr(args, "sample{}_location".format(index))
            kind = getattr(args, "sample{}_type".format(index))
            print("Sample {}: {} [{}]; discovering {} (max-depth {})".format(
                index, name, kind, directory, args.max_depth), flush=True)
            files, _ = list_pico_files(directory, max_depth=args.max_depth)
            print("  ALL {:,} ROOT files; plotting weight = {}".format(
                len(files), "weight/w_nnlo" if kind == "HJ" else "weight"), flush=True)
            print("  Processing events directly (no preliminary file-validation scan)...", flush=True)
            histograms, summary = book_histograms(ROOT, files, "sample" + str(index), sample_type=kind)
            entries = summary["total"]
            print("{}: {:,} baseline events; sum(plotting_weight)={:.8g}; sumw2={:.8g}; "
                  "negative-weight events={:,}".format(
                      name, entries, summary["sumw"], summary["sumw2"], summary["negative_weights"]), flush=True)
            for panel, count in zip(PANELS, summary["counts"]):
                print("  {}: N={:,}, undefined/no-object={:,}".format(panel[0], count, entries - count))
            samples[slot] = histograms
            event_counts[slot] = summary["counts"]
        output = outputs[0]
        draw_plots(ROOT, samples["Central"], samples["Private"], "", output,
                   event_counts=event_counts, panels=PANELS, columns=3, panel_size=(600, 700),
                   reference_lines={"m_ll": 91.2, "m_llgamma": 125.}, legend_text_size=0.040,
                   sample_names=sample_names,
                   ratio_title="Samples / 1" if "Third" in samples else "Sample 2 / 1",
                   third_histograms=samples.get("Third"))
        print("Saved " + str(output))
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "Kinematic plotting failed: {}\n".format(error))


if __name__ == "__main__":
    main()
