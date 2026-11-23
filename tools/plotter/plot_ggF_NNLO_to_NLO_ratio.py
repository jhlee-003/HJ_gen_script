#!/usr/bin/env python3
"""Compare up to three named HJ/ggH Higgs-pT ratios with the Run 3 NNLO lookup.

Required: --sample1-name, --sample1-nnlo-location, --sample1-nlo-location.
Optional: the same three arguments for sample2 and/or sample3, all-or-none.
Use signed weight/w_nnlo for BOTH samples, then normalize each shape.
Select the last mc_id == 25 particle. One comparison PNG; no further cuts.
Standalone script; no dependency on the retired Higgs-pT/w_nnlo plotters.
"""

import argparse
from array import array
import math
import re
from pathlib import Path
from pathlib import PurePosixPath
import shutil
import subprocess
from urllib.parse import urlsplit

# Run 3 lookup pinned from richstu/nano2pico htozgamma_sequoia_v1:
# data/zgamma/GluGluHToZG_NNLO_reweight_run3.json (same as the old 2022 plotter).
# Do not mix Run 2 inputs with this table. No runtime network downloads.
PT_BIN_EDGES = (0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 70.,
                100., 150., 200., 300., 400., 500.)
PT_BIN_LABELS = tuple("{:g}-{:g}".format(low, high)
                      for low, high in zip(PT_BIN_EDGES, PT_BIN_EDGES[1:]))
LOOKUP_WEIGHTS = (
        1.2685261185511314,
        1.0590338886098563,
        0.9718308207358878,
        0.9518462182100818,
        0.9268171153336202,
        0.9251468016571035,
        0.9418127875150621,
        0.9557982359446595,
        0.9967997964299106,
        1.0323538167414477,
        1.0555095800469174,
        1.0745694746526047,
        1.1132983728193748,
        1.200345859778174,
        1.4516191583465505,
        1.6961873040033373,
)

CPP_HELPERS = r"""
#ifndef HJ_PICO_HIGGS_PT_HELPERS
#define HJ_PICO_HIGGS_PT_HELPERS
#include <ROOT/RVec.hxx>
#include <cmath>
#include <limits>
namespace hjpicohiggspt {
struct Values {
    int higgs_count = 0;
    double pt = std::numeric_limits<double>::quiet_NaN();
    bool malformed = false;
};
Values build(const ROOT::VecOps::RVec<int>& ids,
             const ROOT::VecOps::RVec<float>& pts) {
    Values out;
    if (ids.size() != pts.size()) {
        out.malformed = true;
        return out;
    }
    // Match NNLOCorrection's repeated hidx assignment exactly. Do not sort,
    // sum copies, select maximum pT, or fall back to an earlier valid copy.
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (ids[i] == 25) {
            ++out.higgs_count;
            out.pt = pts[i];
        }
    }
    return out;
}
} // namespace hjpicohiggspt
#endif
"""

def lookup_weights():
    return LOOKUP_WEIGHTS


REQUIRED_BRANCHES = ("mc_id", "mc_pt", "weight", "w_nnlo")


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


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for index in (1, 2, 3):
        required = index == 1
        parser.add_argument("--sample{}-name".format(index), required=required, type=nonempty,
                            help="Label for NNLO/NLO comparison {}.".format(index))
        parser.add_argument("--sample{}-nnlo-location".format(index), required=required, type=nonempty,
                            help="HJ baseline-pico directory (NNLO numerator).")
        parser.add_argument("--sample{}-nlo-location".format(index), required=required, type=nonempty,
                            help="ggH baseline-pico directory (NLO denominator).")
    args = parser.parse_args(argv)
    args.comparisons = []
    for index in (1, 2, 3):
        name = getattr(args, "sample{}_name".format(index))
        nnlo = getattr(args, "sample{}_nnlo_location".format(index))
        nlo = getattr(args, "sample{}_nlo_location".format(index))
        if any(value is not None for value in (name, nnlo, nlo)):
            if any(value is None for value in (name, nnlo, nlo)):
                parser.error("Sample {} needs all three: name, nnlo-location and nlo-location.".format(index))
            args.comparisons.append({"index": index, "name": name,
                                     "nnlo_location": nnlo, "nlo_location": nlo})
    return parser, args


def nonempty(value):
    if not value.strip():
        raise argparse.ArgumentTypeError("Names and locations must not be empty.")
    return value.strip()


def output_path(args):
    names = [re.sub(r"[^A-Za-z0-9]+", "_", comparison["name"]).strip("_")[:45]
             or "sample" + str(comparison["index"]) for comparison in args.comparisons]
    return Path(__file__).resolve().parents[2] / "plots" / (
        "ggF_NNLO_to_NLO_ratio_" + "_".join(names) + ".png")




def book_histograms(root, files, sample_name):
    """Fill one signed weight/w_nnlo histogram and all diagnostics in one pass."""
    if not root.gInterpreter.Declare(CPP_HELPERS):
        raise RuntimeError("ROOT could not compile the generator Higgs pT helper.")
    root.TH1.SetDefaultSumw2(True)
    inputs = root.std.vector("string")()
    for filename in files:
        inputs.push_back(filename)
    frame = root.RDataFrame("tree", inputs)
    total = frame.Count()
    weighted = frame.Define(
        "nnlo_plot_weight",
        "std::isfinite(static_cast<double>(weight)) && "
        "std::isfinite(static_cast<double>(w_nnlo)) && w_nnlo != 0. "
        "? static_cast<double>(weight) / static_cast<double>(w_nnlo) "
        ": std::numeric_limits<double>::quiet_NaN()")
    weighted = weighted.Define("nnlo_plot_weight2", "nnlo_plot_weight * nnlo_plot_weight")
    finite = weighted.Filter("std::isfinite(nnlo_plot_weight) && std::isfinite(nnlo_plot_weight2)")
    usable = finite.Count()
    values = finite.Define("nnlo_higgs", "hjpicohiggspt::build(mc_id, mc_pt)")
    malformed = values.Filter("nnlo_higgs.malformed").Count()
    well_formed = values.Filter("!nnlo_higgs.malformed")
    no_higgs = well_formed.Filter("nnlo_higgs.higgs_count == 0").Count()
    selected = well_formed.Filter("nnlo_higgs.higgs_count > 0").Define(
        "nnlo_higgs_pt", "nnlo_higgs.pt")
    selected_count = selected.Count()
    valid = selected.Filter("std::isfinite(nnlo_higgs_pt) && nnlo_higgs_pt >= 0.")
    valid_count = valid.Count()
    above_range = valid.Filter("nnlo_higgs_pt >= " + str(PT_BIN_EDGES[-1])).Count()
    in_range = valid.Filter("nnlo_higgs_pt < " + str(PT_BIN_EDGES[-1]))
    plotted_count = in_range.Count()
    action = in_range.Histo1D(
        (sample_name + "_weighted_higgs_pt", "", len(PT_BIN_LABELS), array("d", PT_BIN_EDGES)),
        "nnlo_higgs_pt", "nnlo_plot_weight")
    entries = int(total.GetValue())
    if entries == 0:
        raise ValueError(sample_name + ": no events.")
    if int(usable.GetValue()) != entries:
        raise ValueError(sample_name + ": invalid weight/w_nnlo (zero/non-finite divisor or non-finite weight).")
    if int(malformed.GetValue()):
        raise ValueError(sample_name + ": different mc_id/mc_pt vector lengths.")
    if int(selected_count.GetValue()) != int(valid_count.GetValue()):
        raise ValueError(sample_name + ": invalid last Higgs pT (negative or non-finite).")
    plotted = int(plotted_count.GetValue())
    if not plotted:
        raise ValueError(sample_name + ": no saved Higgs in the plotted pT range.")
    hist = action.GetValue().Clone()
    hist.SetDirectory(0)
    return [hist], {
        "total": entries, "counts": (plotted,),
        "no_higgs": int(no_higgs.GetValue()), "above_range": int(above_range.GetValue()),
        "sumw": float(hist.Integral()), "event_loops": int(frame.GetNRuns()),
    }


def normalized_fractions(sums, sumw2, sample_name):
    """Propagate each signed bin sum and its shared sample normalization."""
    if any(not math.isfinite(s) for s in sums):
        raise ValueError(sample_name + ": non-finite bin sums.")
    if any(not math.isfinite(q) or q < 0 for q in sumw2):
        raise ValueError(sample_name + ": invalid summed squared weights.")
    total, total_q = math.fsum(sums), math.fsum(sumw2)
    if not math.isfinite(total) or total <= 0:
        raise ValueError(sample_name + ": non-positive total plotting weight.")
    fractions, variances = [], []
    for s, q in zip(sums, sumw2):
        fraction = s / total
        # Var(S_i/S) includes Cov(S_i,S)=sum(w^2)_i; valid for signed weights.
        variance = ((1. - fraction) ** 2 * q
                    + fraction ** 2 * max(0., total_q - q)) / total ** 2
        fractions.append(fraction)
        variances.append(variance)
    return fractions, variances, total


def build_comparison(nlo_sums, nnlo_sums, nlo_sumw2, nnlo_sumw2, weights):
    bins = len(PT_BIN_LABELS)
    if not all(len(values) == bins for values in (
            nlo_sums, nnlo_sums, nlo_sumw2, nnlo_sumw2, weights)):
        raise ValueError("Bin sums, squared weights, and lookup must match all pT intervals.")
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError("Lookup corrections must be finite and positive.")
    fc, vc, sc = normalized_fractions(nlo_sums, nlo_sumw2, "NLO")
    fh, vh, sh = normalized_fractions(nnlo_sums, nnlo_sumw2, "NNLO")
    rows = []
    for i, w in enumerate(weights):
        ratio = error = closure = closure_error = None
        if fc[i] != 0.:
            ratio = fh[i] / fc[i]
            variance = vh[i] / fc[i] ** 2 + fh[i] ** 2 * vc[i] / fc[i] ** 4
            error = math.sqrt(max(0., variance))
            closure, closure_error = ratio / w, error / w
            if not all(math.isfinite(v) for v in (ratio, error, closure, closure_error)):
                raise ValueError("Non-finite weighted ratio in pT bin " + str(i + 1))
        rows.append({
            "low": PT_BIN_EDGES[i], "high": PT_BIN_EDGES[i + 1],
            "nlo_sumw": nlo_sums[i], "nnlo_sumw": nnlo_sums[i],
            "nlo_sumw2": nlo_sumw2[i], "nnlo_sumw2": nnlo_sumw2[i],
            "nlo_fraction": fc[i], "nnlo_fraction": fh[i],
            "ratio": ratio, "error": error, "lookup": w,
            "closure": closure, "closure_error": closure_error,
        })
    return rows, {
        "nlo_sumw": sc, "nnlo_sumw": sh,
        "undefined_bins": sum(row["ratio"] is None for row in rows),
    }


def make_graphs(root, rows, tag="comparison", x_offset=0.):
    observed = root.TGraphErrors()
    observed.SetName(tag + "_nnlo_over_nlo")
    closure = root.TGraphErrors()
    closure.SetName(tag + "_over_lookup")
    for index, row in enumerate(rows, 1):
        if row["ratio"] is None:
            continue
        point = observed.GetN()
        observed.SetPoint(point, float(index) + x_offset, row["ratio"])
        observed.SetPointError(point, 0., row["error"])
        closure.SetPoint(point, float(index) + x_offset, row["closure"])
        closure.SetPointError(point, 0., row["closure_error"])
    return observed, closure


def make_bars(root, rows, field, color, width, offset, limits):
    """Bars occupy a sub-interval of each category; undefined ratios stay gaps.

    Clip only the painted rectangle to the frame. Full values and errors stay
    unchanged in rows/error graphs, including points outside the lower range.
    """
    bars = []
    baseline = min(limits[1], max(limits[0], 0.))
    for index, row in enumerate(rows, 1):
        if row["ratio"] is None:
            continue
        value = min(limits[1], max(limits[0], row[field]))
        center = float(index) + offset
        bar = root.TBox(center - width / 2., min(baseline, value),
                        center + width / 2., max(baseline, value))
        bar.SetFillColorAlpha(color, 0.35)
        bar.SetLineColor(color)
        bar.SetLineWidth(2)
        bars.append(bar)
    return bars


def draw_comparison(root, comparisons, output):
    """Grouped bars for 1-3 normalized HJ/ggH ratios, with errors and lookup."""
    if not 1 <= len(comparisons) <= 3:
        raise ValueError("Provide between one and three comparisons.")
    root.gStyle.SetOptStat(0)
    root.gStyle.SetTextFont(42)
    root.gStyle.SetLegendFont(42)
    root.gStyle.SetCanvasColor(0)
    root.gStyle.SetPadColor(0)
    root.gStyle.SetEndErrorSize(0)
    canvas = root.TCanvas("nnlo_to_nlo_canvas", "", 900, 800)
    top = root.TPad("nnlo_to_nlo_top", "", 0., 0.30, 1., 1.)
    bottom = root.TPad("nnlo_to_nlo_bottom", "", 0., 0., 1., 0.30)
    for pad in (top, bottom):
        pad.SetLeftMargin(0.20)
        pad.SetRightMargin(0.05)
        pad.SetTicks(1, 1)
    top.SetTopMargin(0.08)
    top.SetBottomMargin(0.025)
    bottom.SetTopMargin(0.025)
    bottom.SetBottomMargin(0.55)
    top.Draw()
    bottom.Draw()
    keep = [top, bottom]
    bins = len(PT_BIN_LABELS)
    axes = root.TH1D("nnlo_to_nlo_axes", "", bins, 0.5, bins + 0.5)
    axes.SetDirectory(0)
    for b, label in enumerate(PT_BIN_LABELS, 1):
        axes.GetXaxis().SetBinLabel(b, label)
    axes.SetCanExtend(root.TH1.kNoAxis)
    raw = axes.Clone("raw_nnlo_lookup")
    for b, value in enumerate(lookup_weights(), 1):
        raw.SetBinContent(b, value)
    raw.SetLineColor(root.TColor.GetColor("#D62728"))
    raw.SetLineWidth(3)
    raw.SetFillStyle(0)
    curves = []
    bar_width = 0.78 / len(comparisons)
    for index, (comparison, color) in enumerate(zip(comparisons, ("#5B88CF", "#2CA02C", "#9467BD")), 1):
        color = root.TColor.GetColor(color)
        offset = -0.39 + (index - 0.5) * bar_width
        observed, closure = make_graphs(root, comparison["rows"], "comparison" + str(index), offset)
        for graph in (observed, closure):
            graph.SetLineColor(color)
            graph.SetMarkerColor(color)
            graph.SetMarkerStyle(20)
            graph.SetMarkerSize(0.)
            graph.SetLineWidth(2)
        legend_bar = root.TBox()
        legend_bar.SetFillColorAlpha(color, 0.35)
        legend_bar.SetLineColor(color)
        legend_bar.SetLineWidth(2)
        curves.append({"observed": observed, "closure": closure, "color": color,
                       "offset": offset, "rows": comparison["rows"], "legend": legend_bar})
        keep.extend((observed, closure, legend_bar))

    top.cd()
    rows = [row for comparison in comparisons for row in comparison["rows"]]
    upper = max([1.5] + list(lookup_weights())
                + [row["ratio"] + row["error"] for row in rows if row["ratio"] is not None])
    lower = min([0.] + [row["ratio"] - row["error"] for row in rows if row["ratio"] is not None])
    axes.SetMinimum(lower * 1.15)
    axes.SetMaximum(upper * 1.45)
    axes.GetXaxis().SetLabelSize(0.)
    axes.GetXaxis().SetTitleSize(0.)
    axes.GetYaxis().SetLabelSize(0.049)
    axes.GetYaxis().SetNdivisions(508)
    axes.Draw("AXIS")
    for curve in curves:
        bars = make_bars(root, curve["rows"], "ratio", curve["color"], bar_width,
                         curve["offset"], (axes.GetMinimum(), axes.GetMaximum()))
        for bar in bars:
            bar.Draw("SAME")
        keep.extend(bars)
    unity = root.TLine(0.5, 1., bins + 0.5, 1.)
    unity.SetLineColor(root.kBlack)
    unity.SetLineStyle(2)
    unity.Draw()
    raw.Draw("HIST SAME")
    for curve in curves:
        curve["observed"].Draw("Z SAME")
    keep.append(draw_y_title(root, top, "NNLO / NLO", "nnlo_to_nlo_y", x=0.11))
    legend = root.TLegend(0.24, 0.88 - 0.055 * (len(comparisons) + 1), 0.92, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.040)
    legend.SetMargin(0.09)
    for comparison, curve in zip(comparisons, curves):
        legend.AddEntry(curve["legend"], " " + comparison["name"], "f")
    legend.AddEntry(raw, " w_{NNLO} lookup (Run 3)", "l")
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
    top.RedrawAxis()

    bottom.cd()
    ratio_axes = axes.Clone("nnlo_to_nlo_closure_axes")
    ratio_axes.SetDirectory(0)
    ratio_axes.SetMinimum(0.)
    ratio_axes.SetMaximum(2.)
    ratio_axes.GetXaxis().SetTitle("p_{T}(H_{gen}) interval [GeV]")
    ratio_axes.GetXaxis().SetTitleFont(42)
    ratio_axes.GetXaxis().SetLabelFont(42)
    ratio_axes.GetXaxis().SetTitleSize(0.125)
    ratio_axes.GetXaxis().SetTitleOffset(2.0)
    ratio_axes.GetXaxis().SetLabelSize(0.085)
    ratio_axes.GetXaxis().SetLabelOffset(0.015)
    ratio_axes.GetXaxis().LabelsOption("v")
    ratio_axes.GetYaxis().SetLabelSize(0.108)
    ratio_axes.GetYaxis().SetNdivisions(505)
    ratio_axes.Draw("AXIS")
    for curve in curves:
        bars = make_bars(root, curve["rows"], "closure", curve["color"], bar_width,
                         curve["offset"], (0., 2.))
        for bar in bars:
            bar.Draw("SAME")
        keep.extend(bars)
    line = root.TLine(0.5, 1., bins + 0.5, 1.)
    line.SetLineColor(root.kBlack)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()
    for curve in curves:
        curve["closure"].Draw("Z SAME")
    keep.append(draw_y_title(root, bottom, "Ratio / w_{NNLO}", "nnlo_closure_y", x=0.11))
    bottom.RedrawAxis()
    keep.extend((axes, raw, unity, legend, cms, energy, ratio_axes, line))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(output))
    canvas.Close()
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError("ROOT did not create the comparison PNG: " + str(output))


def print_comparison(name, rows, summary):
    print("\n{}: sum(weight/w_nnlo) in range: NLO={nlo_sumw:.12g}, NNLO={nnlo_sumw:.12g}".format(
        name, **summary), flush=True)
    print("pT [GeV] | sumW(NLO) | sumW(NNLO) | NNLO/NLO +/- stat | lookup | ratio/lookup", flush=True)
    for row in rows:
        if row["ratio"] is None:
            measured, closure = "undefined (zero NLO bin weight)", "undefined"
        else:
            measured = "{:.6f} +/- {:.6f}".format(row["ratio"], row["error"])
            closure = "{:.6f}".format(row["closure"])
        print("{:g}-{:g} | {:.8g} | {:.8g} | {} | {:.12g} | {}".format(
            row["low"], row["high"], row["nlo_sumw"], row["nnlo_sumw"],
            measured, row["lookup"], closure), flush=True)


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
        output = output_path(args)
        if output.exists():
            raise ValueError("Output already exists: {}. Rename it or remove it explicitly before rerunning.".format(output))
        print("Lookup: pinned Run 3 table; 16 generator-Higgs-pT intervals, 0 <= pT < 500 GeV.", flush=True)
        comparisons = []
        for group in args.comparisons:
            sums, squared = {}, {}
            for role in ("nlo", "nnlo"):
                location = group[role + "_location"]
                label = "{}: {} ({})".format(group["name"], role.upper(), "ggH" if role == "nlo" else "HJ")
                print(label + ": discovering " + location, flush=True)
                files, _ = list_pico_files(location)
                print("  ALL {} ROOT files; signed weight/w_nnlo".format(len(files)), flush=True)
                print("  Processing events directly (no preliminary file-validation scan)...", flush=True)
                histograms, stats = book_histograms(ROOT, files, "sample{}_{}".format(group["index"], role))
                hist = histograms[0]
                sums[role] = [float(hist.GetBinContent(b)) for b in range(1, hist.GetNbinsX() + 1)]
                squared[role] = [float(hist.GetBinError(b)) ** 2 for b in range(1, hist.GetNbinsX() + 1)]
                print("{}: total={total:,}; plotted={plotted:,}; no Higgs={no_higgs:,}; "
                      "pT >= 500 GeV omitted={above_range:,}".format(
                          label, plotted=stats["counts"][0], **stats), flush=True)
            rows, summary = build_comparison(sums["nlo"], sums["nnlo"],
                                              squared["nlo"], squared["nnlo"], lookup_weights())
            print_comparison(group["name"], rows, summary)
            comparisons.append({"name": group["name"], "rows": rows, "summary": summary})
        draw_comparison(ROOT, comparisons, output)
        print("Saved " + str(output), flush=True)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "NNLO/NLO comparison failed: {}\n".format(error))


if __name__ == "__main__":
    main()
