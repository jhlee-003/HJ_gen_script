#!/usr/bin/env python3
"""Compare weighted HJ/ggH_qme Higgs-pT shapes with the raw NNLO lookup.

Usage: python3 tools/plotter/plot_ggF_NNLO_over_NLO.py 2022
Use signed weight/w_nnlo for BOTH samples, then normalize each shape.
Select the last mc_id == 25 particle, as in the Higgs pT plotter.
"""

import argparse
from array import array
import math
from pathlib import Path
import subprocess

from plot_ggF_BDT_var import SAMPLE_DIRECTORIES, configured_year, list_pico_files, draw_y_title
from plot_ggF_Higgs_pT import PT_BIN_EDGES, PT_BIN_LABELS, CPP_HELPERS
from plot_ggF_wnnlo import lookup_weights

REQUIRED_BRANCHES = ("mc_id", "mc_pt", "weight", "w_nnlo")


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("year", type=configured_year, metavar="YEAR")
    return parser, parser.parse_args(argv)


def validate_files(root, files):
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
            for name in ("weight", "w_nnlo"):
                leaf = tree.GetLeaf(name)
                if (not leaf or leaf.GetLeafCount() or leaf.GetLenStatic() != 1
                        or leaf.GetTypeName() not in ("Float_t", "Double_t")):
                    raise ValueError("{} must be a floating-point scalar: {}".format(name, filename))
            entries += int(tree.GetEntries())
        finally:
            source.Close()
        if index % 100 == 0 or index == len(files):
            print("Validated {}/{} pico files".format(index, len(files)), flush=True)
    if entries == 0:
        raise ValueError("The pico sample contains no events.")
    return entries


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
    finite = weighted.Filter("std::isfinite(nnlo_plot_weight)")
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


def build_comparison(central_sums, private_sums, central_sumw2, private_sumw2, weights):
    bins = len(PT_BIN_LABELS)
    if not all(len(values) == bins for values in (
            central_sums, private_sums, central_sumw2, private_sumw2, weights)):
        raise ValueError("Bin sums, squared weights, and lookup must match all pT intervals.")
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError("Lookup corrections must be finite and positive.")
    fc, vc, sc = normalized_fractions(central_sums, central_sumw2, "Central")
    fh, vh, sh = normalized_fractions(private_sums, private_sumw2, "Private")
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
            "central_sumw": central_sums[i], "private_sumw": private_sums[i],
            "central_sumw2": central_sumw2[i], "private_sumw2": private_sumw2[i],
            "central_fraction": fc[i], "private_fraction": fh[i],
            "ratio": ratio, "error": error, "lookup": w,
            "closure": closure, "closure_error": closure_error,
        })
    return rows, {
        "central_sumw": sc, "private_sumw": sh,
        "undefined_bins": sum(row["ratio"] is None for row in rows),
    }


def make_graphs(root, rows):
    observed = root.TGraphErrors()
    observed.SetName("hj_over_ggh_weighted_shape")
    closure = root.TGraphErrors()
    closure.SetName("hj_over_ggh_over_raw_lookup")
    for index, row in enumerate(rows, 1):
        if row["ratio"] is None:
            continue
        point = observed.GetN()
        observed.SetPoint(point, float(index), row["ratio"])
        observed.SetPointError(point, 0., row["error"])
        closure.SetPoint(point, float(index), row["closure"])
        closure.SetPointError(point, 0., row["closure_error"])
    return observed, closure


def draw_comparison(root, rows, summary, output):
    """Keep the established CMS overlay/ratio format; show only the raw lookup."""
    root.gStyle.SetOptStat(0)
    root.gStyle.SetTextFont(42)
    root.gStyle.SetLegendFont(42)
    root.gStyle.SetCanvasColor(0)
    root.gStyle.SetPadColor(0)
    root.gStyle.SetEndErrorSize(0)
    canvas = root.TCanvas("nnlo_over_nlo_canvas", "", 900, 800)
    top = root.TPad("nnlo_over_nlo_top", "", 0., 0.30, 1., 1.)
    bottom = root.TPad("nnlo_over_nlo_bottom", "", 0., 0., 1., 0.30)
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
    bins = len(rows)
    axes = root.TH1D("nnlo_over_nlo_axes", "", bins, 0.5, bins + 0.5)
    axes.SetDirectory(0)
    for b, label in enumerate(PT_BIN_LABELS, 1):
        axes.GetXaxis().SetBinLabel(b, label)
    axes.SetCanExtend(root.TH1.kNoAxis)
    raw = axes.Clone("raw_nnlo_lookup")
    for b, row in enumerate(rows, 1):
        raw.SetBinContent(b, row["lookup"])
    raw.SetLineColor(root.TColor.GetColor("#D62728"))
    raw.SetLineWidth(3)
    raw.SetFillStyle(0)
    observed, closure = make_graphs(root, rows)
    observed.SetLineColor(root.TColor.GetColor("#5B88CF"))
    observed.SetMarkerColor(root.TColor.GetColor("#5B88CF"))
    observed.SetMarkerStyle(20)
    observed.SetMarkerSize(0.85)
    observed.SetLineWidth(2)
    closure.SetLineColor(root.kBlack)
    closure.SetMarkerColor(root.kBlack)
    closure.SetMarkerStyle(20)
    closure.SetMarkerSize(0.45)

    top.cd()
    upper = max([1.5] + [row["lookup"] for row in rows]
                + [row["ratio"] + row["error"] for row in rows if row["ratio"] is not None])
    lower = min([0.] + [row["ratio"] - row["error"]
                         for row in rows if row["ratio"] is not None])
    axes.SetMinimum(lower * 1.15)
    axes.SetMaximum(upper * 1.35)
    axes.GetXaxis().SetLabelSize(0.)
    axes.GetXaxis().SetTitleSize(0.)
    axes.GetYaxis().SetLabelSize(0.049)
    axes.GetYaxis().SetNdivisions(508)
    axes.Draw("AXIS")
    unity = root.TLine(0.5, 1., bins + 0.5, 1.)
    unity.SetLineColor(root.kBlack)
    unity.SetLineStyle(2)
    unity.Draw()
    raw.Draw("HIST SAME")
    observed.Draw("P SAME")
    keep.append(draw_y_title(root, top, "NNLO / NLO", "nnlo_over_nlo_y", x=0.11))
    legend = root.TLegend(0.24, 0.72, 0.65, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.040)
    legend.SetMargin(0.11)
    legend.AddEntry(observed, " NNLO/NLO ratio", "pe")
    legend.AddEntry(raw, " w_{NNLO}", "l")
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
    ratio_axes = axes.Clone("nnlo_over_nlo_closure_axes")
    ratio_axes.SetDirectory(0)
    ratio_axes.SetMinimum(0.)
    ratio_axes.SetMaximum(2.)
    ratio_axes.GetXaxis().SetTitle("p_{T}(H_{gen}) interval [GeV]")
    ratio_axes.GetXaxis().SetTitleFont(42)
    ratio_axes.GetXaxis().SetLabelFont(42)
    ratio_axes.GetXaxis().SetTitleSize(0.125)
    ratio_axes.GetXaxis().SetTitleOffset(2.4)
    ratio_axes.GetXaxis().SetLabelSize(0.085)
    ratio_axes.GetXaxis().SetLabelOffset(0.015)
    ratio_axes.GetXaxis().LabelsOption("v")
    ratio_axes.GetYaxis().SetLabelSize(0.108)
    ratio_axes.GetYaxis().SetNdivisions(505)
    ratio_axes.Draw("AXIS")
    line = root.TLine(0.5, 1., bins + 0.5, 1.)
    line.SetLineColor(root.kBlack)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()
    closure.Draw("P SAME")
    keep.append(draw_y_title(root, bottom, "Ratio / w_{NNLO}", "nnlo_closure_y", x=0.11))
    bottom.RedrawAxis()
    keep.extend((axes, raw, observed, closure, unity, legend, cms, energy, ratio_axes, line))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(output))
    canvas.Close()
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError("ROOT did not create the comparison PNG: " + str(output))


def print_comparison(rows, summary):
    print("Sum(weight/w_nnlo): Central={central_sumw:.12g}, Private={private_sumw:.12g}".format(
        **summary), flush=True)
    print("pT [GeV] | sumW(ggH) | sumW(HJ) | HJ/ggH +/- stat | lookup | ratio/lookup", flush=True)
    for row in rows:
        if row["ratio"] is None:
            measured, closure = "undefined (zero ggH bin weight)", "undefined"
        else:
            measured = "{:.6f} +/- {:.6f}".format(row["ratio"], row["error"])
            closure = "{:.6f}".format(row["closure"])
        print("{:g}-{:g} | {:.8g} | {:.8g} | {} | {:.12g} | {}".format(
            row["low"], row["high"], row["central_sumw"], row["private_sumw"],
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
        sums, squared = {}, {}
        for sample, key in (("Central", "central"), ("Private", "private")):
            directory = SAMPLE_DIRECTORIES[args.year][key]
            files, _ = list_pico_files(directory)
            print("{}: ALL {} pico ROOT files in {}; signed weight/w_nnlo".format(
                sample, len(files), directory), flush=True)
            entries = validate_files(ROOT, files)
            histograms, summary = book_histograms(ROOT, files, sample.lower() + "_comparison")
            if summary["total"] != entries:
                raise RuntimeError("Input event count changed while reading " + sample)
            hist = histograms[0]
            sums[sample] = [float(hist.GetBinContent(b)) for b in range(1, hist.GetNbinsX() + 1)]
            squared[sample] = [float(hist.GetBinError(b)) ** 2
                               for b in range(1, hist.GetNbinsX() + 1)]
            print("{}: total={total:,}; plotted={plotted:,}; no Higgs={no_higgs:,}; "
                  "pT >= 500 GeV omitted={above_range:,}".format(
                      sample, plotted=summary["counts"][0], **summary), flush=True)
        rows, comparison = build_comparison(
            sums["Central"], sums["Private"], squared["Central"], squared["Private"],
            lookup_weights(args.year))
        print_comparison(rows, comparison)
        output = Path(__file__).resolve().parents[2] / "plots" / (
            "ggF_NNLO_over_NLO_" + args.year + ".png")
        draw_comparison(ROOT, rows, comparison, output)
        print("Saved " + str(output), flush=True)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "NNLO/NLO comparison failed: {}\n".format(error))


if __name__ == "__main__":
    main()
