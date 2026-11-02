#!/usr/bin/env python3
"""Compare normalized HJ/ggH_qme Higgs-pT fractions with the NNLO lookup.

Usage: python3 tools/plotter/plot_ggF_NNLO_over_NLO.py 2022
Use one count per event and the last mc_id == 25 particle, as in the Higgs pT
plotter. HJ is the candidate NNLO sample; ggH_qme is the NLO reference.
"""

import argparse
import math
from pathlib import Path
import subprocess

from plot_ggF_BDT_var import SAMPLE_DIRECTORIES, configured_year, list_pico_files, draw_y_title
from plot_ggF_Higgs_pT import PT_BIN_EDGES, PT_BIN_LABELS, validate_files, book_histograms
from plot_ggF_wnnlo import lookup_weights


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("year", type=configured_year, metavar="YEAR")
    return parser, parser.parse_args(argv)


def build_comparison(central_counts, private_counts, weights):
    """Return normalized count ratios, binomial delta errors, and table curves."""
    bins = len(PT_BIN_LABELS)
    if not len(central_counts) == len(private_counts) == len(weights) == bins:
        raise ValueError("Counts and lookup must match all Higgs pT intervals.")
    for sample, counts in (("Central", central_counts), ("Private", private_counts)):
        if any(not math.isfinite(n) or n < 0 or n != int(n) for n in counts):
            raise ValueError(sample + ": expected nonnegative integer event counts.")
        if sum(counts) <= 0:
            raise ValueError(sample + ": no events in the plotted pT range.")
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError("Lookup corrections must be finite and positive.")
    nc, nh = sum(central_counts), sum(private_counts)
    mean_weight = math.fsum(n * w for n, w in zip(central_counts, weights)) / nc
    rows = []
    for i, (n, m, w) in enumerate(zip(central_counts, private_counts, weights)):
        fc, fh = n / nc, m / nh
        expected = w / mean_weight
        ratio = error = closure = closure_error = None
        if n:
            ratio = fh / fc
            # Fractions are normalized to their own sample totals: use the
            # binomial variances, then propagate the independent samples.
            variance = (fh * (1. - fh) / nh / (fc * fc)
                        + fh * fh * fc * (1. - fc) / nc / (fc ** 4))
            error = math.sqrt(max(0., variance))
            closure, closure_error = ratio / expected, error / expected
        rows.append({
            "low": PT_BIN_EDGES[i], "high": PT_BIN_EDGES[i + 1],
            "central": int(n), "private": int(m),
            "central_fraction": fc, "private_fraction": fh,
            "ratio": ratio, "error": error, "lookup": w,
            "shape_lookup": expected, "closure": closure, "closure_error": closure_error,
        })
    return rows, {
        "central_total": int(nc), "private_total": int(nh),
        "mean_lookup": mean_weight,
        "undefined_bins": sum(row["ratio"] is None for row in rows),
    }


def make_graphs(root, rows):
    observed = root.TGraphErrors()
    observed.SetName("hj_over_ggh_shape")
    closure = root.TGraphErrors()
    closure.SetName("hj_over_ggh_over_shape_lookup")
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
    """Draw ratio and lookup values directly; never normalize these curves."""
    root.gStyle.SetOptStat(0)
    root.gStyle.SetTextFont(42)
    root.gStyle.SetLegendFont(42)
    root.gStyle.SetEndErrorSize(0)
    canvas = root.TCanvas("nnlo_over_nlo_canvas", "", 1000, 850)
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
    expected = axes.Clone("shape_normalized_nnlo_lookup")
    for b, row in enumerate(rows, 1):
        raw.SetBinContent(b, row["lookup"])
        expected.SetBinContent(b, row["shape_lookup"])
    raw.SetLineColor(root.TColor.GetColor("#D62728"))
    raw.SetLineWidth(3)
    expected.SetLineColor(root.kBlack)
    expected.SetLineStyle(2)
    expected.SetLineWidth(2)
    for hist in (raw, expected):
        hist.SetFillStyle(0)
    observed, closure = make_graphs(root, rows)
    observed.SetLineColor(root.TColor.GetColor("#5B88CF"))
    observed.SetMarkerColor(root.TColor.GetColor("#5B88CF"))
    observed.SetMarkerStyle(20)
    observed.SetMarkerSize(0.85)
    observed.SetLineWidth(2)
    closure.SetLineColor(root.kBlack)
    closure.SetMarkerColor(root.kBlack)
    closure.SetMarkerStyle(20)
    closure.SetMarkerSize(0.75)

    top.cd()
    upper = max([1.5] + [row["lookup"] for row in rows]
                + [row["shape_lookup"] for row in rows]
                + [row["ratio"] + row["error"] for row in rows if row["ratio"] is not None])
    axes.SetMinimum(0.)
    axes.SetMaximum(upper * 1.32)
    axes.GetXaxis().SetLabelSize(0.)
    axes.GetXaxis().SetTitleSize(0.)
    axes.GetYaxis().SetLabelSize(0.045)
    axes.GetYaxis().SetNdivisions(508)
    axes.Draw("AXIS")
    unity = root.TLine(0.5, 1., bins + 0.5, 1.)
    unity.SetLineColor(root.kGray + 2)
    unity.SetLineStyle(3)
    unity.Draw()
    raw.Draw("HIST SAME")
    expected.Draw("HIST SAME")
    observed.Draw("P SAME")
    keep.append(draw_y_title(root, top, "NNLO / NLO", "nnlo_over_nlo_y", x=0.11))
    legend = root.TLegend(0.235, 0.695, 0.65, 0.885)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.034)
    legend.AddEntry(observed, "HJ / ggH_qme (unit-area)", "pe")
    legend.AddEntry(raw, "Lookup w_{NNLO}", "l")
    legend.AddEntry(expected, "Shape-normalized lookup", "l")
    legend.Draw()
    counts = root.TLegend(0.69, 0.745, 0.94, 0.88)
    counts.SetBorderSize(0)
    counts.SetFillStyle(0)
    counts.SetTextSize(0.034)
    counts.SetMargin(0.)
    counts.AddEntry(root.nullptr, "N(ggH)={:,}".format(summary["central_total"]), "")
    counts.AddEntry(root.nullptr, "N(HJ)={:,}".format(summary["private_total"]), "")
    counts.Draw()
    cms = root.TLatex()
    cms.SetNDC(True)
    cms.SetTextSize(0.052)
    cms.DrawLatex(0.195, 0.935, "#bf{CMS}#kern[-0.4]{ }#scale[0.9]{#it{Private Work}}")
    cms.SetTextAlign(31)
    cms.SetTextSize(0.045)
    cms.DrawLatex(0.95, 0.935, "13.6 TeV")
    top.RedrawAxis()

    bottom.cd()
    ratio_axes = axes.Clone("nnlo_over_nlo_closure_axes")
    ratio_axes.SetDirectory(0)
    ratio_axes.SetMinimum(0.)
    ratio_axes.SetMaximum(max(2., max(
        (row["closure"] + row["closure_error"] for row in rows if row["closure"] is not None),
        default=1.) * 1.08))
    ratio_axes.GetXaxis().SetTitle("p_{T}(H_{gen}) interval [GeV]")
    ratio_axes.GetXaxis().SetTitleSize(0.125)
    ratio_axes.GetXaxis().SetTitleOffset(2.4)
    ratio_axes.GetXaxis().SetLabelSize(0.085)
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
    keep.append(draw_y_title(root, bottom, "Ratio / shape lookup", "nnlo_closure_y", x=0.11))
    bottom.RedrawAxis()
    keep.extend((axes, raw, expected, observed, closure, unity, legend, counts, cms, ratio_axes, line))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(output))
    canvas.Close()
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError("ROOT did not create the comparison PNG: " + str(output))


def print_comparison(rows, summary):
    print("Central-average lookup = {:.12g}".format(summary["mean_lookup"]), flush=True)
    print("pT [GeV] | HJ/ggH +/- stat | lookup | shape lookup | ratio/shape lookup", flush=True)
    for row in rows:
        if row["ratio"] is None:
            measured, closure = "undefined (zero ggH)", "undefined"
        else:
            measured = "{:.6f} +/- {:.6f}".format(row["ratio"], row["error"])
            closure = "{:.6f}".format(row["closure"])
        print("{:g}-{:g} | {} | {:.12g} | {:.12g} | {}".format(
            row["low"], row["high"], measured, row["lookup"], row["shape_lookup"], closure),
            flush=True)


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
        counts = {}
        for sample, key in (("Central", "central"), ("Private", "private")):
            directory = SAMPLE_DIRECTORIES[args.year][key]
            files, _ = list_pico_files(directory)
            print("{}: ALL {} pico ROOT files in {}; unit event weights".format(
                sample, len(files), directory), flush=True)
            entries = validate_files(ROOT, files)
            histograms, summary = book_histograms(ROOT, files, sample.lower() + "_comparison")
            if summary["total"] != entries:
                raise RuntimeError("Input event count changed while reading " + sample)
            hist = histograms[0]
            counts[sample] = [int(hist.GetBinContent(b)) for b in range(1, hist.GetNbinsX() + 1)]
            if sum(counts[sample]) != summary["counts"][0]:
                raise RuntimeError("Histogram count disagrees with plotted events for " + sample)
            print("{}: total={total:,}; plotted={plotted:,}; no Higgs={no_higgs:,}; "
                  "pT >= 500 GeV omitted={above_range:,}".format(
                      sample, plotted=summary["counts"][0], **summary), flush=True)
        rows, comparison = build_comparison(counts["Central"], counts["Private"],
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
