#!/usr/bin/env python3
"""Compare the w_nnlo correction itself in central ggF and private HJ picos.

Usage: python3 tools/plotter/plot_ggF_wnnlo.py 2022
Each event contributes ONE, not weight, weight/w_nnlo, or w_nnlo.
Use all configured baseline pico files; apply no further event/object cuts.
"""

import argparse
import math
from pathlib import Path
import subprocess

from plot_ggF_BDT_var import SAMPLE_DIRECTORIES, configured_year, list_pico_files, draw_plots

BINS = 60


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("year", type=configured_year, metavar="YEAR")
    return parser, parser.parse_args(argv)


def validate_files(root, files):
    """Only tree/w_nnlo is needed, in both samples. Never substitute unity."""
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
            branch = tree.GetBranch("w_nnlo")
            leaf = tree.GetLeaf("w_nnlo")
            if not branch or not leaf or leaf.GetLeafCount() or leaf.GetLenStatic() != 1:
                raise ValueError("Missing or non-scalar w_nnlo branch: " + filename)
            if leaf.GetTypeName() not in ("Float_t", "Double_t"):
                raise ValueError("w_nnlo must be a floating-point scalar: " + filename)
            entries += int(tree.GetEntries())
        finally:
            source.Close()
        if index % 100 == 0 or index == len(files):
            print("Validated {}/{} pico files".format(index, len(files)), flush=True)
    if entries == 0:
        raise ValueError("The pico sample contains no events.")
    return entries


def inspect_sample(root, files, sample_name):
    inputs = root.std.vector("string")()
    for filename in files:
        inputs.push_back(filename)
    frame = root.RDataFrame("tree", inputs).Define("hj_nnlo_value", "static_cast<double>(w_nnlo)")
    total = frame.Count()
    finite = frame.Filter("std::isfinite(hj_nnlo_value)")
    count, low, high, mean = (finite.Count(), finite.Min("hj_nnlo_value"),
                              finite.Max("hj_nnlo_value"), finite.Mean("hj_nnlo_value"))
    entries = int(total.GetValue())
    if entries == 0:
        raise ValueError(sample_name + ": no events")
    if int(count.GetValue()) != entries:
        raise ValueError(sample_name + ": non-finite w_nnlo; refusing to silently drop events.")
    return frame, {"total": entries, "min": float(low.GetValue()),
                   "max": float(high.GetValue()), "mean": float(mean.GetValue())}


def common_range(summaries):
    """Cover every correction in both samples, plus the unity reference."""
    low = min(1., *(summary["min"] for summary in summaries))
    high = max(1., *(summary["max"] for summary in summaries))
    padding = max(0.05 * (high - low), 0.05 * max(abs(low), abs(high), 1.))
    result = low - padding, high + padding
    if not all(math.isfinite(value) for value in result):
        raise ValueError("w_nnlo values are too large for a finite plotting range.")
    return result


def book_histogram(frame, sample_name, low, high):
    # No third (weight) argument: every event contributes exactly one.
    action = frame.Histo1D((sample_name + "_wnnlo", "", BINS, low, high), "hj_nnlo_value")
    hist = action.GetValue().Clone()
    hist.SetDirectory(0)
    return hist


def main():
    parser, args = parse_arguments()
    try:
        import ROOT
    except ImportError:
        parser.exit(1, "PyROOT is unavailable. Activate CMSSW, then run with its python3.\n")
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.EnableImplicitMT(2)
    ROOT.TH1.SetDefaultSumw2(True)
    try:
        frames, summaries = {}, {}
        for sample, key in (("Central", "central"), ("Private", "private")):
            directory = SAMPLE_DIRECTORIES[args.year][key]
            print(sample + ": discovering ALL baseline pico files in " + directory, flush=True)
            files, _ = list_pico_files(directory)
            print("{}: {} ROOT files; unit event weights".format(sample, len(files)), flush=True)
            expected = validate_files(ROOT, files)
            frame, summary = inspect_sample(ROOT, files, sample)
            if summary["total"] != expected:
                raise RuntimeError("Input event count changed while reading " + sample)
            frames[sample], summaries[sample] = frame, summary
            print("{}: N={total:,}; w_nnlo min={min:.8g}, mean={mean:.8g}, max={max:.8g}".format(
                sample, **summary), flush=True)
        low, high = common_range(list(summaries.values()))
        print("Shared w_nnlo range: {:.8g} to {:.8g}; {} bins".format(low, high, BINS), flush=True)
        histograms = {sample: book_histogram(frame, sample.lower(), low, high)
                      for sample, frame in frames.items()}
        for sample, hist in histograms.items():
            if hist.GetEntries() != summaries[sample]["total"]:
                raise RuntimeError("Input event count changed while histogramming " + sample)
        output = Path(__file__).resolve().parents[2] / "plots" / ("ggF_wnnlo_" + args.year + ".png")
        draw_plots(ROOT, [histograms["Central"]], [histograms["Private"]], args.year, output,
                   event_counts={sample: (summary["total"],) for sample, summary in summaries.items()},
                   panels=(("w_nnlo", "NNLO correction factor", "w_{NNLO}", BINS, low, high),),
                   columns=1, panel_size=(900, 800), legend_text_size=0.040,
                   reference_lines={"w_nnlo": 1.}, y_title_x=0.11)
        print("Saved " + str(output), flush=True)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "w_nnlo comparison failed: {}\n".format(error))


if __name__ == "__main__":
    main()
