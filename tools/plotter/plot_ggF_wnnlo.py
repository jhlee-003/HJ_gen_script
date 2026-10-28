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
import struct

from plot_ggF_BDT_var import SAMPLE_DIRECTORIES, configured_year, list_pico_files, draw_plots

# Pinned from richstu/nano2pico, tag htozgamma_sequoia_v1:
# data/zgamma/GluGluHToZG_NNLO_reweight_run{2,3}.json.
# Keep all digits; protect both these doubles and their stored float versions.
LOOKUP_WEIGHTS = {
    "run2": (
        1.2754058804973247,
        1.0584550668422257,
        0.9757041089557047,
        0.9442594094230987,
        0.9399967767292806,
        0.934984755937246,
        0.943758946143129,
        0.9633701876518976,
        0.9837051483557692,
        1.0140467126363697,
        1.0505937343714287,
        1.06645850270649,
        1.0980205331522623,
        1.197874444820148,
        1.4473128627697647,
        1.6606953569876508,
    ),
    "run3": (
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
    ),
}


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


def lookup_weights(year):
    run = "run2" if str(year)[:4] in ("2016", "2017", "2018") else "run3"
    return LOOKUP_WEIGHTS[run]


def protected_groups(values):
    """Keep each lookup double and its Float_t representation in the same bin."""
    groups = []
    for value in sorted(set(float(value) for value in values)):
        stored = struct.unpack("f", struct.pack("f", value))[0]
        groups.append((min(value, stored), max(value, stored)))
    if len(groups) < 3 or any(left[1] >= right[0] for left, right in zip(groups, groups[1:])):
        raise ValueError("Need at least three distinct, float-resolvable lookup values.")
    return tuple(groups)


def separates_values(binning, groups, low, high):
    bins, start, stop = binning
    if not start < low <= high < stop:
        return False
    width = (stop - start) / bins
    occupied = set()
    for first, last in groups:
        i = math.floor(bins * (first - start) / (stop - start))
        j = math.floor(bins * (last - start) / (stop - start))
        if i != j or not 0 <= i < bins or i in occupied:
            return False
        # Leave room for floating-point evaluation of bin indices in ROOT.
        margin = min(first - (start + i * width), start + (i + 1) * width - last)
        if margin <= 64 * math.ulp(max(abs(first), abs(last), 1.)):
            return False
        occupied.add(i)
    return True


def minimum_uniform_grid(groups, low, high):
    """Minimize uniform-bin count over width AND origin, without extra packages.

    The order of point phases x modulo width changes only at widths
    (x_j-x_i)/k. Between consecutive such widths, and consecutive point
    phases, bin membership is constant. Search those cells, rather than a
    sampled grid of widths/offsets. Protect Float_t/Double_t envelopes from
    splitting and keep range extrema strictly inside the histogram.
    """
    if not math.isfinite(low) or not math.isfinite(high) or low >= high:
        raise ValueError("Need a finite nonempty w_nnlo range.")
    span = high - low
    padding = max(0.05 * span, 0.05 * max(abs(low), abs(high), 1.))
    start, stop = low - padding, high + padding
    if not math.isfinite(stop - start):
        raise ValueError("w_nnlo values are too large for a finite plotting range.")

    # A guaranteed feasible starting grid: bins narrower than the smallest
    # gap can never combine two groups. Search counts, not only that bound.
    bins = len(groups)
    while not separates_values((bins, start, stop), groups, low, high):
        bins += 1
    best = bins, start, stop
    if bins == len(groups):
        return best

    # Any improvement must cover the span with at most bins-1 bins. Three
    # adjacent groups need two boundaries, limiting their boundary spacing.
    min_width = span / (bins - 1)
    max_width = min(groups[i + 2][0] - groups[i][1]
                    for i in range(len(groups) - 2))
    if min_width >= max_width:
        return best
    points = sorted({low, high, *(value for group in groups for value in group)})
    critical = {min_width, max_width}
    for index, first in enumerate(points):
        for last in points[index + 1:]:
            difference = last - first
            first_k = max(1, math.ceil(difference / max_width))
            last_k = math.floor(difference / min_width)
            for k in range(first_k, last_k + 1):
                width = difference / k
                if min_width < width < max_width:
                    critical.add(width)

    widths = sorted(critical, reverse=True)
    for upper, lower in zip(widths, widths[1:]):
        if upper <= span / (best[0] - 1):
            break  # Smaller widths cannot improve the current bin count.
        width = lower + (upper - lower) / 2.
        phases = sorted({0., width, *(point % width for point in points)})
        for left, right in zip(phases, phases[1:]):
            origin = left + (right - left) / 2.
            first_bin = math.floor((low - origin) / width)
            last_bin = math.floor((high - origin) / width)
            bins = last_bin - first_bin + 1
            if bins >= best[0]:
                continue
            start = origin + first_bin * width
            candidate = bins, start, start + bins * width
            if separates_values(candidate, groups, low, high):
                best = candidate
    return best


def optimal_binning(summaries, year):
    """Cover ALL observed events and separate the year's lookup values + unity."""
    groups = protected_groups((*lookup_weights(year), 1.))
    low = min(groups[0][0], *(summary["min"] for summary in summaries))
    high = max(groups[-1][1], *(summary["max"] for summary in summaries))
    return minimum_uniform_grid(groups, low, high)


def book_histogram(frame, sample_name, binning):
    bins, low, high = binning
    # No third (weight) argument: every event contributes exactly one.
    action = frame.Histo1D((sample_name + "_wnnlo", "", bins, low, high), "hj_nnlo_value")
    hist = action.GetValue().Clone()
    hist.SetDirectory(0)
    if hist.GetBinContent(0) or hist.GetBinContent(hist.GetNbinsX() + 1):
        raise RuntimeError("Optimized w_nnlo range left events outside the histogram.")
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
        binning = optimal_binning(list(summaries.values()), args.year)
        bins, low, high = binning
        print("Optimal uniform w_nnlo grid: {} bins; width={:.12g}; range [{:.12g}, {:.12g})".format(
            bins, (high - low) / bins, low, high), flush=True)
        histograms = {sample: book_histogram(frame, sample.lower(), binning)
                      for sample, frame in frames.items()}
        for sample, hist in histograms.items():
            if hist.GetEntries() != summaries[sample]["total"]:
                raise RuntimeError("Input event count changed while histogramming " + sample)
        output = Path(__file__).resolve().parents[2] / "plots" / ("ggF_wnnlo_" + args.year + ".png")
        draw_plots(ROOT, [histograms["Central"]], [histograms["Private"]], args.year, output,
                   event_counts={sample: (summary["total"],) for sample, summary in summaries.items()},
                   panels=(("w_nnlo", "NNLO correction factor", "w_{NNLO}", bins, low, high),),
                   columns=1, panel_size=(900, 800), legend_text_size=0.040,
                   reference_lines={"w_nnlo": 1.}, y_title_x=0.11)
        print("Saved " + str(output), flush=True)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "w_nnlo comparison failed: {}\n".format(error))


if __name__ == "__main__":
    main()
