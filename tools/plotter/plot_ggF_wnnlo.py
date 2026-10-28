#!/usr/bin/env python3
"""Compare the w_nnlo correction itself in central ggF and private HJ picos.

Usage: python3 tools/plotter/plot_ggF_wnnlo.py 2022
Each event contributes ONE, not weight, weight/w_nnlo, or w_nnlo.
Use all configured baseline pico files; apply no further event/object cuts.
"""

import argparse
from pathlib import Path
import subprocess

from plot_ggF_BDT_var import SAMPLE_DIRECTORIES, configured_year, list_pico_files, draw_plots

# Pinned from richstu/nano2pico, tag htozgamma_sequoia_v1:
# data/zgamma/GluGluHToZG_NNLO_reweight_run{2,3}.json.
# Keep all digits; category matching accepts the stored Float_t representation.
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


def category_values(year):
    """Sixteen lookup values plus unity, sorted by the correction itself."""
    return tuple(sorted(set((*lookup_weights(year), 1.))))


def category_expression(values):
    # Compare at Float_t precision to accept both JSON doubles and their
    # stored pico floats, without snapping unrelated corrections to a bin.
    constants = ", ".join(format(value, ".17g") for value in values)
    return (
        "const double expected[] = {" + constants + "}; "
        "const float value = static_cast<float>(hj_nnlo_value); "
        "for (unsigned int i = 0; i < sizeof(expected)/sizeof(expected[0]); ++i) { "
        "  if (value == static_cast<float>(expected[i])) return static_cast<int>(i + 1); "
        "} return 0;"
    )


def book_histogram(frame, sample_name, year):
    values = category_values(year)
    bins = len(values)
    categorized = frame.Define("hj_nnlo_bin", category_expression(values))
    unknown = categorized.Filter("hj_nnlo_bin == 0")
    unknown_count = unknown.Count()
    unknown_min = unknown.Min("hj_nnlo_value")
    unknown_max = unknown.Max("hj_nnlo_value")
    # No weight argument: one event contributes one to its correction category.
    action = categorized.Filter("hj_nnlo_bin > 0").Histo1D(
        (sample_name + "_wnnlo", "", bins, 0.5, bins + 0.5), "hj_nnlo_bin")
    missing = int(unknown_count.GetValue())
    if missing:
        raise ValueError(
            "{}: {:,} events have unknown w_nnlo values (range {:.17g} to {:.17g}); "
            "expected the lookup corrections or 1, allowing Float_t rounding.".format(
                sample_name, missing, unknown_min.GetValue(), unknown_max.GetValue()))
    hist = action.GetValue().Clone()
    hist.SetDirectory(0)
    for index in range(1, bins + 1):
        hist.GetXaxis().SetBinLabel(index, str(index))
    hist.SetCanExtend(hist.kNoAxis)
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
        values = category_values(args.year)
        bins = len(values)
        unity_bin = values.index(1.) + 1
        print("Numbered w_nnlo categories (ascending; unit event weights):", flush=True)
        for index, value in enumerate(values, 1):
            print("  {:2d}: {!r}{}".format(
                index, value, "  [unity reference]" if value == 1. else ""), flush=True)
        histograms = {sample: book_histogram(frame, sample.lower(), args.year)
                      for sample, frame in frames.items()}
        for sample, hist in histograms.items():
            if hist.GetEntries() != summaries[sample]["total"]:
                raise RuntimeError("Input event count changed while histogramming " + sample)
        output = Path(__file__).resolve().parents[2] / "plots" / ("ggF_wnnlo_" + args.year + ".png")
        draw_plots(ROOT, [histograms["Central"]], [histograms["Private"]], args.year, output,
                   event_counts={sample: (summary["total"],) for sample, summary in summaries.items()},
                   panels=(("w_nnlo", "NNLO correction category", "w_{NNLO} category (ascending)",
                            bins, 0.5, bins + 0.5),),
                   columns=1, panel_size=(900, 800), legend_text_size=0.040,
                   reference_lines={"w_nnlo": float(unity_bin)}, y_title_x=0.11)
        print("Saved " + str(output), flush=True)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "w_nnlo comparison failed: {}\n".format(error))


if __name__ == "__main__":
    main()
