#!/usr/bin/env python3
"""Compare the generator Higgs pT used by nano2pico's NNLO lookup.

Usage: python3 tools/plotter/plot_ggF_Higgs_pT.py 2022
Select the LAST mc_id == 25 array element, not the highest-pT Higgs or an
explicit isLastCopy selection. Every event contributes ONE in both samples.
Read all configured baseline picos without further physics selections.
"""

import argparse
from pathlib import Path
import subprocess

from plot_ggF_BDT_var import SAMPLE_DIRECTORIES, configured_year, list_pico_files, draw_plots

PANELS = (
    ("higgs_pt", "Generator Higgs transverse momentum", "p_{T}(H_{gen}) [GeV]",
     50, 0., 500.),
)
REQUIRED_BRANCHES = ("mc_id", "mc_pt")

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


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("year", type=configured_year, metavar="YEAR")
    return parser, parser.parse_args(argv)


def validate_files(root, files):
    """Require the relevant pico branches in every file, including empty files."""
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
        raise RuntimeError("ROOT could not compile the generator Higgs pT helper.")
    root.TH1.SetDefaultSumw2(True)
    inputs = root.std.vector("string")()
    for filename in files:
        inputs.push_back(filename)
    frame = root.RDataFrame("tree", inputs)
    total = frame.Count()
    values = frame.Define("hj_higgs_values", "hjpicohiggspt::build(mc_id, mc_pt)")
    malformed = values.Filter("hj_higgs_values.malformed").Count()
    well_formed = values.Filter("!hj_higgs_values.malformed")
    no_higgs = well_formed.Filter("hj_higgs_values.higgs_count == 0").Count()
    one_higgs = well_formed.Filter("hj_higgs_values.higgs_count == 1").Count()
    multiple_higgs = well_formed.Filter("hj_higgs_values.higgs_count > 1").Count()
    selected = well_formed.Filter("hj_higgs_values.higgs_count > 0").Define(
        "hj_higgs_pt", "hj_higgs_values.pt")
    valid = selected.Filter("std::isfinite(hj_higgs_pt) && hj_higgs_pt >= 0.",
                            "Valid last Higgs pT")
    selected_count, valid_count = selected.Count(), valid.Count()
    column, _, _, bins, low, high = PANELS[0]
    # No weight column: every event with a valid last Higgs contributes one.
    action = valid.Histo1D((sample_name + "_" + column, "", bins, low, high),
                          "hj_higgs_pt")
    # Book histogram and diagnostics before triggering one event loop per sample.
    total_events = int(total.GetValue())
    if total_events == 0:
        raise ValueError(sample_name + ": no events")
    if int(malformed.GetValue()):
        raise ValueError("{}: {} events have different mc_id/mc_pt vector lengths.".format(
            sample_name, int(malformed.GetValue())))
    if int(valid_count.GetValue()) != int(selected_count.GetValue()):
        raise ValueError(sample_name + ": invalid last Higgs pT (negative or non-finite).")
    count = int(valid_count.GetValue())
    if count == 0:
        raise ValueError(sample_name + ": no saved Higgs; cannot plot generator Higgs pT.")
    hist = action.GetValue().Clone()
    hist.SetDirectory(0)
    return [hist], {
        "total": total_events, "no_higgs": int(no_higgs.GetValue()),
        "one_higgs": int(one_higgs.GetValue()),
        "multiple_higgs": int(multiple_higgs.GetValue()),
        "counts": (count,),
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
            print("{}: ALL {} pico ROOT files in {}; unit event weights".format(
                sample, len(files), directory), flush=True)
            entries = validate_files(ROOT, files)
            histograms, summary = book_histograms(ROOT, files, sample.lower())
            if summary["total"] != entries:
                raise RuntimeError("Input event count changed while reading " + sample)
            print(("{}: N={total:,}; no Higgs={no_higgs:,}, one Higgs={one_higgs:,}, "
                   "multiple Higgs={multiple_higgs:,}").format(sample, **summary), flush=True)
            print("  Plotted N={:,}; unit event weights; events with no saved Higgs omitted".format(
                summary["counts"][0]), flush=True)
            samples[sample], event_counts[sample] = histograms, summary["counts"]
        output = Path(__file__).resolve().parents[2] / "plots" / ("ggF_Higgs_pT_" + args.year + ".png")
        draw_plots(ROOT, samples["Central"], samples["Private"], args.year, output,
                   event_counts=event_counts, panels=PANELS, columns=1,
                   panel_size=(900, 800), legend_text_size=0.040, y_title_x=0.11)
        print("Saved " + str(output), flush=True)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "Generator Higgs pT comparison failed: {}\n".format(error))


if __name__ == "__main__":
    main()
