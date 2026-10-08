#!/usr/bin/env python3
"""Compare six HJ/central ggF kinematic shapes in baseline-selected picos.

Usage: python3 tools/plotter/plot_ggF_kinematics.py 2022
Read all configured baseline pico files; produce one 3x2 comparison PNG.
"""

import argparse
from pathlib import Path
import subprocess

from plot_ggF_BDT_var import SAMPLE_DIRECTORIES, configured_year, list_pico_files, draw_plots

# column, panel title, x-axis label, bin count, lower edge, upper edge
PANELS = (
    ("pt_gamma", "Photon transverse momentum", "p_{T}(#gamma) [GeV]", 50, 0., 240.),
    ("pt_ll", "Dilepton transverse momentum", "p_{T}(#it{l}#it{l}) [GeV]", 50, 0., 240.),
    ("pt_llgamma", "Dilepton-photon transverse momentum", "p_{T}(#it{l}#it{l}#gamma) [GeV]", 50, 0., 240.),
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
        raise RuntimeError("ROOT could not compile the kinematic variable helpers.")
    root.TH1.SetDefaultSumw2(True)
    inputs = root.std.vector("string")()
    for filename in files:
        inputs.push_back(filename)
    frame = root.RDataFrame("tree", inputs)
    total = frame.Count()
    weighted = frame.Define("hj_kin_weight", "static_cast<double>(weight)")
    finite = weighted.Filter("std::isfinite(hj_kin_weight)", "Finite stored pico weight")
    usable = finite.Count()
    sumw = finite.Sum("hj_kin_weight")
    values = finite.Define("hj_kin_values", "hjpicokin::build(photon_pt, ll_pt, ll_m, "
                           "llphoton_pt, llphoton_m, llphoton_iph, llphoton_ill, npv)")
    actions, counts = [], []
    for column, _, _, bins, low, high in PANELS:
        scalar = "hj_kin_" + column
        defined = values.Define(scalar, "hj_kin_values." + column)
        # Validity checks only: do not reapply object or baseline selections.
        defined = defined.Filter("std::isfinite(" + scalar + ")", "Defined " + column)
        counts.append(defined.Count())
        actions.append(defined.Histo1D(
            (sample_name + "_kin_" + column, "", bins, low, high), scalar, "hj_kin_weight"))
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
        output = Path(__file__).resolve().parents[2] / "plots" / ("ggF_kinematics_" + args.year + ".png")
        draw_plots(ROOT, samples["Central"], samples["Private"], args.year, output,
                   event_counts=event_counts, panels=PANELS, columns=3, panel_size=(600, 700),
                   reference_lines={"m_llgamma": 125.}, legend_text_size=0.040)
        print("Saved " + str(output))
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "Kinematic plotting failed: {}\n".format(error))


if __name__ == "__main__":
    main()
