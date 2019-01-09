#!/usr/bin/env python

"""Print main QG plots, comparing systematics samples

"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
import argparse
from array import array
from copy import deepcopy

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgp
import common_utils as cu

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)

# import hunter
# hunter.trace(module='comparator')

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


# Control plot output format
OUTPUT_FMT = "pdf"

TOTAL_LUMI = 35918


def do_plots(nominal_dir, herwig_dir, plot_dir,
             syst_up_dir=None, syst_up_label="",
             syst_down_dir=None, syst_down_label=""):
    # QG variable plots
    pt_bins = qgc.PT_BINS
    # pt_bins = qgc.THEORY_PT_BINS
    var_list = qgc.COMMON_VARS
    var_prepend = ""

    plot_dirname = "Dijet_QG_tighter"

    for ang in var_list:
        if ang.var not in qgc.ANGLE_REBIN_DICT:
            continue

        v = "%s%s_vs_pt" % (var_prepend, ang.var)

        for pt_ind, (start_val, end_val) in enumerate(pt_bins[:]):
            dijet_entries = []
            # Get all plots
            lw = 2

            ####################
            # DIJET REGION
            ####################

            # NOMINAL MC
            h2d_qcd_mc = grab_obj(os.path.join(nominal_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_dirname, v))
            qcd_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                 marker_color=qgc.QCD_COLOUR, marker_style=qgc.QCD_MARKER, marker_size=0,
                                 label=qgc.QCD_Dijet_LABEL + " [nominal]")
            nominal_hist = qgp.get_projection_plot(h2d_qcd_mc, start_val, end_val)
            dijet_entries.append((nominal_hist, qcd_kwargs_mc))

            # HERWIG MC
            if herwig_dir:
                h2d_herwig_qcd_mc = grab_obj(os.path.join(herwig_dir, qgc.QCD_HERWIG_FILENAME), "%s/%s" % (dj_dirname, v))
                colh = qgc.ZB_COLOUR
                qcd_herwig_kwargs_mc = dict(line_color=colh, line_width=lw, fill_color=colh,
                                     marker_color=colh, marker_style=qgc.QCD_MARKER, marker_size=0,
                                     label=qgc.QCD_Dijet_LABEL + " [Herwig++]", subplot=nominal_hist)
                dijet_entries.append((qgp.get_projection_plot(h2d_herwig_qcd_mc, start_val, end_val), qcd_herwig_kwargs_mc))

            # SHIFT UP
            if syst_up_dir:
                col = qgc.QCD_COLOURS[2]
                h2d_qcd_mc2 = grab_obj(os.path.join(syst_up_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_dirname, v))
                qcd_kwargs_mc2 = dict(line_color=col, line_width=lw, fill_color=col,
                                     marker_color=col, marker_style=qgc.QCD_MARKER, marker_size=0,
                                     label=qgc.QCD_Dijet_LABEL + " [" + syst_up_label + "]", subplot=nominal_hist)
                dijet_entries.append((qgp.get_projection_plot(h2d_qcd_mc2, start_val, end_val), qcd_kwargs_mc2))

            # SHIFT DOWN
            if syst_down_dir:
                col2 = qgc.QCD_COLOURS[3]
                h2d_qcd_mc3 = grab_obj(os.path.join(syst_down_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_dirname, v))
                qcd_kwargs_mc3 = dict(line_color=col2, line_width=lw, fill_color=col2,
                                     marker_color=col2, marker_style=qgc.QCD_MARKER, marker_size=0,
                                     label=qgc.QCD_Dijet_LABEL + " [" + syst_down_label + "]", subplot=nominal_hist)
                # h2d_qcd_mc3.Scale(TOTAL_LUMI)
                dijet_entries.append((qgp.get_projection_plot(h2d_qcd_mc3, start_val, end_val), qcd_kwargs_mc3))

            rebin = 5
            v_lower = v.lower()
            if "multiplicity" in v_lower:
                rebin = 2
            elif "flavour" in v_lower or "thrust" in v_lower:
                rebin = 2
            elif "ptd" in v_lower:
                rebin = 5

            xlim = None
            if "width" in v_lower or "ptd" in v_lower:
                xlim = [0, 1]
            elif"thrust" in v_lower:
                xlim = [0, 0.5]
            elif "multiplicity" in v_lower and "ak4" in nominal_dir.lower():
                xlim = [0, 100]
                xlim = [0, 80]
            if xlim:
                xlim[1] = max(xlim[1], nominal_hist.GetBinLowEdge(nominal_hist.GetNbinsX()+1))

            ylim = None
            if "flavour" in v_lower:
                ylim = (0, 1)
            elif "lha" in v_lower:
                ylim = (0, 5)
                ylim = None

            # plot_dir = os.path.join(root_dir, "plots_dy_vs_qcd_mc_vs_data")
            radius, pus = cu.get_jet_config_from_dirname(nominal_dir)
            jet_str = "AK%s PF %s" % (radius.upper(), pus.upper())
            subplot_title = "Variation / nominal"
            subplot_limits = (0.9, 1.1)

            dijet_entries_rebin = []
            this_rebins = qgc.ANGLE_REBIN_DICT[ang.var]
            rebin_hist_norminal = None
            for hist, kwargs in dijet_entries:
                rebin_hist = hist.Rebin(len(this_rebins)-1, hist.GetName()+"Rebin", array('d', this_rebins))
                if not rebin_hist_norminal:
                    rebin_hist_norminal = rebin_hist
                new_kwargs = deepcopy(kwargs)
                if 'subplot' in new_kwargs:
                    new_kwargs['subplot'] = rebin_hist_norminal
                dijet_entries_rebin.append((rebin_hist, new_kwargs))

            qgp.do_comparison_plot(dijet_entries,
                                   "%s/ptBinned/%s_pt%dto%d_dijet.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin,
                                   title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim,
                                   subplot_type='ratio',
                                   subplot_title=subplot_title,
                                   subplot_limits=subplot_limits,
                                   has_data=False)

            qgp.do_comparison_plot(dijet_entries_rebin,
                                   "%s/ptBinned/%s_pt%dto%d_dijet_rebin.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim,
                                   subplot_type='ratio',
                                   subplot_title=subplot_title,
                                   subplot_limits=subplot_limits,
                                   has_data=False)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nominal",
                        help="Directory name for nominal files")
    parser.add_argument("--herwig",
                        help="Directory name for nominal Herwig files for comparison",
                        default=None)
    parser.add_argument("--neutralHadronUp",
                        help="Directory name for neutral hadron scale shift up files",
                        default=None)
    parser.add_argument("--neutralHadronDown",
                        help="Directory name for neutral hadron scale shift down files",
                        default=None)
    parser.add_argument("--photonUp",
                        help="Directory name for photon scale shift up files",
                        default=None)
    parser.add_argument("--photonDown",
                        help="Directory name for photon scale shift down files",
                        default=None)
    parser.add_argument("--jesUp",
                        help="Directory name for jes scale shift up files",
                        default=None)
    parser.add_argument("--jesDown",
                        help="Directory name for jes scale shift down files",
                        default=None)
    parser.add_argument("--jerUp",
                        help="Directory name for jer shift up files",
                        default=None)
    parser.add_argument("--jerDown",
                        help="Directory name for jer shift down files",
                        default=None)
    parser.add_argument("--pileupUp",
                        help="Directory name for pileup scale shift up files",
                        default=None)
    parser.add_argument("--pileupDown",
                        help="Directory name for pileup scale shift down files",
                        default=None)
    parser.add_argument("--outputDir", help="Directory for output file")
    args = parser.parse_args()

    if args.neutralHadronUp and args.neutralHadronDown:
        do_plots(args.nominal, args.herwig,
                 os.path.join(args.outputDir, 'neutralHadron'),
                 args.neutralHadronUp, "Neutral hadron energy scale up",
                 args.neutralHadronDown, "Neutral hadron energy scale down")

    if args.photonUp and args.photonDown:
        do_plots(args.nominal, args.herwig,
                 os.path.join(args.outputDir, 'photon'),
                 args.photonUp, "Photon energy scale up",
                 args.photonDown, "Photon energy scale down")

    if args.jesUp and args.jesDown:
        do_plots(args.nominal, args.herwig,
                 os.path.join(args.outputDir, 'jes'),
                 args.jesUp, "Jet energy scale up",
                 args.jesDown, "Jet energy scale down")

    if args.jerUp and args.jerDown:
        do_plots(args.nominal, args.herwig,
                 os.path.join(args.outputDir, 'jer'),
                 args.jerUp, "Jet energy resolution up",
                 args.jerDown, "Jet energy resolution down")

    if args.pileupUp and args.pileupDown:
        do_plots(args.nominal, args.herwig,
                 os.path.join(args.outputDir, 'pileup'),
                 args.pileupUp, "Pileup up",
                 args.pileupDown, "Pileup down")

