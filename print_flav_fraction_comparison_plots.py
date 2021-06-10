#!/usr/bin/env python

"""Print flavour fraction plots/graphs"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)
import argparse
from array import array
import numpy as np
import warnings

# My stuff
from comparator import Contribution, Plot, ZeroContributions
import qg_common as qgc
import qg_general_plots as qgg
import qg_flavour_plots as qgf
import common_utils as cu

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)

# import hunter
# hunter.trace(module='comparator')

# monkey-patch warning formatter
warnings.formatwarning = cu._formatwarning

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Control plot output format
OUTPUT_FMT = "pdf"

MG_SAMPLE = "MG5+Pythia8"
HPP_SAMPLE = "Herwig++"


def do_flav_nparton_fraction_vs_pt(input_file,
                                   dirname,
                                   flav,
                                   n_partons,
                                   pt_bins,
                                   output_filename,
                                   var_prepend='',
                                   which_jet='1',
                                   title=''):

    tfile = cu.open_root_file(input_file)
    hist_2d_total = tfile.Get("%s/%sjet%s_flavour_vs_pt" % (dirname, var_prepend, which_jet))
    hists_2d = [tfile.Get("%s/%sjet%s_flavour_vs_pt_npartons_%s" % (dirname, var_prepend, which_jet, n_parton))
                for n_parton in n_partons]
    # pdigid is x axis, pt is y axis
    # so we want the 1D projection on the y axis,
    # for the x bin corresponding to this flav
    flav_dict = {
        'd': 1,
        'u': 2,
        's': 3,
        'c': 4,
        'b': 5,
        't': 6,
        'g': 21,
        'unknown': 0,
    }
    flav_bin = hists_2d[0].GetXaxis().FindFixBin(flav_dict[flav])
    projs = [h.ProjectionY(cu.get_unique_str(), flav_bin, flav_bin, 'e')
             for h in hists_2d]
    hist_total = hist_2d_total.ProjectionY(cu.get_unique_str(), flav_bin, flav_bin, 'e')

    # rebin hists
    # add (0, x) bin for lowest bin, since hist has that
    bins = np.array([(0., qgc.PT_BINS_ZPJ_INC_UFLOW[0][0])] + qgc.PT_BINS_ZPJ_INC_UFLOW)
    # convert to only have each edge once, not pairs of edges
    bins = np.concatenate((bins[:,0], bins[-1:,1]))
    hist_total = hist_total.Rebin(len(bins)-1, cu.get_unique_str(), bins)
    rebinned_hists = [h.Rebin(len(bins)-1, cu.get_unique_str(), bins) for h in projs]

    # divide each by total to get fraction
    for h in rebinned_hists:
        h.Divide(hist_total)

    # create Contributions for Plot
    colours = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2]
    marker = cu.Marker()
    entries = [Contribution(h, label="%s parton%s" % (n_parton, 's' if n_parton != '1' else ''),
                            line_width=1, line_color=colours[int(n_parton)-1], line_style=1,
                            fill_color=colours[int(n_parton)-1], fill_style=0,
                            marker_size=1, marker_color=colours[int(n_parton)-1], marker_style=mark,
                            leg_draw_opt=None)
               for n_parton, h, mark in zip(n_partons, rebinned_hists, marker.cycle())
               if h.GetEntries() > 0]

    plot = Plot(entries,
                what='hist',
                xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''),
                ytitle="Fraction",
                title=title,
                xlim=(pt_bins[0][0], pt_bins[-1][1]),
                ylim=[0, 1.4],
                has_data=False)
    plot.default_canvas_size = (600, 600)
    plot.plot("HIST E NOSTACK")
    plot.set_logx(do_more_labels=True, do_exponent=False)
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(0.72)
    plot.legend.SetY2(0.85)
    print(output_filename)
    plot.save(output_filename)


def do_all_flavour_fraction_plots(root_dir,
                                  plot_dir="flav_fractions",
                                  zpj_dirname="ZPlusJets_QG",
                                  dj_cen_dirname="Dijet_QG_central_tighter",
                                  dj_fwd_dirname="Dijet_QG_forward_tighter",
                                  do_nparton_plots=True,
                                  var_prepend="",
                                  is_preliminary=True):
    """Do plots of jet flavour fractions vs pT, eta for both Z+jets and dijets regions"""
    # pt_bins = qgc.PT_BINS_INC_UFLOW
    pt_bins = qgc.PT_BINS_ZPJ

    # Plots of all flavour fractions vs pT for a given sample/selection
    # --------------------------------------------------------------------------
    if zpj_dirname:
        # Z+jets
        zpj_histname = qgf.get_flavour_hist_name(dirname=zpj_dirname, var_prepend=var_prepend, which_jet="1")

        if os.path.isfile(os.path.join(root_dir, qgc.DY_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.DY_FILENAME),
                                          hist_name=zpj_histname,
                                          pt_bins=pt_bins,
                                          title=qgc.ZpJ_LABEL + "\n" + MG_SAMPLE,
                                          output_filename="%s/zpj_flavour_fractions.%s" % (plot_dir, OUTPUT_FMT))

        if os.path.isfile(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.DY_HERWIG_FILENAME),
                                          hist_name=zpj_histname,
                                          pt_bins=pt_bins,
                                          title=qgc.ZpJ_LABEL + "\n" + HPP_SAMPLE,
                                          output_filename="%s/zpj_flavour_fractions_herwigpp.%s" % (plot_dir, OUTPUT_FMT))
    if dj_cen_dirname:
        # Dijets central
        dijet_cen_histname = qgf.get_flavour_hist_name(dirname=dj_cen_dirname, var_prepend=var_prepend, which_jet="1")

        if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                          hist_name=dijet_cen_histname,
                                          pt_bins=pt_bins,
                                          title=qgc.Dijet_CEN_LABEL + "\n" + MG_SAMPLE,
                                          output_filename="%s/dj_flavour_fractions_central_jet.%s" % (plot_dir, OUTPUT_FMT))

        if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME),
                                          hist_name=dijet_cen_histname,
                                          pt_bins=pt_bins,
                                          title=qgc.Dijet_CEN_LABEL + "\n" + HPP_SAMPLE,
                                          output_filename="%s/dj_flavour_fractions_central_jet_herwigpp.%s" % (plot_dir, OUTPUT_FMT))

    if dj_fwd_dirname:
        # Dijets central
        dijet_fwd_histname = qgf.get_flavour_hist_name(dirname=dj_fwd_dirname, var_prepend=var_prepend, which_jet="1")

        if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                          hist_name=dijet_fwd_histname,
                                          pt_bins=pt_bins,
                                          title=qgc.Dijet_FWD_LABEL + "\n" + MG_SAMPLE,
                                          output_filename="%s/dj_flavour_fractions_forward_jet.%s" % (plot_dir, OUTPUT_FMT))

        if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME),
                                          hist_name=dijet_fwd_histname,
                                          pt_bins=pt_bins,
                                          title=qgc.Dijet_FWD_LABEL + "\n" + HPP_SAMPLE,
                                          output_filename="%s/dj_flavour_fractions_forward_jet_herwigpp.%s" % (plot_dir, OUTPUT_FMT))

    # Plots of all flavour fractions vs eta for a given sample/selection
    # --------------------------------------------------------------------------
    # Do for separate pt regions
    pt_regions = [
        {
            'append': "_lowPt",
            'title': '30 < p_{T}^{jet} < 100 GeV',
        },
        {
            'append': "_midPt",
            'title': '100 < p_{T}^{jet} < 250 GeV',
        },
        {
            'append': "_highPt",
            'title': '250 < p_{T}^{jet} < 500 GeV',
        },
        {
            'append': "_highPt2",
            'title': 'p_{T}^{jet} > 500 GeV',
        },
    ]
    end = 1.7
    interval = 0.25
    eta_bins = np.arange(-end, end+interval, interval)
    eta_bins = list(zip(eta_bins[:-1], eta_bins[1:])) # make pairwise bins
    # print(eta_bins)
    try:
        for pt_region in pt_regions:
            eta_title = pt_region['title']
            append = pt_region['append']

            if zpj_dirname:
                # Z+jets
                if os.path.isfile(os.path.join(root_dir, qgc.DY_FILENAME)):
                    qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.DY_FILENAME),
                                                   title=qgc.ZpJ_LABEL +  "\n" + eta_title +  "\n" + MG_SAMPLE,
                                                   dirname=zpj_dirname,
                                                   eta_bins=eta_bins,
                                                   var_prepend=var_prepend,
                                                   which_jet="both",
                                                   append=append,
                                                   output_filename="%s/zpj_flavour_fractions_vs_eta%s.%s" % (plot_dir, append, OUTPUT_FMT))

                if os.path.isfile(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME)):
                    qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.DY_HERWIG_FILENAME),
                                                   title=qgc.ZpJ_LABEL +  "\n" + eta_title + "\n" + HPP_SAMPLE,
                                                   dirname=zpj_dirname,
                                                   eta_bins=eta_bins,
                                                   var_prepend=var_prepend,
                                                   which_jet="both",
                                                   append=append,
                                                   output_filename="%s/zpj_flavour_fractions_vs_eta%s_herwigpp.%s" % (plot_dir, append, OUTPUT_FMT))

            if dj_cen_dirname:
                # Dijets central
                if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
                    qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                                   title=qgc.Dijet_CEN_LABEL +  "\n" + eta_title +  "\n" + MG_SAMPLE,
                                                   dirname=dj_cen_dirname,
                                                   eta_bins=eta_bins,
                                                   var_prepend=var_prepend,
                                                   which_jet="both",
                                                   append=append,
                                                   output_filename="%s/dj_flavour_fractions_central_jet_vs_eta%s.%s" % (plot_dir, append, OUTPUT_FMT))

                if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
                    qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME),
                                                   title=qgc.Dijet_CEN_LABEL +  "\n" + eta_title + "\n" + HPP_SAMPLE,
                                                   dirname=dj_cen_dirname,
                                                   eta_bins=eta_bins,
                                                   var_prepend=var_prepend,
                                                   which_jet="both",
                                                   append=append,
                                                   output_filename="%s/dj_flavour_fractions_central_jet_vs_eta%s_herwigpp.%s" % (plot_dir, append, OUTPUT_FMT))

            if dj_fwd_dirname:
                # Dijets central
                if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
                    qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                                   title=qgc.Dijet_FWD_LABEL +  "\n" + eta_title +  "\n" + MG_SAMPLE,
                                                   dirname=dj_fwd_dirname,
                                                   eta_bins=eta_bins,
                                                   var_prepend=var_prepend,
                                                   which_jet="both",
                                                   append=append,
                                                   output_filename="%s/dj_flavour_fractions_forward_jet_vs_eta%s.%s" % (plot_dir, append, OUTPUT_FMT))

                if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
                    qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME),
                                                   title=qgc.Dijet_FWD_LABEL +  "\n" + eta_title + "\n" + HPP_SAMPLE,
                                                   dirname=dj_fwd_dirname,
                                                   eta_bins=eta_bins,
                                                   var_prepend=var_prepend,
                                                   which_jet="both",
                                                   append=append,
                                                   output_filename="%s/dj_flavour_fractions_forward_jet_vs_eta%s_herwigpp.%s" % (plot_dir, append, OUTPUT_FMT))
    except ZeroContributions as e:
        print("Skipping flav vs eta plots")

    # Compare one flav fraction across samples/selections
    # --------------------------------------------------------------------------
    if dj_fwd_dirname and zpj_dirname:
        dirnames = [dj_cen_dirname, dj_fwd_dirname, zpj_dirname]
        for this_flav in ['g', 'u', 'd', 's', 'c', 'b'][:]:
            input_files = [os.path.join(root_dir, qgc.QCD_FILENAME) if "dijet" in d.lower()
                           else os.path.join(root_dir, qgc.DY_FILENAME)
                           for d in dirnames]

            labels = [qgc.Dijet_CEN_LABEL, qgc.Dijet_FWD_LABEL, qgc.ZpJ_LABEL]

            append = "_paper" if not is_preliminary else ""
            qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                dirnames=dirnames,
                                                pt_bins=pt_bins,
                                                labels=labels,
                                                flav=this_flav,
                                                output_filename="%s/%s_flav_fraction_compare_jet1%s.%s" % (plot_dir, this_flav, append, OUTPUT_FMT),
                                                var_prepend=var_prepend,
                                                which_jet="1",
                                                xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''),
                                                is_preliminary=is_preliminary)

            if do_nparton_plots:
                # Do for various n partons exclusively
                for n_partons in ["1", "2", "3", "4"]:
                    qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                        dirnames=dirnames,
                                                        pt_bins=pt_bins,
                                                        labels=labels,
                                                        flav=this_flav,
                                                        output_filename="%s/%s_flav_fraction_compare_jet1_npartons_%s.%s" % (plot_dir, this_flav, n_partons, OUTPUT_FMT),
                                                        n_partons=n_partons,
                                                        var_prepend=var_prepend,
                                                        which_jet="1",
                                                        title=n_partons+" outgoing parton"+("s" if n_partons != "1" else ""),
                                                        xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''))

    dirnames, labels = [], []
    if dj_cen_dirname:
        dirnames += [dj_cen_dirname, dj_fwd_dirname]
        labels += [qgc.Dijet_CEN_LABEL, qgc.Dijet_FWD_LABEL]
    if zpj_dirname:
        dirnames += [zpj_dirname]
        labels += [qgc.ZpJ_LABEL]

    # Compare MG+Pythia vs Herwig++ for each region, flav
    for dname, label in zip(dirnames, labels):
        input_files = [os.path.join(root_dir, qgc.DY_FILENAME), os.path.join(root_dir, qgc.DY_HERWIG_FILENAME)]
        if 'dijet' in dname.lower():
            input_files = [os.path.join(root_dir, qgc.QCD_FILENAME), os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)]

        sample_labels = [MG_SAMPLE, HPP_SAMPLE]
        sample_short_names = ['mgpythia', 'herwigpp']

        short_name = "dj_central"
        if "forward" in dname.lower():
            short_name = "dj_forward"
        if "zplusjet" in dname.lower():
            short_name = "zpj"

        for this_flav in ['g', 'u', 'd', 's', 'c', 'b'][:]:

            qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                dirnames=[dname for i in input_files],
                                                pt_bins=pt_bins,
                                                labels=sample_labels,
                                                flav=this_flav,
                                                output_filename="%s/%s_%s_flav_fraction_compare_mgpythia_vs_herwigpp_jet1.%s" % (plot_dir, this_flav, short_name, OUTPUT_FMT),
                                                var_prepend=var_prepend,
                                                which_jet="1",
                                                title=label,
                                                xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''))

            if do_nparton_plots:
                # Do for various n partons exclusively
                all_n_partons = ["2", "3", "4"] if 'dijet' in dname.lower() else ["1", "2", "3", "4"]
                for n_partons in all_n_partons:
                    qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                        dirnames=[dname for i in input_files],
                                                        pt_bins=pt_bins,
                                                        labels=sample_labels,
                                                        flav=this_flav,
                                                        output_filename="%s/%s_%s_flav_fraction_compare_mgpythia_vs_herwigpp_jet1_npartons_%s.%s" % (plot_dir, this_flav, short_name, n_partons, OUTPUT_FMT),
                                                        n_partons=n_partons,
                                                        var_prepend=var_prepend,
                                                        which_jet="1",
                                                        title=label + "\n" + n_partons+" outgoing parton"+("s" if n_partons != "1" else ""),
                                                        xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''))

                # Do all n partons on one plot, for brevity
                qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                    dirnames=[dname for i in input_files],
                                                    pt_bins=pt_bins,
                                                    labels=sample_labels,
                                                    flav=this_flav,
                                                    output_filename="%s/%s_%s_flav_fraction_compare_mgpythia_vs_herwigpp_jet1_npartons_all.%s" % (plot_dir, this_flav, short_name, OUTPUT_FMT),
                                                    n_partons=all_n_partons,
                                                    var_prepend=var_prepend,
                                                    which_jet="1",
                                                    title=label + "\n%s outgoing partons" % (", ".join(["%s-" % n for n in all_n_partons])),
                                                    xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''))


                # Do a plot of fractions of n partons for this flavour / input etc
                for input_file, sample_label, sample_short_name in zip(input_files, sample_labels, sample_short_names):
                    do_flav_nparton_fraction_vs_pt(input_file=input_file,
                                                   dirname=dname,
                                                   flav=this_flav,
                                                   n_partons=['2', '3', '4'] if 'dijet' in dname.lower() else ['1', '2', '3', '4'],
                                                   pt_bins=pt_bins,
                                                   var_prepend=var_prepend,
                                                   which_jet='1',
                                                   title=label + "\n" + sample_label + '\n%s flavour jets' % this_flav,
                                                   output_filename="%s/%s_%s_%s_nparton_fraction.%s" % (plot_dir, this_flav, short_name, sample_short_name, OUTPUT_FMT))



def do_flavour_fraction_input_comparison_plots(root_dirs,
                                               labels,
                                               plot_dir="flav_fractions_comparison",
                                               zpj_dirname="ZPlusJets_QG",
                                               dj_cen_dirname="Dijet_QG_central_tighter",
                                               dj_fwd_dirname="Dijet_QG_forward_tighter",
                                               var_prepend=""):
    """Do plots comparing several input dirs """
    # pt_bins = qgc.PT_BINS_INC_UFLOW
    pt_bins = qgc.PT_BINS_ZPJ
    if dj_cen_dirname:
        this_flav = "g"
        dirnames = [dj_cen_dirname]*len(root_dirs)
        # Compare gluon fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/dj_g_flav_fraction_compare_central_jet.%s" % (plot_dir, OUTPUT_FMT),
                                            var_prepend=var_prepend,
                                            which_jet="1",
                                            title=qgc.Dijet_CEN_LABEL + "\nMG+PYTHIA8",
                                            xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''))
        # this_flav = "1-g"
        # dirnames = [dj_cen_dirname]*len(root_dirs)
        # # Compare non-gluon fractions across samples/selections
        # qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
        #                                     dirnames=dirnames,
        #                                     pt_bins=pt_bins,
        #                                     labels=labels,
        #                                     flav=this_flav,
        #                                     output_filename="%s/dj_q_flav_fraction_compare_central_jet.%s" % (plot_dir, OUTPUT_FMT),
        #                                     var_prepend=var_prepend,
        #                                     which_jet="1",
        #                                     xtitle="p_{T}^{jet 1} [GeV]")
    if dj_fwd_dirname:
        this_flav = "g"
        dirnames = [dj_fwd_dirname]*len(root_dirs)
        # Compare gluon fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/dj_g_flav_fraction_compare_forward_jet.%s" % (plot_dir, OUTPUT_FMT),
                                            var_prepend=var_prepend,
                                            which_jet="1",
                                            title=qgc.Dijet_FWD_LABEL + "\nMG+PYTHIA8",
                                            xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''))
        # this_flav = "1-g"
        # dirnames = [dj_fwd_dirname]*len(root_dirs)
        # # Compare non-gluon fractions across samples/selections
        # qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
        #                                     dirnames=dirnames,
        #                                     pt_bins=pt_bins,
        #                                     labels=labels,
        #                                     flav=this_flav,
        #                                     output_filename="%s/dj_q_flav_fraction_compare_forward_jet.%s" % (plot_dir, OUTPUT_FMT),
        #                                     var_prepend=var_prepend,
        #                                     which_jet="1",
        #                                     xtitle="p_{T}^{jet 1} [GeV]")

    if zpj_dirname:
        this_flav = "g"
        dirnames = [zpj_dirname] * len(root_dirs)
        # Compare gluon fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.DY_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/zpj_g_flav_fraction_compare.%s" % (plot_dir, OUTPUT_FMT),
                                            title=qgc.ZpJ_LABEL + "\nMG+PYTHIA8",
                                            xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''),
                                            var_prepend=var_prepend)

        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.DY_HERWIG_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/zpj_hpp_g_flav_fraction_compare.%s" % (plot_dir, OUTPUT_FMT),
                                            title=qgc.ZpJ_LABEL + "\nHERWIG++",
                                            xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''),
                                            var_prepend=var_prepend)

        # this_flav = "1-g"
        # dirnames = [zpj_dirname] * len(root_dirs)
        # # Compare quark fractions across samples/selections
        # qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.DY_FILENAME) for rd in root_dirs],
        #                                     dirnames=dirnames,
        #                                     pt_bins=pt_bins,
        #                                     labels=labels,
        #                                     flav=this_flav,
        #                                     output_filename="%s/zpj_q_flav_fraction_compare.%s" % (plot_dir, OUTPUT_FMT),
        #                                     var_prepend=var_prepend)


def do_flavour_fraction_plots_with_rivet(plot_dir="flav_fractions_with_rivet",
                                         rivet_dijet_files=None,
                                         rivet_dijet_labels=None,
                                         rivet_zpj_files=None,
                                         rivet_zpj_labels=None):
    """Do plots of jet flavour fractions vs pT for both Z+jets and dijets regions,
    also include RIVET inputs"""
    pt_bins = qgc.PT_BINS_ZPJ

    # Plots of all flavour fractions vs pT for a given sample/selection
    # --------------------------------------------------------------------------
    if rivet_zpj_files:
        # Rivet Z+J files
        if len(rivet_zpj_files) != len(rivet_zpj_labels):
            raise ValueError("Expect len(rivet_zpj_files) == len(rivet_zpj_labels)!")

        for radius_ind, radius in enumerate(['AK4', 'AK8'], 1):
            zpj_rivet_histname = qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_ZPJ', region='zpj', radius_ind=radius_ind)

            for fname, label in zip(rivet_zpj_files, rivet_zpj_labels):
                sample_str = os.path.basename(fname).replace(".root", "")
                qgf.do_flavour_fraction_vs_pt(input_file=fname,
                                              hist_name=zpj_rivet_histname,
                                              title=qgc.ZpJ_LABEL + "\n" + label + "\n %s jets" % radius,
                                              pt_bins=pt_bins,
                                              output_filename="%s/zpj_flavour_fractions_%s_%s.%s" % (plot_dir, sample_str, radius, OUTPUT_FMT))

    if rivet_dijet_files:
        # Dijet files
        if len(rivet_dijet_files) != len(rivet_dijet_labels):
            raise ValueError("Expect len(rivet_dijet_files) == len(rivet_dijet_labels)!")

        for radius_ind, radius in enumerate(['AK4', 'AK8'], 1):
            dijet_cen_rivet_histname = qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_DIJET', region='central_jet', radius_ind=radius_ind)
            dijet_fwd_rivet_histname = qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_DIJET', region='forward_jet', radius_ind=radius_ind)

            for fname, label in zip(rivet_dijet_files, rivet_dijet_labels):
                sample_str = os.path.basename(fname).replace(".root", "")
                qgf.do_flavour_fraction_vs_pt(input_file=fname,
                                              hist_name=dijet_cen_rivet_histname,
                                              title=qgc.Dijet_CEN_LABEL + "\n" + label + "\n%s jets" % radius,
                                              pt_bins=pt_bins,
                                              output_filename="%s/dj_flavour_fractions_central_jet_%s_%s.%s" % (plot_dir, sample_str, radius, OUTPUT_FMT))

                qgf.do_flavour_fraction_vs_pt(input_file=fname,
                                              hist_name=dijet_fwd_rivet_histname,
                                              title=qgc.Dijet_FWD_LABEL + "\n" + label + "\n%s jets" % radius,
                                              pt_bins=pt_bins,
                                              output_filename="%s/dj_flavour_fractions_forward_jet_%s_%s.%s" % (plot_dir, sample_str, radius, OUTPUT_FMT))

    # Compare one flav fraction across regions/selections
    # --------------------------------------------------------------------------
    if rivet_dijet_files and rivet_zpj_files:
        # relies on order being same between dijet and z+j files!
        for dj_fname, dj_label, zpj_fname, zpj_label in zip(rivet_dijet_files, rivet_dijet_labels, rivet_zpj_files, rivet_zpj_labels):

            if dj_label != zpj_label:
                warnings.warn("Mismatched labels: %s vs %s" % (dj_label, zpj_label))

            for radius_ind, radius in enumerate(['AK4', 'AK8'], 1):
                for this_flav in ['g', 'u', 'd', 's', 'c', 'b'][:]:
                    input_files = [dj_fname, dj_fname, zpj_fname]

                    hist_names = [
                        qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_DIJET', region='central_jet', radius_ind=radius_ind),
                        qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_DIJET', region='forward_jet', radius_ind=radius_ind),
                        qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_ZPJ', region='zpj', radius_ind=radius_ind)
                    ]

                    labels = [qgc.Dijet_CEN_LABEL, qgc.Dijet_FWD_LABEL, qgc.ZpJ_LABEL]

                    sample_str = os.path.basename(dj_fname).replace(".root", "").replace("QCD_", "")

                    qgf.compare_flavour_fraction_hists_vs_pt(input_files=input_files,
                                                             hist_names=hist_names,
                                                             pt_bins=pt_bins,
                                                             labels=labels,
                                                             flav=this_flav,
                                                             title=dj_label + "\n%s jets" % radius,
                                                             output_filename="%s/%s_flav_fraction_compare_%s_%s.%s" % (plot_dir, this_flav, sample_str, radius, OUTPUT_FMT),
                                                             xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(''))


def compare_flavour_fraction_plots_with_rivet(root_dir,
                                              root_dir_label,
                                              zpj_dirname="ZPlusJets_QG",
                                              dj_cen_dirname="Dijet_QG_central_tighter",
                                              dj_fwd_dirname="Dijet_QG_forward_tighter",
                                              var_prepend="",
                                              plot_dir="flav_fractions_with_rivet_comparison",
                                              rivet_dijet_files=None,
                                              rivet_dijet_labels=None,
                                              rivet_zpj_files=None,
                                              rivet_zpj_labels=None,
                                              is_preliminary=True):
    """Do flav fraction plots comparing rivet & non-rivet inputs"""
    pt_bins = qgc.PT_BINS_ZPJ

    radius_ind = 1 if '_ak4puppi_' in root_dir else 2
    radius_str = 'AK4' if radius_ind == 1 else 'AK8'

    plot_kwargs = dict(
        is_preliminary=is_preliminary,
        is_supplementary=True
    )

    append = "_paper" if not is_preliminary else ""

    # Z+jets
    if zpj_dirname or rivet_zpj_files:
        zpj_histname = qgf.get_flavour_hist_name(dirname=zpj_dirname, var_prepend=var_prepend, which_jet="1")
        zpj_rivet_histname = qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_ZPJ', region='zpj', radius_ind=radius_ind)

        for this_flav in ['g', 'u', 'd', 's', 'c', 'b'][:]:
                non_rivet_files = []
                non_rivet_histnames = []
                non_rivet_labels = []
                if os.path.isfile(os.path.join(root_dir, qgc.DY_FILENAME)):
                    non_rivet_files.append(os.path.join(root_dir, qgc.DY_FILENAME))
                    non_rivet_histnames.append(zpj_histname)
                    non_rivet_labels.append(MG_SAMPLE)
                if os.path.isfile(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME)):
                    non_rivet_files.append(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME))
                    non_rivet_histnames.append(zpj_histname)
                    non_rivet_labels.append(HPP_SAMPLE)

                contribs = qgf.create_contibutions_compare_vs_pt(input_files=[*non_rivet_files, *rivet_zpj_files],
                                                                 hist_names=[*non_rivet_histnames, *[zpj_rivet_histname for r in rivet_zpj_files]],
                                                                 pt_bins=pt_bins,
                                                                 labels=[*non_rivet_labels, *rivet_zpj_labels],
                                                                 flav=this_flav)

                output_filename="%s/zpj_%s_flavour_fractions_compare_rivet_%s%s.%s" % (plot_dir, this_flav, radius_str, append, OUTPUT_FMT)
                qgf.compare_flavour_fraction_hists_vs_pt_from_contribs(contribs,
                                                                       flav=this_flav,
                                                                       output_filename=output_filename,
                                                                       title=qgc.ZpJ_LABEL + "\n%s jets" % radius_str,
                                                                       **plot_kwargs)

    # Central dijet
    if dj_cen_dirname or rivet_dijet_files:
        dj_cen_histname = qgf.get_flavour_hist_name(dirname=dj_cen_dirname, var_prepend=var_prepend, which_jet="1")
        dijet_cen_rivet_histname = qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_DIJET', region='central_jet', radius_ind=radius_ind)

        for this_flav in ['g', 'u', 'd', 's', 'c', 'b'][:]:
                non_rivet_files = []
                non_rivet_histnames = []
                non_rivet_labels = []
                if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
                    non_rivet_files.append(os.path.join(root_dir, qgc.QCD_FILENAME))
                    non_rivet_histnames.append(dj_cen_histname)
                    non_rivet_labels.append(MG_SAMPLE)

                if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
                    non_rivet_files.append(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME))
                    non_rivet_histnames.append(dj_cen_histname)
                    non_rivet_labels.append(HPP_SAMPLE)

                contribs = qgf.create_contibutions_compare_vs_pt(input_files=[*non_rivet_files, *rivet_dijet_files],
                                                                 hist_names=[*non_rivet_histnames, *[dijet_cen_rivet_histname for r in rivet_dijet_files]],
                                                                 pt_bins=pt_bins,
                                                                 labels=[*non_rivet_labels, *rivet_dijet_labels],
                                                                 flav=this_flav)

                output_filename="%s/dj_central_%s_flavour_fractions_compare_rivet_%s%s.%s" % (plot_dir, this_flav, radius_str, append, OUTPUT_FMT)
                qgf.compare_flavour_fraction_hists_vs_pt_from_contribs(contribs,
                                                                       flav=this_flav,
                                                                       output_filename=output_filename,
                                                                       title=qgc.Dijet_CEN_LABEL + "\n%s jets" % radius_str,
                                                                       **plot_kwargs)

    # Forward dijet
    if dj_fwd_dirname and rivet_dijet_files:
        dj_fwd_histname = qgf.get_flavour_hist_name(dirname=dj_fwd_dirname, var_prepend=var_prepend, which_jet="1")
        dijet_fwd_rivet_histname = qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_DIJET', region='forward_jet', radius_ind=radius_ind)

        for this_flav in ['g', 'u', 'd', 's', 'c', 'b'][:]:
                non_rivet_files = []
                non_rivet_histnames = []
                non_rivet_labels = []
                if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
                    non_rivet_files.append(os.path.join(root_dir, qgc.QCD_FILENAME))
                    non_rivet_histnames.append(dj_fwd_histname)
                    non_rivet_labels.append(MG_SAMPLE)
                if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
                    non_rivet_files.append(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME))
                    non_rivet_histnames.append(dj_fwd_histname)
                    non_rivet_labels.append(HPP_SAMPLE)

                contribs = qgf.create_contibutions_compare_vs_pt(input_files=[*non_rivet_files, *rivet_dijet_files],
                                                                 hist_names=[*non_rivet_histnames, *[dijet_fwd_rivet_histname for r in rivet_dijet_files]],
                                                                 pt_bins=pt_bins,
                                                                 labels=[*non_rivet_labels, *rivet_dijet_labels],
                                                                 flav=this_flav)

                output_filename="%s/dj_forward_%s_flavour_fractions_compare_rivet_%s%s.%s" % (plot_dir, this_flav, radius_str, append, OUTPUT_FMT)
                qgf.compare_flavour_fraction_hists_vs_pt_from_contribs(contribs,
                                                                       flav=this_flav,
                                                                       output_filename=output_filename,
                                                                       title=qgc.Dijet_FWD_LABEL + "\n%s jets" % radius_str,
                                                                       **plot_kwargs)

    # Do g ratio in dijet vs ZJet
    if all([rivet_dijet_files, rivet_zpj_files, dj_cen_dirname, zpj_dirname]):
        # Central dijet
        dj_cen_histname = qgf.get_flavour_hist_name(dirname=dj_cen_dirname, var_prepend=var_prepend, which_jet="1")
        dijet_cen_rivet_histname = qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_DIJET', region='central_jet', radius_ind=radius_ind)
        # Z+jets
        zpj_histname = qgf.get_flavour_hist_name(dirname=zpj_dirname, var_prepend=var_prepend, which_jet="1")
        zpj_rivet_histname = qgf.get_rivet_flavour_hist_name(dirname='CMS_2018_PAS_SMP_18_QGX_ZPJ', region='zpj', radius_ind=radius_ind)

        for this_flav in ['g', 'u', 'd', 's', 'c', 'b'][:]:
            non_rivet_files = []
            non_rivet_histnames = []
            non_rivet_labels = []
            if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
                non_rivet_files.append(os.path.join(root_dir, qgc.QCD_FILENAME))
                non_rivet_histnames.append(dj_cen_histname)
                non_rivet_labels.append(MG_SAMPLE)

            if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
                non_rivet_files.append(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME))
                non_rivet_histnames.append(dj_cen_histname)
                non_rivet_labels.append(HPP_SAMPLE)

            dj_cen_contribs = qgf.create_contibutions_compare_vs_pt(input_files=[*non_rivet_files, *rivet_dijet_files],
                                                                    hist_names=[*non_rivet_histnames, *[dijet_cen_rivet_histname for r in rivet_dijet_files]],
                                                                    pt_bins=pt_bins,
                                                                    labels=[*non_rivet_labels, *rivet_dijet_labels],
                                                                    flav=this_flav)

            non_rivet_files = []
            non_rivet_histnames = []
            non_rivet_labels = []
            if os.path.isfile(os.path.join(root_dir, qgc.DY_FILENAME)):
                non_rivet_files.append(os.path.join(root_dir, qgc.DY_FILENAME))
                non_rivet_histnames.append(zpj_histname)
                non_rivet_labels.append(MG_SAMPLE)
            if os.path.isfile(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME)):
                non_rivet_files.append(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME))
                non_rivet_histnames.append(zpj_histname)
                non_rivet_labels.append(HPP_SAMPLE)

            zpj_contribs = qgf.create_contibutions_compare_vs_pt(input_files=[*non_rivet_files, *rivet_zpj_files],
                                                                 hist_names=[*non_rivet_histnames, *[zpj_rivet_histname for r in rivet_zpj_files]],
                                                                 pt_bins=pt_bins,
                                                                 labels=[*non_rivet_labels, *rivet_zpj_labels],
                                                                 flav=this_flav)

            # Now create ratio hists from each contrib pair
            if len(dj_cen_contribs) != len(zpj_contribs):
                raise ValueError('len(dj_cen_contribs) != len(zpj_contribs)')

            ratio_contribs = []
            for dj_contrib, zpj_contrib in zip(dj_cen_contribs, zpj_contribs):
                h_ratio = cu.tgraph_to_th1(dj_contrib.obj)
                h_ratio.Divide(cu.tgraph_to_th1(zpj_contrib.obj))
                ratio_contrib = Contribution(h_ratio,
                                             label=dj_contrib.label,
                                             line_color=dj_contrib.line_color,
                                             line_width=dj_contrib.line_width,
                                             line_style=dj_contrib.line_style,
                                             marker_style=dj_contrib.marker_style,
                                             marker_color=dj_contrib.marker_color,
                                             marker_size=dj_contrib.marker_size,
                                             leg_draw_opt=dj_contrib.leg_draw_opt,)
                ratio_contribs.append(ratio_contrib)

            flav_str = qgf.FLAV_STR_DICT[this_flav]
            ytitle = "%s : %s" % (qgc.Dijet_CEN_LABEL, qgc.ZpJ_LABEL)
            p = Plot(ratio_contribs,
                     what='hist',
                     xtitle="p_{T}^{jet} [GeV]",
                     ytitle=ytitle,
                     title="%s flavour\n%s jets" % (flav_str, radius_str),
                     xlim=(50, 2000),
                     # ylim=(0, 1),
                     has_data=False,
                     **plot_kwargs)
            p.default_canvas_size = (600, 600)

            output_filename="%s/dj_central_vs_zpj_%s_flavour_fractions_compare_rivet_%s%s.%s" % (plot_dir, this_flav, radius_str, append, OUTPUT_FMT)
            try:
                p.plot("E NOSTACK")
                # p.main_pad.SetLeftMargin(0.18)
                # p.get_modifier().GetXaxis().SetTitleOffset(1.4)
                # p.get_modifier().GetXaxis().SetTitleSize(.045)
                # p.get_modifier().GetYaxis().SetTitleSize(.025)
                p.legend.SetX1(0.56)
                p.legend.SetY1(0.65)
                p.legend.SetY2(0.87)
                if len(ratio_contribs) >=4:
                    p.legend.SetY1(0.7)
                    p.legend.SetX1(0.5)
                    p.legend.SetNColumns(2)
                p.set_logx(do_more_labels=True, do_exponent=False)
                print(output_filename)
                p.save(output_filename)
            except ZeroContributions as e:
                warnings.warn("No contributions for %s" % output_filename)
                print(e.message)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dj",
                        help="Do dijet plots",
                        action="store_true")
    parser.add_argument("--zpj",
                        help="Do Z+jet plots",
                        action="store_true")
    parser.add_argument("--gen",
                        help="Use genjet flavour",
                        action="store_true")
    parser.add_argument("--dir",
                        help="Directory to get plot from. Can be used multiple times",
                        action="append")
    parser.add_argument("--dirLabel",
                        help="Label to be given for dir. Must be used in conjunction with --dir, once per entry.",
                        action="append")
    parser.add_argument("--rivetdj",
                        help="Rivet dijet ROOT file (use yoda2root to convert)",
                        action='append')
    parser.add_argument("--rivetdjLabel",
                        help="Rivet dijet ROOT file label. Must be used in conjunction with --rivetdjLabel, once per entry.",
                        action='append')
    parser.add_argument("--rivetzpj",
                        help="Rivet Z+jet ROOT file (use yoda2root to convert)",
                        action='append')
    parser.add_argument("--rivetzpjLabel",
                        help="Rivet Z+jet ROOT file label. Must be used in conjunction with --rivetzpjLabel, once per entry.",
                        action='append')

    parser.add_argument("--paper",
                        action='store_true',
                        help='Don\'t add "Preliminary" to plots')
    args = parser.parse_args()
    print(args)

    if args.rivetdj or args.rivetzpj:
        do_flavour_fraction_plots_with_rivet(plot_dir=os.path.join("flav_fractions_with_rivet"),
                                             rivet_dijet_files=args.rivetdj,
                                             rivet_dijet_labels=args.rivetdjLabel,
                                             rivet_zpj_files=args.rivetzpj,
                                             rivet_zpj_labels=args.rivetzpjLabel)

        if args.dir and (args.dj or args.zpj):
            compare_flavour_fraction_plots_with_rivet(root_dir=args.dir[0],
                                                      root_dir_label=args.dirLabel[0],
                                                      var_prepend="gen" if args.gen else "",
                                                      zpj_dirname="ZPlusJets_QG" if args.zpj else None,
                                                      dj_cen_dirname="Dijet_QG_central_tighter" if args.dj else None,
                                                      dj_fwd_dirname="Dijet_QG_forward_tighter" if args.dj else None,
                                                      plot_dir=os.path.join("flav_fractions_with_rivet_compare"),
                                                      rivet_dijet_files=args.rivetdj,
                                                      rivet_dijet_labels=args.rivetdjLabel,
                                                      rivet_zpj_files=args.rivetzpj,
                                                      rivet_zpj_labels=args.rivetzpjLabel,
                                                      is_preliminary=not args.paper)

    # One set of plots per input
    for adir in args.dir:
        do_all_flavour_fraction_plots(adir,
                                      plot_dir=os.path.join(adir, "flav_fractions%s" % ("_gen" if args.gen else "")),
                                      var_prepend="gen" if args.gen else "",
                                      zpj_dirname="ZPlusJets_QG" if args.zpj else None,
                                      dj_cen_dirname="Dijet_QG_central_tighter" if args.dj else None,
                                      dj_fwd_dirname="Dijet_QG_forward_tighter" if args.dj else None,
                                      is_preliminary=not args.paper)

        do_all_flavour_fraction_plots(adir,
                                      plot_dir=os.path.join(adir, "flav_fractions_gen_selection"),
                                      var_prepend="",
                                      zpj_dirname="ZPlusJets_QG_gen" if args.zpj else None,
                                      dj_cen_dirname="Dijet_QG_gen_central" if args.dj else None,
                                      dj_fwd_dirname="Dijet_QG_gen_forward" if args.dj else None,
                                      is_preliminary=not args.paper)

    if len(args.dir) > 1:
        # Now comparison across all inputs
        do_flavour_fraction_input_comparison_plots(args.dir,
                                                   plot_dir=os.path.join(args.dir[0], "flav_fractions_comparison%s" % ("_gen" if args.gen else "")),
                                                   labels=args.dirLabel,
                                                   var_prepend="gen" if args.gen else "",
                                                   dj_cen_dirname="Dijet_QG_central_tighter" if args.dj else None,
                                                   dj_fwd_dirname="Dijet_QG_forward_tighter" if args.dj else None,
                                                   zpj_dirname="ZPlusJets_QG" if args.zpj else None)
