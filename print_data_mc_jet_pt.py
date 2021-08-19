#!/usr/bin/env python

"""
Do jet pT distributions for paper plots
"""

import os
os.nice(10)
import ROOT
from MyStyle import My_Style
My_Style.cd()
import sys
from array import array
import numpy as np
import math
import argparse

# My stuff
from comparator import Contribution, Plot
import qg_common as qgc
import qg_general_plots as qgg
import common_utils as cu

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.TH1.AddDirectory(False)  # VERY IMPORTANT - somewhere, closing a TFile for exp systs deletes a map...dunno why


# Control output format
OUTPUT_FMT = "pdf"


def compare_bins(histA, histB):
    h1Array = histA.GetXaxis().GetXbins()
    h2Array = histB.GetXaxis().GetXbins()
    if h1Array.fN != h2Array.fN:
        raise ValueError("h1Array.fN != h2Array.fN")

    for i in range(h1Array.fN):
        if not ROOT.TMath.AreEqualRel(h1Array.GetAt(i), h2Array.GetAt(i), 1E-10):
            raise ValueError("not equal: %g %g" % (h1Array.GetAt(i), h2Array.GetAt(i)))

def do_jet_pt_plot(entries,
                   output_filename,
                   xlim=None,
                   ylim=None,
                   title="",
                   subplot_title="",
                   rebin=1,
                   data_first=True,
                   normalise_hists=True,
                   subplot_limits=None,
                   experimental_syst=None,
                   scale_syst=None,
                   pdf_syst=None,
                   total_syst=None,
                   lumi=cu.get_lumi_str(do_dijet=False, do_zpj=True),
                   is_preliminary=True):
    # entries = [ent for ent in entries if ent]

    conts = [Contribution(ent[0], rebin_hist=rebin, **ent[1])  # Don't use noramlise_hist - screws up with unequal binning
             for ent in entries]

    do_legend = len(conts) > 1
    if len(conts) == 0:
        raise RuntimeError("0 contributions for this plot")

    if data_first:
        # Modify data Contribution to have no horizontal lines
        # We'll redraw it at the end with the proper bars ontop
        conts[0].line_width = 0
        conts[0].marker_size = 0.01
        conts[0].update_obj_styling()

        # For subplot to ensure only MC errors drawn, not MC+data
        data_obj = entries[0][0]
        # data_obj = conts[0].obj
        data_no_errors = data_obj.Clone()
        cu.remove_th1_errors(data_no_errors)
        for i, cont in enumerate(conts[1:], 1):
            if cont.subplot == data_obj:
                cont.subplot = data_no_errors
                # compare_bins(cont.obj, cont.subplot)

    plot = Plot(conts,
                what='hist',
                legend=do_legend,
                has_data=data_first,
                xlim=xlim,
                ylim=ylim,
                xtitle="p_{T}^{jet} [GeV]",
                # ytitle="#DeltaN/N" if normalise_hists else "N",
                ytitle="#frac{dN}{dp_{T}} [events / GeV]",
                title=title,
                subplot_title=subplot_title,
                subplot_type='ratio',
                subplot_limits=subplot_limits,
                lumi=lumi,
                is_preliminary=is_preliminary)
    plot.y_padding_max_log = 1E4
    plot.legend.SetX1(0.6)
    plot.legend.SetY1(0.68)
    plot.legend.SetX2(0.98)
    plot.legend.SetY2(0.88)
    plot.left_margin = 0.16
    plot.title_left_offset = 0.05

    if data_first:
        # we'll do the filling of legend ourselves
        plot.do_legend = False

    draw_opt = "NOSTACK L E"
    plot.plot(draw_opt)

    # avoid x title hitting labels
    plot.subplot_container.GetXaxis().SetLabelSize(plot.subplot_container.GetXaxis().GetLabelSize()*0.9)

    plot.set_logx(do_more_labels=True, do_exponent=False)
    plot.set_logy(do_more_labels=False)

    graphs=[]

    # Special case if data first object
    if data_first:
        # Redraw data on top
        data_hist = plot.contributions[0].obj.Clone()
        data_hist.SetLineWidth(entries[0][1].get('line_width', 1))
        data_hist.SetMarkerSize(entries[0][1].get('marker_size', 1))
        plot.main_pad.cd()
        data_draw_opt = "E1 X0 SAME"
        data_hist.Draw(data_draw_opt)

        for e in [p.obj for p in plot.contributions[1:]]+[data_hist]:
           graphs+=[ROOT.TGraphErrors(e)]
           graphs[-1].Draw("P SAME")

        # Create dummy graphs with the same styling to put into the legend
        # Using graphs we can get the correct endings on the TLegend entries (!)
        # Yes a massive faff for something so trivial
        dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        dummy_hist = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        # Add them to the legend and draw it
        dummy_conts = []  # stop premature deletion
        for i, entry in enumerate(entries):
            leg_draw_opt = "LE"
            # check if line_width > 0?
            if i == 0:
                if "X0" in data_draw_opt:
                    leg_draw_opt = "LE"
                if data_hist.GetMarkerSize() > 0:
                    leg_draw_opt += "P"
            # make Contribution just to ease styling methods
            cont = Contribution(dummy_gr.Clone(), leg_draw_opt=leg_draw_opt, **entry[1])
            dummy_conts.append(cont)
            plot.legend.AddEntry(cont.obj, cont.label, cont.leg_draw_opt)

        plot.canvas.cd()
        plot.legend.Draw()

        # Do the subplot uncertainty shading
        data_no_errors = entries[0][0].Clone()
        cu.remove_th1_errors(data_no_errors)

        # Create hists for data with error region for ratio
        # Easiest way to get errors right is to do data (with 0 errors)
        # and divide by data (with errors), as if you had MC = data with 0 error

        # data_stat_ratio = data_no_errors.Clone()
        # data_stat_ratio.Divide(unfolded_hist_bin_stat_errors)
        # data_stat_ratio.SetFillStyle(3354)
        # data_stat_ratio.SetFillColor(self.plot_colours['unfolded_stat_colour'])
        # data_stat_ratio.SetLineWidth(0)
        # data_stat_ratio.SetMarkerSize(0)

        data_total_ratio = data_no_errors.Clone()
        # compare_bins(data_total_ratio, entries[0][0])
        data_total_ratio.Divide(entries[0][0])
        data_total_ratio.SetFillStyle(3545)
        # data_total_ratio.SetFillStyle(3002)
        data_total_ratio.SetFillColor(entries[0][1]['fill_color'])
        data_total_ratio.SetLineWidth(0)
        data_total_ratio.SetMarkerSize(0)

        # now draw the data error shaded area
        # this is a bit hacky - basically draw them on the ratio pad,
        # then redraw the existing hists & line to get them ontop
        # note that we use "same" for all - this is to keep the original axes
        # (we may want to rethink this later?)
        plot.subplot_pad.cd()
        data_draw_opt = "E2 SAME ]["  # need SAME otherwise axis get rescaled
        # data_stat_ratio.Draw(data_draw_opt)
        data_total_ratio.Draw(data_draw_opt)

        # draw small legend for shadings
        x_left = 0.25
        y_bottom = 0.75
        width = 0.63
        height = 0.15
        plot.subplot_leg = ROOT.TLegend(x_left, y_bottom, x_left+width, y_bottom+height)
        # plot.subplot_leg = ROOT.TLegend(0.25, 0.73, 0.9, 0.9)
        # plot.subplot_leg.SetTextSize(0.07)
        plot.subplot_leg.SetFillStyle(0)
        plot.subplot_leg.SetNColumns(2)
        plot.subplot_leg.AddEntry(data_total_ratio, qgc.DATA_STAT_UNC_STR, "F")

        # Do systematic shading
        if experimental_syst:
            experimental_syst.SetFillStyle(3005)
            # experimental_syst.SetFillStyle(3002)
            experimental_syst.SetFillColor(ROOT.kAzure)
            experimental_syst.SetLineWidth(0)
            experimental_syst.SetMarkerSize(0)
            experimental_syst.Draw("SAME 2")

        if scale_syst:
            scale_syst.SetFillStyle(3003)
            # scale_syst.SetFillStyle(3002)
            scale_syst.SetFillColor(ROOT.kOrange)
            scale_syst.SetLineWidth(0)
            scale_syst.SetMarkerSize(0)
            scale_syst.Draw("SAME 2")

        if pdf_syst:
            pdf_syst.SetFillStyle(3002)
            pdf_syst.SetFillColor(ROOT.kMagenta)
            pdf_syst.SetLineWidth(0)
            pdf_syst.SetMarkerSize(0)
            pdf_syst.Draw("SAME 2")

        if total_syst:
            total_syst.SetFillStyle(3003)
            total_syst.SetFillStyle(3354)
            # total_syst.SetFillColor(ROOT.kRed)
            # total_syst.SetFillColor(qgc.QCD_COLOUR)
            total_syst.SetFillColor(entries[1][1]['fill_color'])
            total_syst.SetLineWidth(0)
            total_syst.SetMarkerSize(0)
            total_syst.Draw("SAME 2")

        plot.subplot_container.Draw("SAME " + draw_opt)
        plot.subplot_line.Draw()

        for e in plot.subplot_contributions:
           graphs+=[ROOT.TGraphErrors(e)]
           graphs[-1].Draw("E SAME")
        for e in [data_total_ratio]:
           graphs+=[ROOT.TGraphErrors(e)]
           graphs[-1].Draw("E2 SAME ][")

        if total_syst and not any([experimental_syst, scale_syst, pdf_syst]):
            plot.subplot_leg.SetTextSize(0.085)

        if experimental_syst:
            plot.subplot_leg.AddEntry(experimental_syst, "Exp. syst.", "F")
        if scale_syst:
            plot.subplot_leg.AddEntry(scale_syst, "Scale syst.", "F")
        if pdf_syst:
            plot.subplot_leg.AddEntry(pdf_syst, "PDF syst.", "F")
        if total_syst:
            plot.subplot_leg.AddEntry(total_syst, "MC total unc.", "F")

        plot.subplot_leg.Draw()

        plot.canvas.cd()

    print(output_filename)
    plot.save(output_filename)


def _rebin_scale(hist, binning, normalise=False):
    if binning is not None:
        new_hist = hist.Rebin(len(binning)-1, hist.GetName()+"_rebin" + cu.get_unique_str(), binning)
    else:
        new_hist = hist.Clone(hist.GetName()+"_rebin_scale" + cu.get_unique_str())
    factor = 1.
    if normalise and new_hist.Integral() != 0:
        factor = new_hist.Integral()
    # Urgh this breaks .Integral() somehow
    new_hist.Scale(1./factor, "width")
    return new_hist


def tunfold_to_physical_bins(hist, bins, divide_by_bin_width=True):
    if hist.GetNbinsX() != len(bins)-1:
        raise ValueError("tunfold_to_physical_bins: bins not correct size (%d vs %d)" % (hist.GetNbinsX(), len(bins)-1))

    new_hist = ROOT.TH1D(hist.GetName()+"_relabel" + cu.get_unique_str(), "", len(bins)-1, bins)
    for ix in range(1, hist.GetNbinsX()+1):
        bin_width = bins[ix] - bins[ix-1]
        if not divide_by_bin_width:
            bin_width = 1.
        # print(ix, hist.GetBinError(ix), hist.GetBinError(ix) / bin_width)
        new_hist.SetBinContent(ix, hist.GetBinContent(ix) / bin_width)
        new_hist.SetBinError(ix, hist.GetBinError(ix) / bin_width)
    new_hist.SetEntries(hist.GetEntries())
    return new_hist


angle = [a for a in qgc.COMMON_VARS if a.var in 'jet_LHA'][0]
LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT['ungroomed']
variable_bin_edges_reco = LAMBDA_VAR_DICTS[angle.var]['reco']
variable_bin_edges_gen = LAMBDA_VAR_DICTS[angle.var]['gen']
variable_name = angle.name

pt_bin_edges_gen = qgc.PT_UNFOLD_DICT['signal_gen']
pt_bin_edges_reco = qgc.PT_UNFOLD_DICT['signal_reco']
pt_bin_edges_underflow_gen = qgc.PT_UNFOLD_DICT['underflow_gen']
pt_bin_edges_underflow_reco = qgc.PT_UNFOLD_DICT['underflow_reco']

zpj_append = "_zpj"
pt_bin_edges_zpj_gen = qgc.PT_UNFOLD_DICT['signal%s_gen' % (zpj_append)]
pt_bin_edges_zpj_reco = qgc.PT_UNFOLD_DICT['signal%s_reco' % (zpj_append)]
pt_bin_edges_zpj_underflow_gen = qgc.PT_UNFOLD_DICT['underflow%s_gen' % (zpj_append)]
pt_bin_edges_zpj_underflow_reco = qgc.PT_UNFOLD_DICT['underflow%s_reco' % (zpj_append)]


def create_pt_hist(hist, binning, binning_uflow, pt_bin_edges, pt_bin_edges_uflow, variable_bin_edges):
    """Create 1D pt hist from var & pt 1d hist"""
    # THIS IS A HORRIBLE HACK BECAUSE I DIDNT FILL MY HISTS
    all_pt_bins = list(np.append(pt_bin_edges_uflow[:-1], pt_bin_edges))
    all_pt_bins.append(8000)
    # print(all_pt_bins)
    nbins_pt = len(all_pt_bins)-1
    h_new = ROOT.TH1D("hpt"+cu.get_unique_str(), "", nbins_pt, array('d', all_pt_bins))
    for pt_ind in range(1, h_new.GetNbinsX()+1):
        this_sum = 0
        this_err_sq = 0
        this_pt = all_pt_bins[pt_ind-1]
        # ARGH THIS IS SO FRUSTRATING
        this_binning = binning if this_pt >= pt_bin_edges[0] else binning_uflow
        for var_ind, var in enumerate(variable_bin_edges[:-1]):
            bin_num = this_binning.GetGlobalBinNumber(var, this_pt)
            this_sum += hist.GetBinContent(bin_num)
            this_err_sq += hist.GetBinError(bin_num)**2
        h_new.SetBinContent(pt_ind, this_sum)
        h_new.SetBinError(pt_ind, math.sqrt(this_err_sq))
    return h_new


def do_dijet_pt_plots(workdir,
                      subplot_vs_data=True, # otherwise vs MC
                      show_total_systematics=True,
                      show_grouped_systematics=True,
                      show_individual_systematics=True,
                      is_preliminary=True):

    data_tfile = cu.open_root_file(os.path.join(workdir, qgc.JETHT_ZB_FILENAME))
    mg_tfile = cu.open_root_file(os.path.join(workdir, qgc.QCD_FILENAME))
    # py_tfile = cu.open_root_file(os.path.join(workdir, qgc.QCD_PYTHIA_ONLY_FILENAME))
    hpp_tfile = cu.open_root_file(os.path.join(workdir, qgc.QCD_HERWIG_FILENAME))

    source_dir_systs = os.path.join(workdir, "systematics_files")

    experimental_systematics = [
        {
            "label": "Charged hadron up",
            "tfile": os.path.join(source_dir_systs, 'chargedHadronShiftUp', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure+1,
        },
        {
            "label": "Charged hadron down",
            "tfile": os.path.join(source_dir_systs, 'chargedHadronShiftDown', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure+1,
            "linestyle": 2,
        },
        {
            "label": "Neutral hadron up",
            "tfile": os.path.join(source_dir_systs, 'neutralHadronShiftUp', qgc.QCD_FILENAME),
            "colour": ROOT.kOrange-4,
        },
        {
            "label": "Neutral hadron down",
            "tfile": os.path.join(source_dir_systs, 'neutralHadronShiftDown', qgc.QCD_FILENAME),
            "colour": ROOT.kOrange-4,
            "linestyle": 2,
        },
        {
            "label": "Photon up",
            "tfile": os.path.join(source_dir_systs, 'photonShiftUp', qgc.QCD_FILENAME),
            "colour": ROOT.kMagenta-3,
        },
        {
            "label": "Photon down",
            "tfile": os.path.join(source_dir_systs, 'photonShiftDown', qgc.QCD_FILENAME),
            "colour": ROOT.kMagenta-3,
            "linestyle": 2,
        },
        {
            "label": "JES up",
            "tfile": os.path.join(source_dir_systs, 'jecsmear_directionUp', qgc.QCD_FILENAME),
            "colour": ROOT.kGreen+3,
        },
        {
            "label": "JES down",
            "tfile": os.path.join(source_dir_systs, 'jecsmear_directionDown', qgc.QCD_FILENAME),
            "colour": ROOT.kGreen+3,
            "linestyle": 2,
        },
        {
            "label": "JER up",
            "tfile": os.path.join(source_dir_systs, 'jersmear_directionUp', qgc.QCD_FILENAME),
            "colour": ROOT.kOrange+3,
        },
        {
            "label": "JER down",
            "tfile": os.path.join(source_dir_systs, 'jersmear_directionDown', qgc.QCD_FILENAME),
            "colour": ROOT.kOrange+3,
            "linestyle": 2,
        },
        {
            "label": "Pileup up",
            "tfile": os.path.join(source_dir_systs, 'pileup_directionUp', qgc.QCD_FILENAME),
            "colour": ROOT.kBlue-4,
        },
        {
            "label": "Pileup down",
            "tfile": os.path.join(source_dir_systs, 'pileup_directionDown', qgc.QCD_FILENAME),
            "colour": ROOT.kBlue-4,
            "linestyle": 2,
        },
        # {
        #     "label": "Tracking up",
        #     "tfile": os.path.join(source_dir_systs, 'track_directionUp', qgc.QCD_FILENAME),
        #     "colour": ROOT.kMagenta+3,
        # },
        # {
        #     "label": "Tracking down",
        #     "tfile": os.path.join(source_dir_systs, 'track_directionDown', qgc.QCD_FILENAME),
        #     "colour": ROOT.kMagenta+3,
        #     "linestyle": 2,
        # },
    ]

    scale_systematics = [
        {
            "label": "muR up, muF nominal",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure,
        },
        {
            "label": "muR down, muF nominal",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure+1,
        },
        {
            "label": "muR nominal, muF up",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure+2,
        },
        {
            "label": "muR nominal, muF down",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure+3,
        },
        {
            "label": "muR down, muF down",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure+4,
        },
        {
            "label": "muR up, muF up",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure+5,
        },
    ]

    pdf_systematics = [
        {
            "label": "PDF",  # this is a template entry, used for future
            # this file has the newer PDF hists for reco pt
            "tfile": os.path.join(source_dir_systs, 'PDFvariationsTrue', qgc.QCD_FILENAME),
            # "tfile": os.path.join('/Users/robin/Projects/QGAnalysis/workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts_fixPassGen_jackknife/systematic_files', 'PDFvariationsTrue', qgc.QCD_FILENAME),
            "colour": ROOT.kCyan+2,
            "variations": range(100),  # list of all the variation #s to be used
        },
    ]

    for region_shortname, region_label in [("central", qgc.Dijet_CEN_LABEL), ("forward", qgc.Dijet_FWD_LABEL)]:

        # histname = "Dijet_QG_%s_tighter/jet_pt" % (region_shortname)  # fine equidistant binning
        histname = "Dijet_QG_Unfold_%s_tighter/hist_pt_reco_all" % (region_shortname) # tunfold binning

        lw = 2
        msize = 1.1
        data_line_width = lw
        mc_msize = 1E-3  # small enough to not be seen but not 0 - otherwise the horizontal lines on error bars don't get drawn

        mgpy_label = "MG5+Pythia8"
        hpp_label = "Herwig++"

        # print(mg_tfile)
        data_hist = cu.get_from_tfile(data_tfile, histname)
        mg_hist = cu.get_from_tfile(mg_tfile, histname)
        # py_hist = cu.get_from_tfile(py_tfile, histname)
        hpp_hist = cu.get_from_tfile(hpp_tfile, histname)

        # For use with tunfold binning, which just has bin indices as x values instead of physical values
        all_pt_bins = np.append(qgc.PT_UNFOLD_DICT['underflow_reco'][:-1], qgc.PT_UNFOLD_DICT['signal_reco'])
        all_pt_bins = np.append(all_pt_bins, 8000)  # the overflow bin
        # print(all_pt_bins)
        data_hist = tunfold_to_physical_bins(data_hist, all_pt_bins, divide_by_bin_width=True)
        mg_hist = tunfold_to_physical_bins(mg_hist, all_pt_bins, divide_by_bin_width=True)
        # py_hist = tunfold_to_physical_bins(py_hist, all_pt_bins, divide_by_bin_width=True)
        hpp_hist = tunfold_to_physical_bins(hpp_hist, all_pt_bins, divide_by_bin_width=True)

        # Scale to data
        mg_sf = data_hist.Integral()/mg_hist.Integral()
        mg_hist.Scale(mg_sf)
        # py_hist.Scale(data_hist.Integral()/py_hist.Integral())
        hpp_hist.Scale(data_hist.Integral()/hpp_hist.Integral())

        # Absolute shifted variations
        exp_hist_up, exp_hist_down = None, None
        scale_hist_up, scale_hist_down = None, None
        pdf_hist_up, pdf_hist_down = None, None
        total_hist_up, total_hist_down = None, None

        # Calculate various systematic variations
        if any([show_total_systematics, show_individual_systematics, show_grouped_systematics]):
            # Calculate total experimental systematic uncertainty
            # ------------------------------------------------------------------
            for exp_dict in experimental_systematics:
                if isinstance(exp_dict['tfile'], str):
                    exp_dict['tfile'] = cu.open_root_file(exp_dict['tfile'])
                hist = cu.get_from_tfile(exp_dict['tfile'], histname)
                hist = tunfold_to_physical_bins(hist, all_pt_bins, divide_by_bin_width=True)
                hist.Scale(mg_sf)
                exp_dict['hist'] = hist

            # create envelope of quadrature sum max/min per bin
            exp_hist_up = mg_hist.Clone("exp_hist_up")
            exp_hist_down = mg_hist.Clone("exp_hist_down")
            for ix in range(1, mg_hist.GetNbinsX()+1):
                nominal_bin = mg_hist.GetBinContent(ix)
                all_bin_diffs = [exp_dict['hist'].GetBinContent(ix) - nominal_bin
                                 for exp_dict in experimental_systematics]

                positive_err = math.sqrt(sum([x**2 for x in all_bin_diffs if x > 0]))
                exp_hist_up.SetBinContent(ix, nominal_bin + positive_err)
                exp_hist_up.SetBinError(ix, 0)

                negative_err = math.sqrt(sum([x**2 for x in all_bin_diffs if x < 0]))
                exp_hist_down.SetBinContent(ix, nominal_bin - negative_err)
                exp_hist_down.SetBinError(ix, 0)

                # print(ix, 'nominal', nominal_bin)
                # print(ix, 'positive_err', positive_err, positive_err/nominal_bin)
                # print(ix, 'negative_err', negative_err, negative_err/nominal_bin)


            # Calculate total scale uncertainty
            # ------------------------------------------------------------------
            for scale_dict in scale_systematics:
                if isinstance(scale_dict['tfile'], str):
                    scale_dict['tfile'] = cu.open_root_file(scale_dict['tfile'])
                hist = cu.get_from_tfile(scale_dict['tfile'], histname)
                hist = tunfold_to_physical_bins(hist, all_pt_bins, divide_by_bin_width=True)
                hist.Scale(mg_sf)
                scale_dict['hist'] = hist

            # create envelope of max/min per bin
            scale_hist_up = mg_hist.Clone("scale_hist_up")
            scale_hist_down = mg_hist.Clone("scale_hist_down")
            for ix in range(1, mg_hist.GetNbinsX()+1):
                nominal_bin = mg_hist.GetBinContent(ix)
                all_bin_values = [scale_dict['hist'].GetBinContent(ix)
                                  for scale_dict in scale_systematics]
                scale_hist_up.SetBinContent(ix, max(all_bin_values))
                scale_hist_up.SetBinError(ix, 0)
                scale_hist_down.SetBinContent(ix, min(all_bin_values))
                scale_hist_down.SetBinError(ix, 0)


            # Calculate total PDF uncertainty
            # ------------------------------------------------------------------
            # Create all dicts first
            tfile = pdf_systematics[0]['tfile']
            if isinstance(tfile, str):
                tfile = cu.open_root_file(tfile)

            these_pdf_systematics = []
            num_vars = len(pdf_systematics[0]['variations'])
            for pdf_ind in pdf_systematics[0]['variations']:
                hist = cu.get_from_tfile(tfile, "Dijet_QG_Unfold_%s_tighter/hist_pt_reco_all_PDF_%d" % (region_shortname, pdf_ind))
                hist = tunfold_to_physical_bins(hist, all_pt_bins, divide_by_bin_width=True)
                hist.Scale(mg_sf)
                these_pdf_systematics.append(
                    {
                        "label": "PDF_%d" % (pdf_ind),
                        "hist": hist,
                        "colour": cu.get_colour_seq(pdf_ind, num_vars)
                    })

            # create RMS up/down
            pdf_hist_up = mg_hist.Clone("pdf_hist_up")
            pdf_hist_down = mg_hist.Clone("pdf_hist_down")
            for ix in range(1, mg_hist.GetNbinsX()+1):
                nominal = mg_hist.GetBinContent(ix)
                values = [pdf_dict['hist'].GetBinContent(ix) for pdf_dict in these_pdf_systematics]
                rms = np.std(values, ddof=1)
                pdf_hist_up.SetBinContent(ix, nominal+rms)
                pdf_hist_up.SetBinError(ix, 0)
                pdf_hist_down.SetBinContent(ix, max(nominal-rms, 0))
                pdf_hist_down.SetBinError(ix, 0)
                # print(ix, 'nominal', nominal, 'rms', rms, 'rms/nominal', rms/nominal)
                # if rms > nominal:
                #     print("!!!!!")

            # Calculate total uncertainty
            # ------------------------------------------------------------------
            total_hist_up = mg_hist.Clone('total_hist_up')
            total_hist_down = mg_hist.Clone('total_hist_down')
            for ix in range(1, mg_hist.GetNbinsX()+1):
                nominal = mg_hist.GetBinContent(ix)
                sum_sq = sum([(up_hist.GetBinContent(ix) - nominal)**2
                              for up_hist in [exp_hist_up, scale_hist_up, pdf_hist_up]
                              if up_hist is not None])
                # print(ix, [(up_hist.GetBinContent(ix) - nominal)**2
                #           for up_hist in [exp_hist_up, scale_hist_up, pdf_hist_up]
                #           if up_hist is not None])
                total_hist_up.SetBinContent(ix, nominal+math.sqrt(sum_sq))
                total_hist_up.SetBinError(ix, 0)

                sum_sq = sum([(down_hist.GetBinContent(ix) - nominal)**2
                              for down_hist in [exp_hist_down, scale_hist_down, pdf_hist_down]
                              if down_hist is not None])
                total_hist_down.SetBinContent(ix, nominal - math.sqrt(sum_sq))
                total_hist_down.SetBinError(ix, 0)

        # Create entries for plot
        # --------------------------------------------------------------------------
        ref_hist = data_hist if subplot_vs_data else mg_hist
        entries = [
            # DATA
            [
                data_hist,
                dict(line_color=ROOT.kBlack, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                     marker_color=ROOT.kBlack, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=msize*0.7,
                     label="Data",
                     subplot=None if subplot_vs_data else ref_hist)
            ],

            # MG5+PYTHIA8 MC
            [
                mg_hist,
                dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                     marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                     label=mgpy_label,
                     subplot=ref_hist if subplot_vs_data else None)
            ],

            # PYTHIA-only MC
            # [
            #     py_hist,
            #     dict(line_color=qgc.QCD_COLOURS[2], line_width=lw, fill_color=qgc.QCD_COLOURS[2],
            #          marker_color=qgc.QCD_COLOURS[2], marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=mc_msize,
            #          label="Pythia8",
            #          subplot=ref_hist)
            # ],

            # HERWIG++
            [
                hpp_hist,
                dict(line_color=qgc.HERWIGPP_QCD_COLOUR, line_width=lw, line_style=2, fill_color=qgc.HERWIGPP_QCD_COLOUR,
                     marker_color=qgc.HERWIGPP_QCD_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                     label=hpp_label,
                     subplot=ref_hist)
            ]
        ]

        # ADD EXPERIMENTAL ENTRIES
        # ------------------------------------------------------------------
        if show_individual_systematics:
            # plot individual variations
            for exp_dict in experimental_systematics:
                entries.append([exp_dict['hist'],
                                dict(line_color=exp_dict['colour'], line_width=lw, fill_color=exp_dict['colour'],
                                     marker_color=exp_dict['colour'], marker_size=0, label=exp_dict['label'],
                                     subplot=data_hist)
                               ])

        if show_grouped_systematics:
            # plot up/down boundaries
            exp_col = ROOT.kRed+2
            exp_col2 = ROOT.kRed-2
            entries.extend([
                [
                    exp_hist_up.Clone(),
                    dict(line_color=exp_col, line_width=lw, line_style=1,
                         fill_color=exp_col,
                         marker_color=exp_col, marker_size=0, label="Exp systs up",
                         subplot=data_hist)
                ],
                [
                    exp_hist_down.Clone(),
                    dict(line_color=exp_col2, line_width=lw, line_style=2,
                         fill_color=exp_col2,
                         marker_color=exp_col2, marker_size=0, label="Exp systs down",
                         subplot=data_hist)
                ],
            ])

        n = mg_hist.GetNbinsX()
        bin_width = [mg_hist.GetBinWidth(i)/2. for i in range(1, n+1)]
        x = array('d', [mg_hist.GetBinCenter(i) for i in range(1, n+1)])
        y = array('d', [1 for i in range(1, n+1)])
        exl = array('d', [bin_width[i-1] for i in range(1, n+1)])
        exh = array('d', [bin_width[i-1] for i in range(1, n+1)])

        # exp_hist_up_ratio = exp_hist_up.Clone()
        # exp_hist_down_ratio = exp_hist_down.Clone()
        # exp_hist_up_ratio.Add(mg_hist, -1) # error doesnt matter
        # exp_hist_up_ratio.Divide(mg_hist)
        # exp_hist_down_ratio.Add(mg_hist, -1)
        # exp_hist_down_ratio.Divide(mg_hist)
        # exp_hist_down_ratio.Scale(-1)
        # exp_eyl = array('d', [exp_hist_down_ratio.GetBinContent(i) for i in range(1, n+1)])
        # exp_eyh = array('d', [exp_hist_up_ratio.GetBinContent(i) for i in range(1, n+1)])
        # exp_gr = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, exp_eyl, exp_eyh)

        # ADD SCALE ENTRIES
        # ------------------------------------------------------------------
        if show_individual_systematics:
            for scale_dict in scale_systematics:
                entries.append([scale_dict['hist'],
                                dict(line_color=scale_dict['colour'], line_width=lw, fill_color=scale_dict['colour'],
                                     marker_color=scale_dict['colour'], marker_size=0, label=scale_dict['label'],
                                     subplot=ref_hist)
                               ])

        if show_grouped_systematics:
            scale_col = ROOT.kGreen+2
            scale_col2 = scale_col
            entries.extend([
                [
                    scale_hist_up.Clone(),
                    dict(line_color=scale_col, line_width=lw, line_style=1,
                         fill_color=scale_col,
                         marker_color=scale_col, marker_size=0, label="Scale up",
                         subplot=ref_hist)
                ],
                [
                    scale_hist_down.Clone(),
                    dict(line_color=scale_col2, line_width=lw, line_style=2,
                         fill_color=scale_col2,
                         marker_color=scale_col2, marker_size=0, label="Scale down",
                         subplot=ref_hist)
                ],
            ])

        # scale_hist_up_ratio = scale_hist_up.Clone()
        # scale_hist_down_ratio = scale_hist_down.Clone()
        # scale_hist_up_ratio.Add(mg_hist, -1)
        # scale_hist_up_ratio.Divide(mg_hist)
        # scale_hist_down_ratio.Add(mg_hist, -1)
        # scale_hist_down_ratio.Divide(mg_hist)
        # scale_hist_down_ratio.Scale(-1)
        # scale_eyl = array('d', [scale_hist_down_ratio.GetBinContent(i) for i in range(1, n+1)])
        # scale_eyh = array('d', [scale_hist_up_ratio.GetBinContent(i) for i in range(1, n+1)])
        # scale_gr = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, scale_eyl, scale_eyh)

        # ADD PDF ENTRIES
        # ------------------------------------------------------------------
        if show_individual_systematics:
            for pdf_dict in these_pdf_systematics:
                entries.append(
                    [
                        pdf_dict['hist'],
                        dict(label=pdf_dict['label'],
                             line_color=pdf_dict['colour'],
                             # subplot=mg_hist)
                             subplot=ref_hist)
                    ]
                )

        if show_grouped_systematics:
            pdf_col = ROOT.kOrange-4
            pdf_col2 = pdf_col
            entries.extend([
                [
                    pdf_hist_up.Clone(),
                    dict(line_color=pdf_col, line_width=lw, line_style=1,
                         fill_color=pdf_col,
                         marker_color=pdf_col, marker_size=0, label="PDF up",
                         subplot=ref_hist)
                ],
                [
                    pdf_hist_down.Clone(),
                    dict(line_color=pdf_col2, line_width=lw, line_style=2,
                         fill_color=pdf_col2,
                         marker_color=pdf_col2, marker_size=0, label="PDF down",
                         subplot=ref_hist)
                ],
            ])

        # pdf_hist_up_ratio = pdf_hist_up.Clone()
        # pdf_hist_down_ratio = pdf_hist_down.Clone()
        # pdf_hist_up_ratio.Add(mg_hist, -1)
        # pdf_hist_up_ratio.Divide(mg_hist)
        # pdf_hist_down_ratio.Add(mg_hist, -1)
        # pdf_hist_down_ratio.Divide(mg_hist)
        # pdf_hist_down_ratio.Scale(-1)
        # pdf_eyl = array('d', [pdf_hist_down_ratio.GetBinContent(i) for i in range(1, n+1)])
        # pdf_eyh = array('d', [pdf_hist_up_ratio.GetBinContent(i) for i in range(1, n+1)])
        # pdf_gr = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, pdf_eyl, pdf_eyh)

        # ADD TOTAL
        # ------------------------------------------------------------------
        if show_grouped_systematics:
            total_col = ROOT.kBlack
            total_col2 = total_col
            entries.extend([
                [
                    total_hist_up.Clone(),
                    dict(line_color=total_col, line_width=lw, line_style=1,
                         fill_color=total_col,
                         marker_color=total_col, marker_size=0, label="Total up",
                         subplot=data_hist)
                ],
                [
                    total_hist_down.Clone(),
                    dict(line_color=total_col2, line_width=lw, line_style=2,
                         fill_color=total_col2,
                         marker_color=total_col2, marker_size=0, label="Total down",
                         subplot=data_hist)
                ],
            ])

        if show_total_systematics:
            total_hist_up_ratio = total_hist_up.Clone()
            total_hist_down_ratio = total_hist_down.Clone()
            total_hist_up_ratio.Add(mg_hist, -1)
            total_hist_up_ratio.Divide(mg_hist)
            total_hist_down_ratio.Add(mg_hist, -1)
            total_hist_down_ratio.Divide(mg_hist)
            total_hist_down_ratio.Scale(-1)
            total_eyl = array('d', [total_hist_down_ratio.GetBinContent(i) for i in range(1, n+1)])
            total_eyh = array('d', [total_hist_up_ratio.GetBinContent(i) for i in range(1, n+1)])
            total_gr = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, total_eyl, total_eyh)

        radius, pus = cu.get_jet_config_from_dirname(workdir)
        jet_str = "AK%s" % (radius.upper())
        title = "{jet_algo}\n{region_label}".format(jet_algo=jet_str,
                                                      region_label=region_label)
        subplot_title = "* / %s" % entries[1][1]['label']
        if subplot_vs_data:
            subplot_title = qgc.SIM_DATA_STR
        paper_str = "_paper" if not is_preliminary else ""
        do_jet_pt_plot(entries,
                       output_filename=os.path.join(workdir, "data_mc_jet_pt/Dijet_%s/jet_pt%s.%s" % (region_shortname, paper_str, OUTPUT_FMT)),
                       rebin=1,
                       xlim=(30, qgc.PT_UNFOLD_DICT['signal_gen'][-1]),
                       ylim=(5E-3, 1E14),
                       title=title,
                       subplot_limits=(0, 2.5),
                       subplot_title=subplot_title,
                       data_first=True,
                       normalise_hists=False,
                       # experimental_syst=exp_gr,
                       # scale_syst=scale_gr,
                       # pdf_syst=pdf_gr,
                       experimental_syst=None,
                       scale_syst=None,
                       pdf_syst=None,
                       total_syst=total_gr if show_total_systematics else None,
                       lumi=cu.get_lumi_str(do_dijet=True, do_zpj=False),
                       is_preliminary=is_preliminary)


def do_zpj_pt_plots(workdir,
                    subplot_vs_data=True, # otherwise vs MC
                    show_total_systematics=True,
                    show_grouped_systematics=True,
                    show_individual_systematics=True,
                    is_preliminary=True):

    single_mu_tfile = cu.open_root_file(os.path.join(workdir, qgc.SINGLE_MU_FILENAME))
    mg_dy_tfile = cu.open_root_file(os.path.join(workdir, qgc.DY_FILENAME))
    hpp_dy_tfile = cu.open_root_file(os.path.join(workdir, qgc.DY_HERWIG_LOW_HIGH_PT_FILENAME))

    # histname = "ZPlusJets_QG/jet_pt"  # fine equidistant binning
    histname = "ZPlusJets_QG_Unfold/hist_pt_reco_all"  # tunfold binning

    source_dir_systs = os.path.join(workdir, "systematics_files")

    experimental_systematics = [
        {
            "label": "Charged hadron up",
            "tfile": os.path.join(source_dir_systs, 'chargedHadronShiftUp', qgc.DY_FILENAME),
            "colour": ROOT.kAzure+1,
        },
        {
            "label": "Charged hadron down",
            "tfile": os.path.join(source_dir_systs, 'chargedHadronShiftDown', qgc.DY_FILENAME),
            "colour": ROOT.kAzure+1,
            "linestyle": 2,
        },
        {
            "label": "Neutral hadron up",
            "tfile": os.path.join(source_dir_systs, 'neutralHadronShiftUp', qgc.DY_FILENAME),
            "colour": ROOT.kOrange-4,
        },
        {
            "label": "Neutral hadron down",
            "tfile": os.path.join(source_dir_systs, 'neutralHadronShiftDown', qgc.DY_FILENAME),
            "colour": ROOT.kOrange-4,
            "linestyle": 2,
        },
        {
            "label": "Photon up",
            "tfile": os.path.join(source_dir_systs, 'photonShiftUp', qgc.DY_FILENAME),
            "colour": ROOT.kMagenta-3,
        },
        {
            "label": "Photon down",
            "tfile": os.path.join(source_dir_systs, 'photonShiftDown', qgc.DY_FILENAME),
            "colour": ROOT.kMagenta-3,
            "linestyle": 2,
        },
        {
            "label": "JES up",
            "tfile": os.path.join(source_dir_systs, 'jecsmear_directionUp', qgc.DY_FILENAME),
            "colour": ROOT.kGreen+3,
        },
        {
            "label": "JES down",
            "tfile": os.path.join(source_dir_systs, 'jecsmear_directionDown', qgc.DY_FILENAME),
            "colour": ROOT.kGreen+3,
            "linestyle": 2,
        },
        {
            "label": "JER up",
            "tfile": os.path.join(source_dir_systs, 'jersmear_directionUp', qgc.DY_FILENAME),
            "colour": ROOT.kOrange+3,
        },
        {
            "label": "JER down",
            "tfile": os.path.join(source_dir_systs, 'jersmear_directionDown', qgc.DY_FILENAME),
            "colour": ROOT.kOrange+3,
            "linestyle": 2,
        },
        {
            "label": "Pileup up",
            "tfile": os.path.join(source_dir_systs, 'pileup_directionUp', qgc.DY_FILENAME),
            "colour": ROOT.kBlue-4,
        },
        {
            "label": "Pileup down",
            "tfile": os.path.join(source_dir_systs, 'pileup_directionDown', qgc.DY_FILENAME),
            "colour": ROOT.kBlue-4,
            "linestyle": 2,
        },
        # {
        #     "label": "Tracking up",
        #     "tfile": os.path.join(source_dir_systs, 'track_directionUp', qgc.DY_FILENAME),
        #     "colour": ROOT.kMagenta+3,
        # },
        # {
        #     "label": "Tracking down",
        #     "tfile": os.path.join(source_dir_systs, 'track_directionDown', qgc.DY_FILENAME),
        #     "colour": ROOT.kMagenta+3,
        #     "linestyle": 2,
        # },
    ]

    scale_systematics = [
        {
            "label": "muR up, muF nominal",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.DY_FILENAME),
            "colour": ROOT.kAzure,
        },
        {
            "label": "muR down, muF nominal",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.DY_FILENAME),
            "colour": ROOT.kAzure+1,
        },
        {
            "label": "muR nominal, muF up",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.DY_FILENAME),
            "colour": ROOT.kAzure+2,
        },
        {
            "label": "muR nominal, muF down",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.DY_FILENAME),
            "colour": ROOT.kAzure+3,
        },
        {
            "label": "muR down, muF down",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.DY_FILENAME),
            "colour": ROOT.kAzure+4,
        },
        {
            "label": "muR up, muF up",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.DY_FILENAME),
            "colour": ROOT.kAzure+5,
        },
    ]

    pdf_systematics = [
        {
            "label": "PDF",  # this is a template entry, used for future
            # this file has the newer PDF hists for reco pt
            "tfile": os.path.join(source_dir_systs, 'PDFvariationsTrue', qgc.DY_FILENAME),
            # "tfile": os.path.join('/Users/robin/Projects/QGAnalysis/workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts_fixPassGen_jackknife/systematic_files', 'PDFvariationsTrue', qgc.QCD_FILENAME),
            "colour": ROOT.kCyan+2,
            "variations": range(100),  # list of all the variation #s to be used
        },
    ]

    lw = 2
    msize = 1.1
    data_line_width = lw
    mc_msize = 1E-3  # small enough to not be seen but not 0 - otherwise the horizontal lines on error bars don't get drawn

    mgpy_label = "MG5+Pythia8"
    hpp_label = "Herwig++"
    col_hpp = qgc.DY_COLOURS[2]

    data_hist = cu.get_from_tfile(single_mu_tfile, histname)
    mg_hist = cu.get_from_tfile(mg_dy_tfile, histname)
    hpp_hist = cu.get_from_tfile(hpp_dy_tfile, histname)

    # For use with tunfold binning, which just has bin indices as x values instead of physical values
    all_pt_bins = np.append(qgc.PT_UNFOLD_DICT['underflow_zpj_reco'][:-1], qgc.PT_UNFOLD_DICT['signal_zpj_reco'])
    all_pt_bins = np.append(all_pt_bins, 8000)  # the overflow bin
    # print(all_pt_bins)
    data_hist = tunfold_to_physical_bins(data_hist, all_pt_bins, divide_by_bin_width=True)
    mg_hist = tunfold_to_physical_bins(mg_hist, all_pt_bins, divide_by_bin_width=True)
    hpp_hist = tunfold_to_physical_bins(hpp_hist, all_pt_bins, divide_by_bin_width=True)

    # Scale to data
    mg_sf = data_hist.Integral()/mg_hist.Integral()
    mg_hist.Scale(mg_sf)
    hpp_hist.Scale(data_hist.Integral()/hpp_hist.Integral())

    # Absolute shifted variations
    exp_hist_up, exp_hist_down = None, None
    scale_hist_up, scale_hist_down = None, None
    pdf_hist_up, pdf_hist_down = None, None
    total_hist_up, total_hist_down = None, None

    # Calculate various systematic variations
    if any([show_total_systematics, show_individual_systematics, show_grouped_systematics]):
        # Calculate total experimental systematic uncertainty
        # ------------------------------------------------------------------
        for exp_dict in experimental_systematics:
            if isinstance(exp_dict['tfile'], str):
                exp_dict['tfile'] = cu.open_root_file(exp_dict['tfile'])
            hist = cu.get_from_tfile(exp_dict['tfile'], histname)
            hist = tunfold_to_physical_bins(hist, all_pt_bins, divide_by_bin_width=True)
            hist.Scale(mg_sf)
            exp_dict['hist'] = hist

        # create envelope of quadrature sum max/min per bin
        exp_hist_up = mg_hist.Clone("exp_hist_up")
        exp_hist_down = mg_hist.Clone("exp_hist_down")
        for ix in range(1, mg_hist.GetNbinsX()+1):
            nominal_bin = mg_hist.GetBinContent(ix)
            all_bin_diffs = [exp_dict['hist'].GetBinContent(ix) - nominal_bin
                             for exp_dict in experimental_systematics]

            positive_err = math.sqrt(sum([x**2 for x in all_bin_diffs if x > 0]))
            exp_hist_up.SetBinContent(ix, nominal_bin + positive_err)
            exp_hist_up.SetBinError(ix, 0)

            negative_err = math.sqrt(sum([x**2 for x in all_bin_diffs if x < 0]))
            exp_hist_down.SetBinContent(ix, nominal_bin - negative_err)
            exp_hist_down.SetBinError(ix, 0)

            # print(ix, 'nominal', nominal_bin)
            # print(ix, 'positive_err', positive_err, positive_err/nominal_bin)
            # print(ix, 'negative_err', negative_err, negative_err/nominal_bin)


        # Calculate total scale uncertainty
        # ------------------------------------------------------------------
        for scale_dict in scale_systematics:
            if isinstance(scale_dict['tfile'], str):
                scale_dict['tfile'] = cu.open_root_file(scale_dict['tfile'])
            hist = cu.get_from_tfile(scale_dict['tfile'], histname)
            hist = tunfold_to_physical_bins(hist, all_pt_bins, divide_by_bin_width=True)
            hist.Scale(mg_sf)
            scale_dict['hist'] = hist

        # create envelope of max/min per bin
        scale_hist_up = mg_hist.Clone("scale_hist_up")
        scale_hist_down = mg_hist.Clone("scale_hist_down")
        for ix in range(1, mg_hist.GetNbinsX()+1):
            nominal_bin = mg_hist.GetBinContent(ix)
            all_bin_values = [scale_dict['hist'].GetBinContent(ix)
                              for scale_dict in scale_systematics]
            scale_hist_up.SetBinContent(ix, max(all_bin_values))
            scale_hist_up.SetBinError(ix, 0)
            scale_hist_down.SetBinContent(ix, min(all_bin_values))
            scale_hist_down.SetBinError(ix, 0)


        # Calculate total PDF uncertainty
        # ------------------------------------------------------------------
        # Create all dicts first
        tfile = pdf_systematics[0]['tfile']
        if isinstance(tfile, str):
            tfile = cu.open_root_file(tfile)

        these_pdf_systematics = []
        num_vars = len(pdf_systematics[0]['variations'])
        for pdf_ind in pdf_systematics[0]['variations']:
            hist = cu.get_from_tfile(tfile, "ZPlusJets_QG_Unfold/hist_pt_reco_all_PDF_%d" % (pdf_ind))
            hist = tunfold_to_physical_bins(hist, all_pt_bins, divide_by_bin_width=True)
            hist.Scale(mg_sf)
            these_pdf_systematics.append(
                {
                    "label": "PDF_%d" % (pdf_ind),
                    "hist": hist,
                    "colour": cu.get_colour_seq(pdf_ind, num_vars)
                })

        # create up/down
        pdf_hist_up = mg_hist.Clone("pdf_hist_up")
        pdf_hist_down = mg_hist.Clone("pdf_hist_down")
        for ix in range(1, mg_hist.GetNbinsX()+1):
            nominal = mg_hist.GetBinContent(ix)
            values = [pdf_dict['hist'].GetBinContent(ix) for pdf_dict in these_pdf_systematics]
            rms = np.std(values, ddof=1)
            pdf_hist_up.SetBinContent(ix, nominal+rms)
            pdf_hist_up.SetBinError(ix, 0)
            pdf_hist_down.SetBinContent(ix, max(nominal-rms, 0))
            pdf_hist_down.SetBinError(ix, 0)
            # print(ix, 'nominal', nominal, 'rms', rms, 'rms/nominal', rms/nominal)
            # if rms > nominal:
            #     print("!!!!!")


        # Calculate total uncertainty
        # ------------------------------------------------------------------
        total_hist_up = mg_hist.Clone('total_hist_up')
        total_hist_down = mg_hist.Clone('total_hist_down')
        for ix in range(1, mg_hist.GetNbinsX()+1):
            nominal = mg_hist.GetBinContent(ix)
            sum_sq = sum([(up_hist.GetBinContent(ix) - nominal)**2
                          for up_hist in [exp_hist_up, scale_hist_up, pdf_hist_up]
                          if up_hist is not None])
            # print(ix, [(up_hist.GetBinContent(ix) - nominal)**2
            #           for up_hist in [exp_hist_up, scale_hist_up, pdf_hist_up]
            #           if up_hist is not None])
            total_hist_up.SetBinContent(ix, nominal+math.sqrt(sum_sq))
            total_hist_up.SetBinError(ix, 0)

            sum_sq = sum([(down_hist.GetBinContent(ix) - nominal)**2
                          for down_hist in [exp_hist_down, scale_hist_down, pdf_hist_down]
                          if down_hist is not None])
            total_hist_down.SetBinContent(ix, nominal - math.sqrt(sum_sq))
            total_hist_down.SetBinError(ix, 0)

    # Create entries for plot
    # --------------------------------------------------------------------------
    ref_hist = data_hist if subplot_vs_data else mg_hist
    entries = [
        # SINGLE MU DATA
        [
            data_hist,
            dict(line_color=ROOT.kBlack, line_width=data_line_width, fill_color=qgc.SINGLE_MU_COLOUR,
                 marker_color=ROOT.kBlack, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=msize*0.7,
                 label="Data",
                 subplot=None if subplot_vs_data else ref_hist)
        ],

        # MG5+PYTHIA8 MC
        [
            mg_hist,
            dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                 marker_color=qgc.DY_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                 label=mgpy_label,
                 subplot=ref_hist if subplot_vs_data else None)
        ],

        # HERWIG++
        [
            hpp_hist,
            dict(line_color=col_hpp, line_width=lw, line_style=2, fill_color=col_hpp,
                 marker_color=col_hpp, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                 label=hpp_label,
                 subplot=ref_hist)
        ]
    ]

    # ADD EXPERIMENTAL ENTRIES
    # ------------------------------------------------------------------
    if show_individual_systematics:
        for exp_dict in experimental_systematics:
            entries.append([exp_dict['hist'],
                            dict(line_color=exp_dict['colour'], line_width=1, fill_color=exp_dict['colour'],
                                 marker_color=exp_dict['colour'], marker_size=0, label=exp_dict['label'],
                                 line_style=1 if "up" in exp_dict['label'].lower() else 2,
                                 subplot=ref_hist)
                           ])

    if show_grouped_systematics:
        exp_col = ROOT.kRed+2
        exp_col2 = ROOT.kRed-2
        entries.extend([
            [
                exp_hist_up.Clone(),
                dict(line_color=exp_col, line_width=lw, line_style=1,
                     fill_color=exp_col,
                     marker_color=exp_col, marker_size=0, label="Exp systs up",
                     subplot=ref_hist)
            ],
            [
                exp_hist_down.Clone(),
                dict(line_color=exp_col2, line_width=lw, line_style=2,
                     fill_color=exp_col2,
                     marker_color=exp_col2, marker_size=0, label="Exp systs down",
                     subplot=ref_hist)
            ],
        ])

    n = mg_hist.GetNbinsX()
    bin_width = [mg_hist.GetBinWidth(i)/2. for i in range(1, n+1)]
    x = array('d', [mg_hist.GetBinCenter(i) for i in range(1, n+1)])
    y = array('d', [1 for i in range(1, n+1)])
    exl = array('d', [bin_width[i-1] for i in range(1, n+1)])
    exh = array('d', [bin_width[i-1] for i in range(1, n+1)])

    # exp_hist_up_ratio = exp_hist_up.Clone()
    # exp_hist_down_ratio = exp_hist_down.Clone()
    # exp_hist_up_ratio.Add(mg_hist, -1)
    # exp_hist_up_ratio.Divide(mg_hist)
    # exp_hist_down_ratio.Add(mg_hist, -1)
    # exp_hist_down_ratio.Divide(mg_hist)
    # exp_hist_down_ratio.Scale(-1)
    # exp_eyl = array('d', [exp_hist_down_ratio.GetBinContent(i) for i in range(1, n+1)])
    # exp_eyh = array('d', [exp_hist_up_ratio.GetBinContent(i) for i in range(1, n+1)])
    # exp_gr = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, exp_eyl, exp_eyh)

    # ADD SCALE ENTRIES
    # ------------------------------------------------------------------
    if show_individual_systematics:
        for scale_dict in scale_systematics:
            entries.append([scale_dict['hist'],
                            dict(line_color=scale_dict['colour'], line_width=lw, fill_color=scale_dict['colour'],
                                 marker_color=scale_dict['colour'], marker_size=0, label=scale_dict['label'],
                                 subplot=ref_hist)
                           ])

    if show_grouped_systematics:
        scale_col = ROOT.kGreen+2
        scale_col2 = scale_col
        entries.extend([
            [
                scale_hist_up.Clone(),
                dict(line_color=scale_col, line_width=lw, line_style=1,
                     fill_color=scale_col,
                     marker_color=scale_col, marker_size=0, label="Scale up",
                     subplot=ref_hist)
            ],
            [
                scale_hist_down.Clone(),
                dict(line_color=scale_col2, line_width=lw, line_style=2,
                     fill_color=scale_col2,
                     marker_color=scale_col2, marker_size=0, label="Scale down",
                     subplot=ref_hist)
            ],
        ])

    # scale_hist_up_ratio = scale_hist_up.Clone()
    # scale_hist_down_ratio = scale_hist_down.Clone()
    # scale_hist_up_ratio.Add(mg_hist, -1)
    # scale_hist_up_ratio.Divide(mg_hist)
    # scale_hist_down_ratio.Add(mg_hist, -1)
    # scale_hist_down_ratio.Divide(mg_hist)
    # scale_hist_down_ratio.Scale(-1)
    # scale_eyl = array('d', [scale_hist_down_ratio.GetBinContent(i) for i in range(1, n+1)])
    # scale_eyh = array('d', [scale_hist_up_ratio.GetBinContent(i) for i in range(1, n+1)])
    # scale_gr = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, scale_eyl, scale_eyh)

    # ADD PDF ENTRIES
    # ------------------------------------------------------------------
    if show_individual_systematics:
        for pdf_dict in these_pdf_systematics:
            entries.append(
                [
                    pdf_dict['hist'],
                    dict(label=pdf_dict['label'],
                         line_color=pdf_dict['colour'],
                         subplot=mg_hist)
                ]
            )

    if show_grouped_systematics:
        pdf_col = ROOT.kOrange-4
        pdf_col2 = pdf_col
        entries.extend([
            [
                pdf_hist_up.Clone(),
                dict(line_color=pdf_col, line_width=lw, line_style=1,
                     fill_color=pdf_col,
                     marker_color=pdf_col, marker_size=0, label="PDF up",
                     subplot=ref_hist)
            ],
            [
                pdf_hist_down.Clone(),
                dict(line_color=pdf_col2, line_width=lw, line_style=2,
                     fill_color=pdf_col2,
                     marker_color=pdf_col2, marker_size=0, label="PDF down",
                     subplot=ref_hist)
            ],
        ])

    # pdf_hist_up_ratio = pdf_hist_up.Clone()
    # pdf_hist_down_ratio = pdf_hist_down.Clone()
    # pdf_hist_up_ratio.Add(mg_hist, -1)
    # pdf_hist_up_ratio.Divide(mg_hist)
    # pdf_hist_down_ratio.Add(mg_hist, -1)
    # pdf_hist_down_ratio.Divide(mg_hist)
    # pdf_hist_down_ratio.Scale(-1)
    # pdf_eyl = array('d', [pdf_hist_down_ratio.GetBinContent(i) for i in range(1, n+1)])
    # pdf_eyh = array('d', [pdf_hist_up_ratio.GetBinContent(i) for i in range(1, n+1)])
    # pdf_gr = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, pdf_eyl, pdf_eyh)

    # ADD TOTAL
    # ------------------------------------------------------------------
    if show_grouped_systematics:
        total_col = ROOT.kBlack
        total_col2 = total_col
        entries.extend([
            [
                total_hist_up.Clone(),
                dict(line_color=total_col, line_width=lw, line_style=1,
                     fill_color=total_col,
                     marker_color=total_col, marker_size=0, label="Total up",
                     subplot=ref_hist)
            ],
            [
                total_hist_down.Clone(),
                dict(line_color=total_col2, line_width=lw, line_style=2,
                     fill_color=total_col2,
                     marker_color=total_col2, marker_size=0, label="Total down",
                     subplot=ref_hist)
            ],
        ])

    if show_total_systematics:
        total_hist_up_ratio = total_hist_up.Clone()
        total_hist_down_ratio = total_hist_down.Clone()
        total_hist_up_ratio.Add(mg_hist, -1)
        total_hist_up_ratio.Divide(mg_hist)
        total_hist_down_ratio.Add(mg_hist, -1)
        total_hist_down_ratio.Divide(mg_hist)
        total_hist_down_ratio.Scale(-1)
        total_eyl = array('d', [total_hist_down_ratio.GetBinContent(i) for i in range(1, n+1)])
        total_eyh = array('d', [total_hist_up_ratio.GetBinContent(i) for i in range(1, n+1)])
        total_gr = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, total_eyl, total_eyh)

    # print("data_hist.Integral('width'):", data_hist.Integral("width"))

    radius, pus = cu.get_jet_config_from_dirname(workdir)
    jet_str = "AK%s" % (radius.upper())
    title = "{jet_algo}\n{region_label}".format(jet_algo=jet_str,
                                                  region_label=qgc.ZpJ_LABEL)
    subplot_title = "* / %s" % entries[1][1]['label']
    if subplot_vs_data:
        subplot_title = qgc.SIM_DATA_STR
    paper_str = "_paper" if not is_preliminary else ""
    do_jet_pt_plot(entries,
                   output_filename=os.path.join(workdir, "data_mc_jet_pt/ZPlusJets/jet_pt%s.%s" % (paper_str, OUTPUT_FMT)),
                   rebin=1,
                   xlim=(30, qgc.PT_UNFOLD_DICT['signal_zpj_gen'][-1]),
                   ylim=(5E-3, 1E7),
                   title=title,
                   data_first=True,
                   subplot_title=subplot_title,
                   subplot_limits=(0, 2.5),
                   # subplot_limits=(0.5, 1.5),
                   normalise_hists=False,
                   # experimental_syst=exp_gr,
                   # scale_syst=scale_gr,
                   # pdf_syst=pdf_gr,
                   experimental_syst=None,
                   scale_syst=None,
                   pdf_syst=None,
                   total_syst=total_gr if show_total_systematics else None,
                   lumi=cu.get_lumi_str(do_dijet=False, do_zpj=True),
                   is_preliminary=is_preliminary
                   )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdirs',
                        nargs='+',
                        help='Workdir(s) with ROOT files to process. '
                             'Several dirs can be specified here, separated by a space.')
    parser.add_argument("--paper",
                        action='store_true',
                        help='Don\'t add "Preliminary" to plots')
    args = parser.parse_args()

    for workdir in args.workdirs:
        do_dijet_pt_plots(workdir,
                          show_total_systematics=True,
                          show_grouped_systematics=False,
                          show_individual_systematics=False,
                          is_preliminary=not args.paper)
        do_zpj_pt_plots(workdir,
                        show_total_systematics=True,
                        show_grouped_systematics=False,
                        show_individual_systematics=False,
                        is_preliminary=not args.paper)

    sys.exit(0)
