#!/usr/bin/env python


"""
Do all the unfolding plots: per pT bin, per lambda bin, summary plot
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
import math
from array import array
import pandas as pd

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, unfolder_from_tdir, HistBinChopper, unpack_unfolding_root_file
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import get_dijet_config, get_zpj_config


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


def setup_regions(args):
    regions = []
    if args.doDijetCentral:
        dijet_region_central_dict = get_dijet_config(args.source, central=True, groomed=False)
        regions.append(dijet_region_central_dict)

    if args.doDijetForward:
        dijet_region_forward_dict = get_dijet_config(args.source, central=False, groomed=False)
        regions.append(dijet_region_forward_dict)

    if args.doDijetCentralGroomed:
        dijet_region_central_groomed_dict = get_dijet_config(args.source, central=True, groomed=True)
        regions.append(dijet_region_central_groomed_dict)

    if args.doDijetForwardGroomed:
        dijet_region_forward_groomed_dict = get_dijet_config(args.source, central=False, groomed=True)
        regions.append(dijet_region_forward_groomed_dict)

    if args.doZPJ:
        zpj_region_dict = get_zpj_config(args.source, groomed=False)
        regions.append(zpj_region_dict)

    if args.doZPJGroomed:
        zpj_region_groomed_dict = get_zpj_config(args.source, groomed=True)
        regions.append(zpj_region_groomed_dict)
    return regions


class Setup(object):
    """Loads of common consts, useful for plotting etc"""

    def __init__(self, jet_algo, region, angle, output_dir='.', has_data=False, is_ave_pt_binning=False):
        self.jet_algo = jet_algo
        self.region = region
        self.pt_var_str = "#LT p_{T}^{jet} #GT" if is_ave_pt_binning else "p_{T}^{jet}"
        self.pt_str = self.pt_var_str + " [GeV]"
        self.has_data = has_data
        self.angle = angle
        angle_prepend = "Groomed " if "groomed" in region['name'] else ""
        this_angle_name = angle.name
        if (angle_prepend != ""
            and 'LHA' not in this_angle_name
            and "_{T}" not in this_angle_name
            and "PUPPI" not in this_angle_name):
            # lower case if Groomed..., but be careful of e.g. pTD, LHA
            this_angle_name = this_angle_name[0].lower() + this_angle_name[1:]
        # for plot axis titles
        self.angle_str = "{prepend}{name} ({lambda_str})".format(prepend=angle_prepend,
                                                                 name=this_angle_name,
                                                                 lambda_str=angle.lambda_str)
        # self.particle_title = "Particle-level " + self.angle_str
        self.particle_title = self.angle_str
        self.detector_title = "Detector-level " + self.angle_str
        self.pt_bin_normalised_differential_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_unnormalised_differential_label = "#frac{1}{dN/dp_{T}} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str)  # FIXME
        self.lambda_bin_normalised_differential_label = "#frac{1}{d#sigma/d%s} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)
        self.lambda_bin_unnormalised_differential_label = "#frac{1}{dN/d%s} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)  # FIXME
        self.output_dir = output_dir
        self.output_fmt = 'pdf'
        self.append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name


PLOT_COLOURS = dict(
    gen_colour=ROOT.kRed,
    unfolded_basic_colour=ROOT.kAzure+7,
    unfolded_stat_colour=ROOT.kAzure+7,
    unfolded_total_colour=ROOT.kBlack,
    unfolded_unreg_colour=ROOT.kViolet+2,
    # alt_gen_colour=ROOT.kOrange-3,
    alt_gen_colour=ROOT.kViolet+1,
    alt_unfolded_colour=ROOT.kOrange-3,
    alt_reco_colour=ROOT.kOrange-3,
    # reco_mc_colour=ROOT.kGreen+2,
    # reco_mc_colour=ROOT.kAzure-7,
    # reco_data_colour=ROOT.kRed,
    reco_data_colour=ROOT.kBlack,
    reco_mc_colour=ROOT.kRed+1,
    reco_unfolding_input_colour=ROOT.kRed,
    reco_folded_unfolded_colour=ROOT.kAzure+1,
    reco_folded_mc_truth_colour=ROOT.kGreen+2,
)


def calc_auto_xlim(entries):
    """Figure out x axis range that includes all non-0 bins from all Contributions"""
    if len(entries) == 0:
        return None
    if not isinstance(entries[0], Contribution):
        raise TypeError("calc_auto_xlim: `entries` should be a list of Contributions")

    x_min = 9999999999
    x_max = -9999999999
    for ent in entries:
        obj = ent.obj
        if not isinstance(obj, ROOT.TH1):
            raise TypeError("Cannot handle obj in calc_auto_xlim")
        xax = obj.GetXaxis()
        nbins = obj.GetNbinsX()
        found_min = False
        for i in range(1, nbins+1):
            val = obj.GetBinContent(i)
            if not found_min:
                # find the first non-empty bin
                if val != 0:
                    x_min = min(xax.GetBinLowEdge(i), x_min)
                    found_min = True
            else:
                # find the first empty bin
                if val == 0:
                    x_max = max(xax.GetBinLowEdge(i+1), x_max)
                    break
    if x_max > x_min:
        return (x_min, x_max)
    else:
        # default x max
        return None


# FIXME: generalise this and LambdaBinnedPlotter into one generic BinnedPlotter?
# Although each has different set of plots, so not easy/possible
class GenPtBinnedPlotter(object):
    def __init__(self, setup, bins, hist_bin_chopper, unfolder):
        self.setup = setup
        self.region = setup.region  # just to make life easier
        self.bins = bins
        self.hist_bin_chopper = hist_bin_chopper

        self.line_width = 2
        self.plot_colours = PLOT_COLOURS
        self.pt_bin_plot_args = dict(
            what="hist",
            xtitle=self.setup.particle_title,
            has_data=self.setup.has_data,
            subplot_type='ratio',
            subplot_title="Unfolded / Gen",
            subplot_limits=(0, 2) if self.setup.has_data else (0.75, 1.25),
        )
        self.unfolder = unfolder

    @staticmethod
    def _modify_plot(this_plot):
        this_plot.legend.SetX1(0.6)
        this_plot.legend.SetY1(0.68)
        this_plot.legend.SetX2(0.98)
        this_plot.legend.SetY2(0.9)
        this_plot.left_margin = 0.16

    @staticmethod
    def check_entries(entries, message=""):
        """Check that at least 1 Contribution has something in it, that is non-zero"""
        has_entries = [c.obj.GetEntries() > 0 for c in entries]
        if not any(has_entries):
            if message:
                print("Skipping 0 entries (%s)" % (message))
            else:
                print("Skipping 0 entries")
            return False

        max_bin = max([c.obj.GetMaximum() for c in entries])
        min_bin = min([c.obj.GetMinimum() for c in entries])
        if max_bin == min_bin:
            if message:
                print("Skipping min=max hists (%s)" % (message))
            else:
                print("Skipping min=max hists")
            return False

        return True

    def get_pt_bin_title(self, bin_edge_low, bin_edge_high):
        title = (("{jet_algo}\n"
                  "{region_label} region\n"
                  "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                 .format(
                    jet_algo=self.setup.jet_algo,
                    region_label=self.region['label'],
                    pt_str=self.setup.pt_var_str,
                    bin_edge_low=bin_edge_low,
                    bin_edge_high=bin_edge_high
                ))
        return title

    def plot_unfolded_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin('unfolded_stat_err', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin('unfolded', ibin, binning_scheme='generator')

            # unnormalised version
            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]
            if not self.check_entries(entries, "plot_unfolded_unnormalised %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle="N",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_unnormalised_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]
            if not self.check_entries(entries, "plot_unfolded_normalised_pt_bin %d" % (ibin)):
                return

            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_alt_truth_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin, binning_scheme='generator')
            alt_mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', ibin, binning_scheme='generator')

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(alt_mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['alt_mc_label']),
                             line_color=self.plot_colours['alt_gen_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_colours['alt_gen_colour'], marker_size=0,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]
            if not self.check_entries(entries, "plot_unfolded_with_alt_truth_normalised_pt_bin %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_%s_alt_truth_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_unreg_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin, binning_scheme='generator')
            unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unreg_unfolded', ibin, binning_scheme='generator')

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unreg_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = 0) (total unc.)",
                             line_color=self.plot_colours['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_unreg_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]
            if not self.check_entries(entries, "plot_unfolded_with_unreg_normalised_pt_bin %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_%s_with_unreg_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_alt_response_normalised(self, alt_unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin, binning_scheme='generator')
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            alt_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_unfolded_stat_err', ibin, binning_scheme='generator')
            # alt_mc_gen_hist_gin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', ibin, binning_scheme='generator')

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                # Contribution(alt_mc_gen_hist_bin,
                #              label="Generator (%s)" % (self.region['alt_mc_label']),
                #              line_color=alt_gen_colour, line_width=self.line_width, line_style=2,
                #              marker_color=alt_gen_colour, marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(alt_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (alt_unfolder.tau, self.region['alt_mc_label']),
                             line_color=self.plot_colours['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['alt_unfolded_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]
            if not self.check_entries(entries, "plot_unfolded_with_unreg_normalised_pt_bin %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_%s_alt_response_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_alt_response_truth_normalised(self, alt_unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin, binning_scheme='generator')
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            alt_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_unfolded_stat_err', ibin, binning_scheme='generator')
            alt_mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', ibin, binning_scheme='generator')

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(alt_mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['alt_mc_label']),
                             line_color=self.plot_colours['alt_gen_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_colours['alt_gen_colour'], marker_size=0,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(alt_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (alt_unfolder.tau, self.region['alt_mc_label']),
                             line_color=self.plot_colours['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['alt_unfolded_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]
            if not self.check_entries(entries, "plot_unfolded_with_unreg_normalised_pt_bin %d" % (ibin)):
                return
            this_plot_args = {k:v for k,v in self.pt_bin_plot_args.items()}
            this_plot_args['subplot_title'] = '#splitline{* / Gen}{(%s)}' % (self.region['mc_label'])
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **this_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_%s_alt_response_truth_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_model_systs_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            syst_entries = []
            for syst_dict in self.region['model_systematics']:
                syst_unfolder = syst_dict['unfolder']
                syst_label = syst_dict['label']
                syst_label_no_spaces = cu.no_space_str(syst_dict['label'])

                self.hist_bin_chopper.add_obj('model_syst_%s_unfolded' % (syst_label_no_spaces), syst_unfolder.unfolded)
                self.hist_bin_chopper.add_obj('model_syst_%s_hist_truth' % (syst_label_no_spaces), syst_unfolder.hist_truth)

                syst_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('model_syst_%s_unfolded' % (syst_label_no_spaces), ibin, binning_scheme='generator')
                syst_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('model_syst_%s_hist_truth' % (syst_label_no_spaces), ibin, binning_scheme='generator')

                syst_entries.extend([
                    Contribution(syst_gen_hist_bin,
                                 label="Generator (%s)" % (syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                                 marker_color=syst_dict['colour'], marker_size=0),
                    Contribution(syst_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (stat. err) (%s)" % (syst_unfolder.tau, syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=syst_gen_hist_bin),
                ])

            # add nominal ones last
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin, binning_scheme='generator')

            syst_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])
            if not self.check_entries(syst_entries, "plot_unfolded_with_model_systs_normalised_pt_bin %d" % (ibin)):
                return
            plot = Plot(syst_entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            if len(syst_entries) > 4:
                # plot.legend.SetX1(0.53)
                plot.legend.SetY1(0.65)
                plot.y_padding_max_linear = 1.8
            if len(syst_entries) > 6:
                plot.legend.SetNColumns(2)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_%s_syst_model_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_pdf_systs_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            pdf_entries = []
            for pdf_dict in self.region['pdf_systematics']:
                pdf_unfolder = pdf_dict['unfolder']
                pdf_label = pdf_dict['label']
                pdf_label_no_spaces = cu.no_space_str(pdf_dict['label'])

                self.hist_bin_chopper.add_obj('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), pdf_unfolder.unfolded)
                self.hist_bin_chopper.add_obj('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), pdf_unfolder.hist_truth)

                pdf_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), ibin, binning_scheme='generator')
                pdf_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), ibin, binning_scheme='generator')

                pdf_entries.extend([
                    Contribution(pdf_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (stat. err) (%s)" % (pdf_unfolder.tau, pdf_label),
                                 line_color=pdf_dict['colour'], line_width=1, line_style=1,
                                 marker_color=pdf_dict['colour'], marker_size=0,
                                 subplot=pdf_gen_hist_bin),
                    Contribution(pdf_gen_hist_bin,
                                 label="Generator (%s)" % (pdf_label),
                                 line_color=pdf_dict['colour'], line_width=1, line_style=2,
                                 marker_color=pdf_dict['colour'], marker_size=0),
                ])

            # add nominal ones last
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin, binning_scheme='generator')

            pdf_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])
            if not self.check_entries(pdf_entries, "plot_unfolded_with_model_systs_normalised_pt_bin %d" % (ibin)):
                return
            plot = Plot(pdf_entries,
                       ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            if len(pdf_entries) > 5:
                plot.legend.SetNColumns(2)
            if len(pdf_entries) > 15:
                plot.legend.SetNColumns(3)
            ROOT.gStyle.SetPalette(55)
            plot.plot("NOSTACK E1 PLC PMC")
            plot.save("%s/unfolded_%s_pdf_model_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
            ROOT.gStyle.SetPalette(ROOT.kViridis)

    def plot_uncertainty_shifts_normalised(self):
        """Do plots of fractional uncertainty shifts on *normalised* unfolded distribution"""
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            entries = []
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin, binning_scheme='generator')
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_rsp_err', ibin, binning_scheme='generator')

            for syst_dict in self.region['experimental_systematics']:
                # For each systematic, get the normalised shift and hence fraction
                obj_name = 'syst_shift_%s' % cu.no_space_str(syst_dict['label'])
                syst_unfolded_fraction = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(obj_name, ibin, binning_scheme='generator').Clone()
                syst_unfolded_fraction.Divide(unfolded_hist_bin_total_errors)
                # Set to abs values so can plot log
                for i in range(1, syst_unfolded_fraction.GetNbinsX()+1):
                    syst_unfolded_fraction.SetBinContent(i, abs(syst_unfolded_fraction.GetBinContent(i)))
                    syst_unfolded_fraction.SetBinError(i, 0)
                c = Contribution(syst_unfolded_fraction,
                                 label=syst_dict['label'],
                                 line_color=syst_dict['colour'],
                                 line_style=syst_dict.get('linestyle', 1),
                                 line_width=self.line_width,
                                 marker_size=0,
                                 marker_color=syst_dict['colour'])
                entries.append(c)

            # Now create total & stat. unc.or hists
            h_stat = unfolded_hist_bin_stat_errors.Clone()
            h_syst = unfolded_hist_bin_rsp_errors.Clone()
            h_total = unfolded_hist_bin_total_errors.Clone()
            for i in range(1, h_stat.GetNbinsX()+1):
                if unfolded_hist_bin_total_errors.GetBinContent(i) > 0:
                    h_stat.SetBinContent(i, unfolded_hist_bin_stat_errors.GetBinError(i) / unfolded_hist_bin_total_errors.GetBinContent(i))
                    h_syst.SetBinContent(i, unfolded_hist_bin_rsp_errors.GetBinError(i) / unfolded_hist_bin_total_errors.GetBinContent(i))
                    h_total.SetBinContent(i, unfolded_hist_bin_total_errors.GetBinError(i) / unfolded_hist_bin_total_errors.GetBinContent(i))
                else:
                    h_stat.SetBinContent(i, 0)
                    h_syst.SetBinContent(i, 0)
                    h_total.SetBinContent(i, 0)
                h_stat.SetBinError(i, 0)
                h_syst.SetBinError(i, 0)
                h_total.SetBinError(i, 0)
            c_stat = Contribution(h_stat,
                                 label="Input stats",
                                 line_color=ROOT.kRed,
                                 line_style=3,
                                 line_width=3,
                                 marker_size=0,
                                 marker_color=ROOT.kRed,
                                 )
            c_syst = Contribution(h_syst,
                                 label="Response matrix stats",
                                 line_color=ROOT.kGray+2,
                                 line_style=3,
                                 line_width=3,
                                 marker_size=0,
                                 marker_color=ROOT.kGray+2,
                                 )
            c_tot = Contribution(h_total,
                                 label="Total uncertainty",
                                 line_color=ROOT.kBlack,
                                 line_style=1,
                                 line_width=3,
                                 marker_size=0,
                                 marker_color=ROOT.kBlack,
                                 )
            entries.extend([c_stat, c_syst, c_tot])
            if not self.check_entries(entries, "plot_uncertainty_shifts_normalised %d" % (ibin)):
                return
            plot = Plot(entries,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        xtitle=self.setup.particle_title,
                        ytitle="| Fractional shift on normalised distribution |",
                        has_data=self.setup.has_data)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            plot.legend.SetNColumns(2)
            plot.left_margin = 0.16
            plot.y_padding_max_linear = 1.4
            plot.plot("NOSTACK HIST")
            output_filename = "%s/unfolded_systs_%s_bin_%d.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt)
            plot.save(output_filename)

            plot.y_padding_max_log = 50
            plot.set_logy(do_more_labels=False)
            plot.get_modifier().SetMinimum(1E-4)
            log_filename, ext = os.path.splitext(output_filename)
            plot.save(log_filename+"_log"+ext)

    def plot_unfolded_with_exp_systs_normalised(self):
        """Plot shifted unfolded normalised distributions for each syst"""
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin, binning_scheme='generator')
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_rsp_err', ibin, binning_scheme='generator')

            def _remove_error_bars(h):
                for i in range(1, h.GetNbinsX()+1):
                    h.SetBinError(i, 0)

            entries = []
            for syst_dict in self.region['experimental_systematics']:
                syst_label_no_spaces = cu.no_space_str(syst_dict['label'])
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('syst_shifted_%s_unfolded' % (syst_label_no_spaces), ibin, binning_scheme='generator')
                _remove_error_bars(syst_unfolded_hist_bin)
                c = Contribution(syst_unfolded_hist_bin,
                                 label=syst_dict['label'],
                                 line_color=syst_dict['colour'], line_width=self.line_width,
                                 line_style=2 if 'down' in syst_dict['label'].lower() else 1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=unfolded_hist_bin_stat_errors)
                entries.append(c)

            entries.append(
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], marker_style=20, marker_size=0.75),
            )
            entries.append(
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75),
            )
            if not self.check_entries(entries, "plot_unfolded_with_exp_systs_normalised %d" % ibin):
                return
            plot = Plot(entries,
                        xtitle=self.setup.particle_title,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='Syst / nominal',
                        subplot_limits=(0.75, 1.25))
            self._modify_plot(plot)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            if len(entries) > 5: plot.legend.SetNColumns(2)
            # plot.subplot_limits = (0.9, 1.1)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_syst_variations_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_exp_syst_variation_normalised(self):
        """Plot varation / central value on normalised hists
        (basically the subplot from plot_unfolded_with_exp_systs_normalised)
        """
        def _convert_error_bars_to_error_ratio_hist(h, direction=1):
            """Create hist with bin content = error shift / bin value, 0 error"""
            h_new = h.Clone(h.GetName() + cu.get_unique_str())
            for ibin in range(1, h_new.GetNbinsX()+1):
                if h.GetBinContent(ibin) > 0:
                    h_new.SetBinContent(ibin, 1+(direction*(h.GetBinError(ibin) / h.GetBinContent(ibin))))
                else:
                    h_new.SetBinContent(ibin, 99999 if direction > 0 else -99999)
                h_new.SetBinError(ibin, 0)
            return h_new

        def _convert_syst_shift_to_error_ratio_hist(h_syst, h_nominal):
            """Create h_syst / h_nominal without error bars"""
            h_new = h_syst.Clone(h_syst.GetName() + cu.get_unique_str())
            h_new.Divide(h_nominal)
            for ibin in range(1, h_new.GetNbinsX()+1):
                h_new.SetBinError(ibin, 0)
            return h_new

        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):

            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin, binning_scheme='generator')
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_rsp_err', ibin, binning_scheme='generator')

            entries = []
            for syst_dict, mark in zip(self.region['experimental_systematics'], cu.Marker().cycle(cycle_filling=True)):
                syst_label_no_spaces = cu.no_space_str(syst_dict['label'])
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('syst_shifted_%s_unfolded' % (syst_label_no_spaces), ibin, binning_scheme='generator')
                this_syst_hist = _convert_syst_shift_to_error_ratio_hist(syst_unfolded_hist_bin, unfolded_hist_bin_total_errors)
                c = Contribution(this_syst_hist,
                                 label=syst_dict['label'],
                                 line_color=syst_dict['colour'],
                                 line_width=0,
                                 # line_width=self.line_width,
                                 # line_style=2 if 'down' in syst_dict['label'].lower() else 1,
                                 marker_color=syst_dict['colour'], marker_size=1.25,
                                 marker_style=mark)
                entries.append(c)

            entries.extend([
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_total_errors),
                             label="Total uncertainty",
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_colours['unfolded_total_colour'], marker_style=20, marker_size=0,
                             fill_style=0, fill_color=15),
                # Add in the -ve side, but no label as we don't want it in the legend
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_total_errors, -1),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_colours['unfolded_total_colour'], marker_style=20, marker_size=0,
                             fill_style=0, fill_color=15),
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_stat_errors),
                             label="Input stats",
                             line_color=ROOT.kRed, line_width=self.line_width, line_style=3,
                             marker_color=ROOT.kRed, marker_style=20, marker_size=0,
                             fill_style=0, fill_color=13),
                # Add in the -ve side, but no label as we don't want it in the legend
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_stat_errors, -1),
                             line_color=ROOT.kRed, line_width=self.line_width, line_style=3,
                             marker_color=ROOT.kRed, marker_style=20, marker_size=0,
                             fill_style=0, fill_color=13),
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_rsp_errors),
                             label="Response matrix stats",
                             line_color=ROOT.kGray+2, line_width=self.line_width, line_style=3,
                             marker_color=ROOT.kGray+2, marker_style=20, marker_size=0,
                             fill_style=0, fill_color=13),
                # Add in the -ve side, but no label as we don't want it in the legend
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_rsp_errors, -1),
                             line_color=ROOT.kGray+2, line_width=self.line_width, line_style=3,
                             marker_color=ROOT.kGray+2, marker_style=20, marker_size=0,
                             fill_style=0, fill_color=13),
            ])

            if not self.check_entries(entries, "plot_exp_syst_variation_normalised %d" % ibin):
                return
            xlim = calc_auto_xlim(entries)
            plot = Plot(entries,
                        xtitle=self.setup.particle_title,
                        ytitle='Variation / nominal (on normalised distribution)',
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        ylim=(0.7, 1.5),
                        xlim=xlim,
                        subplot_type=None)
            self._modify_plot(plot)
            plot.default_canvas_size = (800, 700)
            plot.legend.SetX1(0.5)
            plot.legend.SetY1(0.65)
            plot.legend.SetX2(0.93)
            plot.legend.SetY2(0.88)
            if len(entries) > 5: plot.legend.SetNColumns(2)
            plot.plot("NOSTACK P L") # hard to get one that is points for systs, and line for stats
            plot.main_pad.cd()
            if xlim is not None:
                line = ROOT.TLine(xlim[0], 1, xlim[1], 1)
            else:
                xmin = unfolded_hist_bin_stat_errors.GetXaxis().GetBinLowEdge(1)
                xmax = unfolded_hist_bin_stat_errors.GetXaxis().GetBinLowEdge(unfolded_hist_bin_stat_errors.GetNbinsX()+1)
                line = ROOT.TLine(xmin, 1, xmax, 1)
            line.SetLineColor(ROOT.kGray+2)
            line.SetLineColor(ROOT.kBlack)
            line.Draw()
            plot.save("%s/unfolded_syst_variations_vs_nominal_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))


    def plot_detector_normalised(self, alt_detector=None):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_gen_binning_bg_subtracted', ibin, binning_scheme='generator')
            mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_gen_binning_bg_subtracted', ibin, binning_scheme='generator')

            entries = [
                Contribution(mc_hist_bin,
                             label="MC (bg-subtracted) [%s]" % (self.setup.region['mc_label']) if alt_detector else "MC (bg-subtracted)",
                             line_color=self.plot_colours['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_mc_colour'], marker_size=0,
                             subplot=input_hist_bin if alt_detector else None),
            ]
            if alt_detector:
                self.hist_bin_chopper.add_obj("alt_hist_mc_reco_gen_binning_bg_subtracted", alt_detector)
                alt_mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_mc_reco_gen_binning_bg_subtracted', ibin, binning_scheme='generator')
                entries.append(Contribution(alt_mc_hist_bin,
                                            label="MC (bg-subtracted) [%s]" % self.setup.region['alt_mc_label'],
                                            line_color=self.plot_colours['alt_reco_colour'], line_width=self.line_width,
                                            marker_color=self.plot_colours['alt_reco_colour'], marker_size=0,
                                            subplot=input_hist_bin))

            entries.append(
                Contribution(input_hist_bin,
                             label="Data (bg-subtracted)",
                             line_color=self.plot_colours['reco_data_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_data_colour'], marker_style=20, marker_size=0.75,
                             subplot=None if alt_detector else mc_hist_bin)
            )
            if not self.check_entries(entries, "plot_detector_normalised %d" % (ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='MC / Data' if alt_detector else 'Data / MC',
                        subplot_limits=(0.75, 1.25) if alt_detector else (0.75, 1.25)
                        )
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/detector_gen_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))


# ================================================================================


class GenLambdaBinnedPlotter(object):
    def __init__(self, setup, bins, hist_bin_chopper, unfolder):
        self.setup = setup
        self.region = setup.region
        self.bins = bins
        self.hist_bin_chopper = hist_bin_chopper

        self.line_width = 2
        self.plot_colours = PLOT_COLOURS
        self.lambda_bin_plot_args = dict(
            what="hist",
            xtitle=self.setup.pt_str,
            has_data=self.setup.has_data,
            subplot_type='ratio',
            subplot_title="Unfolded / Gen",
            subplot_limits=(0, 2) if self.setup.has_data else (0.75, 1.25),
        )
        self.unfolder = unfolder

    @staticmethod
    def _modify_plot(this_plot):
        this_plot.legend.SetX1(0.6)
        this_plot.legend.SetY1(0.68)
        this_plot.legend.SetX2(0.98)
        this_plot.legend.SetY2(0.9)
        this_plot.left_margin = 0.16
        this_plot.y_padding_max_log = 5000 # space for title

    @staticmethod
    def check_entries(entries, message=""):
        """Check that at least 1 Contribution has something in it"""
        has_entries = [c.obj.GetEntries() > 0 for c in entries]
        if not any(has_entries):
            if message:
                print("Skipping 0 entries (%s)" % (message))
            else:
                print("Skipping 0 entries")
            return False

        max_bin = max([c.obj.GetMaximum() for c in entries])
        min_bin = min([c.obj.GetMinimum() for c in entries])
        if max_bin == min_bin:
            if message:
                print("Skipping min=max hists (%s)" % (message))
            else:
                print("Skipping min=max hists")
            return False

        return True

    def get_lambda_bin_title(self, bin_edge_low, bin_edge_high):
        title = (("{jet_algo}\n"
                  "{region_label} region\n"
                  "{bin_edge_low:g} < {angle_str} < {bin_edge_high:g}")
                 .format(
                    jet_algo=self.setup.jet_algo,
                    region_label=self.region['label'],
                    angle_str=self.setup.angle.name,
                    bin_edge_low=bin_edge_low,
                    bin_edge_high=bin_edge_high
                ))
        return title

    def plot_unfolded_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')

            # unnormalised version
            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]
            if not self.check_entries(entries, "plot_unfolded_unnormalised %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle="N",
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            plot.save("%s/unfolded_unnormalised_%s_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_unreg_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')
            unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unreg_unfolded', ibin, binning_scheme='generator')

            # unnormalised version
            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unreg_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = 0) (total unc.)",
                             line_color=self.plot_colours['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_unreg_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]
            if not self.check_entries(entries, "plot_unfolded_with_unreg_unnormalised %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle="N",
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            plot.save("%s/unfolded_unnormalised_%s_with_unreg_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_alt_response_unnormalised(self, alt_unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            alt_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('alt_unfolded_stat_err', ibin, binning_scheme='generator')
            alt_mc_gen_hist_gin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('alt_hist_truth', ibin, binning_scheme='generator')

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                # Contribution(alt_mc_gen_hist_bin,
                #              label="Generator (%s)" % (self.region['alt_mc_label']),
                #              line_color=alt_gen_colour, line_width=self.line_width, line_style=2,
                #              marker_color=alt_gen_colour, marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(alt_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (alt_unfolder.tau, self.region['alt_mc_label']),
                             line_color=self.plot_colours['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['alt_unfolded_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]
            if not self.check_entries(entries, "plot_unfolded_with_alt_response_unnormalised %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.lambda_bin_unnormalised_differential_label,
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            plot.save("%s/unfolded_unnormalised_%s_alt_response_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_model_systs_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            syst_entries = []
            for syst_dict in self.region['model_systematics']:
                syst_unfolder = syst_dict['unfolder']
                syst_label = syst_dict['label']
                syst_label_no_spaces = cu.no_space_str(syst_dict['label'])

                self.hist_bin_chopper.add_obj('model_syst_%s_unfolded' % (syst_label_no_spaces), syst_unfolder.unfolded)
                self.hist_bin_chopper.add_obj('model_syst_%s_hist_truth' % (syst_label_no_spaces), syst_unfolder.hist_truth)

                syst_unfolded_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('model_syst_%s_unfolded' % (syst_label_no_spaces), ibin, binning_scheme='generator')
                syst_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('model_syst_%s_hist_truth' % (syst_label_no_spaces), ibin, binning_scheme='generator')

                syst_entries.extend([
                    Contribution(syst_gen_hist_bin,
                                 label="Generator (%s)" % (syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                                 marker_color=syst_dict['colour'], marker_size=0),
                    Contribution(syst_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (total unc.) (%s)" % (syst_unfolder.tau, syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=syst_gen_hist_bin),
                ])

            # add nominal ones last
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')

            syst_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])
            if not self.check_entries(syst_entries, "plot_unfolded_with_model_systs_unnormalised %d" % (ibin)):
                return
            plot = Plot(syst_entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.legend.SetX1(0.55)
            plot.legend.SetY1(0.72)
            plot.legend.SetX2(0.97)
            plot.legend.SetY2(0.88)
            if len(syst_entries) > 4:
                # plot.legend.SetX1(0.53)
                plot.legend.SetY1(0.65)
                plot.y_padding_max_linear = 1.8
            if len(syst_entries) > 6:
                plot.legend.SetNColumns(2)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            plot.save("%s/unfolded_unnormalised_%s_syst_model_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_pdf_systs_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            pdf_entries = []
            for pdf_dict in self.region['pdf_systematics']:
                pdf_unfolder = pdf_dict['unfolder']
                pdf_label = pdf_dict['label']
                pdf_label_no_spaces = cu.no_space_str(pdf_dict['label'])

                self.hist_bin_chopper.add_obj('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), pdf_unfolder.unfolded)
                self.hist_bin_chopper.add_obj('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), pdf_unfolder.hist_truth)

                pdf_unfolded_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), ibin, binning_scheme='generator')
                pdf_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), ibin, binning_scheme='generator')

                pdf_entries.extend([
                    Contribution(pdf_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (total unc.) (%s)" % (pdf_unfolder.tau, pdf_label),
                                 line_color=pdf_dict['colour'], line_width=1, line_style=1,
                                 marker_color=pdf_dict['colour'], marker_size=0,
                                 subplot=pdf_gen_hist_bin),
                    Contribution(pdf_gen_hist_bin,
                                 label="Generator (%s)" % (pdf_label),
                                 line_color=pdf_dict['colour'], line_width=1, line_style=2,
                                 marker_color=pdf_dict['colour'], marker_size=0),
                ])

            # add nominal ones last
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')

            pdf_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])
            if not self.check_entries(pdf_entries, "plot_unfolded_with_pdf_systs_unnormalised %d" % (ibin)):
                return
            plot = Plot(pdf_entries,
                        ytitle="N",
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.legend.SetX1(0.55)
            plot.legend.SetY1(0.72)
            plot.legend.SetX2(0.98)
            plot.legend.SetY2(0.88)
            if len(pdf_entries) > 5:
                plot.legend.SetNColumns(2)
            ROOT.gStyle.SetPalette(55)
            plot.plot("NOSTACK E1 PLC PMC")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            plot.save("%s/unfolded_unnormalised_%s_pdf_model_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
            ROOT.gStyle.SetPalette(ROOT.kViridis)


    def plot_uncertainty_shifts_unnormalised(self):
        """Do plots of fractional uncertainty shifts on *unnormalised* unfolded distribution"""
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            entries = []
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_rsp_err', ibin, binning_scheme='generator')

            for syst_dict in self.region['experimental_systematics']:
                # For each systematic, get the normalised shifted distribution for this bin
                # Then calculate the shift wrt nominal result, and hence fraction,
                # then save and plot that
                syst_label_no_spaces = cu.no_space_str(syst_dict['label'])
                self.hist_bin_chopper.add_obj('syst_shifted_%s_unfolded' % syst_label_no_spaces, self.unfolder.systs_shifted[syst_dict['label']])
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('syst_shifted_%s_unfolded' % (syst_label_no_spaces), ibin, binning_scheme='generator')
                syst_unfolded_fraction = syst_unfolded_hist_bin.Clone()
                syst_unfolded_fraction.Add(unfolded_hist_bin_total_errors, -1)
                syst_unfolded_fraction.Divide(unfolded_hist_bin_total_errors)
                # Set to abs values so can plot log
                for i in range(1, syst_unfolded_fraction.GetNbinsX()+1):
                    syst_unfolded_fraction.SetBinContent(i, abs(syst_unfolded_fraction.GetBinContent(i)))
                c = Contribution(syst_unfolded_fraction,
                                 label=syst_dict['label'],
                                 line_color=syst_dict['colour'],
                                 line_style=syst_dict.get('linestyle', 1),
                                 line_width=self.line_width,
                                 marker_size=0,
                                 marker_color=syst_dict['colour'])
                entries.append(c)

            # Now create total & stat. unc.or hists
            h_stat = unfolded_hist_bin_stat_errors.Clone()
            h_syst = unfolded_hist_bin_rsp_errors.Clone()
            h_total = unfolded_hist_bin_total_errors.Clone()
            for i in range(1, h_stat.GetNbinsX()+1):
                if unfolded_hist_bin_total_errors.GetBinContent(i) > 0:
                    h_stat.SetBinContent(i, unfolded_hist_bin_stat_errors.GetBinError(i) / unfolded_hist_bin_total_errors.GetBinContent(i))
                    h_syst.SetBinContent(i, unfolded_hist_bin_rsp_errors.GetBinError(i) / unfolded_hist_bin_total_errors.GetBinContent(i))
                    h_total.SetBinContent(i, unfolded_hist_bin_total_errors.GetBinError(i) / unfolded_hist_bin_total_errors.GetBinContent(i))
                else:
                    h_stat.SetBinContent(i, 0)
                    h_syst.SetBinContent(i, 0)
                    h_total.SetBinContent(i, 0)
                h_stat.SetBinError(i, 0)
                h_syst.SetBinError(i, 0)
                h_total.SetBinError(i, 0)
            c_stat = Contribution(h_stat,
                                 label="Input stats",
                                 line_color=ROOT.kRed,
                                 line_style=3,
                                 line_width=3,
                                 marker_size=0,
                                 marker_color=ROOT.kRed,
                                 )
            c_syst = Contribution(h_syst,
                                 label="Response matrix stats",
                                 line_color=ROOT.kGray+2,
                                 line_style=3,
                                 line_width=3,
                                 marker_size=0,
                                 marker_color=ROOT.kGray+2,
                                 )
            c_tot = Contribution(h_total,
                                 label="Total",
                                 line_color=ROOT.kBlack,
                                 line_style=1,
                                 line_width=3,
                                 marker_size=0,
                                 marker_color=ROOT.kBlack,
                                 )
            entries.extend([c_stat, c_syst, c_tot])
            if not self.check_entries(entries, "plot_uncertainty_shifts_unnormalised %d" % (ibin)):
                return
            plot = Plot(entries,
                        what="hist",
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        xtitle=self.setup.pt_str,
                        ytitle="| Fractional shift on unnormalised distribution |",
                        has_data=self.setup.has_data)
            plot.legend.SetX1(0.55)
            plot.legend.SetY1(0.68)
            plot.legend.SetX2(0.98)
            plot.legend.SetY2(0.88)
            plot.legend.SetNColumns(2)
            plot.left_margin = 0.16
            plot.y_padding_max_linear = 1.4
            plot.plot("NOSTACK HIST")
            plot.set_logx(do_more_labels=False)
            output_filename = "%s/unfolded_unnormalised_systs_%s_lambda_bin_%d.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt)
            plot.save(output_filename)

            plot.y_padding_max_log = 50
            plot.set_logy(do_more_labels=False)
            plot.get_modifier().SetMinimum(1E-4)
            log_filename, ext = os.path.splitext(output_filename)
            plot.save(log_filename+"_log"+ext)

    def plot_unfolded_with_exp_systs_unnormalised(self):
        """Plot shifted unfolded normalised distributions for each syst"""
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')

            entries = []
            for syst_dict in self.region['experimental_systematics']:
                syst_label_no_spaces = cu.no_space_str(syst_dict['label'])
                self.hist_bin_chopper.add_obj('syst_shifted_%s_unfolded' % syst_label_no_spaces, self.unfolder.systs_shifted[syst_dict['label']])
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('syst_shifted_%s_unfolded' % (syst_label_no_spaces), ibin, binning_scheme='generator')
                c = Contribution(syst_unfolded_hist_bin,
                                 label=syst_dict['label'],
                                 line_color=syst_dict['colour'], line_width=self.line_width,
                                 line_style=2 if 'down' in syst_dict['label'].lower() else 1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=unfolded_hist_bin_total_errors)
                entries.append(c)

            entries.append(
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75),
            )
            if not self.check_entries(entries, "plot_unfolded_with_exp_systs_unnormalised %d" % ibin):
                return
            plot = Plot(entries,
                        xtitle=self.setup.pt_str,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        what="hist",
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='Syst / nominal',
                        subplot_limits=(0.75, 1.25))
            self._modify_plot(plot)
            plot.legend.SetX1(0.55)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.98)
            plot.legend.SetY2(0.9)
            if len(entries) > 5: plot.legend.SetNColumns(2)
            # plot.subplot_limits = (0.9, 1.1)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            plot.save("%s/unfolded_unnormalised_syst_variations_%s_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_detector_unnormalised(self, alt_detector=None):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            input_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('input_hist_gen_binning_bg_subtracted', ibin, binning_scheme='generator')
            mc_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_mc_reco_gen_binning_bg_subtracted', ibin, binning_scheme='generator')

            entries = [
                Contribution(mc_hist_bin,
                             label="MC (bg-subtracted) [%s]" % (self.setup.region['mc_label']) if alt_detector else "MC (bg-subtracted)",
                             line_color=self.plot_colours['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_mc_colour'], marker_size=0,
                             subplot=input_hist_bin if alt_detector else None),
            ]
            if alt_detector:
                self.hist_bin_chopper.add_obj("alt_hist_mc_reco_gen_binning_bg_subtracted", alt_detector)
                alt_mc_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('alt_hist_mc_reco_gen_binning_bg_subtracted', ibin, binning_scheme='generator')
                entries.append(Contribution(alt_mc_hist_bin,
                                            label="MC (bg-subtracted) [%s]" % self.setup.region['alt_mc_label'],
                                            line_color=self.plot_colours['alt_reco_colour'], line_width=self.line_width,
                                            marker_color=self.plot_colours['alt_reco_colour'], marker_size=0,
                                            subplot=input_hist_bin))
            entries.append(Contribution(input_hist_bin,
                             label="Data (bg-subtracted)",
                             line_color=self.plot_colours['reco_data_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_data_colour'], marker_style=20, marker_size=0.75,
                             subplot=None if alt_detector else mc_hist_bin)
            )
            if not self.check_entries(entries, "plot_detector_unnormalised %d" % (ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.pt_str,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        what="hist",
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='MC / Data' if alt_detector else 'Data / MC',
                        subplot_limits=(0.75, 1.25) if alt_detector else (0.75, 1.25)
                       )
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            plot.save("%s/detector_unnormalised_gen_binning_bg_subtracted_%s_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))


class RecoPtBinnedPlotter(object):
    def __init__(self, setup, bins, hist_bin_chopper, unfolder):
        self.setup = setup
        self.region = setup.region
        self.bins = bins
        self.hist_bin_chopper = hist_bin_chopper

        self.line_width = 2
        self.plot_colours = PLOT_COLOURS
        self.pt_bin_plot_args = dict(
            what="hist",
            xtitle=self.setup.detector_title,
            has_data=self.setup.has_data,
            subplot_type='ratio',
            subplot_title="Unfolded / Gen",
            subplot_limits=(0, 2) if self.setup.has_data else (0.75, 1.25),
        )
        self.unfolder = unfolder

    @staticmethod
    def _modify_plot(this_plot):
        this_plot.legend.SetX1(0.6)
        this_plot.legend.SetY1(0.68)
        this_plot.legend.SetX2(0.98)
        this_plot.legend.SetY2(0.9)
        this_plot.left_margin = 0.16

    @staticmethod
    def check_entries(entries, message=""):
        """Check that at least 1 Contribution has something in it"""
        has_entries = [c.obj.GetEntries() > 0 for c in entries]
        if not any(has_entries):
            if message:
                print("Skipping 0 entries (%s)" % (message))
            else:
                print("Skipping 0 entries")
            return False

        max_bin = max([c.obj.GetMaximum() for c in entries])
        min_bin = min([c.obj.GetMinimum() for c in entries])
        if max_bin == min_bin:
            return False

        return True

    def get_pt_bin_title(self, bin_edge_low, bin_edge_high):
        title = (("{jet_algo}\n"
                  "{region_label} region\n"
                  "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                 .format(
                    jet_algo=self.setup.jet_algo,
                    region_label=self.region['label'],
                    pt_str=self.setup.pt_var_str,
                    bin_edge_low=bin_edge_low,
                    bin_edge_high=bin_edge_high
                ))
        return title

    def plot_detector_normalised(self, alt_detector=None):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_bg_subtracted', ibin, binning_scheme='detector')
            mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')

            entries = [
                Contribution(mc_hist_bin,
                             label="MC (bg-subtracted) [%s]" % (self.setup.region['mc_label']) if alt_detector else "MC (bg-subtracted)",
                             line_color=self.plot_colours['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_mc_colour'], marker_size=0,
                             subplot=input_hist_bin if alt_detector else None)
            ]
            if alt_detector:
                self.hist_bin_chopper.add_obj("alt_hist_mc_reco_bg_subtracted", alt_detector)
                alt_mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')
                entries.append(Contribution(alt_mc_hist_bin,
                                            label="MC (bg-subtracted) [%s]" % self.setup.region['alt_mc_label'],
                                            line_color=self.plot_colours['alt_reco_colour'], line_width=self.line_width,
                                            marker_color=self.plot_colours['alt_reco_colour'], marker_size=0,
                                            subplot=input_hist_bin if alt_detector else None)
                )
            entries.append(
                Contribution(input_hist_bin,
                             label="Data (bg-subtracted)",
                             line_color=self.plot_colours['reco_data_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_data_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_hist_bin if not alt_detector else None)
            )
            if not self.check_entries(entries, "plot_detector_normalised %d" % (ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='MC / Data' if alt_detector else 'Data / MC',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/detector_reco_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_folded_unfolded_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_bg_subtracted', ibin, binning_scheme='detector')
            folded_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_unfolded', ibin, binning_scheme='detector')

            entries = [
                Contribution(input_hist_bin,
                             label="Unfolding input (bg-subtracted)",
                             line_color=self.plot_colours['reco_unfolding_input_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_unfolding_input_colour'], marker_style=20, marker_size=0),
                Contribution(folded_unfolded_hist_bin,
                             label="Folded unfolded",
                             line_color=self.plot_colours['reco_folded_unfolded_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_folded_unfolded_colour'], marker_style=20, marker_size=0.75,
                             subplot=input_hist_bin),
            ]
            if not self.check_entries(entries, "plot_folded_unfolded_normalised %d" % (ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='Folded / input',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/detector_folded_unfolded_only_data_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_folded_unfolded_with_mc_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_bg_subtracted', ibin, binning_scheme='detector')
            folded_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_unfolded', ibin, binning_scheme='detector')

            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_colours['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_mc_colour'], marker_style=20, marker_size=0),
                Contribution(input_hist_bin,
                             label="Unfolding input (bg-subtracted)",
                             line_color=self.plot_colours['reco_unfolding_input_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_unfolding_input_colour'], marker_style=20, marker_size=0,
                             subplot=mc_reco_hist_bin),
                Contribution(folded_unfolded_hist_bin,
                             label="Folded unfolded",
                             line_color=self.plot_colours['reco_folded_unfolded_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_folded_unfolded_colour'],marker_style=20, marker_size=0.75,
                             subplot=mc_reco_hist_bin),
            ]
            if not self.check_entries(entries, "plot_folded_unfolded_normalised %d" % (ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='* / MC',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/detector_folded_unfolded_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_folded_gen_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')
            folded_mc_truth_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_mc_truth', ibin, binning_scheme='detector')

            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_colours['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_mc_colour'], marker_style=20, marker_size=0),
                Contribution(folded_mc_truth_hist_bin,
                             label="Folded gen",
                             line_color=self.plot_colours['reco_folded_mc_truth_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_folded_mc_truth_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_reco_hist_bin),
            ]
            if not self.check_entries(entries, "plot_folded_gen_normalised %d" % (ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='Folded / MC reco',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/detector_folded_gen_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))


def do_all_plots_per_region_angle(setup, unpack_dict):
    # Note that experimental systs are only different response matrices, and are stored in the main unfolder
    has_exp_systs = len(setup.region['experimental_systematics']) > 0
    has_model_systs = len(setup.region['model_systematics']) > 0
    has_pdf_systs = len(setup.region['pdf_systematics']) > 0

    if has_exp_systs: print("We have experimental systs")
    if has_model_systs: print("We have model systs")
    if has_pdf_systs: print("We have pdf systs")

    unfolder = unpack_dict['unfolder']
    unreg_unfolder = unpack_dict['unreg_unfolder']
    alt_unfolder = unpack_dict['alt_unfolder']
    alt_hist_truth = unpack_dict['alt_hist_truth']
    alt_hist_reco = unpack_dict['alt_hist_reco']
    alt_hist_reco_bg_subtracted = unpack_dict['alt_hist_reco_bg_subtracted']
    alt_hist_reco_bg_subtracted_gen_binning = unpack_dict['alt_hist_reco_bg_subtracted_gen_binning']

    # Big 1D plots to compare things
    hbc = HistBinChopper(generator_binning=unfolder.generator_binning.FindNode("generatordistribution"),
                         detector_binning=unfolder.detector_binning.FindNode("detectordistribution"))
    hbc.add_obj("hist_truth", unfolder.hist_truth)
    hbc.update(unfolder.hist_bin_chopper)   # update the HistBinChopper with the new normalised systematics already produced in unfolder

    # Iterate through pt bins - gen binning
    # ------------------------------------------------------------------
    gen_pt_binned_plotter = GenPtBinnedPlotter(setup=setup,
                                               bins=unfolder.pt_bin_edges_gen,
                                               hist_bin_chopper=hbc,
                                               unfolder=unfolder)
    gen_pt_binned_plotter.plot_unfolded_unnormalised()
    gen_pt_binned_plotter.plot_unfolded_normalised()
    if alt_hist_truth:
        gen_pt_binned_plotter.hist_bin_chopper.add_obj('alt_hist_truth', alt_hist_truth)
        gen_pt_binned_plotter.plot_unfolded_with_alt_truth_normalised()

    if unfolder.tau > 0 and unreg_unfolder:
        gen_pt_binned_plotter.hist_bin_chopper.add_obj("unreg_unfolded", unreg_unfolder.unfolded)
        gen_pt_binned_plotter.plot_unfolded_with_unreg_normalised()

    if alt_unfolder:
        gen_pt_binned_plotter.hist_bin_chopper.add_obj("alt_unfolded_stat_err", alt_unfolder.unfolded_stat_err)
        gen_pt_binned_plotter.plot_unfolded_with_alt_response_normalised(alt_unfolder=alt_unfolder)
        gen_pt_binned_plotter.plot_unfolded_with_alt_response_truth_normalised(alt_unfolder=alt_unfolder)

    if has_exp_systs:
        gen_pt_binned_plotter.plot_uncertainty_shifts_normalised()
        gen_pt_binned_plotter.plot_unfolded_with_exp_systs_normalised()
        gen_pt_binned_plotter.plot_exp_syst_variation_normalised()

    if has_model_systs:
        gen_pt_binned_plotter.plot_unfolded_with_model_systs_normalised()

    if has_pdf_systs:
        gen_pt_binned_plotter.plot_unfolded_with_pdf_systs_normalised()

    # if has_data:
    gen_pt_binned_plotter.hist_bin_chopper.add_obj("input_hist_gen_binning_bg_subtracted", unfolder.input_hist_gen_binning_bg_subtracted)
    gen_pt_binned_plotter.hist_bin_chopper.add_obj("hist_mc_reco_gen_binning_bg_subtracted", unfolder.hist_mc_reco_gen_binning_bg_subtracted)
    gen_pt_binned_plotter.plot_detector_normalised(alt_detector=alt_hist_reco_bg_subtracted_gen_binning)

    # Iterate through lambda bins - gen binning
    # ------------------------------------------------------------------
    lambda_pt_binned_plotter = GenLambdaBinnedPlotter(setup=setup,
                                                      bins=unfolder.variable_bin_edges_gen,
                                                      hist_bin_chopper=hbc,  # this is the same object as gen_pt_binned_plotter, so has all the objects already
                                                      unfolder=unfolder)

    lambda_pt_binned_plotter.plot_unfolded_unnormalised()

    if unfolder.tau > 0 and unreg_unfolder:
        lambda_pt_binned_plotter.plot_unfolded_with_unreg_unnormalised(unreg_unfolder)

    if alt_unfolder:
        lambda_pt_binned_plotter.plot_unfolded_with_alt_response_unnormalised(alt_unfolder=alt_unfolder)

    if has_exp_systs:
        lambda_pt_binned_plotter.plot_uncertainty_shifts_unnormalised()
        lambda_pt_binned_plotter.plot_unfolded_with_exp_systs_unnormalised()

    if has_model_systs:
        lambda_pt_binned_plotter.plot_unfolded_with_model_systs_unnormalised()

    if has_pdf_systs:
        lambda_pt_binned_plotter.plot_unfolded_with_pdf_systs_unnormalised()

    # if has_data:
    lambda_pt_binned_plotter.plot_detector_unnormalised(alt_detector=alt_hist_reco_bg_subtracted_gen_binning)

    # Iterate through pt bins - reco binning
    # ------------------------------------------------------------------
    hbc.add_obj("hist_mc_reco_bg_subtracted", unfolder.hist_mc_reco_bg_subtracted)
    hbc.add_obj("input_hist_bg_subtracted", unfolder.input_hist_bg_subtracted)
    hbc.add_obj("folded_unfolded", unfolder.folded_unfolded)
    hbc.add_obj("folded_mc_truth", unfolder.folded_mc_truth)
    reco_pt_binned_plotter = RecoPtBinnedPlotter(setup=setup,
                                                 bins=unfolder.pt_bin_edges_reco,
                                                 hist_bin_chopper=hbc,
                                                 unfolder=unfolder)
    reco_pt_binned_plotter.plot_detector_normalised(alt_detector=alt_hist_reco_bg_subtracted)
    reco_pt_binned_plotter.plot_folded_unfolded_normalised()
    reco_pt_binned_plotter.plot_folded_unfolded_with_mc_normalised()
    reco_pt_binned_plotter.plot_folded_gen_normalised()

    return hbc


class BigNormalised1DPlotter(object):

    def __init__(self, setup, hist_bin_chopper, unfolder):
        self.setup = setup
        self.hist_bin_chopper = hist_bin_chopper
        self.unfolder = unfolder
        self.pt_bin_edges = unfolder.pt_bin_edges_gen
        self.lambda_bin_edges = unfolder.variable_bin_edges_gen
        self.num_pt_bins = len(self.pt_bin_edges)-1
        self.num_lambda_bins = len(self.lambda_bin_edges)-1
        self.nbins = self.num_pt_bins * self.num_lambda_bins
        self.all_bin_edges = None # gets filled the first time make_big_1d_normalised_hist() is called
        self.plot_with_bin_widths = False # use bin widths in plot instead of uniform bin widths
        self.line_width = 1

        self._cache_1d = {}

    @staticmethod
    def _get_ylim(entries):
        ymax = max([e.obj.GetMaximum() for e in entries])
        ylim = (0, 3*ymax)  # room for labels
        return ylim

    @staticmethod
    def _modify_plot(this_plot):
        this_plot.legend.SetX1(0.6)
        this_plot.legend.SetY1(0.72)
        this_plot.legend.SetX2(0.9)
        this_plot.legend.SetY2(0.88)
        if len(this_plot.contributions) > 4:
            this_plot.legend.SetNColumns(2)
        # this_plot.left_margin = 0.16
        this_plot.left_margin = 0.1
        this_plot.default_canvas_size = (900, 700)
        this_plot.y_padding_max_linear = 2
        this_plot.subplot_pad_height = 0.4
        this_plot.subplot_line_style = 3
        this_plot.subplot_line_width = 1
        this_plot.subplot_line_color = ROOT.kGray+1


    def _plot_pt_bins(self, plot):
        sub_xax = plot.subplot_container.GetXaxis()
        sub_xax.SetTitleOffset(sub_xax.GetTitleOffset()*.9)
        sub_xax.SetTitleOffset(sub_xax.GetTitleOffset()*.9)
        sub_xax.SetLabelOffset(999)
        ndiv = self.num_pt_bins
        plot.container.GetXaxis().SetNdivisions(ndiv, False)
        sub_xax.SetNdivisions(ndiv, False)

        lines = []
        texts = []
        plot.main_pad.cd()
        obj = plot.container.GetHistogram()
        y_low, y_high = obj.GetMinimum(), obj.GetMaximum()
        for ibin, _ in enumerate(self.pt_bin_edges[1:-1], 1):
            pt_bin = self.all_bin_edges[ibin * self.num_lambda_bins]
            plot.main_pad.cd()
            line = ROOT.TLine(pt_bin, y_low, pt_bin, y_high)
            line.SetLineStyle(2)
            line.SetLineColor(14)
            line.Draw()
            lines.append(line)

            # lines on subplot if available
            if isinstance(plot, Plot) and plot.subplot_pad:
                plot.subplot_pad.cd()
                subplot_y_low, subplot_y_high = plot.subplot_container.GetHistogram().GetMinimum(), plot.subplot_container.GetHistogram().GetMaximum()  # BINGO
                line = ROOT.TLine(pt_bin, subplot_y_low, pt_bin, subplot_y_high)
                line.SetLineStyle(2)
                line.SetLineColor(14)
                line.Draw()
                lines.append(line)

        plot.main_pad.cd()

        # add bin texts
        for ibin, (pt_val, pt_val_upper) in enumerate(zip(self.pt_bin_edges[:-1], self.pt_bin_edges[1:])):
            labels_inside_align = 'lower'

            # fixgure out x location
            pt_bin = self.all_bin_edges[ibin * self.num_lambda_bins]
            pt_bin_higher = self.all_bin_edges[(ibin+1) * self.num_lambda_bins]
            pt_bin_interval = pt_bin_higher - pt_bin
            text_x = pt_bin + 0.35*(pt_bin_interval)

            # figure out y location from axes
            axis_range = y_high - y_low
            # assumes logarithmic?
            if ROOT.gPad.GetLogy():
                if y_low <= 0:
                    print("y_low is %f so can't take log" %y_low)
                    return None, None
                    # raise ValueError("y_low is %f so can't take log" %y_low)
                log_axis_range = math.log10(y_high) - math.log10(y_low)
                y_start = math.pow(10, 0.03*log_axis_range + math.log10(y_low))
                y_end = 10*y_low
            else:
                y_start = (0.02*axis_range) + y_low
                y_start = (0.3*axis_range) + y_low
                y_end = 10*y_low

            if labels_inside_align == 'higher':
                if ROOT.gPad.GetLogy():
                    y_start = y_high/10
                    y_end = y_high/1.05
                else:
                    y_start = y_high*0.8
                    y_end = y_high*0.9

            text = ROOT.TPaveText(text_x, y_start, text_x + .5*pt_bin_interval, y_start)
            text.SetBorderSize(0)
            text.SetFillStyle(0)
            contents = '%g < p_{T} < %g GeV' % (pt_val, pt_val_upper)
            # You have to style the thing that is returned by AddText, not the TPaveText obj itself
            # WTAF
            t = text.AddText(contents)  # t is a TText
            t.SetTextColor(14)
            t.SetTextAngle(89)
            t.SetTextSize(0.0275)
            if (isinstance(plot, Plot) and not plot.subplot_pad) or not isinstance(plot, Plot):
                t.SetTextSize(0.02)  # account for the fact that a subplot makes things smaller
            if labels_inside_align == 'lower':
                t.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignTop)
            else:
                t.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignTop)
            text.Draw()
            texts.append(text)

        return lines, texts

    def make_big_1d_normalised_hist(self, name):
        """Make big 1D plot with normalised distribution per pt bin"""
        if self.all_bin_edges is None:
            if self.plot_with_bin_widths:
                bins = []
                for i in range(self.num_pt_bins):
                    offset = i * math.ceil(self.lambda_bin_edges[-1])
                    bins.extend(self.lambda_bin_edges[:-1] + offset)
                    if i == (self.num_pt_bins-1):
                        bins.append(self.lambda_bin_edges[-1:] + offset)
                self.all_bin_edges = array('d', bins)
            else:
                self.all_bin_edges = array('d', list(range(0, self.nbins+1)))

        h_new = ROOT.TH1D(name + "_1d_all_" + cu.get_unique_str(), "", self.nbins, self.all_bin_edges)
        for ibin, _ in enumerate(self.pt_bin_edges[:-1]):
            hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(name, ibin, binning_scheme='generator')

            for lbin in range(1, hist_bin.GetNbinsX()+1):
                global_bin = (ibin * self.num_lambda_bins) + lbin
                h_new.SetBinContent(global_bin, hist_bin.GetBinContent(lbin))
                h_new.SetBinError(global_bin, hist_bin.GetBinError(lbin))

        return h_new

    def get_mc_truth_kwargs(self):
        return dict(
            label="Generator (%s)" % (self.setup.region['mc_label']),
            line_color=PLOT_COLOURS['gen_colour'], line_width=self.line_width,
            marker_color=PLOT_COLOURS['gen_colour'], marker_size=0
            )

    def get_unfolded_total_err_kwargs(self):
        return dict(
            label="Unfolded (#tau = %.3g) (total unc.)\n(%s response matrix)" % (self.unfolder.tau, self.setup.region['mc_label']),
            line_color=PLOT_COLOURS['unfolded_total_colour'], line_width=self.line_width, line_style=1,
            marker_color=PLOT_COLOURS['unfolded_total_colour']
        )

    def get_unfolded_stat_err_kwargs(self):
        return dict(
            label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (self.unfolder.tau, self.setup.region['mc_label']),
            line_color=PLOT_COLOURS['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
            marker_color=PLOT_COLOURS['unfolded_stat_colour'], marker_style=20, marker_size=0
        )

    def get_alt_mc_truth_kwargs(self):
        return dict(
            label="Generator (%s)" % (self.setup.region['alt_mc_label']),
            line_color=PLOT_COLOURS['alt_gen_colour'], line_width=self.line_width, line_style=2,
            marker_color=PLOT_COLOURS['alt_gen_colour'], marker_size=0
        )

    def get_alt_unfolded_stat_err_kwargs(self, alt_unfolder):
        return dict(
            label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (alt_unfolder.tau, self.setup.region['alt_mc_label']),
            line_color=PLOT_COLOURS['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
            marker_color=PLOT_COLOURS['alt_unfolded_colour']
        )

    def get_title(self):
        return (("{jet_algo}\n"
                 "{region_label} region")
                     .format(
                        jet_algo=self.setup.jet_algo,
                        region_label=self.setup.region['label'],
                        pt_str=self.setup.pt_var_str)
               )

    def get_subplot_ylim(self):
        return (0, 2) if self.setup.has_data else (0.7, 1.3)

    def get_big_1d(self, name):
        if self._cache_1d.get(name, None) is None:
            self._cache_1d[name] = self.make_big_1d_normalised_hist(name)
        return self._cache_1d[name]

    def plot_unfolded_truth(self):
        entries = [
            Contribution(self.get_big_1d('hist_truth'),
                         **self.get_mc_truth_kwargs()),
            Contribution(self.get_big_1d('unfolded'),
                         subplot=self.get_big_1d('hist_truth'),
                         **self.get_unfolded_total_err_kwargs()),
            Contribution(self.get_big_1d('unfolded_stat_err'),
                         subplot=self.get_big_1d('hist_truth'),
                         **self.get_unfolded_stat_err_kwargs()),
        ]
        plot = Plot(entries, 'hist',
                    ytitle=self.setup.pt_bin_normalised_differential_label,
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(entries),
                    subplot_type='ratio',
                    subplot_title="Unfolded / Gen",
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_truth_alt_truth(self):
        entries = [
            Contribution(self.get_big_1d('hist_truth'),
                         **self.get_mc_truth_kwargs()),
            Contribution(self.get_big_1d('alt_hist_truth'),
                         subplot=self.get_big_1d('hist_truth'),
                         **self.get_alt_mc_truth_kwargs()),
            Contribution(self.get_big_1d('unfolded'),
                         subplot=self.get_big_1d('hist_truth'),
                         **self.get_unfolded_total_err_kwargs()),
        ]
        plot = Plot(entries, 'hist',
                    ytitle=self.setup.pt_bin_normalised_differential_label,
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(entries),
                    subplot_type='ratio',
                    subplot_title="#splitline{* / Gen}{(%s)}" % (self.setup.region['mc_label']),
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_alt_truth_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_truth_alt_response(self, alt_unfolder):
        entries = [
            Contribution(self.get_big_1d('hist_truth'),
                         **self.get_mc_truth_kwargs()),
            Contribution(self.get_big_1d('unfolded_stat_err'),
                         subplot=self.get_big_1d('hist_truth'),
                         **self.get_unfolded_stat_err_kwargs()),
            Contribution(self.get_big_1d('alt_unfolded_stat_err'),
                         subplot=self.get_big_1d('hist_truth'),
                         **self.get_alt_unfolded_stat_err_kwargs(alt_unfolder)),
        ]
        plot = Plot(entries, 'hist',
                    ytitle=self.setup.pt_bin_normalised_differential_label,
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(entries),
                    subplot_type='ratio',
                    subplot_title="Unfolded / Gen",
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_alt_response_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_truth_alt_response_truth(self, alt_unfolder):
        entries = [
            Contribution(self.get_big_1d('hist_truth'),
                         **self.get_mc_truth_kwargs()),
            Contribution(self.get_big_1d('alt_hist_truth'),
                         subplot=self.get_big_1d('hist_truth'),
                         **self.get_alt_mc_truth_kwargs()),
            Contribution(self.get_big_1d('unfolded'),
                         subplot=self.get_big_1d('hist_truth'),
                         **self.get_unfolded_total_err_kwargs()),
            Contribution(self.get_big_1d('alt_unfolded_stat_err'),
                         subplot=self.get_big_1d('hist_truth'),
                         **self.get_alt_unfolded_stat_err_kwargs(alt_unfolder)),
        ]
        plot = Plot(entries, 'hist',
                    ytitle=self.setup.pt_bin_normalised_differential_label,
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(entries),
                    subplot_type='ratio',
                    subplot_title="#splitline{* / Gen}{(%s)}" % (self.setup.region['mc_label']),
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_alt_response_truth_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_exp_systs(self):
        for syst_dict in self.setup.region['experimental_systematics']:
            entries = [
                        # Contribution(self.get_big_1d('hist_truth'),
                        #  **self.get_mc_truth_kwargs()),
                        Contribution(self.get_big_1d('unfolded_stat_err'),
                                     **dict(self.get_unfolded_stat_err_kwargs(),
                                            line_color=ROOT.kGray+2))
            ]
            syst_label_no_spaces = cu.no_space_str(syst_dict['label'])
            hbc_name = 'syst_shifted_%s_unfolded' % syst_label_no_spaces
            self.hist_bin_chopper.add_obj(hbc_name, self.unfolder.systs_shifted[syst_dict['label']])
            entries.append(Contribution(self.get_big_1d(hbc_name),
                                        label=syst_dict['label'],
                                        line_color=syst_dict['colour'], line_width=self.line_width,
                                        line_style=2 if 'down' in syst_dict['label'].lower() else 1,
                                        marker_color=syst_dict['colour'], marker_size=0,
                                        subplot=self.get_big_1d('unfolded_stat_err')))
            plot = Plot(entries, 'hist',
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_title(),
                        xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                        has_data=self.setup.has_data,
                        ylim=self._get_ylim(entries),
                        subplot_type='ratio',
                        subplot_title="#splitline{Systematic}{/ Nominal}",
                        subplot_limits=self.get_subplot_ylim()
                        )
            self._modify_plot(plot)
            plot.plot("NOSTACK E")
            l, t = self._plot_pt_bins(plot)
            plot.save("%s/unfolded_1d_normalised_exp_syst_%s_%s_divBinWidth.%s" % (self.setup.output_dir, syst_label_no_spaces, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_model_systs(self):
        all_entries = [Contribution(self.get_big_1d('hist_truth'),
                                     **self.get_mc_truth_kwargs()),
                       Contribution(self.get_big_1d('unfolded_stat_err'),
                                    subplot=self.get_big_1d('hist_truth'),
                                    **dict(self.get_unfolded_stat_err_kwargs(),
                                           line_color=ROOT.kGray+2))
                      ]

        for syst_dict in self.setup.region['model_systematics']:
            entries = [
                        Contribution(self.get_big_1d('hist_truth'),
                                     **self.get_mc_truth_kwargs()),
                        Contribution(self.get_big_1d('unfolded_stat_err'),
                                     subplot=self.get_big_1d('hist_truth'),
                                     **dict(self.get_unfolded_stat_err_kwargs(),
                                            line_color=ROOT.kGray+2))
            ]
            syst_unfolder = syst_dict['unfolder']
            syst_label = syst_dict['label']
            syst_label_no_spaces = cu.no_space_str(syst_dict['label'])

            hbc_name_gen = 'model_syst_%s_hist_truth' % (syst_label_no_spaces)
            self.hist_bin_chopper.add_obj(hbc_name_gen, syst_unfolder.hist_truth)

            hbc_name = 'model_syst_%s_unfolded' % (syst_label_no_spaces)
            self.hist_bin_chopper.add_obj(hbc_name, syst_unfolder.unfolded)

            entries.extend([
                Contribution(self.get_big_1d(hbc_name_gen),
                             label="Generator (%s)" % (syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                             marker_color=syst_dict['colour'], marker_size=0),
                Contribution(self.get_big_1d(hbc_name),
                             label="Unfolded (#tau = %.3g) (stat. err) (%s)" % (syst_unfolder.tau, syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                             marker_color=syst_dict['colour'], marker_size=0,
                             subplot=self.get_big_1d(hbc_name_gen)),
            ])
            all_entries.extend(entries[-2:])
            plot = Plot(entries, 'hist',
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_title(),
                        xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                        has_data=self.setup.has_data,
                        ylim=self._get_ylim(entries),
                        subplot_type='ratio',
                        subplot_title="Unfolded / Gen",
                        subplot_limits=self.get_subplot_ylim()
                        )
            self._modify_plot(plot)
            plot.plot("NOSTACK E")
            l, t = self._plot_pt_bins(plot)
            plot.save("%s/unfolded_1d_normalised_model_syst_%s_%s_divBinWidth.%s" % (self.setup.output_dir, syst_label_no_spaces, self.setup.append, self.setup.output_fmt))

        plot = Plot(all_entries, 'hist',
                    ytitle=self.setup.pt_bin_normalised_differential_label,
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(all_entries),
                    subplot_type='ratio',
                    subplot_title="Unfolded / Gen",
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_all_model_syst_%s_%s_divBinWidth.%s" % (self.setup.output_dir, syst_label_no_spaces, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_model_systs_only_scale(self):
        all_entries = [Contribution(self.get_big_1d('hist_truth'),
                                     **self.get_mc_truth_kwargs()),
                       Contribution(self.get_big_1d('unfolded_stat_err'),
                                    subplot=self.get_big_1d('hist_truth'),
                                    **dict(self.get_unfolded_stat_err_kwargs(),
                                           line_color=ROOT.kGray+2))
                      ]

        for syst_dict in self.setup.region['model_systematics']:
            syst_unfolder = syst_dict['unfolder']
            syst_label = syst_dict['label']
            if 'Herwig' in syst_label or 'Pythia' in syst_label:
                continue

            syst_label_no_spaces = cu.no_space_str(syst_dict['label'])

            hbc_name_gen = 'model_syst_%s_hist_truth' % (syst_label_no_spaces)
            self.hist_bin_chopper.add_obj(hbc_name_gen, syst_unfolder.hist_truth)

            hbc_name = 'model_syst_%s_unfolded' % (syst_label_no_spaces)
            self.hist_bin_chopper.add_obj(hbc_name, syst_unfolder.unfolded)

            all_entries.extend([
                Contribution(self.get_big_1d(hbc_name_gen),
                             label="Generator (%s)" % (syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                             marker_color=syst_dict['colour'], marker_size=0),
                Contribution(self.get_big_1d(hbc_name),
                             label="Unfolded (#tau = %.3g) (stat. err) (%s)" % (syst_unfolder.tau, syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                             marker_color=syst_dict['colour'], marker_size=0,
                             subplot=self.get_big_1d(hbc_name_gen)),
            ])

        plot = Plot(all_entries, 'hist',
                    ytitle=self.setup.pt_bin_normalised_differential_label,
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(all_entries),
                    subplot_type='ratio',
                    subplot_title="Unfolded / Gen",
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_scale_model_syst_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_pdf_systs(self):
        entries = [
                    Contribution(self.get_big_1d('hist_truth'),
                                 **self.get_mc_truth_kwargs()),
                    Contribution(self.get_big_1d('unfolded_stat_err'),
                                 subplot=self.get_big_1d('hist_truth'),
                                 **dict(self.get_unfolded_stat_err_kwargs(),
                                        line_color=ROOT.kGray+2))
        ]
        for syst_dict in self.setup.region['pdf_systematics']:
            syst_unfolder = syst_dict['unfolder']
            syst_label = syst_dict['label']
            syst_label_no_spaces = cu.no_space_str(syst_dict['label'])

            hbc_name_gen = 'model_syst_%s_hist_truth' % (syst_label_no_spaces)
            self.hist_bin_chopper.add_obj(hbc_name_gen, syst_unfolder.hist_truth)

            hbc_name = 'model_syst_%s_unfolded' % (syst_label_no_spaces)
            self.hist_bin_chopper.add_obj(hbc_name, syst_unfolder.unfolded)

            entries.extend([
                Contribution(self.get_big_1d(hbc_name_gen),
                             label="Generator (%s)" % (syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                             marker_color=syst_dict['colour'], marker_size=0),
                Contribution(self.get_big_1d(hbc_name),
                             label="Unfolded (#tau = %.3g) (stat. err) (%s)" % (syst_unfolder.tau, syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                             marker_color=syst_dict['colour'], marker_size=0,
                             subplot=self.get_big_1d(hbc_name_gen)),
            ])
        plot = Plot(entries, 'hist',
                    ytitle=self.setup.pt_bin_normalised_differential_label,
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(entries),
                    subplot_type='ratio',
                    subplot_title="Unfolded / Gen",
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.legend.SetNColumns(3)
        plot.legend.SetX1(0.5)
        plot.plot("NOSTACK PLC PMC E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_pdf_syst_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    # def plot_folded_gen(self):
    #     self.hist_bin_chopper.add_obj("hist_mc_reco_bg_subtracted", self.unfolder.hist_mc_reco_bg_subtracted)
    #     mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(, ibin, binning_scheme='detector')

    #     self.hist_bin_chopper.add_obj("folded_mc_truth", self.unfolder.folded_mc_truth)
    #     folded_mc_truth_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_mc_truth', ibin, binning_scheme='detector')

    #     entries = [
    #         Contribution(self.get_big_1d("hist_mc_reco_bg_subtracted"),
    #                      label=""),
    #         Contribution(self.get_big_1d("folded_mc_truth"),
    #                      label="")
    #     ]
    #     plot = Plot(entries, 'hist',
    #                 ytitle=self.setup.pt_bin_normalised_differential_label,
    #                 title=self.get_title(),
    #                 xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
    #                 has_data=self.setup.has_data,
    #                 ylim=self._get_ylim(entries),
    #                 subplot_type='ratio',
    #                 subplot_title="Unfolded / Gen",
    #                 subplot_limits=self.get_subplot_ylim()
    #                 )
    #     self._modify_plot(plot)
    #     plot.plot("NOSTACK E")
    #     l, t = self._plot_pt_bins(plot)
    #     plot.save("%s/folded_gen_1d_normalised_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))


def do_all_big_1d_plots_per_region_angle(setup, unpack_dict, hist_bin_chopper=None):
    unfolder = unpack_dict['unfolder']
    alt_unfolder = unpack_dict['alt_unfolder']
    alt_hist_truth = unpack_dict['alt_hist_truth']

    if not hist_bin_chopper:
        # unreg_unfolder = unpack_dict['unreg_unfolder']
        hist_bin_chopper = HistBinChopper(generator_binning=unfolder.generator_binning.FindNode("generatordistribution"),
                                          detector_binning=unfolder.detector_binning.FindNode("detectordistribution"))
        hist_bin_chopper.add_obj("hist_truth", unfolder.hist_truth)
        hist_bin_chopper.update(unfolder.hist_bin_chopper)
        if alt_unfolder:
            hist_bin_chopper.add_obj("alt_unfolded_stat_err", alt_unfolder.unfolded_stat_err)
        if alt_hist_truth:
            hist_bin_chopper.add_obj("alt_hist_truth", alt_hist_truth)

    has_exp_systs = len(setup.region['experimental_systematics']) > 0
    has_model_systs = len(setup.region['model_systematics']) > 0
    has_pdf_systs = len(setup.region['pdf_systematics']) > 0

    big_plotter = BigNormalised1DPlotter(setup, hist_bin_chopper, unfolder)
    big_plotter.plot_unfolded_truth()

    if alt_hist_truth:
        big_plotter.plot_unfolded_truth_alt_truth()

    if alt_unfolder:
        big_plotter.plot_unfolded_truth_alt_response(alt_unfolder)
        big_plotter.plot_unfolded_truth_alt_response_truth(alt_unfolder)

    if has_exp_systs:
        big_plotter.plot_unfolded_exp_systs()

    if has_model_systs:
        big_plotter.plot_unfolded_model_systs()
        big_plotter.plot_unfolded_model_systs_only_scale()

    if has_pdf_systs:
        big_plotter.plot_unfolded_pdf_systs()


def store_bottom_line_stats(all_chi2_stats, setup, unpack_dict):
    unfolder = unpack_dict['unfolder']

    smeared_chi2, smeared_ndf, smeared_p = unfolder.calculate_chi2(one_hist=unfolder.get_folded_mc_truth(),
                                                                   other_hist=unfolder.hist_mc_reco_bg_subtracted,
                                                                   cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
                                                                   # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_ndarray(),
                                                                   detector_space=True,
                                                                   ignore_underflow_bins=True,
                                                                   debugging_dir=None)
    print('smeared chi2:', smeared_chi2, smeared_ndf, smeared_chi2/smeared_ndf, smeared_p)

    unfolded_chi2, unfolded_ndf, unfolded_p = unfolder.calculate_chi2(one_hist=unfolder.unfolded,
                                                                      other_hist=unfolder.hist_truth,
                                                                      cov_inv_matrix=unfolder.get_vxx_inv_ndarray(),
                                                                      detector_space=False,
                                                                      ignore_underflow_bins=True,
                                                                      debugging_dir=None)
    print('unfolded chi2:', unfolded_chi2, unfolded_ndf, unfolded_chi2/unfolded_ndf, unfolded_p)

    # smeared_alt = unfolder.fold_generator_level(unpack_dict['alt_hist_truth'])
    # smeared_alt_chi2, smeared_alt_ndf, smeared_alt_p =

    all_chi2_stats.append({
            "region": setup.region['label'],
            "is_groomed": "groomed" in setup.region['name'],
            "angle": setup.angle.var,
            # "angle": setup.angle_str,
            "smeared_chi2": smeared_chi2,
            "smeared_ndf": smeared_ndf,
            "smeared_chi2ndf": smeared_chi2/smeared_ndf,
            "smeared_p": smeared_p,
            "unfolded_chi2": unfolded_chi2,
            "unfolded_ndf": unfolded_ndf,
            "unfolded_chi2ndf": unfolded_chi2/unfolded_ndf,
            "unfolded_p": unfolded_p,
    })

def print_chi2_table(df):
    df_sorted = df.sort_values(by=['region', 'angle'])
    print(r'\begin{tabular}{c|c|c|c|c|c|c}')
    print(r"Region & Angle & $\chi_{\text{unfolded}}^2$, $N_{DoF}$, $\chi_{\text{unfolded}}^2 / N_{DoF}$  & $p_{\text{unfolded}}$ & $\chi_{\text{smeared}}^2$, $N_{DoF}$, $\chi_{\text{smeared}}^2 / N_{DoF}$ & $p_{\text{smeared}}$ & $p_{\text{unfolded}} > p_{\text{smeared}}$\\")
    print(r"\hline \\")
    for row in df_sorted.itertuples():
        angle_name = row.angle.replace("jet_", "").replace("_charged", " (charged)")
        if row.is_groomed:
            angle_name = "Groomed " + angle_name
        print("%s & %s & %d / %d = %.3g & %.3g & %d / %d = %.3g & %.3g & %s \\\\" % (row.region, angle_name,
                                                                                     row.unfolded_chi2, row.unfolded_ndf, row.unfolded_chi2ndf, row.unfolded_p,
                                                                                     row.smeared_chi2, row.smeared_ndf, row.smeared_chi2ndf, row.smeared_p,
                                                                                     "V" if row.unfolded_p > row.smeared_p else "X"))
    print(r'\end{tabular}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory (should be the one made by unfolding.py)")
    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT.keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them (default).",
                        default=["all"])
    parser.add_argument("--outputDir",
                        default=None,
                        help='Output directory (default is the source dir)')
    parser.add_argument("--doBinnedPlots",
                        action='store_true',
                        help='Do all lambda/pt binned plots')

    region_group = parser.add_argument_group('Region selection')
    region_group.add_argument("--doAllRegions",
                              action='store_true',
                              help='Do unfolding for all regions (dijet, Z+J, groomed, ungroomed)')
    region_group.add_argument("--doDijetCentral",
                              action='store_true',
                              help='Do unfolding for dijet (central) jets')
    region_group.add_argument("--doDijetForward",
                              action='store_true',
                              help='Do unfolding for dijet (forward) jets')
    region_group.add_argument("--doDijetCentralGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (central) jets')
    region_group.add_argument("--doDijetForwardGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (forward) jets')
    region_group.add_argument("--doZPJ",
                              action='store_true',
                              help='Do unfolding for Z+jet jets')
    region_group.add_argument("--doZPJGroomed",
                              action='store_true',
                              help='Do unfolding for groomed Z+jet jets')

    args = parser.parse_args()

    if args.doAllRegions:
        for x in ['doDijetCentral', 'doDijetForward', 'doDijetCentralGroomed', 'doDijetForwardGroomed', 'doZPJ', 'doZPJGroomed']:
            setattr(args, x, True)

    regions = setup_regions(args)
    if len(regions) == 0:
        raise RuntimeError("Need at least 1 region")

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    if len(angles) == 0:
        raise RuntimeError("Need at least 1 angle")

    jet_algo = "AK4 PF PUPPI"
    if "ak8puppi" in args.source:
        jet_algo = "AK8 PF PUPPI"

    has_data = not ('_MC_all' in args.source or '_MC_split' in args.source)

    all_chi2_stats = []

    # Iterate through regions & variables
    for region in regions:
        region_dir = os.path.join(args.source, region['name'])
        if not os.path.isdir(region_dir):
            print("! Warning ! cannot find region dir", region_dir, '- skipping region')
            continue

        for angle in angles:
            angle_output_dir = "%s/%s/%s" % (args.source, region['name'], angle.var)
            if not os.path.isdir(angle_output_dir):
                print("! Warning ! cannot find angle dir", angle_output_dir, '- skipping angle', angle.var)
                continue

            # TODO: put this in a method / class
            # Check if ROOT file exists
            root_filename = os.path.join(angle_output_dir, "unfolding_result.root")
            if not os.path.isfile(root_filename):
                print("! Warning ! cannot fine unfolding ROOT file", root_filename, ' - skipping angle')

            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
            print("*"*120)
            print("Region/var: %s" % (append))
            print("*"*120)

            input_tfile = cu.TFileCacher(root_filename)  # keep this here otherwise crashes
            unpack_dict = unpack_unfolding_root_file(input_tfile, region, angle)
            # we have
            # - unfolder
            # - unreg_unfolder
            # - alt_unfolder
            # - experimental_systematics[] in unfolder
            # - model_systematics[] in unfolder
            # - pdf[] in unfolder

            # MAKE ALL THE PLOTS
            # ------------------------------------------------------------------
            setup = Setup(jet_algo=jet_algo,
                          region=region,
                          angle=angle,
                          output_dir=angle_output_dir,
                          has_data=has_data)

            hist_bin_chopper = None
            if args.doBinnedPlots:
                hist_bin_chopper = do_all_plots_per_region_angle(setup, unpack_dict)

            # Do a 1D summary plot, with all the normalised plots with bins divided by their width
            # (unlike the standard plot from MyUnfolderPlotter, which is absolute)
            do_all_big_1d_plots_per_region_angle(setup,
                                                 unpack_dict,
                                                 hist_bin_chopper)

            store_bottom_line_stats(all_chi2_stats, setup, unpack_dict)

    df_stats = pd.DataFrame(all_chi2_stats)
    df_stats['region'] = df_stats['region'].astype('category')
    df_stats['angle'] = df_stats['angle'].astype('category')
    print(df_stats.head())
    print_chi2_table(df_stats)

