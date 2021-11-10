#!/usr/bin/env python


"""
Do all the unfolding plots: per pT bin, per lambda bin, summary plot
"""


from __future__ import print_function, division

import os
import sys
import argparse
import math
from array import array
import pandas as pd
from copy import copy, deepcopy
import numpy as np
import scipy
from itertools import chain
from glob import glob
from collections import defaultdict
import inspect
import warnings

# for webpage
from jinja2 import Environment, FileSystemLoader

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, HistBinChopper, unpickle_region, PtVarPerPtBinning
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import setup_regions_from_argparse
import metric_calculators as metrics

My_Style.cd()
pd.set_option('display.max_columns', None)
os.nice(10)

# Use rootpy to throw exceptions on ROOT errors, but need DANGER enabled
# Doesn't work with python 3.8 for now
if not (sys.version_info.major == 3 and sys.version_info.minor >= 8):
    import rootpy
    # import rootpy.logger.magic as M; M.DANGER.enabled = True


# monkey-patch warning formatter
warnings.formatwarning = cu._formatwarning

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()

# When using memory_profiler/mprof, handy to have @profile to mark functions
# But when running normally we want to pass through without manually commenting out,
# so define our own decorator instead that does nothing
# Note that if this module is imported elsewhere,
# profile will be in __builtins__, not directly in locals()
# so we check multiple places
if not any(['profile' in locals(),
            'profile' in globals(),
            (isinstance(globals()['__builtins__'], dict) and 'profile' in globals()['__builtins__']),
            'profile' in dir(__builtins__)]):
    print("I have no memory_profiler @profile decorator in do_unfolding_plots, creating my own instead")

    def profile(func):
        return func


class Setup(object):
    """Loads of common consts, useful for plotting etc"""

    def __init__(self, jet_algo, region, angle, output_dir='.', has_data=False, is_ave_pt_binning=False, is_preliminary=True):
        self.jet_algo = jet_algo
        self.region = region
        do_zpj = "ZPlusJets" in region['name']
        do_dijet = "Dijet" in region['name']
        self.lumi = cu.get_lumi_str(do_dijet=do_dijet, do_zpj=do_zpj)
        self.pt_var_str = "#LT p_{T}^{jet} #GT" if is_ave_pt_binning else "p_{T}^{jet}"
        self.pt_str = self.pt_var_str + " [GeV]"
        self.has_data = has_data
        self.angle = angle
        angle_prepend = "Groomed " if "groomed" in region['name'] else ""
        lower_angle_name = qgc.lower_angle_name(angle)
        # for plot axis titles
        self.angle_str = "{prepend}{name} ({lambda_str})".format(prepend=angle_prepend,
                                                                 name=lower_angle_name if angle_prepend != "" else angle.name,
                                                                 lambda_str=angle.lambda_str)
        # self.particle_title = "Particle-level " + self.angle_str
        self.particle_title = self.angle_str

        angle_prepend = "groomed " if "groomed" in region['name'] else ""
        self.detector_title = "Detector-level {prepend}{name} ({lambda_str})".format(prepend=angle_prepend,
                                                                                     name=lower_angle_name,
                                                                                     lambda_str=angle.lambda_str)
        self.pt_bin_normalised_differential_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % angle.lambda_str
        self.pt_bin_detector_normalised_differential_label = "#frac{1}{dN/dp_{T}} #frac{d^{2}N}{dp_{T} d%s}" % angle.lambda_str
        self.pt_bin_normalised_differential_times_width_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % angle.lambda_str
        self.pt_bin_unnormalised_differential_label = "#frac{1}{dN/dp_{T}} #frac{d^{2}N}{dp_{T} d%s}" % angle.lambda_str  # FIXME
        self.pt_bin_unnormalised_differential_label = "N"
        self.lambda_bin_normalised_differential_label = "#frac{1}{d#sigma/d%s} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)
        self.lambda_bin_unnormalised_differential_label = "#frac{1}{dN/d%s} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)  # FIXME
        self.lambda_bin_unnormalised_differential_label = "N"
        self.output_dir = output_dir
        self.output_fmt = 'pdf'
        self.append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
        self.is_preliminary = is_preliminary

    @property
    def paper_str(self):
        return "_paper" if not self.is_preliminary else ""


PLOT_STYLES = dict(
    # gen_colour=ROOT.kRed,
    gen_colour=qgc.MGPY_QCD_COLOUR,
    gen_marker=cu.Marker.get('square', filled=False),
    gen_marker_size=1.1,

    unfolded_basic_colour=ROOT.kAzure+7,
    unfolded_stat_colour=ROOT.kAzure+7,
    unfolded_total_colour=ROOT.kBlack,
    unfolded_unreg_colour=ROOT.kOrange-3,
    unfolded_marker_size=1.1,

    #alt_gen_colour=ROOT.kOrange-3,
    # alt_gen_colour=ROOT.kViolet+1,
    alt_gen_colour=qgc.HERWIGPP_QCD_COLOUR,
    alt_unfolded_colour=ROOT.kOrange-3,
    alt_unfolded_total_colour=ROOT.kOrange-7,
    alt_reco_colour=ROOT.kViolet+1,
    # alt_reco_colour=ROOT.kOrange-3,
    alt_gen_marker=cu.Marker.get('triangleUp', filled=False),
    alt_gen_marker_size=1.1,

    # reco_mc_colour=ROOT.kGreen+2,
    # reco_mc_colour=ROOT.kAzure-7,
    # reco_data_colour=ROOT.kRed,

    reco_data_colour=ROOT.kBlack,
    reco_mc_colour=ROOT.kRed,
    # reco_mc_colour=ROOT.kMagenta+1,
    reco_unfolding_input_colour=ROOT.kRed+3,
    reco_folded_unfolded_colour=ROOT.kAzure+1,
    reco_folded_mc_truth_colour=ROOT.kGreen+2,
    default_palette=ROOT.kBird,

    rsp_colour=ROOT.kGray+3,
    scale_colour=ROOT.kTeal+3,
    pdf_colour=ROOT.kOrange+4,

    template_colour=ROOT.kAzure+1,
)


def calc_chi2_stats(one_hist, other_hist, cov_matrix, theory_uncertainty=False):
    one_vec, one_err = cu.th1_to_ndarray(one_hist, False)
    # print(one_err)
    other_vec, other_err = cu.th1_to_ndarray(other_hist, False)
    delta = one_vec - other_vec
    if isinstance(cov_matrix, ROOT.TH2):
        v, cov_err = cu.th2_to_ndarray(cov_matrix)
    else:
        v = cov_matrix
    # print("delta:", delta)
    # v = np.diag(np.diag(v))  # turn off correlations
    # print("v:", v)
    if theory_uncertainty:
      #print("delta:", delta)
      #print("v:", v)
      #print("other_err", other_err)
      other_err_m = np.diag(np.square(other_err[0]))
      #print("other_err_m", other_err_m)
      v = np.add(v, other_err_m)
      #print("v", v)
    try:
        v_inv = np.linalg.inv(v)
    except np.linalg.LinAlgError:
        # print("Trying pseudo-inverse instead")
        v_inv = np.linalg.pinv(v, rcond=1E-30)
    inter = v_inv.dot(delta.T)
    # print("parts:", delta * inter.T)
    chi2 = delta.dot(inter)[0][0]
    ndof = delta.shape[1]
    p = 1-scipy.stats.chi2.cdf(chi2, int(ndof))
    return chi2, ndof, p


def _make_thumbnail_canvas(plot):
    orig_w, orig_h = plot.default_canvas_size
    plot.default_canvas_size = (int(orig_w*0.5), int(orig_h*0.5))

#TODO: centralise this - same code in print_detector_unfolded.py!
def scale_ematrix_by_bin_widths(ematrix, widths):
    this_widths = widths.reshape(len(widths), 1)
    return ematrix * this_widths * this_widths.T


def get_correlated_mean_err(hist, ematrix, is_density=True):
    contents, _ = cu.th1_to_ndarray(hist)
    contents = contents.reshape(-1)  # necessary for jax to avoid shape discrepancy
    centers = cu.get_th1_bin_centers(hist)
    cov_matrix, _ = cu.th2_to_ndarray(ematrix)
    if is_density:
        widths = cu.get_th1_bin_widths(hist)
        bin_areas = contents*widths
        # need to scale ematrix by bin areas
        cov_matrix = scale_ematrix_by_bin_widths(cov_matrix, widths)
    else:
        bin_areas = contents
    # print(bin_areas.shape, centers.shape, cov_matrix.shape)
    mean = float(metrics.calc_mean_jax(bin_areas, centers))
    mean_err = float(metrics.calc_mean_correlated_error_jax(bin_areas, centers, cov_matrix))
    return mean, mean_err


def get_correlated_rms_err(hist, ematrix, is_density=True):
    contents, _ = cu.th1_to_ndarray(hist)
    contents = contents.reshape(-1)  # necessary for jax to avoid shape discrepancy
    centers = cu.get_th1_bin_centers(hist)
    cov_matrix, _ = cu.th2_to_ndarray(ematrix)
    if is_density:
        widths = cu.get_th1_bin_widths(hist)
        bin_areas = contents*widths
        # need to scale ematrix by bin areas
        cov_matrix = scale_ematrix_by_bin_widths(cov_matrix, widths)
    else:
        bin_areas = contents
    # print(bin_areas.shape, centers.shape, cov_matrix.shape)
    rms = float(metrics.calc_rms_jax(bin_areas, centers))
    rms_err = float(metrics.calc_rms_correlated_error_jax(bin_areas, centers, cov_matrix))
    return rms, rms_err


def get_uncorrelated_mean_err(hist, is_density=True):
    contents, errors = cu.th1_to_ndarray(hist)
    centers = cu.get_th1_bin_centers(hist)
    if is_density:
        # need to multiply by widths, since the original hist has bin contents divided by width
        widths = cu.get_th1_bin_widths(hist)
        bin_areas = contents*widths
        bin_errors = errors*widths
    else:
        bin_areas = contents
        bin_errors = errors
    # convert to uncertainty arrays
    areas, centers = metrics.hist_values_to_uarray(bin_areas=bin_areas, bin_centers=centers, bin_errors=bin_errors)
    mean_u = metrics.calc_mean_ucert(areas, centers)
    return mean_u.nominal_value, mean_u.std_dev


def get_uncorrelated_rms_err(hist, is_density=True):
    contents, errors = cu.th1_to_ndarray(hist)
    centers = cu.get_th1_bin_centers(hist)
    if is_density:
        # need to multiply by widths, since the original hist has bin contents divided by width
        widths = cu.get_th1_bin_widths(hist)
        bin_areas = contents*widths
        bin_errors = errors*widths
    else:
        bin_areas = contents
        bin_errors = errors
    # convert to uncertainty arrays
    areas, centers = metrics.hist_values_to_uarray(bin_areas=bin_areas, bin_centers=centers, bin_errors=bin_errors)
    rms_u = metrics.calc_rms_ucert(areas, centers)
    return rms_u.nominal_value, rms_u.std_dev


class BinnedPlotter(object):
    def __init__(self):
        self.plots_dict = defaultdict(list)  # use to cache plot names per function name

        self.line_width = 2
        self.plot_styles = PLOT_STYLES

    def save_plot(self, plot, filename):
        plot.save(filename)
        func_name = inspect.currentframe().f_back.f_code.co_name
        self.plots_dict[func_name].append(filename)

    @staticmethod
    def check_entries(entries, message=""):
        """Check that at least 1 Contribution has something in it, that is non-zero"""
        has_entries = [c.obj.GetEntries() > 0 for c in entries]
        if not any(has_entries):
            if message:
                print("Skipping 0 entries (%s)" % message)
            else:
                print("Skipping 0 entries")
            return False

        max_bin = max([c.obj.GetMaximum() for c in entries])
        min_bin = min([c.obj.GetMinimum() for c in entries])
        if max_bin == min_bin:
            if message:
                print("Skipping min=max hists (%s)" % message)
            else:
                print("Skipping min=max hists")
            return False

        return True


class GenPtBinnedPlotter(BinnedPlotter):
    def __init__(self, setup, bins, hist_bin_chopper, unfolder):
        self.setup = setup
        self.region = setup.region  # just to make life easier
        self.bins = bins
        self.hist_bin_chopper = hist_bin_chopper

        self.pt_bin_plot_args = dict(
            what="hist",
            xtitle=self.setup.particle_title,
            has_data=self.setup.has_data,
            # ylim=[0, None],
            subplot_type='ratio',
            subplot_title="* / %s" % (self.region['mc_label']),
            # subplot_limits=(0, 2) if self.setup.has_data else (0.75, 1.25),
            subplot_limits=(0, 3)
        )
        self.unfolder = unfolder
        super().__init__()

    def _modify_plot(self, this_plot):
        if self.setup.output_fmt == "gif":
            _make_thumbnail_canvas(this_plot)
        this_plot.legend.SetX1(0.6)
        this_plot.legend.SetY1(0.7)
        this_plot.legend.SetX2(0.98)
        this_plot.legend.SetY2(0.88)
        this_plot.left_margin = 0.16
        this_plot.y_padding_max_linear = 1.8
        this_plot.lumi = self.setup.lumi

    def _modify_plot_paper(self, this_plot):
        self._modify_plot(this_plot)
        this_plot.is_preliminary = self.setup.is_preliminary

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
                )).replace("1000","1").replace("4000 GeV","4 TeV")
        return title

    def plot_unfolded_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('hist_truth', **hbc_args)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)

            # unnormalised version
            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator",
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % self.unfolder.tau,
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'],  # marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % self.unfolder.tau,
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.subplot_title = "* / Generator"
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            # key = self.hist_bin_chopper._generate_key("unfolded_stat_err",
            #                                           ind=ibin,
            #                                           axis='pt',
            #                                           do_norm=True,
            #                                           do_div_bin_width=True,
            #                                           binning_scheme='generator')
            # print("plot_unfolded_normalised", key)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            # cu.print_th1_bins(unfolded_hist_bin_stat_errors)

            if 'syst_err_Total' in self.hist_bin_chopper.objects:
                key = self.hist_bin_chopper._generate_key("syst_err_Total",
                                                          ind=ibin,
                                                          axis='pt',
                                                          do_norm=True,
                                                          do_div_bin_width=True,
                                                          binning_scheme='generator')

                print("plot_unfolded_normalised", key)
                unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('syst_err_Total', **hbc_args)
            else:
                # print("plot_unfolded_normalised total")
                unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
            # cu.print_th1_bins(unfolded_hist_bin_total_errors)

            for ix in range(1, unfolded_hist_bin_total_errors.GetNbinsX()+1):
                val, err = unfolded_hist_bin_total_errors.GetBinContent(ix), unfolded_hist_bin_total_errors.GetBinError(ix)
                if val < 0:
                    warnings.warn("plot_unfolded_normalised(): pt bin {}, hist bin {} content < {} ± {} ".format(ibin, ix, val, err))

            # get stats
            cov_matrix = self.hist_bin_chopper.get_bin_plot(self.region['unfolder'].total_ematrix_name,
                                                            ind=ibin, axis='pt', do_norm=True, do_div_bin_width=True, binning_scheme='generator')
            unfolded_mean, unfolded_mean_err = get_correlated_mean_err(unfolded_hist_bin_total_errors, cov_matrix, is_density=True)
            unfolded_rms, unfolded_rms_err = get_correlated_rms_err(unfolded_hist_bin_total_errors, cov_matrix, is_density=True)

            generator_mean, generator_mean_err = get_uncorrelated_mean_err(mc_gen_hist_bin, is_density=True)
            generator_rms, generator_rms_err = get_uncorrelated_rms_err(mc_gen_hist_bin, is_density=True)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator\nMean = %.3f #pm %.3f\nRMS = %.3f #pm %.3f" % (generator_mean, generator_mean_err, generator_rms, generator_rms_err),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Data (total unc.)\nMean = %.3f #pm %.3f\nRMS = %.3f #pm %.3f" % (unfolded_mean, unfolded_mean_err, unfolded_rms, unfolded_rms_err),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'],# marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Data (stat. unc.)",
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return

            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        xlim=qgp.calc_auto_xlim(entries[2:3]),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.subplot_title = "* / Generator"
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

            if self.setup.angle.var in ["jet_LHA_charged", "jet_thrust_charged", "jet_width_charged", "jet_thrust", "jet_width"]:
                upper_bin = [x for x in cu.get_bin_edges(mc_gen_hist_bin, 'x') if x < 0.2][-1]
                plot2 = Plot(entries,
                             ytitle=self.setup.pt_bin_normalised_differential_label,
                             title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                             xlim=(0, upper_bin),
                             **self.pt_bin_plot_args)
                self._modify_plot(plot2)
                plot2.subplot_title = "* / Generator"
                plot2.plot("NOSTACK E1")
                # plot2.set_logx(do_exponent=False)
                self.save_plot(plot2, "%s/unfolded_%s_bin_%d_divBinWidth_lowX.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_alt_truth_normalised(self, do_chi2=False, do_zoomed=True):
        data_total_errors_style = dict(label="Data",
                                       line_color=self.plot_styles['unfolded_total_colour'],
                                       line_width=self.line_width,
                                       line_style=1,
                                       marker_color=self.plot_styles['unfolded_total_colour'],
                                       marker_style=cu.Marker.get('circle'),
                                       marker_size=self.plot_styles['unfolded_marker_size'],
                                       leg_draw_opt="LEP")
        data_stat_errors_style = dict(label="Data (stat. unc.)",
                                      line_color=self.plot_styles['unfolded_stat_colour'],
                                      line_width=self.line_width,
                                      line_style=1,
                                      marker_color=self.plot_styles['unfolded_stat_colour'],
                                      marker_style=cu.Marker.get('circle'),
                                      marker_size=0.0001,
                                      leg_draw_opt="EP")  # you need a non-0 marker to get the horizontal bars at the end of errors

        mc_style = dict(label=self.region['mc_label'],
                        line_color=self.plot_styles['gen_colour'],
                        line_width=self.line_width,
                        marker_color=self.plot_styles['gen_colour'],
                        marker_size=self.plot_styles['gen_marker_size'],
                        marker_style=self.plot_styles['gen_marker'],
                        leg_draw_opt="LEP" if self.plot_styles['gen_marker_size'] > 0 else "LE")
        alt_mc_style = dict(label=self.region['alt_mc_label'],
                            line_color=self.plot_styles['alt_gen_colour'],
                            line_width=self.line_width,
                            line_style=1,
                            marker_color=self.plot_styles['alt_gen_colour'],
                            marker_size=self.plot_styles['alt_gen_marker_size'],
                            marker_style=self.plot_styles['alt_gen_marker'],
                            leg_draw_opt="LEP" if self.plot_styles['alt_gen_marker_size'] > 0 else "LE")
        #### HERE IS WHERE ANDREAS ADDED THE THEORY
        theory_style = dict(label="NLO + NLL'+ NP",
                            line_color=ROOT.kRed-9,
                            line_width=self.line_width,
                            line_style=1,
                            marker_color=ROOT.kRed-9,
                            marker_size=0,
                            marker_style=1,
                            fill_color=0,
                            fill_style=0,
                            leg_draw_opt="L")
        theory_style_legend = theory_style.copy()
        theory_style_legend["fill_color"]=ROOT.kRed-9
        theory_style_legend["fill_style"]=3395
        theory_style_legend["leg_draw_opt"]="FL"

        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            alt_mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', **hbc_args)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

            # Create copy of data to go on top of stat unc,
            # but remove vertical error bar so we can see the stat unc
            # Note that you CAN'T set it to 0, otherwise vertical lines connecting
            # bins start being drawn. Instead set it to some super small value.
            unfolded_hist_bin_total_errors_marker_noerror = unfolded_hist_bin_total_errors.Clone()  # clone to avoid restyling the original as well
            for i in range(1, unfolded_hist_bin_total_errors_marker_noerror.GetNbinsX()+1):
                unfolded_hist_bin_total_errors_marker_noerror.SetBinError(i, 1E-100)

            ### GET THEORY CURVE
            bins=np.array(unfolded_hist_bin_total_errors.GetXaxis().GetXbins())
            bins[0]=bins[1]/1.5
            if not "multiplicity" in self.setup.angle.var.lower():
              for hist in [mc_gen_hist_bin,alt_mc_gen_hist_bin,unfolded_hist_bin_stat_errors,unfolded_hist_bin_total_errors,unfolded_hist_bin_total_errors_marker_noerror]:
                hist.SetBins(hist.GetXaxis().GetNbins(),bins)
            
            theory_hist_bin = unfolded_hist_bin_total_errors.Clone()  # clone to avoid restyling the original as well
            theory_hist_bin_error = unfolded_hist_bin_total_errors.Clone()  # clone to avoid restyling the original as well
            theory_hist_bin_upper = unfolded_hist_bin_total_errors.Clone()  # clone to avoid restyling the original as well
            theory_hist_bin_lower = unfolded_hist_bin_total_errors.Clone()  # clone to avoid restyling the original as well
            var_number=["jet_puppiMultiplicity", "jet_pTD","jet_LHA", "jet_width", "jet_thrust", 
            		"jet_puppiMultiplicity_charged", "jet_pTD_charged", "jet_LHA_charged", "jet_width_charged", "jet_thrust_charged"].index(self.setup.angle.var)+1
            if "Z" in self.region['label'] and var_number in [3,4,5,8,9,10]:
              for post in ["_fix",""]:
                theory_file="/nfs/dust/cms/user/hinzmann/qganalysis/CMSSW_10_2_17/src/zjet_angularities-master-to_send_to_CMS-thr_vs_exp_cvt"+post+"/to_send_to_CMS/thr_vs_exp_cvt"+post+"/"
                theory_file+="data_R"+self.setup.jet_algo[-1]+"/"
                theory_file+="RSIG_THR_RES_NP_"
                theory_file+="d"+("1" if "4" in self.setup.jet_algo else "2")+("2" if "groomed" in self.setup.region['name'] else "1")+"-"
                theory_file+="x"+("0" if var_number<10 else "")+str(var_number)+"-"
                theory_file+="y"+("0" if ibin<9 else "")+str(ibin+1)
                theory_file+=".dat"
                if os.path.exists(theory_file): break
                print(theory_file)
              igenbin=0
              for line in open(theory_file).readlines():
                if line.count(".")==5:
                  igenbin+=1
                  theory_hist_bin.SetBinContent(igenbin, float(line.split("\t")[2]))
                  theory_hist_bin.SetBinError(igenbin, 1e-100)
                  theory_hist_bin_error.SetBinContent(igenbin, float(line.split("\t")[2]))
                  theory_hist_bin_error.SetBinError(igenbin, (float(line.split("\t")[3])+float(line.split("\t")[4]))/2.)
                  theory_hist_bin_upper.SetBinContent(igenbin, float(line.split("\t")[2])+float(line.split("\t")[3]))
                  theory_hist_bin_upper.SetBinError(igenbin, 0)
                  theory_hist_bin_lower.SetBinContent(igenbin, float(line.split("\t")[2])-float(line.split("\t")[4]))
                  theory_hist_bin_lower.SetBinError(igenbin, 0)
              assert(igenbin==theory_hist_bin.GetNbinsX())
            
            data_entries = [
                Contribution(unfolded_hist_bin_total_errors, **data_total_errors_style),
                #Contribution(unfolded_hist_bin_stat_errors, **data_stat_errors_style),
                # do data with black marker to get it on top
                #Contribution(unfolded_hist_bin_total_errors_marker_noerror, **data_total_errors_style),
            ]

            # For subplot to ensure only MC errors drawn, not MC+data
            data_no_errors = unfolded_hist_bin_total_errors_marker_noerror.Clone()
            cu.remove_th1_errors(data_no_errors)

            this_mc_style = deepcopy(mc_style)
            this_alt_mc_style = deepcopy(alt_mc_style)
            this_theory_style = deepcopy(theory_style)
            this_theory_style_legend = deepcopy(theory_style_legend)

            max_rel_err = 0.9 # 0.5 if "multiplicity" in self.setup.angle.var.lower() else -1

            # Calculate chi2 between data and MCs if desired
            if do_chi2:
                # print("unfolded_alt_truth bin", ibin)
                ematrix = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.unfolder.total_ematrix_name, **hbc_args)
                print("Correlation,", str(bin_edge_low)+" < PT < "+str(bin_edge_high)+" GeV")
                print(str([[unfolded_hist_bin_total_errors.GetXaxis().GetBinLowEdge(i+1),unfolded_hist_bin_total_errors.GetXaxis().GetBinUpEdge(i+1)]\
                          for i in range(unfolded_hist_bin_total_errors.GetNbinsX())]).replace("[","").replace("]","\n").replace(",",""))
                print(str(cu.th2_to_ndarray(ematrix)[0]).replace("[","").replace("]","").replace(",",""))
                # stats are chi2, ndof, p
                mc_stats = calc_chi2_stats(unfolded_hist_bin_total_errors, mc_gen_hist_bin, ematrix)
                alt_mc_stats = calc_chi2_stats(unfolded_hist_bin_total_errors, alt_mc_gen_hist_bin, ematrix)
                theory_stats = calc_chi2_stats(unfolded_hist_bin_total_errors, theory_hist_bin_error, ematrix, True)
                # print(mc_stats)
                # print(alt_mc_stats)
                nbins = sum([1 for i in range(1, unfolded_hist_bin_total_errors.GetNbinsX()+1)
                             if unfolded_hist_bin_total_errors.GetBinContent(i) != 0 and \
                             unfolded_hist_bin_total_errors.GetBinError(i)/unfolded_hist_bin_total_errors.GetBinContent(i) < max_rel_err])
                # reduced_chi2 = mc_stats[0] / nbins
                # alt_reduced_chi2 = alt_mc_stats[0] / nbins

                n_sig_fig = 2
                chi2_template = "\n#lower[-0.1]{{(#chi^{{2}}/N_{{bins}} = {chi2:g}/{nbins:d})}}"
                this_mc_style['label'] += chi2_template.format(chi2=cu.nsf(mc_stats[0], n_sig_fig), nbins=nbins)
                this_alt_mc_style['label'] += chi2_template.format(chi2=cu.nsf(alt_mc_stats[0], n_sig_fig), nbins=nbins)
                this_theory_style_legend['label'] += chi2_template.format(chi2=cu.nsf(theory_stats[0], n_sig_fig), nbins=nbins)

            mc_entries = [
                Contribution(mc_gen_hist_bin, subplot=data_no_errors, **this_mc_style),
                Contribution(alt_mc_gen_hist_bin, subplot=data_no_errors, **this_alt_mc_style),
            ]
            if "Z" in self.region['label'] and var_number in [3,4,5,8,9,10]:
              mc_entries += [Contribution(theory_hist_bin, subplot=data_no_errors, **this_theory_style) ]

            entries = [
                # Draw MC
                *mc_entries,
                # Draw data after to put on top of MC
                *data_entries
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return

            ymin = 0
            #if np.any(cu.th1_to_ndarray(unfolded_hist_bin_total_errors)[0] < 0):
            #    ymin = None  # let it do its thing and auto calc ymin
            max_x = qgp.calc_auto_xlim(entries[2:3], max_rel_err=max_rel_err)[1]
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        legend=True,
                        xlim=[bins[0],max_x],  # set x lim to where data is non-0 ### ANDREAS changed from 0.5
                        ylim=[ymin, None],
                        **self.pt_bin_plot_args)

            plot.subplot_title = qgc.SIM_DATA_STR
            #print((self.setup.region['name'],self.setup.angle.var,self.setup.jet_algo,ibin))
            if not (self.setup.region['name'],self.setup.angle.var,self.setup.jet_algo,ibin) in [
                ('ZPlusJets', 'jet_puppiMultiplicity', 'AK4', 3),
                ('Dijet_central', 'jet_puppiMultiplicity', 'AK4', 3),
                ('ZPlusJets', 'jet_pTD', 'AK4', 3),
                ('Dijet_central', 'jet_pTD', 'AK4', 3),
                ('ZPlusJets', 'jet_thrust', 'AK4', 3),
                ('Dijet_central', 'jet_thrust', 'AK4', 3),
                ('ZPlusJets', 'jet_width', 'AK4', 3),
                ('Dijet_central', 'jet_width', 'AK4', 3),
                ('ZPlusJets', 'jet_LHA', 'AK4', 3),
                ('Dijet_central', 'jet_LHA', 'AK4', 3),
                ('ZPlusJets', 'jet_LHA', 'AK4', 8),
                ('Dijet_central', 'jet_LHA', 'AK4', 12),
                ('ZPlusJets', 'jet_LHA', 'AK8', 3),
                ('Dijet_central', 'jet_LHA', 'AK8', 3),
                ('ZPlusJets', 'jet_LHA_charged', 'AK4', 3),
                ('Dijet_central', 'jet_LHA_charged', 'AK4', 3),
                ('ZPlusJets_groomed', 'jet_LHA', 'AK4', 3),
                ('Dijet_central_groomed', 'jet_LHA', 'AK4', 3),
                ]:
                plot.is_supplementary = True
            self._modify_plot_paper(plot)

            # disable adding objects to legend & drawing - we'll do it manually
            plot.do_legend = False
            final_reading_factor = 1.18
            plot.legend.SetTextSize(0.03*final_reading_factor)
            if "Z" in self.region['label'] and var_number in [3,4,5,8,9,10]:
              plot.legend.SetY1(0.68)
            else:
              plot.legend.SetY1(0.7)
            plot.legend.SetX1(0.57)
            plot.legend.SetX2(0.93)
            # plot.legend.SetEntrySeparation(0.005)
            draw_opts = "NOSTACK E1 X0"
            plot.plot(draw_opts, draw_opts)
            #plot.container.SetMinimum(plot.container.GetMinimum()/2.)
            #plot.set_logy(do_exponent=False)
            if not "multiplicity" in self.setup.angle.var.lower():
              plot.set_logx(do_exponent=False)
            plot.container.SetMinimum(0)

            plot.main_pad.cd()
            if "Z" in self.region['label'] and var_number in [3,4,5,8,9,10]:
              theory_hist_bin_upper.SetFillStyle(3395)
              theory_hist_bin_upper.SetFillColor(ROOT.kRed-9)
              theory_hist_bin_upper.SetLineWidth(0)
              theory_hist_bin_upper.SetMarkerSize(0)
              theory_hist_bin_upper.Draw("F SAME")
              theory_hist_bin_lower.SetFillStyle(1001)
              theory_hist_bin_lower.SetFillColor(10)
              theory_hist_bin_lower.SetLineWidth(0)
              theory_hist_bin_lower.SetMarkerSize(0)
              theory_hist_bin_lower.Draw("F SAME")
              mc_entries=mc_entries[:-1]+[Contribution(theory_hist_bin_error, subplot=data_no_errors, **this_theory_style_legend) ]
              ROOT.gStyle.SetNdivisions(510,"Y")
              data_no_errors.Draw("AXIS SAME")
              #plot.main_pad.RedrawAxis()
              plot.container.Draw(draw_opts+" A SAME")
            graphs=[]
            for e in plot.contributions:
              graphs+=[ROOT.TGraphErrors(e.obj)]
              graphs[-1].Draw("PZ SAME")
            for e in plot.contributions:
              graphs+=[ROOT.TGraphErrors(e.obj)]
              for p in range(graphs[-1].GetN()):
                graphs[-1].SetPointError(p,0,graphs[-1].GetErrorY(p))
              graphs[-1].Draw("E1 SAME")
            
            dummy_graphs = qgp.do_fancy_legend(chain(data_entries[:2], mc_entries), plot, use_splitline=False)

            plot.canvas.cd()
            plot.legend.Draw()
	    
            # Create hists for data with error region for ratio
            # Easiest way to get errors right is to do data (with 0 errors)
            # and divide by data (with errors), as if you had MC = data with 0 error
            data_stat_ratio = data_no_errors.Clone()
            data_stat_ratio.Divide(unfolded_hist_bin_stat_errors)
            data_stat_ratio.SetFillStyle(3545)
            data_stat_ratio.SetFillColor(self.plot_styles['unfolded_stat_colour'])
            data_stat_ratio.SetLineWidth(0)
            data_stat_ratio.SetMarkerSize(0)

            data_total_ratio = data_no_errors.Clone()
            data_total_ratio.Divide(unfolded_hist_bin_total_errors)
            data_total_ratio.SetFillStyle(3354)
            data_total_ratio.SetFillColor(16)
            data_total_ratio.SetLineWidth(0)
            data_total_ratio.SetMarkerSize(0)

            # now draw the data error shaded area
            # this is a bit hacky - basically draw them on the ratio pad,
            # then redraw the existing hists & line to get them ontop
            # note that we use "same" for all - this is to keep the original axes
            # (we may want to rethink this later?)
            plot.subplot_pad.cd()

            if "Z" in self.region['label'] and var_number in [3,4,5,8,9,10]:
              theory_hist_bin_upper_ratio = theory_hist_bin_upper.Clone()
              theory_hist_bin_upper_ratio.Divide(data_no_errors)
              theory_hist_bin_upper_ratio.Draw("F SAME")
              theory_hist_bin_lower_ratio = theory_hist_bin_lower.Clone()
              theory_hist_bin_lower_ratio.Divide(data_no_errors)
              theory_hist_bin_lower_ratio.Draw("F SAME")
              plot.subplot_pad.RedrawAxis()
              #data_stat_ratio.Draw("AXIS SAME")

            draw_opt = "E2 SAME"
            data_stat_ratio.Draw(draw_opt)
            data_total_ratio.Draw(draw_opt)
            plot.subplot_line.Draw()
            plot.subplot_container.Draw("SAME" + draw_opts)
            for e in [data_stat_ratio,data_total_ratio]:
              graphs+=[ROOT.TGraphErrors(e)]
              graphs[-1].Draw("2P SAME")
            for e in plot.subplot_contributions:
              graphs+=[ROOT.TGraphErrors(e)]
              graphs[-1].Draw("PZ SAME")
            for e in plot.subplot_contributions:
              graphs+=[ROOT.TGraphErrors(e)]
              for p in range(graphs[-1].GetN()):
                graphs[-1].SetPointError(p,0,graphs[-1].GetErrorY(p))
              graphs[-1].Draw("E1 SAME")

            # Add subplot legend
            x_left = 0.25
            y_bottom = 0.78
            width = 0.60
            height = 0.15
            plot.subplot_legend = ROOT.TLegend(x_left, y_bottom, x_left+width, y_bottom+height)
            # plot.subplot_legend = ROOT.TLegend(width, height, width, height)  # automatic placement doesn't work
            # plot.subplot_legend.AddEntry(data_total_ratio, "Total uncert.", "F")
            # plot.subplot_legend.AddEntry(data_stat_ratio, "Stat. uncert.", "F")
            plot.subplot_legend.AddEntry(data_total_ratio, qgc.DATA_TOTAL_UNC_STR, "F")
            plot.subplot_legend.AddEntry(data_stat_ratio, qgc.DATA_STAT_UNC_STR, "F")
            plot.subplot_legend.SetTextSize(0.085)
            plot.subplot_legend.SetFillStyle(0)
            plot.subplot_legend.SetNColumns(2)
            plot.subplot_legend.Draw()
            plot.canvas.cd()

            stp = self.setup
            fname = f"unfolded_{stp.append}_alt_truth_bin_{ibin:d}_divBinWidth{stp.paper_str}.{stp.output_fmt}"
            self.save_plot(plot, os.path.join(stp.output_dir, fname))

            if do_zoomed:
                if self.setup.angle.var in ["jet_thrust_charged", "jet_width_charged", "jet_thrust", "jet_width"]:
                    # plot.ylim = (1E-5)
                    plot.y_padding_max_log = 50
                    plot.y_padding_min_log = 0.5
                    plot.ylim = None
                    plot.set_logy(do_exponent=False, do_more_labels=False)
                    fname = f"unfolded_{stp.append}_alt_truth_bin_{ibin:d}_divBinWidth_logY.{stp.output_fmt}"
                    self.save_plot(plot, os.path.join(stp.output_dir, fname))

                if self.setup.angle.var in ["jet_LHA_charged", "jet_thrust_charged", "jet_width_charged", "jet_thrust", "jet_width"]:
                    bin_edges = cu.get_bin_edges(mc_gen_hist_bin, 'x')
                    # get the bin edge thats smallest between 0.2, and 5th bin
                    bin_lt_lim = [x for x in bin_edges if x < 0.2][-1]
                    upper_bin = min(bin_edges[5], bin_lt_lim)
                    plot2 = Plot(entries,
                                 ytitle=self.setup.pt_bin_normalised_differential_label,
                                 title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                                 xlim=(0, upper_bin),
                                 **self.pt_bin_plot_args)
                    self._modify_plot(plot2)
                    plot2.subplot_title = "* / Generator"
                    plot2.plot("NOSTACK E1")
                    # plot2.set_logx(do_exponent=False)
                    fname = f"unfolded_{stp.append}_alt_truth_bin_{ibin:d}_divBinWidth_lowX.{stp.output_fmt}"
                    self.save_plot(plot2, os.path.join(stp.output_dir, fname))


    def plot_unfolded_with_unreg_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
            unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unreg_unfolded', **hbc_args)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % self.unfolder.tau,
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'],  #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unreg_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = 0) (total unc.)",
                             line_color=self.plot_styles['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_unreg_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        xlim=qgp.calc_auto_xlim(entries),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.subplot_title = "* / Generator"
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_%s_with_unreg_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_unreg_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('hist_truth', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)
            unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unreg_unfolded', **hbc_args)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % self.unfolder.tau,
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'],  #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unreg_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = 0) (total unc.)",
                             line_color=self.plot_styles['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_unreg_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        xlim=qgp.calc_auto_xlim(entries),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.subplot_title = "* / Generator"
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_with_unreg_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_template_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            alt_mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
            # unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unreg_unfolded', **hbc_args)
            template_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('truth_template', **hbc_args)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(alt_mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['alt_mc_label']),
                             line_color=self.plot_styles['alt_gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['alt_gen_colour'], marker_size=0,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % self.unfolder.tau,
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'],  #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                # Contribution(unreg_unfolded_hist_bin_total_errors,
                #              label="Unfolded (#tau = 0) (total unc.)",
                #              line_color=self.plot_styles['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                #              marker_color=self.plot_styles['unfolded_unreg_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                #              subplot=mc_gen_hist_bin),
                Contribution(template_hist_bin_total_errors,
                             label="Template",
                             line_color=self.plot_styles['template_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['template_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),

            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        xlim=qgp.calc_auto_xlim(entries),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.subplot_title = "* / Generator"
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_%s_with_template_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_template_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('hist_truth', **hbc_args)
            alt_mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('alt_hist_truth', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)
            # unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unreg_unfolded', **hbc_args)
            template_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('truth_template', **hbc_args)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(alt_mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['alt_mc_label']),
                             line_color=self.plot_styles['alt_gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['alt_gen_colour'], marker_size=0,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                # Contribution(unreg_unfolded_hist_bin_total_errors,
                #              label="Unfolded (#tau = 0) (total unc.)",
                #              line_color=self.plot_styles['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                #              marker_color=self.plot_styles['unfolded_unreg_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                #              subplot=mc_gen_hist_bin),
                Contribution(template_hist_bin_total_errors,
                             label="Template",
                             line_color=self.plot_styles['template_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['template_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),

            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        xlim=qgp.calc_auto_xlim(entries),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.subplot_title = "* / Generator"
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_with_template_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_alt_response_normalised(self, alt_unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            # Note that we have to use the hist_bin_chopper in alt_unfolder to get
            # the correct total errors on each binned hist, since it is
            # constructed in a normalised way
            alt_unfolded_hist_bin_stat_errors = alt_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            alt_unfolded_hist_bin_total_errors = alt_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                # Contribution(alt_mc_gen_hist_bin,
                #              label="Generator (%s)" % (self.region['alt_mc_label']),
                #              line_color=alt_gen_colour, line_width=self.line_width, line_style=2,
                #              marker_color=alt_gen_colour, marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(alt_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)\n(%s response matrix)" % (alt_unfolder.tau, self.region['alt_mc_label']),
                             line_color=self.plot_styles['alt_unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['alt_unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(alt_unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (alt_unfolder.tau, self.region['alt_mc_label']),
                             line_color=self.plot_styles['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['alt_unfolded_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        xlim=qgp.calc_auto_xlim(entries),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_%s_alt_response_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_alt_response_truth_normalised(self, alt_unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            alt_mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', **hbc_args)
            # Note that we have to use the hist_bin_chopper in alt_unfolder to get
            # the correct total errors on each binned hist, since it is
            # constructed in a normalised way
            alt_unfolded_hist_bin_total_errors = alt_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(alt_mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['alt_mc_label']),
                             line_color=self.plot_styles['alt_gen_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_styles['alt_gen_colour'], marker_size=0,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(alt_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (alt_unfolder.tau, self.region['alt_mc_label']),
                             line_color=self.plot_styles['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['alt_unfolded_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            this_plot_args = {k: v for k, v in self.pt_bin_plot_args.items()}
            this_plot_args['subplot_title'] = '#splitline{* / Gen}{(%s)}' % (self.region['mc_label'])
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        xlim=qgp.calc_auto_xlim(entries),
                        **this_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_%s_alt_response_truth_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_scale_systs_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            syst_entries = []
            hbc_args = dict(ind=ibin, binning_scheme='generator')

            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)

            for syst_dict in self.region['scale_systematics']:
                syst_unfolder = syst_dict['unfolder']
                syst_label = syst_dict['label']

                # Get binned hists from the scale unfolder, since the error bars may have been setup specially
                syst_unfolded_hist_bin = syst_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

                # syst_truth_hist_bin = syst_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)

                syst_entries.extend([
                    Contribution(syst_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (total unc.) (%s)" % (syst_unfolder.tau, syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=mc_gen_hist_bin),  # use this as ref as scale variations are only different response matrices
                    # Contribution(syst_truth_hist_bin,
                    #              label="MC truth (%s)" % (syst_label),
                    #              line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                    #              marker_color=syst_dict['colour'], marker_size=0,
                    #              subplot=syst_truth_hist_bin),
                ])

            # add nominal ones last
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

            syst_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Nominal unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(syst_entries, "%s bin %d" % (func_name, ibin)):
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
            self.save_plot(plot, "%s/unfolded_%s_syst_scale_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_scale_systs_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            syst_entries = []
            hbc_args = dict(ind=ibin, binning_scheme='generator')

            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('hist_truth', **hbc_args)

            for syst_dict in self.region['scale_systematics']:
                syst_unfolder = syst_dict['unfolder']
                syst_label = syst_dict['label']

                # Get binned hists from the scale unfolder, since the error bars may have been setup specially
                syst_unfolded_hist_bin = syst_unfolder.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)

                syst_entries.extend([
                    Contribution(syst_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (total unc.) (%s)" % (syst_unfolder.tau, syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=mc_gen_hist_bin),  # use this as ref as scale variations are only different response matrices
                ])

            # add nominal ones last
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)

            syst_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Nominal unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(syst_entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(syst_entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
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
            self.save_plot(plot, "%s/unfolded_%s_syst_scale_unnormalised_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_model_systs_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            syst_entries = []
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            for syst_dict in self.region['model_systematics']:
                syst_unfolder = syst_dict['unfolder']
                syst_label = syst_dict['label']

                # Get binned hists from the model unfolder, since the error bars may have been setup specially
                syst_unfolded_hist_bin = syst_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
                syst_gen_hist_bin = syst_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)

                syst_entries.extend([
                    Contribution(syst_gen_hist_bin,
                                 label="Generator (%s)" % (syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                                 marker_color=syst_dict['colour'], marker_size=0),
                    Contribution(syst_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (%s)" % (syst_unfolder.tau, syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=syst_gen_hist_bin),
                ])

            # add nominal ones last
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

            syst_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(syst_entries, "%s bin %d" % (func_name, ibin)):
                return
            plot_args = copy(self.pt_bin_plot_args)
            plot_args['subplot_title'] = "Unfolded / Gen"
            plot = Plot(syst_entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **plot_args)
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
            self.save_plot(plot, "%s/unfolded_%s_syst_model_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_model_systs_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            syst_entries = []
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('hist_truth', **hbc_args)
            for syst_dict in self.region['model_systematics']:
                syst_unfolder = syst_dict['unfolder']
                syst_label = syst_dict['label']

                # Get binned hists from the model unfolder, since the error bars may have been setup specially
                syst_unfolded_hist_bin = syst_unfolder.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)
                syst_gen_hist_bin = syst_unfolder.hist_bin_chopper.get_pt_bin_div_bin_width('hist_truth', **hbc_args)

                syst_entries.extend([
                    Contribution(syst_gen_hist_bin,
                                 label="Generator (%s)" % (syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                                 marker_color=syst_dict['colour'], marker_size=0),
                    Contribution(syst_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (%s)" % (syst_unfolder.tau, syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=syst_gen_hist_bin),
                ])

            # add nominal ones last
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)

            syst_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(syst_entries, "%s bin %d" % (func_name, ibin)):
                return
            plot_args = copy(self.pt_bin_plot_args)
            plot_args['subplot_title'] = "Unfolded / Gen"
            plot = Plot(syst_entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **plot_args)
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
            self.save_plot(plot, "%s/unfolded_%s_syst_model_unnormalised_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_pdf_systs_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            pdf_entries = []
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            for pdf_dict in self.region['pdf_systematics']:
                pdf_unfolder = pdf_dict['unfolder']
                pdf_label = pdf_dict['label']
                pdf_label_no_spaces = cu.no_space_str(pdf_dict['label'])

                self.hist_bin_chopper.add_obj('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), pdf_unfolder.unfolded)
                # self.hist_bin_chopper.add_obj('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), pdf_unfolder.hist_truth)

                pdf_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), **hbc_args)
                # pdf_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), **hbc_args)

                pdf_entries.extend([
                    Contribution(pdf_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (stat. unc.) (%s)" % (pdf_unfolder.tau, pdf_label),
                                 line_color=pdf_dict['colour'], line_width=1, line_style=1,
                                 marker_color=pdf_dict['colour'], marker_size=0,
                                 subplot=mc_gen_hist_bin),
                    # Disable this as the PDF variations are in the response matrix, not input being unfolded
                    # so gen wiill always be the same!
                    # Contribution(pdf_gen_hist_bin,
                    #              label="Generator (%s)" % (pdf_label),
                    #              line_color=pdf_dict['colour'], line_width=1, line_style=2,
                    #              marker_color=pdf_dict['colour'], marker_size=0),
                ])

            # add nominal ones last
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

            pdf_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(pdf_entries, "%s bin %d" % (func_name, ibin)):
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
            self.save_plot(plot, "%s/unfolded_%s_pdf_model_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
            ROOT.gStyle.SetPalette(ROOT.kViridis)

    def plot_unfolded_with_pdf_systs_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            pdf_entries = []
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('hist_truth', **hbc_args)
            for pdf_dict in self.region['pdf_systematics']:
                pdf_unfolder = pdf_dict['unfolder']
                pdf_label = pdf_dict['label']
                pdf_label_no_spaces = cu.no_space_str(pdf_dict['label'])

                self.hist_bin_chopper.add_obj('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), pdf_unfolder.unfolded)
                # self.hist_bin_chopper.add_obj('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), pdf_unfolder.hist_truth)

                pdf_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), **hbc_args)
                # pdf_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), **hbc_args)

                pdf_entries.extend([
                    Contribution(pdf_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (stat. unc.) (%s)" % (pdf_unfolder.tau, pdf_label),
                                 line_color=pdf_dict['colour'], line_width=1, line_style=1,
                                 marker_color=pdf_dict['colour'], marker_size=0,
                                 subplot=mc_gen_hist_bin),
                    # Contribution(pdf_gen_hist_bin,
                    #              label="Generator (%s)" % (pdf_label),
                    #              line_color=pdf_dict['colour'], line_width=1, line_style=2,
                    #              marker_color=pdf_dict['colour'], marker_size=0),
                ])

            # add nominal ones last
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)

            pdf_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(pdf_entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(pdf_entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
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
            self.save_plot(plot, "%s/unfolded_%s_pdf_model_unnormalised_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
            ROOT.gStyle.SetPalette(ROOT.kViridis)

    def plot_unfolded_with_jackknife_input_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            entries = []
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            for ind, jk_dict in enumerate(self.region['jackknife_input_variations']):
                jk_unfolder = jk_dict['unfolder']
                jk_label = jk_dict['label']

                # Get binned hists from the jackknife unfolder, since the error bars may have been setup specially
                jk_unfolded_hist_bin = jk_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
                jk_gen_hist_bin = jk_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)

                entries.extend([
                    Contribution(jk_unfolded_hist_bin,
                                 label="Unfolded (%s)" % (jk_label),
                                 line_color=jk_dict['colour'], line_width=self.line_width,
                                 line_style=1 if ind % 2 == 0 else 3,
                                 marker_color=jk_dict['colour'], marker_size=0,
                                 subplot=jk_gen_hist_bin),
                    Contribution(jk_gen_hist_bin,
                                 label="Generator (%s)" % (jk_label),
                                 line_color=jk_dict['colour'], line_width=1, line_style=2,
                                 marker_color=jk_dict['colour'], marker_size=0),
                ])

            # add nominal ones last
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

            entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            if len(entries) > 4:
                # plot.legend.SetX1(0.53)
                plot.legend.SetY1(0.65)
                plot.y_padding_max_linear = 1.8
            if len(entries) > 6:
                plot.legend.SetNColumns(2)
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_%s_jackknife_input_vars_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_jackknife_response_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            entries = []
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            for ind, jk_dict in enumerate(self.region['jackknife_response_variations']):
                jk_unfolder = jk_dict['unfolder']
                jk_label = jk_dict['label']

                # Get binned hists from the jackknife unfolder, since the error bars may have been setup specially
                jk_unfolded_hist_bin = jk_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

                entries.extend([
                    Contribution(jk_unfolded_hist_bin,
                                 label="Unfolded (%s)" % (jk_label),
                                 line_color=jk_dict['colour'], line_width=self.line_width,
                                 line_style=1 if ind % 2 == 0 else 3,
                                 marker_color=jk_dict['colour'], marker_size=0,
                                 subplot=mc_gen_hist_bin),
                ])

            # add nominal ones last
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

            entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            if len(entries) > 4:
                # plot.legend.SetX1(0.53)
                plot.legend.SetY1(0.65)
                plot.y_padding_max_linear = 1.8
            if len(entries) > 6:
                plot.legend.SetNColumns(2)
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/unfolded_%s_jackknife_response_vars_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_jackknife_residuals(self, jackknife_variations):
        """Plot distributions of unfolded / gen for each bin for all jackknife variations

        Hopefully Gaussian with centre on 1
        """
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')

            # hists of all jackknife variations for all bins
            hists = [ROOT.TH1D("jk_ratios_%d" % ind, "%s, %s, bin %d;Unfolded/Gen;N" % (self.setup.angle_str, self.get_pt_bin_title(bin_edge_low, bin_edge_high).replace("\n", ", "), ind), 25, 0.6, 1.4)
                     for ind in range(len(self.unfolder.binning_handler.get_variable_bins(bin_edge_low, 'generator'))-1)]

            for ind, jk_dict in enumerate(jackknife_variations):
                jk_unfolder = jk_dict['unfolder']
                jk_label = jk_dict['label']
                jk_gen_hist_bin = jk_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
                jk_unfolded_hist_bin = jk_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
                ratio = jk_unfolded_hist_bin.Clone()
                ratio.Divide(jk_gen_hist_bin)
                for ix in range(1, ratio.GetNbinsX()+1):
                    hists[ix-1].Fill(ratio.GetBinContent(ix))

            entries = [
                Contribution(h,
                             label="bin %d" % ind,
                             marker_color=cu.get_colour_seq(ind, len(hists)),
                             line_style=1 + (ind % 2),
                             line_width=2,
                             line_color=cu.get_colour_seq(ind, len(hists)))
                for ind, h in enumerate(hists)
            ]
            canv = ROOT.TCanvas(cu.get_unique_str(), "", 1200, 900)
            n = len(hists)
            nx = math.ceil(n / 3)
            ny = 3
            canv.Divide(nx, ny)
            for i, h in enumerate(hists, 1):
                canv.cd(i)
                h.Draw("HIST")
                h.SetMaximum(h.GetMaximum()*1.2)
            canv.SaveAs("%s/jackknife_residuals_%s_bin_%d.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
            # plot = Plot(entries[1:-1], what='hist',
            #             title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
            #             ylim=(0, None),
            #             xtitle="Unfolded / Gen")
            # plot.plot("NOSTACK HIST")
            # plot.legend.SetNColumns(2)
            # self.save_plot(plot, "%s/jackknife_residuals_%s_bin_%d.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_uncertainty_shifts_normalised(self):
        """Do plots of fractional uncertainty shifts on *normalised* unfolded distribution"""
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            entries = []
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.unfolder.rsp_uncert_name, **hbc_args)

            for syst_dict in self.region['experimental_systematics']:
                # For each systematic, get the normalised shift and hence fraction
                this_syst = self.unfolder.get_exp_syst(syst_dict['label'])
                syst_unfolded_fraction = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(this_syst.syst_shift_label, **hbc_args).Clone()
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
                                 label="Input stat.",
                                 line_color=ROOT.kRed,
                                 line_style=3,
                                 line_width=3,
                                 marker_size=0,
                                 marker_color=ROOT.kRed,
                                 )
            c_syst = Contribution(h_syst,
                                 label="Response matrix stat.",
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

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
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
            self.save_plot(plot, output_filename)

            plot.y_padding_max_log = 50
            plot.set_logy(do_more_labels=False)
            plot.get_modifier().SetMinimum(1E-4)
            log_filename, ext = os.path.splitext(output_filename)
            self.save_plot(plot, log_filename+"_log"+ext)

    def plot_unfolded_with_exp_systs_normalised(self):
        """Plot shifted unfolded normalised distributions for each syst"""
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.unfolder.rsp_uncert_name, **hbc_args)

            def _remove_error_bars(h):
                for i in range(1, h.GetNbinsX()+1):
                    h.SetBinError(i, 0)

            entries = [
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75),
            ]

            for syst_dict in self.region['experimental_systematics']:
                this_syst = self.unfolder.get_exp_syst(syst_dict['label'])
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(this_syst.syst_shifted_label, **hbc_args)
                _remove_error_bars(syst_unfolded_hist_bin)
                default_linestyle = 2 if 'down' in syst_dict['label'].lower() else 1
                c = Contribution(syst_unfolded_hist_bin,
                                 label=syst_dict['label'],
                                 line_color=syst_dict['colour'], line_width=self.line_width,
                                 line_style=syst_dict.get('linestyle', default_linestyle),
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=unfolded_hist_bin_stat_errors)
                entries.append(c)


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            xlim = qgp.calc_auto_xlim(entries[0:1])  # only for data
            plot = Plot(entries,
                        xtitle=self.setup.particle_title,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        xlim=xlim,
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
            self.save_plot(plot, "%s/unfolded_syst_variations_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_exp_systs_unnormalised(self):
        """Plot shifted unfolded absolute distributions for each syst"""
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', **hbc_args)
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.unfolder.rsp_uncert_name, **hbc_args)

            def _remove_error_bars(h):
                for i in range(1, h.GetNbinsX()+1):
                    h.SetBinError(i, 0)

            entries = [
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75),
            ]

            for syst_dict in self.region['experimental_systematics']:
                this_syst = self.unfolder.get_exp_syst(syst_dict['label'])
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width(this_syst.syst_shifted_label, **hbc_args)
                _remove_error_bars(syst_unfolded_hist_bin)
                c = Contribution(syst_unfolded_hist_bin,
                                 label=syst_dict['label'],
                                 line_color=syst_dict['colour'], line_width=self.line_width,
                                 line_style=2 if 'down' in syst_dict['label'].lower() else 1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=unfolded_hist_bin_stat_errors)
                entries.append(c)


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            xlim = qgp.calc_auto_xlim(entries[0:1])  # only for data
            plot = Plot(entries,
                        xtitle=self.setup.particle_title,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        xlim=xlim,
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
            self.save_plot(plot, "%s/unfolded_syst_variations_unnormalised_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_syst_fraction_normalised(self):
        """Plot varation / central value on normalised hists
        (basically the subplot from plot_unfolded_with_exp_systs_normalised)
        """
        def _convert_error_bars_to_error_ratio_hist(h, direction=1):
            """Create hist with bin content = (bin value ± bin error) / bin value, 0 error"""
            h_new = h.Clone(h.GetName() + cu.get_unique_str())
            for ix in range(1, h_new.GetNbinsX()+1):
                if h.GetBinContent(ix) > 0:
                    h_new.SetBinContent(ix, 1+(direction*(h.GetBinError(ix) / h.GetBinContent(ix))))
                else:
                    if h.GetBinContent(ix) < 0:
                        h_new.SetBinContent(ix, 1+(direction*(h.GetBinError(ix) / h.GetBinContent(ix))))
                        print("_convert_error_bars_to_error_ratio_hist() warning: bin %d content < 0!" % (ix))
                    else:
                        h_new.SetBinContent(ix, 1)
                h_new.SetBinError(ix, 0)
            return h_new

        def _convert_syst_shift_to_error_ratio_hist(h_syst, h_nominal):
            """Create h_syst / h_nominal without error bars"""
            h_new = h_syst.Clone(h_syst.GetName() + cu.get_unique_str())
            h_new.Divide(h_nominal)
            for ix in range(1, h_new.GetNbinsX()+1):
                if h_nominal.GetBinContent(ix) == 0:
                    h_new.SetBinContent(ix, 1)
                if h_syst.GetBinContent(ix) < 0:
                    print("Warning: _convert_syst_shift_to_error_ratio_hist bin", ix, "of h_syst < 0", h_syst.GetName())
                h_new.SetBinError(ix, 0)
            return h_new

        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.unfolder.rsp_uncert_name, **hbc_args)

            unfolded_hist_bin_no_errors = unfolded_hist_bin_total_errors.Clone()
            cu.remove_th1_errors(unfolded_hist_bin_no_errors)

            input_stats = unfolded_hist_bin_stat_errors.Clone()
            input_stats.Divide(unfolded_hist_bin_no_errors)

            total_err = unfolded_hist_bin_total_errors.Clone()
            total_err.Divide(unfolded_hist_bin_no_errors)

            entries = [
                # TOTAL UNCERT
                Contribution(total_err,
                             label="Total uncertainty",
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=0, line_style=2,
                             marker_color=self.plot_styles['unfolded_total_colour'], marker_style=cu.Marker.get('circle'), marker_size=0,
                             fill_style=3354,
                             # fill_style=3005,
                             fill_color=16),
                # INPUT UNCERT
                Contribution(input_stats,
                             label="Data stat.",
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=0, line_style=3,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0,
                             fill_style=3545,
                             # fill_style=3003,
                             fill_color=self.plot_styles['unfolded_stat_colour']),
            ]

            # Add experimental systs
            for syst_dict, mark in zip(self.region['experimental_systematics'], cu.Marker().cycle(cycle_filling=True)):
                this_syst = self.unfolder.get_exp_syst(syst_dict['label'])
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(this_syst.syst_shifted_label, **hbc_args)
                this_syst_hist = _convert_syst_shift_to_error_ratio_hist(syst_unfolded_hist_bin, unfolded_hist_bin_total_errors)
                is_herwig = "shower" in syst_dict['label'].lower() or "herwig" in syst_dict['label'].lower()
                c = Contribution(this_syst_hist,
                                 label="" if "down" in syst_dict['label'] else syst_dict['label'].replace(" up","").replace("hadronization","hadr."),
                                 line_color=syst_dict['colour'],
                                 leg_draw_opt="L" if is_herwig else "P",
                                 line_width=self.line_width if is_herwig else 0,
                                 # line_width=self.line_width,
                                 line_style=5,
                                 marker_color=syst_dict['colour'],
                                 marker_size=0 if is_herwig else 1.25,
                                 marker_style=mark)
                entries.append(c)

            # Add scale syst
            if self.unfolder.scale_uncert_name in self.hist_bin_chopper.objects:
                scale_hist = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.unfolder.scale_uncert_name, **hbc_args)
                scale_col = self.plot_styles['scale_colour']
                scale_style = dict(line_color=scale_col, line_width=self.line_width, line_style=2,
                                   marker_color=scale_col, marker_style=cu.Marker.get('circle'), marker_size=0,
                                   fill_style=0, fill_color=15)
                entries.extend([
                    Contribution(_convert_error_bars_to_error_ratio_hist(scale_hist),
                                label='#mu scales', leg_draw_opt="L",
                                **scale_style),
                    # add -ve side, no label as we don't want it in legend
                    Contribution(_convert_error_bars_to_error_ratio_hist(scale_hist, -1),
                                **scale_style)
                ])

            # Add pdf syst
            if self.unfolder.pdf_uncert_name in self.hist_bin_chopper.objects:
                pdf_hist = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.unfolder.pdf_uncert_name, **hbc_args)
                pdf_col = self.plot_styles['pdf_colour']
                pdf_style = dict(line_color=pdf_col, line_width=self.line_width, line_style=6,
                                 marker_color=pdf_col, marker_style=cu.Marker.get('circle'), marker_size=0,
                                 fill_style=0, fill_color=15)
                entries.extend([
                    Contribution(_convert_error_bars_to_error_ratio_hist(pdf_hist),
                                label='PDF', leg_draw_opt="L",
                                **pdf_style),
                    # add -ve side, no label as we don't want it in legend
                    Contribution(_convert_error_bars_to_error_ratio_hist(pdf_hist, -1),
                                **pdf_style)
                ])


            rsp_col = self.plot_styles['rsp_colour']
            rsp_style = dict(line_color=rsp_col, line_width=self.line_width, line_style=3,
                             marker_color=rsp_col, marker_style=cu.Marker.get('circle'), marker_size=0,
                             fill_style=0, fill_color=rsp_col)
            entries.extend([
                # RESPONSE UNCERT
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_rsp_errors),
                             label="Response m. stat.", leg_draw_opt="L",
                             **rsp_style),
                # Add in the -ve side, but no label as we don't want it in the legend
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_rsp_errors, -1),
                             **rsp_style),
            ])


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            # find upper limit where there is no data
            xlim = qgp.calc_auto_xlim([unfolded_hist_bin_total_errors])
            ylim = [0.8, 1.45] if "Dijet" in self.setup.region['name'] else [0.3, 1.9]
            # sepcial one for paper
            ylim = [0.85, 1.3] if ("Dijet" in self.setup.region['name'] and "LHA" in self.setup.angle.var and "AK4" in self.setup.jet_algo.upper()) else ylim
            min_total = _convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_total_errors, -1).GetMinimum()
            max_total = _convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_total_errors).GetMaximum()
            if max_total > ylim[1]:
                ylim[1] = max_total*1.1
            if min_total < ylim[0]:
                ylim[0] = min_total*0.9
            ylim[0] = max(-1, ylim[0])
            ylim[1] = min(5, ylim[1])
            plot = Plot(entries,
                        xtitle=self.setup.particle_title,
                        ytitle='Variation / nominal (on normalised distribution)',
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        xlim=xlim,
                        ylim=ylim,
                        subplot_type=None)
            self._modify_plot_paper(plot)
            #plot.default_canvas_size = (800, 700)
            plot.legend.SetX1(0.48)
            plot.legend.SetY1(0.56)
            plot.legend.SetX2(0.93)
            plot.legend.SetY2(0.87)
            #if len(entries) > 10: plot.legend.SetNColumns(2)
            final_reading_factor=1.15
            plot.legend.SetTextSize(0.035*final_reading_factor)
            plot.plot("NOSTACK E2 P L") # hard to get one that is points for systs, and line for stats, and fill for shaded
            plot.main_pad.cd()
            plot.container.GetXaxis().SetTickSize(0.02)
            if xlim is not None:
                line = ROOT.TLine(xlim[0], 1, xlim[1], 1)
            else:
                xmin = unfolded_hist_bin_stat_errors.GetXaxis().GetBinLowEdge(1)
                xmax = unfolded_hist_bin_stat_errors.GetXaxis().GetBinLowEdge(unfolded_hist_bin_stat_errors.GetNbinsX()+1)
                line = ROOT.TLine(xmin, 1, xmax, 1)
            line.SetLineColor(ROOT.kGray+2)
            line.SetLineColor(ROOT.kBlack)
            line.Draw()
            l2=plot.legend.Clone("l2")
            l2.Clear()
            for e in entries:
              if e.label!="" and e.marker_size>0: continue
              if e.label!="" and e.line_width>0 and not "Shower" in e.label: continue
              l2.AddEntry(e.obj,"","P" if e.marker_size>0 else "")
            l2.SetX1(0.51)
            l2.Draw("SAME")
            stp = self.setup
            fname = f'unfolded_syst_variations_vs_nominal_{stp.append}_bin_{ibin:d}_divBinWidth{stp.paper_str}.{stp.output_fmt}'
            self.save_plot(plot, os.path.join(self.setup.output_dir, fname))
            print(os.path.join(self.setup.output_dir, fname))

    def plot_syst_fraction_unnormalised(self):
        """Plot varation / central value on absolute hists
        (basically the subplot from plot_unfolded_with_exp_systs_unnormalised)
        """
        def _convert_error_bars_to_error_ratio_hist(h, direction=1):
            """Create hist with bin content = error shift / bin value, 0 error"""
            h_new = h.Clone(h.GetName() + cu.get_unique_str())
            for ix in range(1, h_new.GetNbinsX()+1):
                if h.GetBinContent(ix) > 0:
                    h_new.SetBinContent(ix, 1+(direction*(h.GetBinError(ix) / h.GetBinContent(ix))))
                else:
                    if h.GetBinContent(ix) < 0:
                        h_new.SetBinContent(ix, 1+(direction*(h.GetBinError(ix) / h.GetBinContent(ix))))
                        print("_convert_error_bars_to_error_ratio_hist() warning: bin %d content < 0!" % (ix))
                    else:
                        h_new.SetBinContent(ix, 1)
                h_new.SetBinError(ix, 0)
            return h_new

        def _convert_syst_shift_to_error_ratio_hist(h_syst, h_nominal):
            """Create h_syst / h_nominal without error bars"""
            h_new = h_syst.Clone(h_syst.GetName() + cu.get_unique_str())
            h_new.Divide(h_nominal)
            for ix in range(1, h_new.GetNbinsX()+1):
                h_new.SetBinError(ix, 0)
            return h_new

        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', **hbc_args)
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.unfolder.rsp_uncert_name, **hbc_args)

            unfolded_hist_bin_no_errors = unfolded_hist_bin_total_errors.Clone()
            cu.remove_th1_errors(unfolded_hist_bin_no_errors)

            input_stats = unfolded_hist_bin_stat_errors.Clone()
            input_stats.Divide(unfolded_hist_bin_no_errors)

            total_err = unfolded_hist_bin_total_errors.Clone()
            total_err.Divide(unfolded_hist_bin_no_errors)

            entries = [
                # TOTAL UNCERT
                Contribution(total_err,
                             label="Total uncertainty",
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=0, line_style=2,
                             marker_color=self.plot_styles['unfolded_total_colour'], marker_style=cu.Marker.get('circle'), marker_size=0,
                             fill_style=3354,
                             # fill_style=3005,
                             fill_color=16),
                # INPUT UNCERT
                Contribution(input_stats,
                             label="Data stat.",
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=0, line_style=3,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0,
                             fill_style=3545,
                             # fill_style=3003,
                             fill_color=self.plot_styles['unfolded_stat_colour']),
            ]

            # Add experimental systs
            for syst_dict, mark in zip(self.region['experimental_systematics'], cu.Marker().cycle(cycle_filling=True)):
                this_syst = self.unfolder.get_exp_syst(syst_dict['label'])
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width(this_syst.syst_shifted_label, **hbc_args)
                this_syst_hist = _convert_syst_shift_to_error_ratio_hist(syst_unfolded_hist_bin, unfolded_hist_bin_total_errors)
                is_herwig = "shower" in syst_dict['label'].lower() or "herwig" in syst_dict['label'].lower()
                c = Contribution(this_syst_hist,
                                 label="" if "down" in syst_dict['label'] else syst_dict['label'].replace(" up","").replace("hadronization","hadr."),
                                 line_color=syst_dict['colour'],
                                 leg_draw_opt="L" if is_herwig else "P",
                                 line_width=self.line_width if is_herwig else 0,
                                 # line_width=self.line_width,
                                 line_style=5,
                                 marker_color=syst_dict['colour'],
                                 marker_size=0 if is_herwig else 1.25,
                                 marker_style=mark)
                entries.append(c)

            # Add scale syst
            if self.unfolder.scale_uncert_name in self.hist_bin_chopper.objects:
                scale_hist = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.unfolder.scale_uncert_name, **hbc_args)
                scale_col = self.plot_styles['scale_colour']
                scale_style = dict(line_color=scale_col, line_width=self.line_width, line_style=2,
                                   marker_color=scale_col, marker_style=cu.Marker.get('circle'), marker_size=0,
                                   fill_style=0, fill_color=15)
                entries.extend([
                    Contribution(_convert_error_bars_to_error_ratio_hist(scale_hist),
                                label='#mu scales', leg_draw_opt="L",
                                **scale_style),
                    # add -ve side, no label as we don't want it in legend
                    Contribution(_convert_error_bars_to_error_ratio_hist(scale_hist, -1),
                                **scale_style)
                ])

            # Add pdf syst
            if self.unfolder.pdf_uncert_name in self.hist_bin_chopper.objects:
                pdf_hist = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.unfolder.pdf_uncert_name, **hbc_args)
                pdf_col = self.plot_styles['pdf_colour']
                pdf_style = dict(line_color=pdf_col, line_width=self.line_width, line_style=6,
                                 marker_color=pdf_col, marker_style=cu.Marker.get('circle'), marker_size=0,
                                 fill_style=0, fill_color=15)
                entries.extend([
                    Contribution(_convert_error_bars_to_error_ratio_hist(pdf_hist),
                                label='PDF', leg_draw_opt="L",
                                **pdf_style),
                    # add -ve side, no label as we don't want it in legend
                    Contribution(_convert_error_bars_to_error_ratio_hist(pdf_hist, -1),
                                **pdf_style)
                ])

            rsp_col = self.plot_styles['rsp_colour']
            rsp_style = dict(line_color=rsp_col, line_width=self.line_width, line_style=3,
                             marker_color=rsp_col, marker_style=cu.Marker.get('circle'), marker_size=0,
                             fill_style=0, fill_color=rsp_col)
            entries.extend([
                # RESPONSE UNCERT
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_rsp_errors),
                             label="Response m. stat.", leg_draw_opt="L",
                             **rsp_style),
                # Add in the -ve side, but no label as we don't want it in the legend
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_rsp_errors, -1),
                             **rsp_style),
            ])


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            # find upper limit where there is no data
            xlim = qgp.calc_auto_xlim([unfolded_hist_bin_total_errors])
            ylim = [0.8, 1.45] if "Dijet" in self.setup.region['name'] else [0.3, 1.9]
            # special one for paper
            ylim = [0.85, 1.3] if ("Dijet" in self.setup.region['name'] and "LHA" in self.setup.angle.var and "AK4" in self.setup.jet_algo.upper()) else ylim
            min_total = _convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_total_errors, -1).GetMinimum()
            max_total = _convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_total_errors).GetMaximum()
            if max_total > ylim[1]:
                ylim[1] = max_total*1.1
            if min_total < ylim[0]:
                ylim[0] = min_total*0.9
            plot = Plot(entries,
                        xtitle=self.setup.particle_title,
                        ytitle='Variation / nominal (on absolute distribution)',
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        xlim=xlim,
                        ylim=ylim,
                        subplot_type=None)
            self._modify_plot_paper(plot)
            #plot.default_canvas_size = (800, 700)
            plot.legend.SetX1(0.48)
            plot.legend.SetY1(0.56)
            plot.legend.SetX2(0.93)
            plot.legend.SetY2(0.87)
            #if len(entries) > 10: plot.legend.SetNColumns(2)
            final_reading_factor=1.15
            plot.legend.SetTextSize(0.035*final_reading_factor)
            plot.plot("NOSTACK E2 P L") # hard to get one that is points for systs, and line for stats
            plot.main_pad.cd()
            plot.container.GetXaxis().SetTickSize(0.02)
            if xlim is not None:
                line = ROOT.TLine(xlim[0], 1, xlim[1], 1)
            else:
                xmin = unfolded_hist_bin_stat_errors.GetXaxis().GetBinLowEdge(1)
                xmax = unfolded_hist_bin_stat_errors.GetXaxis().GetBinLowEdge(unfolded_hist_bin_stat_errors.GetNbinsX()+1)
                line = ROOT.TLine(xmin, 1, xmax, 1)
            line.SetLineColor(ROOT.kGray+2)
            line.SetLineColor(ROOT.kBlack)
            line.Draw()
            l2=plot.legend.Clone("l2")
            l2.Clear()
            for e in entries:
              if e.label!="" and e.marker_size>0: continue
              if e.label!="" and e.line_width>0 and not "Shower" in e.label: continue
              l2.AddEntry(e.obj,"","P" if e.marker_size>0 else "")
            l2.SetX1(0.51)
            l2.Draw("SAME")
            stp = self.setup
            fname = f'unfolded_unnormalised_syst_variations_vs_nominal_{stp.append}_bin_{ibin:d}_divBinWidth{stp.paper_str}.{stp.output_fmt}'
            self.save_plot(plot, os.path.join(self.setup.output_dir, fname))

    def plot_detector_normalised_bg_subtracted(self, alt_detector=None):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_gen_binning_bg_subtracted', **hbc_args)
            mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_gen_binning_bg_subtracted', **hbc_args)

            entries = [
                Contribution(mc_hist_bin,
                             label="MC (bg-subtracted) [%s]" % (self.setup.region['mc_label']) if alt_detector else "MC (bg-subtracted)",
                             line_color=self.plot_styles['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_mc_colour'], marker_size=0,
                             subplot=input_hist_bin if alt_detector else None),
            ]
            if alt_detector:
                self.hist_bin_chopper.add_obj("alt_hist_mc_reco_gen_binning_bg_subtracted", alt_detector)
                alt_mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_mc_reco_gen_binning_bg_subtracted', **hbc_args)
                entries.append(Contribution(alt_mc_hist_bin,
                                            label="MC (bg-subtracted) [%s]" % self.setup.region['alt_mc_label'],
                                            line_color=self.plot_styles['alt_reco_colour'], line_width=self.line_width,
                                            marker_color=self.plot_styles['alt_reco_colour'], marker_size=0,
                                            subplot=input_hist_bin))

            entries.append(
                Contribution(input_hist_bin,
                             label="Data (bg-subtracted)",
                             line_color=self.plot_styles['reco_data_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_data_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=None if alt_detector else mc_hist_bin)
            )

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
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
            self.save_plot(plot, "%s/detector_gen_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_ematrix(self, h2d, variable_bins, title, output_filename):
        """Generic function to plot error matrix for a given pt bin"""
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 700, 600)
        canv.SetTicks(1, 1)
        canv.SetRightMargin(0.15)
        canv.SetTopMargin(0.18)
        h2d.SetTitle("%s;%s;%s" % (title, self.setup.particle_title, self.setup.particle_title))
        h2d.SetTitleSize(0.015, "t")  # the "t" is anything that isn't x, y, z
        # Setup bin labels
        for i, (val_low, val_high) in enumerate(zip(variable_bins[:-1], variable_bins[1:]), 1):
            h2d.GetXaxis().SetBinLabel(i, "%g-%g" % (val_low, val_high))
            h2d.GetYaxis().SetBinLabel(i, "%g-%g" % (val_low, val_high))
        h2d.SetContour(256)
        # enforce scientific notation for now
        default_paint_text_format = ROOT.gStyle.GetPaintTextFormat()
        ROOT.gStyle.SetPaintTextFormat(".2e")
        h2d.Draw("COLZ TEXT")
        h2d.GetYaxis().SetTitleOffset(h2d.GetYaxis().GetTitleOffset()*1.1)
        # Symmetrize z axis and use French flag colours
        cu.symmetrize_h2d_z_limits(h2d)
        cu.set_french_flag_palette()
        canv.SaveAs(output_filename)
        ROOT.gStyle.SetPalette(PLOT_STYLES['default_palette'])
        ROOT.gStyle.SetPaintTextFormat(default_paint_text_format)

    def plot_total_ematrix(self):
        """Plot total normalised error matrix for each gen pt bin"""
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            ematrix = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.unfolder.total_ematrix_name, ibin, binning_scheme='generator')
            title = "#splitline{Total normalised error matrix for %g < %s < %g GeV}{%s region, %s}" % (bin_edge_low, self.setup.pt_var_str, bin_edge_high, self.region['label'], self.setup.jet_algo)
            output_filename = os.path.join(self.setup.output_dir,
                                           "unfolded_ematrix_total_%s_bin_%d_norm_divBinWidth.%s" % (self.setup.append, ibin, self.setup.output_fmt))
            self.plot_ematrix(ematrix,
                              variable_bins=self.unfolder.binning_handler.get_variable_bins(bin_edge_low, binning_scheme='generator'),
                              title=title,
                              output_filename=output_filename)

# ================================================================================


class GenLambdaBinnedPlotter(BinnedPlotter):
    def __init__(self, setup, bins, hist_bin_chopper, unfolder):
        self.setup = setup
        self.region = setup.region
        self.bins = bins
        self.hist_bin_chopper = hist_bin_chopper

        self.line_width = 2
        self.plot_styles = PLOT_STYLES
        self.lambda_bin_plot_args = dict(
            what="hist",
            xtitle=self.setup.pt_str,
            has_data=self.setup.has_data,
            subplot_type='ratio',
            subplot_title="Unfolded / Gen",
            subplot_limits=(0, 2) if self.setup.has_data else (0.75, 1.25),
        )
        self.unfolder = unfolder
        super().__init__()

    def _modify_plot(self, this_plot):
        if self.setup.output_fmt != "pdf":
            _make_thumbnail_canvas(this_plot)
        this_plot.legend.SetX1(0.6)
        this_plot.legend.SetY1(0.68)
        this_plot.legend.SetX2(0.98)
        this_plot.legend.SetY2(0.9)
        this_plot.left_margin = 0.16
        this_plot.y_padding_max_log = 5000 # space for title
        this_plot.lumi = self.setup.lumi

    def get_lambda_bin_title(self, bin_edge_low, bin_edge_high):
        title = (("{jet_algo}\n"
                  "{region_label} region\n"
                  "{bin_edge_low:g} < {angle_str} < {bin_edge_high:g}")
                 .format(
                    jet_algo=self.setup.jet_algo,
                    region_label=self.region['label'],
                    angle_str=self.setup.angle_str,
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
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'],# marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_unreg_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            # unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')
            unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unreg_unfolded', ibin, binning_scheme='generator')

            # unnormalised version
            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'],# marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unreg_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = 0) (total unc.)",
                             line_color=self.plot_styles['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_unreg_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_with_unreg_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_template_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', **hbc_args)
            alt_mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('alt_hist_truth', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', **hbc_args)
            unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unreg_unfolded', **hbc_args)
            template_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('truth_template', **hbc_args)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(alt_mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['alt_mc_label']),
                             line_color=self.plot_styles['alt_gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['alt_gen_colour'], marker_size=0,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unreg_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = 0) (total unc.)",
                             line_color=self.plot_styles['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_unreg_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(template_hist_bin_total_errors,
                             label="Template",
                             line_color=self.plot_styles['template_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['template_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),

            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_with_template_unreg_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_alt_response_unnormalised(self, alt_unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin, binning_scheme='generator')
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            alt_unfolded_hist_bin_total_errors = alt_unfolder.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')
            # alt_mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('alt_hist_truth', ibin, binning_scheme='generator')

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                # Contribution(alt_mc_gen_hist_bin,
                #              label="Generator (%s)" % (self.region['alt_mc_label']),
                #              line_color=alt_gen_colour, line_width=self.line_width, line_style=2,
                #              marker_color=alt_gen_colour, marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (self.unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
                Contribution(alt_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)\n(%s response matrix)" % (alt_unfolder.tau, self.region['alt_mc_label']),
                             line_color=self.plot_styles['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['alt_unfolded_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.lambda_bin_unnormalised_differential_label,
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_alt_response_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_scale_systs_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            syst_entries = []

            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin, binning_scheme='generator')

            for syst_dict in self.region['scale_systematics']:
                syst_unfolder = syst_dict['unfolder']
                syst_label = syst_dict['label']

                syst_unfolded_hist_bin = syst_unfolder.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')

                syst_entries.extend([
                    Contribution(syst_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (total unc.) (%s)" % (syst_unfolder.tau, syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=mc_gen_hist_bin),
                ])

            # add nominal ones last
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')

            syst_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Nominal unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(syst_entries, "%s bin %d" % (func_name, ibin)):
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
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_syst_scale_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_model_systs_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            syst_entries = []
            for syst_dict in self.region['model_systematics']:
                syst_unfolder = syst_dict['unfolder']
                syst_label = syst_dict['label']

                syst_unfolded_hist_bin = syst_unfolder.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')
                syst_gen_hist_bin = syst_unfolder.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin, binning_scheme='generator')

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
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(syst_entries, "%s bin %d" % (func_name, ibin)):
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
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_syst_model_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

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
                             line_color=self.plot_styles['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['gen_colour'], marker_size=0),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], #marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_gen_hist_bin),
            ])

            func_name = cu.get_current_func_name()
            if not self.check_entries(pdf_entries, "%s bin %d" % (func_name, ibin)):
                return
            plot = Plot(pdf_entries,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
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
            self.save_plot(plot, "%s/unfolded_unnormalised_%s_pdf_model_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
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
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width(self.unfolder.rsp_uncert_name, ibin, binning_scheme='generator')

            for syst_dict in self.region['experimental_systematics']:
                # For each systematic, get the normalised shifted distribution for this bin
                # Then calculate the shift wrt nominal result, and hence fraction,
                # then save and plot that

                this_syst = self.unfolder.get_exp_syst(syst_dict['label'])
                self.hist_bin_chopper.add_obj(this_syst.syst_shifted_label, this_syst.syst_shifted)
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width(this_syst.syst_shifted_label, ibin, binning_scheme='generator')
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
                                  label="Input stat.",
                                  line_color=ROOT.kRed,
                                  line_style=3,
                                  line_width=3,
                                  marker_size=0,
                                  marker_color=ROOT.kRed)
            c_syst = Contribution(h_syst,
                                  label="Response matrix stat.",
                                  line_color=ROOT.kGray+2,
                                  line_style=3,
                                  line_width=3,
                                  marker_size=0,
                                  marker_color=ROOT.kGray+2)
            c_tot = Contribution(h_total,
                                 label="Total",
                                 line_color=ROOT.kBlack,
                                 line_style=1,
                                 line_width=3,
                                 marker_size=0,
                                 marker_color=ROOT.kBlack)
            entries.extend([c_stat, c_syst, c_tot])

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
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
            self.save_plot(plot, output_filename)

            plot.y_padding_max_log = 50
            plot.set_logy(do_more_labels=False)
            plot.get_modifier().SetMinimum(1E-4)
            log_filename, ext = os.path.splitext(output_filename)
            self.save_plot(plot, log_filename+"_log"+ext)

    def plot_unfolded_with_exp_systs_unnormalised(self):
        """Plot shifted unfolded normalised distributions for each syst"""
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin, binning_scheme='generator')
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin, binning_scheme='generator')
            unfolded_hist_bin_no_errors = unfolded_hist_bin_total_errors.Clone()
            cu.remove_th1_errors(unfolded_hist_bin_no_errors)

            entries = [
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_total_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat. unc.)" % (self.unfolder.tau),
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75),
            ]
            for syst_dict in self.region['experimental_systematics']:
                this_syst = self.unfolder.get_exp_syst(syst_dict['label'])

                self.hist_bin_chopper.add_obj(this_syst.syst_shifted_label, this_syst.syst_shifted)
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width(this_syst.syst_shifted_label, ibin, binning_scheme='generator').Clone()
                cu.remove_th1_errors(syst_unfolded_hist_bin)
                c = Contribution(syst_unfolded_hist_bin,
                                 label=syst_dict['label'],
                                 line_color=syst_dict['colour'], line_width=self.line_width,
                                 line_style=2 if 'down' in syst_dict['label'].lower() else 1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=unfolded_hist_bin_no_errors)
                entries.append(c)


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
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
            self.save_plot(plot, "%s/unfolded_unnormalised_syst_variations_%s_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_syst_fraction_unnormalised(self):
        """Plot varation / central value on absolute hists
        (basically the subplot from plot_unfolded_with_exp_systs_unnormalised)
        """
        def _convert_error_bars_to_error_ratio_hist(h, direction=1):
            """Create hist with bin content = error shift / bin value, 0 error"""
            h_new = h.Clone(h.GetName() + cu.get_unique_str())
            for ix in range(1, h_new.GetNbinsX()+1):
                if h.GetBinContent(ix) > 0:
                    h_new.SetBinContent(ix, 1+(direction*(h.GetBinError(ix) / h.GetBinContent(ix))))
                else:
                    if h.GetBinContent(ix) < 0:
                        h_new.SetBinContent(ix, 1+(direction*(h.GetBinError(ix) / h.GetBinContent(ix))))
                        print("_convert_error_bars_to_error_ratio_hist() warning: bin %d content < 0!" % (ix))
                    else:
                        h_new.SetBinContent(ix, 1)
                h_new.SetBinError(ix, 0)
            return h_new

        def _convert_syst_shift_to_error_ratio_hist(h_syst, h_nominal):
            """Create h_syst / h_nominal without error bars"""
            h_new = h_syst.Clone(h_syst.GetName() + cu.get_unique_str())
            h_new.Divide(h_nominal)
            for ix in range(1, h_new.GetNbinsX()+1):
                h_new.SetBinError(ix, 0)
            return h_new

        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            # Get total for this bin
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', **hbc_args)
            # Get stat. unc. from input for this bin
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width(self.unfolder.stat_uncert_name, **hbc_args)
            # Get stat. unc. from response matrix for this bin
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width(self.unfolder.rsp_uncert_name, **hbc_args)

            unfolded_hist_bin_no_errors = unfolded_hist_bin_total_errors.Clone()
            cu.remove_th1_errors(unfolded_hist_bin_no_errors)

            input_stats = unfolded_hist_bin_stat_errors.Clone()
            input_stats.Divide(unfolded_hist_bin_no_errors)

            total_err = unfolded_hist_bin_total_errors.Clone()
            total_err.Divide(unfolded_hist_bin_no_errors)

            entries = [
                # TOTAL UNCERT
                Contribution(total_err,
                             label="Total uncertainty",
                             line_color=self.plot_styles['unfolded_total_colour'], line_width=0, line_style=2,
                             marker_color=self.plot_styles['unfolded_total_colour'], marker_style=cu.Marker.get('circle'), marker_size=0,
                             fill_style=3354,
                             # fill_style=3005,
                             fill_color=16),
                # INPUT UNCERT
                Contribution(input_stats,
                             label="Data stat.",
                             line_color=self.plot_styles['unfolded_stat_colour'], line_width=0, line_style=3,
                             marker_color=self.plot_styles['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0,
                             fill_style=3545,
                             # fill_style=3003,
                             fill_color=self.plot_styles['unfolded_stat_colour']),
            ]

            # Add experimental systs
            for syst_dict, mark in zip(self.region['experimental_systematics'], cu.Marker().cycle(cycle_filling=True)):
                this_syst = self.unfolder.get_exp_syst(syst_dict['label'])
                syst_unfolded_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width(this_syst.syst_shifted_label, **hbc_args)
                this_syst_hist = _convert_syst_shift_to_error_ratio_hist(syst_unfolded_hist_bin, unfolded_hist_bin_total_errors)
                is_herwig = "shower" in syst_dict['label'].lower() or "herwig" in syst_dict['label'].lower()
                c = Contribution(this_syst_hist,
                                 label="" if "down" in syst_dict['label'] else syst_dict['label'].replace(" up",""),
                                 line_color=syst_dict['colour'],
                                 leg_draw_opt="L" if is_herwig else "P",
                                 line_width=self.line_width if is_herwig else 0,
                                 # line_width=self.line_width,
                                 line_style=5,
                                 marker_color=syst_dict['colour'],
                                 marker_size=0 if is_herwig else 1.25,
                                 marker_style=mark)
                entries.append(c)

            # Add scale syst
            if self.unfolder.scale_uncert_name in self.hist_bin_chopper.objects:
                scale_hist = self.hist_bin_chopper.get_lambda_bin_div_bin_width(self.unfolder.scale_uncert_name, **hbc_args)
                scale_col = self.plot_styles['scale_colour']
                scale_style = dict(line_color=scale_col, line_width=self.line_width, line_style=2,
                                   marker_color=scale_col, marker_style=cu.Marker.get('circle'), marker_size=0,
                                   fill_style=0, fill_color=15)
                entries.extend([
                    Contribution(_convert_error_bars_to_error_ratio_hist(scale_hist),
                                label='#mu scales', leg_draw_opt="L",
                                **scale_style),
                    # add -ve side, no label as we don't want it in legend
                    Contribution(_convert_error_bars_to_error_ratio_hist(scale_hist, -1),
                                **scale_style)
                ])

            # Add pdf syst
            if self.unfolder.pdf_uncert_name in self.hist_bin_chopper.objects:
                pdf_hist = self.hist_bin_chopper.get_lambda_bin_div_bin_width(self.unfolder.pdf_uncert_name, **hbc_args)
                pdf_col = self.plot_styles['pdf_colour']
                pdf_style = dict(line_color=pdf_col, line_width=self.line_width, line_style=6,
                                 marker_color=pdf_col, marker_style=cu.Marker.get('circle'), marker_size=0,
                                 fill_style=0, fill_color=15)
                entries.extend([
                    Contribution(_convert_error_bars_to_error_ratio_hist(pdf_hist),
                                label='PDF', leg_draw_opt="L",
                                **pdf_style),
                    # add -ve side, no label as we don't want it in legend
                    Contribution(_convert_error_bars_to_error_ratio_hist(pdf_hist, -1),
                                **pdf_style)
                ])


            rsp_col = self.plot_styles['rsp_colour']
            rsp_style = dict(line_color=rsp_col, line_width=self.line_width, line_style=3,
                             marker_color=rsp_col, marker_style=cu.Marker.get('circle'), marker_size=0,
                             fill_style=0, fill_color=rsp_col)
            entries.extend([
                # RESPONSE UNCERT
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_rsp_errors),
                             label="Response matrix stat.", leg_draw_opt="L",
                             **rsp_style),
                # Add in the -ve side, but no label as we don't want it in the legend
                Contribution(_convert_error_bars_to_error_ratio_hist(unfolded_hist_bin_rsp_errors, -1),
                             **rsp_style),
            ])


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return
            xlim = qgp.calc_auto_xlim(entries)
            ylim = [0.7, 1.5] if "Dijet" in self.setup.region['name'] else [0.3, 1.9]
            min_total = entries[-1].obj.GetMinimum()
            max_total = entries[-2].obj.GetMaximum()
            if max_total > ylim[1]:
                ylim[1] = max_total*1.1
            if min_total < ylim[0]:
                ylim[0] = min_total*0.9
            plot = Plot(entries,
                        xtitle=self.setup.pt_str,
                        ytitle='Variation / nominal (on absolute distribution)',
                        what="hist",
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        xlim=xlim,
                        ylim=ylim,
                        subplot_type=None)
            self._modify_plot(plot)
            #plot.default_canvas_size = (800, 700)
            plot.legend.SetX1(0.48)
            plot.legend.SetY1(0.58)
            plot.legend.SetX2(0.93)
            plot.legend.SetY2(0.87)
            #if len(entries) > 10: plot.legend.SetNColumns(2)
            plot.legend.SetTextSize(0.03)
            plot.plot("NOSTACK E2 P L") # hard to get one that is points for systs, and line for stats
            plot.set_logx(do_more_labels=False)
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
            self.save_plot(plot, "%s/unfolded_unnormalised_syst_variations_vs_nominal_%s_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_detector_unnormalised(self, alt_detector=None):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            input_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('input_hist_gen_binning_bg_subtracted', ibin, binning_scheme='generator')
            mc_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_mc_reco_gen_binning_bg_subtracted', ibin, binning_scheme='generator')

            entries = [
                Contribution(mc_hist_bin,
                             label="MC (bg-subtracted) [%s]" % (self.setup.region['mc_label']) if alt_detector else "MC (bg-subtracted)",
                             line_color=self.plot_styles['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_mc_colour'], marker_size=0,
                             subplot=input_hist_bin if alt_detector else None),
            ]
            if alt_detector:
                self.hist_bin_chopper.add_obj("alt_hist_mc_reco_gen_binning_bg_subtracted", alt_detector)
                alt_mc_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('alt_hist_mc_reco_gen_binning_bg_subtracted', ibin, binning_scheme='generator')
                entries.append(Contribution(alt_mc_hist_bin,
                                            label="MC (bg-subtracted) [%s]" % self.setup.region['alt_mc_label'],
                                            line_color=self.plot_styles['alt_reco_colour'], line_width=self.line_width,
                                            marker_color=self.plot_styles['alt_reco_colour'], marker_size=0,
                                            subplot=input_hist_bin))
            entries.append(Contribution(input_hist_bin,
                             label="Data (bg-subtracted)",
                             line_color=self.plot_styles['reco_data_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_data_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=None if alt_detector else mc_hist_bin)
            )

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
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
            self.save_plot(plot, "%s/detector_unnormalised_gen_binning_bg_subtracted_%s_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))


class RecoPtBinnedPlotter(BinnedPlotter):
    def __init__(self, setup, bins, hist_bin_chopper, unfolder):
        self.setup = setup
        self.region = setup.region
        self.bins = bins
        self.hist_bin_chopper = hist_bin_chopper

        self.line_width = 2
        self.plot_styles = PLOT_STYLES
        self.pt_bin_plot_args = dict(
            what="hist",
            xtitle=self.setup.detector_title,
            has_data=self.setup.has_data,
            subplot_type='ratio',
            subplot_title="Unfolded / Gen",
            subplot_limits=(0, 2) if self.setup.has_data else (0.75, 1.25),
        )
        self.unfolder = unfolder
        super().__init__()

    def _modify_plot(self, this_plot):
        if self.setup.output_fmt != "pdf":
            _make_thumbnail_canvas(this_plot)
        this_plot.legend.SetX1(0.6)
        this_plot.legend.SetY1(0.68)
        this_plot.legend.SetX2(0.98)
        this_plot.legend.SetY2(0.88)
        this_plot.left_margin = 0.16
        this_plot.lumi = self.setup.lumi

    def _modify_plot_paper(self, this_plot):
        self._modify_plot(this_plot)
        this_plot.is_preliminary = self.setup.is_preliminary
        this_plot.is_supplementary = True

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
        data_total_errors_style = dict(label="Data",
                                       line_color=self.plot_styles['unfolded_total_colour'],
                                       line_width=self.line_width,
                                       line_style=1,
                                       marker_color=self.plot_styles['unfolded_total_colour'],
                                       marker_style=cu.Marker.get('circle'),
                                       marker_size=self.plot_styles['unfolded_marker_size'],
                                       leg_draw_opt="LEP")
        mc_style = dict(label=self.region['mc_label'],
                        line_color=self.plot_styles['gen_colour'],
                        line_width=self.line_width,
                        marker_color=self.plot_styles['gen_colour'],
                        marker_size=self.plot_styles['gen_marker_size'],
                        marker_style=self.plot_styles['gen_marker'],
                        leg_draw_opt="LEP" if self.plot_styles['gen_marker_size'] > 0 else "LE")
        alt_mc_style = dict(label=self.region['alt_mc_label'],
                            line_color=self.plot_styles['alt_gen_colour'],
                            line_width=self.line_width,
                            line_style=1,
                            marker_color=self.plot_styles['alt_gen_colour'],
                            marker_size=self.plot_styles['alt_gen_marker_size'],
                            marker_style=self.plot_styles['alt_gen_marker'],
                            leg_draw_opt="LEP" if self.plot_styles['alt_gen_marker_size'] > 0 else "LE")

        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='detector')
            data_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist', **hbc_args)
            mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco', **hbc_args)

            # For subplot to ensure only MC errors drawn, not MC+data
            data_no_errors = data_hist_bin.Clone()
            cu.remove_th1_errors(data_no_errors)

            entries = [
                Contribution(mc_hist_bin, subplot=data_no_errors, **mc_style),
            ]
            if alt_detector:
                name = "alt_hist_mc_reco"
                self.hist_bin_chopper.add_obj(name, alt_detector)
                alt_mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(name, **hbc_args)
                entries.append(Contribution(alt_mc_hist_bin,
                                            subplot=data_no_errors,
                                            **alt_mc_style)
                )
            # Put data last so it gets drawn on top - we do legend ordering ourselves later anyway
            entries.append(
                Contribution(data_hist_bin, **data_total_errors_style),
            )

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                continue
            max_rel_err = 0.5 if "multiplicity" in self.setup.angle.var.lower() else -1
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_detector_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        xlim=qgp.calc_auto_xlim(entries[0:1], max_rel_err=max_rel_err),  # reduce x axis to where reference prediction is non-0
                        subplot_type='ratio',
                        subplot_title=qgc.SIM_DATA_STR,
                        subplot_limits=(0.5, 2.5) if self.setup.has_data else (0.75, 1.25))
            self._modify_plot_paper(plot)

            # disable adding objects to legend & drawing - we'll do it manually
            plot.do_legend = False
            # plot.legend.SetY1(0.6)
            # plot.legend.SetX1(0.59)
            # plot.legend.SetX2(0.95)
            subplot_draw_opts = "NOSTACK E1"
            plot.plot("NOSTACK E1", subplot_draw_opts)

            dummy_graphs = qgp.do_fancy_legend(chain(entries[-1:], entries[:-1]), plot, use_splitline=False)

            plot.canvas.cd()
            plot.legend.Draw()

            # Create hists for data with error region for ratio
            # Easiest way to get errors right is to do data (with 0 errors)
            # and divide by data (with errors), as if you had MC = data with 0 error
            # data_stat_ratio = data_no_errors.Clone()
            # data_stat_ratio.Divide(unfolded_hist_bin_stat_errors)
            # data_stat_ratio.SetFillStyle(3545)
            # data_stat_ratio.SetFillColor(self.plot_styles['unfolded_stat_colour'])
            # data_stat_ratio.SetLineWidth(0)
            # data_stat_ratio.SetMarkerSize(0)

            data_total_ratio = data_no_errors.Clone()
            data_total_ratio.Divide(data_hist_bin)
            data_total_ratio.SetFillStyle(3354)
            data_total_ratio.SetFillColor(16)
            data_total_ratio.SetLineWidth(0)
            data_total_ratio.SetMarkerSize(0)

            # now draw the data error shaded area
            # this is a bit hacky - basically draw them on the ratio pad,
            # then redraw the existing hists & line to get them ontop
            # note that we use "same" for all - this is to keep the original axes
            # (we may want to rethink this later?)
            plot.subplot_pad.cd()
            draw_opt = "E2 SAME"
            data_total_ratio.Draw(draw_opt)
            plot.subplot_container.Draw("SAME" + subplot_draw_opts)
            plot.subplot_line.Draw()

            # Add subplot legend
            x_left = 0.25
            y_bottom = 0.77
            width = 0.55
            height = 0.13
            plot.subplot_legend = ROOT.TLegend(x_left, y_bottom, x_left+width, y_bottom+height)
            plot.subplot_legend.AddEntry(data_total_ratio, qgc.DATA_STAT_UNC_STR, "F")
            plot.subplot_legend.SetFillStyle(0)
            plot.subplot_legend.SetTextSize(0.085)
            plot.subplot_legend.Draw()

            plot.canvas.cd()

            stp = self.setup
            fname = f"detector_reco_binning_{self.setup.append}_bin_{ibin:d}_divBinWidth{self.setup.paper_str}.{self.setup.output_fmt}"
            self.save_plot(plot, os.path.join(stp.output_dir, fname))

    def plot_detector_normalised_bg_subtracted(self, alt_detector=None):
        data_total_errors_style = dict(label="Data (bg-subtracted)",
                                       line_color=self.plot_styles['reco_data_colour'], line_width=self.line_width, line_style=1,
                                       marker_color=self.plot_styles['reco_data_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75)

        mc_style = dict(label="%s (bg-subtracted)" % self.region['mc_label'],
                        line_color=self.plot_styles['reco_mc_colour'], line_width=self.line_width,
                        marker_color=self.plot_styles['reco_mc_colour'], marker_size=0)

        alt_mc_style = dict(label="%s (bg-subtracted)" % self.region['alt_mc_label'],
                            line_color=self.plot_styles['alt_reco_colour'], line_width=self.line_width, #line_style=2,
                            marker_color=self.plot_styles['alt_reco_colour'], marker_size=0)

        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='detector')
            data_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_bg_subtracted', **hbc_args)
            mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', **hbc_args)

            # For subplot to ensure only MC errors drawn, not MC+data
            data_no_errors = data_hist_bin.Clone()
            cu.remove_th1_errors(data_no_errors)

            entries = [
                Contribution(mc_hist_bin, subplot=data_no_errors, **mc_style),
            ]
            if alt_detector:
                name = "alt_hist_mc_reco_bg_subtracted"
                self.hist_bin_chopper.add_obj(name, alt_detector)
                alt_mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(name, **hbc_args)
                entries.append(Contribution(alt_mc_hist_bin,
                                            subplot=data_no_errors,
                                            **alt_mc_style)
                )
            # Put data last so it gets drawn on top - we do legend ordering ourselves later anyway
            entries.append(
                Contribution(data_hist_bin, **data_total_errors_style),
            )

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_detector_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='Simulation / data',
                        subplot_limits=(0, 2))
            self._modify_plot(plot)

            # disable adding objects to legend & drawing - we'll do it manually
            plot.do_legend = False
            subplot_draw_opts = "NOSTACK E1"
            plot.plot("NOSTACK E1", subplot_draw_opts)

            # Create dummy graphs with the same styling to put into the legend
            dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
            dummy_data = Contribution(dummy_gr.Clone(), leg_draw_opt="LEP", **data_total_errors_style)
            dummy_mc = Contribution(dummy_gr.Clone(), leg_draw_opt="LE", **mc_style)
            dummy_alt_mc = Contribution(dummy_gr.Clone(), leg_draw_opt="LE", **alt_mc_style)
            # Add them to the legend and draw it
            for cont in [dummy_data, dummy_mc, dummy_alt_mc]:
                plot.legend.AddEntry(cont.obj, cont.label, cont.leg_draw_opt)
            plot.canvas.cd()
            plot.legend.Draw()

            # Create hists for data with error region for ratio
            # Easiest way to get errors right is to do data (with 0 errors)
            # and divide by data (with errors), as if you had MC = data with 0 error
            # data_stat_ratio = data_no_errors.Clone()
            # data_stat_ratio.Divide(unfolded_hist_bin_stat_errors)
            # data_stat_ratio.SetFillStyle(3545)
            # data_stat_ratio.SetFillColor(self.plot_styles['unfolded_stat_colour'])
            # data_stat_ratio.SetLineWidth(0)
            # data_stat_ratio.SetMarkerSize(0)

            data_total_ratio = data_no_errors.Clone()
            data_total_ratio.Divide(data_hist_bin)
            data_total_ratio.SetFillStyle(3354)
            data_total_ratio.SetFillColor(16)
            data_total_ratio.SetLineWidth(0)
            data_total_ratio.SetMarkerSize(0)

            # now draw the data error shaded area
            # this is a bit hacky - basically draw them on the ratio pad,
            # then redraw the existing hists & line to get them ontop
            # note that we use "same" for all - this is to keep the original axes
            # (we may want to rethink this later?)
            plot.subplot_pad.cd()
            draw_opt = "E2 SAME"
            data_total_ratio.Draw(draw_opt)
            plot.subplot_container.Draw("SAME" + subplot_draw_opts)
            plot.subplot_line.Draw()
            plot.canvas.cd()

            self.save_plot(plot, "%s/detector_reco_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_folded_unfolded_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_bg_subtracted', ibin, binning_scheme='detector')
            folded_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_unfolded', ibin, binning_scheme='detector')

            entries = [
                Contribution(input_hist_bin,
                             label="Unfolding input (bg-subtracted)",
                             line_color=self.plot_styles['reco_unfolding_input_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_unfolding_input_colour'], marker_style=cu.Marker.get('circle'), marker_size=0),
                Contribution(folded_unfolded_hist_bin,
                             label="Folded unfolded",
                             line_color=self.plot_styles['reco_folded_unfolded_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_folded_unfolded_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=input_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_detector_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='Folded / input',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/detector_folded_unfolded_only_data_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_folded_unfolded_with_mc_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_bg_subtracted', ibin, binning_scheme='detector')
            folded_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_unfolded', ibin, binning_scheme='detector')

            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_styles['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_mc_colour'], marker_style=cu.Marker.get('circle'), marker_size=0),
                Contribution(input_hist_bin,
                             label="Unfolding input (bg-subtracted)",
                             line_color=self.plot_styles['reco_unfolding_input_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_unfolding_input_colour'], marker_style=cu.Marker.get('circle'), marker_size=0,
                             subplot=mc_reco_hist_bin),
                Contribution(folded_unfolded_hist_bin,
                             label="Folded unfolded",
                             line_color=self.plot_styles['reco_folded_unfolded_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_folded_unfolded_colour'],marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_reco_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_detector_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='* / MC',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/detector_folded_unfolded_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_folded_gen_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')
            folded_mc_truth_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_mc_truth', ibin, binning_scheme='detector')

            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_styles['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_mc_colour'], marker_style=cu.Marker.get('circle'), marker_size=0),
                Contribution(folded_mc_truth_hist_bin,
                             label="Folded gen",
                             line_color=self.plot_styles['reco_folded_mc_truth_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_folded_mc_truth_colour'], marker_style=cu.Marker.get('circle'), marker_size=0.75,
                             subplot=mc_reco_hist_bin),
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_detector_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='Folded / MC reco',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/detector_folded_gen_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_detector_with_model_systs_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # add nominal
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')
            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_styles['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_mc_colour'], marker_style=cu.Marker.get('circle'), marker_size=0),
            ]
            # add model variations
            for model_dict in self.region['model_systematics']:
                model_unfolder = model_dict['unfolder']
                model_label = model_dict['label']
                model_label_no_spaces = cu.no_space_str(model_dict['label'])
                name = 'hist_mc_reco_bg_subtracted_model_%s' % model_label_no_spaces
                self.hist_bin_chopper.add_obj(name, model_unfolder.input_hist_bg_subtracted)
                entries.append(
                    Contribution(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(name, ibin, binning_scheme='detector'),
                                 label="%s (bg-subtracted)" % (model_label),
                                 line_color=model_dict['colour'], line_width=1, line_style=1,
                                 marker_color=model_dict['colour'], marker_size=0,
                                 subplot=mc_reco_hist_bin)
                )


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_detector_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='#splitline{Variation /}{Nominal}',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            if len(entries) > 5:
                plot.legend.SetNColumns(2)
            if len(entries) > 15:
                plot.legend.SetNColumns(3)
            ROOT.gStyle.SetPalette(55)
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/detector_reco_binning_bg_subtracted_model_systs_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
            ROOT.gStyle.SetPalette(ROOT.kViridis)

    def plot_detector_with_model_systs_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # add nominal
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')
            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_styles['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_mc_colour'], marker_style=cu.Marker.get('circle'), marker_size=0),
            ]
            # add model variations
            for model_dict in self.region['model_systematics']:
                model_unfolder = model_dict['unfolder']
                model_label = model_dict['label']
                model_label_no_spaces = cu.no_space_str(model_dict['label'])
                name = 'hist_mc_reco_bg_subtracted_model_%s' % model_label_no_spaces
                self.hist_bin_chopper.add_obj(name, model_unfolder.input_hist_bg_subtracted)
                entries.append(
                    Contribution(self.hist_bin_chopper.get_pt_bin_div_bin_width(name, ibin, binning_scheme='detector'),
                                 label="%s (bg-subtracted)" % (model_label),
                                 line_color=model_dict['colour'], line_width=1, line_style=1,
                                 marker_color=model_dict['colour'], marker_size=0,
                                 subplot=mc_reco_hist_bin)
                )


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='#splitline{Variation /}{Nominal}',
                        subplot_limits=(0.5, 1.5),)
            self._modify_plot(plot)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            if len(entries) > 5:
                plot.legend.SetNColumns(2)
            if len(entries) > 15:
                plot.legend.SetNColumns(3)
            ROOT.gStyle.SetPalette(55)
            plot.plot("NOSTACK E1")
            self.save_plot(plot, "%s/detector_reco_binning_bg_subtracted_model_systs_unnormalised_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
            ROOT.gStyle.SetPalette(ROOT.kViridis)

    def plot_detector_with_pdf_systs_normalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # add nominal
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')
            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_styles['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_mc_colour'], marker_style=cu.Marker.get('circle'), marker_size=0),
            ]
            # add PDF variations
            for pdf_dict in self.region['pdf_systematics']:
                pdf_unfolder = pdf_dict['unfolder']
                pdf_label = pdf_dict['label']
                pdf_label_no_spaces = cu.no_space_str(pdf_dict['label'])
                name = 'hist_mc_reco_bg_subtracted_PDF_%s' % pdf_label_no_spaces
                self.hist_bin_chopper.add_obj(name, pdf_unfolder.input_hist_bg_subtracted)
                entries.append(
                    Contribution(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(name, ibin, binning_scheme='detector'),
                                 label="%s" % (pdf_label),
                                 line_color=pdf_dict['colour'], line_width=1, line_style=1,
                                 marker_color=pdf_dict['colour'], marker_size=0,
                                 subplot=mc_reco_hist_bin)
                )


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_detector_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='#splitline{Variation /}{Nominal}',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            if len(entries) > 5:
                plot.legend.SetNColumns(2)
            if len(entries) > 15:
                plot.legend.SetNColumns(3)
            ROOT.gStyle.SetPalette(55)
            plot.plot("NOSTACK PLC PMC E1")
            self.save_plot(plot, "%s/detector_reco_binning_bg_subtracted_pdf_systs_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
            ROOT.gStyle.SetPalette(ROOT.kViridis)

    def plot_detector_with_pdf_systs_unnormalised(self):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # add nominal
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('hist_mc_reco_bg_subtracted', ibin, binning_scheme='detector')
            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_styles['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_styles['reco_mc_colour'], marker_style=cu.Marker.get('circle'), marker_size=0),
            ]
            # add PDF variations
            for pdf_dict in self.region['pdf_systematics']:
                pdf_unfolder = pdf_dict['unfolder']
                pdf_label = pdf_dict['label']
                pdf_label_no_spaces = cu.no_space_str(pdf_dict['label'])
                name = 'hist_mc_reco_bg_subtracted_PDF_%s' % pdf_label_no_spaces
                self.hist_bin_chopper.add_obj(name, pdf_unfolder.input_hist_bg_subtracted)
                entries.append(
                    Contribution(self.hist_bin_chopper.get_pt_bin_div_bin_width(name, ibin, binning_scheme='detector'),
                                 label="%s" % (pdf_label),
                                 line_color=pdf_dict['colour'], line_width=1, line_style=1,
                                 marker_color=pdf_dict['colour'], marker_size=0,
                                 subplot=mc_reco_hist_bin)
                )


            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_unnormalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='#splitline{Variation /}{Nominal}',
                        subplot_limits=(0.5, 1.5),)
            self._modify_plot(plot)
            plot.legend.SetX1(0.53)
            plot.legend.SetY1(0.7)
            plot.legend.SetX2(0.96)
            plot.legend.SetY2(0.88)
            if len(entries) > 5:
                plot.legend.SetNColumns(2)
            if len(entries) > 15:
                plot.legend.SetNColumns(3)
            ROOT.gStyle.SetPalette(55)
            plot.plot("NOSTACK PLC PMC E1")
            self.save_plot(plot, "%s/detector_reco_binning_bg_subtracted_pdf_systs_unnormalised_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))
            ROOT.gStyle.SetPalette(ROOT.kViridis)


@profile
def do_binned_plots_per_region_angle(setup, do_binned_gen_pt, do_binned_gen_lambda, do_binned_reco_pt, only_paper_plots=False):
    """Do individual binned plots, can select which binning(s) to plot over"""
    region = setup.region
    # Note that experimental systs are only different response matrices, and are stored in the main unfolder
    has_exp_systs = len(region['experimental_systematics']) > 0
    has_scale_systs = len(region['scale_systematics']) > 0
    has_model_systs = len(region['model_systematics']) > 0
    has_pdf_systs = len(region['pdf_systematics']) > 0
    has_jackknife_input_vars = len(region['jackknife_input_variations']) > 0
    has_jackknife_response_vars = len(region['jackknife_response_variations']) > 0

    if has_exp_systs: print("We have experimental systs")
    if has_model_systs: print("We have model systs")
    if has_scale_systs: print("We have scale systs")
    if has_pdf_systs: print("We have pdf systs")
    if has_jackknife_input_vars: print("We have jackknife input variations")
    if has_jackknife_response_vars: print("We have jackknife response variations")

    unfolder = region['unfolder']
    unreg_unfolder = region.get('unreg_unfolder', None)
    alt_unfolder = region.get('alt_unfolder', None)
    alt_hist_truth = region.get('alt_hist_mc_gen', None)
    alt_hist_reco = region.get('alt_hist_mc_reco', None)
    alt_hist_reco_bg_subtracted = region.get('alt_hist_mc_reco_bg_subtracted', None)
    alt_hist_reco_bg_subtracted_gen_binning = region.get('alt_hist_mc_reco_bg_subtracted_gen_binning', None)

    # hbc = HistBinChopper(generator_binning=unfolder.generator_binning.FindNode("generatordistribution"),
    #                      detector_binning=unfolder.detector_binning.FindNode("detectordistribution"))
    # hbc.add_obj("hist_truth", unfolder.hist_truth)
    # hbc.add_obj('unfolded', unfolder.get_output())
    # hbc.add_obj('unfolded_stat_err', unfolder.get_unfolded_with_ematrix_stat())
    # hbc.update(unfolder.hist_bin_chopper)   # update the HistBinChopper with the new normalised systematics already produced in unfolder
    hbc = unfolder.hist_bin_chopper

    if do_binned_gen_pt:
        # Iterate through pt bins - gen binning
        # ------------------------------------------------------------------
        print("Doing GenPtBinnedPlotter...")
        pt_bins = unfolder.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=True)
        gen_pt_binned_plotter = GenPtBinnedPlotter(setup=setup,
                                                   bins=pt_bins,
                                                   hist_bin_chopper=hbc,
                                                   unfolder=unfolder)
        if alt_hist_truth:
            print("...doing alt truth")
            gen_pt_binned_plotter.hist_bin_chopper.add_obj('alt_hist_truth', alt_hist_truth)
            gen_pt_binned_plotter.plot_unfolded_with_alt_truth_normalised(do_chi2=True, do_zoomed=True)

        if not only_paper_plots:
            gen_pt_binned_plotter.plot_unfolded_normalised()
            gen_pt_binned_plotter.plot_unfolded_unnormalised()
            gen_pt_binned_plotter.plot_total_ematrix()

            if unfolder.tau > 0 and unreg_unfolder:
                print("...doing unregularised vs regularised")
                gen_pt_binned_plotter.hist_bin_chopper.add_obj("unreg_unfolded", unreg_unfolder.unfolded)
                gen_pt_binned_plotter.hist_bin_chopper.add_obj("truth_template", unreg_unfolder.truth_template)
                gen_pt_binned_plotter.plot_unfolded_with_unreg_normalised()
                gen_pt_binned_plotter.plot_unfolded_with_unreg_unnormalised()
                gen_pt_binned_plotter.plot_unfolded_with_template_normalised()
                gen_pt_binned_plotter.plot_unfolded_with_template_unnormalised()

            if alt_unfolder:
                print("...doing alt unfolder")
                gen_pt_binned_plotter.plot_unfolded_with_alt_response_normalised(alt_unfolder=alt_unfolder)
                gen_pt_binned_plotter.plot_unfolded_with_alt_response_truth_normalised(alt_unfolder=alt_unfolder)

            if has_exp_systs:
                print("...doing exp systs")
                # gen_pt_binned_plotter.plot_uncertainty_shifts_normalised()
                gen_pt_binned_plotter.plot_unfolded_with_exp_systs_normalised()
                gen_pt_binned_plotter.plot_unfolded_with_exp_systs_unnormalised()

            if has_scale_systs:
                print("...doing scale systs")
                gen_pt_binned_plotter.plot_unfolded_with_scale_systs_normalised()
                gen_pt_binned_plotter.plot_unfolded_with_scale_systs_unnormalised()
                # Do a set of individual plots for these scale variations
                # for syst_dict in region['scale_systematics']:
                #     print(".......", syst_dict['label'])
                #     this_setup = copy(setup)
                #     this_setup.output_dir = os.path.join(setup.output_dir, "scaleSyst_"+cu.no_space_str(syst_dict['label']))
                #     syst_gen_pt_binned_plotter = GenPtBinnedPlotter(setup=this_setup,
                #                                                     bins=unfolder.pt_bin_edges_gen,
                #                                                     hist_bin_chopper=syst_dict['unfolder'].hist_bin_chopper,
                #                                                     unfolder=syst_dict['unfolder'])
                #     syst_gen_pt_binned_plotter.plot_unfolded_unnormalised()
                #     syst_gen_pt_binned_plotter.plot_unfolded_normalised()

            if has_model_systs:
                print("...doing model systs")
                gen_pt_binned_plotter.plot_unfolded_with_model_systs_normalised()
                gen_pt_binned_plotter.plot_unfolded_with_model_systs_unnormalised()
                # Do a set of individual plots for these model variations
                for syst_dict in region['model_systematics']:
                    print(".......", syst_dict['label'])
                    this_setup = copy(setup)
                    this_setup.output_dir = os.path.join(setup.output_dir, "modelSyst_"+cu.no_space_str(syst_dict['label']))
                    this_setup.region['mc_label'] = syst_dict['label']  # since that's the thing being unfolded
                    syst_gen_pt_binned_plotter = GenPtBinnedPlotter(setup=this_setup,
                                                                    bins=pt_bins,
                                                                    hist_bin_chopper=syst_dict['unfolder'].hist_bin_chopper,
                                                                    unfolder=syst_dict['unfolder'])
                    syst_gen_pt_binned_plotter.plot_unfolded_unnormalised()
                    syst_gen_pt_binned_plotter.plot_unfolded_normalised()
                    if syst_dict['unfolder'].tau > 0:
                        # plot with template since template here is different to nominal
                        syst_dict['unfolder'].hist_bin_chopper.add_obj("truth_template", syst_dict['unfolder'].truth_template)
                        syst_dict['unfolder'].hist_bin_chopper.add_obj("alt_hist_truth", alt_hist_truth)
                        syst_gen_pt_binned_plotter.plot_unfolded_with_template_normalised()
                        syst_gen_pt_binned_plotter.plot_unfolded_with_template_unnormalised()

            # if has_pdf_systs:
            #     print("...doing pdf systs")
            #     gen_pt_binned_plotter.plot_unfolded_with_pdf_systs_normalised()
            #     gen_pt_binned_plotter.plot_unfolded_with_pdf_systs_unnormalised()

                # # Do a set of individual plots for these PDF variations
                # for syst_dict in region['pdf_systematics']:
                #     print(".......", syst_dict['label'])
                #     this_setup = copy(setup)
                #     this_setup.output_dir = os.path.join(setup.output_dir, "pdfSyst", cu.no_space_str(syst_dict['label']))
                #     syst_gen_pt_binned_plotter = GenPtBinnedPlotter(setup=this_setup,
                #                                                     bins=unfolder.pt_bin_edges_gen,
                #                                                     hist_bin_chopper=syst_dict['unfolder'].hist_bin_chopper,
                #                                                     unfolder=syst_dict['unfolder'])
                #     syst_gen_pt_binned_plotter.plot_unfolded_unnormalised()
                #     syst_gen_pt_binned_plotter.plot_unfolded_normalised()

            if has_jackknife_input_vars:
                print("...doing jackknife input variations")
                gen_pt_binned_plotter.plot_unfolded_with_jackknife_input_normalised()
                gen_pt_binned_plotter.plot_jackknife_residuals(region['jackknife_input_variations'])
                # Do a set of individual plots for these jackknife variations
                for jk_dict in region['jackknife_input_variations']:
                    print(".......", jk_dict['label'])
                    this_setup = copy(setup)
                    this_setup.output_dir = os.path.join(setup.output_dir, "jackknife_input", jk_dict['label'])
                    jk_gen_pt_binned_plotter = GenPtBinnedPlotter(setup=this_setup,
                                                                  bins=pt_bins,
                                                                  hist_bin_chopper=jk_dict['unfolder'].hist_bin_chopper,
                                                                  unfolder=jk_dict['unfolder'])
                    jk_gen_pt_binned_plotter.plot_unfolded_unnormalised()
                    jk_gen_pt_binned_plotter.plot_unfolded_normalised()
                    if jk_dict['unfolder'].tau > 0:
                        jk_dict['unfolder'].hist_bin_chopper.add_obj("truth_template", jk_dict['unfolder'].truth_template)
                        jk_dict['unfolder'].hist_bin_chopper.add_obj("alt_hist_truth", alt_hist_truth)
                        jk_gen_pt_binned_plotter.plot_unfolded_with_template_normalised()
                        jk_gen_pt_binned_plotter.plot_unfolded_with_template_unnormalised()

            if has_jackknife_response_vars:
                print("...doing jackknife response variations")
                gen_pt_binned_plotter.plot_unfolded_with_jackknife_response_normalised()
                gen_pt_binned_plotter.plot_jackknife_residuals(region['jackknife_response_variations'])
                # Do a set of individual plots for these jackknife variations
                for jk_dict in region['jackknife_response_variations']:
                    print(".......", jk_dict['label'])
                    this_setup = copy(setup)
                    this_setup.output_dir = os.path.join(setup.output_dir, "jackknife_response", jk_dict['label'])
                    jk_gen_pt_binned_plotter = GenPtBinnedPlotter(setup=this_setup,
                                                                  bins=pt_bins,
                                                                  hist_bin_chopper=jk_dict['unfolder'].hist_bin_chopper,
                                                                  unfolder=jk_dict['unfolder'])
                    jk_gen_pt_binned_plotter.plot_unfolded_unnormalised()
                    jk_gen_pt_binned_plotter.plot_unfolded_normalised()

        print("...doing uncert fraction")
        gen_pt_binned_plotter.plot_syst_fraction_normalised()

        # FIXME delete this next time, should be done in unfolding.py
        # if has_scale_systs:
        #     unfolder.create_scale_syst_uncertainty_per_pt_bin(region['scale_systematics'])
        # if has_pdf_systs:
        #     unfolder.create_pdf_syst_uncertainty_per_pt_bin(region['pdf_systematics'])
        # unfolder.setup_absolute_results_per_pt_bin()

        gen_pt_binned_plotter.plot_syst_fraction_unnormalised()

        # # if has_data:
        # print("...doing detector-level")
        # gen_pt_binned_plotter.hist_bin_chopper.add_obj("input_hist_gen_binning_bg_subtracted", unfolder.input_hist_gen_binning_bg_subtracted)
        # gen_pt_binned_plotter.hist_bin_chopper.add_obj("hist_mc_reco_gen_binning_bg_subtracted", unfolder.hist_mc_reco_gen_binning_bg_subtracted)
        # gen_pt_binned_plotter.plot_detector_normalised_bg_subtracted(alt_detector=alt_hist_reco_bg_subtracted_gen_binning)

    if do_binned_gen_lambda:
        # Iterate through lambda bins - gen binning
        # ------------------------------------------------------------------
        print("Doing GenLambdaBinnedPlotter...")
        if isinstance(unfolder.binning_handler.get_binning_scheme('generator'), PtVarPerPtBinning):
            warnings.warn(cu.pcolors.WARNING +
                          "Using PtVarPerPtBinning, no guarantee of same variable bining, "
                          "skipping gen lambda binned plots" + cu.pcolors.ENDC)
        else:
            gen_lambda_binned_plotter = GenLambdaBinnedPlotter(setup=setup,
                                                               bins=unfolder.variable_bin_edges_gen,
                                                               hist_bin_chopper=hbc,
                                                               # this is the same object as gen_pt_binned_plotter, so has all the objects already
                                                               unfolder=unfolder)

            if not only_paper_plots:
                gen_lambda_binned_plotter.plot_unfolded_unnormalised()

                if unfolder.tau > 0 and unreg_unfolder:
                    print("...doing unregularised vs regularised")
                    gen_lambda_binned_plotter.hist_bin_chopper.add_obj("unreg_unfolded", unreg_unfolder.unfolded)
                    gen_lambda_binned_plotter.hist_bin_chopper.add_obj("truth_template", unreg_unfolder.truth_template)
                    gen_lambda_binned_plotter.plot_unfolded_with_unreg_unnormalised()
                    gen_lambda_binned_plotter.plot_unfolded_with_template_unnormalised()

                if alt_unfolder:
                    print("...doing alt unfolder")
                    gen_lambda_binned_plotter.plot_unfolded_with_alt_response_unnormalised(alt_unfolder=alt_unfolder)

                if has_exp_systs:
                    # gen_lambda_binned_plotter.plot_uncertainty_shifts_unnormalised()
                    print("...doing exp systs")
                    gen_lambda_binned_plotter.plot_unfolded_with_exp_systs_unnormalised()

                if has_scale_systs:
                    print("...doing scale systs")
                    gen_lambda_binned_plotter.plot_unfolded_with_scale_systs_unnormalised()

                if has_model_systs:
                    print("...doing model systs")
                    gen_lambda_binned_plotter.plot_unfolded_with_model_systs_unnormalised()
                    # Do a set of individual plots for these model variations
                    for syst_dict in region['model_systematics']:
                        print(".......", syst_dict['label'])
                        this_setup = copy(setup)
                        this_setup.output_dir = os.path.join(setup.output_dir,
                                                             "modelSyst_" + cu.no_space_str(syst_dict['label']))
                        syst_gen_lambda_binned_plotter = GenLambdaBinnedPlotter(setup=this_setup,
                                                                                bins=unfolder.variable_bin_edges_gen,
                                                                                hist_bin_chopper=syst_dict[
                                                                                    'unfolder'].hist_bin_chopper,
                                                                                unfolder=syst_dict['unfolder'])
                        syst_gen_lambda_binned_plotter.plot_unfolded_unnormalised()
                        if syst_dict['unfolder'].tau > 0:
                            syst_gen_lambda_binned_plotter.hist_bin_chopper.add_obj("unreg_unfolded",
                                                                                    syst_dict['unfolder'].unfolded)
                            syst_gen_lambda_binned_plotter.hist_bin_chopper.add_obj("truth_template",
                                                                                    syst_dict['unfolder'].truth_template)
                            syst_gen_lambda_binned_plotter.plot_unfolded_with_unreg_unnormalised()
                            syst_gen_lambda_binned_plotter.plot_unfolded_with_template_unnormalised()

                # if has_pdf_systs:
                #     print("...doing pdf systs")
                #     gen_lambda_binned_plotter.plot_unfolded_with_pdf_systs_unnormalised()
                #     # Do a set of individual plots for these PDF variations
                #     for syst_dict in region['pdf_systematics']:
                #         print(".......", syst_dict['label'])
                #         this_setup = copy(setup)
                #         this_setup.output_dir = os.path.join(setup.output_dir, "pdfSyst", cu.no_space_str(syst_dict['label']))
                #         syst_gen_lambda_binned_plotter = GenLambdaBinnedPlotter(setup=this_setup,
                #                                                                 bins=unfolder.pt_bin_edges_gen,
                #                                                                 hist_bin_chopper=syst_dict['unfolder'].hist_bin_chopper,
                #                                                                 unfolder=syst_dict['unfolder'])
                #         syst_gen_lambda_binned_plotter.plot_unfolded_unnormalised()
                #         syst_gen_lambda_binned_plotter.plot_unfolded_normalised()

                if has_jackknife_input_vars:
                    # Do a set of individual plots for these jackknife variations
                    for jk_dict in region['jackknife_input_variations']:
                        print(".......", jk_dict['label'])
                        this_setup = copy(setup)
                        this_setup.output_dir = os.path.join(setup.output_dir, "jackknife_input", jk_dict['label'])
                        jk_gen_lambda_binned_plotter = GenLambdaBinnedPlotter(setup=this_setup,
                                                                              bins=unfolder.pt_bin_edges_gen,
                                                                              hist_bin_chopper=jk_dict[
                                                                                  'unfolder'].hist_bin_chopper,
                                                                              unfolder=jk_dict['unfolder'])
                        jk_gen_lambda_binned_plotter.plot_unfolded_unnormalised()
                        if jk_dict['unfolder'].tau > 0:
                            # plot with template since template different to nominal template
                            jk_gen_lambda_binned_plotter.hist_bin_chopper.add_obj("unreg_unfolded",
                                                                                  jk_dict['unfolder'].unfolded)
                            jk_gen_lambda_binned_plotter.hist_bin_chopper.add_obj("truth_template",
                                                                                  jk_dict['unfolder'].truth_template)
                            jk_gen_lambda_binned_plotter.plot_unfolded_with_unreg_unnormalised()
                            jk_gen_lambda_binned_plotter.plot_unfolded_with_template_unnormalised()

                if has_jackknife_response_vars:
                    # Do a set of individual plots for these jackknife variations
                    for jk_dict in region['jackknife_response_variations']:
                        print(".......", jk_dict['label'])
                        this_setup = copy(setup)
                        this_setup.output_dir = os.path.join(setup.output_dir, "jackknife_response", jk_dict['label'])
                        jk_gen_lambda_binned_plotter = GenLambdaBinnedPlotter(setup=this_setup,
                                                                              bins=unfolder.pt_bin_edges_gen,
                                                                              hist_bin_chopper=jk_dict[
                                                                                  'unfolder'].hist_bin_chopper,
                                                                              unfolder=jk_dict['unfolder'])
                        jk_gen_lambda_binned_plotter.plot_unfolded_unnormalised()



                # if has_data:
                print("...doing detector-level")
                gen_lambda_binned_plotter.hist_bin_chopper.add_obj("input_hist_gen_binning_bg_subtracted",
                                                                   unfolder.input_hist_gen_binning_bg_subtracted)
                gen_lambda_binned_plotter.hist_bin_chopper.add_obj("hist_mc_reco_gen_binning_bg_subtracted",
                                                                   unfolder.hist_mc_reco_gen_binning_bg_subtracted)
                gen_lambda_binned_plotter.plot_detector_unnormalised(alt_detector=alt_hist_reco_bg_subtracted_gen_binning)

            # bit of a hack, should do it in main unfolding.py script
            # if has_scale_systs:
            #     unfolder.create_scale_syst_uncertainty_per_lambda_bin(region['scale_systematics'])
            # if has_pdf_systs:
            #     unfolder.create_pdf_syst_uncertainty_per_lambda_bin(region['pdf_systematics'])
            # unfolder.setup_absolute_results_per_lambda_bin()
            print("...doing uncert fraction")
            gen_lambda_binned_plotter.plot_syst_fraction_unnormalised()

    if do_binned_reco_pt:
        # Iterate through pt bins - reco binning
        # ------------------------------------------------------------------
        print("Doing RecoPtBinnedPlotter...")
        hbc.add_obj("hist_mc_reco", unfolder.hist_mc_reco)
        hbc.add_obj("hist_mc_reco_bg_subtracted", unfolder.hist_mc_reco_bg_subtracted)
        hbc.add_obj("input_hist_bg_subtracted", unfolder.input_hist_bg_subtracted)
        hbc.add_obj("input_hist", unfolder.input_hist)
        hbc.add_obj("folded_unfolded", unfolder.folded_unfolded)
        hbc.add_obj("folded_mc_truth", unfolder.folded_mc_truth)
        pt_bins = unfolder.binning_handler.get_pt_bins(binning_scheme='detector', is_signal_region=True)
        reco_pt_binned_plotter = RecoPtBinnedPlotter(setup=setup,
                                                     bins=pt_bins,
                                                     hist_bin_chopper=hbc,
                                                     unfolder=unfolder)
        reco_pt_binned_plotter.plot_detector_normalised(alt_detector=alt_hist_reco)
        # reco_pt_binned_plotter.plot_detector_normalised_bg_subtracted(alt_detector=alt_hist_reco_bg_subtracted)
        if not only_paper_plots:
            reco_pt_binned_plotter.plot_folded_unfolded_normalised()
            reco_pt_binned_plotter.plot_folded_unfolded_with_mc_normalised()
            reco_pt_binned_plotter.plot_folded_gen_normalised()

        # if has_model_systs:
        #     reco_pt_binned_plotter.plot_detector_with_model_systs_normalised()
        #     reco_pt_binned_plotter.plot_detector_with_model_systs_unnormalised()

        # if has_pdf_systs:
        #     reco_pt_binned_plotter.plot_detector_with_pdf_systs_normalised()
        #     reco_pt_binned_plotter.plot_detector_with_pdf_systs_unnormalised()

    return hbc


class BigNormalised1DPlotter(object):

    # TODO: Some of this seems very overlapped with MyUnfolderPlotter...

    def __init__(self, setup, hist_bin_chopper, plot_with_bin_widths=False):
        self.setup = setup
        self.hist_bin_chopper = hist_bin_chopper
        self.unfolder = setup.region['unfolder']
        self.binning_handler = self.unfolder.binning_handler
        self.plot_with_bin_widths = plot_with_bin_widths  # use bin widths in plot instead of uniform bin widths

        # cache all these things since they don't (shouldn't!) change
        # but the user should more often use the methods below, since it takes
        # care if you need to change binning scheme
        self.pt_bin_edges_generator = self.binning_handler.get_pt_bins('generator', is_signal_region=True)
        self.num_pt_bins_generator = len(self.pt_bin_edges_generator) - 1
        self.all_bin_edges_generator = self._calc_bin_edges('generator')
        self.nbins_generator = len(self.all_bin_edges_generator) - 1

        self.pt_bin_edges_detector = self.binning_handler.get_pt_bins('detector', is_signal_region=True)
        self.num_pt_bins_detector = len(self.pt_bin_edges_detector) - 1
        self.all_bin_edges_detector = self._calc_bin_edges('detector')
        self.nbins_detector = len(self.all_bin_edges_detector) - 1

        self.line_width = 1

        self._cache_1d = {}

    # urgh, is there a better way?
    def pt_bin_edges(self, binning_scheme):
        return self.pt_bin_edges_generator if binning_scheme == 'generator' else self.pt_bin_edges_detector

    def num_pt_bins(self, binning_scheme):
        return self.num_pt_bins_generator if binning_scheme == 'generator' else self.num_pt_bins_detector

    def nbins(self, binning_scheme):
        return self.nbins_generator if binning_scheme == 'generator' else self.nbins_detector

    def all_bin_edges(self, binning_scheme):
        return self.all_bin_edges_generator if binning_scheme == 'generator' else self.all_bin_edges_detector

    def _calc_bin_edges(self, binning_scheme='generator'):
        """Construct bin edges for big TH1s

        If self.plot_with_bin_widths, returns properly sized bins
        according to lambda bins widths.
        Otherwise just returns evenly-spaced integers for uniform bin width
        """
        pt_bin_edges = self.pt_bin_edges(binning_scheme)
        num_pt_bins = self.num_pt_bins(binning_scheme)
        bins = []
        for ibin_pt, pt in enumerate(pt_bin_edges[:-1]):
            lambda_bin_edges = self.binning_handler.get_variable_bins(pt=pt,
                                                                      binning_scheme=binning_scheme)
            if self.plot_with_bin_widths:
                offset = ibin_pt * math.ceil(lambda_bin_edges[-1])
                bins.extend(lambda_bin_edges[:-1] + offset)
                if ibin_pt == (num_pt_bins - 1):
                    bins.append(lambda_bin_edges[-1] + offset)
            else:
                start = 0 if len(bins) == 0 else (bins[-1] + 1)
                end = start + len(lambda_bin_edges) - 1
                # -1 to account for first bin having lower edge of last bin
                bins.extend(list(range(start, end)))
        return array('d', bins)

    @staticmethod
    def _get_ylim(entries):
        ymax = max([e.obj.GetMaximum() for e in entries])
        ylim = (0, 3*ymax)  # room for labels
        return ylim

    def _modify_plot(self, this_plot):
        this_plot.legend.SetX1(0.72)
        this_plot.legend.SetY1(0.72)
        this_plot.legend.SetX2(0.9)
        this_plot.legend.SetY2(0.88)
        if len(this_plot.contributions) > 4:
            this_plot.legend.SetNColumns(2)
        # this_plot.left_margin = 0.16
        this_plot.left_margin = 0.07
        this_plot.right_margin = 0.038
        this_plot.default_canvas_size = (800, 500)
        this_plot.y_padding_max_linear = 2
        this_plot.subplot_pad_height = 0.4
        this_plot.subplot_line_style = 2
        this_plot.subplot_line_width = 1
        this_plot.subplot_line_color = ROOT.kGray+2
        # Good if many vertical binning lines
        this_plot.legend.SetFillColorAlpha(ROOT.kWhite, 0.9)
        this_plot.legend.SetFillStyle(1001)
        this_plot.lumi = self.setup.lumi

    def _plot_pt_bins(self, plot, binning_scheme='generator'):
        """Plot vertical pt bin lines & text labels"""
        sub_xax = plot.subplot_container.GetXaxis()
        sub_xax.SetTitleOffset(sub_xax.GetTitleOffset()*.9)
        sub_xax.SetTitleOffset(sub_xax.GetTitleOffset()*.9)
        sub_xax.SetLabelOffset(999)
        ndiv = self.num_pt_bins(binning_scheme)
        plot.container.GetXaxis().SetNdivisions(ndiv, False)
        sub_xax.SetNdivisions(ndiv, False)

        lines = []
        texts = []
        plot.main_pad.cd()
        obj = plot.container.GetHistogram()
        y_low, y_high = obj.GetMinimum(), obj.GetMaximum()

        # iterate over each pt bin
        lower_bin = 1  # ROOT bin numbers. bin num for the start of this pt bin
        for ibin, pt_val in enumerate(self.pt_bin_edges(binning_scheme)[:-1]):
            # Get Lambda bins for this pt bin
            lambda_bins = self.binning_handler.get_variable_bins(pt=pt_val,
                                                                 binning_scheme=binning_scheme)
            num_bins = len(lambda_bins) - 1
            next_lower_bin = lower_bin + num_bins  # bin num for the start of next pt bin

            # Do pt bin boundary line for end of this pt bin
            if ibin != self.num_pt_bins(binning_scheme):
                pt_bin = plot.get_modifier().GetXaxis().GetBinLowEdge(next_lower_bin)
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

            # Do bin label texts
            labels_inside_align = 'lower'

            # Get physical values on x axis of pt bin boundaries
            pt_bin_lower = plot.get_modifier().GetXaxis().GetBinLowEdge(lower_bin)
            pt_bin_higher = plot.get_modifier().GetXaxis().GetBinLowEdge(next_lower_bin)
            pt_bin_interval = pt_bin_higher - pt_bin_lower
            text_x = pt_bin_lower + 0.35*(pt_bin_interval)

            # figure out y location from axes
            axis_range = y_high - y_low
            # assumes logarithmic?
            if ROOT.gPad.GetLogy():
                if y_low <= 0:
                    print("y_low is %f so can't take log" % y_low)
                    return None, None
                    # raise ValueError("y_low is %f so can't take log" %y_low)
                log_axis_range = math.log10(y_high) - math.log10(y_low)
                y_start = math.pow(10, 0.03*log_axis_range + math.log10(y_low))
            else:
                y_start = (0.3*axis_range) + y_low

            if labels_inside_align == 'higher':
                if ROOT.gPad.GetLogy():
                    y_start = y_high/10
                else:
                    y_start = y_high*0.8

            text = ROOT.TPaveText(text_x, y_start, text_x + .5*pt_bin_interval, y_start)  # end y determined by text size
            text.SetBorderSize(0)
            text.SetFillStyle(0)
            contents = '%g < p_{T} < %g GeV' % (pt_val, self.pt_bin_edges(binning_scheme)[ibin+1])
            # You have to style the thing that is returned by AddText, not the TPaveText obj itself
            # WTAF
            t = text.AddText(contents)  # t is a TText
            t.SetTextColor(14)
            t.SetTextAngle(89)
            # t.SetTextSize(0.0275)
            t.SetTextSize(0.0375)
            if (isinstance(plot, Plot) and not plot.subplot_pad) or not isinstance(plot, Plot):
                t.SetTextSize(0.02)  # account for the fact that a subplot makes things smaller
            if labels_inside_align == 'lower':
                t.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignTop)
            else:
                t.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignTop)
            text.Draw()
            texts.append(text)

            lower_bin = next_lower_bin

        return lines, texts

    def make_big_1d_normalised_hist(self, name, binning_scheme='generator', hist_bin_chopper=None, do_div_bin_width=True):
        """Make big 1D plot with normalised distribution per pt bin"""
        hist_bin_chopper = hist_bin_chopper or self.hist_bin_chopper
        h_new = ROOT.TH1D(name + "_1d_all_" + cu.get_unique_str(), "", self.nbins(binning_scheme), self.all_bin_edges(binning_scheme))
        global_bin = 1
        for ibin, _ in enumerate(self.pt_bin_edges(binning_scheme)[:-1]):
            # hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width(name, ibin, binning_scheme=binning_scheme)
            hist_bin = hist_bin_chopper.get_bin_plot(name, ibin,
                                                     axis='pt',
                                                     do_norm=True,
                                                     do_div_bin_width=do_div_bin_width,
                                                     binning_scheme=binning_scheme)

            for lbin in range(1, hist_bin.GetNbinsX()+1):
                h_new.SetBinContent(global_bin, hist_bin.GetBinContent(lbin))
                h_new.SetBinError(global_bin, hist_bin.GetBinError(lbin))
                global_bin += 1

        return h_new

    def get_mc_truth_kwargs(self):
        return dict(
            label=self.setup.region['mc_label'],
            line_color=PLOT_STYLES['gen_colour'], line_width=self.line_width,
            marker_color=PLOT_STYLES['gen_colour'], marker_size=0
            )

    def get_unfolded_total_err_kwargs(self):
        return dict(
            label="Data" if self.setup.has_data else "Unfolded MC",
            line_color=PLOT_STYLES['unfolded_total_colour'], line_width=self.line_width, line_style=1,
            marker_color=PLOT_STYLES['unfolded_total_colour']
        )

    def get_unfolded_stat_err_kwargs(self):
        src_type = "data" if self.setup.has_data else "MC"
        tau_str = "" if self.unfolder.tau == 0 else "(#tau = %.3g) " % self.unfolder.tau
        return dict(
            label="Unfolded %s %s(%s response matrix)" % (src_type, tau_str, self.setup.region['mc_label']),
            line_color=PLOT_STYLES['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
            marker_color=PLOT_STYLES['unfolded_stat_colour'], marker_style=cu.Marker.get('circle'), marker_size=0
        )

    def get_alt_mc_truth_kwargs(self):
        return dict(
            label=self.setup.region['alt_mc_label'],
            line_color=PLOT_STYLES['alt_gen_colour'], line_width=self.line_width, line_style=1,  # hard to see other line styles
            marker_color=PLOT_STYLES['alt_gen_colour'], marker_size=0  # hard to make small enough markers
        )

    def get_alt_unfolded_stat_err_kwargs(self, alt_unfolder):
        tau_str = "" if alt_unfolder.tau == 0 else "(#tau = %.3g) " % alt_unfolder.tau
        return dict(
            label="Unfolded MC %s(%s response matrix)" % (tau_str, self.setup.region['alt_mc_label']),
            line_color=PLOT_STYLES['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
            marker_color=PLOT_STYLES['alt_unfolded_colour']
        )

    def get_title(self):
        return (("{jet_algo}\n"
                 "{region_label} region")
                     .format(
                        jet_algo=self.setup.jet_algo,
                        region_label=self.setup.region['label'],
                        pt_str=self.setup.pt_var_str)
               )

    def get_ytitle(self):
        return self.setup.pt_bin_normalised_differential_label if self.plot_with_bin_widths else self.setup.pt_bin_normalised_differential_times_width_label

    def get_subplot_ylim(self):
        return (0, 2) if self.setup.has_data else (0.7, 1.3)

    def get_big_1d(self, name, binning_scheme='generator', hist_bin_chopper=None, hist_bin_chopper_key=None, do_div_bin_width=True):
        """Getter for big 1D hist - use cached version if available, otherwise create and store in cache

        By default uses the main unfolder.hist_bin_chopper to construct the hists
        using histogram cached name "name". Can also specify an alternate
        hist_bin_chopper, as well as a key for it (default is same as name arg)
        This allows you to do e.g. name="alt_unfolded", hist_bin_chopper_key="unfolded"
        to store the unfolded hist from alt_unfolder in this cache as "alt_unfolded"
        """
        key = name + '_' + binning_scheme + ("_divBinWidth" if do_div_bin_width else "")
        if self._cache_1d.get(key, None) is None:
            hist_bin_chopper_key = hist_bin_chopper_key or name
            self._cache_1d[key] = self.make_big_1d_normalised_hist(hist_bin_chopper_key,
                                                                   binning_scheme=binning_scheme,
                                                                   hist_bin_chopper=hist_bin_chopper,
                                                                   do_div_bin_width=do_div_bin_width)
        return self._cache_1d[key]

    def plot_unfolded_truth(self):
        entries = [
            Contribution(self.get_big_1d('hist_truth', 'generator'),
                         **self.get_mc_truth_kwargs()),
            Contribution(self.get_big_1d('unfolded', 'generator'),
                         subplot=self.get_big_1d('hist_truth', 'generator'),
                         **self.get_unfolded_total_err_kwargs()),
            Contribution(self.get_big_1d('unfolded_stat_err', 'generator'),
                         subplot=self.get_big_1d('hist_truth', 'generator'),
                         **self.get_unfolded_stat_err_kwargs()),
        ]
        plot = Plot(entries, 'hist',
                    ytitle=self.get_ytitle(),
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
        plot.legend.SetX1(0.4)
        plot.legend.SetY1(0.8)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_reco_input(self):
        entries = []
        self.hist_bin_chopper.add_obj("input_hist_bg_subtracted", self.unfolder.input_hist_bg_subtracted)
        self.hist_bin_chopper.add_obj("hist_mc_reco_bg_subtracted", self.unfolder.hist_mc_reco_bg_subtracted)
        entries.append(
            Contribution(self.get_big_1d("input_hist_bg_subtracted", 'detector'),
                         label="Data bg-subtracted [detector-level]",
                         line_color=ROOT.kBlack, line_width=0,
                         marker_color=ROOT.kBlack, marker_size=0.6, marker_style=cu.Marker.get('circle'),
                         normalise_hist=False, subplot=self.get_big_1d("hist_mc_reco_bg_subtracted", 'detector'),
                         subplot_line_width=1),
        )

        entries.append(
            Contribution(self.get_big_1d("hist_mc_reco_bg_subtracted", 'detector'),
                         label="MC bg-subtracted [detector-level]",
                         line_color=ROOT.kAzure+4, line_width=1,
                         marker_color=ROOT.kAzure+4, marker_size=0,
                         normalise_hist=False),
        )

        plot = Plot(entries, 'hist',
                    ytitle=self.get_ytitle(),
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(entries),
                    subplot_type='ratio',
                    subplot_title="Data / MC",
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot, binning_scheme='detector')
        plot.save("%s/detector_1d_normalised_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_truth_alt_truth(self):
        data = self.get_big_1d('unfolded', 'generator')
        data_no_errors = data.Clone()
        cu.remove_th1_errors(data_no_errors)
        entries = [
            Contribution(self.get_big_1d('alt_hist_truth', 'generator'),
                         subplot=data_no_errors,
                         **self.get_alt_mc_truth_kwargs()),
            Contribution(self.get_big_1d('hist_truth', 'generator'),
                         subplot=data_no_errors,
                         **self.get_mc_truth_kwargs()),
            Contribution(data,
                         **self.get_unfolded_total_err_kwargs()),
        ]
        plot = Plot(entries, 'hist',
                    ytitle=self.get_ytitle(),
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(entries),
                    subplot_type='ratio',
                    # subplot_title="#splitline{* / Gen}{(%s)}" % (self.setup.region['mc_label']),
                    subplot_title="Simulation / data",
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.reverse_legend = True
        plot.legend.SetNColumns(3)
        plot.legend.SetX1(0.4)
        plot.legend.SetY1(0.8)
        subplot_draw_opts = "NOSTACK E"
        plot.plot("NOSTACK E", subplot_draw_opts)

        # on subplot, plot shaded region for data uncertainties
        data_total_ratio = data_no_errors.Clone()
        data_total_ratio.Divide(data)

        data_total_ratio.SetFillStyle(3354)
        # data_total_ratio.SetFillColor(PLOT_STYLES['unfolded_total_colour'])
        data_total_ratio.SetFillColor(16) # not black, too dark
        data_total_ratio.SetLineWidth(0)
        data_total_ratio.SetMarkerSize(0)

        # now draw the data error shaded area
        # this is a bit hacky - basically draw them on the ratio pad,
        # then redraw the existing hists & line to get them ontop
        # note that we use "same" for all - this is to keep the original axes
        # (we may want to rethink this later?)
        plot.subplot_pad.cd()
        draw_opt = "E2 SAME"
        plot.subplot_line.Draw()
        data_total_ratio.Draw(draw_opt)
        plot.subplot_container.Draw("SAME" + subplot_draw_opts)

        l, t = self._plot_pt_bins(plot)

        plot.save("%s/unfolded_1d_normalised_alt_truth_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_truth_alt_response(self, alt_unfolder):
        entries = [
            Contribution(self.get_big_1d('hist_truth', 'generator'),
                         **self.get_mc_truth_kwargs()),
            Contribution(self.get_big_1d('unfolded_stat_err', 'generator'),
                         subplot=self.get_big_1d('hist_truth', 'generator'),
                         **self.get_unfolded_stat_err_kwargs()),
            Contribution(self.get_big_1d('alt_unfolded_stat_err', 'generator',
                                         hist_bin_chopper=alt_unfolder.hist_bin_chopper,
                                         hist_bin_chopper_key='unfolded_stat_err'),
                         subplot=self.get_big_1d('hist_truth', 'generator'),
                         **self.get_alt_unfolded_stat_err_kwargs(alt_unfolder)),
        ]
        plot = Plot(entries, 'hist',
                    ytitle=self.get_ytitle(),
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(entries),
                    subplot_type='ratio',
                    subplot_title="* / %s" % (self.setup.region['mc_label']),
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.legend.SetNColumns(2)
        plot.legend.SetX1(0.35)
        plot.legend.SetY1(0.8)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_alt_response_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_truth_alt_response_truth(self, alt_unfolder):
        entries = [
            Contribution(self.get_big_1d('hist_truth', 'generator'),
                         **self.get_mc_truth_kwargs()),
            Contribution(self.get_big_1d('alt_hist_truth', 'generator'),
                         subplot=self.get_big_1d('hist_truth', 'generator'),
                         **self.get_alt_mc_truth_kwargs()),
            Contribution(self.get_big_1d('unfolded_stat_err', 'generator'),
                         subplot=self.get_big_1d('hist_truth', 'generator'),
                         **self.get_unfolded_stat_err_kwargs()),
            Contribution(self.get_big_1d('alt_unfolded_stat_err', 'generator',
                                         hist_bin_chopper=alt_unfolder.hist_bin_chopper,
                                         hist_bin_chopper_key='unfolded_stat_err'),
                         subplot=self.get_big_1d('hist_truth', 'generator'),
                         **self.get_alt_unfolded_stat_err_kwargs(alt_unfolder)),
        ]
        plot = Plot(entries, 'hist',
                    ytitle=self.get_ytitle(),
                    title=self.get_title(),
                    xtitle=self.setup.angle_str + ', per %s bin' % (self.setup.pt_var_str),
                    has_data=self.setup.has_data,
                    ylim=self._get_ylim(entries),
                    subplot_type='ratio',
                    subplot_title="* / %s" % (self.setup.region['mc_label']),
                    subplot_limits=self.get_subplot_ylim()
                    )
        self._modify_plot(plot)
        plot.legend.SetNColumns(2)
        plot.legend.SetX1(0.35)
        plot.legend.SetY1(0.8)
        plot.plot("NOSTACK E")
        l, t = self._plot_pt_bins(plot)
        plot.save("%s/unfolded_1d_normalised_alt_response_truth_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_exp_systs(self):
        for syst_dict in self.setup.region['experimental_systematics']:
            entries = [
                        # Contribution(self.get_big_1d('hist_truth', 'generator'),
                        #  **self.get_mc_truth_kwargs()),
                        Contribution(self.get_big_1d('unfolded_stat_err', 'generator'),
                                     **dict(self.get_unfolded_stat_err_kwargs(),
                                            line_color=ROOT.kGray+2))
            ]

            this_syst = self.unfolder.get_exp_syst(syst_dict['label'])

            # syst_label_no_spaces = cu.no_space_str()
            # hbc_name = 'syst_shifted_%s_unfolded' % syst_label_no_spaces
            self.hist_bin_chopper.add_obj(this_syst.syst_shifted_label, this_syst.syst_shifted)
            entries.append(Contribution(self.get_big_1d(this_syst.syst_shifted_label, 'generator'),
                                        label=syst_dict['label'],
                                        line_color=syst_dict['colour'], line_width=self.line_width,
                                        line_style=2 if 'down' in syst_dict['label'].lower() else 1,
                                        marker_color=syst_dict['colour'], marker_size=0,
                                        subplot=self.get_big_1d('unfolded_stat_err', 'generator')))
            plot = Plot(entries, 'hist',
                        ytitle=self.get_ytitle(),
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
            plot.save("%s/unfolded_1d_normalised_exp_syst_%s_%s_divBinWidth.%s" % (self.setup.output_dir, this_syst.label_no_spaces, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_scale_systs(self):
        all_entries = [Contribution(self.get_big_1d('hist_truth', 'generator'),
                                     **self.get_mc_truth_kwargs()),
                       Contribution(self.get_big_1d('unfolded_stat_err', 'generator'),
                                    subplot=self.get_big_1d('hist_truth', 'generator'),
                                    **dict(self.get_unfolded_stat_err_kwargs(),
                                           line_color=ROOT.kGray+2))
                      ]

        for syst_dict in self.setup.region['scale_systematics']:
            syst_unfolder = syst_dict['unfolder']
            syst_label = syst_dict['label']

            syst_label_no_spaces = cu.no_space_str(syst_dict['label'])

            hbc_name = 'scale_syst_%s_unfolded' % (syst_label_no_spaces)
            self.hist_bin_chopper.add_obj(hbc_name, syst_unfolder.unfolded)

            all_entries.extend([
                Contribution(self.get_big_1d(hbc_name, 'generator'),
                             label="Unfolded (#tau = %.3g) (stat. unc.) (%s)" % (syst_unfolder.tau, syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                             marker_color=syst_dict['colour'], marker_size=0,
                             subplot=self.get_big_1d('hist_truth', 'generator')),
            ])

        plot = Plot(all_entries, 'hist',
                    ytitle=self.get_ytitle(),
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
        plot.save("%s/unfolded_1d_normalised_scale_scale_syst_%s_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, self.setup.output_fmt))

    def plot_unfolded_model_systs(self):
        all_entries = [Contribution(self.get_big_1d('hist_truth', 'generator'),
                                    line_style=2,
                                    **self.get_mc_truth_kwargs()),
                       Contribution(self.get_big_1d('unfolded_stat_err', 'generator'),
                                    subplot=self.get_big_1d('hist_truth', 'generator'),
                                    **dict(self.get_unfolded_stat_err_kwargs(),
                                           line_color=ROOT.kGray+2))
                      ]

        for syst_dict in self.setup.region['model_systematics']:
            entries = [
                        Contribution(self.get_big_1d('hist_truth', 'generator'),
                                     line_style=2,
                                     **self.get_mc_truth_kwargs()),
                        Contribution(self.get_big_1d('unfolded_stat_err', 'generator'),
                                     subplot=self.get_big_1d('hist_truth', 'generator'),
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
                Contribution(self.get_big_1d(hbc_name_gen, 'generator'),
                             label="Generator (%s)" % (syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                             marker_color=syst_dict['colour'], marker_size=0),
                Contribution(self.get_big_1d(hbc_name, 'generator'),
                             label="Unfolded (#tau = %.3g) (stat. unc.) (%s)" % (syst_unfolder.tau, syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                             marker_color=syst_dict['colour'], marker_size=0,
                             subplot=self.get_big_1d(hbc_name_gen, 'generator')),
            ])
            all_entries.extend(entries[-2:])
            plot = Plot(entries, 'hist',
                        ytitle=self.get_ytitle(),
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
                    ytitle=self.get_ytitle(),
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
        all_entries = [Contribution(self.get_big_1d('hist_truth', 'generator'),
                                    line_style=2,
                                    **self.get_mc_truth_kwargs()),
                       Contribution(self.get_big_1d('unfolded_stat_err', 'generator'),
                                    subplot=self.get_big_1d('hist_truth', 'generator'),
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
                Contribution(self.get_big_1d(hbc_name_gen, 'generator'),
                             label="Generator (%s)" % (syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                             marker_color=syst_dict['colour'], marker_size=0),
                Contribution(self.get_big_1d(hbc_name, 'generator'),
                             label="Unfolded (#tau = %.3g) (stat. unc.) (%s)" % (syst_unfolder.tau, syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                             marker_color=syst_dict['colour'], marker_size=0,
                             subplot=self.get_big_1d(hbc_name_gen, 'generator')),
            ])

        plot = Plot(all_entries, 'hist',
                    ytitle=self.get_ytitle(),
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
                    Contribution(self.get_big_1d('hist_truth', 'generator'),
                                 **self.get_mc_truth_kwargs()),
                    Contribution(self.get_big_1d('unfolded_stat_err', 'generator'),
                                 subplot=self.get_big_1d('hist_truth', 'generator'),
                                 **dict(self.get_unfolded_stat_err_kwargs(),
                                        line_color=ROOT.kGray+2))
        ]
        for syst_dict in self.setup.region['pdf_systematics']:
            syst_unfolder = syst_dict['unfolder']
            syst_label = syst_dict['label']
            syst_label_no_spaces = cu.no_space_str(syst_dict['label'])

            hbc_name_gen = 'pdf_syst_%s_hist_truth' % (syst_label_no_spaces)
            self.hist_bin_chopper.add_obj(hbc_name_gen, syst_unfolder.hist_truth)

            hbc_name = 'pdf_syst_%s_unfolded' % (syst_label_no_spaces)
            self.hist_bin_chopper.add_obj(hbc_name, syst_unfolder.unfolded)

            entries.extend([
                Contribution(self.get_big_1d(hbc_name_gen, 'generator'),
                             label="Generator (%s)" % (syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                             marker_color=syst_dict['colour'], marker_size=0),
                Contribution(self.get_big_1d(hbc_name, 'generator'),
                             label="Unfolded (#tau = %.3g) (stat. unc.) (%s)" % (syst_unfolder.tau, syst_label),
                             line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                             marker_color=syst_dict['colour'], marker_size=0,
                             subplot=self.get_big_1d(hbc_name_gen, 'generator')),
            ])
        plot = Plot(entries, 'hist',
                    ytitle=self.get_ytitle(),
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
    #         Contribution(self.get_big_1d("hist_mc_reco_bg_subtracted", 'generator'),
    #                      label=""),
    #         Contribution(self.get_big_1d("folded_mc_truth", 'generator'),
    #                      label="")
    #     ]
    #     plot = Plot(entries, 'hist',
    #                 ytitle=self.get_ytitle(),
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


@profile
def do_all_big_normalised_1d_plots_per_region_angle(setup, hist_bin_chopper=None):
    """Do various big 1D plots collating all normalised distributions"""
    region = setup.region
    unfolder = region['unfolder']
    alt_unfolder = region.get('alt_unfolder', None)
    alt_hist_truth = region.get('alt_hist_mc_gen', None)

    if not hist_bin_chopper:
        # unreg_unfolder = unpack_dict['unreg_unfolder']
        # hist_bin_chopper = HistBinChopper(generator_binning=unfolder.generator_binning.FindNode("generatordistribution"),
        #                                   detector_binning=unfolder.detector_binning.FindNode("detectordistribution"))
        # hist_bin_chopper.add_obj("hist_truth", unfolder.hist_truth)
        hist_bin_chopper = unfolder.hist_bin_chopper
        if alt_unfolder:
            hist_bin_chopper.add_obj("alt_unfolded_stat_err", alt_unfolder.unfolded_stat_err)
        if alt_hist_truth:
            hist_bin_chopper.add_obj("alt_hist_truth", alt_hist_truth)

    has_exp_systs = len(region['experimental_systematics']) > 0
    has_scale_systs = len(region['scale_systematics']) > 0
    has_model_systs = len(region['model_systematics']) > 0
    has_pdf_systs = len(region['pdf_systematics']) > 0

    if has_exp_systs: print("We have experimental systs")
    if has_model_systs: print("We have model systs")
    if has_scale_systs: print("We have scale systs")
    if has_pdf_systs: print("We have pdf systs")
    # if has_jackknife_input_vars: print("We have jackknife input variations")
    # if has_jackknife_response_vars: print("We have jackknife response variations")

    big_plotter = BigNormalised1DPlotter(setup, hist_bin_chopper, plot_with_bin_widths=False)

    print("...doing standard big 1D plots")
    big_plotter.plot_unfolded_truth()
    big_plotter.plot_reco_input()

    if alt_hist_truth:
        print("...doing alt truth big 1D plots")
        big_plotter.plot_unfolded_truth_alt_truth()

    if alt_unfolder:
        print("...doing alt unfolder big 1D plots")
        big_plotter.plot_unfolded_truth_alt_response(alt_unfolder)
        big_plotter.plot_unfolded_truth_alt_response_truth(alt_unfolder)

    # if has_exp_systs:
    #     print("...doing exp systs big 1D plots")
    #     big_plotter.plot_unfolded_exp_systs()

    # if has_scale_systs:
    #     print("...doing scale systs big 1D plots")
    #     big_plotter.plot_unfolded_scale_systs()

    # if has_model_systs:
    #     print("...doing model systs big 1D plots")
    #     big_plotter.plot_unfolded_model_systs()
    #     big_plotter.plot_unfolded_model_systs_only_scale()

    # if has_pdf_systs:
    #     print("...doing pdf systs big 1D plots")
    #     big_plotter.plot_unfolded_pdf_systs()


@profile
def do_all_big_absolute_1d_plots_per_region_angle(setup):
    """Do various big absolute 1d plots"""
    region = setup.region
    unfolder = region['unfolder']
    alt_unfolder = region.get('alt_unfolder', None)
    alt_hist_truth = region.get('alt_hist_mc_gen', None)

    has_exp_systs = len(region['experimental_systematics']) > 0
    has_scale_systs = len(region['scale_systematics']) > 0
    has_model_systs = len(region['model_systematics']) > 0
    has_pdf_systs = len(region['pdf_systematics']) > 0

    if has_exp_systs: print("We have experimental systs")
    if has_model_systs: print("We have model systs")
    if has_scale_systs: print("We have scale systs")
    if has_pdf_systs: print("We have pdf systs")
    # if has_jackknife_input_vars: print("We have jackknife input variations")
    # if has_jackknife_response_vars: print("We have jackknife response variations")

    unfolder_plotter = MyUnfolderPlotter(unfolder,
                                         is_data=setup.has_data,
                                         lumi=setup.lumi)

    print("...doing standard big absolute1D plots")

    title = "%s\n%s region, %s" % (setup.jet_algo, region['label'], setup.angle_str)
    # put it in the over-arching directory as a quick check
    unfolder_plotter.draw_unfolded_1d(output_dir=os.path.dirname(os.path.dirname(setup.output_dir)),
                                      append=setup.append,
                                      title=title)

    # put another (without arrows) in usual plot dir
    unfolder_plotter.draw_unfolded_1d(output_dir=setup.output_dir,
                                      mark_negatives=False,
                                      append=setup.append,
                                      title=title)

    # put another (without arrows) + alternate in usual plot dir
    if alt_unfolder is not None:
        unfolder_plotter.draw_unfolded_1d(output_dir=setup.output_dir,
                                          mark_negatives=False,
                                          other_contributions=[
                                              Contribution(alt_unfolder.get_output(),
                                                           label="Alt MC",
                                                           line_color=ROOT.kGreen+2,
                                                           marker_color=ROOT.kGreen+2,
                                                           subplot=unfolder.hist_truth)
                                          ],
                                          append=setup.append + "_altMC",
                                          title=title)

    if alt_hist_truth is not None:
        unfolder_plotter.draw_unfolded_1d(output_dir=setup.output_dir,
                                          mark_negatives=False,
                                          other_contributions=[
                                              Contribution(alt_hist_truth,
                                                           label="Alt MC",
                                                           line_color=ROOT.kGreen+2,
                                                           marker_color=ROOT.kGreen+2,
                                                           subplot=unfolder.hist_truth)
                                          ],
                                          append=setup.append + "_altMCTruth",
                                          title=title)

    # reco using detector binning
    unfolder_plotter.draw_detector_1d(do_reco_mc=True,
                                      do_reco_data=setup.has_data,
                                      do_reco_bg=True,
                                      output_dir=setup.output_dir,
                                      append=setup.append,
                                      title=title)

    # same plot but with bg-subtracted reco (incl fakes)
    unfolder_plotter.draw_detector_1d(do_reco_data_bg_sub=setup.has_data,
                                      do_reco_bg=True,
                                      do_reco_mc_bg_sub=True,
                                      output_dir=setup.output_dir,
                                      append='bg_fakes_subtracted_%s' % setup.append,
                                      title=title)

    # with alt MC
    if alt_hist_truth is not None:
        unfolder_plotter.draw_detector_1d(do_reco_data_bg_sub=setup.has_data,
                                          do_reco_bg=False,
                                          do_reco_mc_bg_sub=True,
                                          output_dir=setup.output_dir,
                                          other_contributions=[
                                                  Contribution(region['alt_hist_mc_reco_bg_subtracted'],
                                                               label="Alt MC bg-subtracted [detector-level]",
                                                               line_color=ROOT.kGreen+2,
                                                               marker_color=ROOT.kGreen+2,
                                                               subplot=None)
                                              ],
                                          append='bg_fakes_subtracted_altMC_%s' % setup.append,
                                          title=title)


    # reco using gen binning
    unfolder_plotter.draw_generator_1d(do_reco_data=setup.has_data,
                                       do_reco_data_bg_sub=False,
                                       do_reco_bg=False,
                                       do_reco_mc=True,
                                       do_reco_mc_bg_sub=False,
                                       do_truth_mc=True,
                                       output_dir=setup.output_dir,
                                       append=setup.append,
                                       title=title)

    # same but with bg-subtracted generator-binning
    unfolder_plotter.draw_generator_1d(do_reco_data=False,
                                       do_reco_data_bg_sub=setup.has_data,
                                       do_reco_bg=True,
                                       do_reco_mc=False,
                                       do_reco_mc_bg_sub=True,
                                       do_truth_mc=False,
                                       output_dir=setup.output_dir,
                                       append='bg_fakes_subtracted_%s' % setup.append,
                                       title=title)

    plot_args = dict(output_dir=setup.output_dir, append=setup.append)

    unfolder_plotter.draw_rel_uncertainty_1d(draw_detector=True,
                                             draw_unfolded=True,
                                             title=title,
                                             **plot_args)

    # Folded unfolded & folded truth
    unfolder_plotter.draw_truth_unfolded_folded(title=title, **plot_args)

    rsp_title = "Response matrix, %s, %s region, %s" % (setup.jet_algo, region['label'], setup.angle_str)
    unfolder_plotter.draw_response_matrix(title=rsp_title, **plot_args)

    prob_title = ("#splitline{Probability matrix, %s, %s region, %s}{Condition number: #sigma_{max} / #sigma_{min} = %.3g / %.3g = %g}"
                % (setup.jet_algo, region['label'], setup.angle_str, unfolder.sigma_max, unfolder.sigma_min, unfolder.condition_number))
    unfolder_plotter.draw_probability_matrix(title=prob_title, **plot_args)

    if has_scale_systs:
        # detector scale variations - useless as all unfold same input
        # scale_contributions = [
        #     Contribution(mdict['unfolder'].hist_mc_reco,
        #                  label=mdict['label'],
        #                  line_color=mdict['colour'], line_style=1, line_width=1,
        #                  marker_color=mdict['colour'], marker_size=0, marker_style=21,
        #                  subplot=unfolder.hist_mc_reco)
        #     for mdict in region['scale_systematics']
        # ]
        # unfolder_plotter.draw_detector_1d(do_reco_mc=True,
        #                                   output_dir=setup.output_dir,
        #                                   append='scale_systs_%s' % setup.append,
        #                                   title=title,
        #                                   other_contributions=scale_contributions,
        #                                   subplot_title='#splitline{Variation /}{nominal}')
        # unfolded scale variations
        scale_contributions = [
            Contribution(mdict['unfolder'].get_unfolded_with_ematrix_stat(),
                         label=mdict['label'],
                         line_color=mdict['colour'], line_style=1, line_width=1,
                         marker_color=mdict['colour'], marker_size=0, marker_style=21,
                         subplot=unfolder.get_unfolded_with_ematrix_stat())
            for mdict in region['scale_systematics']
        ]
        unfolder_plotter.draw_unfolded_1d(do_unfolded=True,
                                          do_gen=False,
                                          output_dir=setup.output_dir,
                                          append='scale_systs_%s' % setup.append,
                                          title=title,
                                          other_contributions=scale_contributions,
                                          subplot_limits=(0.975, 1.025),
                                          mark_negatives=False,
                                          subplot_title='#splitline{Variation /}{nominal}')
        # truth scale variations
        # useless as all use same truth as same input
        # scale_contributions = [
        #     Contribution(mdict['unfolder'].hist_truth,
        #                  label=mdict['label'],
        #                  line_color=mdict['colour'], line_style=1, line_width=1,
        #                  marker_color=mdict['colour'], marker_size=0, marker_style=21,
        #                  subplot=unfolder.hist_truth)
        #     for mdict in region['scale_systematics']
        # ]
        # unfolder_plotter.draw_generator_1d(do_truth_mc=True,
        #                                    output_dir=setup.output_dir,
        #                                    append='scale_systs_%s' % setup.append,
        #                                    title=title,
        #                                    other_contributions=scale_contributions,
        #                                    subplot_title='#splitline{Variation /}{nominal}')


@profile
def get_bottom_line_stats(setup):
    """Construct dict of bottom-line (i.e. chi2) stats for this region/angle combo"""
    unfolder = setup.region['unfolder']

    # ABOSULTE VERSION
    # Not sure if this makes sense? Since subject to overall normalisation issues
    # But we can't do it on the Jacobian-transformed normalised version,
    # since the Jacobian is singular (i.e. non-invertible)

    # TODO: convert detector space to gen binning to match NdF, improve chi2

    # do smeared chi2
    # folded_mc_truth_gen_binning = unfolder.get_ndarray_signal_region_no_overflow(
    #                                   unfolder.convert_reco_binned_hist_to_gen_binned(
    #                                       unfolder.get_folded_mc_truth()
    #                                   ),
    #                                   xbinning='generator',
    #                                   ybinning=None
    #                               )
    # folded_alt_mc_truth_gen_binning = unfolder.get_ndarray_signal_region_no_overflow(
    #                                       unfolder.convert_reco_binned_hist_to_gen_binned(
    #                                         unfolder.fold_generator_level(setup.region['alt_hist_mc_gen'])
    #                                       ),
    #                                       xbinning='generator',
    #                                       ybinning=None
    #                                   )
    # input_hist_gen_binning_bg_subtracted_signal = unfolder.get_ndarray_signal_region_no_overflow(
    #                                                   unfolder.input_hist_gen_binning_bg_subtracted,
    #                                                   xbinning='generator',
    #                                                   ybinning=None
    #                                               )

    vyy = unfolder.make_diag_cov_hist_from_errors(
              h1d=unfolder.input_hist_gen_binning_bg_subtracted,
              inverse=False
          )
    # vyy_gen_binning = unfolder.get_ndarray_signal_region_no_overflow(
    #                       vyy,
    #                       xbinning='generator',
    #                       ybinning='generator'
    #                   )

    vyy_inv = unfolder.make_diag_cov_hist_from_errors(
                  h1d=unfolder.input_hist_gen_binning_bg_subtracted,
                  inverse=True
              )
    # vyy_inv_gen_binning = unfolder.get_ndarray_signal_region_no_overflow(
    #                           vyy_inv,
    #                           xbinning='generator',
    #                           ybinning='generator'
    #                       )

    folded_mc_truth_gen_binning = unfolder.convert_reco_binned_hist_to_gen_binned(unfolder.get_folded_mc_truth())
    folded_alt_mc_truth_gen_binning = unfolder.convert_reco_binned_hist_to_gen_binned(unfolder.fold_generator_level(setup.region['alt_hist_mc_gen']))
    input_hist_gen_binning_bg_subtracted_signal = unfolder.input_hist_gen_binning_bg_subtracted
    vyy_gen_binning = vyy
    vyy_inv_gen_binning = vyy_inv

    smeared_chi2, smeared_ndf, smeared_p = unfolder.calculate_chi2(one_hist=folded_mc_truth_gen_binning,
                                                                   # one_hist=unfolder.get_folded_mc_truth(),
                                                                   # other_hist=unfolder.input_hist_bg_subtracted,
                                                                   other_hist=input_hist_gen_binning_bg_subtracted_signal,
                                                                   # cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
                                                                   # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_subtraction_ndarray(),
                                                                   cov_inv_matrix=vyy_inv_gen_binning,
                                                                   cov_matrix=vyy_gen_binning,
                                                                   detector_space=False,
                                                                   ignore_underflow_bins=False,
                                                                   has_underflow=False,
                                                                   debugging_dir=os.path.join(setup.output_dir, "smeared_chi2_gen_binning_signal"))
                                                                   # debugging_dir=None)
    print('smeared chi2, mdf, chi2/ndf, p:', smeared_chi2, smeared_ndf, smeared_chi2/smeared_ndf, smeared_p)

    # do for alt model
    smeared_alt_chi2, smeared_alt_ndf, smeared_alt_p = unfolder.calculate_chi2(one_hist=folded_alt_mc_truth_gen_binning,
                                                                               # one_hist=unfolder.fold_generator_level(setup.region['alt_hist_mc_gen']),
                                                                               # other_hist=unfolder.input_hist_bg_subtracted,
                                                                               other_hist=input_hist_gen_binning_bg_subtracted_signal,
                                                                               # cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
                                                                               # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_subtraction_ndarray(),
                                                                               cov_inv_matrix=vyy_inv_gen_binning,
                                                                               detector_space=False,
                                                                               ignore_underflow_bins=False,
                                                                               has_underflow=False,
                                                                               debugging_dir=os.path.join(setup.output_dir, "smeared_alt_chi2_gen_binning_signal"))
                                                                               # debugging_dir=None)
    print('smeared alt chi2, mdf, chi2/ndf, p:', smeared_alt_chi2, smeared_alt_ndf, smeared_alt_chi2/smeared_alt_ndf, smeared_alt_p)

    # unfolded_signal = unfolder.get_ndarray_signal_region_no_overflow(
    #                       unfolder.unfolded,
    #                       xbinning='generator',
    #                       ybinning=None
    #                   )
    # hist_truth_signal = unfolder.get_ndarray_signal_region_no_overflow(
    #                         unfolder.hist_truth,
    #                         xbinning='generator',
    #                         ybinning=None
    #                     )
    # hist_alt_truth_signal = unfolder.get_ndarray_signal_region_no_overflow(
    #                         setup.region['alt_hist_mc_gen'],
    #                         xbinning='generator',
    #                         ybinning=None
    #                     )
    # vxx_total_signal = unfolder.get_ndarray_signal_region_no_overflow(
    #                        # unfolder.get_ematrix_tunfold_total(),
    #                        unfolder.get_vxx_ndarray(),  # tunfold's version
    #                        xbinning='generator',
    #                        ybinning='generator'
    #                    )

    # vxx_total_signal_inv = np.linalg.pinv(vxx_total_signal, 1E-200)  # no need to cut down to signal region

    # tunfold's version
    # vxx_total_signal_inv = unfolder.get_ndarray_signal_region_no_overflow(
    #                           unfolder.get_vxx_inv_ndarray(),
    #                           xbinning='generator',
    #                           ybinning='generator'
    #                       )

    unfolded_signal = unfolder.unfolded
    hist_truth_signal = unfolder.hist_truth
    hist_alt_truth_signal = setup.region['alt_hist_mc_gen']
    vxx_total_signal = unfolder.get_ematrix_tunfold_total()
    vxx_total_signal_inv = unfolder.get_vxx_inv_ndarray()

    # do unfolded chi2
    unfolded_chi2, unfolded_ndf, unfolded_p = unfolder.calculate_chi2(one_hist=unfolded_signal,
                                                                      other_hist=hist_truth_signal,
                                                                      cov_inv_matrix=vxx_total_signal_inv,
                                                                      cov_matrix=vxx_total_signal,
                                                                      detector_space=False,
                                                                      ignore_underflow_bins=False,
                                                                      has_underflow=False,
                                                                      debugging_dir=os.path.join(setup.output_dir, "unfolded_chi2_signal_full_inv"))
                                                                      # debugging_dir=None)
    print('unfolded chi2, mdf, chi2/ndf, p:', unfolded_chi2, unfolded_ndf, unfolded_chi2/unfolded_ndf, unfolded_p)

    # do for alt model
    unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_p = unfolder.calculate_chi2(one_hist=unfolded_signal,
                                                                                  other_hist=hist_alt_truth_signal,
                                                                                  cov_inv_matrix=vxx_total_signal_inv,
                                                                                  detector_space=False,
                                                                                  ignore_underflow_bins=False,
                                                                                  has_underflow=False,
                                                                                  debugging_dir=os.path.join(setup.output_dir, "unfolded_alt_chi2_signal_full_inv"))
                                                                                  # debugging_dir=None)
    print('unfolded alt chi2, mdf, chi2/ndf, p:', unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_chi2/unfolded_alt_ndf, unfolded_alt_p)

    return {
                "region": setup.region['label'],
                "is_groomed": "groomed" in setup.region['name'],
                "angle": setup.angle.var,

                "smeared_chi2": smeared_chi2,
                "smeared_ndf": smeared_ndf,
                "smeared_chi2_ov_ndf": smeared_chi2/smeared_ndf,
                "smeared_p": smeared_p,

                "smeared_alt_chi2": smeared_alt_chi2,
                "smeared_alt_ndf": smeared_alt_ndf,
                "smeared_alt_chi2_ov_ndf": smeared_alt_chi2/smeared_alt_ndf,
                "smeared_alt_p": smeared_alt_p,

                "unfolded_chi2": unfolded_chi2,
                "unfolded_ndf": unfolded_ndf,
                "unfolded_chi2_ov_ndf": unfolded_chi2/unfolded_ndf,
                "unfolded_p": unfolded_p,

                "unfolded_alt_chi2": unfolded_alt_chi2,
                "unfolded_alt_ndf": unfolded_alt_ndf,
                "unfolded_alt_chi2_ov_ndf": unfolded_alt_chi2/unfolded_alt_ndf,
                "unfolded_alt_p": unfolded_alt_p,
            }


def print_chi2_table(df):
    df_sorted = df.sort_values(by=['region', 'angle'])
    # print(r'\begin{tabular}{c|c|c|c|c|c|c|c|c|c|c|c}')
    print(r'\begin{tabular}{c|c|c|c|c|c|c|c}')
    print(("Region & "
           "Angle & "
           r"$\chi_{\text{unfolded}}^2 / N_{DoF}$ & "
           # r"$p_{\text{unfolded}}$ & "
           r"$\chi_{\text{unfolded, alt. model}}^2 / N_{DoF}$ & "
           # r"$p_{\text{unfolded, alt. model}}$ & "
           r"$\chi_{\text{smeared}}^2 / N_{DoF}$ & "
           # r"$p_{\text{smeared}}$ & "
           r"$\chi_{\text{smeared, alt. model}}^2 / N_{DoF}$ & "
           # r"$p_{\text{smeared, alt. model}}$ & "
           r"$\chi_{\text{unfolded}}^2 < \chi_{\text{smeared}}^2$ & "
           r"$\chi_{\text{unfolded, alt. model}}^2 < \chi_{\text{smeared, alt. model}}^2$ \\"))
    print(r"\hline \\")
    for row in df_sorted.itertuples():
        angle_name = row.angle.replace("jet_", "").replace("_charged", " (charged)")
        if row.is_groomed:
            angle_name = "Groomed " + angle_name
        # print(dir(row))
        print(("{region} & "
               "{angle_name} & "
               "{unfolded_chi2:.3e} / {unfolded_ndf} = {unfolded_chi2_ov_ndf:.3g} & "
               # "{unfolded_p:.3e} & "
               "{unfolded_alt_chi2:.3e} / {unfolded_alt_ndf} = {unfolded_alt_chi2_ov_ndf:.3g} & "
               # "{unfolded_alt_p:.3e} & "
               "{smeared_chi2:.3e} / {smeared_ndf} = {smeared_chi2_ov_ndf:.3g} & "
               # "{smeared_p:.3e} & "
               "{smeared_alt_chi2:.3e} / {smeared_alt_ndf} = {smeared_alt_chi2_ov_ndf:.3g} & "
               # "{smeared_alt_p:.3e} & "
               "{result} & "
               "{result_alt} \\\\").format(angle_name=angle_name,
                                           result=r"$\checkmark$" if row.unfolded_chi2 < row.smeared_chi2 else r"$\times$",
                                           result_alt=r"$\checkmark$" if row.unfolded_alt_chi2 < row.smeared_alt_chi2 else r"$\times$",
                                           **row._asdict()))
    print(r'\end{tabular}')


def slim_region(region):
    """Remove excess bits from Region dict to save memory"""
    if 'pdf_systematics' in region:
        del region['pdf_systematics']
    if 'scale_systematics' in region:
        del region['scale_systematics']
    if 'model_systematics' in region:
        del region['model_systematics']


class PlotWebpageMaker(object):

    def __init__(self, webpage_dir):
        self.webpage_dir = webpage_dir
        self.thumb_dir = os.path.join(webpage_dir, 'thumbnails')
        self.pdf_dir = os.path.join(webpage_dir, 'images')
        self.plot_setups = []  # store info each time we call make_webpage_plots()

    def make_webpage_plots(self, setup, do_plotting=True):
        """Create all thumbnails & images for webpage"""
        thumb_setup = Setup(jet_algo=setup.jet_algo,
                            region=setup.region,
                            angle=setup.angle,
                            output_dir=self.thumb_dir,
                            has_data=setup.has_data)
        if not do_plotting:
            slim_region(thumb_setup.region)
        thumb_setup.output_fmt = "gif"
        if do_plotting:
            do_binned_plots_per_region_angle(setup=thumb_setup,
                                             do_binned_gen_pt=True,
                                             do_binned_gen_lambda=False,
                                             do_binned_reco_pt=True,
                                             only_paper_plots=True)

        pdf_setup = Setup(jet_algo=setup.jet_algo,
                          region=setup.region,
                          angle=setup.angle,
                          output_dir=self.pdf_dir,
                          has_data=setup.has_data)
        if do_plotting:
            do_binned_plots_per_region_angle(setup=pdf_setup,
                                             do_binned_gen_pt=True,
                                             do_binned_gen_lambda=False,
                                             do_binned_reco_pt=True,
                                             only_paper_plots=True)

        self.plot_setups.append([thumb_setup, pdf_setup])

    def make_webpage(self, source_dir=""):
        """Make the HTML pages"""
        env = Environment(loader=FileSystemLoader('.'), trim_blocks=True, lstrip_blocks=True)
        template = env.get_template("unfolding_webpage_template.html")

        unfolded_file = os.path.join(self.webpage_dir, "index.html")
        uncertainty_file = os.path.join(self.webpage_dir, "index_uncertainty.html")
        detector_file = os.path.join(self.webpage_dir, "index_detector.html")

        menu_items = [
            {'href': os.path.basename(unfolded_file), 'text': "Unfolded differential plots"},
            {'href': os.path.basename(uncertainty_file), 'text': "Unfolded uncertainty plots"},
            {'href': os.path.basename(detector_file), 'text': "Detector differential plots"},
        ]

        def _generate_title(setup):
            """Geenrate title for a given setup (jet algo, region, angle)"""
            groom_str = " Groomed" if "groomed" in setup.region['name'] else ''
            return '%s %s%s %s' % (setup.jet_algo, setup.region['label'], groom_str, setup.angle.name)

        def _generate_id(setup):
            """Geenrate ID for a given setup (jet algo, region, angle)"""
            return '%s_%s_%s' % (setup.jet_algo, setup.region['name'], setup.angle.var)

        def _get_pdf_filename(thumb_setup, pdf_setup, thumb_filename):
            """Convert thumbnail filename to PDF filename"""
            pdf_base = os.path.basename(thumb_filename).replace("."+thumb_setup.output_fmt, "."+pdf_setup.output_fmt)
            return os.path.relpath(os.path.join(pdf_setup.output_dir, pdf_base),
                                   self.webpage_dir)

        # make differential unfolded page
        thumb_item_groups = []
        for thumb_setup, pdf_setup in self.plot_setups:
            thumbnails = [
                {
                    'src': thumb_fname,
                    'href': _get_pdf_filename(thumb_setup, pdf_setup, thumb_fname),
                    'alt': ''
                }
                for thumb_fname
                in glob(thumb_setup.output_dir + "/unfolded_%s_alt_truth_bin_*_divBinWidth.%s" % (pdf_setup.append, thumb_setup.output_fmt))
                if "lowX" not in thumb_fname and "logY" not in thumb_fname
            ]
            this_dict = {
                'title': _generate_title(thumb_setup),
                'id': _generate_id(thumb_setup),
                'thumbnails': thumbnails
            }
            thumb_item_groups.append(this_dict)
            print('Adding', this_dict['title'])

        template_dict = dict(
            title=source_dir + "</h2><h2>Unfolded plots",
            menu=menu_items,
            thumbnail_groups=thumb_item_groups
        )
        output = template.render(template_dict)
        with open(unfolded_file, 'w') as f:
            f.write(output)
        print("Webpage written to", unfolded_file)

        # make systematics page
        thumb_item_groups = []
        for thumb_setup, pdf_setup in self.plot_setups:
            thumbnails = [
                {
                    'src': thumb_fname,
                    'href': _get_pdf_filename(thumb_setup, pdf_setup, thumb_fname),
                    'alt': ''
                }
                for thumb_fname
                in glob(thumb_setup.output_dir + "/unfolded_syst_variations_vs_nominal_%s_bin_*_divBinWidth.%s" % (pdf_setup.append, thumb_setup.output_fmt))
            ]
            this_dict = {
                'title': _generate_title(thumb_setup),
                'id': _generate_id(thumb_setup),
                'thumbnails': thumbnails
            }
            thumb_item_groups.append(this_dict)

        template_dict = dict(
            title=source_dir + "</h2><h2>Unfolded uncertainty plots",
            menu=menu_items,
            thumbnail_groups=thumb_item_groups
        )
        output = template.render(template_dict)
        with open(uncertainty_file, 'w') as f:
            f.write(output)
        print("Webpage written to", uncertainty_file)

        # make differential detector page
        thumb_item_groups = []
        for thumb_setup, pdf_setup in self.plot_setups:
            thumbnails = [
                {
                    'src': thumb_fname,
                    'href': _get_pdf_filename(thumb_setup, pdf_setup, thumb_fname),
                    'alt': ''
                }
                for thumb_fname
                in glob(thumb_setup.output_dir + "/detector_reco_binning_%s_bin_*_divBinWidth.%s" % (pdf_setup.append, thumb_setup.output_fmt))
            ]
            this_dict = {
                'title': _generate_title(thumb_setup),
                'id': _generate_id(thumb_setup),
                'thumbnails': thumbnails
            }
            thumb_item_groups.append(this_dict)

        template_dict = dict(
            title=source_dir + "</h2><h2>Detector plots",
            menu=menu_items,
            thumbnail_groups=thumb_item_groups
        )
        output = template.render(template_dict)
        with open(detector_file, 'w') as f:
            f.write(output)
        print("Webpage written to", detector_file)


@profile
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory (should be the one made by unfolding.py)")
    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT['ungroomed'].keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them (default).",
                        default=["all"])
    parser.add_argument("--doBinnedPlots",
                        action='store_true',
                        help='Do all lambda/pt binned plots')
    parser.add_argument("--doBinnedPlotsGenPt",
                        action='store_true',
                        help='Do lambda plots, binned by gen pT')
    parser.add_argument("--doBinnedPlotsGenLambda",
                        action='store_true',
                        help='Do pT plots, binned by gen lambda')
    parser.add_argument("--doBinnedPlotsRecoPt",
                        action='store_true',
                        help='Do lambda plots, binned by reco pT')
    parser.add_argument("--doBigNormed1DPlots",
                        action='store_true',
                        help='Do big 1D plots (all normalised hists on one big axis)')
    parser.add_argument("--doBigAbs1DPlots",
                        action='store_true',
                        help='Do big absolute 1D plots')
    parser.add_argument("--doBottomLineTest",
                        action='store_true',
                        help='Do bottom line test')

    parser.add_argument("--onlyPaperPlots",
                        action='store_true',
                        help='Only do paper plots (applies to --doBinnedPlotsGenPt, etc)')

    parser.add_argument("--paper",
                        action='store_true',
                        help='Don\'t add "Preliminary" to plots')

    parser.add_argument("--webpage",
                        action='store_true',
                        help='Do webpage of unfolded results & systematics')
    parser.add_argument("--webpagePlots",
                        action='store_true',
                        help='Do plots for webpage')

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

    regions = setup_regions_from_argparse(args)
    if len(regions) == 0:
        raise RuntimeError("Need at least 1 region")

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    if len(angles) == 0:
        raise RuntimeError("Need at least 1 angle")

    if args.doBinnedPlots:
        for x in ['doBinnedPlotsRecoPt', 'doBinnedPlotsGenLambda', 'doBinnedPlotsGenPt']:
            setattr(args, x, True)

    # jet_algo = "AK4 PF PUPPI"
    jet_algo = "AK4"
    if "ak8" in args.source:
        # jet_algo = "AK8 PF PUPPI"
        jet_algo = "AK8"

    has_data = not ('_MC_all' in args.source or '_MC_split' in args.source)

    webpage_dir = os.path.join(args.source, 'webpages')
    webpage_maker = PlotWebpageMaker(webpage_dir=webpage_dir)

    all_chi2_stats = []

    num_all_iterations = len(regions) * len(angles)
    counter = 1

    # Iterate through regions & variables
    for region in regions:
        prof_start_region()

        region_dir = os.path.join(args.source, region['name'])
        if not os.path.isdir(region_dir):
            print("! Warning ! cannot find region dir", region_dir, '- skipping region')
            continue

        for angle in angles:
            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
            print("*"*120)
            print("Algo/Region/Var: %s %s %s (%d/%d)" % (jet_algo, region['name'], angle.var, counter, num_all_iterations))
            print("*"*120)
            counter += 1

            angle_output_dir = "%s/%s/%s" % (args.source, region['name'], angle.var)
            if not os.path.isdir(angle_output_dir):
                print("! Warning ! cannot find angle dir", angle_output_dir, '- skipping angle', angle.var)
                continue

            prof_unpickle_angle()

            # Get region dict from pickle file
            pickle_filename = os.path.join(angle_output_dir, "unfolding_result.pkl")
            unpickled_region = unpickle_region(pickle_filename)
            if unpickled_region is None:
                continue

            # check
            if region['name'] != unpickled_region['name']:
                raise RuntimeError("Mismatch region name")

            # use region dict from unpickling
            # don't use update(), mega slow
            this_region = unpickled_region
            this_region['label'] = region['label']

            prof_start_binned_plots()

            # MAKE ALL THE PLOTS
            # ------------------------------------------------------------------
            setup = Setup(jet_algo=jet_algo,
                          region=this_region,
                          angle=angle,
                          output_dir=angle_output_dir,
                          has_data=has_data,
                          is_preliminary=not args.paper)

            hist_bin_chopper = this_region['unfolder'].hist_bin_chopper

            if any([args.doBinnedPlotsGenPt, args.doBinnedPlotsGenLambda, args.doBinnedPlotsRecoPt, args.onlyPaperPlots]):
                print("...........................................................")
                print(" Doing lots of plots per pt/lambda bin...")
                print("...........................................................")
                hist_bin_chopper = do_binned_plots_per_region_angle(setup,
                                                                    do_binned_gen_pt=args.doBinnedPlotsGenPt or args.onlyPaperPlots,
                                                                    do_binned_gen_lambda=args.doBinnedPlotsGenLambda or args.onlyPaperPlots,
                                                                    do_binned_reco_pt=args.doBinnedPlotsRecoPt or args.onlyPaperPlots,
                                                                    only_paper_plots=args.onlyPaperPlots)
            # prof_done_binned_plots()

            # print("Setup obj:")
            # cu.print_dict_item_sizes(setup.__dict__)
            # print("region PDF obj:")
            # cu.print_dict_item_sizes(setup.region['pdf_systematics'][0])
            # print("region scale obj:")
            # cu.print_dict_item_sizes(setup.region['scale_systematics'][0])

            if args.doBigNormed1DPlots:
                print("...........................................................")
                print(" Doing big normed 1D plots...")
                print("...........................................................")
                # Do a 1D summary plot, with all the normalised plots with bins divided by their width
                # (unlike the standard plot from MyUnfolderPlotter, which is absolute)
                do_all_big_normalised_1d_plots_per_region_angle(setup, hist_bin_chopper)

            # prof_done_big_normed_plots()

            if args.doBigAbs1DPlots:
                print("...........................................................")
                print(" Doing big absolute 1D plots...")
                print("...........................................................")
                # Do standard 1D absolute plots
                do_all_big_absolute_1d_plots_per_region_angle(setup)

            # prof_done_big_abs_plots()

            if args.doBottomLineTest:
                all_chi2_stats.append(get_bottom_line_stats(setup))


            if args.webpage:
                webpage_maker.make_webpage_plots(setup, do_plotting=False)
            if args.webpagePlots:
                webpage_maker.make_webpage_plots(setup, do_plotting=True)

            prof_done_chi2()

            # cleanup object, as it uses loads of memory
            if num_all_iterations > 1:
                print("...tidying up...")
                del hist_bin_chopper
                del unpickled_region
                del setup

            prof_done_cleanup()

    if args.webpage:
        prof_start_webpage()
        webpage_maker.make_webpage(source_dir=args.source)
        prof_done_webpage()

    if len(all_chi2_stats) > 0:
        df_stats = pd.DataFrame(all_chi2_stats)
        df_stats['region'] = df_stats['region'].astype('category')
        df_stats['angle'] = df_stats['angle'].astype('category')
        print(df_stats.head())
        print_chi2_table(df_stats)


@profile
def prof_start_region():
    pass

@profile
def prof_unpickle_angle():
    pass

@profile
def prof_start_binned_plots():
    pass

@profile
def prof_done_binned_plots():
    pass

@profile
def prof_done_binned_plots():
    pass

@profile
def prof_done_big_normed_plots():
    pass

@profile
def prof_done_big_abs_plots():
    pass

@profile
def prof_done_chi2():
    pass

@profile
def prof_done_cleanup():
    pass

@profile
def prof_start_webpage():
    pass

@profile
def prof_done_webpage():
    pass


if __name__ == "__main__":
    main()
