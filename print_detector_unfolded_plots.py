#!/usr/bin/env python


"""
Plot detector & unfolding distributions together for data only

Does dijet/zpj separately, and together

"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
from copy import copy, deepcopy
from itertools import chain
import numpy as np
from pprint import pprint
from array import array
from math import floor, log10


import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, HistBinChopper, unpickle_region
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import get_dijet_config, get_zpj_config
# from do_unfolding_plots import PLOT_STYLES
import metric_calculators as metrics

# Doesn't work with python 3.8 for now
if not (sys.version_info.major == 3 and sys.version_info.minor >= 8):
    import rootpy
    # import rootpy.logger.magic as M; M.DANGER.enabled = True

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


# Define own linestyle for smaller plots
# linestyle 2 (dashed) has too big dashes
ROOT.gStyle.SetLineStyleString(22, "8 6")
ROOT.gStyle.SetLineStyleString(23, "5 10")


def scale_ematrix_by_bin_widths(ematrix, widths):
    this_widths = widths.reshape(len(widths), 1)
    return ematrix * this_widths * this_widths.T


class DijetZPJGenPtBinnedPlotter(object):
    """Class to do plots of dijet and/or ZPJ with generator pt binning"""

    def __init__(self, jet_algo, angle, bins, dijet_region, zpj_region, output_dir, is_groomed, is_preliminary):
        self.jet_algo = jet_algo
        self.has_data = has_data
        self.angle = angle
        self.dijet_region = dijet_region
        if self.dijet_region is not None:
            self.dijet_hbc = self.dijet_region['unfolder'].hist_bin_chopper
            self.dijet_hbc.add_obj("input_hist_gen_binning_bg_subtracted", dijet_region['unfolder'].input_hist_gen_binning_bg_subtracted)
        self.zpj_region = zpj_region
        if self.zpj_region is not None:
            self.zpj_hbc = self.zpj_region['unfolder'].hist_bin_chopper
            self.zpj_hbc.add_obj("input_hist_gen_binning_bg_subtracted", zpj_region['unfolder'].input_hist_gen_binning_bg_subtracted)
        self.bins = bins
        self.is_groomed = is_groomed

        # for plot axis titles
        self.pt_var_str = "p_{T}^{jet}"
        self.pt_str = self.pt_var_str + " [GeV]"

        angle_prepend = "Groomed " if self.is_groomed else ""
        lower_angle_name = qgc.lower_angle_name(angle)
        self.angle_str = "{prepend}{name} ({lambda_str})".format(prepend=angle_prepend,
                                                                 name=lower_angle_name,
                                                                 lambda_str=angle.lambda_str)
        # self.particle_title = "Particle-level " + self.angle_str
        self.particle_title = self.angle_str

        self.detector_title = "Detector-level {prepend}{name} ({lambda_str})".format(prepend=angle_prepend.lower(),  # lower as no longer first
                                                                                     name=lower_angle_name,
                                                                                     lambda_str=angle.lambda_str)
        self.pt_bin_normalised_differential_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_detector_normalised_differential_label = "#frac{1}{dN/dp_{T}} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_normalised_differential_times_width_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_unnormalised_differential_label = "#frac{1}{dN/dp_{T}} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str)  # FIXME
        self.pt_bin_unnormalised_differential_label = "N"
        self.lambda_bin_normalised_differential_label = "#frac{1}{d#sigma/d%s} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)
        self.lambda_bin_unnormalised_differential_label = "#frac{1}{dN/d%s} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)  # FIXME
        self.lambda_bin_unnormalised_differential_label = "N"
        self.output_dir = output_dir
        self.output_fmt = 'pdf'
        angle_prepend = "groomed_" if self.is_groomed else ""
        self.append = "%s%s" % (angle_prepend, angle.var)  # common str to put on filenames, etc
        self.is_preliminary = is_preliminary
        self.paper_str = "paper" if not self.is_preliminary else ""

        self.plot_colours = dict(
            dijet_colour=ROOT.kBlack,
            zpj_colour=ROOT.kRed
        )
        self.line_width = 2
        self.line_style_detector = 22


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
                  "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                 .format(
                    jet_algo=self.jet_algo,
                    pt_str=self.pt_var_str,
                    bin_edge_low=bin_edge_low,
                    bin_edge_high=bin_edge_high
                ))
        return title

    @staticmethod
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

    @staticmethod
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
        mean = float(metrics.calc_mean_jax(bin_areas, centers))
        mean_err = float(metrics.calc_mean_correlated_error_jax(bin_areas, centers, cov_matrix))
        return mean, mean_err

    def plot_detector_unfolded(self, do_dijet=True, do_zpj=True):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, axis='pt', do_norm=True, do_div_bin_width=True, binning_scheme='generator')
            print("----------pt bin", ibin, "=", bin_edge_low, "-", bin_edge_high)

            def determine_num_dp(values):
                """Determine correct number of decimal places
                so that the smallest value has 1 digit

                If smallest > 1, then just returns 0
                """
                smallest_value = min([v for v in values if v > 0])
                print('determine_num_dp', smallest_value)
                if smallest_value > 1:
                    return 0
                # use log10 to figure out exponent, then floor it
                # e.g. log10(5.5E-3) = -2.25...., floor(-2.25...) = -3
                n_dp = abs(floor(log10(smallest_value)))
                return n_dp

            def setup_beta_function(name):
                fit = ROOT.TF1(name, "[2]*TMath::BetaDist(x,[0],[1])", 0, 1)
                fit.SetParameter(0, 3)
                fit.SetParLimits(0, 0, 100)
                fit.SetParameter(1, 5)
                fit.SetParLimits(1, 0, 100)
                fit.SetParameter(2, 0.05)
                fit.SetParLimits(2, 0, 10)
                # fit.SetParameter(3, 0)
                # fit.SetParLimits(3, 0, 10)
                # fit.SetParameter(4, 1)
                # fit.SetParLimits(4, 0.1, 10)
                return fit

            fit_opts = "S"

            # Weird setup here: basically need all the values going into legend first,
            # to be able to determine the number of decimal points
            # Then we can actually create the Contributions and plot afterwards
            dijet_detector_hist, dijet_unfolded_hist = None, None
            dijet_detector_mean, dijet_detector_mean_err = None, None
            dijet_unfolded_mean, dijet_unfolded_mean_err = None, None
            zpj_detector_hist, zpj_unfolded_hist = None, None
            zpj_detector_mean, zpj_detector_mean_err = None, None
            zpj_unfolded_mean, zpj_unfolded_mean_err = None, None
            dijet_entries, zpj_entries = [], []
            errors = []

            if do_dijet:
                # get detector-level data
                dijet_detector_hist = self.dijet_hbc.get_bin_plot('input_hist_gen_binning_bg_subtracted', **hbc_args)
                dijet_detector_mean, dijet_detector_mean_err = self.get_uncorrelated_mean_err(dijet_detector_hist, is_density=True)
                errors.append(dijet_detector_mean_err)

                # fit with beta function
                # dijet_detector_fit = setup_beta_function("beta_fit_dijet_detector")
                # fit_result = dijet_detector_hist.Fit(dijet_detector_fit, fit_opts, "")
                # fit_result.Print()

                # get unfolded data
                dijet_unfolded_hist = self.dijet_hbc.get_bin_plot("unfolded", **hbc_args)
                dijet_cov_matrix = self.dijet_hbc.get_bin_plot(self.dijet_region['unfolder'].total_ematrix_name, **hbc_args)
                dijet_unfolded_mean, dijet_unfolded_mean_err = self.get_correlated_mean_err(dijet_unfolded_hist, dijet_cov_matrix, is_density=True)
                errors.append(dijet_unfolded_mean_err)

                # dijet_unfolded_fit = setup_beta_function("beta_fit_dijet_unfolded")
                # fit_result = dijet_unfolded_hist.Fit(dijet_unfolded_fit, fit_opts, "")
                # fit_result.Print()

            if do_zpj:
                # get detector-level data
                zpj_detector_hist = self.zpj_hbc.get_bin_plot('input_hist_gen_binning_bg_subtracted', **hbc_args)
                zpj_detector_mean, zpj_detector_mean_err = self.get_uncorrelated_mean_err(zpj_detector_hist, is_density=True)
                errors.append(zpj_detector_mean_err)

                # zpj_detector_fit = setup_beta_function("beta_fit_zpj_detector")
                # fit_result = zpj_detector_hist.Fit(zpj_detector_fit, fit_opts, "")
                # fit_result.Print()

                # get unfolded data
                zpj_unfolded_hist = self.zpj_hbc.get_bin_plot("unfolded", **hbc_args)
                zpj_cov_matrix = self.zpj_hbc.get_bin_plot(self.zpj_region['unfolder'].total_ematrix_name, **hbc_args)
                zpj_unfolded_mean, zpj_unfolded_mean_err = self.get_correlated_mean_err(zpj_unfolded_hist, zpj_cov_matrix, is_density=True)
                errors.append(zpj_unfolded_mean_err)

                # zpj_unfolded_fit = setup_beta_function("beta_fit_zpj_unfolded")
                # fit_result = zpj_unfolded_hist.Fit(zpj_unfolded_fit, fit_opts, "")
                # fit_result.Print()

            # n_dp = determine_num_dp(errors)
            n_dp = 3

            # kerning necessary as it puts massive space around #pm
            # but the first #kern doesn't affect the space after #pm as much (?!),
            # so I have to add another one with harder kerning
            # we use %.(n_dp)f as our float format str to ensure the correct number of dp are shown (and not rounded off)
            stat_template = 'Mean = {:.%df}#kern[-0.2dx]{{ #pm}}#kern[-0.5dx]{{ }}{}' % (n_dp)
            err_template = "{:.%df}" % n_dp

            def _stat_label(mean, err, dp):
                err_str = err_template.format(err)
                # if the error is so small that it would display as 0.000,
                # i.e. < 0.0005, instead show < 0.001
                if err < 5*10**(-dp-1):
                    err_str = "#lower[-0.09dy]{<}#kern[-0.75dx]{ }" + err_template.format(1*10**(-dp))
                return stat_template.format(round(mean, dp), err_str)

            if do_dijet:
                dijet_entries.append(Contribution(dijet_detector_hist,
                                                  label='Detector-level\n%s' % (_stat_label(dijet_detector_mean, dijet_detector_mean_err, n_dp)),
                                                  line_color=self.plot_colours['dijet_colour'],
                                                  line_width=self.line_width,
                                                  line_style=self.line_style_detector,
                                                  marker_color=self.plot_colours['dijet_colour'],
                                                  marker_style=cu.Marker.get('circle', filled=False),
                                                  marker_size=1.2,
                                                  subplot=zpj_detector_hist if do_zpj else None))
                dijet_entries.append(Contribution(dijet_unfolded_hist,
                                                  label='Particle-level\n%s' % (_stat_label(dijet_unfolded_mean, dijet_unfolded_mean_err, n_dp)),
                                                  line_color=self.plot_colours['dijet_colour'],
                                                  line_width=self.line_width,
                                                  line_style=1,
                                                  marker_color=self.plot_colours['dijet_colour'],
                                                  marker_style=cu.Marker.get('circle', filled=True),
                                                  marker_size=1.2,
                                                  subplot=zpj_unfolded_hist if do_zpj else None))

            if do_zpj:
                zpj_entries.append(Contribution(zpj_detector_hist,
                                            label='Detector-level\n%s' % (_stat_label(zpj_detector_mean, zpj_detector_mean_err, n_dp)),
                                            line_color=self.plot_colours['zpj_colour'],
                                            line_width=self.line_width,
                                            line_style=self.line_style_detector,
                                            marker_color=self.plot_colours['zpj_colour'],
                                            marker_style=cu.Marker.get('square', filled=False),
                                            marker_size=1.2))
                zpj_entries.append(Contribution(zpj_unfolded_hist,
                                            label='Particle-level\n%s' % (_stat_label(zpj_unfolded_mean, zpj_unfolded_mean_err, n_dp)),
                                            line_color=self.plot_colours['zpj_colour'],
                                            line_width=self.line_width,
                                            line_style=1,
                                            marker_color=self.plot_colours['zpj_colour'],
                                            marker_style=cu.Marker.get('square', filled=True),
                                            marker_size=1.2))

            all_entries = list(chain(dijet_entries, zpj_entries))
            plot = Plot(all_entries,
                        ytitle=self.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        legend=True,
                        xlim=qgp.calc_auto_xlim(all_entries),  # set x lim to where data is non-0
                        what="hist",
                        xtitle=self.particle_title,
                        has_data=self.has_data,
                        ylim=[0, None],
                        is_preliminary=self.is_preliminary
                        )
            # plot.default_canvas_size = (600, 600)
            # plot.left_margin = 0.2
            # plot.left_title_offset_fudge_factor = 8
            plot.y_padding_max_linear = 1.8
            plot.top_margin = 0.07
            plot.title_start_y = 0.888
            plot.cms_text_y = 0.94
            plot.lumi = cu.get_lumi_str(do_dijet=do_dijet, do_zpj=do_zpj)
            if do_zpj and do_dijet:
                plot.subplot_type = 'ratio'
                plot.subplot_title = 'Dijet / Z+jet'
                plot.subplot_limits = (0.25, 2.75)

            # disable adding objects to legend & drawing - we'll do it manually
            # since we want proper error bar
            plot.do_legend = False
            plot.splitline_legend = False
            # plot.legend.SetFillColor(ROOT.kRed)
            # plot.legend.SetFillStyle(1001)
            plot.plot("NOSTACK E1 X0")
            # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.5)
            plot.main_pad.cd()

            for e in [plot.contributions[0].obj,plot.contributions[2].obj]:
              e.Draw("HIST SAME")
            for e in [plot.contributions[1].obj,plot.contributions[3].obj]:
              e.Draw("L SAME")

            plot.canvas.cd()

            # unfolded_fit = ROOT.TF1("beta_fit_dijet_unfolded", "[2]*TMath::BetaDist(x,[0],[1])", 0, 1)
            # unfolded_fit.SetParameter(0, 3)
            # unfolded_fit.SetParLimits(0, 0, 100)
            # unfolded_fit.SetParameter(1, 5)
            # unfolded_fit.SetParLimits(1, 0, 100)
            # unfolded_fit.SetParameter(2, .1)
            # unfolded_fit.SetParLimits(2, 0, 1000)
            # # fit_result = unfolded_hist.Fit(unfolded_fit, "EMSR", "", 0, 1)
            # # fit_result.Print()
            # unfolded_fit.SetLineColor(ROOT.kRed)
            # unfolded_fit.SetLineWidth(2)
            # plot.main_pad.cd()
            # unfolded_fit.Draw("SAME")

            # Create dummy graphs with the same styling to put into the legend
            dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
            dummies = []  # to stop garbage collection
            dummy_entries = []  # to stop garbage collection
            label_height = 0.03
            legend_height = 0.12
            legend_x1 = 0.54
            legend_x2 = 0.75 # this doesn't really control width - legend_text_size mainly does
            label_left_offset = 0.01
            final_reading_factor=1.2
            label_text_size = 0.032*final_reading_factor
            label_top = plot.title_start_y
            legend_text_size = 0.028*final_reading_factor
            inter_region_offset = 0.025
            if do_dijet:
                dijet_legend = plot.legend.Clone()
                dijet_legend.SetX1(legend_x1)
                dijet_legend.SetX2(legend_x2)
                dijet_legend.SetY1(label_top-label_height-legend_height)
                dijet_legend.SetY2(label_top-label_height+0.01)
                dijet_legend.SetTextSize(legend_text_size)
                # Add text with region label
                dijet_pt = ROOT.TPaveText(legend_x1-label_left_offset, label_top-label_height, legend_x2-label_left_offset, label_top, "NDC NB")
                dijet_pt.SetFillStyle(0)
                dijet_pt.SetBorderSize(0)
                text = dijet_pt.AddText(qgc.Dijet_CEN_LABEL)
                text.SetTextAlign(11)
                text.SetTextFont(62)
                text.SetTextSize(label_text_size)
                dummies.append(dijet_pt)
                dummies.append(text)
                dijet_pt.Draw()
                for cont in dijet_entries:
                    this_dummy_entry = dummy_gr.Clone()
                    cont.update_obj_styling(this_dummy_entry)
                    if "\n" in cont.label:
                        parts = cont.label.split("\n")
                        dijet_legend.AddEntry(this_dummy_entry, "#splitline{%s}{%s}" % (parts[0], parts[1]), "LEP")
                    else:
                        dijet_legend.AddEntry(this_dummy_entry, cont.label, "LEP")
                    dummy_entries.append(this_dummy_entry)  # to avoid garbage collection
                dijet_legend.Draw()
                dummies.append(dijet_legend)
                # setup for Z+J
                label_top -= label_height+legend_height+inter_region_offset

            if do_zpj:
                zpj_legend = plot.legend.Clone()
                zpj_legend.SetX1(legend_x1)
                zpj_legend.SetX2(legend_x2)
                zpj_legend.SetY1(label_top-label_height-legend_height)
                zpj_legend.SetY2(label_top-label_height+0.01)
                zpj_legend.SetTextSize(legend_text_size)
                # Add text with region label
                zpj_pt = ROOT.TPaveText(legend_x1-label_left_offset, label_top-label_height, legend_x2-label_left_offset, label_top, "NDC NB")
                zpj_pt.SetFillStyle(0)
                zpj_pt.SetBorderSize(0)
                text = zpj_pt.AddText(qgc.ZpJ_LABEL)
                text.SetTextAlign(11)
                text.SetTextFont(62)
                text.SetTextSize(label_text_size)
                dummies.append(zpj_pt)
                dummies.append(text)
                zpj_pt.Draw()
                for cont in zpj_entries:
                    this_dummy_entry = dummy_gr.Clone()
                    cont.update_obj_styling(this_dummy_entry)
                    if "\n" in cont.label:
                        parts = cont.label.split("\n")
                        zpj_legend.AddEntry(this_dummy_entry, "#splitline{%s}{%s}" % (parts[0], parts[1]), "LEP")
                    else:
                        zpj_legend.AddEntry(this_dummy_entry, cont.label, "LEP")
                    dummy_entries.append(this_dummy_entry)
                zpj_legend.Draw()
                dummies.append(zpj_legend)

            # Add legend to ratio plot
            plot.subplot_pad.cd()

            for e in [plot.subplot_contributions[0]]:
              e.Draw("HIST SAME")
            for e in [plot.subplot_contributions[1]]:
              e.Draw("L SAME")


            plot.subplot_leg = ROOT.TLegend(0.3, 0.73, 0.9, 0.9)
            plot.subplot_leg.SetTextSize(0.07)
            plot.subplot_leg.SetFillStyle(0)
            plot.subplot_leg.SetNColumns(2)
            plot.subplot_leg.AddEntry(dummy_entries[0], "Detector-level", "LEP")
            plot.subplot_leg.AddEntry(dummy_entries[1], "Particle-level", "LEP")
            plot.subplot_leg.Draw()

            plot.canvas.cd()



            parts = ['detector_unfolded',
                     'dijet' if do_dijet else None,
                     'zpj' if do_zpj else None,
                     self.append,
                     'bin_%d' % ibin,
                     'divBinWidth',
                     f'{self.paper_str}.{self.output_fmt}']
            filename = '_'.join([x for x in parts if x])
            plot.save("%s/%s" % (self.output_dir, filename))
            print("%s/%s" % (self.output_dir, filename))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory (should be the one made by unfolding.py)")
    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT['ungroomed'].keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them (default).",
                        default=["all"])

    parser.add_argument("--paper",
                        action='store_true',
                        help='Don\'t add "Preliminary" to plots')
    args = parser.parse_args()

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    if len(angles) == 0:
        raise RuntimeError("Need at least 1 angle")

    jet_algo = "AK4"
    if "ak8puppi" in args.source:
        jet_algo = "AK8"

    has_data = not ('_MC_all' in args.source or '_MC_split' in args.source)

    grooming_status = [False, True][:1]
    num_all_iterations = len(grooming_status) * len(angles) * 2 # 2 for dijet & Z+Jet in each angle/grooming variation
    counter = 1

    # Iterate through (un)groomed, variables
    for do_grooming in grooming_status:

        for angle in angles:

            unpickled_regions = []

            # get all regions for this plot
            dijet_region = get_dijet_config(args.source, central=True, groomed=do_grooming)
            zpj_region = get_zpj_config(args.source, groomed=do_grooming)

            def _get_unpicked_region_dict(region):
                region_dir = os.path.join(args.source, region['name'])
                if not os.path.isdir(region_dir):
                    raise RuntimeError("! Warning ! cannot find region dir", region_dir, '- skipping region')

                global counter
                append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
                print("*"*120)
                print("Algo/Region/Var: %s %s %s (%d/%d)" % (jet_algo, region['name'], angle.var, counter, num_all_iterations))
                print("*"*120)
                counter += 1

                angle_output_dir = "%s/%s/%s" % (args.source, region['name'], angle.var)
                if not os.path.isdir(angle_output_dir):
                    raise RuntimeError("! Warning ! cannot find angle dir", angle_output_dir, '- skipping angle', angle.var)

                # Get region dict from pickle file
                pickle_filename = os.path.join(angle_output_dir, "unfolding_result.pkl")
                unpickled_region = unpickle_region(pickle_filename)
                if unpickled_region is None:
                    raise RuntimeError("Unpickling for %s gave None" % append)

                # check
                if region['name'] != unpickled_region['name']:
                    raise RuntimeError("Mismatch region name")

                return unpickled_region

            dijet_region = _get_unpicked_region_dict(dijet_region)
            zpj_region = _get_unpicked_region_dict(zpj_region)

            # Now make plots
            groom_str = "_groomed" if do_grooming else ""
            output_dir = "%s/detector_unfolded_dijet_zpj%s/%s" % (args.source, groom_str, angle.var)
            binned_plotter = DijetZPJGenPtBinnedPlotter(jet_algo=jet_algo,
                                                        angle=angle,
                                                        bins=dijet_region['unfolder'].pt_bin_edges_gen[:-3], # stop early as Z+J missing some
                                                        dijet_region=dijet_region,
                                                        zpj_region=zpj_region,
                                                        output_dir=output_dir,
                                                        is_groomed=do_grooming,
                                                        is_preliminary=not args.paper)

            # dijet-only
            # binned_plotter.plot_detector_unfolded(do_dijet=True, do_zpj=False)

            # ZPJ-only
            # binned_plotter.plot_detector_unfolded(do_dijet=False, do_zpj=True)

            # dijet+ZPJ
            binned_plotter.plot_detector_unfolded(do_dijet=True, do_zpj=True)


