#!/usr/bin/env python


"""
Do all the unfolding plots: per pT bin, per lambda bin, summary plot
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, unfolder_from_tdir
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import get_dijet_config, get_zpj_config


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")
ROOT.gStyle.SetHistTopMargin(0.)


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


def unpack_unfolding_root_file(input_tfile, region, angle):
    input_tdir_name = "%s/%s" % (region['name'], angle.var)
    input_tdir = input_tfile.Get(input_tdir_name)
    cu.check_root_obj(input_tdir)
    unfolder = unfolder_from_tdir(input_tdir)
    print("Loaded main unfolder")

    list_of_obj = cu.get_list_of_element_names(input_tdir)

    # Get unregularised unfolder, if available
    unreg_tdir = [x for x in list_of_obj if x.startswith("unreg_unfolder")]
    unreg_unfolder = None
    if len(unreg_tdir) == 1:
        unreg_unfolder = unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, unreg_tdir[0])))
        print("Loaded comparison unregularised unfolder")

    # Get alternate response object, if it exists
    alt_tdir = [x for x in list_of_obj if x.startswith("alt_response_")]
    alt_unfolder = None
    alt_unfolder_name = None
    alt_hist_truth = None
    if len(alt_tdir)  == 1:
        alt_unfolder = unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, alt_tdir[0])))
        region['alt_unfolder'] = alt_unfolder
        alt_unfolder_name = alt_tdir[0].replace("alt_response_", "").replace("_", " ")
        if region['alt_mc_label'] != alt_unfolder_name:
            raise RuntimeError("Bad unpacking of alt response unfolder: expected %s, got %s" % (region['alt_mc_label'], alt_unfolder_name))
        print("Loaded alt unfolder")
        alt_hist_truth = input_tfile.Get(os.path.join(input_tdir_name, alt_tdir[0], "alt_hist_mc_gen"))
    if len(alt_tdir) > 1:
        raise RuntimeError(">1 alt_response?! %s" % (alt_tdir))

    # Get model systs
    print(list_of_obj)
    model_tdirs = [x for x in list_of_obj if x.startswith("modelSyst_")]
    if len(model_tdirs) > 0:
        for model_tdir_name in model_tdirs:
            syst_name = model_tdir_name.replace("modelSyst_", "").replace("_", " ")
            this_one = [x for x in region['model_systematics'] if x['label'] == syst_name]
            if len(this_one) == 0:
                print("No entry for model systematic", syst_name, "- skipping")
                continue
            # TODO: check it agrees with region dict?
            this_one[0]['unfolder'] = unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, model_tdir_name)))
            print("Loaded", len(model_tdirs), "model systematic unfolders")
    # remove entries without an unfolder
    region['model_systematics'] = [k for k in region['model_systematics']
                                 if k.get('unfolder', None) is not None]

    # Get PDF systs
    # For some reason, this is done as a list instead of dict
    pdf_tdirs = [x for x in list_of_obj if x.startswith("pdfSyst_")]
    if len(pdf_tdirs) > 0:
        # Remove original, construct all other
        region['pdf_systematics'] = []

        for pdf_tdir_name in pdf_tdirs:
            pdf_name = pdf_tdir_name.replace("pdfSyst_", "")
            region['pdf_systematics'].append({
                'label': pdf_name,
                'unfolder': unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, pdf_tdir_name))),
                'colour': ROOT.kCyan+2,
            })
        print("Loaded", len(pdf_tdirs), "PDF systematic unfolders")
    # remove entries without an unfolder
    region['pdf_systematics'] = [k for k in region['pdf_systematics']
                                 if k.get('unfolder', None) is not None]

    return dict(
        unfolder=unfolder,
        alt_unfolder=alt_unfolder,
        alt_hist_truth=alt_hist_truth
    )


class HistBinChopper(object):
    """Get histogram for pt or variable bin, and cache it in dict, so can be used later"""

    def __init__(self, unfolder):
        self.unfolder = unfolder
        self.objects = {}
        self._cache = {}

    def add_obj(self, name, obj):
        # TODO: allow overwrite?
        self.objects[name] = obj

    def get_pt_bin(self, name, ind, binning_scheme='generator'):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_pt_bin_%d_%s" % (ind, binning_scheme)
        if this_name not in self._cache:
            self._cache[this_name] = self.unfolder.get_var_hist_pt_binned(self.objects[name], ind, binning_scheme)
        return self._cache[this_name]

    def get_pt_bin_div_bin_width(self, name, ind, binning_scheme='generator'):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_pt_bin_%d_%s_divBinWidth" % (ind, binning_scheme)
        if this_name not in self._cache:
            self._cache[this_name] = qgp.hist_divide_bin_width(self.get_pt_bin(name, ind, binning_scheme))
        return self._cache[this_name]

    def get_pt_bin_normed_div_bin_width(self, name, ind, binning_scheme='generator'):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_pt_bin_%d_%s_norm_divBinWidth" % (ind, binning_scheme)
        if this_name not in self._cache:
            self._cache[this_name] = qgp.normalise_hist_divide_bin_width(self.get_pt_bin(name, ind, binning_scheme))
        return self._cache[this_name]

    def get_lambda_bin(self, name, ind, binning_scheme='generator'):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_lambda_bin_%d_%s" % (ind, binning_scheme)
        if this_name not in self._cache:
            self._cache[this_name] = self.unfolder.get_pt_hist_var_binned(self.objects[name], ind, binning_scheme=binning_scheme)
        return self._cache[this_name]

    def get_lambda_bin_div_bin_width(self, name, ind, binning_scheme='generator'):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_lambda_bin_%d_%s_divBinWidth" % (ind, binning_scheme)
        if this_name not in self._cache:
            self._cache[this_name] = qgp.hist_divide_bin_width(self.get_lambda_bin(name, ind, binning_scheme))
        return self._cache[this_name]

    def get_lambda_bin_normed_div_bin_width(self, name, ind, binning_scheme='generator'):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_lambda_bin_%d_%s_norm_divBinWidth" % (ind, binning_scheme)
        if this_name not in self._cache:
            self._cache[this_name] = qgp.normalise_hist_divide_bin_width(self.get_lambda_bin(name, ind, binning_scheme))
        return self._cache[this_name]


class Setup(object):
    """Loads of common consts, useful for plotting etc"""

    def __init__(self, jet_algo, region, angle, output_dir='.', has_data=False, is_ave_pt_binning=False):
        self.jet_algo = jet_algo
        self.region = region
        self.pt_str = "#LT p_{T}^{jet} #GT" if is_ave_pt_binning else "p_{T}^{jet}"
        self.has_data = has_data
        self.angle = angle
        angle_prepend = "groomed " if "groomed" in region['name'] else ""
        this_angle_name = angle.name
        if (angle_prepend != ""
            and this_angle_name != 'LHA'
            and "_{T}" not in this_angle_name
            and "PUPPI" not in this_angle_name):
            # lower case if Groomed..., but be careful of e.g. pTD, LHA
            this_angle_name = this_angle_name[0].lower() + this_angle_name[1:]
        # for plot axis titles
        self.angle_str = "{prepend}{name} ({lambda_str})".format(prepend=angle_prepend,
                                                                 name=this_angle_name,
                                                                 lambda_str=angle.lambda_str)
        self.particle_title = "Particle-level " + self.angle_str
        self.detector_title = "Detector-level " + self.angle_str
        self.pt_bin_normalised_differential_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str)
        self.lambda_bin_normalised_differential_label = "#frac{1}{d#sigma/d%s}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)
        self.output_dir = output_dir
        self.output_fmt = 'pdf'
        self.append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name

PLOT_COLOURS = dict(
    gen_colour=ROOT.kRed,
    unfolded_basic_colour=ROOT.kAzure+7,
    unfolded_stat_colour=ROOT.kAzure+7,
    unfolded_total_colour=ROOT.kBlack,
    unfolded_unreg_colour=ROOT.kViolet+2,
    alt_gen_colour=ROOT.kViolet+1,
    alt_unfolded_colour=ROOT.kBlue-4,
    # reco_mc_colour=ROOT.kGreen+2,
    # reco_mc_colour=ROOT.kAzure-7,
    reco_data_colour=ROOT.kRed,
    reco_mc_colour=ROOT.kRed+3,
    reco_unfolding_input_colour=ROOT.kRed,
    reco_folded_unfolded_colour=ROOT.kAzure+1,
    reco_folded_mc_truth_colour=ROOT.kGreen+2,
)

# FIXME: generalise thise and LambdaBinnedPlotter into one generic BinnedPlotter?
# Although each has different set of plots, so not easy/possible
class GenPtBinnedPlotter(object):
    def __init__(self, setup, bins, hist_bin_chopper):
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
            subplot_limits=(0.75, 1.25),
        )

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
        return True

    def get_pt_bin_title(self, bin_edge_low, bin_edge_high):
        title = (("{jet_algo}\n"
                  "{region_label} region\n"
                  "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                 .format(
                    jet_algo=self.setup.jet_algo,
                    region_label=self.region['label'],
                    pt_str=self.setup.pt_str,
                    bin_edge_low=bin_edge_low,
                    bin_edge_high=bin_edge_high
                ))
        return title

    def plot_unfolded_unnormalised(self, unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', ibin)

            # unnormalised version
            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
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

    def plot_unfolded_normalised(self, unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
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

    def plot_unfolded_with_alt_truth_normalised(self, unfolder, alt_truth):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin)

            self.hist_bin_chopper.add_obj('alt_hist_truth', alt_truth)
            alt_mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', ibin)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(alt_mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['alt_mc_label']),
                             line_color=self.plot_colours['alt_gen_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_colours['alt_gen_colour'], marker_size=0,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
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

    def plot_unfolded_with_unreg_normalised(self, unfolder, unreg_unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin)

            self.hist_bin_chopper.add_obj("unreg_unfolded", unreg_unfolder.unfolded)
            unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unreg_unfolded', ibin)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
                Contribution(unreg_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = 0) (total err)",
                             line_color=self.plot_colours['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_unreg_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
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

    def plot_unfolded_with_alt_response_normalised(self, unfolder, alt_unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # TODO: should this be inside or outside this func?
            self.hist_bin_chopper.add_obj("alt_unfolded", alt_unfolder.unfolded)
            self.hist_bin_chopper.add_obj("alt_hist_truth", alt_unfolder.hist_truth)

            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin)
            alt_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_unfolded', ibin)
            alt_mc_gen_hist_gin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', ibin)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0,
                             normalise_hist=False),
                # Contribution(alt_mc_gen_hist_bin,
                #              label="Generator (%s)" % (self.region['alt_mc_label']),
                #              line_color=alt_gen_colour, line_width=self.line_width, line_style=2,
                #              marker_color=alt_gen_colour, marker_size=0,
                #              normalise_hist=False),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)\n(%s response matrix)" % (unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
                Contribution(alt_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)\n(%s response matrix)" % (alt_unfolder.tau, self.region['alt_mc_label']),
                             line_color=self.plot_colours['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['alt_unfolded_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
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

    def plot_unfolded_with_alt_response_truth_normalised(self, unfolder, alt_unfolder, alt_truth):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # TODO: should this be inside or outside this func?
            self.hist_bin_chopper.add_obj("alt_unfolded", alt_unfolder.unfolded)
            self.hist_bin_chopper.add_obj("alt_truth", alt_truth)

            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin)
            alt_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_unfolded', ibin)
            alt_mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_truth', ibin)

            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(alt_mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['alt_mc_label']),
                             line_color=self.plot_colours['alt_gen_colour'], line_width=self.line_width, line_style=2,
                             marker_color=self.plot_colours['alt_gen_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)\n(%s response matrix)" % (unfolder.tau, self.region['mc_label']),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
                Contribution(alt_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)\n(%s response matrix)" % (alt_unfolder.tau, self.region['alt_mc_label']),
                             line_color=self.plot_colours['alt_unfolded_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['alt_unfolded_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
            ]
            if not self.check_entries(entries, "plot_unfolded_with_unreg_normalised_pt_bin %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_%s_alt_response_truth_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_model_systs_normalised(self, unfolder, model_systs):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            syst_entries = []
            for syst_dict in model_systs:
                syst_unfolder = syst_dict['unfolder']
                syst_label = syst_dict['label']
                syst_label_no_spaces = syst_dict['label'].replace(" ", "_")

                self.hist_bin_chopper.add_obj('model_syst_%s_unfolded' % (syst_label_no_spaces), syst_unfolder.unfolded)
                self.hist_bin_chopper.add_obj('model_syst_%s_hist_truth' % (syst_label_no_spaces), syst_unfolder.hist_truth)

                syst_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('model_syst_%s_unfolded' % (syst_label_no_spaces), ibin)
                syst_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('model_syst_%s_hist_truth' % (syst_label_no_spaces), ibin)

                syst_entries.extend([
                    Contribution(syst_unfolded_hist_bin,
                                 label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_unfolder.tau, syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=1,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 subplot=syst_gen_hist_bin,
                                 normalise_hist=True),
                    Contribution(syst_gen_hist_bin,
                                 label="Generator (%s)" % (syst_label),
                                 line_color=syst_dict['colour'], line_width=self.line_width, line_style=2,
                                 marker_color=syst_dict['colour'], marker_size=0,
                                 normalise_hist=True),
                ])

            # add nominal ones last
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin)

            syst_entries.extend([
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
            ])
            if not self.check_entries(syst_entries, "plot_unfolded_with_model_systs_normalised_pt_bin %d" % (ibin)):
                return
            plot = Plot(syst_entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        **self.pt_bin_plot_args)
            plot.legend.SetX1(0.55)
            plot.legend.SetY1(0.72)
            plot.legend.SetX2(0.98)
            plot.legend.SetY2(0.88)
            if len(syst_entries) > 5:
                plot.legend.SetNColumns(2)
            plot.plot("NOSTACK E1")
            plot.save("%s/unfolded_%s_syst_model_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_pdf_systs_normalised(self, unfolder, pdf_systs):
        pdf_entries = []
        for pdf_dict in pdf_systs:
            pdf_unfolder = pdf_dict['unfolder']
            pdf_label = pdf_dict['label']
            pdf_label_no_spaces = pdf_dict['label'].replace(" ", "_")

            self.hist_bin_chopper.add_obj('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), pdf_unfolder.unfolded)
            self.hist_bin_chopper.add_obj('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), pdf_unfolder.hist_truth)

            pdf_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), ibin)
            pdf_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), ibin)

            pdf_entries.extend([
                Contribution(pdf_unfolded_hist_bin,
                             label="Unfolded (#tau = %.3g) (total err) (%s)" % (pdf_unfolder.tau, pdf_label),
                             line_color=pdf_dict['colour'], line_width=self.line_width, line_style=1,
                             marker_color=pdf_dict['colour'], marker_size=0,
                             subplot=pdf_gen_hist_bin,
                             normalise_hist=True),
                Contribution(pdf_gen_hist_bin,
                             label="Generator (%s)" % (pdf_label),
                             line_color=pdf_dict['colour'], line_width=self.line_width, line_style=2,
                             marker_color=pdf_dict['colour'], marker_size=0,
                             normalise_hist=True),
            ])

        # add nominal ones last
        mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
        unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin)

        pdf_entries.extend([
            Contribution(mc_gen_hist_bin,
                         label="Generator (%s)" % (self.region['mc_label']),
                         line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                         marker_color=self.plot_colours['gen_colour'], marker_size=0,
                         normalise_hist=False),
            Contribution(unfolded_hist_bin_total_errors,
                         label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                         line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                         marker_color=self.plot_colours['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                         subplot=mc_gen_hist_bin,
                         normalise_hist=False),
        ])
        if not self.check_entries(pdf_entries, "plot_unfolded_with_model_systs_normalised_pt_bin %d" % (ibin)):
            return
        plot = Plot(pdf_entries,
                    xtitle=self.setup.particle_title,
                    ytitle=self.setup.pt_bin_normalised_differential_label,
                    title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                    has_data=self.setup.has_data,
                    **self.pt_bin_plot_args)
        plot.legend.SetX1(0.55)
        plot.legend.SetY1(0.72)
        plot.legend.SetX2(0.98)
        plot.legend.SetY2(0.88)
        if len(pdf_entries) > 5:
            plot.legend.SetNColumns(2)
        plot.plot("NOSTACK E1")
        plot.save("%s/unfolded_%s_pdf_model_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_detector_normalised(self, unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            self.hist_bin_chopper.add_obj("input_hist_gen_binning_bg_subtracted", unfolder.input_hist_gen_binning_bg_subtracted)
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_gen_binning_bg_subtracted', ibin)
            self.hist_bin_chopper.add_obj("hist_mc_reco_gen_binning_bg_subtracted", unfolder.hist_mc_reco_gen_binning_bg_subtracted)
            mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_gen_binning_bg_subtracted', ibin)

            entries = [
                Contribution(mc_hist_bin,
                             label="MC (bg-subtracted)",
                             line_color=self.plot_colours['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_mc_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(input_hist_bin,
                             label="Data (bg-subtracted)",
                             line_color=self.plot_colours['reco_data_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_data_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_hist_bin,
                             normalise_hist=False),
            ]
            if not self.check_entries(entries, "plot_detector_normalised %d" % (ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='Data / MC',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/detector_gen_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))


# ================================================================================


class GenLambdaBinnedPlotter(object):
    def __init__(self, setup, bins, hist_bin_chopper):
        self.setup = setup
        self.region = setup.region
        self.bins = bins
        self.hist_bin_chopper = hist_bin_chopper

        self.line_width = 2
        self.plot_colours = PLOT_COLOURS
        self.lambda_bin_plot_args = dict(
            what="hist",
            xtitle="%s [GeV]" % self.setup.pt_str,
            has_data=self.setup.has_data,
            subplot_type='ratio',
            subplot_title="Unfolded / Gen",
            subplot_limits=(0.75, 1.25),
        )

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

    def plot_unfolded_unnormalised(self, unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin)

            # unnormalised version
            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_stat_errors,
                             label="Unfolded (#tau = %.3g) (stat err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_stat_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
            ]
            if not self.check_entries(entries, "plot_unfolded_unnormalised %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle="N",
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.y_padding_max_log = 5000  # space for title
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            plot.save("%s/unfolded_unnormalised_%s_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_unfolded_with_unreg_unnormalised(self, unfolder, unreg_unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            mc_gen_hist_bin = self.hist_bin_chopper.get_lambda_bin_div_bin_width('hist_truth', ibin)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', ibin)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', ibin)

            self.hist_bin_chopper.add_obj("unreg_unfolded", unreg_unfolder.unfolded)
            unreg_unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unreg_unfolded', ibin)

            # unnormalised version
            entries = [
                Contribution(mc_gen_hist_bin,
                             label="Generator (%s)" % (self.region['mc_label']),
                             line_color=self.plot_colours['gen_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['gen_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                             line_color=self.plot_colours['unfolded_total_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
                Contribution(unreg_unfolded_hist_bin_total_errors,
                             label="Unfolded (#tau = 0) (total err)",
                             line_color=self.plot_colours['unfolded_unreg_colour'], line_width=self.line_width, line_style=1,
                             marker_color=self.plot_colours['unfolded_unreg_colour'], #marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=False),
            ]
            if not self.check_entries(entries, "plot_unfolded_with_unreg_unnormalised %d" % (ibin)):
                return
            plot = Plot(entries,
                        ytitle="N",
                        title=self.get_lambda_bin_title(bin_edge_low, bin_edge_high),
                        **self.lambda_bin_plot_args)
            self._modify_plot(plot)
            plot.y_padding_max_log = 5000  # space for title
            plot.plot("NOSTACK E1")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            plot.save("%s/unfolded_unnormalised_%s_with_unreg_lambda_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))


class RecoPtBinnedPlotter(object):
    def __init__(self, setup, bins, hist_bin_chopper):
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
            subplot_limits=(0.75, 1.25),
        )

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
        return True

    def get_pt_bin_title(self, bin_edge_low, bin_edge_high):
        title = (("{jet_algo}\n"
                  "{region_label} region\n"
                  "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                 .format(
                    jet_algo=self.setup.jet_algo,
                    region_label=self.region['label'],
                    pt_str=self.setup.pt_str,
                    bin_edge_low=bin_edge_low,
                    bin_edge_high=bin_edge_high
                ))
        return title

    def plot_detector_normalised(self, unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            self.hist_bin_chopper.add_obj("input_hist_bg_subtracted", unfolder.input_hist_bg_subtracted)
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_bg_subtracted', ibin, 'detector')
            self.hist_bin_chopper.add_obj("hist_mc_reco_bg_subtracted", unfolder.hist_mc_reco_bg_subtracted)
            mc_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, 'detector')

            entries = [
                Contribution(mc_hist_bin,
                             label="MC (bg-subtracted)",
                             line_color=self.plot_colours['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_mc_colour'], marker_size=0,
                             normalise_hist=False),
                Contribution(input_hist_bin,
                             label="Data (bg-subtracted)",
                             line_color=self.plot_colours['reco_data_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_data_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_hist_bin,
                             normalise_hist=False),
            ]
            if not self.check_entries(entries, "plot_detector_normalised %d" % (ibin)):
                continue
            plot = Plot(entries,
                        xtitle=self.setup.detector_title,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        what="hist",
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        has_data=self.setup.has_data,
                        subplot_type='ratio',
                        subplot_title='Data / MC',
                        subplot_limits=(0.75, 1.25),)
            self._modify_plot(plot)
            plot.plot("NOSTACK E1")
            plot.save("%s/detector_reco_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (self.setup.output_dir, self.setup.append, ibin, self.setup.output_fmt))

    def plot_folded_unfolded_normalised(self, unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            self.hist_bin_chopper.add_obj("input_hist_bg_subtracted", unfolder.input_hist_bg_subtracted)
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_bg_subtracted', ibin, 'detector')
            self.hist_bin_chopper.add_obj("folded_unfolded", unfolder.folded_unfolded)
            folded_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_unfolded', ibin, 'detector')

            entries = [
                Contribution(input_hist_bin,
                             label="Unfolding input (bg-subtracted)",
                             line_color=self.plot_colours['reco_unfolding_input_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_unfolding_input_colour'], marker_style=20, marker_size=0,
                             normalise_hist=False),
                Contribution(folded_unfolded_hist_bin,
                             label="Folded unfolded",
                             line_color=self.plot_colours['reco_folded_unfolded_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_folded_unfolded_colour'], marker_style=20, marker_size=0.75,
                             subplot=input_hist_bin,
                             normalise_hist=False),
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

    def plot_folded_unfolded_with_mc_normalised(self, unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            self.hist_bin_chopper.add_obj("hist_mc_reco_bg_subtracted", unfolder.hist_mc_reco_bg_subtracted)
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, 'detector')
            self.hist_bin_chopper.add_obj("input_hist_bg_subtracted", unfolder.input_hist_bg_subtracted)
            input_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('input_hist_bg_subtracted', ibin, 'detector')
            self.hist_bin_chopper.add_obj("folded_unfolded", unfolder.folded_unfolded)
            folded_unfolded_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_unfolded', ibin, 'detector')

            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_colours['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_mc_colour'], marker_style=20, marker_size=0,
                             normalise_hist=False),
                Contribution(input_hist_bin,
                             label="Unfolding input (bg-subtracted)",
                             line_color=self.plot_colours['reco_unfolding_input_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_unfolding_input_colour'], marker_style=20, marker_size=0,
                             subplot=mc_reco_hist_bin,
                             normalise_hist=False),
                Contribution(folded_unfolded_hist_bin,
                             label="Folded unfolded",
                             line_color=self.plot_colours['reco_folded_unfolded_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_folded_unfolded_colour'],marker_style=20, marker_size=0.75,
                             subplot=mc_reco_hist_bin,
                             normalise_hist=False),
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

    def plot_folded_gen_normalised(self, unfolder):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            self.hist_bin_chopper.add_obj("hist_mc_reco_bg_subtracted", unfolder.hist_mc_reco_bg_subtracted)
            mc_reco_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_mc_reco_bg_subtracted', ibin, 'detector')
            self.hist_bin_chopper.add_obj("folded_mc_truth", unfolder.folded_mc_truth)
            folded_mc_truth_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('folded_mc_truth', ibin, 'detector')

            entries = [
                Contribution(mc_reco_hist_bin,
                             label="MC (reco, bg-subtracted)",
                             line_color=self.plot_colours['reco_mc_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_mc_colour'], marker_style=20, marker_size=0,
                             normalise_hist=False),
                Contribution(folded_mc_truth_hist_bin,
                             label="Folded gen",
                             line_color=self.plot_colours['reco_folded_mc_truth_colour'], line_width=self.line_width,
                             marker_color=self.plot_colours['reco_folded_mc_truth_colour'], marker_style=20, marker_size=0.75,
                             subplot=mc_reco_hist_bin,
                             normalise_hist=False),
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory (should be the one made by unfolding.py")
    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT.keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them (default).",
                        default=["all"])
    parser.add_argument("--outputDir",
                        default=None,
                        help='Output directory (default is the source dir')

    region_group = parser.add_argument_group('Region selection')
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

    regions = setup_regions(args)

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    jet_algo = "AK4 PF PUPPI"
    if "ak8puppi" in args.source:
        jet_algo = "AK8 PF PUPPI"

    has_data = not ('_MC_all' in args.source or '_MC_split' in args.source)

    # Store all things for final summary plots
    all_binned_hists = []

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
            unfolder = unpack_dict['unfolder']
            alt_unfolder = unpack_dict['alt_unfolder']
            alt_hist_truth = unpack_dict['alt_hist_truth']

            # Note that experimental systs are only different response matrices, and are stored in the main unfolder
            # we have
            # - unfolder
            # - unreg_unfolder
            # - alt_unfolder
            # - model_systematics[name] unfolder
            # - pdf[ind] unfolder

            # MAKE ALL THE PLOTS
            # ------------------------------------------------------------------

            # Big 1D plots to compare things
            hbc = HistBinChopper(unfolder)
            hbc.add_obj("unfolded", unfolder.unfolded)
            hbc.add_obj("unfolded_stat_err", unfolder.unfolded_stat_err)
            hbc.add_obj("hist_truth", unfolder.hist_truth)

            setup = Setup(jet_algo=jet_algo,
                          region=region,
                          angle=angle,
                          output_dir=angle_output_dir,
                          has_data=has_data)

            # Iterate through pt bins - gen binning
            gen_pt_binned_plotter = GenPtBinnedPlotter(setup=setup,
                                                       bins=unfolder.pt_bin_edges_gen,
                                                       hist_bin_chopper=hbc)
            gen_pt_binned_plotter.plot_unfolded_unnormalised(unfolder)
            gen_pt_binned_plotter.plot_unfolded_normalised(unfolder)
            if alt_hist_truth:
                gen_pt_binned_plotter.plot_unfolded_with_alt_truth_normalised(unfolder=unfolder,
                                                                              alt_truth=alt_hist_truth)

            if unfolder.tau > 0 and unreg_unfolder:
                gen_pt_binned_plotter.plot_unfolded_with_unreg_normalised(unfolder=unfolder,
                                                                          unreg_unfolder=unreg_unfolder)

            if alt_unfolder:
                gen_pt_binned_plotter.plot_unfolded_with_alt_response_normalised(unfolder=unfolder,
                                                                                 alt_unfolder=alt_unfolder)
                gen_pt_binned_plotter.plot_unfolded_with_alt_response_truth_normalised(unfolder=unfolder,
                                                                                       alt_unfolder=alt_unfolder,
                                                                                       alt_truth=alt_hist_truth)

            if len(region['model_systematics']) > 0:
                print(region['model_systematics'])
                gen_pt_binned_plotter.plot_unfolded_with_model_systs_normalised(unfolder=unfolder,
                                                                                model_systs=region['model_systematics'])

            # if has_data:
            gen_pt_binned_plotter.plot_detector_normalised(unfolder)

            # Iterate through lambda bins - gen binning
            lambda_pt_binned_plotter = GenLambdaBinnedPlotter(setup=setup,
                                                              bins=unfolder.variable_bin_edges_gen,
                                                              hist_bin_chopper=hbc)
            lambda_pt_binned_plotter.plot_unfolded_unnormalised(unfolder)

            if unfolder.tau > 0 and unreg_unfolder:
                lambda_pt_binned_plotter.plot_unfolded_with_unreg_unnormalised(unfolder, unreg_unfolder)

            # Iterate through pt bins - reco binning
            reco_pt_binned_plotter = RecoPtBinnedPlotter(setup=setup,
                                                         bins=unfolder.pt_bin_edges_reco,
                                                         hist_bin_chopper=hbc)
            reco_pt_binned_plotter.plot_detector_normalised(unfolder)
            reco_pt_binned_plotter.plot_folded_unfolded_normalised(unfolder)
            reco_pt_binned_plotter.plot_folded_unfolded_with_mc_normalised(unfolder)
            reco_pt_binned_plotter.plot_folded_gen_normalised(unfolder)

            # Iterate through lambda bins - reco binning
            #
