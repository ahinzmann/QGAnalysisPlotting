"""Common functions for doing QG flavour info plots"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import bisect
import numpy as np
import os
from array import array
import warnings
from itertools import cycle

# My stuff
from comparator import Contribution, Plot, grab_obj, ZeroContributions
import common_utils as cu
from qg_common import *
import qg_general_plots as qgp
from qg_general_plots import get_projection_plot

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)

# monkey-patch warning formatter
warnings.formatwarning = cu._formatwarning

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


FLAV_STR_DICT = {
    'u': 'Up quark',
    'd': 'Down quark',
    'c': 'Charm quark',
    's': 'Strange quark',
    'b': 'Bottom quark',
    't': 'Top quark',
    'g': 'Gluon',
    '1-g': 'Non-gluon',
}

class FlavourInfoCache(object):
    """Cache flavour efficiencies, such that getting efficiencies from TFiles
    only needs to be done once per set of input args"""

    def __init__(self, is_rivet_src=False):
        self._cache = {}
        self.is_rivet_src = is_rivet_src

    def get_flavour_efficiencies(self, **kwargs):
        tuple_args = tuple(sorted(kwargs.items))
        print(tuple_args)
        if tuple_args in self._cache:
            return self._cache[tuple_args]

        if self.is_rivet_src:
            eff = get_flavour_efficiencies_rivet(**kwargs)
            self._cache[tuple_args] = eff
        else:
            eff = get_flavour_efficiencies(**kwargs)
            self._cache[tuple_args] = eff


def get_jet_str(thing):
    return "jet"
    if "gen" in thing.lower():
        return "genjet"
    else:
        return "jet"


def get_flavour_hist_name(dirname, var_prepend="", which_jet="both", metric="pt"):
    jet_str = "" if which_jet == "both" else str(which_jet)
    return "%s/%sjet%s_flavour_vs_%s" % (dirname, var_prepend, jet_str, metric)


def get_rivet_flavour_hist_name(dirname, region, radius_ind):
    return "%s/%s_flav_vs_pt_radius%s" % (dirname, region, radius_ind)


def get_flavour_efficiencies(input_file, hist_name, bins):
    """Get dict of flav : [TEfficiency for specified bins] for a given input file & hist name."""
    h2d_flav = grab_obj(input_file, hist_name)

    flav_dict = {'d': [], 'u': [], 's': [], 'c': [], 'b': [] ,'t': [], 'g': [], 'unknown': [], 'total': []}

    for (pt_min, pt_max) in bins:
        h_flav = get_projection_plot(h2d_flav, pt_min, pt_max, cut_axis='y')

        total_err = array('d', [-1.])
        total2 = h_flav.IntegralAndError(1, h_flav.GetNbinsX()+1, total_err)
        total_err = total_err[0]

        unknown_num = h_flav.GetBinContent(1)
        d_num = h_flav.GetBinContent(2)
        u_num = h_flav.GetBinContent(3)
        s_num = h_flav.GetBinContent(4)
        c_num = h_flav.GetBinContent(5)
        b_num = h_flav.GetBinContent(6)
        t_num = h_flav.GetBinContent(7)
        g_num = h_flav.GetBinContent(22)

        total = d_num+u_num+s_num+c_num+b_num+t_num+g_num+unknown_num
        # print(total, total2)
        if not cu.same_floats(total, total2):
            # raise RuntimeError("totals dont match: %.9g vs %.9g" % (total, total2))
            warnings.warn("totals dont match: %.9g vs %.9g" % (total, total2))
        if total <= 0:
            warnings.warn("Total = 0")
        if total2 <= 0:
            warnings.warn("Total2 = 0")

        unknown_err = h_flav.GetBinError(1)
        d_err = h_flav.GetBinError(2)
        u_err = h_flav.GetBinError(3)
        s_err = h_flav.GetBinError(4)
        c_err = h_flav.GetBinError(5)
        b_err = h_flav.GetBinError(6)
        t_err = h_flav.GetBinError(7)
        g_err = h_flav.GetBinError(22)

        # print(pt_min, pt_max, b_num, b_err, total2, total_err)
        # print(pt_min, pt_max, g_num, g_err, total2, total_err)

        flav_dict['unknown'].append([unknown_num, unknown_err])
        flav_dict['d'].append([d_num, d_err])
        flav_dict['u'].append([u_num, u_err])
        flav_dict['s'].append([s_num, s_err])
        flav_dict['c'].append([c_num, c_err])
        flav_dict['b'].append([b_num, b_err])
        flav_dict['t'].append([t_num, t_err])
        flav_dict['g'].append([g_num, g_err])
        flav_dict['total'].append([total, total_err])


    flav_eff_dict = {
        'd': None,
        'u': None,
        's': None,
        'c': None,
        'b': None,
        't': None,
        'g': None,
        'unknown': None,
    }

    bin_edges = [p[0] for p in bins]
    bin_edges += bins[-1][1:]
    bin_edges = array('d', bin_edges)

    h_total = ROOT.TH1D("h_total", ";p_{T} [GeV];N", len(bins), bin_edges)
    for ind, val in enumerate(flav_dict['total'], 1):
        h_total.SetBinContent(ind, val[0])
        h_total.SetBinError(ind, val[1])

    for flav_key in [k for k in flav_dict.keys() if k != 'total']:
        h_flav = ROOT.TH1D("h_%s" % flav_key, ";p_{T} [GeV];N", len(bins), bin_edges)
        for ind, val in enumerate(flav_dict[flav_key], 1):
            h_flav.SetBinContent(ind, val[0])
            h_flav.SetBinError(ind, val[1])

        if (h_flav.GetNbinsX() != h_total.GetNbinsX()):
            raise ValueError("flav %s has diff # bins to total (%d vs %d)" % flav_key, h_flav.GetNbinsX(), h_total.GetNbinsX())

        eff_flav = ROOT.TEfficiency(h_flav, h_total)
        eff_flav.SetStatisticOption(1)  # normal approx
        flav_eff_dict[flav_key] = eff_flav

    return flav_eff_dict


def compare_flavour_fractions_vs_pt(input_files, 
                                    dirnames, 
                                    pt_bins, 
                                    labels, 
                                    flav, 
                                    output_filename, 
                                    title="", 
                                    var_prepend="", 
                                    which_jet="both", 
                                    xtitle="p_{T}^{jet} [GeV]", 
                                    n_partons='all', 
                                    is_preliminary=True):
    """Plot a specified flavour fraction vs pT for several sources.
    Each entry in input_files, dirnames, and labels corresponds to one line
    n_partons can be a str, 'all', '1', etc, or a list of str to include

    TODO: fix this - bit stupid input format
    """
    bin_centers = [0.5*(x[0]+x[1]) for x in pt_bins]
    bin_widths = [0.5*(x[1]-x[0]) for x in pt_bins]

    if isinstance(n_partons, str):
        n_partons = [n_partons]

    contribs = []
    for n_parton_ind, n_parton in enumerate(n_partons):
        metric = 'pt'
        if n_parton.lower() != 'all':
            metric = 'pt_npartons_%s' % n_parton
        info = [get_flavour_efficiencies(ifile,
                                         bins=pt_bins,
                                         hist_name=get_flavour_hist_name(
                                                 dirname=dname,
                                                 var_prepend=var_prepend,
                                                 which_jet=(which_jet if "Dijet" in dname else "both"),
                                                 metric=metric)
                                         )
                for ifile, dname in zip(input_files, dirnames)]
        N = len(bin_centers)

        colours = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2]

        for i, fdict in enumerate(info):
            if flav in ['u', 'd', 's', 'c', 'b', 't', 'g']:
                obj = fdict[flav].CreateGraph()
            else:
                raise RuntimeError("Robin broke 1-X functionality")
                obj = ROOT.TGraphErrors(N, np.array(bin_centers), 1.-np.array(fdict[flav.replace("1-", '')]), np.array(bin_widths), np.zeros(N))
            if obj.GetN() == 0:
                continue
            n_parton_str = "" if n_parton == "all" else " (%s-parton)" % n_parton
            c = Contribution(obj,
                             label="%s%s" % (labels[i], n_parton_str),
                             line_color=colours[i]+n_parton_ind, line_width=1,
                             line_style=n_parton_ind+1,
                             marker_style=20+i, marker_color=colours[i]+n_parton_ind, marker_size=1,
                             leg_draw_opt="LP")
            contribs.append(c)

    flav_str = FLAV_STR_DICT[flav]
    ytitle = "Fraction of %s %ss" % (flav_str.lower(), get_jet_str(''))
    p = Plot(contribs,
             what='graph',
             xtitle=xtitle,
             ytitle=ytitle,
             title=title,
             xlim=(pt_bins[0][0], pt_bins[-1][1]),
             ylim=(0, 1),
             has_data=False,
             is_preliminary=is_preliminary)
    p.default_canvas_size = (600, 600)
    try:
        p.plot("APZ")
        p.main_pad.SetBottomMargin(0.16)
        p.get_modifier().GetXaxis().SetTitleOffset(1.4)
        p.get_modifier().GetXaxis().SetTitleSize(.045)
        p.legend.SetX1(0.5)
        p.legend.SetY1(0.65)
        p.legend.SetY2(0.87)
        p.set_logx(do_more_labels=True, do_exponent=False)
        p.save(output_filename)
    except ZeroContributions:
        pass


def create_contibutions_compare_vs_pt(input_files, hist_names, pt_bins, labels, flav, **contrib_kwargs):
    bin_centers = [0.5*(x[0]+x[1]) for x in pt_bins]
    bin_widths = [0.5*(x[1]-x[0]) for x in pt_bins]

    contribs = []
    info = [get_flavour_efficiencies(ifile, bins=pt_bins, hist_name=hname)
            for ifile, hname in zip(input_files, hist_names)]
    N = len(bin_centers)

    colours = cycle([ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange-3, ROOT.kViolet+6])

    for i, (fdict, label, col) in enumerate(zip(info, labels, colours)):
        if flav in ['u', 'd', 's', 'c', 'b', 't', 'g']:
            obj = fdict[flav].CreateGraph()
        else:
            raise RuntimeError("Robin broke 1-X functionality")
            obj = ROOT.TGraphErrors(N, np.array(bin_centers), 1.-np.array(fdict[flav.replace("1-", '')]), np.array(bin_widths), np.zeros(N))
        if obj.GetN() == 0:
            continue
        c = Contribution(obj,
                         label=label,
                         line_color=col, line_width=1,
                         line_style=1,
                         marker_style=20+i, marker_color=col, marker_size=1,
                         leg_draw_opt="LP",
                         **contrib_kwargs)
        contribs.append(c)
    return contribs


def compare_flavour_fraction_hists_vs_pt_from_contribs(contribs, flav, output_filename, title="", xtitle="p_{T}^{jet} [GeV]", **plot_kwargs):
    """Plot a specified flavour fraction vs pT for several sources.

    TODO: use this one more often - compare_flavour_fractions_vs_pt() is almost identical but has to deal with npartons
    """
    flav_str = FLAV_STR_DICT[flav]
    ytitle = "Fraction of %s %ss" % (flav_str.lower(), get_jet_str(''))
    p = Plot(contribs,
             what='graph',
             xtitle=xtitle,
             ytitle=ytitle,
             title=title,
             xlim=(50, 2000),
             ylim=(0, 1),
             has_data=False,
             **plot_kwargs)
    p.default_canvas_size = (600, 600)
    try:
        p.plot("AP")
        p.main_pad.SetBottomMargin(0.16)
        p.get_modifier().GetXaxis().SetTitleOffset(1.4)
        p.get_modifier().GetXaxis().SetTitleSize(.045)
        p.legend.SetX1(0.56)
        p.legend.SetY1(0.65)
        p.legend.SetY2(0.87)
        if len(contribs) >=4:
            p.legend.SetY1(0.7)
            p.legend.SetX1(0.5)
            p.legend.SetNColumns(2)
        p.set_logx(do_more_labels=True, do_exponent=False)
        p.save(output_filename)
    except ZeroContributions:
        warnings.warn("No contributions for %s" % output_filename)


def compare_flavour_fraction_hists_vs_pt(input_files, hist_names, pt_bins, labels, flav, output_filename, title="", xtitle="p_{T}^{jet} [GeV]", is_preliminary=True):
    """Plot a specified flavour fraction vs pT for several sources.
    Each entry in input_files, dirnames, and labels corresponds to one line
    n_partons can be a str, 'all', '1', etc, or a list of str to include

    TODO: use this one more often - compare_flavour_fractions_vs_pt() is almost identical but has to deal with npartons
    """
    bin_centers = [0.5*(x[0]+x[1]) for x in pt_bins]
    bin_widths = [0.5*(x[1]-x[0]) for x in pt_bins]

    contribs = []
    info = [get_flavour_efficiencies(ifile, bins=pt_bins, hist_name=hname)
            for ifile, hname in zip(input_files, hist_names)]
    N = len(bin_centers)

    colours = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2]

    for i, fdict in enumerate(info):
        if flav in ['u', 'd', 's', 'c', 'b', 't', 'g']:
            obj = fdict[flav].CreateGraph()
        else:
            raise RuntimeError("Robin broke 1-X functionality")
            obj = ROOT.TGraphErrors(N, np.array(bin_centers), 1.-np.array(fdict[flav.replace("1-", '')]), np.array(bin_widths), np.zeros(N))
        if obj.GetN() == 0:
            continue
        c = Contribution(obj,
                         label=labels[i],
                         line_color=colours[i], line_width=1,
                         line_style=1,
                         marker_style=20+i, marker_color=colours[i], marker_size=1,
                         leg_draw_opt="LP")
        contribs.append(c)

    flav_str = FLAV_STR_DICT[flav]
    ytitle = "Fraction of %s %ss" % (flav_str.lower(), get_jet_str(''))
    p = Plot(contribs,
             what='graph',
             xtitle=xtitle,
             ytitle=ytitle,
             title=title,
             xlim=(50, 2000),
             ylim=(0, 1),
             has_data=False,
             is_preliminary=is_preliminary)
    p.default_canvas_size = (600, 600)
    try:
        p.plot("AP")
        p.main_pad.SetBottomMargin(0.16)
        p.get_modifier().GetXaxis().SetTitleOffset(1.4)
        p.get_modifier().GetXaxis().SetTitleSize(.045)
        p.legend.SetX1(0.56)
        p.legend.SetY1(0.65)
        p.legend.SetY2(0.87)
        p.set_logx(do_more_labels=True, do_exponent=False)
        p.save(output_filename)
    except ZeroContributions:
        pass

def do_flavour_fraction_vs_pt(input_file, hist_name, pt_bins, output_filename, title=""):
    """Plot all flavour fractions vs PT for one input file & hist_name in the ROOT file"""
    info = get_flavour_efficiencies(input_file, bins=pt_bins, hist_name=hist_name)

    leg_draw_opt = "LP"
    plot_u = Contribution(info['u'].CreateGraph(), label="Up", line_color=ROOT.kRed, marker_color=ROOT.kRed, marker_style=20, leg_draw_opt=leg_draw_opt)
    plot_d = Contribution(info['d'].CreateGraph(), label="Down", line_color=ROOT.kBlue, marker_color=ROOT.kBlue, marker_style=21, leg_draw_opt=leg_draw_opt)
    plot_s = Contribution(info['s'].CreateGraph(), label="Strange", line_color=ROOT.kBlack, marker_color=ROOT.kBlack, marker_style=22, leg_draw_opt=leg_draw_opt)
    plot_c = Contribution(info['c'].CreateGraph(), label="Charm", line_color=ROOT.kGreen+2, marker_color=ROOT.kGreen+2, marker_style=23, leg_draw_opt=leg_draw_opt)
    plot_b = Contribution(info['b'].CreateGraph(), label="Bottom", line_color=ROOT.kOrange-3, marker_color=ROOT.kOrange-3, marker_style=33, leg_draw_opt=leg_draw_opt)
    plot_g = Contribution(info['g'].CreateGraph(), label="Gluon", line_color=ROOT.kViolet, marker_color=ROOT.kViolet, marker_style=29, leg_draw_opt=leg_draw_opt)
    plot_unknown = Contribution(info['unknown'].CreateGraph(), label="Unknown", line_color=ROOT.kGray+1, marker_color=ROOT.kGray+1, marker_style=26, leg_draw_opt=leg_draw_opt)

    p_flav = Plot([plot_d, plot_u, plot_s, plot_c, plot_b, plot_g, plot_unknown],
                  what='graph',
                  xtitle="p_{T}^{%s} [GeV]" % get_jet_str(''),
                  ytitle="Fraction",
                  title=title,
                  xlim=(pt_bins[0][0], pt_bins[-1][1]),
                  ylim=[0, 1],
                  has_data=False)
    p_flav.default_canvas_size = (600, 600)
    p_flav.plot("AP")
    p_flav.main_pad.SetBottomMargin(0.16)
    p_flav.get_modifier().GetXaxis().SetTitleOffset(1.4)
    p_flav.get_modifier().GetXaxis().SetTitleSize(.045)
    p_flav.set_logx(do_more_labels=True, do_exponent=False)
    p_flav.legend.SetX1(0.55)
    p_flav.legend.SetY1(0.72)
    p_flav.legend.SetY2(0.85)
    p_flav.legend.SetNColumns(2)
    p_flav.save(output_filename)


def do_flavour_fraction_vs_eta(input_file, dirname, eta_bins, output_filename, title="", var_prepend="", which_jet="both", append=""):
    """Plot all flavour fractions vs eta for one input file & dirname in the ROOT file"""
    info = get_flavour_efficiencies(input_file,
                                    bins=eta_bins,
                                    hist_name=get_flavour_hist_name(
                                        dirname=dirname,
                                        var_prepend=var_prepend,
                                        which_jet="both", # since no jet/jet1/jet2 in hist name
                                        metric='eta'+append)
                                    )

    leg_draw_opt = "LP"
    plot_u = Contribution(info['u'].CreateGraph(), label="Up", line_color=ROOT.kRed, marker_color=ROOT.kRed, marker_style=20, leg_draw_opt=leg_draw_opt)
    plot_d = Contribution(info['d'].CreateGraph(), label="Down", line_color=ROOT.kBlue, marker_color=ROOT.kBlue, marker_style=21, leg_draw_opt=leg_draw_opt)
    plot_s = Contribution(info['s'].CreateGraph(), label="Strange", line_color=ROOT.kBlack, marker_color=ROOT.kBlack, marker_style=22, leg_draw_opt=leg_draw_opt)
    plot_c = Contribution(info['c'].CreateGraph(), label="Charm", line_color=ROOT.kGreen+2, marker_color=ROOT.kGreen+2, marker_style=23, leg_draw_opt=leg_draw_opt)
    plot_b = Contribution(info['b'].CreateGraph(), label="Bottom", line_color=ROOT.kOrange-3, marker_color=ROOT.kOrange-3, marker_style=33, leg_draw_opt=leg_draw_opt)
    plot_g = Contribution(info['g'].CreateGraph(), label="Gluon", line_color=ROOT.kViolet, marker_color=ROOT.kViolet, marker_style=29, leg_draw_opt=leg_draw_opt)
    plot_unknown = Contribution(info['unknown'].CreateGraph(), label="Unknown", line_color=ROOT.kGray+1, marker_color=ROOT.kGray+1, marker_style=26, leg_draw_opt=leg_draw_opt)

    p_flav = Plot([plot_d, plot_u, plot_s, plot_c, plot_b, plot_g, plot_unknown],
                  what='graph',
                  xtitle="y^{%s}" % get_jet_str(''),
                  ytitle="Fraction",
                  title=title,
                  xlim=(eta_bins[0][0], eta_bins[-1][1]),
                  ylim=[0, 1],
                  has_data=False)
    p_flav.default_canvas_size = (600, 600)
    p_flav.plot("AP")
    p_flav.main_pad.SetBottomMargin(0.16)
    p_flav.get_modifier().GetXaxis().SetTitleOffset(1.4)
    p_flav.get_modifier().GetXaxis().SetTitleSize(.045)
    p_flav.legend.SetX1(0.55)
    p_flav.legend.SetY1(0.72)
    p_flav.legend.SetY2(0.85)
    p_flav.legend.SetNColumns(2)
    p_flav.save(output_filename)
