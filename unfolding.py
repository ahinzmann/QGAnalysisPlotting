#!/usr/bin/env python


"""TUnfold it all

Thanks to Ashley, Dennis
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
from array import array
import numpy as np
import math
import distutils
from distutils import util
from itertools import product

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from unfolding_classes import TauScanner, LCurveScanner, MyUnfolder, MyUnfolderPlotter


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")
ROOT.gStyle.SetHistTopMargin(0.)


# Control plot output format
OUTPUT_FMT = "pdf"


def calculate_chi2(hist_test, hist_truth, hist_covariance=None):
    """Calculate the chi2 between 2 histograms, given a covariance matrix

    = (hist1 - hist2)T * cov^-1 * (hist1 - hist2)

    Parameters
    ----------
    hist1 : list, np.array
        List of bin values for one histogram
    hist2 : list, np.array
        List of bin values for other histogram
    cov : np.array
        Covariance matrix

    Returns
    -------
    float
    """
    diff = hist1 - hist2
    diff = diff.reshape(len(diff), 1)
    # for now, hack the inversion, since we know it's diagonal
    inv_cov = np.zeros_like(cov)
    # if True:
    #     for i, e in enumerate(range(cov.shape[0])):
    #         cov_entry = cov[i][i]
    #         if cov_entry != 0:
    #             inv_cov[i][i] = 1./cov_entry
    #         else:
    #             inv_cov[i][i] = 0
    # else:
    inv_cov = np.linalg.inv(cov)
    # print('inv_cov', inv_cov)
    part = np.dot(inv_cov, diff)
    # print('part', part)
    result = np.dot(diff.T, part)
    # print('result', result)
    return result[0][0]


def create_hist_with_errors(hist, err_matrix):
    hnew = hist.Clone(cu.get_unique_str())
    nbins = hist.GetNbinsX()
    for i in range(1, nbins+1):
        err = math.sqrt(err_matrix.GetBinContent(i, i))
        hnew.SetBinError(i, err)
    return hnew


def make_hist_from_diagonal_errors(h2d, do_sqrt=True):
    nbins = h2d.GetNbinsX()
    hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, 0, nbins)
    for i in range(1, nbins+1):
        err = h2d.GetBinContent(i, i)
        if do_sqrt and err > 0:
            err = math.sqrt(err)
        hnew.SetBinError(i, err)
    return hnew


def update_hist_bin_content(h_orig, h_to_be_updated):
    if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
        raise RuntimeError("Need same # x bins, %d vs %s" % (h_orig.GetNbinsX(), h_to_be_updated.GetNbinsX()))
    for i in range(0, h_orig.GetNbinsX()+2):
        h_to_be_updated.SetBinContent(i, h_orig.GetBinContent(i))
        # h_to_be_updated.SetBinError(i, h_orig.GetBinError(i))


def update_hist_bin_error(h_orig, h_to_be_updated):
    if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
        raise RuntimeError("Need same # x bins, %d vs %s" % (h_orig.GetNbinsX(), h_to_be_updated.GetNbinsX()))
    for i in range(0, h_orig.GetNbinsX()+2):
        h_to_be_updated.SetBinError(i, h_orig.GetBinError(i))


def draw_projection_comparison(h_orig, h_projection, title, xtitle, output_filename, print_bin_comparison=True):
    """Draw 2 hists, h_orig the original, and h_projection the projection of a 2D hist"""

    # Check integrals
    int_orig = h_orig.Integral()
    int_proj = h_projection.Integral()
    if abs(int_orig - int_proj)/int_orig > 0.01:
        print("draw_projection_comparison: different integrals: %f vs %f" % (int_orig, int_proj))

    # Check bin-by-bin
    if print_bin_comparison:
        for i in range(1, h_orig.GetNbinsX()+1):
            value_orig = h_orig.GetBinContent(i)
            value_proj = h_projection.GetBinContent(i)
            if value_orig == 0 and value_proj == 0:
                continue
            rel_diff = abs((value_orig - value_proj)/max(abs(value_orig), abs(value_proj)))
            if rel_diff > 1E-3:
                print("draw_projection_comparison: bin %s has different contents: %f vs %f (rel diff %f)" % (i, value_orig, value_proj, rel_diff))

    entries = [
        Contribution(h_orig, label="1D hist",
                     line_color=ROOT.kBlue, line_width=1,
                     marker_color=ROOT.kBlue, marker_size=0,
                     normalise_hist=False),
        Contribution(h_projection, label="Response map projection",
                     line_color=ROOT.kRed, line_width=1,
                     marker_color=ROOT.kRed, marker_size=0,
                     normalise_hist=False,
                     subplot=h_orig),
    ]
    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle=xtitle,
                ytitle="N",
                subplot_type='ratio',
                subplot_title='Projection / 1D',
                subplot_limits=(0.999, 1.001))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HIST")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [h_orig, h_projection])
    plot.container.SetMaximum(ymax * 10)
    # plot.container.SetMinimum(1E-8)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX2NDC(0.85)
    plot.save(output_filename)


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


def plot_uncertainty_shifts(total_hist, stat_hist, syst_shifts, systs, output_filename, title, angle_str):
    """
    Parameters
    ----------
    total_hist : TH1
        Distribution with total uncertainty
    stat_hist : TH1
        Distribution with stat-only uncertainties
    syst_shifts : list[TH1]
        Distributions with 1-sigma shift from systematics
    systs : list[dict]
        Dicts describing each syst
    output_filename : str
        Output plot filename
    title : str
        Title to put on plot
    angle_str : str
        Angle name for x axis

    Returns
    -------
    None
        If all hists are empty
    """
    entries = []
    hists = []
    for h, syst_dict in zip(syst_shifts, systs):
        h_fraction = h.Clone()
        h_fraction.Divide(total_hist)
        hists.append(h_fraction)
        for i in range(1, h_fraction.GetNbinsX()+1):
            h_fraction.SetBinContent(i, abs(h_fraction.GetBinContent(i)))
        c = Contribution(h_fraction,
                         label=syst_dict['label'],
                         line_color=syst_dict['colour'],
                         line_style=syst_dict.get('linestyle', 1),
                         line_width=2,
                         marker_size=0,
                         marker_color=syst_dict['colour'],
                         )
        entries.append(c)

    # Add systematic
    h_syst = stat_hist.Clone()
    h_total = total_hist.Clone()
    for i in range(1, h_syst.GetNbinsX()+1):
        if total_hist.GetBinContent(i) > 0:
            h_syst.SetBinContent(i, stat_hist.GetBinError(i) / total_hist.GetBinContent(i))
            h_total.SetBinContent(i, total_hist.GetBinError(i) / total_hist.GetBinContent(i))
        else:
            h_syst.SetBinContent(i, 0)
            h_total.SetBinContent(i, 0)
        h_syst.SetBinError(i, 0)
        h_total.SetBinError(i, 0)
    c_stat = Contribution(h_syst,
                         label="Stat.",
                         line_color=ROOT.kRed,
                         line_style=3,
                         line_width=3,
                         marker_size=0,
                         marker_color=ROOT.kRed,
                         )
    entries.append(c_stat)
    c_tot = Contribution(h_total,
                         label="Total",
                         line_color=ROOT.kBlack,
                         line_style=1,
                         line_width=3,
                         marker_size=0,
                         marker_color=ROOT.kBlack,
                         )
    entries.append(c_tot)

    if not check_entries(entries, "systematic shifts"):
        return
    plot = Plot(entries,
                what="hist",
                title=title,
                xtitle="Particle-level "+angle_str,
                ytitle="| Fractional shift |")
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(0.68)
    plot.legend.SetX2(0.98)
    plot.legend.SetY2(0.88)
    plot.legend.SetNColumns(2)
    plot.y_padding_max_linear = 1.4
    plot.plot("NOSTACK HIST")
    plot.save(output_filename)

    plot.y_padding_max_log = 50
    plot.set_logy()
    log_filename, ext = os.path.splitext(output_filename)
    plot.save(log_filename+"_log"+ext)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory with ROOT files")

    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT.keys()),
                        nargs='+',
                        help="Lambda angles to unfold")

    parser.add_argument("--doSummaryPlot",
                        type=lambda x:bool(distutils.util.strtobool(x)),
                        default=False,
                        help='Do summary plot')

    parser.add_argument("--outputDir",
                        default='',
                        help='Output directory')

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

    regularization_group = parser.add_argument_group('Regularization options')
    regularization_group.add_argument("--regularize",
                                      choices=['None', 'tau', 'L'],
                                      default='None',
                                      help='Regularization scheme')
    regularization_group.add_argument("--regularizeAxis",
                                      choices=['both', 'pt', 'angle'],
                                      default='both',
                                      help='Axis to regularize')
    regularization_group.add_argument("--nScan",
                                      type=int,
                                      default=100,
                                      help='Number of scan points for regularization')
    regularization_group.add_argument("--biasFactor",
                                      type=float,
                                      default=0,
                                      help='Bias factor for regularization')

    mc_group = parser.add_argument_group('MC input options')
    mc_group.add_argument("--MCinput",
                          type=lambda x:bool(distutils.util.strtobool(x)),
                          default=False,
                          help='Unfold MC instead of data')

    mc_group.add_argument("--MCsplit",
                          type=lambda x:bool(distutils.util.strtobool(x)),
                          default=False,
                          help='Split MC between response & 1D reco, good for testing procedure')

    bg_group = parser.add_argument_group("Backgrounds options")
    bg_group.add_argument("--subtractBackgrounds",
                          type=lambda x:bool(distutils.util.strtobool(x)),
                          default=False,
                          help='Subtract true backgrounds (e.g. ttbar)')

    syst_group = parser.add_argument_group('Systematics options')
    syst_group.add_argument("--doExperimentalSysts",
                            type=lambda x:bool(distutils.util.strtobool(x)),
                            default=False,
                            help='Do experimental systematics (i.e. those that modify response matrix)')

    syst_group.add_argument("--doModelSysts",
                            type=lambda x:bool(distutils.util.strtobool(x)),
                            default=False,
                            help='Do model systematics (i.e. those that modify input to be unfolded)')

    syst_group.add_argument("--doPDFSysts",
                            type=lambda x:bool(distutils.util.strtobool(x)),
                            default=False,
                            help='Do pdf systematics (may be slow!)')

    syst_group.add_argument("--useAltResponse",
                            type=lambda x:bool(distutils.util.strtobool(x)),
                            default=False,
                            help='Use alternate response matrix to unfold')

    args = parser.parse_args()
    print(args)

    if not any([args.doDijetCentral, args.doDijetForward, args.doDijetCentralGroomed, args.doDijetForwardGroomed, args.doZPJ, args.doZPJGroomed]):
        raise RuntimeError("You need to specify at least one signal region e.g. --doDijetCentral")

    if not args.MCinput and args.doModelSysts:
        raise RuntimeError("You cannot do both model systs and run on data")

    if args.MCinput and args.subtractBackgrounds:
        print("")
        print("!!!! Cannot subtract backgrounds while using MC input, ignoring for now")
        print("")
        args.subtractBackgrounds = False

    # if args.useAltResponse and args.doExperimentalSysts:
    #     args.doExperimentalSysts = False
    #     print("You cannot use both --useAltResponse and --doExperimentalSysts: disabling doExperimentalSysts")

    # # TODO handle both?
    # if args.useAltResponse and args.doModelSysts:
    #     args.doExperimentalSysts = False
    #     print("You cannot use both --useAltResponse and --doModelSysts: disabling doModelSysts")

    # Setup files and regions to unfold
    # --------------------------------------------------------------------------
    src_dir = args.source
    src_dir_systs = os.path.join(src_dir, "systematics_files")

    regions = []

    if any([args.doDijetCentral, args.doDijetForward, args.doDijetCentralGroomed, args.doDijetForwardGroomed]):
        # FOR DIJET:
        input_mc_qcd_mgpythia_tfile = os.path.join(src_dir, qgc.QCD_FILENAME)
        input_mc_qcd_pythia_tfile = cu.open_root_file(os.path.join(src_dir, qgc.QCD_PYTHIA_ONLY_FILENAME))
        input_mc_qcd_herwig_tfile = os.path.join(src_dir, qgc.QCD_HERWIG_FILENAME)
        input_mc_qcd_herwig_tfile_reweight = os.path.join(src_dir, "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD_PtReweight.root")

        input_jetht_tfile = os.path.join(src_dir, qgc.JETHT_ZB_FILENAME)

        # actually these are all pretty similar...
        tau_limits_central = {
            'jet_puppiMultiplicity': (1E-9, 1E-6),
            'jet_pTD': (1E-12, 1E-8),
            'jet_LHA': (1E-11, 1E-8) if args.regularizeAxis == 'angle' else (1E-13, 1E-10),
            'jet_width': (1E-12, 1E-8),
            'jet_thrust': (1E-12, 1E-8),
            'jet_puppiMultiplicity_charged': (1E-12, 1E-8),
            'jet_pTD_charged': (1E-12, 1E-8),
            'jet_LHA_charged': (1E-10, 1E-8),
            'jet_width_charged': (1E-13, 1E-8),
            'jet_thrust_charged': (1E-12, 1E-9),
        }

        tau_limits_central_groomed = {
            'jet_puppiMultiplicity': (1E-9, 1E-6),
            'jet_pTD': (1E-12, 1E-8),
            'jet_LHA': (1E-11, 1E-8) if args.regularizeAxis == 'angle' else (1E-13, 1E-10),
            'jet_width': (1E-12, 1E-8),
            'jet_thrust': (1E-12, 1E-8),
            'jet_puppiMultiplicity_charged': (1E-12, 1E-8),
            'jet_pTD_charged': (1E-12, 1E-8),
            'jet_LHA_charged': (1E-10, 1E-8),
            'jet_width_charged': (1E-13, 1E-8),
            'jet_thrust_charged': (1E-12, 1E-9),
        }

        tau_limits_forward = {
            'jet_puppiMultiplicity': (1E-9, 1E-6),
            'jet_pTD': (1E-12, 1E-8),
            'jet_LHA': (1E-11, 1E-8) if args.regularizeAxis == 'angle' else (1E-13, 1E-10),
            'jet_width': (1E-12, 1E-8),
            'jet_thrust': (1E-12, 1E-8),
            'jet_puppiMultiplicity_charged': (1E-12, 1E-8),
            'jet_pTD_charged': (1E-12, 1E-8),
            'jet_LHA_charged': (1E-10, 1E-8),
            'jet_width_charged': (1E-13, 1E-8),
            'jet_thrust_charged': (1E-12, 1E-9),
        }

        tau_limits_forward_groomed = {
            'jet_puppiMultiplicity': (1E-9, 1E-6),
            'jet_pTD': (1E-12, 1E-8),
            'jet_LHA': (1E-11, 1E-8) if args.regularizeAxis == 'angle' else (1E-13, 1E-10),
            'jet_width': (1E-12, 1E-8),
            'jet_thrust': (1E-12, 1E-8),
            'jet_puppiMultiplicity_charged': (1E-12, 1E-8),
            'jet_pTD_charged': (1E-12, 1E-8),
            'jet_LHA_charged': (1E-10, 1E-8),
            'jet_width_charged': (1E-13, 1E-8),
            'jet_thrust_charged': (1E-12, 1E-9),
        }

        dijet_region_dict_template = {
            "name": "Dijet",
            "dirname": "Dijet_QG_Unfold_central_tighter",
            "label": "Dijet",
            "data_tfile": input_jetht_tfile,
            "mc_tfile": input_mc_qcd_mgpythia_tfile,
            "mc_label": "MG+Pythia8",
            # "mc_tfile": input_mc_qcd_herwig_tfile,
            # "mc_label": "Herwig++",
            "alt_mc_tfile": input_mc_qcd_herwig_tfile,
            "alt_mc_label": "Herwig++",
            # "alt_mc_tfile": input_mc_qcd_herwig_tfile_reweight,
            # "alt_mc_label": "Herwig++ (p_{T} reweight)",
            "tau_limits": tau_limits_central,
            "experimental_systematics": [
                {
                    "label": "Neutral hadron up",
                    "tfile": os.path.join(src_dir_systs, 'neutralHadronShiftUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kOrange-3,
                },
                {
                    "label": "Neutral hadron down",
                    "tfile": os.path.join(src_dir_systs, 'neutralHadronShiftDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kOrange-3,
                    "linestyle": 2,
                },
                {
                    "label": "Photon up",
                    "tfile": os.path.join(src_dir_systs, 'photonShiftUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kMagenta-3,
                },
                {
                    "label": "Photon down",
                    "tfile": os.path.join(src_dir_systs, 'photonShiftDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kMagenta-3,
                    "linestyle": 2,
                },
                {
                    "label": "JEC up",
                    "tfile": os.path.join(src_dir_systs, 'jecsmear_directionUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kGreen+2,
                },
                {
                    "label": "JEC down",
                    "tfile": os.path.join(src_dir_systs, 'jecsmear_directionDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kGreen+2,
                    "linestyle": 2,
                },
                {
                    "label": "JER up",
                    "tfile": os.path.join(src_dir_systs, 'jersmear_directionUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kOrange+3,
                },
                {
                    "label": "JER down",
                    "tfile": os.path.join(src_dir_systs, 'jersmear_directionDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kOrange+3,
                    "linestyle": 2,
                },
                {
                    "label": "Pileup up",
                    "tfile": os.path.join(src_dir_systs, 'pileup_directionUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kBlue-4,
                },
                {
                    "label": "Pileup down",
                    "tfile": os.path.join(src_dir_systs, 'pileup_directionDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kBlue-4,
                    "linestyle": 2,
                },
            ],
            "model_systematics": [
                {
                    "label": "muR up, muF nominal",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure,
                },
                {
                    "label": "muR down, muF nominal",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+1,
                },
                {
                    "label": "muR nominal, muF up",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+2,
                },
                {
                    "label": "muR nominal, muF down",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+3,
                },
                {
                    "label": "muR down, muF down",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+4,
                },
                {
                    "label": "muR up, muF up",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+5,
                },
                {
                    "label": "Herwig++",
                    "tfile": input_mc_qcd_herwig_tfile,
                    "colour": ROOT.kOrange-3,
                },
            ],
            "pdf_systematics": [
                {
                    "label": "PDF",  # this is a tempalte entry, used for future
                    "tfile": os.path.join(src_dir_systs, 'PDFvariationsTrue', qgc.QCD_FILENAME),
                    "colour": ROOT.kCyan+2,
                    "variations": range(100),  # list of all the variation #s to be used
                },
            ]
        }

        if args.doDijetCentral:
            dijet_region_central_dict = dijet_region_dict_template.copy()
            dijet_region_central_dict['tau_limits'] = tau_limits_central
            dijet_region_central_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter'
            dijet_region_central_dict['label'] = 'Dijet central'
            dijet_region_central_dict['name'] = 'Dijet_central'
            regions.append(dijet_region_central_dict)

        if args.doDijetForward:
            dijet_region_forward_dict = dijet_region_dict_template.copy()
            dijet_region_forward_dict['tau_limits'] = tau_limits_forward
            dijet_region_forward_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter'
            dijet_region_forward_dict['label'] = 'Dijet forward'
            dijet_region_forward_dict['name'] = 'Dijet_forward'
            regions.append(dijet_region_forward_dict)

        if args.doDijetCentralGroomed:
            dijet_region_central_groomed_dict = dijet_region_dict_template.copy()
            dijet_region_central_groomed_dict['tau_limits'] = tau_limits_central_groomed
            dijet_region_central_groomed_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter_groomed'
            dijet_region_central_groomed_dict['label'] = 'Dijet central'
            dijet_region_central_groomed_dict['name'] = 'Dijet_central_groomed'
            regions.append(dijet_region_central_groomed_dict)

        if args.doDijetForwardGroomed:
            dijet_region_forward_groomed_dict = dijet_region_dict_template.copy()
            dijet_region_forward_groomed_dict['tau_limits'] = tau_limits_forward_groomed
            dijet_region_forward_groomed_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter_groomed'
            dijet_region_forward_groomed_dict['label'] = 'Dijet forward'
            dijet_region_forward_groomed_dict['name'] = 'Dijet_forward_groomed'
            regions.append(dijet_region_forward_groomed_dict)

    if any([args.doZPJ, args.doZPJGroomed]):
        # FOR Z+JETS:
        input_mc_dy_mgpythia_tfile = os.path.join(src_dir, qgc.DY_FILENAME)
        input_mc_dy_mgherwig_tfile = os.path.join(src_dir, qgc.DY_MG_HERWIG_FILENAME)
        input_mc_dy_herwig_tfile = os.path.join(src_dir, qgc.DY_HERWIG_FILENAME)

        input_singlemu_tfile = os.path.join(src_dir, qgc.SINGLE_MU_FILENAME)

        tau_limits = {
            'jet_puppiMultiplicity': (1E-5, 1E-2),
            'jet_pTD': (1E-6, 1E-4),
            'jet_LHA': (1E-5, 1E-3),
            'jet_width': (1E-5, 1E-2),
            'jet_thrust': (1E-6, 1E-2),
            'jet_puppiMultiplicity_charged': (1E-6, 1E-2),
            'jet_pTD_charged': (1E-6, 1E-2),
            'jet_LHA_charged': (1E-5, 1E-2),
            'jet_width_charged': (1E-6, 1E-3),
            'jet_thrust_charged': (1E-8, 1E-5),
        }

        tau_limits_groomed = {
            'jet_puppiMultiplicity': (1E-5, 1E-2),
            'jet_pTD': (1E-7, 1E-3),
            'jet_LHA': (1E-5, 1E-3),
            'jet_width': (1E-5, 1E-2),
            'jet_thrust': (1E-6, 1E-2),
            'jet_puppiMultiplicity_charged': (1E-5, 1E-2),
            'jet_pTD_charged': (1E-6, 1E-2),
            'jet_LHA_charged': (1E-5, 1E-2),
            'jet_width_charged': (1E-6, 1E-3),
            'jet_thrust_charged': (1E-7, 1E-5),
        }

        zpj_region_dict = {
            "name": "ZPlusJets",
            "dirname": "ZPlusJets_QG_Unfold",
            "label": "Z+jets",
            "data_tfile": input_singlemu_tfile,
            "mc_tfile": input_mc_dy_mgpythia_tfile,
            "mc_label": "MG+Pythia8",
            "alt_mc_tfile": input_mc_dy_mgherwig_tfile,
            "alt_mc_label": "MG+Herwig++",
            "tau_limits": None,
            "backgrounds": [
                {
                    "name": "t#bar{t}",
                    "tfile": os.path.join(src_dir, qgc.TTBAR_FILENAME),
                    "rate_unc": 1.,
                },
                {
                    "name": "WW",
                    "tfile": os.path.join(src_dir, qgc.WW_FILENAME),
                    "rate_unc": 1.
                },
                {
                    "name": "WZ",
                    "tfile": os.path.join(src_dir, qgc.WZ_FILENAME),
                    "rate_unc": 1.
                },
                {
                    "name": "ZZ",
                    "tfile": os.path.join(src_dir, qgc.ZZ_FILENAME),
                    "rate_unc": 1.
                },
            ],
            "experimental_systematics": [
                {
                    "label": "Luminosity up",
                    "tfile": None,
                    "factor": 0.025,
                    "colour": ROOT.kCyan,
                },
                {
                    "label": "Luminosity down",
                    "tfile": None,
                    "factor": -0.025,
                    "colour": ROOT.kCyan,
                    'linestyle': 2,
                },
            ],
            "model_systematics": [
                {
                    "label": "muR up, muF nominal",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure,
                },
                {
                    "label": "muR down, muF nominal",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+1,
                },
                {
                    "label": "muR nominal, muF up",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+2,
                },
                {
                    "label": "muR nominal, muF down",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+3,
                },
                {
                    "label": "muR down, muF down",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+4,
                },
                {
                    "label": "muR up, muF up",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+5,
                },
            ],
            "pdf_systematics": [
                {
                    "label": "PDF",
                    "tfile": os.path.join(src_dir_systs, 'PDFvariationsTrue', qgc.DY_FILENAME),
                    "colour": ROOT.kCyan+2,
                    "variations": range(100),
                },
            ]
        }

        if args.doZPJ:
            zpj_region_dict['tau_limits'] = tau_limits
            regions.append(zpj_region_dict)

        if args.doZPJGroomed:
            zpj_region_groomed_dict = zpj_region_dict.copy()
            zpj_region_groomed_dict['tau_limits'] = tau_limits_groomed
            zpj_region_groomed_dict['dirname'] = 'ZPlusJets_QG_Unfold_groomed'
            zpj_region_groomed_dict['name'] = 'ZPlusJets_groomed'
            zpj_region_groomed_dict['label'] = 'Z+jets'
            regions.append(zpj_region_groomed_dict)

    # Setup various options
    # --------------------------------------------------------------------------

    REGULARIZE = args.regularize

    # Run with MC input instead of data
    MC_INPUT = args.MCinput
    mc_append = "_MC" if MC_INPUT else ""

    # If True, use part of MC for response matrix, and separate part for unfolding
    # as independent test
    # Should always be false for data
    if not args.MCinput:
        print("Ignoring your --MCsplit setting as running over data")

    MC_SPLIT = args.MCsplit if args.MCinput else False
    if MC_INPUT:
        mc_append += "_split" if MC_SPLIT else "_all"

    SUBTRACT_FAKES = True  # this should alwys be True
    sub_append = "_subFakes" if SUBTRACT_FAKES else ""

    append = ""

    if args.subtractBackgrounds:
        append += "_subBkg"

    if args.doExperimentalSysts:
        append += "_experimentalSyst"

    if args.doModelSysts:
        append += "_modelSystScale"

    if args.doPDFSysts:
        append += "_pdfSyst"

    if args.useAltResponse:
        append += "_altResponse"

    bias_str = ""
    if args.biasFactor != 0:
        bias_str = "_biasFactor%g" % args.biasFactor
        bias_str = bias_str.replace(".", "p")

    reg_axis_str = ""
    if REGULARIZE != "None":
        if args.regularizeAxis == 'pt':
            reg_axis_str = '_onlyRegPt'
        elif args.regularizeAxis == 'angle':
            reg_axis_str = '_onlyRegAngle'

    area_constraint = ROOT.TUnfold.kEConstraintArea
    area_constraint = ROOT.TUnfold.kEConstraintNone
    area_constraint_str = "Area" if area_constraint == ROOT.TUnfold.kEConstraintArea else "None"

    regularize_str = "regularize%s%s%s" % (str(REGULARIZE).capitalize(), bias_str, reg_axis_str)

    str_parts = dict(
        regularize_str=regularize_str,
        mc_append=mc_append,
        area=area_constraint_str,
        append=append,
        sub_append=sub_append,
    )
    output_dir = os.path.join(src_dir, "unfolding_{regularize_str}{mc_append}{sub_append}_densityModeBinWidth_constraint{area}{append}_signalRegionOnly_noHerwigPtReweight".format(**str_parts))

    if args.outputDir:
        output_dir = args.outputDir
    cu.check_dir_exists_create(output_dir)

    jet_algo = "AK4 PF PUPPI"
    if "ak8puppi" in src_dir:
        jet_algo = "AK8 PF PUPPI"

    angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]
    print(angles)

    print("Running TUnfold version", ROOT.TUnfold.GetTUnfoldVersion())

    LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT
    if 'target0p5' in src_dir:
        LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT_TARGET0p5
    elif 'target0p6' in src_dir:
        LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT_TARGET0p6


    # Do unfolding per signal region
    # --------------------------------------------------------------------------
    for region in regions[:]:

        # Setup pt bins
        # -------------
        # need different ones for Z+Jets region
        is_zpj = "ZPlusJets" in region['name']

        zpj_append = "_zpj" if is_zpj else ""

        pt_bin_edges_gen = qgc.PT_UNFOLD_DICT['signal%s_gen' % (zpj_append)]
        pt_bin_edges_reco = qgc.PT_UNFOLD_DICT['signal%s_reco' % (zpj_append)]

        pt_bin_edges_underflow_gen = qgc.PT_UNFOLD_DICT['underflow%s_gen' % (zpj_append)]
        pt_bin_edges_underflow_reco = qgc.PT_UNFOLD_DICT['underflow%s_reco' % (zpj_append)]

        # new_tdir = region['name']
        # output_tfile.mkdir(new_tdir)
        # region_tdir = output_tfile.Get(new_tdir)
        # region_tdir.cd()

        # Modify systematics as necessary
        # ----------------------------------------------------------------------

        # Remove the lumi one if we have no backgrounds, or the user has not said to remove backgrounds
        region['experimental_systematics'] = [syst_dict for syst_dict in region['experimental_systematics']
                                              if not ('lumi' in syst_dict['label'].lower()
                                                       and (len(region.get('backgrounds', [])) == 0
                                                            or not args.subtractBackgrounds))]

        # Do 1D unfolding of pt
        # ----------------------------------------------------------------------
        append = "%s_pt" % (region['name'])  # common str to put on filenames, etc
        # print("*"*80)
        # print("Region/var: %s" % (append))
        # print("*"*80)

        # hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_pt_reco_all" % (region['dirname']))
        mc_hname_append = "split" if MC_SPLIT else "all"
        if isinstance(region['mc_tfile'], str):
            region['mc_tfile'] = cu.open_root_file(region['mc_tfile'])
        # hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_pt_reco_%s" % (region['dirname'], mc_hname_append))
        # hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_pt_truth_%s" % (region['dirname'], mc_hname_append))
        hist_mc_gen_pt = cu.get_from_tfile(region['mc_tfile'], "%s/hist_pt_truth_%s" % (region['dirname'], mc_hname_append))
        # hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_pt_GenReco_%s" % (region['dirname'], mc_hname_append))
        # TODO!

        # Remake gen hist with physical bins & save to file
        all_pt_bins_gen = np.concatenate((pt_bin_edges_underflow_gen[:-1], pt_bin_edges_gen))
        hist_mc_gen_pt_physical = ROOT.TH1F("mc_gen_pt", ";p_{T}^{jet} [GeV];N", len(all_pt_bins_gen)-1, array('d', all_pt_bins_gen))
        update_hist_bin_content(h_orig=hist_mc_gen_pt, h_to_be_updated=hist_mc_gen_pt_physical)
        update_hist_bin_error(h_orig=hist_mc_gen_pt, h_to_be_updated=hist_mc_gen_pt_physical)

        # Do unfolding for each angle
        # ----------------------------------------------------------------------
        for angle in angles:
            angle_prepend = "groomed " if "groomed" in region['name'] else ""
            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
            this_angle_name = angle.name
            if (angle_prepend != ""
                and this_angle_name != 'LHA'
                and "_{T}" not in this_angle_name
                and "PUPPI" not in this_angle_name):
                # lower case if Groomed..., but be careful of e.g. pTD, LHA
                this_angle_name = this_angle_name[0].lower() + this_angle_name[1:]
            angle_str = "%s%s (%s)" % (angle_prepend, this_angle_name, angle.lambda_str)

            print("*"*120)
            print("Region/var: %s" % (append))
            print("*"*120)

            # put plots in subdir, to avoid overcrowding
            this_output_dir = "%s/%s/%s" % (output_dir, region['name'], angle.var)
            cu.check_dir_exists_create(this_output_dir)

            # Save hists etc to ROOT file for access later
            output_tfile = ROOT.TFile("%s/unfolding_result.root" % (this_output_dir), "RECREATE")

            new_tdir = "%s/%s" % (region['name'], angle.var)
            output_tfile.mkdir(new_tdir)
            this_tdir = output_tfile.Get(new_tdir)
            this_tdir.cd()
            this_tdir.WriteTObject(hist_mc_gen_pt_physical, "mc_gen_pt")


            # Setup MyUnfolder object to do unfolding etc
            # -------------------------------------------
            angle_bin_edges_reco = LAMBDA_VAR_DICTS[angle.var]['reco']
            angle_bin_edges_gen = LAMBDA_VAR_DICTS[angle.var]['gen']
            angle_shortname = angle.var.replace("jet_", "")

            if not isinstance(region['data_tfile'], ROOT.TFile):
                region['data_tfile'] = cu.open_root_file(region['data_tfile'])
            hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_all" % (region['dirname'], angle_shortname))
            mc_hname_append = "split" if MC_SPLIT else "all"
            hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_%s_GenReco_%s" % (region['dirname'], angle_shortname, mc_hname_append))

            # Need to scale if using H++ as input
            # hist_mc_gen_reco_map.Scale(1E8)
            # hist_mc_gen.Scale(1E8)
            # hist_mc_reco.Scale(1E8)

            # Actual distribution to be unfolded
            reco_1d = hist_mc_reco.Clone() if MC_INPUT else hist_data_reco

            hist_fakes_reco_fraction = None

            if SUBTRACT_FAKES:
                # to construct our "fakes" template, we use the ratio as predicted by MC, and apply it to data
                # this way we ensure we don't have -ve values, and avoid any issue with cross sections
                hist_mc_fakes_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                # hist_mc_fakes_reco.Scale(1E8)
                hist_fakes_reco = hist_mc_fakes_reco.Clone("hist_%s_fakes" % angle_shortname)
                hist_fakes_reco.Divide(hist_mc_reco)

                # plot fake fraction before multiplting by 'data'
                hist_fakes_reco_fraction = hist_fakes_reco.Clone("hist_fakes_reco_fraction")
                hist_fakes_reco.Multiply(reco_1d)

                # background-subtracted reco hists, only for plotting purposes, not for TUnfold (that does background subtraction internally)
                reco_1d_bg_subtracted = reco_1d.Clone()
                reco_1d_bg_subtracted.Add(hist_fakes_reco, -1)

                chosen_bin = 15
                print("1D reco input without background subtraction:", reco_1d.GetBinContent(chosen_bin))
                print("1D reco input with background subtraction:", reco_1d_bg_subtracted.GetBinContent(chosen_bin))
                print("1D reco input fakes:", hist_fakes_reco.GetBinContent(chosen_bin))

                hist_data_reco_bg_subtracted = hist_data_reco.Clone(hist_data_reco.GetName() + "_bgrSubtracted")
                hist_data_reco_bg_subtracted.Add(hist_fakes_reco, -1)

                hist_mc_reco_bg_subtracted = hist_mc_reco.Clone(hist_mc_reco.GetName() + "_bgrSubtracted")
                hist_mc_reco_bg_subtracted.Add(hist_mc_fakes_reco, -1)  # should this be hist_fakes_reco? depends on what we want to see...

            mc_hname_append = "_split" if MC_SPLIT else ""  # FIXME consistency in unfold hist module!
            hist_data_reco_gen_binning = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname'], angle_shortname))
            hist_mc_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_gen_binning%s" % (region['dirname'], angle_shortname, mc_hname_append))

            # hist_data_reco_gen_binning.Scale(1e8)
            # hist_mc_reco_gen_binning.Scale(1e8)

            # Actual distribution to be unfolded, but with gen binning
            reco_1d_gen_binning = hist_mc_reco_gen_binning.Clone() if MC_INPUT else hist_data_reco_gen_binning

            hist_fakes_reco_gen_binning = None
            if SUBTRACT_FAKES:
                mc_hname_append = "_split" if MC_SPLIT else ""  # FIXME consistency in unfold hist module!
                hist_mc_fakes_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_gen_binning%s" % (region['dirname'], angle_shortname, mc_hname_append))
                # create template as above, but with gen binning
                hist_fakes_reco_gen_binning = hist_mc_fakes_reco_gen_binning.Clone("hist_%s_fakes_gen_binning" % angle_shortname)
                hist_fakes_reco_gen_binning.Divide(hist_mc_reco_gen_binning)
                hist_fakes_reco_gen_binning.Multiply(reco_1d_gen_binning)

                # background-subtracted reco hists, only for plotting purposes, not for TUnfold (that does background subtraction internally)
                reco_1d_gen_binning_bg_subtracted = reco_1d_gen_binning.Clone()
                reco_1d_gen_binning_bg_subtracted.Add(hist_fakes_reco_gen_binning, -1)

                hist_data_reco_gen_binning_bg_subtracted = hist_data_reco_gen_binning.Clone(hist_data_reco_gen_binning.GetName() + "_bgrSubtracted")
                hist_data_reco_gen_binning_bg_subtracted.Add(hist_fakes_reco_gen_binning, -1)

                hist_mc_reco_gen_binning_bg_subtracted = hist_mc_reco_gen_binning.Clone(hist_mc_reco_gen_binning.GetName() + "_bgrSubtracted")
                hist_mc_reco_gen_binning_bg_subtracted.Add(hist_mc_fakes_reco_gen_binning, -1)  # should this be hist_fakes_reco_gen_binning? depends on what we want to see...

            # Setup unfolder object
            # ---------------------
            variable_name = "%s%s" % (angle_prepend, angle.name)
            axis_steering = '*[B]'
            if args.regularizeAxis == 'pt':
                axis_steering = 'pt[B];%s[N]' % variable_name
            elif args.regularizeAxis == 'angle':
                axis_steering = 'pt[N];%s[B]' % variable_name

            unfolder = MyUnfolder(response_map=hist_mc_gen_reco_map,
                                  variable_bin_edges_reco=angle_bin_edges_reco,
                                  variable_bin_edges_gen=angle_bin_edges_gen,
                                  variable_name=variable_name,
                                  pt_bin_edges_reco=pt_bin_edges_reco,
                                  pt_bin_edges_gen=pt_bin_edges_gen,
                                  pt_bin_edges_underflow_reco=pt_bin_edges_underflow_reco,
                                  pt_bin_edges_underflow_gen=pt_bin_edges_underflow_gen,
                                  orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                                  constraintMode=area_constraint,
                                  regMode=ROOT.TUnfold.kRegModeCurvature,
                                  densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth, # important as we have varying bin sizes!
                                  distribution='generatordistribution',  # the one to use for actual final regularisation/unfolding
                                  axisSteering=axis_steering)

            unfolder.save_binning(txt_filename="%s/binning_scheme.txt" % (this_output_dir), print_xml=False)

            unfolder_plotter = MyUnfolderPlotter(unfolder)
            plot_args = dict(output_dir=this_output_dir, append=append)

            # Set what is to be unfolded
            # ------------------------------------------------------------------
            unfolder.set_input(input_hist=reco_1d,
                               input_hist_gen_binning=reco_1d_gen_binning,
                               hist_truth=hist_mc_gen,
                               hist_mc_reco=hist_mc_reco,
                               hist_mc_reco_bg_subtracted=hist_mc_reco_bg_subtracted,
                               hist_mc_reco_gen_binning=hist_mc_reco_gen_binning,
                               hist_mc_reco_gen_binning_bg_subtracted=hist_mc_reco_gen_binning_bg_subtracted,
                               bias_factor=args.biasFactor)

            # Add systematic errors as different response matrices
            # ------------------------------------------------------------------
            if args.doExperimentalSysts:
                chosen_rsp_bin = (18, 18)
                print("nominal response bin content for", chosen_rsp_bin, hist_mc_gen_reco_map.GetBinContent(*chosen_rsp_bin))

                for syst_ind, syst_dict in enumerate(region['experimental_systematics']):
                    print("Adding systematic:", syst_dict['label'])
                    if 'factor' in syst_dict:
                        # special case for e.g. lumi - we construct a reponse hist, and add it using relative mode
                        rel_map = unfolder.response_map.Clone(syst_dict['label']+"Map")
                        for xbin, ybin in product(range(1, rel_map.GetNbinsX()+1), range(1, rel_map.GetNbinsY()+1)):
                            rel_map.SetBinContent(xbin, ybin, syst_dict['factor'])
                            rel_map.SetBinError(xbin, ybin, 0)
                        # unfolder.tunfolder.AddSysError(rel_map, syst_dict['label'], unfolder.orientation, ROOT.TUnfoldDensity.kSysErrModeRelative)
                        unfolder.add_sys_error(rel_map, syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeRelative)
                    else:
                        if not isinstance(syst_dict['tfile'], ROOT.TFile):
                            syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                        map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                        print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                        # unfolder.tunfolder.AddSysError(map_syst, syst_dict['label'], unfolder.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)
                        unfolder.add_sys_error(map_syst, syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)

            # Subtract fakes (treat as background)
            # ------------------------------------------------------------------
            if SUBTRACT_FAKES:
                unfolder.subtract_background(hist_fakes_reco, "Signal fakes", scale=1., scale_err=0.0)
                unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "Signal fakes", scale=1., scale_err=0.0)

            # Subtract actual backgrounds if necessary
            # ------------------------------------------------------------------
            background_reco_1d = None
            background_gen_1d = None
            if "backgrounds" in region and args.subtractBackgrounds:
                for bg_ind, bg_dict in enumerate(region['backgrounds']):
                    print("Subtracting", bg_dict['name'], 'background')
                    if not isinstance(bg_dict['tfile'], ROOT.TFile):
                        bg_dict['tfile'] = cu.open_root_file(bg_dict['tfile'])

                    mc_hname_append = "split" if MC_SPLIT else "all"
                    bg_hist = cu.get_from_tfile(bg_dict['tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    bg_hist_gen = cu.get_from_tfile(bg_dict['tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    bg_dict['hist'] = bg_hist
                    bg_dict['hist_gen'] = bg_hist

                    # keep one big hist
                    if not background_reco_1d:
                        background_reco_1d = bg_hist.Clone()
                        background_gen_1d = bg_hist_gen.Clone()
                    else:
                        background_reco_1d.Add(bg_hist, bg_dict.get('rate', 1.))
                        background_gen_1d.Add(bg_hist_gen, bg_dict.get('rate', 1.))

                    unfolder.subtract_background(hist=bg_hist,
                                                 name=bg_dict['name'],
                                                 scale=bg_dict.get('rate', 1.),
                                                 scale_err=bg_dict.get('rate_unc', 0.))

            title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_background_fractions(title=title, **plot_args)

            # Do any regularization
            # ---------------------
            unfolder.print_condition_number()

            tau = 0
            scan_mode = ROOT.TUnfoldDensity.kEScanTauRhoAvgSys
            scan_distribution = unfolder.distribution
            if REGULARIZE == "L":
                print("Regularizing with ScanLcurve, please be patient...")
                l_scanner = LCurveScanner()
                tau = l_scanner.scan_L(tunfolder=unfolder.tunfolder,
                                       n_scan=args.nScan,
                                       tau_min=region['tau_limits'][angle.var][0],
                                       tau_max=region['tau_limits'][angle.var][1])
                print("Found tau:", tau)
                l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (this_output_dir, append, OUTPUT_FMT))
                l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (this_output_dir, append, OUTPUT_FMT))
                l_scanner.save_to_tfile(this_tdir)

            elif REGULARIZE == "tau":
                print("Regularizing with ScanTau, please be patient...")
                tau_scanner = TauScanner()
                tau = tau_scanner.scan_tau(tunfolder=unfolder.tunfolder,
                                           n_scan=args.nScan,
                                           tau_min=region['tau_limits'][angle.var][0],
                                           tau_max=region['tau_limits'][angle.var][1],
                                           scan_mode=scan_mode,
                                           distribution=scan_distribution,
                                           axis_steering=unfolder.axisSteering)
                print("Found tau:", tau)
                tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (this_output_dir, append, OUTPUT_FMT))
                tau_scanner.save_to_tfile(this_tdir)

            if REGULARIZE != "None":
                title = "L matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_L_matrix(title=title, **plot_args)
                title = "L^{T}L matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_L_matrix_squared(title=title, **plot_args)
                title = "L * (x - bias vector), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_Lx_minus_bias(title=title, **plot_args)

            # Do unfolding!
            # ---------------------
            unfolder.do_unfolding(tau)
            unfolded_1d = unfolder.get_output(hist_name="unfolded_1d")
            chosen_bin = 18
            print("Bin %d:" % chosen_bin, unfolded_1d.GetBinContent(chosen_bin))
            print("original uncert:", unfolded_1d.GetBinError(chosen_bin))
            unfolder._post_process()

            # Get various error matrices
            # ------------------------------------------------------------------
            # stat errors only - do before or after systematics?
            print("stat uncert:", unfolder.unfolded_stat_err.GetBinError(chosen_bin))
            print("new uncert:", unfolder.unfolded.GetBinError(chosen_bin))

            # hist1, err1 = cu.th1_to_arr(unfolded_1d)
            # hist2, err2 = cu.th1_to_arr(hist_mc_gen)
            # cov, cov_err = cu.th2_to_arr(ematrix_total)
            # chi2 = calculate_chi2(hist1, hist2, cov)
            # print("my chi2 =", chi2)

            # Get shifts due to systematics
            # ------------------------------------------------------------------
            systematic_shift_hists = []
            if args.doExperimentalSysts:
                for syst_dict in region['experimental_systematics']:
                    h_syst = unfolder.get_delta_sys_shift(syst_label=syst_dict['label'])
                    systematic_shift_hists.append(h_syst)

            # Draw big 1D distributions
            # ------------------------------------------------------------------
            title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_unfolded_1d(is_data=not args.MCinput, output_dir=this_output_dir, append=append, title=title)

            # reco using detector binning
            unfolder_plotter.draw_detector_1d(do_reco_mc=True,
                                              do_reco_data=not MC_INPUT,
                                              output_dir=this_output_dir,
                                              append=append,
                                              title=title)

            # reco using gen binning
            unfolder_plotter.draw_generator_1d(do_reco_data=not MC_INPUT,
                                               do_reco_data_bg_sub=False,
                                               do_reco_bg=False,
                                               do_reco_mc=True,
                                               do_reco_mc_bg_sub=False,
                                               do_truth_mc=True,
                                               output_dir=this_output_dir,
                                               append=append,
                                               title=title)

            if SUBTRACT_FAKES:
                # same plot but with bg-subtracted reco (incl fakes)
                unfolder_plotter.draw_detector_1d(do_reco_data_bg_sub=not MC_INPUT,
                                                  do_reco_bg=SUBTRACT_FAKES,
                                                  do_reco_mc_bg_sub=True,
                                                  output_dir=this_output_dir,
                                                  append='bg_fakes_subtracted_%s' % append,
                                                  title=title)

                # same but with generator-binning
                unfolder_plotter.draw_generator_1d(do_reco_data=False,
                                                   do_reco_data_bg_sub=not MC_INPUT,
                                                   do_reco_bg=True,
                                                   do_reco_mc=False,
                                                   do_reco_mc_bg_sub=True,
                                                   do_truth_mc=True,
                                                   output_dir=this_output_dir,
                                                   append='bg_fakes_subtracted_%s' % append,
                                                   title=title)

            # Draw collapsed distributions
            # ---------------------
            # TODO!

            # Draw projections of response matrix vs 1D hist to check normalisation OK
            # Only makes sense if the same MC events go into matrix & 1D plot
            # ------------------------------------------------------------------
            if not MC_SPLIT:
                proj_reco = hist_mc_gen_reco_map.ProjectionY("proj_reco_%s" % (append))

                proj_gen = hist_mc_gen_reco_map.ProjectionX("proj_gen_%s" % (append))
                draw_projection_comparison(unfolder.hist_truth, proj_gen,
                                           title="%s\n%s region" % (jet_algo, region['label']),
                                           xtitle="%s, Generator binning" % (angle_str),
                                           output_filename="%s/projection_gen_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

                print("projection reco #bins:", proj_reco.GetNbinsX())
                print("response map # bins x:", hist_mc_gen_reco_map.GetNbinsX())
                print("response map # bins y:", hist_mc_gen_reco_map.GetNbinsY())
                if SUBTRACT_FAKES:
                    print("reco bg subtracted #bins:", hist_mc_reco_bg_subtracted.GetNbinsX())
                    print(proj_reco.GetNbinsX())
                    # Do the same but with backgrounds subtracted from the 1D
                    draw_projection_comparison(hist_mc_reco_bg_subtracted, proj_reco,
                                               title="%s\n%s region" % (jet_algo, region['label']),
                                               xtitle="%s, Detector binning" % (angle_str),
                                               output_filename="%s/projection_reco_bg_subtracted_%s.%s" % (this_output_dir, append, OUTPUT_FMT))
                else:
                    draw_projection_comparison(hist_mc_reco, proj_reco,
                           title="%s\n%s region" % (jet_algo, region['label']),
                           xtitle="%s, Detector binning" % (angle_str),
                           output_filename="%s/projection_reco_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                           print_bin_comparison=False)

            # Draw matrices
            # ------------------------------------------------------------------
            unfolder_plotter.plot_bias_vector(title=title, **plot_args)

            title = "Response matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_response_matrix(title=title, **plot_args)

            # title = "Covariance matrix, %s region, %s" % (region['label'], angle_str)
            # unfolder_plotter.draw_covariance_matrix(title=title, **plot_args)

            title = ("#splitline{Probability matrix, %s, %s region, %s}{Condition number: #sigma_{max} / #sigma_{min} = %.3g / %.3g = %g}"
                        % (jet_algo, region['label'], angle_str, unfolder.sigma_max, unfolder.sigma_min, unfolder.condition_number))
            unfolder_plotter.draw_probability_matrix(title=title, **plot_args)

            title = "Correlation matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_correlation_matrix(title=title, draw_values=True, **plot_args)
            unfolder_plotter.draw_correlation_matrix(title=title, draw_values=False, **plot_args)

            title = "Error matrix (input), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_input(title=title, **plot_args)

            title = "Error matrix (sys uncorr), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_sys_uncorr(title=title, **plot_args)

            title = "Error matrix (total), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_total(title=title, **plot_args)

            # Do forward-folding to check unfolding
            # ------------------------------------------------------------------
            # Do it on the unfolded result
            title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_unfolded_folded(title=title, **plot_args)

            if MC_INPUT:
                # Folded MC truth
                title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_truth_folded(title=title, **plot_args)

            # Save everything to TFile
            unfolder.save_to_tfile(this_tdir)

            # ------------------------------------------------------------------
            # UNFOLDING WITH ALTERNATIVE RESPONSE MATRIX
            # ------------------------------------------------------------------
            alt_hist_mc_gen = None  # mc truth of the generator used to make reponse matrix
            alt_unfolder = None
            if args.useAltResponse:
                print("*" * 80)
                print("*** Unfolding with alternate response matrix ***")
                print("*" * 80)

                if not isinstance(region['alt_mc_tfile'], ROOT.TFile):
                    region['alt_mc_tfile'] = cu.open_root_file(region['alt_mc_tfile'])
                alt_hist_mc_gen = cu.get_from_tfile(region['alt_mc_tfile'], "%s/hist_%s_truth_all" % (region['dirname'], angle_shortname))
                hist_mc_gen_reco_map_alt = cu.get_from_tfile(region['alt_mc_tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                alt_unfolder = MyUnfolder(response_map=hist_mc_gen_reco_map_alt,
                                          variable_bin_edges_reco=unfolder.variable_bin_edges_reco,
                                          variable_bin_edges_gen=unfolder.variable_bin_edges_gen,
                                          variable_name=unfolder.variable_name,
                                          pt_bin_edges_reco=unfolder.pt_bin_edges_reco,
                                          pt_bin_edges_gen=unfolder.pt_bin_edges_gen,
                                          pt_bin_edges_underflow_reco=unfolder.pt_bin_edges_underflow_reco,
                                          pt_bin_edges_underflow_gen=unfolder.pt_bin_edges_underflow_gen,
                                          orientation=unfolder.orientation,
                                          constraintMode=unfolder.constraintMode,
                                          regMode=unfolder.regMode,
                                          densityFlags=unfolder.densityFlags,
                                          distribution=unfolder.distribution,
                                          axisSteering=unfolder.axisSteering)

                this_tdir.cd()
                alt_tdir = this_tdir.mkdir("alt_response_%s" % region['alt_mc_label'].replace(" ", "-"))
                alt_tdir.cd()

                alt_unfolder_plotter = MyUnfolderPlotter(alt_unfolder)
                alt_output_dir = this_output_dir+"/altResponse"
                alt_plot_args = dict(output_dir=alt_output_dir,
                                     append=append)

                # Only plot response matrix
                # --------------------------------------------------------------
                title = ("#splitline{Probability matrix, %s region, %s, %s}{Condition number: #sigma_{max} / #sigma_{min} = %.3g / %.3g = %g}"
                            % (region['label'], angle_str, region['alt_mc_label'], unfolder.sigma_max, unfolder.sigma_min, unfolder.condition_number))
                alt_unfolder_plotter.draw_probability_matrix(title=title, **alt_plot_args)

                # Set what is to be unfolded - same as main unfolder
                # --------------------------------------------------------------
                alt_unfolder.set_input(input_hist=reco_1d,
                                       input_hist_gen_binning=reco_1d_gen_binning,
                                       hist_truth=unfolder.hist_truth.Clone(),
                                       hist_mc_reco=unfolder.hist_mc_reco.Clone(),
                                       hist_mc_reco_bg_subtracted=unfolder.hist_mc_reco_bg_subtracted.Clone(),
                                       hist_mc_reco_gen_binning=unfolder.hist_mc_reco_gen_binning.Clone(),
                                       hist_mc_reco_gen_binning_bg_subtracted=unfolder.hist_mc_reco_gen_binning_bg_subtracted.Clone(),
                                       bias_factor=args.biasFactor)

                # Subtract fakes (treat as background)
                # --------------------------------------------------------------
                if SUBTRACT_FAKES:
                    alt_unfolder.subtract_background(hist_fakes_reco, "fakes")

                # Do any regularization
                # --------------------------------------------------------------
                alt_tau = 0
                if REGULARIZE == "L":
                    print("Regularizing alternative with ScanL, please be patient...")
                    alt_L_scanner = LCurveScanner()
                    alt_tau = alt_l_scanner.scan_L(tunfolder=alt_unfolder.tunfolder,
                                               n_scan=args.nScan,
                                               tau_min=region['tau_limits'][angle.var][0],
                                               tau_max=region['tau_limits'][angle.var][1])
                    print("Found tau:", alt_tau)
                    alt_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_alt_%s.%s" % (alt_output_dir, unfolder.variable_name, OUTPUT_FMT))
                    alt_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_alt_%s.%s" % (alt_output_dir, unfolder.variable_name, OUTPUT_FMT))
                    alt_l_scanner.save_to_tfile(alt_tdir)

                elif REGULARIZE == "tau":
                    print("Regularizing alternative with ScanTau, please be patient...")
                    alt_tau_scanner = TauScanner()
                    alt_tau = alt_tau_scanner.scan_tau(tunfolder=alt_unfolder.tunfolder,
                                                       n_scan=args.nScan,
                                                       tau_min=region['tau_limits'][angle.var][0],
                                                       tau_max=region['tau_limits'][angle.var][1],
                                                       scan_mode=scan_mode,
                                                       distribution=scan_distribution,
                                                       axis_steering=alt_unfolder.axisSteering)
                    print("Found tau for alt matrix:", alt_tau)
                    alt_tau_scanner.plot_scan_tau(output_filename="%s/scantau_alt_%s.%s" % (alt_output_dir, alt_unfolder.variable_name, OUTPUT_FMT))
                    alt_tau_scanner.save_to_tfile(alt_tdir)

                if REGULARIZE != "None":
                    title = "L matrix, %s region, %s, alt. response (%s)" % (region['label'], angle_str, region['alt_mc_label'])
                    alt_unfolder_plotter.draw_L_matrix(title=title, **alt_plot_args)
                    title = "L^{T}L matrix, %s region, %s, alt. response (%s)" % (region['label'], angle_str, region['alt_mc_label'])
                    alt_unfolder_plotter.draw_L_matrix_squared(title=title, **alt_plot_args)
                    title = "L * (x - bias vector), %s region, %s,  alt. response (%s)" % (region['label'], angle_str, region['alt_mc_label'])
                    alt_unfolder_plotter.draw_Lx_minus_bias(title=title, **alt_plot_args)

                # Do unfolding!
                # --------------------------------------------------------------
                alt_unfolder.do_unfolding(alt_tau)
                alt_unfolded_1d = alt_unfolder.get_output(hist_name="alt_unfolded_1d")
                print("Bin %d:" % chosen_bin, alt_unfolded_1d.GetBinContent(chosen_bin))
                print("original uncert:", alt_unfolded_1d.GetBinError(chosen_bin))
                alt_unfolder._post_process()

                alt_title = "%s\n%s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_unfolded_1d(is_data=not args.MCinput, title=alt_title, **alt_plot_args)

                # Save important stuff to TFile
                # --------------------------------------------------------------
                alt_unfolder.save_to_tfile(alt_tdir)

            # ------------------------------------------------------------------
            # MODEL INPUT VARIATIONS
            # ------------------------------------------------------------------
            # For each model variation, we unfold using the same settings as
            # the nominal one, just changing the input 1D hist
            if args.doModelSysts:
                for ind, syst_dict in enumerate(region['model_systematics']):
                    syst_label = syst_dict['label']
                    syst_label_no_spaces = syst_dict['label'].replace(", ", "_").replace(" ", "_").replace("{", "").replace("}", "")

                    print("*" * 80)
                    print("*** Unfolding with alternate input:", syst_label, "(%d/%d) ***" % (ind+1, len(region['model_systematics'])))
                    print("*" * 80)

                    is_herwig = "Herwig" in syst_label

                    mc_hname_append = "split" if MC_SPLIT else "all"
                    if is_herwig:
                        # use all the stats!
                        mc_hname_append = "all"
                    if not isinstance(syst_dict['tfile'], ROOT.TFile):
                        syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                    hist_syst_reco = cu.get_from_tfile(syst_dict['tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    hist_syst_gen = cu.get_from_tfile(syst_dict['tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))

                    syst_unfolder = MyUnfolder(response_map=unfolder.response_map,
                                               variable_bin_edges_reco=unfolder.variable_bin_edges_reco,
                                               variable_bin_edges_gen=unfolder.variable_bin_edges_gen,
                                               variable_name=unfolder.variable_name,
                                               pt_bin_edges_reco=unfolder.pt_bin_edges_reco,
                                               pt_bin_edges_gen=unfolder.pt_bin_edges_gen,
                                               pt_bin_edges_underflow_reco=unfolder.pt_bin_edges_underflow_reco,
                                               pt_bin_edges_underflow_gen=unfolder.pt_bin_edges_underflow_gen,
                                               orientation=unfolder.orientation,
                                               constraintMode=unfolder.constraintMode,
                                               regMode=unfolder.regMode,
                                               densityFlags=unfolder.densityFlags,
                                               distribution=unfolder.distribution,
                                               axisSteering=unfolder.axisSteering)

                    this_tdir.cd()
                    syst_tdir_name = "modelSyst_"+syst_label_no_spaces
                    syst_tdir = this_tdir.mkdir(syst_tdir_name)
                    syst_tdir.cd()

                    syst_output_dir = this_output_dir+"/modelSyst_"+syst_label_no_spaces
                    syst_unfolder_plotter = MyUnfolderPlotter(syst_unfolder)
                    syst_plot_args = dict(output_dir=syst_output_dir,
                                          append=append)

                    # if is_herwig:
                        # herwig_sf = 1E6  # should've done this earlier
                        # hist_syst_reco.Scale(herwig_sf)
                        # hist_syst_gen.Scale(herwig_sf)
                        # SetEpsMatrix ensures rank properly calculated when inverting
                        # "rank of matrix E 55 expect 170"
                        # syst_unfolder.tunfolder.SetEpsMatrix(1E-18)

                    # because we only care about shape, not overall normalisation
                    # (which can artificially increase/decrease errors)
                    # we normalise to the nominal integral
                    # Note that we use the scaling from gen level, to take
                    # into account any reco-dependent efficiencies
                    # TODO: is this right?
                    sf = hist_mc_gen.Integral() / hist_syst_gen.Integral()
                    hist_syst_reco.Scale(sf)
                    hist_syst_gen.Scale(sf)

                    if SUBTRACT_FAKES:
                        # Use the background template from the nominal MC
                        # (since we're only testing different input shapes,
                        # and our bkg estimate is always from MC)
                        hist_fakes_syst = hist_fakes_reco_fraction.Clone("hist_fakes_syst_%s" % syst_label_no_spaces)
                        hist_fakes_syst.Multiply(hist_syst_reco)

                    hist_mc_reco_bg_subtracted = hist_syst_reco.Clone()
                    hist_mc_reco_bg_subtracted.Add(hist_fakes_syst, -1)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    syst_unfolder.set_input(input_hist=hist_syst_reco,
                                            hist_truth=hist_syst_gen,
                                            hist_mc_reco=hist_syst_reco,
                                            hist_mc_reco_bg_subtracted=hist_mc_reco_bg_subtracted,
                                            bias_factor=args.biasFactor)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    if SUBTRACT_FAKES:
                        syst_unfolder.subtract_background(hist_fakes_syst, "fakes")

                    syst_unfolder_plotter.draw_detector_1d(do_reco_data=False,
                                                           do_reco_data_bg_sub=False,
                                                           do_reco_bg=True,
                                                           do_reco_mc=False,
                                                           do_reco_mc_bg_sub=True,
                                                           output_dir=syst_plot_args['output_dir'],
                                                           append='bg_fakes_subtracted_%s' % append,
                                                           title="%s region, %s, %s" % (region['label'], angle_str, syst_label))

                    # Add systematic errors as different response matrices
                    # ----------------------------------------------------
                    if args.doExperimentalSysts:
                        chosen_rsp_bin = (18, 18)
                        print("nominal response bin content for", chosen_rsp_bin, syst_unfolder.response_map.GetBinContent(*chosen_rsp_bin))
                        for exp_dict in region['experimental_systematics']:
                            if not isinstance(exp_dict['tfile'], ROOT.TFile):
                                exp_dict['tfile'] = cu.open_root_file(exp_dict['tfile'])
                            map_syst = cu.get_from_tfile(exp_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                            print("Adding systematic:", exp_dict['label'])
                            print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                            # syst_unfolder.tunfolder.AddSysError(map_syst, exp_dict['label'], syst_unfolder.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)
                            syst_unfolder.add_sys_error(map_syst, exp_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)


                    # Do any regularization
                    # --------------------------------------------------------------
                    syst_tau = 0
                    if REGULARIZE == "L":
                        print("Regularizing systematic model with ScanL, please be patient...")
                        syst_l_scanner = LCurveScanner()
                        syst_tau = syst_l_scanner.scan_L(tunfolder=syst_unfolder.tunfolder,
                                                         n_scan=args.nScan,
                                                         tau_min=region['tau_limits'][angle.var][0],
                                                         tau_max=region['tau_limits'][angle.var][1])
                        print("Found tau:", syst_tau)
                        syst_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_syst_%s_%s.%s" % (syst_output_dir, syst_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                        syst_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_syst_%s_%s.%s" % (syst_output_dir, syst_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                        syst_l_scanner.save_to_tfile(syst_tdir)

                    elif REGULARIZE == "tau":
                        print("Regularizing systematic model with ScanTau, please be patient...")
                        syst_tau_scanner = TauScanner()
                        syst_tau = syst_tau_scanner.scan_tau(tunfolder=syst_unfolder.tunfolder,
                                                             n_scan=args.nScan,
                                                             tau_min=region['tau_limits'][angle.var][0],
                                                             tau_max=region['tau_limits'][angle.var][1],
                                                             scan_mode=scan_mode,
                                                             distribution=scan_distribution,
                                                             axis_steering=syst_unfolder.axisSteering)
                        print("Found tau for syst matrix:", syst_tau)
                        syst_tau_scanner.plot_scan_tau(output_filename="%s/scantau_syst_%s_%s.%s" % (syst_output_dir, syst_label_no_spaces, syst_unfolder.variable_name, OUTPUT_FMT))
                        syst_tau_scanner.save_to_tfile(syst_tdir)

                    region['model_systematics'][ind]['tau'] = syst_tau

                    if REGULARIZE != "None":
                        title = "L matrix, %s region, %s,\n%s" % (region['label'], angle_str, syst_label)
                        syst_unfolder_plotter.draw_L_matrix(title=title, **syst_plot_args)
                        title = "L^{T}L matrix, %s region, %s,\n%s" % (region['label'], angle_str, syst_label)
                        syst_unfolder_plotter.draw_L_matrix_squared(title=title, **syst_plot_args)
                        title = "L * (x - bias vector), %s region, %s,\n%s" % (region['label'], angle_str, syst_label)
                        syst_unfolder_plotter.draw_Lx_minus_bias(title=title, **syst_plot_args)

                    # Do unfolding!
                    # --------------------------------------------------------------
                    syst_unfolder.do_unfolding(syst_tau)
                    syst_unfolded_1d = syst_unfolder.get_output(hist_name="syst_%s_unfolded_1d" % (syst_label_no_spaces))
                    print("Bin %d:" % (chosen_bin), syst_unfolded_1d.GetBinContent(chosen_bin))
                    print("original uncert:", syst_unfolded_1d.GetBinError(chosen_bin))
                    syst_unfolder._post_process()

                    region['model_systematics'][ind]['unfolded_1d'] = syst_unfolder.unfolded
                    region['model_systematics'][ind]['gen_1d'] = syst_unfolder.hist_truth

                    syst_title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, syst_label)
                    syst_unfolder_plotter.draw_unfolded_1d(is_data=not args.MCinput, title=syst_title, **syst_plot_args)

                    # Save important stuff to TFile
                    # --------------------------------------------------------------
                    syst_unfolder.save_to_tfile(syst_tdir)

            # ------------------------------------------------------------------
            # DO PDF VARIATIONS
            # ------------------------------------------------------------------
            if args.doPDFSysts:
                # first construct all new systematic variations
                original_pdf_dict = region['pdf_systematics'][0]
                original_pdf_dict['label'] = '_PDF_template'  # initial _ to ignore it later on
                tfile = original_pdf_dict['tfile']
                if not isinstance(tfile, ROOT.TFile):
                    tfile = cu.open_root_file(tfile)

                region['pdf_systematics']  = []
                for pdf_ind in original_pdf_dict['variations']:
                    mc_hname_append = "split" if MC_SPLIT else "all"
                    mc_hname_append = "all"
                    region['pdf_systematics'].append(
                        {
                            "label": "PDF_%d" % (pdf_ind),
                            "hist_reco": cu.get_from_tfile(tfile, "%s/hist_%s_reco_%s_PDF_%d" % (region['dirname'], angle_shortname, mc_hname_append, pdf_ind)),
                            "hist_gen": cu.get_from_tfile(tfile, "%s/hist_%s_gen_%s_PDF_%d" % (region['dirname'], angle_shortname, mc_hname_append, pdf_ind)),
                            "colour": ROOT.kCyan+2,
                        })


                # Now run over all variations like for model systs
                for ind, pdf_dict in enumerate(region['pdf_systematics']):
                    pdf_label = pdf_dict['label']
                    pdf_label_no_spaces = pdf_label.replace(" ", "_")

                    if pdf_label.startswith("_"):
                        continue

                    print("*" * 80)
                    print("*** Unfolding with alternate input:", pdf_label, "(%d/%d) ***" % (ind+1, len(region['pdf_systematics'])))
                    print("*" * 80)

                    hist_pdf_reco = pdf_dict['hist_reco']
                    hist_pdf_gen = pdf_dict['hist_gen']

                    pdf_unfolder = MyUnfolder(response_map=unfolder.response_map,
                                              variable_bin_edges_reco=unfolder.variable_bin_edges_reco,
                                              variable_bin_edges_gen=unfolder.variable_bin_edges_gen,
                                              variable_name=unfolder.variable_name,
                                              pt_bin_edges_reco=unfolder.pt_bin_edges_reco,
                                              pt_bin_edges_gen=unfolder.pt_bin_edges_gen,
                                              pt_bin_edges_underflow_reco=unfolder.pt_bin_edges_underflow_reco,
                                              pt_bin_edges_underflow_gen=unfolder.pt_bin_edges_underflow_gen,
                                              orientation=unfolder.orientation,
                                              constraintMode=unfolder.constraintMode,
                                              regMode=unfolder.regMode,
                                              densityFlags=unfolder.densityFlags,
                                              distribution=unfolder.distribution,
                                              axisSteering=unfolder.axisSteering)

                    this_tdir.cd()
                    pdf_tdir_name = "pdfSyst_"+pdf_label_no_spaces
                    pdf_tdir = this_tdir.mkdir(pdf_tdir_name)
                    pdf_tdir.cd()

                    pdf_unfolder_plotter = MyUnfolderPlotter(pdf_unfolder)
                    pdf_output_dir = this_output_dir+"/pdfSyst/"+pdf_label_no_spaces,
                    pdf_plot_args = dict(output_dir=pdf_output_dir,
                                         append=append)

                    if SUBTRACT_FAKES:
                        # Use the background template from the nominal MC
                        # (since we're only testing different input shapes,
                        # and our bkg estimate is always from MC)
                        hist_fakes_pdf = hist_fakes_reco_fraction.Clone("hist_fakes_pdf_%s" % pdf_label_no_spaces)
                        hist_fakes_pdf.Multiply(hist_pdf_reco)

                    hist_pdf_reco_bg_subtracted = hist_pdf_reco.Clone()
                    hist_pdf_reco_bg_subtracted.Add(hist_fakes_pdf, -1)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    pdf_unfolder.set_input(input_hist=hist_pdf_reco,
                                           hist_truth=hist_pdf_gen,
                                           hist_mc_reco=hist_pdf_reco,
                                           hist_mc_reco_bg_subtracted=hist_pdf_reco_bg_subtracted,
                                           bias_factor=args.biasFactor)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    if SUBTRACT_FAKES:
                        pdf_unfolder.subtract_background(hist_fakes_pdf, "fakes")

                    pdf_unfolder_plotter.draw_detector_1d(do_reco_data=False,
                                                          do_reco_data_bg_sub=False,
                                                          do_reco_bg=True,
                                                          do_reco_mc=False,
                                                          do_reco_mc_bg_sub=True,
                                                          output_dir=pdf_plot_args['output_dir'],
                                                          append='bg_fakes_subtracted_%s' % append,
                                                          title="%s region, %s, %s" % (region['label'], angle_str, pdf_label))

                    # Add systematic errors as different response matrices
                    # ----------------------------------------------------
                    if args.doExperimentalSysts:
                        chosen_rsp_bin = (18, 18)
                        print("nominal response bin content for", chosen_rsp_bin, pdf_unfolder.response_map.GetBinContent(*chosen_rsp_bin))
                        for exp_dict in region['experimental_systematics']:
                            print("Adding systematic:", exp_dict['label'])
                            if 'factor' in exp_dict:
                                # special case for e.g. lumi - we construct a reponse hist, and add it using relative mode
                                rel_map = pdf_unfolder.response_map.Clone(exp_dict['label']+"MapPDF")
                                for xbin, ybin in product(range(1, rel_map.GetNbinsX()+1), range(1, rel_map.GetNbinsY()+1)):
                                    rel_map.SetBinContent(xbin, ybin, exp_dict['factor'])
                                    rel_map.SetBinError(xbin, ybin, 0)
                                # pdf_unfolder.tunfolder.AddSysError(rel_map, exp_dict['label'], ROOT.TUnfold.kHistMapOutputHoriz, ROOT.TUnfoldDensity.kSysErrModeRelative)
                                pdf_unfolder.add_sys_error(rel_map, exp_dict['label'], ROOT.TUnfoldDensity.kSysErrModeRelative)
                            else:
                                if not isinstance(exp_dict['tfile'], ROOT.TFile):
                                    exp_dict['tfile'] = cu.open_root_file(exp_dict['tfile'])
                                map_syst = cu.get_from_tfile(exp_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                                print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                                # pdf_unfolder.tunfolder.AddSysError(map_syst, exp_dict['label'], pdf_unfolder.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)
                                pdf_unfolder.add_sys_error(map_syst, exp_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)

                    # Do any regularization
                    # --------------------------------------------------------------
                    syst_tau = 0
                    if REGULARIZE == "L":
                        print("Regularizing systematic model with ScanL, please be patient...")
                        syst_l_scanner = LCurveScanner()
                        syst_tau = syst_l_scanner.scan_L(tunfolder=pdf_unfolder.tunfolder,
                                                         n_scan=args.nScan,
                                                         tau_min=region['tau_limits'][angle.var][0],
                                                         tau_max=region['tau_limits'][angle.var][1])
                        print("Found tau:", syst_tau)
                        syst_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_syst_%s_%s.%s" % (pdf_output_dir, pdf_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                        syst_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_syst_%s_%s.%s" % (pdf_output_dir, pdf_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                        syst_l_scanner.save_to_tfile(pdf_tdir)

                    elif REGULARIZE == "tau":
                        print("Regularizing systematic model with ScanTau, please be patient...")
                        syst_tau_scanner = TauScanner()
                        syst_tau = syst_tau_scanner.scan_tau(tunfolder=pdf_unfolder.tunfolder,
                                                             n_scan=args.nScan,
                                                             tau_min=region['tau_limits'][angle.var][0],
                                                             tau_max=region['tau_limits'][angle.var][1],
                                                             scan_mode=scan_mode,
                                                             distribution=scan_distribution,
                                                             axis_steering=pdf_unfolder.axisSteering)
                        print("Found tau for syst matrix:", syst_tau)
                        syst_tau_scanner.plot_scan_tau(output_filename="%s/scantau_syst_%s_%s.%s" % (pdf_output_dir, pdf_label_no_spaces, pdf_unfolder.variable_name, OUTPUT_FMT))
                        syst_tau_scanner.save_to_tfile(pdf_tdir)

                    region['pdf_systematics'][ind]['tau'] = syst_tau

                    # Do unfolding!
                    # --------------------------------------------------------------
                    pdf_unfolder.do_unfolding(syst_tau)
                    pdf_unfolded_1d = pdf_unfolder.get_output(hist_name="syst_%s_unfolded_1d" % (pdf_label_no_spaces))
                    print("Bin %d:" % (chosen_bin), pdf_unfolded_1d.GetBinContent(chosen_bin))
                    print("original uncert:", pdf_unfolded_1d.GetBinError(chosen_bin))
                    pdf_unfolder._post_process()

                    region['pdf_systematics'][ind]['unfolded_1d'] = pdf_unfolder.unfolded
                    region['pdf_systematics'][ind]['gen_1d'] = pdf_unfolder.hist_truth

                    pdf_title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, pdf_label)
                    pdf_unfolder_plotter.draw_unfolded_1d(is_data=not args.MCinput, title=pdf_title, **pdf_plot_args)

                    # Save important stuff to TFile
                    # ----------------------------------------------------------
                    pdf_unfolder.save_to_tfile(pdf_tdir)

            # ------------------------------------------------------------------
            # PLOTTING LOTS OF THINGS
            # ------------------------------------------------------------------

            # Some common plotting vars
            # ------------------------------------------------------------------
            detector_title = "Detector-level " + angle_str
            particle_title = "Particle-level " + angle_str
            normalised_differential_label = "#frac{1}{#sigma} #frac{d^{2}#sigma}{dp_{T} d%s}" % angle.lambda_str
            summary_1d_entries = []  # for final summary plot

            # Draw individual pt bin plots - GEN binning
            # ------------------------------------------------------------------
            for ibin_pt in range(0, len(unfolder.pt_bin_edges_gen)-1):

                this_pt_bin_tdir = this_tdir.mkdir("gen_bin_%d" % (ibin_pt))
                print("Individual gen bin", ibin_pt, "=", unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1])

                # Produce 1D hists for this pt bin
                # --------------------------------------------------------------
                # Unfolded hists
                mc_gen_hist_bin = unfolder.get_var_hist_pt_binned(unfolder.hist_truth, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(mc_gen_hist_bin, "mc_gen_hist_bin")

                unfolded_hist_bin_stat_errors = unfolder.get_var_hist_pt_binned(unfolder.unfolded_stat_err, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin_stat_errors, "unfolded_hist_bin_stat_errors")

                unfolded_hist_bin_total_errors = unfolder.get_var_hist_pt_binned(unfolder.unfolded, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin_total_errors, "unfolded_hist_bin_total_errors") # same as unfolded_1d?

                # Reco level but with gen binning
                # For 'data'
                reco_hist_bin_gen_binning = unfolder.get_var_hist_pt_binned(unfolder.input_hist_gen_binning, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(reco_hist_bin_gen_binning, "reco_hist_bin_gen_binning")

                if SUBTRACT_FAKES:
                    reco_hist_bg_subtracted_bin_gen_binning = unfolder.get_var_hist_pt_binned(unfolder.input_hist_gen_binning_bg_subtracted, ibin_pt, binning_scheme="generator")
                    this_pt_bin_tdir.WriteTObject(reco_hist_bg_subtracted_bin_gen_binning, "reco_hist_bg_subtracted_bin_gen_binning")

                # For MC
                mc_reco_hist_bin_gen_binning = unfolder.get_var_hist_pt_binned(unfolder.hist_mc_reco_gen_binning, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(mc_reco_hist_bin_gen_binning, "mc_reco_hist_bin_gen_binning")

                if SUBTRACT_FAKES:
                    mc_reco_hist_bg_subtracted_bin_gen_binning = unfolder.get_var_hist_pt_binned(unfolder.hist_mc_reco_gen_binning_bg_subtracted, ibin_pt, binning_scheme="generator")
                    this_pt_bin_tdir.WriteTObject(mc_reco_hist_bg_subtracted_bin_gen_binning, "mc_reco_hist_bg_subtracted_bin_gen_binning")


                # Make lots of plots
                # ------------------------------------------------------------
                lw = 2
                # common hist settings
                title = "%s\n%s region\n%g < p_{T}^{jet} < %g GeV" % (jet_algo, region['label'], unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1])
                if "ptavebinning" in src_dir.lower():
                    title = "%s\n%s region\n%g < #LT p_{T}^{jet} #GT < %g GeV" % (jet_algo, region['label'], unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1])
                common_hist_args = dict(
                    what="hist",
                    title=title,
                    subplot_type='ratio',
                    # subplot_limits=(0.5, 1.5),
                    subplot_limits=(0.75, 1.25),
                )
                subplot_title = "Unfolded / Gen"

                # PLOT UNFOLDED DATA
                # --------------------------------------------------------------
                gen_colour = ROOT.kRed
                unfolded_basic_colour = ROOT.kAzure+7
                unfolded_stat_colour = ROOT.kAzure+7
                unfolded_total_colour = ROOT.kBlack

                def _modify_plot(this_plot):
                    this_plot.legend.SetX1(0.6)
                    this_plot.legend.SetY1(0.68)
                    this_plot.legend.SetX2(0.98)
                    this_plot.legend.SetY2(0.9)
                    this_plot.left_margin = 0.16

                # unnormalised version
                entries = [
                    Contribution(mc_gen_hist_bin,
                                 label="Generator",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_total_errors,
                                 label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_gen_hist_bin,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_stat_errors,
                                 label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 subplot=mc_gen_hist_bin,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=particle_title,
                            ytitle="N",
                            subplot_title='Unfolded / gen',
                            **common_hist_args)
                _modify_plot(plot)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_unnormalised_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # now normalise each plot to unity
                # Note that this modifies e.g. mc_gen_hist_bin, so from this point
                # onwards it will be normalised to unity
                entries = [
                    Contribution(mc_gen_hist_bin,
                                 label="Generator (MG+Pythia8)",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin_total_errors,
                                 label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_gen_hist_bin,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin_stat_errors,
                                 label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=particle_title,
                            ytitle=normalised_differential_label,
                            subplot_title=subplot_title,
                            **common_hist_args)
                _modify_plot(plot)
                plot.plot("NOSTACK E1")
                # plot.save("%s/unfolded_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Do a version where divided by bin width
                # Note that these hists are already normalised to 1!
                # Do not use normalise_hist!
                mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(mc_gen_hist_bin)
                unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)
                entries = [
                    Contribution(mc_gen_hist_bin_div_bin_width,
                                 label="Generator (MG+Pythia8)",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_total_errors_div_bin_width,
                                 label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_gen_hist_bin_div_bin_width,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                 label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=particle_title,
                            ytitle=normalised_differential_label,
                            subplot_title=subplot_title,
                            **common_hist_args)
                _modify_plot(plot)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                summary_1d_entries.append([
                    (mc_gen_hist_bin_div_bin_width,
                        dict(label=region['mc_label'],
                             line_color=gen_colour, line_width=lw,
                             marker_color=gen_colour, marker_size=0,
                             normalise_hist=True)),  # generator
                    (unfolded_hist_bin_total_errors_div_bin_width,
                        dict(label="Unfolded\n($\\tau = %.3g$)" % (tau),
                             line_color=unfolded_total_colour, line_width=lw, line_style=3,
                             marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=True)),  # unfolded data
                ])

                # Unfolded plots with alternate response matrix results as well
                # --------------------------------------------------------------
                # And alternate MC gen level to compare
                # (i.e. is the difference between unfolded results < diff at gen level?)
                if args.useAltResponse:
                    alt_mc_gen_hist_bin = alt_unfolder.get_var_hist_pt_binned(alt_hist_mc_gen, ibin_pt, binning_scheme="generator")
                    alt_unfolded_hist_bin_total_errors = alt_unfolder.get_var_hist_pt_binned(alt_unfolder.unfolded, ibin_pt, binning_scheme="generator")
                    this_pt_bin_tdir.WriteTObject(alt_unfolded_hist_bin_total_errors, "alt_unfolded_hist_bin_total_errors")

                    alt_gen_colour = ROOT.kBlack
                    alt_colour = ROOT.kBlue+1

                    entries = [
                        Contribution(mc_gen_hist_bin,
                                     label="Generator (%s)" % (region['mc_label']),
                                     line_color=gen_colour, line_width=lw,
                                     marker_color=gen_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(alt_mc_gen_hist_bin,
                                     label="Generator (%s)" % (region['alt_mc_label']),
                                     line_color=alt_gen_colour, line_width=lw,
                                     marker_color=alt_gen_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin,
                                     normalise_hist=True),
                        Contribution(unfolded_hist_bin_stat_errors,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (tau, region['mc_label']),
                                     line_color=unfolded_stat_colour, line_width=lw, line_style=2,
                                     marker_color=unfolded_stat_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin,
                                     normalise_hist=True),
                        Contribution(alt_unfolded_hist_bin_total_errors,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (alt_unfolder.tau, region['alt_mc_label']),
                                     line_color=alt_colour, line_width=lw, line_style=3,
                                     marker_color=alt_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=particle_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Unfolded / Gen}{(%s)}' % (region['mc_label']),
                                **common_hist_args)
                    plot.legend.SetX1(0.55)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.88)
                    plot.plot("NOSTACK E1")
                    # plot.save("%s/unfolded_%s_alt_response_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    # Do a version where divided by bin width
                    # Note that inputs are already normalised to 1
                    # Do not use normalise_hist!
                    alt_mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(alt_mc_gen_hist_bin)
                    alt_unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(alt_unfolded_hist_bin_total_errors)
                    unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                    # unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)

                    entries = [
                        Contribution(mc_gen_hist_bin_div_bin_width,
                                     label="Generator (%s)" % (region['mc_label']),
                                     line_color=gen_colour, line_width=lw,
                                     marker_color=gen_colour, marker_size=0,
                                     normalise_hist=False),
                        Contribution(alt_mc_gen_hist_bin_div_bin_width,
                                     label="Generator (%s)" % (region['alt_mc_label']),
                                     line_color=alt_gen_colour, line_width=lw,
                                     marker_color=alt_gen_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin_div_bin_width,
                                     normalise_hist=False),
                        Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (tau, region['mc_label']),
                                     line_color=unfolded_stat_colour, line_width=lw, line_style=2,
                                     marker_color=unfolded_stat_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin_div_bin_width,
                                     normalise_hist=False),
                        Contribution(alt_unfolded_hist_bin_total_errors_div_bin_width,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (alt_unfolder.tau, region['alt_mc_label']),
                                     line_color=alt_colour, line_width=lw, line_style=3,
                                     marker_color=alt_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin_div_bin_width,
                                     normalise_hist=False),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=particle_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Unfolded / Gen}{(%s)}' % (region['mc_label']),
                                **common_hist_args)
                    plot.legend.SetX1(0.55)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.88)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/unfolded_%s_alt_response_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Unfolded plots with variations in input model systematics plotted
                # --------------------------------------------------------------
                if args.doModelSysts:
                    syst_entries = []
                    syst_entries_div_bin_width = []
                    for syst_dict in region['model_systematics']:
                        syst_label = syst_dict['label']
                        syst_label_no_spaces = syst_dict['label'].replace(" ", "_")
                        syst_tau = syst_dict['tau']
                        syst_unfolded_1d = syst_dict.get('unfolded_1d', None)
                        if not syst_unfolded_1d:
                            continue
                        syst_unfolded_hist_bin_total_errors = unfolder.get_var_hist_pt_binned(syst_unfolded_1d, ibin_pt, binning_scheme="generator")
                        this_pt_bin_tdir.WriteTObject(syst_unfolded_hist_bin_total_errors, "syst_%s_unfolded_hist_bin_total_errors" % (syst_label_no_spaces))

                        syst_gen_1d = syst_dict.get('gen_1d', None)
                        if not syst_gen_1d:
                            continue
                        syst_gen_1d_bin = unfolder.get_var_hist_pt_binned(syst_gen_1d, ibin_pt, binning_scheme="generator")
                        this_pt_bin_tdir.WriteTObject(syst_gen_1d_bin, "syst_%s_gen_hist_bin" % (syst_label_no_spaces))

                        syst_entries.extend([
                            Contribution(syst_unfolded_hist_bin_total_errors,
                                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=syst_gen_1d_bin,
                                         normalise_hist=True),
                            Contribution(syst_gen_1d_bin,
                                         label="Generator (%s)" % (syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=2,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         normalise_hist=True),
                        ])
                        # already normalised to 1
                        # do not use normalise_hist!
                        syst_unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(syst_unfolded_hist_bin_total_errors)
                        syst_gen_1d_bin_div_bin_width = qgp.hist_divide_bin_width(syst_gen_1d_bin)
                        syst_entries_div_bin_width.extend([
                            Contribution(syst_unfolded_hist_bin_total_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=syst_gen_1d_bin_div_bin_width,
                                         normalise_hist=False),
                            Contribution(syst_gen_1d_bin_div_bin_width,
                                         label="Generator (%s)" % (syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=2,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         normalise_hist=False),
                        ])

                    if len(syst_entries):
                        entries = [
                            Contribution(mc_gen_hist_bin,
                                         label="Generator (%s)" % (region['mc_label']),
                                         line_color=gen_colour, line_width=lw,
                                         marker_color=gen_colour, marker_size=0,
                                         normalise_hist=True),
                            Contribution(unfolded_hist_bin_total_errors,
                                         label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                         line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                         subplot=mc_gen_hist_bin,
                                         normalise_hist=True),
                            Contribution(unfolded_hist_bin_stat_errors,
                                         label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                         line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_stat_colour, marker_size=0,
                                         normalise_hist=True),
                        ]
                        entries.extend(syst_entries)
                        if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                            continue
                        plot = Plot(entries,
                                    xtitle=particle_title,
                                    ytitle=normalised_differential_label,
                                    subplot_title=subplot_title,
                                    **common_hist_args)
                        plot.legend.SetX1(0.55)
                        plot.legend.SetY1(0.72)
                        plot.legend.SetX2(0.98)
                        plot.legend.SetY2(0.88)
                        plot.legend.SetNColumns(3)
                        plot.plot("NOSTACK E1")
                        # plot.save("%s/unfolded_%s_syst_model_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                        # Do a version where divided by bin width
                        mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(mc_gen_hist_bin)
                        unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                        unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)
                        entries_div_bin_width = [
                            Contribution(mc_gen_hist_bin_div_bin_width,
                                         label="Generator (%s)" % (region['mc_label']),
                                         line_color=gen_colour, line_width=lw,
                                         marker_color=gen_colour, marker_size=0,
                                         normalise_hist=False),
                            Contribution(unfolded_hist_bin_total_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                         line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                         subplot=mc_gen_hist_bin_div_bin_width,
                                         normalise_hist=False),
                            Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                         line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_stat_colour, marker_size=0,
                                         normalise_hist=False),
                        ]
                        entries_div_bin_width.extend(syst_entries_div_bin_width)
                        if not check_entries(entries_div_bin_width, "%s %d" % (append, ibin_pt)):
                            continue
                        plot = Plot(entries_div_bin_width,
                                    xtitle=particle_title,
                                    ytitle=normalised_differential_label,
                                    subplot_title=subplot_title,
                                    **common_hist_args)
                        plot.legend.SetX1(0.55)
                        plot.legend.SetY1(0.72)
                        plot.legend.SetX2(0.98)
                        plot.legend.SetY2(0.88)
                        plot.legend.SetNColumns(3)
                        plot.plot("NOSTACK E1")
                        plot.save("%s/unfolded_%s_syst_model_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Unfolded plots with PDF variations
                # --------------------------------------------------------------
                if args.doPDFSysts:
                    syst_entries = []
                    syst_entries_div_bin_width = []
                    for syst_dict in region['pdf_systematics']:
                        syst_label = syst_dict['label']
                        if syst_label.startswith("_"):
                            continue
                        syst_label_no_spaces = syst_dict['label'].replace(" ", "_")
                        syst_tau = syst_dict['tau']
                        syst_unfolded_1d = syst_dict.get('unfolded_1d', None)
                        if not syst_unfolded_1d:
                            continue
                        syst_unfolded_hist_bin_total_errors = unfolder.get_var_hist_pt_binned(syst_unfolded_1d, ibin_pt, binning_scheme="generator")
                        this_pt_bin_tdir.WriteTObject(syst_unfolded_hist_bin_total_errors, "syst_%s_unfolded_hist_bin_total_errors" % (syst_label_no_spaces))

                        syst_gen_1d = syst_dict.get('gen_1d', None)
                        if not syst_gen_1d:
                            continue
                        syst_gen_1d_bin = unfolder.get_var_hist_pt_binned(syst_gen_1d, ibin_pt, binning_scheme="generator")
                        this_pt_bin_tdir.WriteTObject(syst_gen_1d_bin, "syst_%s_gen_hist_bin" % (syst_label_no_spaces))

                        syst_entries.extend([
                            Contribution(syst_unfolded_hist_bin_total_errors,
                                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=syst_gen_1d_bin,
                                         normalise_hist=True),
                            Contribution(syst_gen_1d_bin,
                                         label="Generator (%s)" % (syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=2,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         normalise_hist=True),
                        ])
                        # already normalised to 1
                        # do not use normalise_hist!
                        syst_unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(syst_unfolded_hist_bin_total_errors)
                        syst_gen_1d_bin_div_bin_width = qgp.hist_divide_bin_width(syst_gen_1d_bin)
                        syst_entries_div_bin_width.extend([
                            Contribution(syst_unfolded_hist_bin_total_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=syst_gen_1d_bin_div_bin_width,
                                         normalise_hist=False),
                            Contribution(syst_gen_1d_bin_div_bin_width,
                                         label="Generator (%s)" % (syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=2,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         normalise_hist=False),
                        ])

                    if len(syst_entries):
                        entries = [
                            Contribution(mc_gen_hist_bin,
                                         label="Generator (%s)" % (region['mc_label']),
                                         line_color=gen_colour, line_width=lw,
                                         marker_color=gen_colour, marker_size=0,
                                         normalise_hist=True),
                            Contribution(unfolded_hist_bin_total_errors,
                                         label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                         line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                         subplot=mc_gen_hist_bin,
                                         normalise_hist=True),
                            Contribution(unfolded_hist_bin_stat_errors,
                                         label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                         line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_stat_colour, marker_size=0,
                                         normalise_hist=True),
                        ]
                        entries.extend(syst_entries)
                        if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                            continue
                        plot = Plot(entries,
                                    xtitle=particle_title,
                                    ytitle=normalised_differential_label,
                                    subplot_title=subplot_title,
                                    **common_hist_args)
                        plot.legend.SetX1(0.55)
                        plot.legend.SetY1(0.72)
                        plot.legend.SetX2(0.98)
                        plot.legend.SetY2(0.88)
                        plot.legend.SetNColumns(3)
                        plot.plot("NOSTACK E1")
                        # plot.save("%s/unfolded_%s_syst_model_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                        # Do a version where divided by bin width
                        mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(mc_gen_hist_bin)
                        unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                        unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)
                        entries_div_bin_width = [
                            Contribution(mc_gen_hist_bin_div_bin_width,
                                         label="Generator (%s)" % (region['mc_label']),
                                         line_color=gen_colour, line_width=lw,
                                         marker_color=gen_colour, marker_size=0,
                                         normalise_hist=False),
                            Contribution(unfolded_hist_bin_total_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                         line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                         subplot=mc_gen_hist_bin_div_bin_width,
                                         normalise_hist=False),
                            Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                         line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_stat_colour, marker_size=0,
                                         normalise_hist=False),
                        ]
                        entries_div_bin_width.extend(syst_entries_div_bin_width)
                        if not check_entries(entries_div_bin_width, "%s %d" % (append, ibin_pt)):
                            continue
                        plot = Plot(entries_div_bin_width,
                                    xtitle=particle_title,
                                    ytitle=normalised_differential_label,
                                    subplot_title=subplot_title,
                                    **common_hist_args)
                        plot.legend.SetX1(0.55)
                        plot.legend.SetY1(0.72)
                        plot.legend.SetX2(0.98)
                        plot.legend.SetY2(0.88)
                        plot.legend.SetNColumns(3)
                        plot.plot("NOSTACK E1")
                        plot.save("%s/unfolded_%s_pdf_model_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))


                # --------------------------------------------------------------
                # PLOT UNCERTAINTY SHIFTS
                # --------------------------------------------------------------
                systematic_shift_hists_bin = [unfolder.get_var_hist_pt_binned(h, ibin_pt, binning_scheme='generator')
                                              for h in systematic_shift_hists]
                unfolded_stat_error_bin = unfolder.get_var_hist_pt_binned(unfolder.unfolded_stat_err, ibin_pt, binning_scheme="generator")
                unfolded_total_error_bin =  unfolder.get_var_hist_pt_binned(unfolder.unfolded, ibin_pt, binning_scheme="generator")
                plot_uncertainty_shifts(total_hist=unfolded_total_error_bin,
                                        stat_hist=unfolded_stat_error_bin,
                                        syst_shifts=systematic_shift_hists_bin,
                                        systs=region['experimental_systematics'],
                                        output_filename='%s/unfolded_systs_%s_bin_%d.%s' % (this_output_dir, append, ibin_pt, OUTPUT_FMT),
                                        title=title,
                                        angle_str=angle_str)

                # --------------------------------------------------------------
                # PLOT RECO-LEVEL DISTRIBUTIONS (with gen binning)
                # --------------------------------------------------------------
                # Reco level, generator-binning
                reco_mc_colour = ROOT.kGreen+2
                reco_data_colour = ROOT.kRed
                entries = [
                    Contribution(mc_reco_hist_bin_gen_binning,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(reco_hist_bin_gen_binning,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_gen_binning,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                # plot.save("%s/detector_gen_binning_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Again but divided by bin width
                mc_reco_hist_bin_gen_binning_div_bin_width = qgp.hist_divide_bin_width(mc_reco_hist_bin_gen_binning)
                reco_hist_bin_gen_binning_div_bin_width = qgp.hist_divide_bin_width(reco_hist_bin_gen_binning)
                entries = [
                    Contribution(mc_reco_hist_bin_gen_binning_div_bin_width,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(reco_hist_bin_gen_binning_div_bin_width,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_gen_binning_div_bin_width,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_gen_binning_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                if SUBTRACT_FAKES:
                    # Same but background-subtracted
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_gen_binning,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(reco_hist_bg_subtracted_bin_gen_binning,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_gen_binning,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    # plot.save("%s/detector_gen_binning_bg_subtracted_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    mc_reco_hist_bg_subtracted_bin_gen_binning_div_bin_width = qgp.hist_divide_bin_width(mc_reco_hist_bg_subtracted_bin_gen_binning)
                    reco_hist_bg_subtracted_bin_gen_binning_div_bin_width = qgp.hist_divide_bin_width(reco_hist_bg_subtracted_bin_gen_binning)
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_gen_binning_div_bin_width,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=False),
                        Contribution(reco_hist_bg_subtracted_bin_gen_binning_div_bin_width,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_gen_binning_div_bin_width,
                                     normalise_hist=False),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/detector_gen_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))


            # Draw individual pt bin plots - RECO binning
            # ------------------------------------------------------------------
            for ibin_pt in range(0, len(unfolder.pt_bin_edges_reco)-1):
                this_pt_bin_tdir = this_tdir.mkdir("reco_bin_%d" % (ibin_pt))
                print("Individual detector bin", ibin_pt, "=", unfolder.pt_bin_edges_reco[ibin_pt], unfolder.pt_bin_edges_reco[ibin_pt+1])

                # Get 1D hists
                # The folded unfolded result
                folded_unfolded_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(unfolder.get_folded_unfolded(), ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(folded_unfolded_hist_bin_reco_binning, "folded_unfolded_hist_bin_reco_binning")

                # here this is the thing to be unfolded, should be data or MC depending on MC_INPUT flag
                reco_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(unfolder.input_hist, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(reco_hist_bin_reco_binning, "reco_hist_bin_reco_binning")

                if SUBTRACT_FAKES:
                    reco_hist_bg_subtracted_bin_reco_binning = unfolder.get_var_hist_pt_binned(unfolder.input_hist_bg_subtracted, ibin_pt, binning_scheme="detector")
                    this_pt_bin_tdir.WriteTObject(reco_hist_bg_subtracted_bin_reco_binning, "reco_hist_bg_subtracted_bin_reco_binning")

                # Get MC hists
                mc_reco_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(unfolder.hist_mc_reco, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(mc_reco_hist_bin_reco_binning, "mc_reco_hist_bin_reco_binning")

                if SUBTRACT_FAKES:
                    mc_reco_hist_bg_subtracted_bin_reco_binning = unfolder.get_var_hist_pt_binned(unfolder.hist_mc_reco_bg_subtracted, ibin_pt, binning_scheme="detector")
                    this_pt_bin_tdir.WriteTObject(mc_reco_hist_bg_subtracted_bin_reco_binning, "mc_reco_hist_bg_subtracted_bin_reco_binning")

                # Do the plots
                # --------------------------------------------------------------
                # common hist settings
                lw = 2
                title = "%s\n%s region\n%g < p_{T}^{Jet} < %g GeV" % (jet_algo, region['label'], unfolder.pt_bin_edges_reco[ibin_pt], unfolder.pt_bin_edges_reco[ibin_pt+1])
                common_hist_args = dict(
                    what="hist",
                    title=title,
                    subplot_type='ratio',
                    # subplot_limits=(0.5, 1.5),
                    subplot_limits=(0.75, 1.25),
                )

                # Reco only, detector-binning
                entries = [
                    Contribution(mc_reco_hist_bin_reco_binning,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(reco_hist_bin_reco_binning,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_reco_binning,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                # plot.save("%s/detector_reco_binning_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # same but divided by bin width
                mc_reco_hist_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(mc_reco_hist_bin_reco_binning)
                reco_hist_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(reco_hist_bin_reco_binning)
                entries = [
                    Contribution(mc_reco_hist_bin_reco_binning_div_bin_width,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(reco_hist_bin_reco_binning_div_bin_width,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_reco_binning_div_bin_width,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_reco_binning_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                if SUBTRACT_FAKES:
                    # Same but background-subtracted
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_reco_binning,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(reco_hist_bg_subtracted_bin_reco_binning,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_reco_binning,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    # plot.save("%s/detector_reco_binning_bg_subtracted_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    # again but divided by bin width
                    mc_reco_hist_bg_subtracted_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(mc_reco_hist_bg_subtracted_bin_reco_binning)
                    reco_hist_bg_subtracted_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(reco_hist_bg_subtracted_bin_reco_binning)
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_reco_binning_div_bin_width,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=False),
                        Contribution(reco_hist_bg_subtracted_bin_reco_binning_div_bin_width,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_reco_binning_div_bin_width,
                                     normalise_hist=False),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/detector_reco_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # PLOT FOLDED UNFOLDED DATA
                # --------------------------------------------------------------
                # Reco + folded, detector binning
                # FIXME: shouldn't this use RECO gen bin
                reco_unfolding_input_colour = ROOT.kRed
                reco_folded_unfolded_colour = ROOT.kAzure+1
                reco_folded_truth_colour = ROOT.kGreen+2
                reco_mc_colour = ROOT.kGreen

                mc_reco_bin_hist = mc_reco_hist_bg_subtracted_bin_reco_binning if SUBTRACT_FAKES else mc_reco_hist_bin_reco_binning
                reco_bin_hist = reco_hist_bg_subtracted_bin_reco_binning if SUBTRACT_FAKES else reco_hist_bin_reco_binning

                # Folded, but only comparing data with data to check it is sane
                entries = [
                    Contribution(reco_bin_hist,
                                 label="Unfolding input (bg-subtracted)" if SUBTRACT_FAKES else "Unfolding input",
                                 line_color=reco_unfolding_input_colour, line_width=lw,
                                 marker_color=reco_unfolding_input_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(folded_unfolded_hist_bin_reco_binning,
                                 label="Folded unfolded result (#tau = %.3g)" % (tau),
                                 line_color=reco_folded_unfolded_colour, line_width=lw,
                                 marker_color=reco_folded_unfolded_colour, marker_size=0,
                                 subplot=reco_bin_hist,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Folded unfolded / reco',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                # plot.save("%s/detector_folded_unfolded_only_data_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Same but divided by bin width
                # Do not normalise again!
                # data + folded + MC
                mc_reco_bin_hist_div_bin_width = qgp.hist_divide_bin_width(mc_reco_bin_hist)
                reco_bin_hist_div_bin_width = qgp.hist_divide_bin_width(reco_bin_hist)
                folded_unfolded_hist_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(folded_unfolded_hist_bin_reco_binning)
                entries = [
                    Contribution(mc_reco_bin_hist_div_bin_width,
                                 label="MC (reco, bg-subtracted)" if SUBTRACT_FAKES else "MC (reco)",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(reco_bin_hist_div_bin_width,
                                 label="Unfolding input (bg-subtracted)" if SUBTRACT_FAKES else "Unfolding input",
                                 line_color=reco_unfolding_input_colour, line_width=lw,
                                 marker_color=reco_unfolding_input_colour, marker_size=0,
                                 subplot=mc_reco_bin_hist_div_bin_width,
                                 normalise_hist=False),
                    Contribution(folded_unfolded_hist_bin_reco_binning_div_bin_width,
                                 label="Folded unfolded data (#tau = %.3g)" % (tau),
                                 line_color=reco_folded_unfolded_colour, line_width=lw,
                                 marker_color=reco_folded_unfolded_colour, marker_size=0,
                                 subplot=mc_reco_bin_hist_div_bin_width,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.56)
                plot.legend.SetY1(0.72)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                # plot.save("%s/detector_folded_unfolded_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # data + folded
                entries = [
                    Contribution(reco_bin_hist_div_bin_width,
                                 label="Unfolding input (bg-subtracted)" if SUBTRACT_FAKES else "Unfolding input",
                                 line_color=reco_unfolding_input_colour, line_width=lw,
                                 marker_color=reco_unfolding_input_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(folded_unfolded_hist_bin_reco_binning_div_bin_width,
                                 label="Folded unfolded data (#tau = %.3g)" % (tau),
                                 line_color=reco_folded_unfolded_colour, line_width=lw,
                                 marker_color=reco_folded_unfolded_colour, marker_size=0,
                                 subplot=reco_bin_hist_div_bin_width,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='#splitline{Folded unfolded}{/ reco}',
                            **common_hist_args)
                plot.legend.SetX1(0.56)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_folded_unfolded_only_data_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # PLOT FOLDED GEN MC
                # --------------------------------------------------------------
                if MC_INPUT and unfolder.get_folded_mc_truth():
                    # Get 1D hists
                    # The folded result
                    folded_gen_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(unfolder.get_folded_mc_truth(), ibin_pt, binning_scheme="detector")
                    this_pt_bin_tdir.WriteTObject(folded_gen_hist_bin_reco_binning, "folded_gen_hist_bin_reco_binning")

                    # Folded gen, comparing to original reco
                    entries = [
                        Contribution(mc_reco_bin_hist,
                                     label="Unfolding input (bg-subtracted)" if SUBTRACT_FAKES else "Unfolding input",
                                     line_color=reco_unfolding_input_colour, line_width=lw,
                                     marker_color=reco_unfolding_input_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(folded_gen_hist_bin_reco_binning,
                                     label="Folded particle-level MC",
                                     line_color=reco_folded_truth_colour, line_width=lw,
                                     marker_color=reco_folded_truth_colour, marker_size=0,
                                     subplot=reco_bin_hist,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Folded gen}{/ reco}',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    # plot.save("%s/detector_folded_gen_only_data_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    # Same but divided by bin width
                    # Do not normalise again!
                    mc_reco_bin_hist_div_bin_width = qgp.hist_divide_bin_width(mc_reco_bin_hist)
                    folded_gen_hist_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(folded_gen_hist_bin_reco_binning)
                    entries = [
                        Contribution(mc_reco_bin_hist_div_bin_width,
                                     label="Unfolding input (bg-subtracted)" if SUBTRACT_FAKES else "Unfolding input",
                                     line_color=reco_unfolding_input_colour, line_width=lw,
                                     marker_color=reco_unfolding_input_colour, marker_size=0,
                                     normalise_hist=False),
                        Contribution(folded_gen_hist_bin_reco_binning_div_bin_width,
                                     label="Folded particle-level MC",
                                     line_color=reco_folded_truth_colour, line_width=lw,
                                     marker_color=reco_folded_truth_colour, marker_size=0,
                                     subplot=mc_reco_bin_hist_div_bin_width,
                                     normalise_hist=False),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Folded gen}{/ reco}',
                                **common_hist_args)
                    plot.legend.SetX1(0.56)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/detector_folded_gen_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

            # DO SUMMARY PLOT
            # ------------------------------------------------------------------
            if args.doSummaryPlot:
                marker = ""
                if "_" in angle.name or "^" in angle.name:
                    marker = "$"
                var_label = "Particle-level " + marker + angle.name + marker + " ($%s$)" % angle.lambda_str
                v = "%s_vs_pt" % (angle.var)
                bins = [(pt_bin_edges_gen[i], pt_bin_edges_gen[i+1]) for i in range(len(pt_bin_edges_gen)-1)]
                print("Making summary plot from pt bins:", bins)
                xlim = (50, 614) if "ZPlusJets" in region['name'] else (50, 2000)
                region_label = region['label'].replace("Dijet", "dijet")  # to ensure correct capitalisation
                qgp.do_mean_rms_summary_plot(summary_1d_entries, bins,
                                             "%s/%s_box_dijet_mpl.%s" % (this_output_dir, v, OUTPUT_FMT),
                                             var_label=var_label,
                                             xlim=xlim,
                                             region_title=region_label)

    print("Saved hists to", output_tfile.GetName())
