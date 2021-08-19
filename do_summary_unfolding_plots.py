#!/usr/bin/env python


"""
Do all the unfolding plots: per pT bin, per lambda bin, summary plot

Input is HDF5 file made by extract_unfolding_summary_stats.py
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
import pandas as pd
import numpy as np
from itertools import product
from array import array
from copy import copy
import warnings

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from unfolding_config import get_dijet_config, get_zpj_config
import rivet_naming as rn
import metric_calculators as metrics
from extract_rivet_summary_stats import dataframe_yoda_key, convert_df_types

# monkey-patch warning formatter
warnings.formatwarning = cu._formatwarning

ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()

# import hunter
# hunter.trace(module='my_unfolder', action=hunter.CallPrinter)

# Define own linestyle for smaller plots
# linestyle 2 (dashed) has too big dashes
ROOT.gStyle.SetLineStyleString(22, "8 6")
ROOT.gStyle.SetLineStyleString(23, "5 10")

COMMON_STYLE_DICT = {
    "line_width": 2,
    "dijet_cen_color": ROOT.kBlack,
    "dijet_fwd_color": ROOT.kRed,
    "zpj_color": ROOT.kBlue,
    "jet_pt_units_str": "Jet p_{T} [GeV]",

    "data_line_style": 1,
    "data_color": ROOT.kBlack,
    "data_marker_style": cu.Marker.get("circle", filled=True),

    "marker_size": 1.3,

    # "mc_line_style": 22,
    "mc_line_style": 1,
    # format: 3ijk,
    # i=distance between lines,
    # j=angle between 0 and 90 degrees (5 = not drawn),
    # k=angle between 90 and 180 degrees (5 = not drawn)
    "mc_fill_style": 3445,
    "mc_color": qgc.MGPY_QCD_COLOUR,
    "mc_marker_style": cu.Marker.get('square', filled=False),

    # "mc_alt_line_style": 22,
    "mc_alt_line_style": 1,
    "mc_alt_fill_style": 3454,
    # "mc_alt_color": ROOT.kAzure+1,
    "mc_alt_color": qgc.HERWIGPP_QCD_COLOUR,
    "mc_alt_marker_style": cu.Marker.get('triangleUp', filled=False),
}


# some pre-defined sample dicts for YODA files
# TODO: better just to add as specific commandline args?
SAMPLE_STYLE_DICTS = {
    dataframe_yoda_key("Sherpa"): {
        "color": ROOT.kGreen+2,
        "label": "Sherpa LO",
        "marker_style": cu.Marker.get('triangleDown', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1,
    },

    dataframe_yoda_key("Sherpa LO"): {
        "color": ROOT.kGreen+2,
        "label": "Sherpa LO",
        "marker_style": cu.Marker.get('triangleDown', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1,
    },

    dataframe_yoda_key("Sherpa LO+jet"): {
        "color": ROOT.kViolet-6,
        "label": "Sherpa LO+jet",
        "marker_style": cu.Marker.get('triangleUp', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1,
    },

    dataframe_yoda_key("Herwig7 CH3"): {
        "color": ROOT.kOrange-3,
        "label": "Herwig7 CH3",
        "marker_style": cu.Marker.get('star', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1.2, # these shapes always come out small
    },

    dataframe_yoda_key("Herwig7 CH3 alphaS=0.136"): {
        "color": ROOT.kRed-2,
        "label": "Herwig7 CH3\n#alpha_{S} = 0.136",
        "marker_style": cu.Marker.get('diamond', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1.2, # these shapes always come out small
    },

    dataframe_yoda_key("Pythia8 CP5"): {
        # "color": ROOT.kMagenta-7,
        # "color": ROOT.kViolet-1,
        "color": ROOT.kRed-4,
        "label": "Pythia8 CP5",
        "marker_style": cu.Marker.get('crossX', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1.2, # these shapes always come out small
    },

    dataframe_yoda_key("Pythia8 CP2"): {
        # "color": ROOT.kMagenta-7,
        # "color": ROOT.kViolet-1,
        "color": ROOT.kAzure+1,
        "label": "Pythia8 CP2",
        "marker_style": cu.Marker.get('doubleDiamond', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1.2, # these shapes always come out small
    },
}


def create_angle_label(angle, do_groomed=False):
    angle_prepend = "Groomed " if do_groomed else ""
    if "charged" in angle.var:
        if angle_prepend != "":
            angle_prepend += "c"
        else:
            angle_prepend += "C"
        angle_prepend += "harged-only "

    # for plot axis titles
    angle_str = "{prepend}{lambda_str}".format(prepend=angle_prepend,
                                                 lambda_str=angle.lambda_str)
    return angle_str


class SummaryPlotter(object):
    """Do lots of summary plots"""

    def __init__(self, jet_algos, angles, pt_bins_dijet, pt_bins_zpj, df, output_dir, has_data, is_preliminary, is_supplementary, only_yoda_data):
        if len(jet_algos) == 0:
            raise RuntimeError("jet_algos is empty")
        self.jet_algos = jet_algos
        if len(angles) == 0:
            raise RuntimeError("angles is empty")
        self.angles = angles
        self.pt_bins_dijet = pt_bins_dijet
        self.pt_bins_zpj = pt_bins_zpj
        self.min_pt_bin_ind_ak4 = 0
        self.min_pt_bin_ind_ak8 = 2  # lose first 2 bins of AK8 as MC not proper
        self.df = df
        self.output_fmt = 'pdf'
        self.output_dir = output_dir
        self.has_data = has_data
        self.is_preliminary = is_preliminary
        self.is_supplementary = is_supplementary
        # use kerning to avoid splitline taking up too much space
        # lower the whole thing a little to avoid clashing with hashed bit in plots with ratio
        self.mc_label = '#lower[0.1]{#splitline{MG5+Pythia8}{#lower[-0.15]{CUETP8M1}}}'
        self.mc_label_short = 'MG5+Pythia8'
        self.alt_mc_label = 'Herwig++'
        self.other_samples = []
        self.only_yoda_data = only_yoda_data
        self.filename_append = "_onlyYodaData" if self.only_yoda_data else ""
        if not self.is_preliminary:
            self.filename_append += "_paper"

    def add_sample(self, key, style_dict):
        self.other_samples.append({
            "key": key, # for accesing the variable in dataframe
            "style_dict": style_dict
        })

    def min_pt_bin_ind(self, jet_algo):
        if "ak4" in jet_algo['name'].lower():
            return self.min_pt_bin_ind_ak4
        else:
            return self.min_pt_bin_ind_ak8

    @staticmethod
    def data_to_hist(data, data_err, bins):
        if len(data) != len(bins)-1:
            print("data:", data)
            print("bins:", bins)
            raise RuntimeError("len(data) != len(bins)-1")
        h = ROOT.TH1D("h_"+cu.get_unique_str(), "", len(bins)-1, array('d', bins))
        for ind, (y, err) in enumerate(zip(data, data_err), 1):
            h.SetBinContent(ind, y)
            h.SetBinError(ind, err)
        return h

    @staticmethod
    def _generate_filename_prefix(do_dijet_cen, do_dijet_fwd, do_zpj):
        this_str = ""
        if do_dijet_cen:
            this_str += "dijet_cen_"
        if do_dijet_fwd:
            this_str += "dijet_fwd_"
        if do_zpj:
            this_str += "zpj_"
        return this_str


    def plot_dijet_zpj_metric_vs_pt_one_angle_one_jet(self, metric, angle, jet_algo, do_groomed, output_dir, do_zpj=True, do_dijet_cen=True, do_dijet_fwd=True, min_pt_bin_ind=0):
        """Do plot of lambda metric vs pt, for any combo of the dijet cen/fwd and zpj regions, for a given angle/jet algo/grooming"""
        if metric not in ["mean", "rms", "delta"]:
            raise ValueError("metric must be one of 'mean', 'rms', 'delta'")

        only_one_region = sum([do_dijet_cen, do_dijet_fwd, do_zpj]) == 1

        df = self.df
        mask = ((df['angle'] == angle.var) & (df['jet_algo'] == jet_algo['name']) & (df['isgroomed'] == do_groomed) & (df['pt_bin'] >= min_pt_bin_ind))
        if not mask.any():
            return

        # create hists for dijet central
        dijet_cen_col = COMMON_STYLE_DICT['dijet_cen_color']
        other_samples_dijet_central_hists = []
        dijet_central_hist_no_errors = None
        dijet_central_hist_truth = None
        dijet_central_hist_alt_truth = None
        dijet_central_hist = None
        dijet_central_hist_ratio_error = None
        if do_dijet_cen:
            region_name = 'Dijet_central'
            if do_groomed:
                region_name += "_groomed"
            dijet_central_data = df[mask & (df['region'] == region_name)]

            if not self.only_yoda_data:
                dijet_central_hist_truth = self.data_to_hist(dijet_central_data['%s_truth' % metric], dijet_central_data['%s_err_truth' % metric], self.pt_bins_dijet[min_pt_bin_ind:])
                dijet_central_hist_alt_truth = self.data_to_hist(dijet_central_data['%s_alt_truth' % metric], dijet_central_data['%s_err_alt_truth' % metric], self.pt_bins_dijet[min_pt_bin_ind:])

            for sample in self.other_samples:
                hist = self.data_to_hist(dijet_central_data['%s_%s' % (metric, sample['key'])], dijet_central_data['%s_err_%s' % (metric, sample['key'])], self.pt_bins_dijet[min_pt_bin_ind:])
                other_samples_dijet_central_hists.append(hist)

            if metric != 'delta':
                dijet_central_hist = self.data_to_hist(dijet_central_data[metric], dijet_central_data['%s_err' % metric], self.pt_bins_dijet[min_pt_bin_ind:])
                # Create copy with 0 error bars, for the ratio subplot
                # The data will have its own error region
                dijet_central_hist_no_errors = dijet_central_hist.Clone()
                cu.remove_th1_errors(dijet_central_hist_no_errors)
                # Create hists for data with error reigon for ratio
                # Easiest way to get errors right is to do data (with 0 errors)
                # and divide by data (with errors), as if you had MC = data with 0 error
                dijet_central_hist_ratio_error = dijet_central_hist_no_errors.Clone()
                dijet_central_hist_ratio_error.Divide(dijet_central_hist)
                dijet_central_hist_ratio_error.SetFillStyle(3354)
                dijet_central_hist_ratio_error.SetFillColor(16)
                dijet_central_hist_ratio_error.SetLineWidth(0)
                dijet_central_hist_ratio_error.SetMarkerSize(0)

        # create hists for dijet forward
        dijet_fwd_col = COMMON_STYLE_DICT['dijet_fwd_color']
        other_samples_dijet_forward_hists = []
        dijet_forward_hist_no_errors = None
        dijet_forward_hist_truth = None
        dijet_forward_hist_alt_truth = None
        dijet_forward_hist = None
        dijet_forward_hist_ratio_error = None
        if do_dijet_fwd:
            region_name = 'Dijet_forward'
            if do_groomed:
                region_name += "_groomed"
            dijet_forward_data = df[mask & (df['region'] == region_name)]

            if not self.only_yoda_data:
                dijet_forward_hist_truth = self.data_to_hist(dijet_forward_data['%s_truth' % metric], dijet_forward_data['%s_err_truth' % metric], self.pt_bins_dijet[min_pt_bin_ind:])
                dijet_forward_hist_alt_truth = self.data_to_hist(dijet_forward_data['%s_alt_truth' % metric], dijet_forward_data['%s_err_alt_truth' % metric], self.pt_bins_dijet[min_pt_bin_ind:])

            for sample in self.other_samples:
                hist = self.data_to_hist(dijet_forward_data['%s_%s' % (metric, sample['key'])], dijet_forward_data['%s_err_%s' % (metric, sample['key'])], self.pt_bins_dijet[min_pt_bin_ind:])
                other_samples_dijet_forward_hists.append(hist)

            if metric != 'delta':
                dijet_forward_hist = self.data_to_hist(dijet_forward_data[metric], dijet_forward_data['%s_err' % metric], self.pt_bins_dijet[min_pt_bin_ind:])
                dijet_forward_hist_no_errors = dijet_forward_hist.Clone()
                cu.remove_th1_errors(dijet_forward_hist.Clone())

                dijet_forward_hist_ratio_error = dijet_forward_hist_no_errors.Clone()
                dijet_forward_hist_ratio_error.Divide(dijet_forward_hist)
                dijet_forward_hist_ratio_error.SetFillStyle(3354)
                dijet_forward_hist_ratio_error.SetFillColor(16)
                dijet_forward_hist_ratio_error.SetLineWidth(0)
                dijet_forward_hist_ratio_error.SetMarkerSize(0)

        # create hists for Z+jets
        zpj_col = COMMON_STYLE_DICT['zpj_color']
        other_samples_zpj_hists = []
        zpj_hist_no_errors = None
        zpj_hist_truth = None
        zpj_hist_alt_truth = None
        zpj_hist_theory = None
        zpj_hist_sherpa = None
        zpj_hist = None
        zpj_hist_ratio_error = None
        if do_zpj:
            region_name = 'ZPlusJets'
            if do_groomed:
                region_name += "_groomed"
            # drop last pt bin as massive error
            pt_bins = self.pt_bins_zpj[:]
            zpj_data = df[mask & (df['region'] == region_name) & (df['pt_bin'] < len(pt_bins)-1) & (df['pt_bin'] >= min_pt_bin_ind)]
            # the -1 is because the last entry in pt_bins is the upper edge of the last bin
            if not self.only_yoda_data:
                zpj_hist_truth = self.data_to_hist(zpj_data['%s_truth' % metric], zpj_data['%s_err_truth' % metric], pt_bins[min_pt_bin_ind:])
                zpj_hist_alt_truth = self.data_to_hist(zpj_data['%s_alt_truth' % metric], zpj_data['%s_err_alt_truth' % metric], pt_bins[min_pt_bin_ind:])
                #### HERE IS WHERE ANDREAS ADDED THE THEORY
                ### GET THEORY CURVE
                zpj_hist_theory = zpj_hist_truth.Clone()
                zpj_hist_theory_error = zpj_hist_truth.Clone()
                zpj_hist_theory_upper = zpj_hist_truth.Clone()
                zpj_hist_theory_lower = zpj_hist_truth.Clone()
                zpj_hist_sherpa_error = zpj_hist_truth.Clone()
                var_number=["jet_puppiMultiplicity", "jet_pTD","jet_LHA", "jet_width", "jet_thrust", 
            		"jet_puppiMultiplicity_charged", "jet_pTD_charged", "jet_LHA_charged", "jet_width_charged", "jet_thrust_charged"].index(angle.var)+1
                if var_number in [3,4,5,8,9,10]:
                 angle_output_dir = "%s/%s/%s" % ("/nfs/dust/cms/user/hinzmann/qganalysis/CMSSW_10_2_17/src/trees/unfolded_ak"+jet_algo['label'][-1]+"puppi", region_name, angle.var)
                 root_filename = os.path.join(angle_output_dir, "unfolding_result_slim.root")
                 import uproot
                 uproot_file = uproot.open(root_filename)
                 unfolding_dict = unpack_slim_unfolding_root_file_uproot(uproot_file, region_name, angle.var, pt_bins[min_pt_bin_ind:])
                 for ibin in range(len(pt_bins[min_pt_bin_ind:])-1):
                  errors=np.ndarray(1)
                  areas=np.ndarray(1)
                  centers=np.ndarray(1)
                  errors=np.delete(errors,0)
                  areas=np.delete(areas,0)
                  centers=np.delete(centers,0)
                  #print(min_pt_bin_ind, ibin)
                  for post in ["_fix",""]:
                    theory_file="/nfs/dust/cms/user/hinzmann/qganalysis/CMSSW_10_2_17/src/zjet_angularities-master-to_send_to_CMS-thr_vs_exp_cvt"+post+"/to_send_to_CMS/thr_vs_exp_cvt"+post+"/"
                    theory_file+="data_R"+jet_algo['label'][-1]+"/"
                    theory_file+="RSIG_THR_RES_NP_"
                    theory_file+="d"+("1" if "4" in jet_algo['label'] else "2")+("2" if do_groomed else "1")+"-"
                    theory_file+="x"+("0" if var_number<10 else "")+str(var_number)+"-"
                    theory_file+="y"+("0" if ibin<9 else "")+str(ibin+1)
                    theory_file+=".dat"
                    if os.path.exists(theory_file): break
                  #print(theory_file)
                  igenbin=0
                  for line in open(theory_file).readlines():
                    if line.count(".")==5:
                      igenbin+=1
                      width=unfolding_dict['truth_hists'][ibin].edges[igenbin]-unfolding_dict['truth_hists'][ibin].edges[igenbin-1]
                      errors=np.append(errors,float(line.split("\t")[3])*width)
                      areas=np.append(areas,float(line.split("\t")[2])*width)
                      centers=np.append(centers,unfolding_dict['truth_hists'][ibin].edges[igenbin-1]+width/2.)
                      assert(float(line.split("\t")[0])==unfolding_dict['truth_hists'][ibin].edges[igenbin-1])
                  mean = metrics.calc_mean_jax(areas, centers)
                  err = metrics.calc_mean_uncorrelated_error_jax(areas, centers, errors)
                  zpj_hist_theory.SetBinContent(ibin+1,mean)
                  zpj_hist_theory.SetBinError(ibin+1,1e-100)
                  zpj_hist_theory_error.SetBinContent(ibin+1,mean)
                  zpj_hist_theory_error.SetBinError(ibin+1,err)
                  zpj_hist_theory_upper.SetBinContent(ibin+1,mean+err)
                  zpj_hist_theory_upper.SetBinError(ibin+1,0)
                  zpj_hist_theory_lower.SetBinContent(ibin+1,mean-err)
                  zpj_hist_theory_lower.SetBinError(ibin+1,0)

                  for post in ["_fix",""]:
                    sherpa_file="/nfs/dust/cms/user/hinzmann/qganalysis/CMSSW_10_2_17/src/zjet_angularities-master-to_send_to_CMS-thr_vs_exp_cvt"+post+"/to_send_to_CMS/thr_vs_exp_cvt"+post+"/"
                    sherpa_file+="data_R"+jet_algo['label'][-1]+"/"
                    sherpa_file+="SIG_SRP_NP_"
                    sherpa_file+="d"+("1" if "4" in jet_algo['label'] else "2")+("2" if do_groomed else "1")+"-"
                    sherpa_file+="x"+("0" if var_number<10 else "")+str(var_number)+"-"
                    sherpa_file+="y"+("0" if ibin<9 else "")+str(ibin+1)
                    sherpa_file+=".dat"
                    if os.path.exists(sherpa_file): break
                  #print(sherpa_file)
                  igenbin=0
                  for line in open(sherpa_file).readlines():
                    if line.count(".")==5:
                      igenbin+=1
                      width=unfolding_dict['truth_hists'][ibin].edges[igenbin]-unfolding_dict['truth_hists'][ibin].edges[igenbin-1]
                      errors=np.append(errors,float(line.split("\t")[3])*width)
                      areas=np.append(areas,float(line.split("\t")[2])*width)
                      centers=np.append(centers,unfolding_dict['truth_hists'][ibin].edges[igenbin-1]+width/2.)
                      assert(float(line.split("\t")[0])==unfolding_dict['truth_hists'][ibin].edges[igenbin-1])
                  mean = metrics.calc_mean_jax(areas, centers)
                  err = metrics.calc_mean_uncorrelated_error_jax(areas, centers, errors)
                  zpj_hist_sherpa_error.SetBinContent(ibin+1,mean)
                  zpj_hist_sherpa_error.SetBinError(ibin+1,err)
            
            for sample in self.other_samples:
                key = sample['key']
                hist = self.data_to_hist(zpj_data['%s_%s' % (metric, key)], zpj_data['%s_err_%s' % (metric, key)], pt_bins[min_pt_bin_ind:])
                other_samples_zpj_hists.append(hist)

            if metric != 'delta':
                zpj_hist = self.data_to_hist(zpj_data[metric], zpj_data['%s_err' % metric], pt_bins[min_pt_bin_ind:])
                # 0 error bar hists & ratio = 1 hists for subplot
                zpj_hist_no_errors = zpj_hist.Clone()
                cu.remove_th1_errors(zpj_hist_no_errors)
                zpj_hist_ratio_error = zpj_hist_no_errors.Clone()
                zpj_hist_ratio_error.Divide(zpj_hist)
                zpj_hist_ratio_error.SetFillStyle(3354)
                zpj_hist_ratio_error.SetFillColor(16)
                zpj_hist_ratio_error.SetLineWidth(0)
                zpj_hist_ratio_error.SetMarkerSize(0)

        m_size = COMMON_STYLE_DICT['marker_size']
        lw = COMMON_STYLE_DICT['line_width']
        entries = []

        # Create dummy graphs with the same styling to put into the legend,
        # since graphs get error bar ends, but hists don't
        dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        dummy_entries = []
        # NB Spaces in legend labels are important for padding

        mc_label = self.mc_label_short

        if len(self.other_samples) > 0:
            # also write tune if > 1 Pythia
            mc_label = self.mc_label

        if not self.only_yoda_data:
            # Add nominal MC
            if do_dijet_cen:
                cont_args = dict(label=mc_label if only_one_region else '#splitline{ %s  }{ [%s]}' % (qgc.Dijet_CEN_LABEL.replace(" region", ""), self.mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "LEP",
                                 subplot=dijet_central_hist_no_errors)
                entries.append(Contribution(dijet_central_hist_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_dijet_fwd:
                cont_args = dict(label=mc_label if only_one_region else '#splitline{ %s  }{ [%s]}' % (qgc.Dijet_FWD_LABEL.replace(" region", ""), self.mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "LEP",
                                 subplot=dijet_forward_hist_no_errors)
                entries.append(Contribution(dijet_forward_hist_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_zpj:
                cont_args = dict(label=mc_label if only_one_region else '#splitline{ %s  }{ [%s]}' % (qgc.ZpJ_LABEL.replace(" region", ""), self.mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "LEP",
                                 subplot=zpj_hist_no_errors)
                entries.append(Contribution(zpj_hist_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            # add alt MC
            if do_dijet_cen:
                cont_args = dict(label=self.alt_mc_label if only_one_region else '#splitline{ %s  }{ [%s]}' % (qgc.Dijet_CEN_LABEL.replace(" region", ""), self.alt_mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_alt_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "LEP",
                                 subplot=dijet_central_hist_no_errors)
                entries.append(Contribution(dijet_central_hist_alt_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_dijet_fwd:
                cont_args = dict(label=self.alt_mc_label if only_one_region else '#splitline{ %s  }{ [%s]}' % (qgc.Dijet_FWD_LABEL.replace(" region", ""), self.alt_mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_alt_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "LEP",
                                 subplot=dijet_forward_hist_no_errors)
                entries.append(Contribution(dijet_forward_hist_alt_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_zpj:
                cont_args = dict(label=self.alt_mc_label if only_one_region else '#splitline{ %s  }{ [%s]}' % (qgc.ZpJ_LABEL.replace(" region", ""), self.alt_mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_alt_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "LEP",
                                 subplot=zpj_hist_no_errors)
                entries.append(Contribution(zpj_hist_alt_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

        # Add other samples: dijet cen
        if do_dijet_cen:
            marker = cu.Marker(shape='triangleDown')
            for sample, hist, mark in zip(self.other_samples, other_samples_dijet_central_hists, marker.cycle()):
                style_dict = sample['style_dict']
                color = style_dict.get('color', COMMON_STYLE_DICT['mc_alt_color'])
                fancy_label = style_dict.get("label", sample['key'])
                marker_size = style_dict.get('marker_size', m_size)
                cont_args = dict(label=fancy_label if only_one_region else '#splitline{ %s }{ [%s]}' % (qgc.Dijet_CEN_LABEL.replace(" region", ""), fancy_label),
                                 line_color=color,
                                 line_width=lw,
                                 line_style=style_dict.get('line_style', COMMON_STYLE_DICT['mc_alt_line_style']),
                                 marker_color=color,
                                 marker_style=style_dict.get('marker_style', mark),
                                 marker_size=marker_size,
                                 leg_draw_opt="LE" if marker_size == 0 else "LEP",
                                 subplot=dijet_central_hist_no_errors)
                entries.append(Contribution(hist, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

        # Add other samples: dijet fwd
        if do_dijet_fwd:
            marker = cu.Marker(shape='triangleDown')
            for sample, hist, mark in zip(self.other_samples, other_samples_dijet_forward_hists, marker.cycle()):
                style_dict = sample['style_dict']
                color = style_dict.get('color', COMMON_STYLE_DICT['mc_alt_color'])
                fancy_label = style_dict.get("label", sample['key'])
                marker_size = style_dict.get('marker_size', m_size)
                cont_args = dict(label=fancy_label if only_one_region else '#splitline{ %s }{ [%s]}' % (qgc.Dijet_FWD_LABEL.replace(" region", ""),fancy_label),
                                 line_color=color,
                                 line_width=lw,
                                 line_style=style_dict.get('line_style', COMMON_STYLE_DICT['mc_alt_line_style']),
                                 marker_color=color,
                                 marker_style=style_dict.get('marker_style', mark),
                                 marker_size=marker_size,
                                 leg_draw_opt="LE" if marker_size == 0 else "LEP",
                                 subplot=dijet_forward_hist_no_errors)
                entries.append(Contribution(hist, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

        # Add other samples: Z+jet
        if do_zpj:
            marker = cu.Marker(shape='triangleDown')
            for sample, hist, mark in zip(self.other_samples, other_samples_zpj_hists, marker.cycle()):
                style_dict = sample['style_dict']
                color = style_dict.get('color', COMMON_STYLE_DICT['mc_alt_color'])
                fancy_label = style_dict.get("label", sample['key'])
                marker_size = style_dict.get('marker_size', m_size)
                cont_args = dict(label=fancy_label if only_one_region else '#splitline{ %s}{ [%s]}' % (qgc.ZpJ_LABEL.replace(" region", ""), fancy_label),
                                 line_color=color,
                                 line_width=lw,
                                 line_style=style_dict.get('line_style', COMMON_STYLE_DICT['mc_alt_line_style']),
                                 marker_color=color,
                                 marker_style=style_dict.get('marker_style', mark),
                                 marker_size=marker_size,
                                 leg_draw_opt="LE" if marker_size == 0 else "LEP",
                                 subplot=zpj_hist_no_errors)
                entries.append(Contribution(hist, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

        if not self.only_yoda_data:
            if do_zpj:
                if var_number in [3,4,5,8,9,10] and metric != 'delta':
                  cont_args = dict(label="Sherpa NLO+jet",
                                   line_color=ROOT.kGray,
                                   line_width=lw,
                                   line_style=1,
                                   marker_color=ROOT.kGray,
                                   marker_style=27,
                                   marker_size=m_size,
                                   fill_color=ROOT.kGray,
                                   fill_style=0,
                                   leg_draw_opt="LEP",
                                   subplot=zpj_hist_no_errors)
                  if self.is_supplementary:                 
                    entries.append(Contribution(zpj_hist_sherpa_error, **cont_args))
                    dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

                  cont_args = dict(label="NLO + NLL'+ NP",
                                   line_color=ROOT.kRed-9,
                                   line_width=lw,
                                   line_style=1,
                                   marker_color=ROOT.kRed-9,
                                   marker_style=1,
                                   marker_size=0,
                                   fill_color=ROOT.kRed-9,
                                   fill_style=0,
                                   leg_draw_opt="L",
                                   subplot=zpj_hist_no_errors)
                  theory_style_legend = cont_args.copy()
                  theory_style_legend["fill_color"]=ROOT.kRed-9
                  theory_style_legend["fill_style"]=3395
                  theory_style_legend["leg_draw_opt"]="FL"
                  entries.append(Contribution(zpj_hist_theory, **cont_args))
                  dummy_entries.append(Contribution(dummy_gr.Clone(), **theory_style_legend))
                  
        # Add data last so it gets drawn ontop
        if metric != 'delta':
            if do_dijet_cen:
                cont_args = dict(label='Data' if only_one_region else ' '+qgc.Dijet_CEN_LABEL.replace(" region", ""),
                                 leg_draw_opt="LEP",
                                 line_color=ROOT.kBlack if only_one_region else dijet_cen_col,
                                 line_width=lw,
                                 marker_color=ROOT.kBlack if only_one_region else dijet_cen_col,
                                 marker_style=COMMON_STYLE_DICT['data_marker_style'],
                                 marker_size=m_size)
                entries.append(Contribution(dijet_central_hist, **cont_args))
                dummy_entries.insert(0, Contribution(dummy_gr.Clone(), **cont_args))

            if do_dijet_fwd:
                cont_args = dict(label='Data' if only_one_region else ' '+qgc.Dijet_FWD_LABEL.replace(" region", ""),
                                 leg_draw_opt="LEP",
                                 line_color=ROOT.kBlack if only_one_region else dijet_fwd_col,
                                 line_width=lw,
                                 marker_color=ROOT.kBlack if only_one_region else dijet_fwd_col,
                                 marker_style=COMMON_STYLE_DICT['data_marker_style'],
                                 marker_size=m_size)
                entries.append(Contribution(dijet_forward_hist, **cont_args))
                dummy_entries.insert(0, Contribution(dummy_gr.Clone(), **cont_args))

            if do_zpj:
                cont_args = dict(label='Data' if only_one_region else ' '+qgc.ZpJ_LABEL.replace(" region", ""),
                                 leg_draw_opt="LEP",
                                 line_color=ROOT.kBlack if only_one_region else zpj_col,
                                 line_width=lw,
                                 marker_color=ROOT.kBlack if only_one_region else zpj_col,
                                 marker_style=COMMON_STYLE_DICT['data_marker_style'],
                                 marker_size=m_size)
                entries.append(Contribution(zpj_hist, **cont_args))
                dummy_entries.insert(0, Contribution(dummy_gr.Clone(), **cont_args))


        # for plot axis titles
        angle_str = ""
        if metric == "mean":
            angle_str = "#LT%s#GT" % create_angle_label(angle, do_groomed)
        elif metric == "rms":
            angle_str = "RMS %s" % create_angle_label(angle, do_groomed)
        elif metric == "delta":
            angle_str = "#Delta (Simulation, Data), %s" % create_angle_label(angle, do_groomed)

        region_str = ""
        if only_one_region:
            # only 1 signal region, put it on the title
            region_str = "\n"
            if do_dijet_cen:
                region_str += qgc.Dijet_CEN_LABEL
            elif do_dijet_fwd:
                region_str += qgc.Dijet_FWD_LABEL
            elif do_zpj:
                region_str += qgc.ZpJ_LABEL

        # determine y range of main pad
        # look for global min/max across all contributions
        nbins = entries[-1].obj.GetNbinsX()

        # For the paper, we want LHA to have all the same limits for onlyYodaData and onlyDataOldMC
        is_paper = (metric == "mean" and angle.var == "jet_LHA" and jet_algo['name'] == "ak4puppi" and not do_groomed)
        use_entries = entries[-1:] if is_paper else entries
        
        h_min, min_bin_num, h_max, max_bin_num = cu.get_min_max_bin_contents_multiple_hists([c.obj for c in use_entries])
        # print('h_min, min_bin_num, h_max, max_bin_num:', h_min, min_bin_num, h_max, max_bin_num)
        h_range = h_max - h_min
        lower_padding = h_range*0.1
        upper_padding = h_range*0.1
        new_h_max = h_max+(4*upper_padding)
        new_h_min = h_min-(3*lower_padding)

        if not is_paper:
            # If the any bin overlaps with the legend region (upper 50% end of x axis)
            # then increase range so plot sits in lower half or so
            check_x_bin = int(0.75*nbins) # any bin > this is checked to see if it might overlap w/legend

            # if not ish_min_all, min_bin_num_all, h_max_all, max_bin_num_all = cu.get_min_max_bin_contents_multiple_hists([c.obj for c in entries])
            # new_h_min = h_min_all-lower_padding
            
            hits_legend = max_bin_num >= check_x_bin  # quick: check if maximum overlaps legend
            hits_legend = False
            frac_height = 0.4
            max_h_check = -1
            max_ind_check = -1
            if not hits_legend:
                # check if any bin, not just the maximum, hits the legend
                # legend is location is based on a rough fraction of the current y range
                y_limit = frac_height*(new_h_max-new_h_min) + new_h_min
                for ent in entries:
                    hist = ent.obj
                    for i in range(check_x_bin, nbins+1):
                        new_max = hist.GetBinContent(i) + hist.GetBinError(i)
                        if new_max > y_limit and new_max > max_h_check:
                            max_h_check = new_max
                            hits_legend = True
                            max_ind_check = i
                            # break

            # if metric == "mean" and angle.var == "jet_LHA_charged" and jet_algo['name'] == "ak8puppi" and not do_groomed:
            if metric == "mean" and angle.var == "jet_pTD" and jet_algo['name'] == "ak8puppi" and not do_groomed:
                print('h_min, min_bin_num, h_max, max_bin_num:', h_min, min_bin_num, h_max, max_bin_num)
                # print('h_min_all, min_bin_num_all, h_max_all, max_bin_num_all:', h_min_all, min_bin_num_all, h_max_all, max_bin_num_all)
                # print('0.6*entries[-1].obj.GetNbinsX():', 0.6*entries[-1].obj.GetNbinsX())
                print("max_h_check, max_ind_check:", max_h_check, max_ind_check)
                print("hits_legend:", hits_legend)

            if hits_legend:
                buffer_factor = 0.9
                new_h_max = (max_h_check / (frac_height*buffer_factor)) - (new_h_min * (1 - (frac_height*buffer_factor)) / (frac_height*buffer_factor))
                # print("new_h_max", new_h_max)
                # add more upper headroom, so plot essentially only covers lower half
                # new_h_max += 1*upper_padding
                # new_h_max = h_max + 0.5*upper_padding
                # new_h_max += (new_h_max-new_h_min)
            # else:
                # print("Not hit legend")
                # print('h_min, min_bin_num, h_max, max_bin_num:', h_min, min_bin_num, h_max, max_bin_num)
                # print(new_h_min, new_h_max)

        ylim = (max(0, new_h_min), new_h_max)
        plot = Plot(entries,
                    what='hist',
                    xtitle=COMMON_STYLE_DICT['jet_pt_units_str'],
                    ytitle=angle_str,
                    title="%s%s" % (jet_algo['label'], region_str),
                    # ylim=(0, h_max*1.75),
                    # ylim=(h_min*0.75, h_max*1.5),
                    ylim=ylim,
                    has_data=self.has_data,
                    is_preliminary=self.is_preliminary,
                    is_supplementary=self.is_supplementary,
                    subplot_type='ratio' if metric != 'delta' else None,
                    subplot_title=qgc.SIM_DATA_STR,
                    subplot_limits=(0.7, 1.4)# if self.has_data else (0.9, 1.1)
                    )
        plot.lumi = cu.get_lumi_str(do_dijet=any([do_dijet_fwd, do_dijet_cen]),
                                    do_zpj=do_zpj)

        # plot.default_canvas_size = (700, 600)
        # plot.title_start_y = 0.85
        plot.title_left_offset = 0.05
        plot.title_font_size = 0.035

        plot.legend.SetX1(0.55)
        plot.legend.SetX2(0.9)

        leg_y2 = plot.title_start_y + 0.01  # magic extra bit
        plot.legend.SetY2(leg_y2)
        plot.legend.SetTextSize(0.035)
        final_reading_factor=1.2
        plot.legend.SetTextSize(0.03*final_reading_factor)
        plot.legend.SetEntrySeparation(0.01)

        # scale height with number of entries
        if len(entries) < 8:
            plot.legend.SetY1(leg_y2 - (len(entries) * (0.27/6))) # magic good separation
        else:
            plot.legend.SetY1(leg_y2 - (len(entries) * (0.20/6))) # magic good separation

        # plot.legend.SetBorderSize(1)
        # plot.legend.SetLineColor(ROOT.kBlack)

        #if len(entries) > 8:
        #    plot.legend.SetNColumns(3)
        #    plot.legend.SetTextSize(0.02)
        # plot.legend.SetTextAlign(13)

        plot.left_margin = 0.16
        plot.subplot_line_style = 1
        # plot.y_padding_max_linear = 1.9

        plot.do_legend = False  # do it ourselves manually
        subplot_draw_opts = "NOSTACK E1"
        plot.plot("NOSTACK E1", subplot_draw_opts)
        # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        plot.set_logx(do_more_labels=False)

        plot.main_pad.cd()
        if do_zpj and not self.only_yoda_data and var_number in [3,4,5,8,9,10] and metric != 'delta':
          zpj_hist_theory_upper.SetFillStyle(3395)
          zpj_hist_theory_upper.SetFillColor(ROOT.kRed-9)
          zpj_hist_theory_upper.SetLineWidth(0)
          zpj_hist_theory_upper.SetMarkerSize(0)
          zpj_hist_theory_upper.Draw("F SAME")
          zpj_hist_theory_lower.SetFillStyle(1001)
          zpj_hist_theory_lower.SetFillColor(10)
          zpj_hist_theory_lower.SetLineWidth(0)
          zpj_hist_theory_lower.SetMarkerSize(0)
          zpj_hist_theory_lower.Draw("F SAME")
          ROOT.gStyle.SetNdivisions(510,"Y")
          zpj_hist_no_errors.Draw("AXIS SAME")
          plot.container.Draw("NOSTACK E1 A SAME")
        graphs=[]
        for e in plot.contributions:
          graphs+=[ROOT.TGraphErrors(e.obj)]
          graphs[-1].Draw("P SAME")

        #print("NAME: ",angle.var.replace("jet_",""),region_name,jet_algo['label'])
        #datag=plot.contributions[-1].obj
        #for b in range(datag.GetNbinsX()):
        #  print(str(datag.GetXaxis().GetBinCenter(b+1))+" "+str(datag.GetXaxis().GetBinCenter(b+1)-datag.GetXaxis().GetBinLowEdge(b+1))+" "+str(datag.GetXaxis().GetBinUpEdge(b+1)-datag.GetXaxis().GetBinCenter(b+1))+" "+str(datag.GetBinContent(b+1))+" "+str(datag.GetBinError(b+1))+" "+str(datag.GetBinError(b+1)))

        # Calculate automatic subplot limits, accounting for the range of values,
        # and allowing for the subplot legend
        #if len(plot.subplot_contributions) > 0:
        #    min_max_all = [cu.get_min_max_bin_contents(c) for c in plot.subplot_contributions]
        #    min_ratio = min([m[0] for m in min_max_all])
        #    max_ratio = max([m[2] for m in min_max_all])
        #    ratio_range = max_ratio - min_ratio
        #    # add some padding, fraction of range
        #    ratio_padding = 0.1*ratio_range
        #    min_ratio -= 2*ratio_padding  # bit more padding on lower limit as upper has legend space already
        #    max_ratio += ratio_padding
        #    new_ratio_range = max_ratio - min_ratio
        #    # now add on a fraction of this range to accommodate subplot legend
        #    new_ratio_range += 0.6*new_ratio_range
        #    plot.subplot_container.SetMinimum(min_ratio)
        #    plot.subplot_container.SetMaximum(min_ratio + new_ratio_range)

        # Do legend manually with graphs to get the right bars on the ends (sigh)
        for cont in dummy_entries:
            this_label = cont.label
            if '\n' not in this_label:
                plot.legend.AddEntry(cont.obj, this_label, cont.leg_draw_opt)
            else:
                # use #splitline instead of \n,
                # since the latter with > 1 columns will go in the wrong place
                parts = this_label.split("\n")
                if len(parts) > 2:
                    raise RuntimeError("Cannot handle > newlines in legend yet - requires nested #splitline")
                plot.legend.AddEntry(cont.obj, "#splitline{%s}{%s}" % (parts[0], parts[1]), cont.leg_draw_opt)

        plot.canvas.cd()
        plot.legend.Draw()

        if metric != 'delta':
            # now draw the data error shaded area
            # this is a bit hacky - basically draw them on the ratio pad,
            # then redraw the existing hists & line to get them ontop
            # note that we use "same" for all - this is to keep the original axes
            # (we may want to rethink this later?)
            plot.subplot_pad.cd()
            plot.subplot_legend = ROOT.TLegend(0.25, 0.75, 0.65, 0.9)
            plot.subplot_legend.SetFillStyle(0)
            plot.subplot_legend.SetTextSize(0.085)
            
            if do_zpj and not self.only_yoda_data and var_number in [3,4,5,8,9,10] and metric != 'delta':
              zpj_hist_theory_upper_ratio = zpj_hist_theory_upper.Clone()
              zpj_hist_theory_upper_ratio.Divide(zpj_hist_no_errors)
              zpj_hist_theory_upper_ratio.Draw("F SAME")
              zpj_hist_theory_lower_ratio = zpj_hist_theory_lower.Clone()
              zpj_hist_theory_lower_ratio.Divide(zpj_hist_no_errors)
              zpj_hist_theory_lower_ratio.Draw("F SAME")
              #zpj_hist_no_errors.Draw("AXIS SAME")
              plot.subplot_pad.RedrawAxis()
            
            draw_opt = "E2 SAME"
            if do_dijet_cen:
                dijet_central_hist_ratio_error.Draw(draw_opt)
                ratiohist=[dijet_central_hist_ratio_error]
                # plot.subplot_legend.AddEntry(dijet_central_hist_ratio_error, "Total data uncert%s" % (". (central)" if do_dijet_fwd else "ainty"), "F")
                plot.subplot_legend.AddEntry(dijet_central_hist_ratio_error, "Total data uncert. (central)" if do_dijet_fwd else qgc.DATA_TOTAL_UNC_STR, "F")
                if do_dijet_fwd:
                    plot.subplot_legend.SetNColumns(2)
                    plot.subplot_legend.SetX2(0.8)
            if do_dijet_fwd:
                dijet_forward_hist_ratio_error.Draw(draw_opt)
                ratiohist=[dijet_forward_hist_ratio_error]
                # plot.subplot_legend.AddEntry(dijet_forward_hist_ratio_error, "Total data uncert%s" % (". (forward)" if do_dijet_cen else "ainty"), "F")
                plot.subplot_legend.AddEntry(dijet_forward_hist_ratio_error, "Total data uncert (forward)" if do_dijet_cen else qgc.DATA_TOTAL_UNC_STR, "F")
            if do_zpj:
                # plot.subplot_legend.AddEntry(zpj_hist_ratio_error, "Total data uncertainty", "F")
                ratiohist=[zpj_hist_ratio_error]
                plot.subplot_legend.AddEntry(zpj_hist_ratio_error, qgc.DATA_TOTAL_UNC_STR, "F")
                zpj_hist_ratio_error.Draw(draw_opt)
            plot.subplot_line.Draw()
            # draw hists after line otherwise ugly overlap
            plot.subplot_container.Draw("SAME" + subplot_draw_opts)

            for e in ratiohist:
              graphs+=[ROOT.TGraphErrors(e)]
              graphs[-1].Draw("2P SAME")
            for e in plot.subplot_contributions:
              graphs+=[ROOT.TGraphErrors(e)]
              graphs[-1].Draw("P SAME")

            plot.subplot_legend.Draw()
            plot.canvas.cd()

        prefix = self._generate_filename_prefix(do_dijet_cen, do_dijet_fwd, do_zpj)
        groomed_str = '_groomed' if do_groomed else ''
        plot.canvas.Update()
        plot.save("{output_dir}/{prefix}{metric}_vs_pt_{angle_var}{groomed_str}_{algo}{append}.{fmt}"
                    .format(
                        output_dir=output_dir,
                        prefix=prefix,
                        metric=metric,
                        angle_var=angle.var,
                        groomed_str=groomed_str,
                        algo=jet_algo['name'],
                        fmt=self.output_fmt,
                        append=self.filename_append))
        #print("{output_dir}/{prefix}{metric}_vs_pt_{angle_var}{groomed_str}_{algo}{append}.{fmt}"
        #            .format(
        #                output_dir=output_dir,
        #                prefix=prefix,
        #                metric=metric,
        #                angle_var=angle.var,
        #                groomed_str=groomed_str,
        #                algo=jet_algo['name'],
        #                fmt=self.output_fmt,
        #                append=self.filename_append))

    def plot_dijet_zpj_means_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_dijet_zpj_means_vs_pt_all' % self.output_dir,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))

    def plot_dijet_means_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd, cen, fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            # self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_means_vs_pt_all' % self.output_dir, do_zpj=False, min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_dijet_cen_means_vs_pt_all' % self.output_dir,
                                                               do_zpj=False,
                                                               do_dijet_fwd=False,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_dijet_fwd_means_vs_pt_all' % self.output_dir,
                                                               do_zpj=False,
                                                               do_dijet_cen=False,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))

    def plot_zpj_means_vs_pt_all(self):
        """Plot mean vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_zpj_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            # if 'AK8' not in jet_algo['label']:
            #     continue
            # if angle.var != "jet_pTD":
            #     continue
            # if groomed:
            #     continue
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_zpj_means_vs_pt_all' % self.output_dir,
                                                               do_dijet_cen=False,
                                                               do_dijet_fwd=False,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))

    def plot_dijet_zpj_rms_vs_pt_all(self):
        """Plot RMS vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_dijet_zpj_rms_vs_pt_all' % self.output_dir,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))

    def plot_dijet_rms_vs_pt_all(self):
        """Plot RMS vs pt for dijet (cen+fwd, cen, fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            # self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_rms_vs_pt_all' % self.output_dir, do_zpj=False, min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_dijet_cen_rms_vs_pt_all' % self.output_dir,
                                                               do_zpj=False,
                                                               do_dijet_fwd=False,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_dijet_fwd_rms_vs_pt_all' % self.output_dir,
                                                               do_zpj=False,
                                                               do_dijet_cen=False,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))

    def plot_zpj_rms_vs_pt_all(self):
        """Plot RMS vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_zpj_rms_vs_pt_all' % self.output_dir,
                                                               do_dijet_cen=False,
                                                               do_dijet_fwd=False,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))

    def plot_dijet_zpj_delta_vs_pt_all(self):
        """Plot delta vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_delta_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_dijet_zpj_delta_vs_pt_all' % self.output_dir,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))

    def plot_dijet_delta_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd, cen, fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_delta_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            # self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_delta_vs_pt_all' % self.output_dir, do_zpj=False, min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_dijet_cen_delta_vs_pt_all' % self.output_dir,
                                                               do_zpj=False,
                                                               do_dijet_fwd=False,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_dijet_fwd_delta_vs_pt_all' % self.output_dir,
                                                               do_zpj=False,
                                                               do_dijet_cen=False,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))

    def plot_zpj_delta_vs_pt_all(self):
        """Plot mean vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_zpj_delta_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta',
                                                               angle,
                                                               jet_algo,
                                                               do_groomed=groomed,
                                                               output_dir='%s/plot_zpj_delta_vs_pt_all' % self.output_dir,
                                                               do_dijet_cen=False,
                                                               do_dijet_fwd=False,
                                                               min_pt_bin_ind=self.min_pt_bin_ind(jet_algo))

    @staticmethod
    def _make_hist_from_values(value_error_pairs, bins=None, title="", name="", bin_names=None):
        name = name or cu.get_unique_str()
        if bins is not None:
            nbins = len(bins)-1
            if len(value_error_pairs) != nbins:
                return ValueError("len(value) != len(bins)-1")
            h = ROOT.TH1D(name, title, bins)
        else:
            nbins = len(value_error_pairs)
            h = ROOT.TH1D(name, title, nbins, 0, nbins)
        if bin_names is not None:
            if len(bin_names) != nbins:
                raise ValueError("len(bin_names) != len(value_error_pairs)")
        for i, (v, e) in enumerate(value_error_pairs, 1):
            if bin_names is not None:
                h.GetXaxis().SetBinLabel(i, bin_names[i-1])
            h.SetBinContent(i, v)
            h.SetBinError(i, e)
        # h.SetDirectory(0)
        return h

    @staticmethod
    def _style_hist(hist,
                    line_style=None,
                    color=None,
                    fill_style=None,
                    marker_size=None,
                    marker_style=None,
                    **kwargs):
        if line_style is not None:
            hist.SetLineStyle(line_style)
        if fill_style is not None:
            hist.SetFillStyle(fill_style)
        if color is not None:
            hist.SetLineColor(color)
            hist.SetFillColor(color)
            hist.SetMarkerColor(color)
        if marker_size is not None:
            hist.SetMarkerSize(marker_size)
        if marker_style is not None:
            hist.SetMarkerStyle(marker_style)

    def _style_data_hist(self, hist):
        self._style_hist(hist,
                         color=COMMON_STYLE_DICT['data_color'],
                         line_style=COMMON_STYLE_DICT['data_line_style'],
                         marker_size=COMMON_STYLE_DICT['marker_size'],
                         marker_style=COMMON_STYLE_DICT['data_marker_style'])

    def _style_mc_hist(self, hist):
        self._style_hist(hist,
                         line_style=COMMON_STYLE_DICT['mc_line_style'],
                         color=COMMON_STYLE_DICT['mc_color'],
                         fill_style=COMMON_STYLE_DICT['mc_fill_style'],
                         marker_size=COMMON_STYLE_DICT['marker_size'],
                         marker_style=COMMON_STYLE_DICT['mc_marker_style'])

    def _style_alt_mc_hist(self, hist):
        self._style_hist(hist,
                         line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                         color=COMMON_STYLE_DICT['mc_alt_color'],
                         fill_style=COMMON_STYLE_DICT['mc_alt_fill_style'],
                         marker_size=COMMON_STYLE_DICT['marker_size'],
                         marker_style=COMMON_STYLE_DICT['mc_alt_marker_style'])

    def _style_theory_hist(self, hist):
        self._style_hist(hist,
                         line_style=1,
                         color=ROOT.kRed-9,
                         fill_style=3395,
                         fill_color=ROOT.kRed-9,
                         marker_size=0,
                         marker_style=1)

    @staticmethod
    def calc_hists_max_min(hists):
        y_up, y_down = -999999, 999999
        for h in hists:
            y_up = max(max([h.GetBinContent(i) + h.GetBinError(i) for i in range(1, h.GetNbinsX()+1)]), y_up)
            y_down = min(min([h.GetBinContent(i) - h.GetBinError(i) for i in range(1, h.GetNbinsX()+1)]), y_down)
        return y_up, y_down

    def calc_auto_ylim(self, hists, up_padding_frac=0.2, down_padding_frac=0.2):
        y_up, y_down = self.calc_hists_max_min(hists)
        y_range = y_up - y_down
        down_padding = down_padding_frac * y_range
        up_padding = up_padding_frac * y_range
        return y_up + up_padding, y_down - down_padding

    def construct_mean_rms_hist_groups(self, selections, add_theory=True):
        """Construct mean & RMS hists for plots

        See plot_mean_rms_bins_summary() docstring about `selections` arg.

        Each of the mean/RMS returns are lists. Each item corresponds to the
        respective item in the outer level of the `selections` arg.
        (e.g. each angle).
        Each of the hist list items is itself a list, corresponding to the hists
        for [data, nominal MC, alt MC, other samples, ...].

        Also does styling of hists (as easier here)
        """
        mean_hists, rms_hists = [], []
        for selection_group in selections:

            mean_entries_data = []
            mean_entries_mc = []
            mean_entries_alt_mc = []
            mean_entries_theory = []
            # each top-level index is a sample
            mean_entries_other_samples = [[] for _ in self.other_samples]  # for loop needed to create indepedent lists

            rms_entries_data = []
            rms_entries_mc = []
            rms_entries_alt_mc = []
            rms_entries_other_samples = [[] for _ in self.other_samples]

            bin_names = []

            for query, label in selection_group['selections']:
                # print(query)
                results = self.df.query(query)
                if len(results.index) != 1:
                    print(query)
                    print(results)
                    raise ValueError("Got != 1 results, check query")

                mean_entries_data.append([results['mean'].item(), results['mean_err'].item()])

                rms_entries_data.append([results['rms'].item(), results['rms_err'].item()])

                if not self.only_yoda_data:
                    mean_entries_mc.append([results['mean_truth'].item(), results['mean_err_truth'].item()])
                    mean_entries_alt_mc.append([results['mean_alt_truth'].item(), results['mean_err_alt_truth'].item()])
                    ### new theory
                    var_number=["jet_puppiMultiplicity", "jet_pTD","jet_LHA", "jet_width", "jet_thrust", 
                            "jet_puppiMultiplicity_charged", "jet_pTD_charged", "jet_LHA_charged", "jet_width_charged", "jet_thrust_charged"].index(query.split('"')[-2])+1
                    ibin=int(query.split('&')[1].split("=")[-1])
                    if "ZPlusJets" in query and var_number in [3,4,5,8,9,10]:
                     angle_output_dir = "%s/%s/%s" % ("/nfs/dust/cms/user/hinzmann/qganalysis/CMSSW_10_2_17/src/trees/unfolded_ak"+query.split('"')[1][2]+"puppi", query.split('"')[-4], query.split('"')[-2])
                     root_filename = os.path.join(angle_output_dir, "unfolding_result_slim.root")
                     import uproot
                     uproot_file = uproot.open(root_filename)
                     #print(root_filename, query.split('"')[-4], query.split('"')[-2], self.pt_bins_zpj[(0 if "4" in query.split('"')[1] else 2):])
                     unfolding_dict = unpack_slim_unfolding_root_file_uproot(uproot_file, query.split('"')[-4], query.split('"')[-2], self.pt_bins_zpj[(0 if "4" in query.split('"')[1] else 2):])
                     errors=np.ndarray(1)
                     areas=np.ndarray(1)
                     centers=np.ndarray(1)
                     errors=np.delete(errors,0)
                     areas=np.delete(areas,0)
                     centers=np.delete(centers,0)
                     for post in ["_fix",""]:
                       theory_file="/nfs/dust/cms/user/hinzmann/qganalysis/CMSSW_10_2_17/src/zjet_angularities-master-to_send_to_CMS-thr_vs_exp_cvt"+post+"/to_send_to_CMS/thr_vs_exp_cvt"+post+"/"
                       theory_file+="data_R"+query.split('"')[1][2]+"/"
                       theory_file+="RSIG_THR_RES_NP_"
                       theory_file+="d"+("1" if "4" in query.split('"')[1] else "2")+("2" if (not "~isgroomed" in query) else "1")+"-"
                       theory_file+="x"+("0" if var_number<10 else "")+str(var_number)+"-"
                       theory_file+="y"+("0" if ibin<9 else "")+str(ibin+1)
                       theory_file+=".dat"
                       if os.path.exists(theory_file): break
                     #print(theory_file)
                     igenbin=0
                     for line in open(theory_file).readlines():
                       if line.count(".")==5:
                         igenbin+=1
                         width=unfolding_dict['truth_hists'][ibin].edges[igenbin]-unfolding_dict['truth_hists'][ibin].edges[igenbin-1]
                         errors=np.append(errors,float(line.split("\t")[3])*width)
                         areas=np.append(areas,float(line.split("\t")[2])*width)
                         centers=np.append(centers,unfolding_dict['truth_hists'][ibin].edges[igenbin-1]+width/2.)
                         assert(float(line.split("\t")[0])==unfolding_dict['truth_hists'][ibin].edges[igenbin-1])
                     mean = metrics.calc_mean_jax(areas, centers)
                     err = metrics.calc_mean_uncorrelated_error_jax(areas, centers, errors)
                     mean_entries_theory.append([mean,err])
                    else:
                     mean_entries_theory.append([-999,0])

                    rms_entries_mc.append([results['rms_truth'].item(), results['rms_err_truth'].item()])
                    rms_entries_alt_mc.append([results['rms_alt_truth'].item(), results['rms_err_alt_truth'].item()])

                for ind, sample in enumerate(self.other_samples):
                    key = sample['key']
                    mean_entries_other_samples[ind].append([results['mean_%s' % key].item(), results['mean_err_%s' % key].item()])
                    rms_entries_other_samples[ind].append([results['rms_%s' % key].item(), results['rms_err_%s' % key].item()])

                bin_names.append(label)

            # check for dodgy values
            def _check_values(values, label):
                arr = np.array(values)
                if np.any(arr == 0):
                    print("0 in", label)
                if np.all(arr == 0):
                    print("0 for all", label)
                if np.any(arr < 0):
                    print("<0 in", label)
                if np.all(arr < 0):
                    print("<0 for all", label)
                if np.isnan(arr).any():
                    print("nan in", label)
                    print(arr)
                if np.isinf(arr).any():
                    print("inf in", label)
                    print(arr)

            _check_values(mean_entries_data, "mean_entries_data")
            # set bin names here, easier
            hist_mean_data = self._make_hist_from_values(mean_entries_data, bin_names=bin_names)
            self._style_data_hist(hist_mean_data)

            hists = [hist_mean_data]

            if not self.only_yoda_data:
                _check_values(mean_entries_mc, "mean_entries_mc")
                hist_mean_mc = self._make_hist_from_values(mean_entries_mc, bin_names=bin_names)
                self._style_mc_hist(hist_mean_mc)
                hists.append(hist_mean_mc)

                _check_values(mean_entries_alt_mc, "mean_entries_alt_mc")
                hist_mean_alt_mc = self._make_hist_from_values(mean_entries_alt_mc, bin_names=bin_names)
                self._style_alt_mc_hist(hist_mean_alt_mc)
                hists.append(hist_mean_alt_mc)

                _check_values(mean_entries_theory, "mean_entries_theory")
                hist_mean_theory = self._make_hist_from_values(mean_entries_theory, bin_names=bin_names)
                self._style_theory_hist(hist_mean_theory)
                if add_theory:
                  hists.append(hist_mean_theory)

            for sample, entries in zip(self.other_samples, mean_entries_other_samples):
                _check_values(entries, "mean_entries_%s" % sample['key'])
                hist = self._make_hist_from_values(entries, bin_names=bin_names)
                self._style_hist(hist, **sample['style_dict'])
                hists.append(hist)

            mean_hists.append(hists)

            _check_values(rms_entries_data, "rms_entries_data")
            hist_rms_data = self._make_hist_from_values(rms_entries_data, bin_names=bin_names)
            self._style_data_hist(hist_rms_data)

            hists = [hist_rms_data]

            if not self.only_yoda_data:
                _check_values(rms_entries_mc, "rms_entries_mc")
                hist_rms_mc = self._make_hist_from_values(rms_entries_mc, bin_names=bin_names)
                self._style_mc_hist(hist_rms_mc)
                hists.append(hist_rms_mc)

                _check_values(rms_entries_alt_mc, "rms_entries_alt_mc")
                hist_rms_alt_mc = self._make_hist_from_values(rms_entries_alt_mc, bin_names=bin_names)
                self._style_alt_mc_hist(hist_rms_alt_mc)
                hists.append(hist_rms_alt_mc)

            for sample, entries in zip(self.other_samples, rms_entries_other_samples):
                _check_values(entries, "rms_entries_%s" % sample['key'])
                hist = self._make_hist_from_values(entries, bin_names=bin_names)
                self._style_hist(hist, **sample['style_dict'])
                hists.append(hist)

            rms_hists.append(hists)

        return mean_hists, rms_hists

    def plot_mean_rms_bins_summary(self, selections, output_file, legend_header=None, ylims_upper=None, ylims_lower=None):
        """Make plot of mean & RMS for various selections, showing data, mc, alt mc, etc """
        print("Plotting mean_rms_bins_summary for", legend_header)
        # Get data for plots
        mean_hists, rms_hists = self.construct_mean_rms_hist_groups(selections)

        self.plot_two_row_bins_summary(selection_groups=selections,
                                       upper_row_hist_groups=mean_hists,
                                       lower_row_hist_groups=rms_hists,
                                       output_file=output_file,
                                       upper_row_label="#LT#lambda^{#kappa}_{#beta}#GT",
                                       lower_row_label="RMS #lambda^{#kappa}_{#beta}",
                                       legend_header=legend_header,
                                       label_every_bin=False,
                                       ylims_upper=ylims_upper,
                                       ylims_lower=ylims_lower)

    def construct_delta_hist_groups(self, selections):
        """Construct delta hists for plots (nominal MC, & alt MC)

        See plot_delta_bins_summary() docstring about `selections` arg.

        Returns a list. Each item corresponds to the respective item in the
        outer level of the `selections` arg. (e.g. each angle).

        Each of the hist list items is itself a list, corresponding to the hists
        for [nominal MC, alt MC].

        Also does styling of hists (as easier here)
        """
        delta_hists = []
        for selection_group in selections:

            entries_nominal = []
            entries_alt = []

            bin_names = []

            for query, label in selection_group['selections']:
                # print(query)
                results = self.df.query(query)
                if len(results.index) != 1:
                    print(query)
                    print(results)
                    raise ValueError("Got != 1 results, check query")

                # if any([np.isnan(v) for v in [results['delta_truth'], results['delta_err_truth']]]):
                #     raise ValueError("Got a nan in delta_truth, delta_err_truth")
                # if any([np.isnan(v) for v in [results['delta_alt_truth'], results['delta_err_alt_truth']]]):
                #     raise ValueError("Got a nan in delta_truth, delta_err_truth")
                entries_nominal.append([results['delta_truth'], results['delta_err_truth']])
                entries_alt.append([results['delta_alt_truth'], results['delta_err_alt_truth']])

                # FIXME: add other_samples

                if "\n" in label:
                    parts = label.split("\n")
                    if len(parts)>2:
                        raise RuntimeError("too many \\n")
                    label = "#splitline{%s}{%s}" % (parts[0], parts[1])
                bin_names.append(label)

            hist_delta_nominal = self._make_hist_from_values(entries_nominal, bin_names=bin_names)
            self._style_mc_hist(hist_delta_nominal)

            hist_delta_alt = self._make_hist_from_values(entries_alt, bin_names=bin_names)
            self._style_alt_mc_hist(hist_delta_alt)

            delta_hists.append([hist_delta_nominal, hist_delta_alt])

        return delta_hists


    def plot_delta_bins_summary(self, selections, output_file, legend_header=None):
        """"""
        print("Plotting delta_bins_summary for", legend_header)
        # Get data for plots
        delta_hists = self.construct_delta_hist_groups(selections)

        gc_stash = [] # to stop stuff being deleted

        # Setup canvas and pads
        # Setup same dimension etc as mean & RMS plot?
        canvas = ROOT.TCanvas("c_"+cu.get_unique_str(), "", 1200, 400)
        canvas.SetBottomMargin(0.0)
        canvas.SetTopMargin(0.0)
        canvas.SetLeftMargin(0.0)
        canvas.SetRightMargin(0.0)
        canvas.cd()

        # For each selection group, create a pad
        pads = []
        n_pads = len(selections)

        # gap between right end of plots and edge of canvas, used for legend
        right_margin = 0.15
        # pad_left_titles_gap = 0.01 # gap between pad_left_titles and all plots
        pad_to_pad_gap = 0.005  # gap between plot pad columns
        # how far in from the left the first plotting pad starts. used for y axis title
        left_margin = 0.04
        # figure out width per pad - get total width available, then divide by number of pads
        pad_width = (1 - left_margin - right_margin - (pad_to_pad_gap*(n_pads - 1))) / n_pads

        pad_offset_bottom = 0.01 # spacing between bottom of pad and canvas edge
        pad_offset_top = 0.08  # spacing between top of pad and canvas edge - for CMS and lumi text

        # per-pad margins: these determine where the hist axes lie,
        # and are fractions of the **pad** width/height, not the global canvas
        pad_right_margin = 0.02
        pad_left_margin = 0.25

        # bottom margin includes space for x axis labels
        pad_bottom_margin = 0.49

        # extra bit to add to the top margins of lower and upper pads
        # to ensure y axis numbers don't get cut off
        pad_top_margin = 0.018

        for isel, selection_group in enumerate(selections):
            canvas.cd()
            pad_start_x = left_margin + (isel*pad_to_pad_gap) + (isel*pad_width)
            pad_end_x = pad_start_x + pad_width
            # Create pad
            pad = ROOT.TPad(cu.get_unique_str(), "", pad_start_x, pad_offset_bottom, pad_end_x, 1-pad_offset_top)
            ROOT.SetOwnership(pad, False)
            pad.SetFillColor(isel+2)
            pad.SetFillStyle(3004)
            pad.SetFillStyle(4000)
            pad.SetTopMargin(pad_top_margin)
            pad.SetBottomMargin(pad_bottom_margin)
            pad.SetRightMargin(pad_right_margin)
            pad.SetLeftMargin(pad_left_margin)
            pad.SetTicks(1, 1)
            pad.Draw()
            pads.append(pad)

        # Now draw the histograms
        for isel, (selection_group, pad, hist_group) in enumerate(zip(selections, pads, delta_hists)):
            pad.cd()
            # FIXME add other_samples
            hist_group[0].Draw("E1")
            hist_group[1].Draw("E1 SAME")

            draw_hist = hist_group[0]

            xax = draw_hist.GetXaxis()
            xax.CenterLabels()
            xax.LabelsOption("v")

            factor = 1.5
            label_size_fudge = 1.6
            label_offset_fudge = 4
            tick_fudge = 4

            xax.SetTickLength(xax.GetTickLength()*factor)
            xax.SetLabelSize(xax.GetLabelSize()*factor*2)
            xax.SetLabelOffset(xax.GetLabelOffset()*factor)

            yax = draw_hist.GetYaxis()
            yax.SetLabelSize(yax.GetLabelSize()*factor*label_size_fudge)
            yax.SetLabelOffset(yax.GetLabelOffset()*factor*label_offset_fudge)
            yax.SetTickLength(yax.GetTickLength()*factor*tick_fudge)  # extra bit of fudging, probably becasue pads are sligthly different sizes
            primary = 4 # 4 is the magic primary number, not 3
            secondary = 5
            n_divisions_opts = [primary, secondary, 0, True]
            yax.SetNdivisions(*n_divisions_opts)

            # Set range using bin contents + error bars
            y_up, y_down = self.calc_hists_max_min(hist_group)
            y_range = y_up - y_down
            down_padding = 0.2 * y_range
            up_padding = 0.22 * y_range
            draw_hist.SetMinimum(y_down - down_padding)
            draw_hist.SetMaximum(y_up + up_padding)

            # Draw variable name
            var_latex = ROOT.TLatex()
            var_latex.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignBottom)
            var_latex.SetTextSize(0.1)
            # var_latex.SetTextFont(42)
            # these are relative to the RMS pad! not the canvas
            # put in middle of the plot (if 0.5, woud look off-centre)
            var_x = 0.5*(1-pad_right_margin-pad_left_margin) + pad_left_margin
            var_y = 0.03  # do by eye
            var_latex.DrawLatexNDC(var_x, var_y, selection_group['label'])
            gc_stash.append(var_latex)

        canvas.cd()

        # Add legend
        leg_y_top = 1 - pad_offset_top
        leg_left = pads[-1].GetAbsXlowNDC() + pads[-1].GetAbsWNDC() + pad_to_pad_gap
        leg_right = 1-0.02
        leg_y_bottom = leg_y_top-pads[0].GetAbsHNDC()+(pad_bottom_margin*pads[0].GetAbsHNDC()) # want the pad to correspond to the plotable area
        leg_pad = ROOT.TPad("leg_pad_"+cu.get_unique_str(), "", leg_left, leg_y_bottom, leg_right, leg_y_top)
        ROOT.SetOwnership(leg_pad, False)  # important! otherwise seg fault
        leg_pad.SetFillStyle(4000)
        # leg_pad.SetFillColor(ROOT.kYellow)
        # leg_pad.SetFillStyle(3004)
        leg_pad.SetLeftMargin(0)
        leg_pad.SetRightMargin(0)
        leg_pad.SetTopMargin(0)
        leg_pad.SetBottomMargin(0)
        leg_pad.Draw()
        leg_pad.cd()
        gc_stash.append(leg_pad)
        leg = ROOT.TLegend(0., pads[0].GetBottomMargin(), 1, 1)

        pt = None
        if legend_header:
            # Add title to legend
            # Add ability to do multiple lines by splitting on \n
            # Assumes first line most important, so bolded
            # Dont account for blank lines
            num_header_lines = len([x for x in legend_header.split("\n") if len(x) > 0])
            line_height = 0.15
            offset = num_header_lines * line_height
            # move legend down by the height of the new TPaveText
            leg.SetY1(leg.GetY1()-offset)
            leg.SetY2(leg.GetY2()-offset)
            pt = ROOT.TPaveText(leg.GetX1(), leg.GetY2(), leg.GetX2(), leg_y_top, "NDC NB")
            pt.SetFillStyle(0)
            pt.SetBorderSize(0)
            for line_ind, line in enumerate(legend_header.split("\n")):
                text = pt.AddText(line)
                text.SetTextAlign(11)
                if line_ind == 0:
                    text.SetTextFont(62)
                    text.SetTextSize(0.1)
                else:
                    text.SetTextFont(42)
                    text.SetTextSize(0.09)
            pt.Draw()

        # Replace legend markers with graph to get correct error bar endings
        # Yes this is ridiculous
        dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        dummy_mc = dummy_gr.Clone()
        self._style_mc_hist(dummy_mc)
        dummy_alt_mc = dummy_gr.Clone()
        self._style_alt_mc_hist(dummy_alt_mc)
        leg.AddEntry(dummy_mc, self.mc_label, "EL")
        leg.AddEntry(dummy_alt_mc, self.alt_mc_label, "EL")

        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        # leg.SetFillColor(ROOT.kGreen)
        # leg.SetFillStyle(3004)
        leg.SetTextSize(0.1)
        leg.SetTextAlign(12)
        leg.SetEntrySeparation(0.08)
        leg.Draw()

        canvas.cd()

        # Add delta text
        text_width = left_margin
        text_x = (0.5 * left_margin) - (0.5*text_width)
        # let ROOT center it, just make the box the height of the axis
        text_y_end = 1 - pad_offset_top - (pad_top_margin*pads[0].GetAbsHNDC())
        text_y = 1 - pad_offset_top - pads[0].GetAbsHNDC() + (pads[0].GetAbsHNDC()*pad_bottom_margin)
        delta_text = ROOT.TPaveText(text_x, text_y, text_x+text_width, text_y_end, "NDC NB")
        delta_text.SetFillStyle(0)
        delta_text.SetBorderSize(0)
        t_delta = delta_text.AddText("#Delta (Simulation, Data)")
        t_delta.SetTextAngle(90)
        text_size = 0.06
        t_delta.SetTextSize(text_size)
        delta_text.Draw()

        # Add CMS text
        cms_latex = ROOT.TLatex()
        cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        cms_latex.SetTextFont(42)
        extraOverCmsTextSize  = 0.76
        cms_latex.SetTextSize(0.05/extraOverCmsTextSize)
        # Get the text sitting just above the axes of the mean plot
        # Axes end inside the mean pad at (1-top_margin), but this has
        # to be scaled to canvas NDC
        # Then add a little extra spacing ontop to separate text from axes line
        latex_height = 1 - pad_offset_top - (pads[0].GetAbsHNDC() * pads[0].GetTopMargin()) + 0.02
        # Want it to start at the left edge of the first plot
        start_x = left_margin + (pad_width*pad_left_margin)
        cms_latex.DrawLatexNDC(start_x, latex_height, "#font[61]{CMS}")

        cms_latex = ROOT.TLatex()
        cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        cms_latex.SetTextFont(42)
        cms_latex.SetTextSize(0.05)
        # Get the text sitting just above the axes of the mean plot
        # Axes end inside the mean pad at (1-top_margin), but this has
        # to be scaled to canvas NDC
        # Then add a little extra spacing ontop to separate text from axes line
        latex_height = 1 - pad_offset_top - (pads[0].GetAbsHNDC() * pads[0].GetTopMargin()) + 0.02

        # Want it to start at the left edge of the first plot
        start_x = left_margin + (pad_width*pad_left_margin) + 0.03
        # # start_x = 100
        is_supplementary = True
        if self.is_preliminary:
            if self.has_data:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[52]{ Preliminary}")
            else:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[52]{ Simulation Preliminary}")
        else:
            if is_supplementary:
                if self.has_data:
                    cms_latex.DrawLatexNDC(start_x, latex_height, "#font[52]{ Supplementary}")
                else:
                    cms_latex.DrawLatexNDC(start_x, latex_height, "#font[52]{ Simulation Supplementary}")
            else:
                if self.has_data:
                    cms_latex.DrawLatexNDC(start_x, latex_height, "")
                else:
                    cms_latex.DrawLatexNDC(start_x, latex_height, "#font[52]{ Simulation}")
        cms_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
        # Get the lumi text aligned to right edge of axes
        # i.e. 1-pad_right_margin, but remember to scale by pad width
        end_x = 1 - right_margin - (pads[0].GetAbsWNDC() * pad_right_margin)
        # end_x = 0.985  # to match last pad
        # Figure out why physical region(s) are in plot, to get correct lumi
        do_zpj, do_dijet = False, False
        for sg in selections:
            if any(["ZPlusJets" in x[0] for x in sg['selections']]):
                do_zpj = True
            if any(["Dijet" in x[0] for x in sg['selections']]):
                do_dijet = True
        lumi = cu.get_lumi_str(do_dijet=do_dijet, do_zpj=do_zpj)
        cms_latex.DrawLatexNDC(end_x, latex_height, " %s fb^{-1} (13 TeV)" % lumi)
        gc_stash.append(cms_latex)

        canvas.Update()
        cu.check_dir_exists_create(os.path.dirname(os.path.abspath(output_file)))
        canvas.SaveAs(output_file)


    def construct_q_vs_g_hist_groups(self, gluon_selections, quark_selections, add_theory=True):
        """Summary

        Parameters
        ----------
        gluon_selections : TYPE
            Description
        quark_selections : TYPE
            Description

        Returns
        -------
        TYPE
            Description

        Raises
        ------
        ValueError
            Description
        """
        if len(gluon_selections) != len(quark_selections):
            raise ValueError("Require len(gluon_selections) == len(quark_selections)")

        gluon_mean_hists, gluon_rms_hists = self.construct_mean_rms_hist_groups(gluon_selections, add_theory)
        quark_mean_hists, quark_rms_hists = self.construct_mean_rms_hist_groups(quark_selections, add_theory)

        mean_hists, rms_hists = [], []
        for sel_ind in range(len(gluon_mean_hists)):
            # for each selection group (i.e. angle),
            # we have a list of hists [data, mc, alt mc, sample A, ...]
            # We need to construct copies with gluon / quark

            # Mean
            this_gluon_hists = gluon_mean_hists[sel_ind]
            this_quark_hists = quark_mean_hists[sel_ind]
            this_ratio_hists = []
            for g_hist, q_hist in zip(this_gluon_hists, this_quark_hists):
                ratio_hist = g_hist.Clone(g_hist.GetName() + "_qg_ratio")
                ratio_hist.Divide(q_hist)
                this_ratio_hists.append(ratio_hist)
                # no need to worry about styling, carried over from originals
            mean_hists.append(this_ratio_hists)

            # RMS
            this_gluon_hists = gluon_rms_hists[sel_ind]
            this_quark_hists = quark_rms_hists[sel_ind]
            this_ratio_hists = []
            for g_hist, q_hist in zip(this_gluon_hists, this_quark_hists):
                ratio_hist = g_hist.Clone(g_hist.GetName() + "_qg_ratio")
                ratio_hist.Divide(q_hist)
                this_ratio_hists.append(ratio_hist)
                # no need to worry about styling, carried over from originals
            rms_hists.append(this_ratio_hists)
        return mean_hists, rms_hists

    def construct_q_vs_g_mc_vs_data_hist_groups(self, ratio_hist_groups, data_index=0):
        """Summary

        Parameters
        ----------
        ratio_hist_groups : TYPE
            Description
        data_index : int, optional
            Description
        """
        new_hists = []
        for hist_group in ratio_hist_groups:
            this_group = []
            ref_hist = hist_group[data_index]
            for ind, h in enumerate(hist_group):
                if ind == data_index:
                    continue
                h_new = h.Clone(h.GetName() + "_vs_data")
                h_new.Divide(ref_hist)
                this_group.append(h_new)
            new_hists.append(this_group)
        return new_hists

    def plot_q_vs_g_bins_summary(self, quark_selections, gluon_selections, output_file, legend_header=None, ylims_upper=None, ylims_lower=None):
        """Plot gluon/quark summary stats for various selections, showing data, MCs.

        Also show MC / Data for the various simulations.
        """
        print("plotting q_vs_g_bins_summary")
        mean_q_vs_g_ratio_hists, rms_q_vs_g_ratio_hists = self.construct_q_vs_g_hist_groups(gluon_selections=gluon_selections,
                                                                                            quark_selections=quark_selections, add_theory=False)
        mean_q_vs_g_data_vs_mc_ratio_hists = self.construct_q_vs_g_mc_vs_data_hist_groups(mean_q_vs_g_ratio_hists)
        self.plot_two_row_bins_summary(selection_groups=quark_selections[0:1]+gluon_selections[1:], # merge the two to get correct identification of doing dijet and Z+J e.g. for lumi
                                       upper_row_hist_groups=mean_q_vs_g_ratio_hists,
                                       lower_row_hist_groups=mean_q_vs_g_data_vs_mc_ratio_hists,
                                       # upper_row_label="#splitline{g-enriched mean /}{q-enriched mean}",  # don't use #frac, numerator too close to line, and kerning doesn't work
                                       # upper_row_label="#frac{g-enriched mean}{q-enriched mean}",  # don't use #frac, numerator too close to line, and kerning doesn't work
                                       # upper_row_label="#frac{g-enriched #LT #lambda #GT}{q-enriched #LT #lambda #GT}",  # don't use #frac, numerator too close to line, and kerning doesn't work
                                       # upper_row_label="#splitline{g-enriched #LT #lambda #GT /}{q-enriched #LT #lambda #GT}",  # don't use #frac, numerator too close to line, and kerning doesn't work
                                       # upper_row_label="#splitline{g-enriched #LT #lambda^{#kappa}_{#beta} #GT /}{q-enriched #LT #lambda^{#kappa}_{#beta} #GT}",  # don't use #frac, numerator too close to line, and kerning doesn't work
                                       upper_row_label="#frac{g-enriched #LT#lambda^{#kappa}_{#beta}#GT}{#lower[0.08]{q-enriched #LT#lambda^{#kappa}_{#beta}#GT}}",  # don't use #frac, numerator too close to line, and kerning doesn't work
                                       # upper_row_label="#splitline{#LT g-enriched #GT /}{#LT q-enriched #GT}",
                                       # upper_row_label="#frac{#lower[-0.3]{#LT g-enriched #GT}}{#lower[0.15]{#LT q-enriched #GT}}",
                                       # upper_row_label="#lower[0.2]{#frac{#lower[0.]{#LT g-enriched #GT}}{#lower[0.15]{#LT q-enriched #GT}}}",
                                       # upper_row_label="#lower[0.2]{#frac{#lower[0.]{#LT g-enriched #GT}}{#lower[0.15]{#LT q-enriched #GT}}}",
                                       # lower_row_label="#splitline{Simulation /}{     Data}",
                                       lower_row_label="#frac{#lower[-0.08]{Simulation}}{#lower[0.15]{Data}}",
                                       output_file=output_file,
                                       legend_header=legend_header,
                                       label_every_bin=False,
                                       lower_row_is_ratio=True,
                                       ylims_upper=ylims_upper,
                                       ylims_lower=ylims_lower, add_theory=False)


    def plot_q_g_mean_bins_summary(self, quark_selections, gluon_selections, output_file, legend_header=None, ylims=None):
        """Plot means for quark & gluon separately for various selections, showing data, MCs

        Like plot_mean_rms_bins_summary, but instead of RMS, bottom row is also mean,
        but with different selection
        """
        print("Plot q g mean summary")

        if len(gluon_selections) != len(quark_selections):
            raise ValueError("Require len(gluon_selections) == len(quark_selections)")

        gluon_mean_hists, _ = self.construct_mean_rms_hist_groups(gluon_selections)
        quark_mean_hists, _ = self.construct_mean_rms_hist_groups(quark_selections)
        self.plot_two_row_bins_summary(selection_groups=quark_selections[0:1]+gluon_selections[1:], # merge the two to get correct identification of doing dijet and Z+J e.g. for lumi
                                       upper_row_hist_groups=gluon_mean_hists,
                                       lower_row_hist_groups=quark_mean_hists,
                                       # upper_row_label="#splitline{Gluon-enriched}{         mean}",  # manually do spacing, since ROOT only left-aligns splitline
                                       # lower_row_label="#splitline{Quark-enriched}{         mean}",
                                       # upper_row_label=" #LT Gluon-enriched #GT",  # manually do spacing, since ROOT only left-aligns splitline
                                       # lower_row_label="#LT Quark-enriched #GT",
                                       upper_row_label=" g-enriched #LT#lambda^{#kappa}_{#beta}#GT",  # manually do spacing, since ROOT only left-aligns splitline
                                       lower_row_label="q-enriched #LT#lambda^{#kappa}_{#beta}#GT",
                                       output_file=output_file,
                                       legend_header=legend_header,
                                       label_every_bin=False,
                                       upper_lower_same_ylim=True,
                                       ylims_upper=ylims,
                                       ylims_lower=ylims)


    @staticmethod
    def construct_data_ratio_hists(data_hist, other_hists, data_is_denominator=True):
        """Construct ratio histograms of other/data (if data_is_denominator),
        or data/other

        Parameters
        ----------
        data_hist : TH1
            Description
        other_hists : list[TH1]
            Description
        data_is_denominator : bool, optional
            If true, do other/data. Otherwise data/other
        """
        new_hists = []
        for h in other_hists:
            h_new = h.Clone()
            if data_is_denominator:
                h_new.Divide(data_hist)
            else:
                h_new.Divide(data_hist, h)
            new_hists.append(h_new)
        return new_hists

    def plot_two_row_bins_summary(self,
                                  selection_groups,
                                  upper_row_hist_groups,
                                  lower_row_hist_groups,
                                  output_file,
                                  legend_header=None,
                                  upper_row_label="",
                                  lower_row_label="",
                                  label_every_bin=True,
                                  lower_row_is_ratio=False,
                                  upper_lower_same_ylim=False,
                                  ylims_upper=None,
                                  ylims_lower=None, add_theory=True):
        """Plot 2 row summary plot showing values from choice bins

        `selection_groups` is a multi-level list
        The first level is of dicts,
            {'label': str, 'selections': list}
        where each dict represents an angle (i.e. its own column).

        The 'selections' key of that dict then provides a list of tuples,
        where each is a bin in the histograms that make up each plot column.

        Each tuple is of the form (query_str, bin_label_str),
        where the query_str gets passed to the dataframe to retrieve a single entry.
        If 0 or >1 entries are found, this raises a ValueError.
        From that entry, we get mean ± error, and RMS ± error.
        That bin in the histogram is then labelled with bin_label_str.

        In this way, we build up the contents of the histograms for each angle.


        The *_row_hist_groups are a bit complicated: each is a list, where
        each entry corresponds to a column group.
        Within that list are the different sample entries, e.g. [data, MC, ...]
        Data is assumed to be first.

        The plotting part is complicated, since we need 2 pads (upper & lower)
        for each column. These need to overlap, since otherwise it will chop
        off the y axes labels. However, we want the top x axis of the lower pad,
        and the lower x axis of the upper pad, to align vertically. So we have
        to be careful with the exact pad & margin sizes.

        Note that the margins basically determine where the axes go in the pad,
        and are a fraction of the pad width/height. Axis labels must go in the
        pad margin - nothing can be drawn outside the pad itself.

        We also manually put in the y axis titles, and the angle names, as
        individual text elements, just so we can control their position exactly.

        Parameters
        ----------
        selection_groups : TYPE
            Description
        upper_row_hist_groups : list[list[ROOT.TH1]]
            TH1s to plot on upper row
        lower_row_hist_groups : list[list[ROOT.TH1]]
            TH1s to plot on lower row
        output_file : str
            Output filename
        legend_header : str, optional
            Optional title to put on legend
        upper_row_label : str, optional
            Y label for upper row
        lower_row_label : str, optional
            Y label for lower row
        label_every_bin : bool, optional
            Description
        lower_row_is_ratio : bool, optional
            True if lower row should be made into ratio plot of simulation / data
            Where data is the first entry in each hist group
        upper_lower_same_ylim : bool, optional
            Enforce same ylimit for upper and lower plots, per column
        ylims_upper : list[(float, float)], optional
            List of y axis limits to use on upper row
            One entry per column (i.e. variable), of a tuple of (lower, upper)
        ylims_lower : list[(float, float)], optional
            List of y axis limits to use on lower row
            Ignored if upper_lower_same_ylim = True and ylims_upper is set

        TODO: move into own class to be more customisable?

        Raises
        ------
        RuntimeError
            Description
        ValueError
            Description
        """
        if lower_row_is_ratio and not self.has_data:
            raise ValueError("Cannot do lower_row_is_ratio if no data")

        gc_stash = [] # to stop stuff being deleted

        # Setup canvas and pads
        canvas = ROOT.TCanvas("c_"+cu.get_unique_str(), "", 1600, 600)
        canvas.SetBottomMargin(0.0)
        canvas.SetTopMargin(0.0)
        canvas.SetLeftMargin(0.0)
        canvas.SetRightMargin(0.0)
        canvas.cd()

        # For each selection group, create 2 pads, one for gluon, one for quark
        upper_pads, lower_pads = [], []
        n_pads = len(selection_groups)

        def _is_multiline(text):
            return any(f in text for f in ['frac', 'splitline'])

        # gap between right end of plots and edge of canvas, used for legend
        final_reading_factor=1.15
        right_margin = 0.19*final_reading_factor
        # pad_left_titles_gap = 0.01 # gap between pad_left_titles and all plots
        pad_to_pad_gap = 0.0025  # gap between plot pad columns
        # how far in from the left the first plotting pad starts. used for y axis title
        left_margin = 0.055 if any(_is_multiline(t) for t in [upper_row_label, lower_row_label]) else 0.04
        # figure out width per pad - get total width available, then divide by number of pads
        pad_width = (1 - left_margin - right_margin - (pad_to_pad_gap*(n_pads - 1))) / n_pads

        pad_offset_bottom = 0.01 # spacing between bottom of RMS pad and canvas edge
        pad_offset_top = 0.08  # spacing between top of mean pad and canvas edge - for CMS and lumi text

        # per-pad margins: these determine where the hist axes lie,
        # and are fractions of the **pad** width/height, not the global canvas
        pad_right_margin = 0.04
        pad_left_margin = 0.2

        # bottom margin includes space for x axis labels & and text at the bottom
        # note that we apply it BOTH to the upper and lower pads,
        # even though the labels are on the lower pads, since we need the pads
        # to be exactly the same size, in order to get man things the same size,
        # e.g. ticks, labels, hatching, all of which depend on pad size,
        # and not histogram axes size!
        pad_bottom_margin = 0.49 if label_every_bin else 0.3

        # extra bit to add to the top margins of lower and upper pads
        # to ensure y axis numbers don't get cut off
        pad_top_margin = 0.012

        # pad height is constrained by the available space
        # (i.e. after pad_offset_top and pad_offset_bottom),
        # and the fact that we want the pads to overlap exactly by both the
        # top and bottom margins, to ensure that the x axes align vertically
        pad_height = (1 - pad_offset_top - pad_offset_bottom) / (2 - pad_bottom_margin - pad_top_margin)

        for isel in range(n_pads):
            canvas.cd()
            pad_start_x = left_margin + (isel*pad_to_pad_gap) + (isel*pad_width)
            pad_end_x = pad_start_x + pad_width
            # Create pad for upper hist - upper half of this column
            pad_upper = ROOT.TPad(cu.get_unique_str(), "", pad_start_x, 1-pad_offset_top-pad_height, pad_end_x, 1-pad_offset_top)
            ROOT.SetOwnership(pad_upper, False)
            pad_upper.SetFillColor(isel+2)
            pad_upper.SetFillStyle(3004)
            pad_upper.SetFillStyle(4000)
            pad_upper.SetTopMargin(pad_top_margin)
            pad_upper.SetBottomMargin(pad_bottom_margin)
            pad_upper.SetRightMargin(pad_right_margin)
            pad_upper.SetLeftMargin(pad_left_margin)
            pad_upper.SetTicks(1, 1)
            pad_upper.Draw()
            upper_pads.append(pad_upper)

            canvas.cd()
            # Create pad for lower hist - lower half of this column
            pad_lower = ROOT.TPad(cu.get_unique_str(), "", pad_start_x, pad_offset_bottom, pad_end_x, pad_offset_bottom+pad_height)
            ROOT.SetOwnership(pad_lower, False)
            pad_lower.SetFillColor(isel+2)
            pad_lower.SetFillStyle(3003)
            pad_lower.SetFillStyle(4000)
            pad_lower.SetTopMargin(pad_top_margin)
            pad_lower.SetBottomMargin(pad_bottom_margin)
            pad_lower.SetRightMargin(pad_right_margin)
            pad_lower.SetLeftMargin(pad_left_margin)
            pad_lower.SetTicks(1, 1)
            pad_lower.Draw()
            lower_pads.append(pad_lower)

        def _get_nbins(selection_group):
            return len(selection_group['selections'])

        if not label_every_bin:
            if len(set([_get_nbins(selection_groups[i]) for i in range(n_pads)])) > 1:
                raise RuntimeError("Cannot disable label_every_bin, since differing number of bins")

        key_bin_names = ['(%d)' % (i+1) for i in range(_get_nbins(selection_groups[0]))]

        # Construct lower row of hists if doing ratio
        # Note that we remove data error bars - those get added as a separate hashed region
        if lower_row_is_ratio:
            lower_row_hist_groups = []
            for hist_group in upper_row_hist_groups:
                # always assumes data is [0]
                data_hist = hist_group[0]
                data_no_errors = data_hist.Clone()
                cu.remove_th1_errors(data_no_errors)
                new_hist_group = self.construct_data_ratio_hists(data_no_errors, hist_group[1:])
                lower_row_hist_groups.append(new_hist_group)

        data_total_ratio = None  # for legend
        hists_extra=[]

        # Now draw the histograms
        for isel, (selection_group, upper_pad, upper_hist_group, lower_pad, lower_hist_group) \
            in enumerate(zip(selection_groups, upper_pads, upper_row_hist_groups, lower_pads, lower_row_hist_groups)):

            # DO UPPER ROW
            # ------------
            upper_pad.cd()

            # Draw hists - reverse order, so data on top
            for ind, hist in enumerate(upper_hist_group[::-1]):
                draw_opt = "E1"
                if ind != 0:
                    draw_opt += " SAME"
                hist.Draw(draw_opt)

            upper_draw_hist = upper_hist_group[-1]
            # remove x axis label
            xax = upper_draw_hist.GetXaxis()
            xax.SetLabelSize(0)
            xax.SetTitleSize(0)

            factor = 1.5
            xax.SetTickLength(xax.GetTickLength()*factor)

            yax = upper_draw_hist.GetYaxis()
            label_size_fudge = 1.6
            yax.SetLabelSize(yax.GetLabelSize()*factor*label_size_fudge)
            label_offset_fudge = 4
            yax.SetLabelOffset(yax.GetLabelOffset()*factor*label_offset_fudge)
            tick_fudge = 2
            yax.SetTickLength(yax.GetTickLength()*factor*tick_fudge)

            primary = 4 # 4 is the magic primary number, not 3
            secondary = 5
            # n_divisions = 510  # default
            n_divisions_opts = [primary, secondary, 0, True]
            yax.SetNdivisions(*n_divisions_opts)

            if ylims_upper is not None:
                y_down, y_up = ylims_upper[isel]
            else:
                if upper_lower_same_ylim:
                    # Use all hists
                    y_up, y_down = self.calc_auto_ylim(upper_hist_group + lower_hist_group)
                else:
                    # Set range using bin contents + error bars
                    y_up, y_down = self.calc_auto_ylim(upper_hist_group + lower_hist_group)

            upper_draw_hist.SetMinimum(y_down)
            upper_draw_hist.SetMaximum(y_up)


            # DO LOWER ROW
            # ------------
            lower_pad.cd()
            draw_opt = "E1"
            for ind, hist in enumerate(lower_hist_group[::-1][:1]):
                if hist.GetFillStyle()==3395:
                    hist_upper=hist.Clone(hist.GetTitle()+"upper")
                    hists_extra+=[hist_upper]
                    for b in range(hist.GetNbinsX()):
                      if hist.GetBinContent(b+1)>0:
                       hist_upper.SetBinContent(b+1,hist.GetBinContent(b+1)+hist.GetBinError(b+1)*(1 if hist.GetBinError(b+1)>1e-10 else 1e100))
                       hist_upper.SetBinError(b+1,0)
                    hist_upper.SetFillStyle(3395)
                    hist_upper.SetFillColor(ROOT.kRed-9)
                    hist_upper.SetLineColor(0)
                    hist_lower=hist.Clone(hist.GetTitle()+"lower")
                    hists_extra+=[hist_lower]
                    for b in range(hist.GetNbinsX()):
                      if hist.GetBinContent(b+1)>0:
                       hist_lower.SetBinContent(b+1,hist.GetBinContent(b+1)-hist.GetBinError(b+1)*(1 if hist.GetBinError(b+1)>1e-10 else 1e100))
                       hist_lower.SetBinError(b+1,0)
                    hist_lower.SetFillStyle(1001)
                    hist_lower.SetFillColor(10)
                    hist_lower.SetLineColor(0)
                    for b in range(hist.GetNbinsX()):
                       hist.SetBinError(b+1,hist.GetBinError(b+1)*(1e-100 if hist.GetBinError(b+1)>1e-10 else 1))
                    hist.Draw(draw_opt)
                    hist_upper.Draw("F SAME")
                    hist_lower.Draw("F SAME")
                    hist.Draw(draw_opt + " SAME")
                    hist.Draw("AXIS SAME")
                else:
                    hist.Draw(draw_opt)
            for ind, hist in enumerate(lower_hist_group[::-1][1:]):
                hist.Draw(draw_opt + " SAME")

            lower_draw_hist = lower_hist_group[-1]
            xax = lower_draw_hist.GetXaxis()

            if lower_row_is_ratio:
                # need to have drawn hists already. now draw other things,
                # then redraw hists on top

                # draw hashed area for data uncertainty
                # Easiest way to get errors right is to do data (with 0 errors)
                # and divide by data (with errors), as if you had MC = data with 0 error
                data_hist = upper_hist_group[0]
                data_no_errors = data_hist.Clone()
                cu.remove_th1_errors(data_no_errors)
                data_total_ratio = data_no_errors.Clone()
                data_total_ratio.Divide(data_hist)
                data_total_ratio.SetFillStyle(3754)
                data_total_ratio.SetFillColor(ROOT.kGray+1)
                data_total_ratio.SetLineWidth(0)
                data_total_ratio.SetMarkerSize(0)
                data_total_ratio.Draw("E2 SAME")
                gc_stash.append(data_total_ratio)

                # create line at 1
                ax_min, ax_max = xax.GetXmin(), xax.GetXmax()
                line = ROOT.TLine(ax_min, 1, ax_max, 1)
                line.SetLineWidth(2)
                line.SetLineStyle(2)
                line.SetLineColor(ROOT.kBlack)
                line.Draw()
                gc_stash.append(line)

                # now draw all the other hists
                for ind, hist in enumerate(lower_hist_group):
                    hist.Draw(draw_opt + "SAME")

            xax.CenterLabels()
            xax.LabelsOption("v")

            xax.SetTickLength(xax.GetTickLength()*factor)
            xax.SetLabelSize(xax.GetLabelSize()*factor*2)
            xax.SetLabelOffset(xax.GetLabelOffset()*factor)

            yax = lower_draw_hist.GetYaxis()
            yax.SetLabelSize(yax.GetLabelSize()*factor*label_size_fudge)
            yax.SetLabelOffset(yax.GetLabelOffset()*factor*label_offset_fudge)
            yax.SetTickLength(yax.GetTickLength()*factor*tick_fudge)  # extra bit of fudging, probably becasue pads are sligthly different sizes
            yax.SetNdivisions(*n_divisions_opts)

            # Instead of labelling every bin, replace with numbers,
            # and add key to plot
            if not label_every_bin:
                for i in range(1, lower_draw_hist.GetNbinsX()+1):
                    xax.SetBinLabel(i, key_bin_names[i-1])
                xax.LabelsOption("h")
                xax.SetLabelSize(xax.GetLabelSize()*1.5)
                xax.SetLabelOffset(xax.GetLabelOffset()*2)

            if ylims_lower is not None:
                y_down, y_up = ylims_lower[isel]
            else:
                if upper_lower_same_ylim:
                    if ylims_upper is not None:
                        y_down, y_up = ylims_upper[isel]
                else:
                    y_up, y_down = self.calc_auto_ylim(lower_hist_group)

            lower_draw_hist.SetMinimum(y_down)
            lower_draw_hist.SetMaximum(y_up)

            if lower_row_is_ratio:
                # avoid awkward bit where 1 is right at the top, or not even shown
                y_range = y_up - y_down
                min_ratio_y_max = 1. + y_range*0.1
                if lower_draw_hist.GetMaximum() < min_ratio_y_max:
                    lower_draw_hist.SetMaximum(min_ratio_y_max)

            # Draw variable name
            var_latex = ROOT.TLatex()
            var_latex.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignBottom)
            var_latex.SetTextSize(0.13)
            # var_latex.SetTextFont(42)
            # these are relative to the RMS pad! not the canvas
            # put in middle of the plot (if 0.5, woud look off-centre)
            var_x = 0.5*(1-pad_right_margin-pad_left_margin) + pad_left_margin
            var_y = 0.03  # do by eye
            if not label_every_bin:
                var_y = 0.05
            var_latex.DrawLatexNDC(var_x, var_y, selection_group['label'])
            gc_stash.append(var_latex)

        canvas.cd()

        # Add legend
        # Put legend_header + legend in own TPad
        leg_y_top = 1-pad_offset_top+0.03
        leg_left = upper_pads[-1].GetAbsXlowNDC() + upper_pads[-1].GetAbsWNDC()
        leg_right = 1-0.005
        leg_y_bottom = leg_y_top-(1.5*upper_pads[0].GetAbsHNDC())
        leg_pad = ROOT.TPad("leg_pad_"+cu.get_unique_str(), "", leg_left, leg_y_bottom, leg_right, leg_y_top)
        ROOT.SetOwnership(leg_pad, False)  # important! otherwise seg fault
        # leg_pad.SetFillColor(ROOT.kGreen)
        # leg_pad.SetFillStyle(3004)
        # leg_pad.SetFillStyle(4000)
        leg_pad.SetLeftMargin(0)
        leg_pad.SetRightMargin(0)
        leg_pad.SetTopMargin(0)
        leg_pad.SetBottomMargin(0)
        leg_pad.Draw()
        leg_pad.cd()
        gc_stash.append(leg_pad)

        n_leg_entries = len(self.other_samples) # other samples
        if not self.only_yoda_data:
            n_leg_entries += 2  # mc, alt mc
        if self.has_data:
            n_leg_entries += 1
        if lower_row_is_ratio:
            n_leg_entries += 1 # for data hashing

        # to figure out legend height, account for any #splitline or \n in labels
        multiline_extra = 2.5
        for label in [self.mc_label, self.alt_mc_label]:
            if '#splitline' in label or '\n' in label:
                n_leg_entries += multiline_extra
        for sample in self.other_samples:
            label = sample['style_dict'].get('label', sample['key'])
            if '#splitline' in label or '\n' in label:
                n_leg_entries += multiline_extra

        leg_entry_spacing = 0.07
        leg = ROOT.TLegend(0., max(0.5,0.98 - (leg_entry_spacing*n_leg_entries)), 1, 0.98) # relative to leg_pad

        # Replace legend markers with graph to get correct error bar endings
        # Yes this is ridiculous
        dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        dummy_data = dummy_gr.Clone()
        self._style_data_hist(dummy_data)
        dummy_mc = dummy_gr.Clone()
        self._style_mc_hist(dummy_mc)
        dummy_alt_mc = dummy_gr.Clone()
        self._style_alt_mc_hist(dummy_alt_mc)
        dummy_theory = dummy_gr.Clone()
        self._style_theory_hist(dummy_theory)

        def _add_entry(obj, label, option):
            parts = label.split("\n")
            if len(parts) == 1:
                leg.AddEntry(obj, label, option)
            elif len(parts) == 2:
                leg.AddEntry(obj, "#splitline{%s}{%s}" % (parts[0], parts[1]), option)
            else:
                raise RuntimeError("Cannot handle > 1 \n in legend label")

            # for line_ind, line in enumerate(label.split("\n")):
            #     if line_ind == 0:
            #         leg.AddEntry(obj, line, option)
            #     else:
            #         leg.AddEntry(0, line, "")

        # TODO: tie together labels and hists from upper_hists
        if self.has_data:
            _add_entry(dummy_data, "Data", "LEP")
            if lower_row_is_ratio:
                _add_entry(data_total_ratio, "Data uncertainty", "F")
        if not self.only_yoda_data:
            _add_entry(dummy_mc, self.mc_label, "LEP")
            _add_entry(dummy_alt_mc, self.alt_mc_label, "LEP")
            if add_theory:
              _add_entry(dummy_theory, "NLO+NLL'+NP", "FL")

        dummy_other = []  # keep reference
        for sample in self.other_samples:
            dummy_sample = dummy_gr.Clone()
            self._style_hist(dummy_sample, **sample['style_dict'])
            dummy_other.append(dummy_sample)
            sample_label = sample['style_dict'].get('label', sample['key'])
            _add_entry(dummy_sample, sample_label, "LEP")

        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        # leg.SetFillColor(ROOT.kBlue)
        # leg.SetFillStyle(3004)
        leg.SetTextSize(0.1)
        leg.SetEntrySeparation(leg.GetEntrySeparation()*2.5)
        if lower_row_is_ratio:
            leg.SetEntrySeparation(leg.GetEntrySeparation()*1.2)

        leg.SetMargin(0.12)
        leg.SetTextAlign(12)
        leg.Draw()

        canvas.cd()

        # Add key for bins below legend
        if not label_every_bin:
            # the bottom is the offset, plus a little bit
            # the top of the text is the margin (scaled according to pad height),
            # minus a bit for labels

            # this is the bottom of the legend in canvas co-ords,
            # note that leg.GetY1 is relative to leg_pad height, so needs scaling
            # to convert to canvas coords
            leg_bottom = leg_y_bottom + (0.48 * (leg_y_top-leg_y_bottom))
            leg_bottom -= 0.03  # offset
            # set bottom of key_text to align with bottom axis of plots,
            # plus some fudging since it doesn't quite align
            # make X2 off the canvas, since it tends to crop it in a bit
            key_text = ROOT.TPaveText(leg_left, pad_offset_bottom+(pad_bottom_margin*pad_height)-0.03-0.05, 1.005, leg_bottom, "NB NDC")
            key_text.SetFillStyle(4000)
            key_text.SetFillColor(ROOT.kWhite)
            # key_text.SetFillStyle(1001)
            # key_text.SetFillColor(ROOT.kYellow)
            key_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
            # if you don't set text size, it scales it automatically
            final_reading_factor=1.15
            key_text.SetTextSize(0.047*final_reading_factor)
            key_text.SetMargin(0.01)

            # generate all the key items
            text_items = []
            for key_bin, selection in zip(key_bin_names, selection_groups[0]['selections']):
                key = selection[1]
                key_str = "%s " % key_bin
                inset = (len(key_str)+2) * ' '
                # inset = ''
                if '\n' not in key and '#splitline' not in key:
                    text_items.append("%s%s" % (key_str, key))
                elif '\n' in key:
                    parts = key.split('\n')
                    text_items.append("%s%s" % (key_str, parts[0]))
                    text_items.extend(["%s%s" % (inset, p) for p in parts[1:]])

            for t in text_items:
                key_text.AddText(t)

            key_text.Draw()
            gc_stash.append(key_text)

        # Add upper/lower y labels
        text_width = left_margin

        def _calc_text_x(is_multiline=False):
            x = (0.55 * left_margin) - (0.5*text_width)
            if is_multiline:
                x += 0.005
            return x

        text_x = _calc_text_x(is_multiline=_is_multiline(upper_row_label))
        # let ROOT center it, just make the box the height of the axis
        text_y_end = 1 - pad_offset_top - (pad_top_margin*upper_pads[0].GetAbsHNDC())
        text_y = 1 - pad_offset_top - upper_pads[0].GetAbsHNDC() + (upper_pads[0].GetAbsHNDC()*pad_bottom_margin)
        upper_text = ROOT.TPaveText(text_x, text_y, text_x+text_width, text_y_end, "NDC NB")
        upper_text.SetFillStyle(0)
        upper_text.SetBorderSize(0)
        upper_text.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignCenter)
        text = upper_text.AddText(upper_row_label)
        text.SetTextAngle(90)
        text_size = 0.05
        text.SetTextSize(text_size)
        upper_text.Draw()

        text_x = _calc_text_x(is_multiline=_is_multiline(lower_row_label))
        text_y = (lower_pads[0].GetAbsHNDC() * pad_bottom_margin) + pad_offset_bottom
        text_y_end = (lower_pads[0].GetAbsHNDC() * (1-pad_top_margin)) + pad_offset_bottom
        lower_text = ROOT.TPaveText(text_x, text_y, text_x+text_width, text_y_end, "NDC NB")
        lower_text.SetFillStyle(0)
        lower_text.SetBorderSize(0)
        lower_text.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignCenter)
        l_text = lower_text.AddText(lower_row_label)
        l_text.SetTextAngle(90)
        l_text.SetTextSize(text_size)
        lower_text.Draw()

        # Add CMS text
        cms_latex = ROOT.TLatex()
        cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        cms_latex.SetTextFont(42)
        extraOverCmsTextSize  = 0.76
        cms_latex.SetTextSize(0.045/extraOverCmsTextSize)
        # Get the text sitting just above the axes of the mean plot
        # Axes end inside the mean pad at (1-top_margin), but this has
        # to be scaled to canvas NDC
        # Then add a little extra spacing ontop to separate text from axes line
        latex_height = 1 - pad_offset_top - (upper_pads[0].GetAbsHNDC() * upper_pads[0].GetTopMargin()) + 0.02
        # Want it to start at the left edge of the first plot
        start_x = left_margin + (pad_width*pad_left_margin)
        cms_latex.DrawLatexNDC(start_x, latex_height, "#font[61]{CMS}")

        cms_latex = ROOT.TLatex()
        cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        cms_latex.SetTextFont(42)
        cms_latex.SetTextSize(0.045)
        start_x = left_margin + (pad_width*pad_left_margin) + 0.1
        # # start_x = 100
        if self.is_preliminary:
            if self.has_data:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[52]{ Preliminary}")
            else:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[52]{ Preliminary Simulation}")
        else:
            if self.has_data:
                cms_latex.DrawLatexNDC(start_x, latex_height, "")
            else:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[52]{ Simulation}")
        cms_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
        # # Get the lumi text aligned to right edge of axes
        # # i.e. 1-pad_right_margin, but remember to scale by pad width
        end_x = 1 - right_margin - (upper_pads[0].GetAbsWNDC() * pad_right_margin)
        # end_x = 0.985  # to match legend
        # Figure out why physical region(s) are in plot, to get correct lumi
        do_zpj, do_dijet = False, False
        for sg in selection_groups:
            if any(["ZPlusJets" in x[0] for x in sg['selections']]):
                do_zpj = True
            if any(["Dijet" in x[0] for x in sg['selections']]):
                do_dijet = True
        lumi = cu.get_lumi_str(do_dijet=do_dijet, do_zpj=do_zpj)
        cms_latex.DrawLatexNDC(end_x, latex_height, " %s fb^{-1} (13 TeV)" % lumi)
        gc_stash.append(cms_latex)

        canvas.Update()
        cu.check_dir_exists_create(os.path.dirname(os.path.abspath(output_file)))
        canvas.SaveAs(output_file)
        print(output_file)


def unpack_slim_unfolding_root_file_uproot(input_tfile, region_name, angle_name, pt_bins):
    tdir = "%s/%s" % (region_name, angle_name)
    indices = range(len(pt_bins)-1 if isinstance(pt_bins[0], float) else len(pt_bins))  # adjust if pairs of bins, or individual bin edges

    unfolding_stat_err_hists = [
        input_tfile["%s/unfolded_stat_err_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]
    unfolding_total_err_hists = [
        input_tfile["%s/unfolded_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]
    unfolding_total_err_ematrices = [
        input_tfile["%s/unfolded_total_ematrix_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]
    truth_hists = [
        input_tfile["%s/hist_truth_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]
    alt_truth_hists = [
        input_tfile["%s/alt_hist_truth_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]

    return dict(
        unfolding_stat_err_hists=unfolding_stat_err_hists,
        unfolding_total_err_hists=unfolding_total_err_hists,
        unfolding_total_ematrices=unfolding_total_err_ematrices,
        truth_hists=truth_hists,
        alt_truth_hists=alt_truth_hists,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5input",
                        help="Read analysis data from H5 input file (from running extract_unfolding_summary_stats.py)")
    parser.add_argument("--h5inputRivet",
                        help="Read RIVET data from H5 input file (from running extract_rivet_summary_stats.py)")
    parser.add_argument("--outputDir",
                        default=None,
                        help='Output directory (default is the source dir)')
    parser.add_argument("--onlyYodaData",
                        action='store_true',
                        help='Only plot data & Yoda inputs (ignore MG+Py, H++)')
    parser.add_argument("--onlyDataOldMC",
                        action='store_true',
                        help='Only plot data & non-Yoda MCs (i.e. only MG+Py, H++)')
    parser.add_argument("--doMetricVsPt",
                        action='store_true',
                        help='Plot metric vs pT (e.g. mean vs pT)')
    parser.add_argument("--doSummaryBins",
                        action='store_true',
                        help='Plot summary bin plots')
    parser.add_argument("--paper",
                        action='store_true',
                        help='Don\'t add "Preliminary" to plots')
    parser.add_argument("--supplementary",
                        action='store_true',
                        help='Add "Supplementary" to plots')
    args = parser.parse_args()

    # Get input data
    if not any([args.h5input, args.h5inputRivet]):
        raise RuntimeError("Need one of --h5input, --h5inputRivet")

    if not any([args.doMetricVsPt, args.doSummaryBins]):
        raise RuntimeError("You should do at least one of --doMetricVsPt, --doSummaryBins")

    if args.onlyYodaData and not args.h5inputRivet:
        raise RuntimeError("--onlyYodaData requires --h5inputRivet")

    if args.onlyYodaData and args.onlyDataOldMC:
        raise RuntimeError("--onlyYodaData and --onlyDataOldMC are mutually exclusive, can only use one of them")

    if not args.outputDir:
        if args.h5input:
            args.outputDir = os.path.dirname(os.path.abspath(args.h5input))
        else:
            args.outputDir = os.getcwd()
        print("Setting output dir to", args.outputDir)

    # ----------------------------------------------------------------------
    # Read in data from h5 files
    # -----------------------------------------------------------------------
    print("Reading in unfolding data from existing HDF5 file...")
    if not os.path.isfile(args.h5input):
        raise IOError("Cannot find --h5input file")

    with pd.HDFStore(args.h5input) as store:
        df = store['df']
    print(df.head())
    print("# entries:", len(df.index))

    yoda_labels = []
    if args.h5inputRivet and not args.onlyDataOldMC:
        print("Reading in RIVET data from existing HDF5 file...")
        with pd.HDFStore(args.h5inputRivet) as store:
            df_rivet = store['df']
        print(df_rivet.head())
        print("# Rivet entries:", len(df_rivet.index))
        # Figure out YODA entries from column names
        mean_columns = [c.replace("mean_err_", '') for c in df_rivet.columns if 'mean_err_' in c]
        print(mean_columns)
        df = pd.merge(df, df_rivet, how='outer')
        # sort manually, but check the ones we want are actually in the dataframe
        only_these_yoda_labels = ['Pythia8_CP2', 'Pythia8_CP5', 'Herwig7_CH3', 'Sherpa_LO', 'Sherpa_LO+jet'][:]
        for mc in mean_columns:
            if only_these_yoda_labels:
                if mc in only_these_yoda_labels:
                    yoda_labels.append(mc)
                else:
                    print("Skipping column", mc)

        # yoda_labels_ideal = ['Pythia8_CP2', 'Pythia8_CP5', 'Herwig7_CH3', 'Sherpa']
        # for yl in yoda_labels_ideal:
        #     if yl not in mean_columns:
        #         warnings.warn("Missing yoda input %s from rivet file" % yl)
        #     else:
        #         yoda_labels.append(yl)
        print("Setting yoda_labels to", yoda_labels)

    convert_df_types(df)
    print(df.columns)
    print(df.head())
    print(df.dtypes)
    print(len(df.index), 'entries in dataframe')

    # Filter only regions/algos/angles in the dataframe, since it could have
    # been modified earlier
    all_jet_algos = [
        {'label': 'AK4', 'name': 'ak4puppi'},
        {'label': 'AK8', 'name': 'ak8puppi'}
    ]
    jet_algos = [j for j in all_jet_algos if j['name'] in df['jet_algo'].unique()]
    # print("Plotting jet_algos:", jet_algos)

    all_regions = [
        get_dijet_config('', central=True, groomed=False),
        get_dijet_config('', central=False, groomed=False),
        get_dijet_config('', central=True, groomed=True),
        get_dijet_config('', central=False, groomed=True),
        get_zpj_config('', groomed=False),
        get_zpj_config('', groomed=True),
    ]
    regions = [r for r in all_regions if r['name'] in df['region'].unique()]
    # print("Plotting regions:", regions)

    angles = [a for a in qgc.COMMON_VARS if a.var in df['angle'].unique()]
    # print("Plotting angles:", angles)
    charged_only_angles = [a for a in angles if "charged" in a.var]
    charged_and_neutral_angles = [a for a in angles if "charged" not in a.var]
    # charged_and_neutral_angles[0] = charged_and_neutral_angles[1]

    # --------------------------------------------------------------------------
    # Do all the plotting
    # --------------------------------------------------------------------------
    print("*"*120)
    print("* Plotting time!")
    print("*"*120)
    plotter = SummaryPlotter(jet_algos,
                             angles,
                             qgc.PT_UNFOLD_DICT['signal_gen'],
                             qgc.PT_UNFOLD_DICT['signal_zpj_gen'],
                             df,
                             args.outputDir,
                             has_data=True,
                             is_preliminary=not args.paper,
                             is_supplementary=args.supplementary,
                             only_yoda_data=args.onlyYodaData)

    has_dijet = any(["Dijet" in r['name'] for r in regions])
    has_zpj = any(["ZPlusJet" in r['name'] for r in regions])

    # Add extra samples we got from dataframe
    for ylabel in yoda_labels:
        ylabel = dataframe_yoda_key(ylabel)
        if ylabel not in SAMPLE_STYLE_DICTS:
            print("No entry found in SAMPLE_STYLE_DICTS for", ylabel, ", using defaults")
        plotter.add_sample(key=ylabel,
                           style_dict=SAMPLE_STYLE_DICTS.get(ylabel, dict()))

    # For summary plots
    low_pt = 120
    high_pt = 614
    high_pt = 800
    high_pt = 1000

    ak4_str = "AK4"
    ak8_str = "AK8"

    if (not args.h5inputRivet and not args.onlyYodaData) or args.onlyDataOldMC:
        plotter.filename_append = "_onlyDataNomMC%s" % plotter.filename_append

    filename_append = plotter.filename_append

    g_selections, q_selections = None, None

    # charged_only_template = "#splitline{{{jet_str}, {pt_str}}}{{      Charged-only}}"
    # groomed_template = "#splitline{{{jet_str}, {pt_str}}}{{         Groomed}}"

    normal_template = "{jet_str}, {pt_str}"
    charged_only_template = "{jet_str}, {pt_str},\ncharged"
    groomed_template = "{jet_str}, {pt_str},\ngroomed"

    # gev_template = "p_{{T}}#in  [{:g}, {:g}] GeV"
    # tev_template = "p_{{T}}#in  [{:g}, {:g}] TeV"

    gev_template = "[{:g}, {:g}] GeV"
    tev_template = "[{:g}, {:g}] TeV"

    g_legend_header = "Gluon-enriched jets:\n" + qgc.Dijet_CEN_LABEL

    if has_dijet:
        if args.doMetricVsPt:
            plotter.plot_dijet_means_vs_pt_all()
            plotter.plot_dijet_rms_vs_pt_all()
            plotter.plot_dijet_delta_vs_pt_all()

        pt_bins = qgc.PT_UNFOLD_DICT['signal_gen']
        low_pt_bin = np.where(pt_bins == low_pt)[0][0]
        low_pt_str = gev_template.format(low_pt, pt_bins[low_pt_bin+1])

        high_pt_bin = np.where(pt_bins == high_pt)[0][0]
        high_pt_str = tev_template.format(high_pt/1000, pt_bins[high_pt_bin+1]/1000)

        if args.doSummaryBins:
            # DIJET CENTRAL/GLUON
            # ---------------------------------------------
            # Create selection queries for mean/RMS/delta summary plots
            g_selections = []
            for angle in charged_and_neutral_angles:
                this_angle_str = "%s (%s)" % (angle.name, angle.lambda_str)
                # this_angle_str = "%s" % (angle.lambda_str)
                this_selection = [
                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s"' % (low_pt_bin, angle.var),
                        normal_template.format(jet_str=ak4_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s"' % (high_pt_bin, angle.var),
                        normal_template.format(jet_str=ak4_str, pt_str=high_pt_str)),

                    ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s"' % (low_pt_bin, angle.var),
                        normal_template.format(jet_str=ak8_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                        charged_only_template.format(jet_str=ak4_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="Dijet_central_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                        groomed_template.format(jet_str=ak4_str, pt_str=low_pt_str)),
                ]
                g_selections.append({'label': this_angle_str, 'selections': this_selection})

            dijet_central_legend_header = qgc.Dijet_CEN_LABEL
            plotter.plot_mean_rms_bins_summary(
                selections=g_selections,
                legend_header=dijet_central_legend_header,
                output_file=os.path.join(args.outputDir, "dijet_central_mean_rms_summary%s.pdf" % (filename_append))
            )

            plotter.plot_delta_bins_summary(
                selections=g_selections,
                legend_header=dijet_central_legend_header,
                output_file=os.path.join(args.outputDir, "dijet_central_delta_summary%s.pdf" % (filename_append))
            )

            plotter.plot_mean_rms_bins_summary(
                selections=g_selections,
                legend_header=g_legend_header,
                output_file=os.path.join(args.outputDir, "gluon_mean_rms_summary%s.pdf" % (filename_append))
            )

            plotter.plot_delta_bins_summary(
                selections=g_selections,
                legend_header=g_legend_header,
                output_file=os.path.join(args.outputDir, "gluon_delta_summary%s.pdf" % (filename_append))
            )

            # DIJET FORWARD
            # ---------------------------------------------

            # Create selection queries for mean/RMS/delta summary plots
            selections = []
            for angle in charged_and_neutral_angles:
                this_angle_str = "%s (%s)" % (angle.name, angle.lambda_str)
                # this_angle_str = "%s" % (angle.lambda_str)
                this_selection = [
                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s"' % (low_pt_bin, angle.var),
                        normal_template.format(jet_str=ak4_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s"' % (high_pt_bin, angle.var),
                        normal_template.format(jet_str=ak4_str, pt_str=high_pt_str)),

                    ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s"' % (low_pt_bin, angle.var),
                        normal_template.format(jet_str=ak8_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                        charged_only_template.format(jet_str=ak4_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="Dijet_forward_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                        groomed_template.format(jet_str=ak4_str, pt_str=low_pt_str)),
                ]
                selections.append({'label': this_angle_str, 'selections': this_selection})

            dijet_fwd_legend_header = qgc.Dijet_FWD_LABEL
            plotter.plot_mean_rms_bins_summary(
                selections=selections,
                legend_header=dijet_fwd_legend_header,
                output_file=os.path.join(args.outputDir, "dijet_forward_mean_rms_summary%s.pdf" % (filename_append))
            )

            plotter.plot_delta_bins_summary(
                selections=selections,
                legend_header=dijet_fwd_legend_header,
                output_file=os.path.join(args.outputDir, "dijet_forward_delta_summary%s.pdf" % (filename_append))
            )

    if has_zpj:
        if args.doMetricVsPt:
            plotter.plot_zpj_means_vs_pt_all()
            plotter.plot_zpj_rms_vs_pt_all()
            plotter.plot_zpj_delta_vs_pt_all()

        pt_bins = qgc.PT_UNFOLD_DICT['signal_zpj_gen']
        low_pt_bin = np.where(pt_bins == low_pt)[0][0]
        low_pt_str = gev_template.format(low_pt, pt_bins[low_pt_bin+1])

        high_pt = 326
        high_pt_bin = np.where(pt_bins == high_pt)[0][0]
        high_pt_str = gev_template.format(high_pt, pt_bins[high_pt_bin+1])

        if args.doSummaryBins:
            selections = []
            for angle in charged_and_neutral_angles:
                this_angle_str = "%s (%s)" % (angle.name, angle.lambda_str)
                # this_angle_str = "%s" % (angle.lambda_str)
                this_selection = [
                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                        normal_template.format(jet_str=ak4_str, pt_str=low_pt_str)),

                    # ignore high pt bin as not useful - same composition as dijet but fewer stats
                    # ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (high_pt_bin, angle.var),
                    #     normal_template.format(jet_str=ak4_str, pt_str=high_pt_str)),
                    #     # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=high_pt_str)),

                    ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                        normal_template.format(jet_str=ak8_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                        charged_only_template.format(jet_str=ak4_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="ZPlusJets_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                        groomed_template.format(jet_str=ak4_str, pt_str=low_pt_str)),
                ]
                selections.append({'label': this_angle_str, 'selections': this_selection})

            zpj_legend_header = qgc.ZpJ_LABEL
            plotter.plot_mean_rms_bins_summary(
                selections=selections,
                legend_header=zpj_legend_header,
                output_file=os.path.join(args.outputDir, "zpj_mean_rms_summary%s.pdf" % (filename_append))
            )

            plotter.plot_delta_bins_summary(
                selections=selections,
                legend_header=zpj_legend_header,
                output_file=os.path.join(args.outputDir, "zpj_delta_summary%s.pdf" % (filename_append))
            )

    if has_dijet and has_zpj and args.doSummaryBins:
        low_pt = 120
        pt_bins = qgc.PT_UNFOLD_DICT['signal_gen']
        low_pt_bin = np.where(pt_bins == low_pt)[0][0]
        low_pt_str = gev_template.format(low_pt, pt_bins[low_pt_bin+1])

        high_pt = 1000
        high_pt_bin = np.where(pt_bins == high_pt)[0][0]
        high_pt_str = tev_template.format(high_pt/1000, pt_bins[high_pt_bin+1]/1000)

        # QUARK
        # ---------------------------------------------
        # Create selection queries for mean/RMS/delta summary plots
        q_selections = []
        for angle in charged_and_neutral_angles:
            this_angle_str = "%s (%s)" % (angle.name, angle.lambda_str)
            # this_angle_str = "%s" % (angle.lambda_str)
            this_selection = [
                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                    normal_template.format(jet_str=ak4_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s"' % (high_pt_bin, angle.var),
                    normal_template.format(jet_str=ak4_str, pt_str=high_pt_str)),

                ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                    normal_template.format(jet_str=ak8_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                    charged_only_template.format(jet_str=ak4_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="ZPlusJets_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                    groomed_template.format(jet_str=ak4_str, pt_str=low_pt_str)),
            ]
            q_selections.append({'label': this_angle_str, 'selections': this_selection})

        q_legend_header = "Quark-enriched jets:\n%s:\n%s\n%s:\n%s" % (low_pt_str, qgc.ZpJ_LABEL, high_pt_str, qgc.Dijet_FWD_LABEL)
        plotter.plot_mean_rms_bins_summary(
            selections=q_selections,
            legend_header=q_legend_header,
            output_file=os.path.join(args.outputDir, "quark_mean_rms_summary%s.pdf" % (filename_append))
        )

        plotter.plot_delta_bins_summary(
            selections=q_selections,
            legend_header=q_legend_header,
            output_file=os.path.join(args.outputDir, "quark_delta_summary%s.pdf" % (filename_append))
        )

        qg_legend_header = g_legend_header + "\n" + q_legend_header

        # hard-coded y axis limits to ensure e.g. data w/ mg_py, h++ is same as data + other MC
        # I CBA to try and loop over all, figure it out, etc
        # ...re-evaluate as necessary
        ylims_g_q_enriched_mean = [
            (0.14, 0.42),
            (0.07, 0.24),
            (0.02, 0.16),
            (0, 55),
            (0.06, 0.38)
        ]

        plotter.plot_q_g_mean_bins_summary(
            gluon_selections=g_selections,
            quark_selections=q_selections,
            legend_header=qg_legend_header,
            output_file=os.path.join(args.outputDir, "quark_gluon_mean_summary%s.pdf" % (filename_append)),
            ylims=ylims_g_q_enriched_mean,
        )

        ylims_g_q_ratio = [
            (0.98, 1.39),
            (0.95, 1.57),
            (0.98, 1.65),
            (0.96, 1.45),
            (0.55, 1.05)
        ]
        ylims_g_q_data_mc_ratio = [
            (0.98, 1.14),
            (0.96, 1.19),
            (0.94, 1.35),
            (0.95, 1.18),
            (0.77, 1.06)
        ]
        plotter.plot_q_vs_g_bins_summary(
            gluon_selections=g_selections,
            quark_selections=q_selections,
            legend_header=qg_legend_header,
            output_file=os.path.join(args.outputDir, "quark_gluon_data_mc_mean_summary%s.pdf" % (filename_append)),
            ylims_upper=ylims_g_q_ratio,
            ylims_lower=ylims_g_q_data_mc_ratio
        )

