#!/usr/bin/env python

"""produce plots comparing different kinematics between different flavours"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product
import numpy as np
np.seterr(all='raise')
from array import array
from uuid import uuid1

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgg
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

# Control output format
OUTPUT_FMT = "pdf"


def get_list_of_obj(directory):
    key_list = directory.GetListOfKeys()
    return [x.GetName() for x in key_list]


def do_1D_plot(hists, output_filename, components_styles_dicts=None, 
               draw_opts="NOSTACK HISTE", do_ratio=True, logy=False,
               normalise_hists=True, title=""):
    hists = [h.Clone(h.GetName() + str(uuid1())) for h in hists]
    contributions = [Contribution(hist, normalise_hist=normalise_hists, **csd)
                     for hist, csd in zip(hists, components_styles_dicts)]

    if len(contributions) == 0:
        return
    
    # Ignore if all empty objs
    total_entries = sum(c.obj.GetEntries() for c in contributions)
    if total_entries == 0:
        return
    
    min_val = min([h.GetMinimum(0) for h in hists])
    max_val = max([h.GetMaximum() for h in hists])
    if logy:
        ylim = [0.5*min_val, 10*max_val]
    else:
        # ylim = [0.5*min_val, 1.5*max_val]
        ylim = [0, 1.5*max_val]
    p = Plot(contributions, what='hist', 
             ytitle="p.d.f." if normalise_hists else "N", 
             title=title,
             subplot_type="ratio" if do_ratio else None, 
             subplot_title="#splitline{Ratio wrt}{%s}" % contributions[0].label,
             subplot=contributions[0], ylim=ylim)
    p.legend.SetX1(0.55)
    p.legend.SetX2(0.95)
    p.legend.SetY1(0.7)
    p.legend.SetY2(0.85)
    p.plot(draw_opts)
    p.main_pad.cd()

    if logy:
        p.main_pad.SetLogy(1)
                
    # p.save(os.path.join(output_dir, obj_name+".%s" % (OUTPUT_FMT)))
    p.save(output_filename)


def do_all_1D_projection_plots_in_dir(directories, output_dir, components_styles_dicts=None, 
                                      draw_opts="NOSTACK HISTE", do_ratio=True, 
                                      normalise_hists=True,
                                      jet_config_str="",
                                      signal_mask=None):
    """
    Given a set of TDirs, loop over all 2D hists, do projection hists for chosen bins, and plot all TDir contributions on a canvas for comparison.

    components_styles_dicts should be a list of style dicts, one for each directory/contribution to a plot. 
    This should include label for this component.

    """
    if len(directories) != len(components_styles_dicts):
        raise IndexError("Need same number of style dicts and directories")

    list_of_obj = [get_list_of_obj(d) for d in directories]
    print(list_of_obj[0])

    # check all have same list of plots
    if not all(x == list_of_obj[0] for x in list_of_obj):
        print("Different number of object in the TDirectorys")

    common_list_of_obj = set(list_of_obj[0])
    for l in list_of_obj[1:]:
        common_list_of_obj = common_list_of_obj & set(l)


    # pt_bins = [(20, 40), (40, 60), (60, 80), (100, 120), (160, 200), (260, 300), (500, 600), (1000, 2000)]
    pt_bins = [(20, 40), (40, 60), (90, 100), (100, 120), (160, 200), (280, 300), (300, 350), (400, 500), (500, 600), (1000, 2000)]

    for obj_name in common_list_of_obj:
        if "flav" in obj_name:
            continue
        if obj_name in [
                        # 'eta_jet1_vs_eta_jet2', 
                        'phi_jet1_vs_pt_jet1', 'phi_jet2_vs_pt_jet1', 
                        # 'reliso_mu1_vs_pt_jet1', 'reliso_mu2_vs_pt_jet1', 
                        # 'dphi_mumu_jet1_vs_pt_jet1', 
                        # 'dphi_mumu_vs_pt_jet1', 
                        'pt_jet1_z_pt_jet2_z_ratio']:
            continue
        
        objs = [d.Get(obj_name).Clone(obj_name + str(uuid1())) for d in directories]
        
        # Special case for Dijets summing JetHT and ZeroBias
        if len(objs) == 3:
            objs[0].Scale(35860)  # QCD scaling
            objs[2].Scale(1235009.27580634)  # ZB scaling
            objs[1].Add(objs[2])
            objs = objs[:2]
            components_styles_dicts = components_styles_dicts[:2]

        # Ignore TH1s
        if not isinstance(objs[0], (ROOT.TH2F, ROOT.TH2D, ROOT.TH2I)):
            do_1D_plot(objs, components_styles_dicts=components_styles_dicts, 
                       draw_opts=draw_opts, do_ratio=do_ratio, normalise_hists=normalise_hists, logy=True,
                       title=jet_config_str,
                       output_filename=os.path.join(output_dir, obj_name+".%s" % (OUTPUT_FMT)))
        else:

            for pt_min, pt_max in pt_bins:
                # print(pt_min, pt_max)
                rebin = 1
                # exceptions...why didn't I pick the same number of bins...
                do_not_rebin = any([
                    "n_jets" in obj_name, 
                    "n_mu" in obj_name, 
                    "met_sig" in obj_name, 
                    obj_name.startswith('dphi_mumu'), 
                    obj_name.startswith('pt_jet3_frac'), 
                    obj_name.startswith('pt_jet1_jet2_ratio'),
                    obj_name.startswith('pt_jet1_z_ratio'),
                    obj_name.startswith('pt_jet2_z_ratio'),
                    obj_name.startswith('dphi_jet1_z_vs_pt_jet1'),
                    # obj_name.startswith('m_jj'),
                ])
                if not do_not_rebin:
                    if objs[0].GetNbinsX() % 5 == 0 and objs[0].GetNbinsX() >= 100:
                        rebin = 5
                    elif objs[0].GetNbinsX() % 4 == 0 and objs[0].GetNbinsX() >= 80:
                        rebin = 2
                    elif objs[0].GetNbinsX() % 3 == 0 and objs[0].GetNbinsX() >= 60:
                        rebin = 3
                if obj_name.startswith("m_jj"):
                    rebin = 2
                if "reliso" in obj_name:
                    rebin = 2
                hists = [qgg.get_projection_plot(ob, pt_min, pt_max) for ob in objs]

                do_1D_plot(hists, components_styles_dicts=components_styles_dicts, 
                           draw_opts=draw_opts, do_ratio=do_ratio, 
                           normalise_hists=normalise_hists, logy=False,
                           title="#splitline{%s}{%d < p_{T}^{jet 1} < %d GeV}" % (jet_config_str, pt_min, pt_max),
                           output_filename=os.path.join(output_dir, obj_name+"_pt%dto%d.%s" % (pt_min, pt_max, OUTPUT_FMT)))


def do_dijet_distributions(root_dir):
    """Do plots comparing different different inputs in dijet region"""
    root_files = [qgc.QCD_FILENAME, qgc.JETHT_FILENAME, qgc.ZB_FILENAME][:]
    root_files = [qgc.QCD_PYTHIA_ONLY_FILENAME, qgc.JETHT_FILENAME, qgc.ZB_FILENAME][:]
    root_files = [cu.open_root_file(os.path.join(root_dir, r)) for r in root_files]
    
    directories = [cu.get_from_file(rf, "Dijet") for rf in root_files]
    mc_col = ROOT.kAzure
    data_col = ROOT.kRed
    zb_col = ROOT.kGreen+2
    csd = [
        {"label": "QCD PY8 MC", "line_color": mc_col, "fill_color": mc_col, "marker_color": mc_col, "marker_style": 20, "fill_style": 0},
        {"label": "JetHT + ZeroBias Data", "line_color": data_col, "fill_color": data_col, "marker_color": data_col, "marker_style": 22, "fill_style": 0},
        {"label": "ZeroBias Data", "line_color": zb_col, "fill_color": zb_col, "marker_color": zb_col, "marker_style": 23, "fill_style": 0},
    ]
    jet_config_str = qgc.extract_jet_config(os.path.basename(os.path.abspath(root_dir)))

    # Compare shapes
    do_all_1D_projection_plots_in_dir(directories=directories, 
                                      output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_normalised"),
                                      components_styles_dicts=csd,
                                      jet_config_str=jet_config_str)
                                      
    # Compare yields
    do_all_1D_projection_plots_in_dir(directories=directories, 
                                      output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_absolute"),
                                      components_styles_dicts=csd,
                                      jet_config_str=jet_config_str,
                                      normalise_hists=False)


def do_zpj_distributions(root_dir):
    """Do plots comparing different different inputs in dijet region"""
    root_files = [qgc.DY_FILENAME, qgc.SINGLE_MU_FILENAME]
    root_files = [cu.open_root_file(os.path.join(root_dir, r)) for r in root_files]
    
    directories = [cu.get_from_file(rf, "ZPlusJets") for rf in root_files]
    mc_col = ROOT.kAzure
    data_col = ROOT.kRed
    csd = [
        {"label": "DYJetsToLL MC", "line_color": mc_col, "fill_color": mc_col, "marker_color": mc_col, "marker_style": 20, "fill_style": 0},
        {"label": "SingleMu Data", "line_color": data_col, "fill_color": data_col, "marker_color": data_col, "marker_style": 21, "fill_style": 0},
    ]

    # Compare shapes
    do_all_1D_projection_plots_in_dir(directories=directories, 
                                      output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_normalised"),
                                      components_styles_dicts=csd,
                                      filter_noisy=False)
                                      
    # Compare yields
    do_all_1D_projection_plots_in_dir(directories=directories, 
                                      output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_absolute"),
                                      components_styles_dicts=csd,
                                      normalise_hists=False,
                                      filter_noisy=False)


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
        # if (os.path.isfile(os.path.join(workdir, qgc.QCD_FILENAME)) and 
        #     os.path.isfile(os.path.join(workdir, qgc.JETHT_FILENAME)) and
        #     os.path.isfile(os.path.join(workdir, qgc.ZB_FILENAME))):
            do_dijet_distributions(workdir)

        # if (os.path.isfile(os.path.join(workdir, qgc.DY_FILENAME)) and 
        #     os.path.isfile(os.path.join(workdir, qgc.SINGLE_MU_FILENAME))):
        #     do_zpj_distributions(workdir)

    sys.exit(0)
