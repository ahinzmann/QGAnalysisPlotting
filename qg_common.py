"""Common vars for doing QG plots"""

from __future__ import print_function, division

import numpy as np
import argparse
from collections import namedtuple

# Use ints for ROOT colours as importing ROOT seems OTT

# Purple & blue look ~ same to red-green colourblind people
# (V.Dark blue & light purple do look a bit diff)
# Red and blue does look different
# DY_COLOUR =  880 # kViolet
DY_COLOUR =  628
QCD_COLOUR = 867
QCD_CEN_COLOUR = 867
QCD_FWD_COLOUR = 868

# kBlue = 600
MGPY_QCD_COLOUR = 601
# kViolet = 880
HERWIGPP_QCD_COLOUR = 800+3 #880-2

# SINGLE_MU_COLOUR =  884 # kViolet+3
SINGLE_MU_COLOUR =  635
JETHT_COLOUR = 600 # kBlue
JETHT_CEN_COLOUR = 601 # kBlue
JETHT_FWD_COLOUR = 600 # kBlue
ZB_COLOUR = 416+1 # kGreen+1

HERWIGPP_DY_COLOUR = 797

# You should use markers as well for colourblindness
DY_MARKER = 20
DY_MARKER = 'circle'

QCD_MARKER = 22
QCD_MARKER = 'triangleUp'

QCD_CEN_MARKER = 22
QCD_CEN_MARKER = 'triangleUp'

QCD_FWD_MARKER = 23
QCD_FWD_MARKER = 'triangleDown'

# When comparing e.g. PU bins
# FIXME: update for colourblindness
# DY_COLOURS = [880, 881, 884, 873]
DY_COLOURS = [880, 881, 797, 902]
DY_COLOURS *= 3
QCD_COLOURS = [867, 600, 853, 882, 835, 592, 883]
QCD_COLOURS *= 3  # for lots

PT_BINS = [(80, 100), (100, 200), (400, 500), (1000, 2000), (80, 2000)]
PT_BINS = [(80, 100), (100, 120), (200, 250), (400, 500), (1000, 2000), (80, 2000)]
THEORY_PT_BINS = [(100, 200), (400, 500), (1000, 2000), (80, 2000)]
THEORY_PT_BINS = [(100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 800), (800, 1000), (1000, 2000), (80, 2000)]
# THEORY_PT_BINS = [(20, 40), (40, 60), (60, 80), (80, 100), (100, 120),
THEORY_PT_BINS = [(20, 25), (40, 50), (50, 55), (60, 70), (70, 80), (80, 90), (90, 100), (100, 120),
                  (120, 160), (160, 200), (200, 220), (280, 300), (300, 350), (400, 410),
                  (400, 500), (500, 600), (500, 550), (600, 800), (800, 1000),
                  (1000, 1400), (1400, 2000)]
# PT_BINS = [
# (0, 29), (29,38), (38,50), (50,65), (65,88), (88,120), (120,150), (150,186), (186,254),
# (254,326), (326,408), (408,481), (481,614), (614,800), (800,1000), (1000,1300),
# (1300,1700), (1700,2200), (2200,3000), (3000,4000), (4000,5000), (5000,10000),
# ]

PT_BINS = [
(15, 23), (23, 30), (30, 38), (38, 50), (50, 65), (65, 88), (88, 120), (120, 150),
(150, 186), (186, 254), (254, 326), (326, 408), (408, 481), (481, 614), (614, 800),
(800, 1000), (1000, 2000), (2000, 6500)
]

PT_BINS_ZPJ = [
(15, 23), (23, 30), (30, 38), (38, 50), (50, 65), (65, 88), (88, 120), (120, 150),
(150, 186), (186, 254), (254, 326), (326, 408), (408, 481), (481, 614), (614, 800),
(800, 6500)
]

# PT_BINS = [
# (29,38), (38,50), (50,65), (65,88), (88,114), (114,145), (145,180), (180,254),
# (254,318), (318,400), (400,500), (500,625), (625,800), (800,1000), (1000,1300),
# (1300,1700), (1700,2200), (2200,3000), (3000,4000), (4000,5000), (5000,10000),
# ]
#

PT_BINS = [(30,50),
(50, 65), (65, 88), (88, 120), (120, 150), (150, 186),
(186, 254), (254, 326), (326, 408), (408, 481), (481, 614),
(614, 800), (800, 1000), (1000, 4000)
][1:]

PT_BINS_ZPJ = [
(50, 65), (65, 88), (88, 120), (120, 150), (150, 186),
(186, 254), (254, 326), (326, 408), (408, 481), (481, 614),
(614, 800), (800, 2000)
]

PT_BINS_INC_UFLOW = [(30, 38), (38, 50)]
PT_BINS_INC_UFLOW += PT_BINS

PT_BINS_ZPJ_INC_UFLOW = [(30, 38), (38, 50)]
PT_BINS_ZPJ_INC_UFLOW += PT_BINS_ZPJ


def construct_fine_binning(coarse_bin_edges):
    fine_bin_edges = []
    for x, y in zip(coarse_bin_edges[:-1], coarse_bin_edges[1:]):
        fine_bin_edges.append(x)
        fine_bin_edges.append(0.5*(x+y))
    fine_bin_edges.append(coarse_bin_edges[-1])
    return np.array(fine_bin_edges, dtype='d')


def construct_all_fine_binnings(var_dict):
    # Construct fine binning from splitting coarser bins
    # Coarser bins dervied from determine_lambda_binning.py
    for group_name, group_dict in var_dict.items():
        for angle_name, angle_dict in group_dict.items():
            if angle_dict['reco'] is None:
                angle_dict['reco'] = construct_fine_binning(angle_dict['gen'])
            var_dict[group_name][angle_name] = angle_dict


Angle = namedtuple("Angle", ['var', 'kappa', 'beta', 'name', "lambda_str", "colour", "mathmode"])
COMMON_VARS_WITH_FLAV = [
    # The order here is important: determines that in summary plots
    # Angle("jet_multiplicity", 0, 0, "Multiplicity", "#lambda_{0}^{0}", 2),
    Angle('jet_LHA', 1, 0.5, "LHA", "#lambda_{0.5}^{1}", 600, r"\text{LHA}"),
    Angle('jet_width', 1, 1, "Width", "#lambda_{1}^{1}", 861, r"\text{Width}"),
    Angle('jet_thrust', 1, 2, "Thrust", "#lambda_{2}^{1}", 617, r"\text{Thrust}"),
    Angle("jet_puppiMultiplicity", 0, 0, "Multiplicity", "#lambda_{0}^{0}", 2, r"\text{Multiplicity}"),
    Angle('jet_pTD', 2, 0, "(p_{T}^{D})^{2}", "#lambda_{0}^{2}", 418, r"(p_{T}^{D})^{2}"),
    # charged-only constit
    # Angle("jet_multiplicity_charged", 0, 0, "Multiplicity (charged)", "#lambda_{0}^{0}", 2),
    Angle('jet_LHA_charged', 1, 0.5, "LHA (charged-only)", "#lambda_{0.5}^{1}", 600, r"\text{LHA (charged-only)}"),
    Angle('jet_width_charged', 1, 1, "Width (charged-only)", "#lambda_{1}^{1}", 861, r"\text{Width (charged-only)}"),
    Angle('jet_thrust_charged', 1, 2, "Thrust (charged-only)", "#lambda_{2}^{1}", 617, r"\text{Thrust (charged-only)}"),
    Angle("jet_puppiMultiplicity_charged", 0, 0, "Multiplicity (charged-only)", "#lambda_{0}^{0}", 2, r"\text{Multiplicity (charged-only)}"),
    Angle('jet_pTD_charged', 2, 0, "(p_{T}^{D})^{2} (charged-only)", "#lambda_{0}^{2}", 418, r"(p_{T}^{D})^{2}~\text{(charged-only)}"),
    # Angle('jet_multiplicity', 0, 0, "Multiplicity", "#lambda_{0}^{0}", 2),  # leave last as we don't want it in our RIVET numbering
    Angle('jet_flavour', 0, 0, "Flavour", "PDGID", 7, r"\text{Flavour}"),
    Angle("jet_genParton_flavour", 0, 0, "Flavour", "PDGID", 8, r"\text{Flavour}")
]

COMMON_VARS = COMMON_VARS_WITH_FLAV[:-2]


def lower_angle_name(angle):
    lower_angle_name = angle.name
    if ('LHA' not in lower_angle_name
        and "_{T}" not in lower_angle_name
        and "PUPPI" not in lower_angle_name):
        lower_angle_name = lower_angle_name[0].lower() + lower_angle_name[1:]
    return lower_angle_name


ANGLE_REBIN_DICT = {
    "jet_puppiMultiplicity": [0.0, 6.0, 9.0, 12.0, 18.0, 150.0],
    'jet_multiplicity': [0.0, 6.0, 9.0, 12.0, 18.0, 150.0],
    'jet_pTD': [0.0, 0.1, 0.13, 0.17, 0.23, 0.33, 0.51, 0.85, 1.0],
    'jet_LHA': [0.0, 0.29, 0.37, 0.44, 0.5, 0.56, 0.62, 0.68, 0.75, 1.0],
    'jet_width': [0.0, 0.12, 0.18, 0.24, 0.3, 0.36, 0.43, 0.51, 1.0],
    'jet_thrust': [0.0, 0.04, 0.08, 0.12, 0.17, 0.24, 0.33, 1.0],
}

PT_UNFOLD_DICT = {
    # "signal_gen": np.array([50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 6500], dtype='d'),

    "underflow_gen": np.array([15, 30, 38, 50], dtype='d'),
    "signal_gen": np.array([50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 4000], dtype='d'),

    # "signal_gen": np.array([65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 4000], dtype='d'),
    # "underflow_gen": np.array([30, 38, 50, 65], dtype='d'),

    # "signal_zpj_gen": np.array([50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 2000], dtype='d'),
    "signal_zpj_gen": np.array([50, 65, 88, 120, 150, 186, 254, 326, 408, 1500], dtype='d'),
    # "signal_zpj_gen": np.array([50, 65, 88, 120, 150, 186, 254, 408, 614, 2000], dtype='d'),
    "underflow_zpj_gen": np.array([15, 30, 38, 50], dtype='d'),
}


for pt_name in list(PT_UNFOLD_DICT.keys()):
    new_name = pt_name.replace("_gen", "_reco")
    PT_UNFOLD_DICT[new_name] = construct_fine_binning(PT_UNFOLD_DICT[pt_name])

# remove bins below reco pt cut
reco_jet_pt = 30
PT_UNFOLD_DICT['underflow_reco'] = PT_UNFOLD_DICT['underflow_reco'][PT_UNFOLD_DICT['underflow_reco'] >= reco_jet_pt]
PT_UNFOLD_DICT['underflow_zpj_reco'] = PT_UNFOLD_DICT['underflow_zpj_reco'][PT_UNFOLD_DICT['underflow_zpj_reco'] >= reco_jet_pt]


# new version of WTA, cen+fwd dijet, target 0.5
# from ./workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts/determine_lambda_binning_target0p5/binning_uhh2.AnalysisModuleRunner.MC.MC_QCD.txt
VAR_UNFOLD_DICT_TARGET0p5 = {
    'jet_puppiMultiplicity': {
        # 'gen': np.array([0.0, 10, 15, 20, 27, 50, 75, 100, 150], dtype='d'),
        'gen': np.array([0.0, 10, 15, 20, 30, 50, 75, 100, 150], dtype='d'),
        # 'gen': np.array([0.0, 10, 26, 50, 76, 150 ], dtype='d'),
        # 'gen': np.array([0,  4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 50,  60,  70,  80,  90, 100, 110, 120, 130, 140, 150], dtype='d'),  # super fine
        'reco': None
    },
    'jet_pTD': {
        # 'gen': np.array([0.0, 0.07, 0.1, 0.15, 0.24, 0.45, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.07, 0.1, 0.15, 0.24, 1.0], dtype='d'),
        # 'gen': np.array([0.0, 0.07, 0.1, 0.15, 0.22, 1.0], dtype='d'), # widen last bin
        # 'gen': np.array([0.0, 0.04, 0.08, 0.12, 0.14, 0.18, 0.22, 0.26, 0.32, 1.0], dtype='d'),  # super fine
        'reco': None
    },
    'jet_LHA': {
        # 'gen': np.array([0.0, 0.14, 0.22, 0.29, 0.35, 0.42, 0.49, 0.56, 0.64, 0.75, 1.0], dtype='d'),
        # 'gen': np.array([0.0, 0.14, 0.22, 0.29, 0.35, 0.42, 0.49, 0.56, 0.64, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.14, 0.21, 0.28, 0.34, 0.4, 0.47, 0.54, 0.61, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width': {
        # 'gen': np.array([0.0, 0.09, 0.145, 0.205, 0.28, 0.36, 0.445, 0.545, 1.0 ], dtype='d'),
        'gen': np.array([ 0.0, 0.09, 0.145, 0.205, 0.28, 0.355, 0.43, 0.515, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust': {
        # 'gen': np.array([0.0, 0.04, 0.08, 0.145, 0.225, 0.32, 0.445, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.04, 0.08, 0.14, 0.195, 0.25, 1.0], dtype='d'),
        'reco': None
    },

    'jet_puppiMultiplicity_charged': {
        # 'gen': np.array([0.0, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 18.0, 21.0, 25.0, 32.0, 91, 150.0], dtype='d'),
        'gen': np.array([0.0, 3.0, 5.0, 7.0, 10.0, 13.0, 16.0, 20.0, 30, 50, 75, 100, 150.0], dtype='d'),
        # 'gen': np.array([0, 4, 6, 8, 10, 14, 18, 26, 50, 76, 150], dtype='d'), # wide
        # 'gen': np.array([0,  4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 50,  60,  70,  80,  90, 100, 110, 120, 130, 140, 150 ], dtype='d'), # super fine binning
        'reco': None
    },
    'jet_pTD_charged': {
        # 'gen': np.array([0.0, 0.08, 0.1, 0.12, 0.14, 0.17, 0.2, 0.24, 0.29, 0.34, 0.4, 0.47, 0.55, 0.65, 0.76, 0.89, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.07, 0.09, 0.11, 0.14, 0.18, 0.23, 0.3, 0.39, 0.51, 0.64, 1.0], dtype='d'),
        # 'gen': np.array([0.0, 0.09, 0.11, 0.14, 0.18, 0.23, 0.3, 0.39, 0.51, 0.64, 1.0], dtype='d'),
        # 'gen': np.array([    0.0, 0.04, 0.08, 0.12, 0.14, 0.18, 0.22, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.6, 1.0], dtype='d'), # super fine
        'reco': None
    },
    'jet_LHA_charged': {
        # 'gen': np.array([0.0, 0.05, 0.09, 0.12, 0.15, 0.18, 0.22, 0.26, 0.31, 0.36, 0.42, 0.49, 0.56, 0.64, 0.75, 1.0], dtype='d'),
        # 'gen': np.array([0.0, 0.05, 0.09, 0.12, 0.15, 0.18, 0.22, 0.26, 0.31, 0.36, 0.42, 0.49, 0.56, 0.64, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.05, 0.09, 0.14, 0.19, 0.25, 0.32, 0.39, 0.47, 0.55, 0.63, 0.73, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width_charged': {
        # 'gen': np.array([0.0, 0.015, 0.025, 0.035, 0.05, 0.065, 0.085, 0.105, 0.13, 0.16, 0.195, 0.235, 0.28, 0.335, 0.4, 0.48, 0.585, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.01, 0.025, 0.045, 0.075, 0.11, 0.155, 0.21, 0.275, 0.34, 0.415, 0.495, 0.595, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust_charged': {
        # 'gen': np.array([0.0, 0.005, 0.01, 0.015, 0.025, 0.035, 0.05, 0.065, 0.085, 0.115, 0.15, 0.2, 0.265, 0.355, 0.48, 0.725, 1.], dtype='d'),
        'gen': np.array([0.0, 0.005, 0.015, 0.03, 0.06, 0.105, 0.16, 0.23, 0.315, 0.425, 1.0], dtype='d'),
        'reco': None
    },
}

# new version of WTA, cen+fwd dijet, target 0.5, pt>0 (except multiplicity), softdrop, fix lambda, charged clustering
# target 0.65 + spike smoothing for groomed charged
# from workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged/
VAR_UNFOLD_DICT_TARGET0p5 = {
    "ungroomed": {
        'jet_puppiMultiplicity': {
            'gen': np.array([-0.5, 9.5, 15.5, 21.5, 29.5, 39.5, 59.5, 99.5, 149.5], dtype='d'),
            'reco': None
        },
        'jet_pTD': {
            'gen': np.array([0.0, 0.06, 0.09, 0.13, 0.19, 0.3, 1.0], dtype='d'),
            'reco': None
        },
        'jet_LHA': {
            'gen': np.array([0.0, 0.17, 0.25, 0.32, 0.38, 0.45, 0.52, 0.59, 1.0], dtype='d'),
            'reco': None
        },
        'jet_width': {
            'gen': np.array([ 0.0, 0.105, 0.165, 0.23, 0.305, 0.38, 0.46, 0.55, 1.0], dtype='d'),
            'reco': None
        },
        'jet_thrust': {
            'gen': np.array([0.0, 0.05, 0.09, 0.15, 0.205, 1], dtype='d'),
            'reco': None
        },

        'jet_puppiMultiplicity_charged': {
            'gen': np.array([-0.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 17.5, 21.5, 25.5, 29.5, 35.5, 41.5, 59.5, 99.5], dtype='d'),
            'reco': None
        },
        'jet_pTD_charged': {
            'gen': np.array([0, 0.07, 0.09, 0.12, 0.15, 0.19, 0.24, 0.31, 0.4, 0.52, 0.69, 1], dtype='d'),
            'reco': None
        },
        'jet_LHA_charged': {
            'gen': np.array([0, 0.06, 0.11, 0.15, 0.19, 0.23, 0.27, 0.31, 0.35, 0.39, 0.44, 0.49, 0.54, 0.6, 1], dtype='d'),
            'reco': None
        },
        'jet_width_charged': {
            'gen': np.array([0, 0.0225, 0.04, 0.0575, 0.0775, 0.1, 0.125, 0.152, 0.185, 0.22, 0.26, 0.307, 0.362, 0.425, 0.497, 1], dtype='d'),
            'reco': None
        },
        'jet_thrust_charged': {
            'gen': np.array([0, 0.005, 0.0125, 0.0225, 0.035, 0.05, 0.07, 0.0925, 0.12, 0.152, 0.188, 0.228, 1], dtype='d'),
            'reco': None
        },
    },

    "groomed": {
        'jet_puppiMultiplicity': {
            'gen': np.array([-0.5, 7.5, 13.5, 19.5, 29.5, 39.5, 49.5, 75.5, 99.5, 149.5], dtype='d'),
            'reco': None
        },
        'jet_pTD': {
            'gen': np.array([0, 0.06, 0.09, 0.14, 0.22, 1], dtype='d'),
            'reco': None
        },
        'jet_LHA': {
            'gen': np.array([0, 0.1, 0.18, 0.26, 0.34, 0.42, 0.5, 0.57, 0.64, 1.], dtype='d'),
            'reco': None
        },
        'jet_width': {
            'gen': np.array([0, 0.02, 0.05, 0.095, 0.147, 0.225, 0.307, 0.388, 0.468, 1], dtype='d'),
            'reco': None
        },
        'jet_thrust': {
            'gen': np.array([0, 0.0025, 0.01, 0.025, 0.06, 0.12, 0.177, 0.23, 1], dtype='d'),
            'reco': None
        },

        'jet_puppiMultiplicity_charged': {
            'gen': np.array([-0.5, 3.5, 5.5, 9.5, 13.5, 17.5, 21.5, 27.5, 35.5, 49.5, 99.5], dtype='d'),
            'reco': None
        },
        'jet_pTD_charged': {
            'gen': np.array([0, 0.07, 0.09, 0.12, 0.16, 0.21, 0.27, 0.35, 0.44, 0.54, 0.66, 1], dtype='d'),
            'reco': None
        },
        'jet_LHA_charged': {
            'gen': np.array([0, 0.06, 0.09, 0.12, 0.15, 0.19, 0.23, 0.27, 0.32, 0.37, 0.42, 0.48, 0.54, 0.6, 1], dtype='d'),
            'reco': None
        },
        'jet_width_charged': {
            'gen': np.array([0, 0.0125, 0.0225, 0.035, 0.05, 0.07, 0.095, 0.128, 0.17, 0.225, 0.29, 0.365, 0.45, 1], dtype='d'),
            'reco': None
        },
        'jet_thrust_charged': {
            'gen': np.array([0, 0.0025, 0.005, 0.0075, 0.0125, 0.02, 0.0325, 0.05, 0.0775, 0.115, 0.16, 0.21, 1], dtype='d'),
            'reco': None
        },
    }
}


# construct_all_fine_binnings(VAR_UNFOLD_DICT_TARGET0p6)
# VAR_UNFOLD_DICT = VAR_UNFOLD_DICT_TARGET0p6

construct_all_fine_binnings(VAR_UNFOLD_DICT_TARGET0p5)
VAR_UNFOLD_DICT = VAR_UNFOLD_DICT_TARGET0p5

# Common labels for legends etc
ZpJ_LABEL = "Z+jet region"
DY_ZpJ_LABEL = "DY+jets MC, Z+jet region"
DY_ZpJ_GROOMED_LABEL = "DY+jets MC, Z+jet region\n(groomed)"
DY_ZpJ_GEN_LABEL = "DY+jets MC, Z+jet region\n(GenJets)"
DY_ZpJ_QFLAV_LABEL = "DY+jets MC, Z+jet region\n(uds-matched)"
DY_ZpJ_GFLAV_LABEL = "DY+jets MC, Z+jet region\n(g-matched)"

DY_Dijet_LABEL = "DY+jets MC, Dijet region"
DY_Dijet_GEN_LABEL = "DY+jets MC, Dijet region\n(GenJets)"
DY_Dijet_QFLAV_LABEL = "DY+jets MC, Dijet region\n(uds-matched)"
DY_Dijet_GFLAV_LABEL = "DY+jets MC, Dijet region\n(g-matched)"

QCD_ZpJ_LABEL = "QCD MC, Z+jet region"
QCD_ZpJ_GEN_LABEL = "QCD MC, Z+jet region\n(GenJets)"
QCD_ZpJ_QFLAV_LABEL = "QCD MC, Z+jet region\n(uds-matched)"
QCD_ZpJ_GFLAV_LABEL = "QCD MC, Z+jet region\n(g-matched)"

Dijet_LABEL = "Dijet region"
Dijet_GROOMED_LABEL = "Dijet region (groomed)"
QCD_Dijet_LABEL = "QCD MC, Dijet region"
QCD_Dijet_GEN_LABEL = "QCD MC, Dijet region\n(GenJets)"
QCD_Dijet_QFLAV_LABEL = "QCD MC, Dijet region\n(uds-matched)"
QCD_Dijet_GFLAV_LABEL = "QCD MC, Dijet region\n(g-matched)"

Dijet_CEN_LABEL = "Central dijet region"
QCD_Dijet_CEN_LABEL = "QCD MC, Central dijet region"
QCD_Dijet_CEN_GROOMED_LABEL = "QCD MC, Central dijet region (groomed)"
QCD_Dijet_CEN_GEN_LABEL = "QCD MC, Central dijet region\n(GenJets)"
QCD_Dijet_CEN_QFLAV_LABEL = "QCD MC, Central dijet region\n(uds-matched)"
QCD_Dijet_CEN_GFLAV_LABEL = "QCD MC, Central dijet region\n(g-matched)"

Dijet_FWD_LABEL = "Forward dijet region"
QCD_Dijet_FWD_LABEL = "QCD MC, Forward dijet region"
QCD_Dijet_FWD_GROOMED_LABEL = "QCD MC, Forward dijet region (groomed)"
QCD_Dijet_FWD_GEN_LABEL = "QCD MC, Forward dijet region\n(GenJets)"
QCD_Dijet_FWD_QFLAV_LABEL = "QCD MC, Forward dijet region\n(uds-matched)"
QCD_Dijet_FWD_GFLAV_LABEL = "QCD MC, Forward dijet region\n(g-matched)"

SINGLE_MU_LABEL = "SingleMu, Z+jet region"
JETHT_LABEL = "JetHT, Dijet region"
ZB_LABEL = "ZeroBias, Dijet region"
JETHT_ZB_LABEL = "JetHT+ZeroBias, Dijet region"

JETHT_CEN_LABEL = "JetHT, Central dijet region"
ZB_CEN_LABEL = "ZeroBias, Central dijet region"
JETHT_ZB_CEN_LABEL = "JetHT+ZeroBias, Central dijet region"

JETHT_FWD_LABEL = "JetHT, Forward dijet region"
ZB_FWD_LABEL = "ZeroBias, Forward dijet region"
JETHT_ZB_FWD_LABEL = "JetHT+ZeroBias, Forward dijet region"

# Dirs in ROOT files
ZPJ_RECOJET_RDIR = "ZPlusJets_QG"
DJ_RECOJET_RDIR = "Dijet_QG_tighter"
ZPJ_GENJET_RDIR = "ZPlusJets_genjet"
DJ_GENJET_RDIR = "Dijet_genjet"

# Common filenames
DY_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root"
DY_INCL_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_HT-INCL.root"
QCD_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_QCD.root"
QCD_PYTHIA_ONLY_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root"
QCD_FLAT_PYTHIA_ONLY_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_QCD_PYTHIA_FLAT.root"
ZPJ_ALL_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-All.root"

QCD_HERWIG_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root"
QCD_HERWIG_PTREWEIGHT_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD_PtReweight.root"

DY_HERWIG_INCL_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_Incl.root"
DY_HERWIG_LOW_PT_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_Incl_PartonKtMax300.root"
DY_HERWIG_HIGH_PT_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_JetKtMin170_PartonKtMin300.root"
DY_HERWIG_LOW_HIGH_PT_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_merged_PartonKtMin300.root"

DY_HERWIG_FILENAME = DY_HERWIG_LOW_HIGH_PT_FILENAME

DY_MG_HERWIG_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_MG_HERWIG_DYJetsToLL.root"
DY_AMCATNLO_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_AMCATNLO_DYJetsToLL.root"

SINGLE_MU_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root"
JETHT_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT.root"
ZB_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias.root"
JETHT_ZB_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias.root"

TTBAR_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_TTBAR.root"
WW_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_WW.root"
WZ_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_WZ.root"
ZZ_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_ZZ.root"

ZEROBIAS_RUNB_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias_RunB.root"
ZEROBIAS_RUNC_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias_RunC.root"
ZEROBIAS_RUND_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias_RunD.root"
ZEROBIAS_RUNE_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias_RunE.root"
ZEROBIAS_RUNF_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias_RunF.root"
ZEROBIAS_RUNG_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias_RunG.root"
ZEROBIAS_RUNH_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias_RunH.root"

JETHT_RUNB_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunB.root"
JETHT_RUNC_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunC.root"
JETHT_RUND_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunD.root"
JETHT_RUNE_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunE.root"
JETHT_RUNF_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunF.root"
JETHT_RUNG_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunG.root"
JETHT_RUNH_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunH.root"



DATA_TOTAL_UNC_STR = "Data total unc."
DATA_STAT_UNC_STR = "Data stat. unc."

SIM_DATA_STR = "Sim. / data"


def get_parser():
    """Return a parser to loop over several input dirs"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdirs',
                        nargs='+',
                        help='Workdir(s) with ROOT files to process. '
                        'Several dirs can be specified here, separated by a space.')
    parser.add_argument("-o", "--output",
                        help="Directory to put output plot dirs into",
                        default=".")
    parser.add_argument("--title",
                        help="Optional label to put on plot")
    return parser


def extract_jet_config(workdir):
    """Return jet algo/PUS from dir name"""
    import re
    res = re.search(r"(ak|ca)([0-9]+)(chs|puppi)?", workdir, flags=re.IGNORECASE)
    if res:
        algo = res.groups()[0]
        r = res.groups()[1]
        pus = ""
        if len(res.groups()) == 3:
            pus = res.groups()[2]
        return "%s%s %s" % (algo.upper(), r, pus.upper())
    else:
        print("Couldn't decypher jet algo")
        return ""
