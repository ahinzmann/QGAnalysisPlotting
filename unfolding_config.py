"""Common region dict generating function for unfolding"""


import os

import ROOT
from MyStyle import My_Style
My_Style.cd()

import qg_common as qgc

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


# FOR DIJET:
def get_dijet_config(source_dir, central=True, groomed=False):
    source_dir_systs = os.path.join(source_dir, "systematics_files")

    input_mc_qcd_mgpythia_tfile = os.path.join(source_dir, qgc.QCD_FILENAME)
    # input_mc_qcd_mgpythia_tfile = os.path.join(source_dir, "uhh2.AnalysisModuleRunner.MC.MC_QCD_jackknife10.root")
    input_mc_qcd_pythia_tfile = os.path.join(source_dir, qgc.QCD_PYTHIA_ONLY_FILENAME)
    input_mc_qcd_flat_pythia_tfile = os.path.join(source_dir, qgc.QCD_FLAT_PYTHIA_ONLY_FILENAME)
    input_mc_qcd_herwig_tfile = os.path.join(source_dir, qgc.QCD_HERWIG_FILENAME)
    input_mc_qcd_herwig_tfile_reweight = os.path.join(source_dir, qgc.QCD_HERWIG_PTREWEIGHT_FILENAME)
    
    input_mc_dy_mgpythia_tfile = os.path.join(source_dir, qgc.DY_FILENAME)
    input_mc_dy_herwig_tfile = os.path.join(source_dir, qgc.DY_HERWIG_FILENAME)

    input_jetht_tfile = os.path.join(source_dir, qgc.JETHT_ZB_FILENAME)

    dijet_region_dict_template = {
        "name": "Dijet",
        "dirname": "Dijet_QG_Unfold_central_tighter",
        "label": "Dijet",
        "data_tfile": input_jetht_tfile,
        "mc_tfile": input_mc_qcd_mgpythia_tfile,
        "mc_label": "MG5+Pythia8",
        
        "dirname_otherProc": "ZPlusJets_QG_Unfold",
        "mc_otherProc_tfile": input_mc_dy_mgpythia_tfile,
        "mc_otherProc_label": "MG5+Pythia8 (DY)",

        # "mc_tfile": input_mc_qcd_herwig_tfile,
        # "mc_label": "Herwig++",
        # "mc_tfile": input_mc_qcd_pythia_tfile,
        # "mc_label": "Pythia8",
        # "mc_tfile": input_mc_qcd_flat_pythia_tfile,
        # "mc_label": "Pythia8 (Flat)",
        "alt_mc_tfile": input_mc_qcd_herwig_tfile,
        "alt_mc_label": "Herwig++",
        # "alt_mc_tfile": input_mc_qcd_pythia_tfile,
        # "alt_mc_label": "Pythia8",
        # "alt_mc_tfile": input_mc_qcd_herwig_tfile_reweight,
        # "alt_mc_label": "Herwig++ (reweighted)",

        "alt_mc_otherProc_tfile": input_mc_dy_herwig_tfile,
        "alt_mc_otherProc_label": "Herwig++ (DY)",

        "tau_limits": None,  # user should set this
        "unreg_unfolder": None,  # set later if regularisation used
        "experimental_systematics": [
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
            {
                "label": "Shower & hadr.",
                "tfile": input_mc_qcd_herwig_tfile,
                # "colour": ROOT.kGreen-3,
                "colour": qgc.HERWIGPP_QCD_COLOUR,
            },
            # {
            #     "label": "Herwig++ (reweighted)",
            #     "tfile": input_mc_qcd_herwig_tfile_reweight,
            #     "colour": ROOT.kOrange+4,
            # },
            # {
            #     "label": "NNPDF30_lo_as_0118",
            #     "tfile": os.path.join(source_dir, "uhh2.AnalysisModuleRunner.MC.MC_QCD_NNPDF30_lo_as_0118.root"),
            #     "colour": 835,
            # },
            # {
            #     "label": "HERAPDF20_LO_EIG",
            #     "tfile": os.path.join(source_dir, "uhh2.AnalysisModuleRunner.MC.MC_QCD_HERAPDF20_LO_EIG.root"),
            #     "colour": 592,
            # },
            # {
            #     "label": "CT14lo",
            #     "tfile": os.path.join(source_dir, "uhh2.AnalysisModuleRunner.MC.MC_QCD_CT14lo.root"),
            #     "colour": 883,
            # },
        ],

        "scale_systematics": [
            {
                "label": "muR up, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure,
                "unfolder": None,
            },
            {
                "label": "muR down, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+1,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+2,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+3,
                "unfolder": None,
            },
            {
                "label": "muR down, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+4,
                "unfolder": None,
            },
            {
                "label": "muR up, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+5,
                "unfolder": None,
            },
        ],

        "model_systematics": [
            {
                "label": "muR up, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure,
                "unfolder": None,
            },
            {
                "label": "muR down, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+1,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+2,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+3,
                "unfolder": None,
            },
            {
                "label": "muR down, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+4,
                "unfolder": None,
            },
            {
                "label": "muR up, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.QCD_FILENAME),
                "colour": ROOT.kAzure+5,
                "unfolder": None,
            },
            {
                "label": "Herwig++",
                "tfile": input_mc_qcd_herwig_tfile,
                "colour": ROOT.kGreen-3,
                "unfolder": None,
            },
            # {
            #     "label": "Herwig++ (reweighted)",  # don't put in p_{T} as braces cause havoc
            #     "tfile": input_mc_qcd_herwig_tfile_reweight,
            #     "colour": ROOT.kOrange+4,
            #     "unfolder": None,
            # },
            {
                "label": "Pythia8",
                "tfile": input_mc_qcd_pythia_tfile,
                "colour": ROOT.kMagenta+3,
                "unfolder": None,
            },
        ],

        "pdf_systematics": [
            {
                "label": "PDF",  # this is a template entry, used for future
                "tfile": os.path.join(source_dir_systs, 'PDFvariationsTrue', qgc.QCD_FILENAME),
                "colour": ROOT.kCyan+2,
                "unfolder": None,
                "variations": range(100),  # list of all the variation #s to be used
            },
        ],

        "jackknife_input_variations": [
            {
                "label": "jackknife",  # this is a template entry, used for future
                "tfile": input_mc_qcd_mgpythia_tfile,
                "colour": ROOT.kCyan+2,
                "unfolder": None,
                "variations": range(25),  # list of all the variation #s to be used
            },
        ],

        "jackknife_response_variations": [
            {
                "label": "jackknife",  # this is a template entry, used for future
                "tfile": input_mc_qcd_mgpythia_tfile,
                "colour": ROOT.kCyan+2,
                "unfolder": None,
                "variations": range(25),  # list of all the variation #s to be used
            },
        ]

    }

    if central and not groomed:
        this_dict = dijet_region_dict_template.copy()
        this_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter'
        this_dict['label'] = qgc.Dijet_CEN_LABEL.replace(" region", "")
        this_dict['name'] = 'Dijet_central'
        return this_dict

    elif not central and not groomed:
        this_dict = dijet_region_dict_template.copy()
        this_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter'
        this_dict['label'] = qgc.Dijet_FWD_LABEL.replace(" region", "")
        this_dict['name'] = 'Dijet_forward'
        return this_dict

    elif central and groomed:
        this_dict = dijet_region_dict_template.copy()
        this_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter_groomed'
        this_dict['label'] = qgc.Dijet_CEN_LABEL.replace(" region", "")
        this_dict['name'] = 'Dijet_central_groomed'
        this_dict["dirname_otherProc"] = "ZPlusJets_QG_Unfold_groomed"
        return this_dict

    elif not central and groomed:
        this_dict = dijet_region_dict_template.copy()
        this_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter_groomed'
        this_dict['label'] = qgc.Dijet_FWD_LABEL.replace(" region", "")
        this_dict['name'] = 'Dijet_forward_groomed'
        this_dict["dirname_otherProc"] = "ZPlusJets_QG_Unfold_groomed"
        return this_dict


def get_zpj_config(source_dir, groomed=False):
    source_dir_systs = os.path.join(source_dir, "systematics_files")

    input_mc_dy_mgpythia_tfile = os.path.join(source_dir, qgc.DY_FILENAME)  # incl + HT binned
    input_mc_dy_mgpythia_incl_tfile = os.path.join(source_dir, qgc.DY_INCL_FILENAME)  # only inclusive
    input_mc_dy_mgherwig_tfile = os.path.join(source_dir, qgc.DY_MG_HERWIG_FILENAME)
    input_mc_dy_herwig_tfile = os.path.join(source_dir, qgc.DY_HERWIG_FILENAME)  # inclusive + my high pt sample
    
    input_mc_qcd_mgpythia_tfile = os.path.join(source_dir, qgc.QCD_FILENAME)  # HT binned
    input_mc_qcd_herwig_tfile = os.path.join(source_dir, qgc.QCD_HERWIG_FILENAME)

    input_singlemu_tfile = os.path.join(source_dir, qgc.SINGLE_MU_FILENAME)

    zpj_region_dict = {
        "name": "ZPlusJets",
        "dirname": "ZPlusJets_QG_Unfold",
        "label": qgc.ZpJ_LABEL.replace(" region", ""),
        "data_tfile": input_singlemu_tfile,

        "mc_tfile": input_mc_dy_mgpythia_tfile,
        "mc_label": "MG5+Pythia8",

        # "mc_tfile": input_mc_dy_mgpythia_incl_tfile,
        # "mc_label": "MG5+Pythia8 (inclusive)",

        "dirname_otherProc": "Dijet_QG_Unfold_central_tighter",
        "mc_otherProc_tfile": input_mc_qcd_mgpythia_tfile,
        "mc_otherProc_label": "MG5+Pythia8 (QCD)",

        "alt_mc_tfile": input_mc_dy_herwig_tfile,
        "alt_mc_label": "Herwig++",

        "alt_mc_otherProc_tfile": input_mc_qcd_herwig_tfile,
        "alt_mc_otherProc_label": "Herwig++ (QCD)",

        # "alt_mc_tfile": input_mc_dy_mgherwig_tfile,
        # "alt_mc_label": "MG+Herwig++",

        "tau_limits": None,
        "backgrounds": [
            {
                "name": "t#bar{t}",
                "tfile": os.path.join(source_dir, qgc.TTBAR_FILENAME),
                "rate_unc": 1.,
            },
            {
                "name": "WW",
                "tfile": os.path.join(source_dir, qgc.WW_FILENAME),
                "rate_unc": 1.
            },
            {
                "name": "WZ",
                "tfile": os.path.join(source_dir, qgc.WZ_FILENAME),
                "rate_unc": 1.
            },
            {
                "name": "ZZ",
                "tfile": os.path.join(source_dir, qgc.ZZ_FILENAME),
                "rate_unc": 1.
            },
        ],

        "experimental_systematics": [
            # {
            #     "label": "Luminosity up",
            #     "tfile": None,
            #     "factor": 0.025,
            #     "colour": ROOT.kCyan,
            # },
            # {
            #     "label": "Luminosity down",
            #     "tfile": None,
            #     "factor": -0.025,
            #     "colour": ROOT.kCyan,
            #     'linestyle': 2,
            # },
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
                "colour": ROOT.kOrange-3,
            },
            {
                "label": "Neutral hadron down",
                "tfile": os.path.join(source_dir_systs, 'neutralHadronShiftDown', qgc.DY_FILENAME),
                "colour": ROOT.kOrange-3,
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
            # {
            #     "label": "MG+Herwig++",
            #     "tfile": input_mc_dy_mgherwig_tfile,
            #     "colour": ROOT.kOrange-3,
            # },
            {
                "label": "Shower & hadr.",
                "tfile": input_mc_dy_herwig_tfile,
                # "colour": ROOT.kGreen-3,
                "colour": qgc.HERWIGPP_QCD_COLOUR,
                # "colour": 797,
            },
        ],

        "scale_systematics": [
            {
                "label": "muR up, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.DY_FILENAME),
                "colour": ROOT.kAzure,
                "unfolder": None,
            },
            {
                "label": "muR down, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+1,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+2,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+3,
                "unfolder": None,
            },
            {
                "label": "muR down, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+4,
                "unfolder": None,
            },
            {
                "label": "muR up, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+5,
                "unfolder": None,
            },
        ],

        "model_systematics": [
            {
                "label": "muR up, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.DY_FILENAME),
                "colour": ROOT.kAzure,
                "unfolder": None,
            },
            {
                "label": "muR down, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+1,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+2,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+3,
                "unfolder": None,
            },
            {
                "label": "muR down, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+4,
                "unfolder": None,
            },
            {
                "label": "muR up, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+5,
                "unfolder": None,
            },
            # {
            #     "label": "MG+Herwig++",
            #     "tfile": input_mc_dy_mgherwig_tfile,
            #     "colour": ROOT.kOrange-3,
            #     "unfolder": None,
            # },
            {
                "label": "Herwig++",
                "tfile": input_mc_dy_herwig_tfile,
                "colour": ROOT.kGreen-3,
                "unfolder": None,
            },
        ],

        "pdf_systematics": [
            {
                "label": "PDF",
                "tfile": os.path.join(source_dir_systs, 'PDFvariationsTrue', qgc.DY_FILENAME),
                "colour": ROOT.kCyan+2,
                "unfolder": None,
                "variations": range(100),
            },
        ],

        "jackknife_input_variations": [
            {
                "label": "jackknife",  # this is a template entry, used for future
                "tfile": input_mc_dy_mgpythia_tfile,
                "colour": ROOT.kCyan+2,
                "unfolder": None,
                "variations": range(25),  # list of all the variation #s to be used
            },
        ],

        "jackknife_response_variations": [
            {
                "label": "jackknife",  # this is a template entry, used for future
                "tfile": input_mc_dy_mgpythia_tfile,
                "colour": ROOT.kCyan+2,
                "unfolder": None,
                "variations": range(25),  # list of all the variation #s to be used
            },
        ]
    }

    if not groomed:
        return zpj_region_dict.copy()
    else:
        this_dict = zpj_region_dict.copy()
        this_dict['dirname'] = 'ZPlusJets_QG_Unfold_groomed'
        this_dict['name'] = 'ZPlusJets_groomed'
        this_dict['label'] = qgc.ZpJ_LABEL.replace(" region", "")
        this_dict["dirname_otherProc"] = "Dijet_QG_Unfold_central_tighter_groomed"
        return this_dict



def setup_regions_from_argparse(args):
    """Utility method for argparser, to setup list of region(s) based on args

    Assumes args object is setup like:
        parser.add_argument("source",
                            help="Source directory with ROOT files")
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
    """
    return setup_regions(do_dijet_central=args.doDijetCentral,
                         do_dijet_forward=args.doDijetForward,
                         do_dijet_central_groomed=args.doDijetCentralGroomed,
                         do_dijet_forward_groomed=args.doDijetForwardGroomed,
                         do_zpj=args.doZPJ,
                         do_zpj_groomed=args.doZPJGroomed,
                         source=args.source)


def setup_regions(do_dijet_central,
                  do_dijet_forward,
                  do_dijet_central_groomed,
                  do_dijet_forward_groomed,
                  do_zpj,
                  do_zpj_groomed,
                  source):
    """Setup list of region(s)"""
    regions = []
    if do_dijet_central:
        regions.append(get_dijet_config(source, central=True, groomed=False))

    if do_dijet_forward:
        regions.append(get_dijet_config(source, central=False, groomed=False))

    if do_dijet_central_groomed:
        regions.append(get_dijet_config(source, central=True, groomed=True))

    if do_dijet_forward_groomed:
        regions.append(get_dijet_config(source, central=False, groomed=True))

    if do_zpj:
        regions.append(get_zpj_config(source, groomed=False))

    if do_zpj_groomed:
        regions.append(get_zpj_config(source, groomed=True))

    return regions
