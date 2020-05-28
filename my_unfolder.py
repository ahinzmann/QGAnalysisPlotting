"""
Main class to handle unfolding
"""


from __future__ import print_function, division

import os
import sys
from array import array
import numpy as np
import math
from itertools import chain
import scipy
from scipy import stats
import inspect
import warnings
import pickle
import gzip

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder_plotter import MyUnfolderPlotter


# This doesn't seem to work...sigh
np.set_printoptions(edgeitems=3,infstr='Infinity',
                    linewidth=75, nanstr='nan', precision=8,
                    suppress=False, threshold=1000, formatter=None)

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


# Load my derived class
with open("MyTUnfoldDensity.cpp") as f:
    code = f.read()
    ROOT.gInterpreter.ProcessLine(code)

# Load my derived class
# ROOT.gInterpreter.ProcessLine(".L MyTUnfoldDensity.cc+")
# ROOT.gSystem.Load("MyTUnfoldDensity_cc.so")


class MyUnfolder(ROOT.MyTUnfoldDensity):
    """Main class to handle unfolding input/outputs, all the associated objects

    Derived from MyTUnfoldDensity to get access to protected methods/vars
    """
    # for saving/loading
    _simple_attr = [
        # input & truth histograms
        "input_hist",
        "input_hist_bg_subtracted",
        "input_hist_gen_binning",
        "input_hist_gen_binning_bg_subtracted",
        "hist_truth",
        "hist_mc_reco",
        "hist_mc_reco_bg_subtracted",
        "hist_mc_reco_gen_binning",
        "hist_mc_reco_gen_binning_bg_subtracted",
        # save other matrices
        "rhoij_total",
        "probability_matrix",
        # save error matrices
        "ematrix_input",
        "ematrix_stat_response",
        "ematrix_stat",
        "ematrix_tau",
        "ematrix_total",
        # save folded things
        "folded_unfolded",
        "folded_mc_truth",
        # unfolded
        "unfolded",
        "unfolded_stat_err",
        "unfolded_rsp_err",
        # inverse cov for chi2 tests
        "vyy_inv_tmatrix",
        "vxx_inv_th2",
        "vyy_inv_no_bg_th2",
        # cov matr
        "vyy_no_bg_th2",
    ]

    def __init__(self,
                 response_map,  # 2D GEN-RECO heatmap
                 variable_bin_edges_reco, # 'variable' refers to e.g. ptD, LHA
                 variable_bin_edges_gen, # reco for detector binnig, gen for generator (final) binning
                 variable_name,
                 pt_bin_edges_reco,
                 pt_bin_edges_gen,
                 pt_bin_edges_underflow_reco,
                 pt_bin_edges_underflow_gen,
                 orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                 constraintMode=ROOT.TUnfold.kEConstraintArea,
                 regMode=ROOT.TUnfold.kRegModeCurvature,
                 densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth,
                 distribution='generatordistribution',
                 axisSteering='*[b]'):

        # print("Calling __init__ with args:", locals())

        self.response_map = response_map
        if self.response_map is not None:
            # check no uflow
            if orientation == ROOT.TUnfold.kHistMapOutputHoriz:
                if response_map.GetBinContent(0, 0) != 0 or response_map.GetBinContent(0, 1) != 0:
                    raise RuntimeError("Your response_map has entries in 0th gen bin - this means you've got unintended underflow!")
            elif orientation == ROOT.TUnfold.kHistMapOutputVert:
                if response_map.GetBinContent(0, 1) != 0 or response_map.GetBinContent(1, 0) != 0:
                    raise RuntimeError("Your response_map has entries in 0th gen bin - this means you've got unintended underflow!")
        self.variable_name = variable_name
        self.variable_name_safe = cu.no_space_str(variable_name)

        self.variable_bin_edges_reco = variable_bin_edges_reco
        self.nbins_variable_reco = len(variable_bin_edges_reco)-1 if variable_bin_edges_reco is not None else 0
        self.variable_bin_edges_gen = variable_bin_edges_gen
        self.nbins_variable_gen = len(variable_bin_edges_gen)-1 if variable_bin_edges_gen is not None else 0

        self.pt_bin_edges_reco = pt_bin_edges_reco
        self.nbins_pt_reco = len(pt_bin_edges_reco)-1 if pt_bin_edges_reco is not None else 0
        self.pt_bin_edges_gen = pt_bin_edges_gen
        self.nbins_pt_gen = len(pt_bin_edges_gen)-1 if pt_bin_edges_gen is not None else 0

        self.pt_bin_edges_underflow_reco = pt_bin_edges_underflow_reco
        self.nbins_pt_underflow_reco = len(pt_bin_edges_underflow_reco)-1 if pt_bin_edges_underflow_reco is not None else 0
        self.pt_bin_edges_underflow_gen = pt_bin_edges_underflow_gen
        self.nbins_pt_underflow_gen = len(pt_bin_edges_underflow_gen)-1 if pt_bin_edges_underflow_gen is not None else 0

        # Binning setup here MUST match how it was setup in making the input files, otherwise
        # you will have untold pain and suffering!
        # TODO read in from XML
        var_uf, var_of = False, True
        pt_uf, pt_of = False, True  # handle pt under/over flow ourselves
        self.detector_binning = ROOT.TUnfoldBinning("detector")

        self.detector_distribution_underflow = self.detector_binning.AddBinning("detectordistribution_underflow")
        if self.variable_bin_edges_reco is not None:
            self.detector_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        if self.pt_bin_edges_underflow_reco is not None:
            self.detector_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_reco, self.pt_bin_edges_underflow_reco, False, False)

        self.detector_distribution = self.detector_binning.AddBinning("detectordistribution")
        if self.variable_bin_edges_reco is not None:
            self.detector_distribution.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        if self.pt_bin_edges_reco is not None:
            self.detector_distribution.AddAxis("pt", self.nbins_pt_reco, self.pt_bin_edges_reco, False, pt_of)


        self.generator_binning = ROOT.TUnfoldBinning("generator")

        self.generator_distribution_underflow = self.generator_binning.AddBinning("generatordistribution_underflow")
        if self.variable_bin_edges_gen is not None:
            self.generator_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        if self.pt_bin_edges_underflow_gen is not None:
            self.generator_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_gen, self.pt_bin_edges_underflow_gen, pt_uf, False)

        self.generator_distribution = self.generator_binning.AddBinning("generatordistribution")
        if self.variable_bin_edges_gen is not None:
            self.generator_distribution.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        if self.pt_bin_edges_gen is not None:
            self.generator_distribution.AddAxis("pt", self.nbins_pt_gen, self.pt_bin_edges_gen, False, pt_of)

        self.orientation = orientation
        self.constraintMode = constraintMode
        self.regMode = regMode
        self.densityFlags = densityFlags
        self.distribution = distribution
        self.axisSteering = axisSteering

        tunf_args = [
            self.response_map,
            self.orientation,
            self.regMode,
            self.constraintMode,
            self.densityFlags,
            self.generator_binning,
            self.detector_binning,
            # hmm these take preference over whatever is use for scantau?
            self.distribution,
            self.axisSteering
        ]
        # Ensure all necessary arguments are there, or invoke blank ctor
        # Be careful - some will be 0 intentionally but evaulate False
        # Thus need to do is not None instead
        if all([a is not None for a in tunf_args]):
            ROOT.TUnfoldDensity.__init__(self, *tunf_args)
        else:
            ROOT.TUnfoldDensity.__init__(self)

        # for things like get_probability_matrix()...but doesn't seem to do anything?!
        # this is only used when output_distribution_name = 'generatordistribution': then it makes a TH2
        # since it can map to it
        # Otherwise it makes a TH1
        self.use_axis_binning = False
        # self.probability_ndarray = self.response_matrix_to_probability_array(self.response_map)
        # self.probability_ndarray, _ = cu.th2_to_ndarray(self.get_probability_matrix(), oflow_x=False, oflow_y=False)
        # self.probability_ndarray, _ = None, None

        # hists that will be assigned later
        # TODO: change to properties? although still need to cache somewhere

        # reco, to be unfolded (data or MC)
        self.input_hist = None
        self.input_hist_bg_subtracted = None

        # reco, but using gen binning
        self.input_hist_gen_binning = None
        self.input_hist_gen_binning_bg_subtracted = None

        # fakes
        # for total bg, use get_total_background(),
        # or get_total_background_gen_binning()

        # reco MC
        self.hist_mc_reco = None
        self.hist_mc_reco_bg_subtracted = None

        # reco MC, but using gen binning
        self.hist_mc_reco_gen_binning = None
        self.hist_mc_reco_gen_binning_bg_subtracted = None

        # generator-level MC truth
        self.hist_truth = None  # gen truth

        self.tau = 0  # to be set by user later, via TauScanner or LCurveScanner
        self.backgrounds = {}  # gets filled with subtract_background()
        self.backgrounds_gen_binning = {}  # gets filled with subtract_background_gen_binning()

        self.unfolded = None  # set in get_output(), total error
        self.unfolded_stat_err = None  # set in get_unfolded_with_ematrix_stat()
        self.unfolded_rsp_err = None  # set in get_unfolded_with_ematrix_response()

        self.exp_systs = []  # list of ExpSystematic objects to hold experimental systematics

        # use "generator" for signal + underflow region
        # "generatordistribution" only for ???
        self.output_distribution_name = "generator"

        self.folded_unfolded = None  # set in get_folded_unfolded()
        self.folded_unfolded_tunfold = None  # set in get_folded_unfolded()
        self.folded_mc_truth = None  # set in get_folded_mc_truth()

        # For producing normalised distributions
        self.hist_bin_chopper = HistBinChopper(generator_binning=self.generator_binning.FindNode("generatordistribution"),
                                               detector_binning=self.detector_binning.FindNode("detectordistribution"))

        # For setting/getting various uncerts from HistBinChopper
        self.stat_uncert_name = 'unfolded_stat_err'
        self.stat_ematrix_name = "stat_ematrix"

        self.rsp_uncert_name = 'unfolded_rsp_err'
        self.rsp_ematrix_name = "rsp_ematrix"

        self.scale_uncert_name = "scale_uncert"
        self.scale_uncert_ematrix_name = "scale_uncert_ematrix"

        self.pdf_uncert_name = "pdf_uncert"
        self.pdf_uncert_ematrix_name = "pdf_uncert_ematrix"

        self.total_ematrix_name = "total_ematrix"

    @staticmethod
    def construct_tunfold_binning(variable_bin_edges_reco,
                                  variable_bin_edges_gen,
                                  variable_name,
                                  pt_bin_edges_reco,
                                  pt_bin_edges_gen,
                                  pt_bin_edges_underflow_reco,
                                  pt_bin_edges_underflow_gen):
        """Setup TUnfoldBinning objects

        TODO: integrate this into __init__ properly!
        """
        variable_bin_edges_reco = variable_bin_edges_reco
        nbins_variable_reco = len(variable_bin_edges_reco)-1 if variable_bin_edges_reco is not None else 0
        variable_bin_edges_gen = variable_bin_edges_gen
        nbins_variable_gen = len(variable_bin_edges_gen)-1 if variable_bin_edges_gen is not None else 0

        pt_bin_edges_reco = pt_bin_edges_reco
        nbins_pt_reco = len(pt_bin_edges_reco)-1 if pt_bin_edges_reco is not None else 0
        pt_bin_edges_gen = pt_bin_edges_gen
        nbins_pt_gen = len(pt_bin_edges_gen)-1 if pt_bin_edges_gen is not None else 0

        pt_bin_edges_underflow_reco = pt_bin_edges_underflow_reco
        nbins_pt_underflow_reco = len(pt_bin_edges_underflow_reco)-1 if pt_bin_edges_underflow_reco is not None else 0
        pt_bin_edges_underflow_gen = pt_bin_edges_underflow_gen
        nbins_pt_underflow_gen = len(pt_bin_edges_underflow_gen)-1 if pt_bin_edges_underflow_gen is not None else 0

        # Binning setup here MUST match how it was setup in making the input files, otherwise
        # you will have untold pain and suffering!
        # TODO read in from XML
        var_uf, var_of = False, True
        pt_uf, pt_of = False, True  # handle pt under/over flow ourselves
        detector_binning = ROOT.TUnfoldBinning("detector")

        detector_distribution_underflow = detector_binning.AddBinning("detectordistribution_underflow")
        if variable_bin_edges_reco is not None:
            detector_distribution_underflow.AddAxis(variable_name, nbins_variable_reco, variable_bin_edges_reco, var_uf, var_of)
        if pt_bin_edges_underflow_reco is not None:
            detector_distribution_underflow.AddAxis("pt", nbins_pt_underflow_reco, pt_bin_edges_underflow_reco, False, False)

        detector_distribution = detector_binning.AddBinning("detectordistribution")
        if variable_bin_edges_reco is not None:
            detector_distribution.AddAxis(variable_name, nbins_variable_reco, variable_bin_edges_reco, var_uf, var_of)
        if pt_bin_edges_reco is not None:
            detector_distribution.AddAxis("pt", nbins_pt_reco, pt_bin_edges_reco, False, pt_of)


        generator_binning = ROOT.TUnfoldBinning("generator")

        generator_distribution_underflow = generator_binning.AddBinning("generatordistribution_underflow")
        if variable_bin_edges_gen is not None:
            generator_distribution_underflow.AddAxis(variable_name, nbins_variable_gen, variable_bin_edges_gen, var_uf, var_of)
        if pt_bin_edges_underflow_gen is not None:
            generator_distribution_underflow.AddAxis("pt", nbins_pt_underflow_gen, pt_bin_edges_underflow_gen, pt_uf, False)

        generator_distribution = generator_binning.AddBinning("generatordistribution")
        if variable_bin_edges_gen is not None:
            generator_distribution.AddAxis(variable_name, nbins_variable_gen, variable_bin_edges_gen, var_uf, var_of)
        if pt_bin_edges_gen is not None:
            generator_distribution.AddAxis("pt", nbins_pt_gen, pt_bin_edges_gen, False, pt_of)

        return generator_binning, detector_binning


    def save_binning(self, print_xml=True, txt_filename=None):
        """Save binning scheme to txt and/or print XML to screen"""
        if txt_filename:
            with open(txt_filename, 'w') as of:
                of.write("GEN BINNING\n")
                of.write("--------------------\n")
                for name, region in [
                    ("generator_distribution", self.generator_distribution),
                    ("generator_distribution_underflow", self.generator_distribution_underflow)]:
                    of.write("%s: bins %d - %d\n" % (name, region.GetStartBin(), region.GetEndBin()))
                of.write("\nDETECTOR BINNING\n")
                of.write("--------------------\n")
                for name, region in [
                    ("detector_distribution", self.detector_distribution),
                    ("detector_distribution_underflow", self.detector_distribution_underflow)]:
                    of.write("%s: bins %d - %d\n" % (name, region.GetStartBin(), region.GetEndBin()))

        if print_xml:
            # don't know how to create a ofstream in python :( best we can do is ROOT.cout
            ROOT.TUnfoldBinningXML.ExportXML(self.detector_binning, ROOT.cout, True, False)
            ROOT.TUnfoldBinningXML.ExportXML(self.generator_binning, ROOT.cout, False, True)

    @staticmethod
    def _check_save_to_tfile(tfile, obj, name=None):
        if obj is not None:  # be careful of int = 0!
            # for "simple" types, need to do some trickery
            if isinstance(obj, type(np.array([1.]))):
                if len(obj.shape) > 1:
                    raise RuntimeError("cannot save numpy unless 1D")
                n = ROOT.TVectorD(len(obj))
                for ind, val in enumerate(obj):
                    n[ind] = val
                tfile.WriteTObject(n, name)
            elif isinstance(obj, (bool, int, float)):
                n = ROOT.TVectorD(1)  # the only class that works for these types
                n[0] = obj
                tfile.WriteTObject(n, name)
            elif isinstance(obj, str):
                n = ROOT.TNamed(name, obj)
                tfile.WriteTObject(n, name)
            else:
                if name is None:
                    name = obj.GetName()
                tfile.WriteTObject(obj, name)
        else:
            print("Not saving", name, "as does not exist")

    # def save_to_tfile(self, tfile):
    #     """Save important stuff to TFile/TDirectory"""
    #     self._check_save_to_tfile(tfile, self.detector_binning, "detector_binning")
    #     self._check_save_to_tfile(tfile, self.generator_binning, "generator_binning")
    #     self._check_save_to_tfile(tfile, self.response_map, "response_map")
    #     self._check_save_to_tfile(tfile, self.orientation, "orientation")
    #     self._check_save_to_tfile(tfile, self.constraintMode, "constraintMode")
    #     self._check_save_to_tfile(tfile, self.regMode, "regMode")
    #     self._check_save_to_tfile(tfile, self.densityFlags, "densityFlags")
    #     self._check_save_to_tfile(tfile, self.distribution, "distribution")
    #     self._check_save_to_tfile(tfile, self.axisSteering, "axisSteering")

    #     self._check_save_to_tfile(tfile, self.variable_bin_edges_reco, "variable_bin_edges_reco")
    #     self._check_save_to_tfile(tfile, self.variable_bin_edges_gen, "variable_bin_edges_gen")
    #     self._check_save_to_tfile(tfile, self.variable_name, "variable_name")

    #     self._check_save_to_tfile(tfile, self.pt_bin_edges_reco, "pt_bin_edges_reco")
    #     self._check_save_to_tfile(tfile, self.pt_bin_edges_gen, "pt_bin_edges_gen")

    #     self._check_save_to_tfile(tfile, self.pt_bin_edges_underflow_reco, "pt_bin_edges_underflow_reco")
    #     self._check_save_to_tfile(tfile, self.pt_bin_edges_underflow_gen, "pt_bin_edges_underflow_gen")

    #     # handle most of the simple hists
    #     for name in MyUnfolder._simple_attr:
    #         self._check_save_to_tfile(tfile, getattr(self, name, None), name)

    #     # save all backgrounds (incl fakes)
    #     for name, hist in self.backgrounds.items():
    #         self._check_save_to_tfile(tfile, hist, "background_reco_binning_%s" % cu.no_space_str(name))
    #     for name, hist in self.backgrounds_gen_binning.items():
    #         self._check_save_to_tfile(tfile, hist, "background_gen_binning_%s" % cu.no_space_str(name))

    #     # save systematic response matrices
    #     for name, syst_map in self.syst_maps.items():
    #         self._check_save_to_tfile(tfile, syst_map, "syst_map_%s" % cu.no_space_str(name))

    #     # save systematic error matrices
    #     for name, syst_ematrix in self.syst_ematrices.items():
    #         self._check_save_to_tfile(tfile, syst_ematrix, "syst_ematrix_%s" % cu.no_space_str(name))

    #     # save systematic shifts
    #     for name, syst_shift in self.syst_shifts.items():
    #         self._check_save_to_tfile(tfile, syst_shift, "syst_shift_%s" % cu.no_space_str(name))

    #     # save systematic shifted hists (yes this is a bit wasteful)
    #     for name, syst_shift in self.systs_shifted.items():
    #         self._check_save_to_tfile(tfile, syst_shift, "syst_shifted_unfolded_%s" % cu.no_space_str(name))

    #     # Write the TUnfoldDensity object
    #     tfile.cd()
        # super(ROOT.MyTUnfoldDensity, self).Write()

    def save_unfolded_binned_hists_to_tfile(self, tfile):
        if self.stat_uncert_name not in self.hist_bin_chopper.objects:
            print("Cannot save unfolded_stat_err binned as not in HistBinChopper")
            return

        if "unfolded" not in self.hist_bin_chopper.objects:
            print("Cannot save unfolded binned as not in HistBinChopper")
            return
        try:
            for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
                tfile.WriteTObject(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.stat_uncert_name, ibin_pt, 'generator'), "unfolded_stat_err_norm_divBinWidth_%d" % (ibin_pt))
                tfile.WriteTObject(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded", ibin_pt, 'generator'), "unfolded_norm_divBinWidth_%d" % (ibin_pt))
                tfile.WriteTObject(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width("hist_truth", ibin_pt, 'generator'), "hist_truth_norm_divBinWidth_%d" % (ibin_pt))
                tfile.WriteTObject(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width("alt_hist_truth", ibin_pt, 'generator'), "alt_hist_truth_norm_divBinWidth_%d" % (ibin_pt))
                tfile.WriteTObject(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.stat_ematrix_name, ibin_pt, 'generator'), 'unfolded_stat_ematrix_norm_divBinWidth_%d' % (ibin_pt))
                tfile.WriteTObject(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.rsp_ematrix_name, ibin_pt, 'generator'), 'unfolded_rsp_ematrix_norm_divBinWidth_%d' % (ibin_pt))
                tfile.WriteTObject(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.total_ematrix_name, ibin_pt, 'generator'), 'unfolded_total_ematrix_norm_divBinWidth_%d' % (ibin_pt))
        except KeyError:
            pass

    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        def _del_state(key):
            if key in state:
                del state[key]
        # Remove the large entries from being pickled
        # Normally they can be reconstructed from ROOT objects intstead
        _del_state('_probability_ndarray')
        _del_state('response_map_normed_by_detector_pt')
        _del_state('vyy_inv_ndarray')
        _del_state('vxx_inv_ndarray')
        return state

    # def __new__(cls, *args, **kwargs):
    #     print("__new__ called with:")
    #     print("cls:", cls)
    #     print("args:", args)
    #     print("kwargs:", kwargs)
    #     obj = super(MyUnfolder, cls).__new__(cls)
    #     return obj

    def __reduce__(self):
        # This is the magic method needed to pickle this class
        # For some reason, this class has issue with __new__?
        # We pass the attributes needed for the ctor from __getstate__
        # using the getfullargspec() method
        # If we didn't do that, we'd need to return an empty tuple as the 2nd
        # arg, and allow __init__ to take an empty ctor
        # (i.e. default every arg to None)

        # print("Calling __reduce__")
        argspec = inspect.getargspec(MyUnfolder.__init__)
        # print(argspec)
        state = self.__getstate__()
        # print("state:", state)
        this_args = [state.get(a, None) for a in argspec.args if a != "self"]
        # print("this_args:", this_args)
        return (self.__class__, tuple(this_args), state)

    # SETUP INPUT, BACKGROUNDS
    # --------------------------------------------------------------------------
    def set_input(self,
                  input_hist,
                  input_hist_gen_binning=None,
                  hist_truth=None,
                  hist_mc_reco=None,
                  hist_mc_reco_bg_subtracted=None,
                  hist_mc_reco_gen_binning=None,
                  hist_mc_reco_gen_binning_bg_subtracted=None,
                  bias_factor=0):
        """Set hist to be unfolded, plus other basic hists

        Also allow other args to be passed to TUnfoldSys::SetInput
        """
        self.input_hist = input_hist.Clone()
        self.input_hist_bg_subtracted = input_hist.Clone() if input_hist else None
        self.input_hist_gen_binning = input_hist_gen_binning.Clone() if input_hist_gen_binning else None
        self.input_hist_gen_binning_bg_subtracted = input_hist_gen_binning.Clone() if input_hist_gen_binning else None
        self.hist_truth = hist_truth
        self.hist_mc_reco = hist_mc_reco
        self.hist_mc_reco_bg_subtracted = hist_mc_reco_bg_subtracted
        self.hist_mc_reco_gen_binning = hist_mc_reco_gen_binning
        self.hist_mc_reco_gen_binning_bg_subtracted = hist_mc_reco_gen_binning_bg_subtracted
        self.SetInput(input_hist, bias_factor)

    def subtract_background(self, hist, name, scale=1.0, scale_err=0.0):
        """Subtract background source from input hist"""
        # Save into dict of components - needed? since TUnfoldDensity does this as well
        self.backgrounds[name] = hist.Clone()
        self.backgrounds[name].Scale(scale)
        # Also save total input subtracted
        self.input_hist_bg_subtracted.Add(hist, -1*scale)
        # sanity check that none end up < 0
        for ix in range(1, self.input_hist_bg_subtracted.GetNbinsX()+1):
            val = self.input_hist_bg_subtracted.GetBinContent(ix)
            if val < 0:
                raise ValueError("self.input_hist_bg_subtracted bin %d has contents <0: %g" % (ix, val))
        self.SubtractBackground(hist.Clone(), name, scale, scale_err)

    def subtract_background_gen_binning(self, hist, name, scale=1.0, scale_err=0.0):
        """Subtract background source with gen binning from input hist

        NB doesn't affect TUnfold, only for own book keeping
        """
        # Save into dict of components
        self.backgrounds_gen_binning[name] = hist.Clone()
        self.backgrounds_gen_binning[name].Scale(scale)
        # Also save total input subtracted
        self.input_hist_gen_binning_bg_subtracted.Add(hist, -1*scale)

    def get_total_background(self):
        """Get total cumulative background"""
        total_bg_hist = None
        for name, hist in self.backgrounds.items():
            if total_bg_hist is None:
                total_bg_hist = hist.Clone()
            else:
                total_bg_hist.Add(hist)
        return total_bg_hist

    def get_total_background_gen_binning(self):
        """Get total cumulative background with generator binning"""
        total_bg_hist = None
        for name, hist in self.backgrounds_gen_binning.items():
            if total_bg_hist is None:
                total_bg_hist = hist.Clone()
            else:
                total_bg_hist.Add(hist)
        return total_bg_hist

    def do_unfolding(self, tau):
        """Carry out unfolding with given regularisastion parameter"""
        print(">>> Unfolding with tau =", tau)
        self.tau = tau
        self.DoUnfold(tau)

    def calculate_pt_bin_factors(self, which):
        """Calculate bin factors to account for falling distributions when regularising

        NB done according to signal region - excludes underflow!

        which : str
            Choose which histogram to use for integrals, 'gen' or 'unfolded'
        """

        # Tricky - need to get counts in spectrum, ideally with gen binning
        # FIXME use input_hist instead with detector binning?
        which = which.lower()
        if which not in ['unfolded', 'gen']:
            raise ArgumentError("calculate_pt_bin_factors: 'which' arg should be 'unfolded' or 'gen'")

        hist = None
        if which == 'gen':
            if self.input_hist_gen_binning_bg_subtracted is None:
                raise RuntimeError("Need input_hist_gen_binning_bg_subtracted to be able to do calculate_pt_bin_factors")
            hist = self.input_hist_gen_binning_bg_subtracted

        if which == 'unfolded':
            if self.unfolded is None:
                raise RuntimeError("Need unfolded to be able to do calculate_pt_bin_factors")
            hist = self.unfolded

        # Get integral of 1st pt bin in signal region
        gen_node = self.generator_binning.FindNode('generatordistribution')
        first_var = self.variable_bin_edges_gen[0]
        last_var = self.variable_bin_edges_gen[-1]
        pt_val = self.pt_bin_edges_gen[0]
        start_bin = gen_node.GetGlobalBinNumber(first_var+0.00001, pt_val+0.001)
        end_bin = gen_node.GetGlobalBinNumber(last_var-0.00001, pt_val+0.001)
        first_bin_integral = hist.Integral(start_bin, end_bin)  # ROOTs integral is inclusive of last bin

        bin_factors = {}
        # Add 1s for the 1st pt bin
        for var in self.variable_bin_edges_gen[:-1]:
            this_bin = gen_node.GetGlobalBinNumber(var+0.00001, pt_val+0.001)
            bin_factors[this_bin] = 1

        # Iterate through pt bins, figure out integral, scale according to first bin
        for pt_val in self.pt_bin_edges_gen[1:-1]:
            start_bin = gen_node.GetGlobalBinNumber(first_var+0.00001, pt_val+0.001)
            end_bin = gen_node.GetGlobalBinNumber(last_var-0.00001, pt_val+0.001)
            integral = hist.Integral(start_bin, end_bin)
            sf = first_bin_integral / integral

            # Store bin factor for each lambda bin
            for var in self.variable_bin_edges_gen[:-1]:
                this_bin = gen_node.GetGlobalBinNumber(var+0.00001, pt_val+0.001)
                bin_factors[this_bin] = sf

        return bin_factors

    def get_gen_bin_widths(self):
        results = {}
        # loop through the real axes, convert to global bin number, store width
        # This is because TUnfoldBinning has no *public* method to go the other
        # way...ToAxisBins() is *protected* FFS
        for lambda_ind, lambda_var in enumerate(self.variable_bin_edges_gen[:-1]):
            # underflow region
            for pt_ind, pt_val in enumerate(self.pt_bin_edges_underflow_gen[:-1]):
                global_bin = self.generator_distribution_underflow.GetGlobalBinNumber(lambda_var+0.00000001, pt_val+0.0000001)
                lambda_width = self.variable_bin_edges_gen[lambda_ind+1] - self.variable_bin_edges_gen[lambda_ind]
                pt_width = self.pt_bin_edges_underflow_gen[pt_ind+1] - self.pt_bin_edges_underflow_gen[pt_ind]
                results[global_bin] = (lambda_width, pt_width)
            # signal region
            for pt_ind, pt_val in enumerate(self.pt_bin_edges_gen[:-1]):
                global_bin = self.generator_distribution.GetGlobalBinNumber(lambda_var+0.00000001, pt_val+0.0000001)
                lambda_width = self.variable_bin_edges_gen[lambda_ind+1] - self.variable_bin_edges_gen[lambda_ind]
                pt_width = self.pt_bin_edges_gen[pt_ind+1] - self.pt_bin_edges_gen[pt_ind]
                results[global_bin] = (lambda_width, pt_width)
        return results

    # HANDLE SYSTEMATIC UNCERTAINTIES
    # --------------------------------------------------------------------------
    def get_all_exp_syst_labels(self):
        return [x.label for x in self.exp_systs]

    def get_exp_syst(self, label):
        items = [x for x in self.exp_systs if x.label == label]
        if len(items) == 0:
            raise ValueError("Found no exp systematic with label '%s', only have: %s" % (label, [x.label for x in self.exp_systs]))
        if len(items) > 1:
            raise ValueError("Found >1 exp systematic with label '%s': %s" % (label, [x.label for x in items]))
        return items[0]

    def add_sys_error(self, map_syst, label, mode):
        """Add systematic error via response map, arguments as per AddSysError()"""
        this_syst = ExpSystematic(label=label, syst_map=map_syst)
        self.AddSysError(map_syst, label, self.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)
        self.exp_systs.append(this_syst)

    def get_delta_sys_shift(self, syst_label):
        """Get shift in result due to a particular systeamtic

        Label must be same as used to add it in add_sys_error()
        """
        # TODO: syst_label -> ExpSystematic obj
        this_syst = self.get_exp_syst(syst_label)
        if this_syst.syst_shift is None:
            hist = self.GetDeltaSysSource(this_syst.label,
                                          this_syst.syst_shift_label,  # name given to hist obj
                                          "",
                                          self.output_distribution_name, # must be the same as what's used in get_output
                                          self.axisSteering,
                                          self.use_axis_binning)
            this_syst.syst_shift = hist  # cache shifts
        return this_syst.syst_shift

    def get_syst_shifted_hist(self, syst_label, unfolded=None):
        """
        Can specify unfolded hist, default is the one with all errors
        """
        # TODO: syst_label -> ExpSystematic obj
        this_syst = self.get_exp_syst(syst_label)
        if this_syst.syst_shifted is None:
            hist_shift = self.get_delta_sys_shift(syst_label).Clone(this_syst.syst_shifted_label)
            unfolded = unfolded or self.get_unfolded_with_ematrix_stat()
            hist_shift.Add(unfolded)  # TODO what about errors?
            this_syst.syst_shifted = hist_shift
        return this_syst.syst_shifted

    # POST-UNFOLDING FUNCTIONS
    # --------------------------------------------------------------------------
    def get_output(self, hist_name='unfolded', max_chi2_ndf=1000):
        """Get 1D unfolded histogram covering all bins"""
        if getattr(self, 'unfolded', None) is None:
            print("Ndf:", self.GetNdf())
            self.Ndf = self.GetNdf()
            print("Npar:", self.GetNpar())
            self.Npar = self.GetNpar()
            print("chi2sys:", self.GetChi2Sys())
            self.chi2sys = self.GetChi2Sys()
            print("chi2A:", self.GetChi2A())
            self.chi2A = self.GetChi2A()
            print("chi2L:", self.GetChi2L())
            self.chi2L = self.GetChi2L()

            self.unfolded = self.GetOutput(hist_name, "", self.output_distribution_name, "*[]", self.use_axis_binning)
            print("self.unfolded is a", type(self.unfolded), "and has", self.unfolded.GetNbinsX(), "x bins")

            # Sanity check incase unfolding went really wrong - e.g. if rank
            # of matrix E is way smaller than what is expected
            # May mean SetEpsMatrix() is needed to help inversion
            if max_chi2_ndf > 0 and self.chi2sys / self.Ndf > max_chi2_ndf:
                raise RuntimeError("chi2sys / Ndf > %d - unfolding is rubbish! Maybe SetEpsMatrix()? if you get 'rank of matrix E X expect Y' warning" % (max_chi2_ndf))

        return self.unfolded

    def _post_process(self):
        """Do some standard things & store various things that are done after unfolding"""
        self.get_ematrix_input()
        self.get_ematrix_stat_response()
        self.get_ematrix_stat()
        self.get_ematrix_tau()
        self.get_ematrix_total()
        self.get_rhoij_total()
        self.get_probability_matrix()
        self.update_unfolded_with_ematrix_total()
        self.get_unfolded_with_ematrix_stat()
        self.get_unfolded_with_ematrix_response()
        self.get_folded_unfolded()
        self.get_folded_mc_truth()
        for exp_syst in self.exp_systs:
            # setup all the internal maps
            self.get_delta_sys_shift(exp_syst.label)
            self.get_ematrix_syst(exp_syst.label)
            self.get_syst_shifted_hist(exp_syst.label)

        self.hist_bin_chopper.add_obj('hist_truth', self.hist_truth)
        self.hist_bin_chopper.add_obj('unfolded', self.get_output())
        self.hist_bin_chopper.add_obj(self.stat_uncert_name, self.get_unfolded_with_ematrix_stat())
        self.hist_bin_chopper.add_obj(self.rsp_uncert_name, self.get_unfolded_with_ematrix_response())

    @staticmethod
    def make_hist_from_diagonal_errors(h2d, do_sqrt=True):
        """Make 1D hist, with errors set to diagonal elements from h2d

        Can be TH2 or numpy.ndarray, cos we have to use both
        Yes that is majorly wack
        """
        if isinstance(h2d, ROOT.TH2):
            nbins = h2d.GetNbinsX()
            hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, 0, nbins)
            for i in range(1, nbins+1):
                err = h2d.GetBinContent(i, i)
                if do_sqrt and err > 0:
                    err = math.sqrt(err)
                hnew.SetBinError(i, err)
            return hnew
        elif isinstance(h2d, np.ndarray):
            nbins = h2d.shape[0]
            hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, 0, nbins)
            for i in range(1, nbins+1):
                err = h2d[i-1, i-1]
                if do_sqrt and err > 0:
                    err = math.sqrt(err)
                hnew.SetBinError(i, err)
            return hnew

    @staticmethod
    def make_diag_cov_hist_from_errors(h1d, do_squaring=True, inverse=False):
        """Make diagonal TH2 from errors on TH1.

        Assumes off-diag = 0.

        Can also do inverse by doing 1/(err^2)
        """
        nbins = h1d.GetNbinsX()
        bin_edges = array('d', [h1d.GetBinLowEdge(i) for i in range(1, nbins+2)])
        h = ROOT.TH2D(cu.get_unique_str(), "", nbins, bin_edges, nbins, bin_edges)
        for i in range(1, nbins+1):
            err = h1d.GetBinError(i)
            if do_squaring:
                err *= err
            if inverse:
                if err != 0:
                    err = 1/err
                else:
                    err = 0
            h.SetBinContent(i, i, err)
        return h

    @staticmethod
    def update_hist_bin_error(h_orig, h_to_be_updated):
        """Change the errors in h_to_be_updated to those from h_orig"""
        if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
            raise RuntimeError("Need same # x bins, %d vs %s" % (h_orig.GetNbinsX(), h_to_be_updated.GetNbinsX()))
        for i in range(0, h_orig.GetNbinsX()+2):
            h_to_be_updated.SetBinError(i, h_orig.GetBinError(i))

    def update_unfolded_with_ematrix_total(self):
        """Update unfolded hist with total errors from total error matrix"""
        error_total_1d = self.make_hist_from_diagonal_errors(self.get_ematrix_total(), do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
        self.update_hist_bin_error(h_orig=error_total_1d, h_to_be_updated=self.get_output())

    def get_unfolded_with_ematrix_stat(self):
        """Make copy of unfolded, but only stat errors"""
        if getattr(self, 'unfolded_stat_err', None) is None:
            error_1d = self.make_hist_from_diagonal_errors(self.get_ematrix_stat(), do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
            self.unfolded_stat_err = self.get_output().Clone("unfolded_stat_err")
            self.update_hist_bin_error(h_orig=error_1d, h_to_be_updated=self.unfolded_stat_err)
        return self.unfolded_stat_err

    def get_unfolded_with_ematrix_response(self):
        """Create unfolded with error bars from response matrix uncertainties"""
        if getattr(self, 'unfolded_rsp_err', None) is None:
            error_stat_response = self.make_hist_from_diagonal_errors(self.get_ematrix_stat_response(), do_sqrt=True) # note that bin contents need to be correct, otherwise won't normalise correctly
            self.unfolded_rsp_err = self.get_output().Clone("unfolded_rsp_err")
            self.update_hist_bin_error(h_orig=error_stat_response, h_to_be_updated=self.unfolded_rsp_err)
        return self.unfolded_rsp_err

    def get_bias_vector(self):
        if getattr(self, "bias_vector", None) is None:
            self.bias_vector = self.GetBias("bias_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.bias_vector

    def get_probability_matrix(self):
        if getattr(self, "probability_matrix", None) is None:
            self.probability_matrix = self.GetProbabilityMatrix("prob_matrix_"+cu.get_unique_str(), "", self.use_axis_binning)
        return self.probability_matrix

    @property
    def probability_ndarray(self):
        cached_attr_name = '_probability_ndarray'
        if not hasattr(self, cached_attr_name):
            arr, _ = cu.th2_to_ndarray(self.get_probability_matrix())
            setattr(self, cached_attr_name, arr)
        return getattr(self, cached_attr_name)

    def get_rhoij_total(self):
        if getattr(self, "rhoij_total", None) is None:
            self.rhoij_total = self.GetRhoIJtotal("rhoij_total_"+cu.get_unique_str(), "", self.output_distribution_name, "*[]", self.use_axis_binning)
        return self.rhoij_total

    # LOTS OF COVARIANCE MATRIX FUNCTIONS
    # --------------------------------------------------------------------------
    def get_ematrix_input(self):
        """Get error matrix due to statistics from thing being unfolded"""
        if getattr(self, "ematrix_input", None) is None:
            self.ematrix_input = self.GetEmatrixInput("ematrix_input_"+cu.get_unique_str(), "", self.output_distribution_name, "*[]", self.use_axis_binning)
        return self.ematrix_input

    def get_ematrix_stat_response(self):
        """Statistical uncertainty error matrix from response matrix, should be considered a systematic uncert"""
        if getattr(self, "ematrix_stat_response", None) is None:
            self.ematrix_stat_response = self.GetEmatrixSysUncorr("ematrix_stat_response_"+cu.get_unique_str(), "", self.output_distribution_name, "*[]", self.use_axis_binning)
        return self.ematrix_stat_response

    @property
    def ematrix_stat_response_ndarray(self):
        cached_attr_name = '_ematrix_stat_response_ndarray'
        if not hasattr(self, cached_attr_name):
            arr, _ = cu.th2_to_ndarray(self.get_ematrix_stat_response())
            setattr(self, cached_attr_name, arr)
        return getattr(self, cached_attr_name)

    def get_ematrix_total(self):
        """Total error matrix, from stat+systs"""
        if getattr(self, "ematrix_total", None) is None:
            self.ematrix_total = self.GetEmatrixTotal("ematrix_total_"+cu.get_unique_str(), "", self.output_distribution_name, "*[]", self.use_axis_binning)
            print("ematrix_total is:", type(self.ematrix_total), "with #xbins:", self.ematrix_total.GetNbinsX())
        return self.ematrix_total

    @property
    def ematrix_total_ndarray(self):
        cached_attr_name = '_ematrix_total_ndarray'
        if not hasattr(self, cached_attr_name):
            arr, _ = cu.th2_to_ndarray(self.get_ematrix_total())
            setattr(self, cached_attr_name, arr)
        return getattr(self, cached_attr_name)

    def get_ematrix_stat(self):
        """Get total statitical error matrix (from input being unfolded + background sources, including fakes)"""
        if getattr(self, 'ematrix_stat', None) is None:
            # Have to manually create hist first, awkward
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            # Gnahhhhhhh
            self.ematrix_stat = this_binning.CreateErrorMatrixHistogram("ematrix_stat_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
            self.GetEmatrix(self.ematrix_stat)
        return self.ematrix_stat

    @property
    def ematrix_stat_ndarray(self):
        cached_attr_name = '_ematrix_stat_ndarray'
        if not hasattr(self, cached_attr_name):
            arr, _ = cu.th2_to_ndarray(self.get_ematrix_stat())
            setattr(self, cached_attr_name, arr)
        return getattr(self, cached_attr_name)

    def get_ematrix_syst(self, syst_label):
        """Get error matrix from a systematic source"""

        # TODO: syst_label -> ExpSystematic obj

        this_syst = self.get_exp_syst(syst_label)

        if this_syst.syst_ematrix is None:
            # query the variables inside TUnfolder itself
            syst_source_names = [x.GetName() for x in self.GetSysSources()]
            if syst_label in syst_source_names:
                # Have to manually create hist first, awkward
                this_binning = self.generator_binning.FindNode('generator')
                # I cannot figure out how to make the int** object for bin_map
                # So we are trusting that the default args for title and axisSteering are correct
                # Gnahhhhhhh
                hist = this_binning.CreateErrorMatrixHistogram(this_syst.syst_ematrix_label, self.use_axis_binning) #, bin_map, "", "*[]")
                self.GetEmatrixSysSource(hist, syst_label)
                this_syst.syst_ematrix = hist
            else:
                # TODO: make it ourself from deltaX
                pass
        return this_syst.syst_ematrix

    # property or getter?
    @property
    def ematrix_syst_ndarray(self, syst_label):
        cached_attr_name = '_ematrix_syst_%s_ndarray' % (cu.no_space_str(syst_label))
        if not hasattr(self, cached_attr_name):
            arr, _ = cu.th2_to_ndarray(self.get_ematrix_syst(syst_label))
            setattr(self, cached_attr_name, arr)
        return getattr(self, cached_attr_name)

    def get_ematrix_tau(self):
        """Get error matrix due to regularisation uncertainty"""
        if getattr(self, 'ematrix_tau', None) is None:
            # Have to manually create hist first, awkward
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            # Gnahhhhhhh
            self.ematrix_tau = this_binning.CreateErrorMatrixHistogram("ematrix_tau_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
            self.GetEmatrixSysTau(self.ematrix_tau)
        return self.ematrix_tau

    def get_ematrix_total_inv(self):
        """Total error matrix inverted, from stat+systs"""
        return self.InvertMSparseSymmPos(self.GetSummedErrorMatrixXX(), False)  # TMatrixDSparse

    def get_vyy_inv_ndarray(self):
        """Get inverse of V_yy (ie input). Note done after BG-subtraction"""
        if getattr(self, 'vyy_inv_ndarray', None) is None:
            if getattr(self, 'vyy_inv_tmatrix', None) is None:
                self.vyy_inv_tmatrix = self.GetVyyInv()
            self.vyy_inv_ndarray = cu.tmatrixdsparse_to_ndarray(self.vyy_inv_tmatrix)
        return self.vyy_inv_ndarray

    def get_vyy_no_bg_th2(self):
        """Same as get_vyy_inv_ndarray() but before BG-subtraction"""
        if getattr(self, 'vyy_no_bg_th2', None) is None:
            self.vyy_no_bg_th2 = self.make_diag_cov_hist_from_errors(self.input_hist, inverse=False)
        return self.vyy_no_bg_th2

    def get_vyy_no_bg_ndarray(self):
        if getattr(self, 'vyy_no_bg_ndarray', None) is None:
            self.vyy_no_bg_ndarray, _ = cu.th2_to_ndarray(self.get_vyy_no_bg_th2())
        return self.vyy_no_bg_ndarray

    def get_vyy_inv_no_bg_th2(self):
        if getattr(self, 'vyy_inv_no_bg_th2', None) is None:
            self.vyy_inv_no_bg_th2 = self.make_diag_cov_hist_from_errors(self.input_hist, inverse=True)
            # this_binning = self.detector_binning.FindNode('detector')
            # self.vyy_inv_no_bg_th2 = this_binning.CreateErrorMatrixHistogram("ematrix_vyyinv_no_bg_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
            # for i in range(1, self.input_hist.GetNbinsX()+1):
            #     bin_err = self.input_hist.GetBinError(i)
            #     new_err = 1./(bin_err*bin_err) if bin_err != 0 else 0
            #     self.vyy_inv_no_bg_th2.SetBinContent(i, i, new_err)
        return self.vyy_inv_no_bg_th2

    def get_vyy_inv_no_bg_ndarray(self):
        if getattr(self, 'vyy_inv_no_bg_ndarray', None) is None:
            self.vyy_inv_no_bg_ndarray, _ = cu.th2_to_ndarray(self.get_vyy_inv_no_bg_th2())
        return self.vyy_inv_no_bg_ndarray

    def get_vxx_inv_th2(self):
        if getattr(self, 'vxx_inv_th2', None) is None:
            # Have to manually create hist first, awkward
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            # Gnahhhhhhh
            self.vxx_inv_th2 = this_binning.CreateErrorMatrixHistogram("ematrix_vxxinv_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
            self.ErrorMatrixToHist(self.vxx_inv_th2, self.GetVxxInv())
        return self.vxx_inv_th2

    def get_vxx_inv_ndarray(self):
        if getattr(self, 'vxx_inv_ndarray', None) is None:
            self.vxx_inv_ndarray, _ = cu.th2_to_ndarray(self.get_vxx_inv_th2())
        return self.vxx_inv_ndarray

    def get_vyy_ndarray(self):
        """Get V_yy (ie input). Note done after BG-subtraction"""
        if getattr(self, 'vyy_ndarray', None) is None:
            if getattr(self, 'vyy_tmatrix', None) is None:
                self.vyy_tmatrix = self.GetVyy()
            self.vyy_ndarray = cu.tmatrixdsparse_to_ndarray(self.vyy_tmatrix)
        return self.vyy_ndarray

    @staticmethod
    def calculate_singular_max_min(matrix):
        """Calculate max & min singular values condition number as per StatsComm guidelines

        These are found using TDecompSVD.
        (we ignore the builtin condition() method as it calculates it differently)
        """
        if matrix.GetNcols() > matrix.GetNrows():
            raise RuntimeError("Condition number only for matrix where # rows >= # cols")

        svd = ROOT.TDecompSVD(matrix)
        sig = svd.GetSig()  # by construction, ordered descending
        sigma_max = sig[0]
        sigma_min = max(0, sig[sig.GetNrows()-1])
        print("sigma_max:", sigma_max, "sigma_min:", sigma_min)
        if sigma_min == 0:
            print("sigma_min > 0:", min([x for x in sig if x>0]))
        return sigma_max, sigma_min

    def print_condition_number(self):
        """Store & print response matrix condition number and some advice

        Defined as sigma_max / max(0, sigma_min), where sigma_{max/min} are the
        largest/smallest singular values.
        These are also stored for later usage if needed (since expensive to calc)
        """
        if getattr(self, 'condition_number', None) is None:
            sigma_max, sigma_min = self.calculate_singular_max_min(cu.th2_to_tmatrixd(self.get_probability_matrix()))
            if sigma_min == 0:
                # avoid DivisionError
                print("Minmum singular value = 0, condition number = Infinity")
                num = np.inf
            else:
                num = sigma_max / sigma_min
            self.sigma_max = sigma_max
            self.sigma_min = sigma_min
            self.condition_number = num

        print("Condition number:", self.condition_number)
        if self.condition_number < 50:
            print(" - You probably shouldn't regularize this")
        elif self.condition_number > 1E5:
            print(" - You probably should regularize this")
        else:
            print(" - You probably should look into regularization")

    def get_response_normed_by_detector_pt(self):
        if getattr(self, 'response_map_normed_by_detector_pt', None) is None:
            normed_response_map = self.response_map.Clone("response_map_normed_by_detector_pt")
            sums = []
            bins = list(self.pt_bin_edges_underflow_reco) + list(self.pt_bin_edges_reco)[1:]
            # print(bins)
            # First get the sum over all bins corresponding to this reco pT bin
            # Means summing over a number of reco lambda bins, and all generator bins
            for ibin_pt, (pt, pt_next) in enumerate(zip(bins[:-1], bins[1:])):
                # print(ibin_pt, pt, pt_next)
                this_sum = 0
                var = self.variable_bin_edges_reco[0] * 1.00000001
                # urgh this is horrible, but crashes if you use self.detector_binning.GetGlobalBinNumber() whyyyy
                # need separate binning obj for each value, else it misses bins
                binning = self.detector_distribution_underflow if pt in self.pt_bin_edges_underflow_reco else self.detector_distribution
                binning_next = self.detector_distribution_underflow if pt_next in self.pt_bin_edges_underflow_reco else self.detector_distribution
                global_bin_reco = binning.GetGlobalBinNumber(var, pt*1.00001)
                global_bin_reco_next = binning_next.GetGlobalBinNumber(var, pt_next*1.00001)
                # print(ibin_pt, pt, pt_next, global_bin_reco, global_bin_reco_next)
                for i_reco in range(global_bin_reco, global_bin_reco_next):
                    for global_bin_gen in range(self.response_map.GetNbinsX()+1):
                        this_sum += self.response_map.GetBinContent(global_bin_gen, i_reco)
                sums.append(this_sum)
            # print(sums)

            # Now set the new bin contents and error bars by scaling using these sums
            for ibin_pt, (pt, pt_next) in enumerate(zip(bins[:-1], bins[1:])):
                var = self.variable_bin_edges_reco[0] * 1.00000001
                # urgh this is horrible, but crashes if you use self.detector_binning.GetGlobalBinNumber() whyyyy
                binning = self.detector_distribution_underflow if pt in self.pt_bin_edges_underflow_reco else self.detector_distribution
                binning_next = self.detector_distribution_underflow if pt_next in self.pt_bin_edges_underflow_reco else self.detector_distribution
                global_bin_reco = binning.GetGlobalBinNumber(var, pt*1.00001)
                global_bin_reco_next = binning_next.GetGlobalBinNumber(var, pt_next*1.00001)
                for i_reco in range(global_bin_reco, global_bin_reco_next):
                    for global_bin_gen in range(self.response_map.GetNbinsX()+1):
                        factor = 1./sums[ibin_pt] if sums[ibin_pt] != 0 else 0
                        # print("new content", global_bin_gen, i_reco, self.response_map.GetBinContent(global_bin_gen, i_reco) * factor)
                        normed_response_map.SetBinContent(global_bin_gen, i_reco, self.response_map.GetBinContent(global_bin_gen, i_reco) * factor)
                        normed_response_map.SetBinError(global_bin_gen, i_reco, self.response_map.GetBinError(global_bin_gen, i_reco) * factor)

            self.response_map_normed_by_detector_pt = normed_response_map
        return self.response_map_normed_by_detector_pt

    # METHODS FOR FORWARD-FOLDING & CHI2 TESTS
    # --------------------------------------------------------------------------
    def get_folded_unfolded(self):
        # don't use getfoldedoutput, because it doesn't have the updated errors from the total error matrix
        # so we'll have to do it ourselves
        # 1. Make unfolded hist into TVector/TMatrix

        # 2. Make response 2d hist into matrix

        # 3. Multiply the two, convert to TH1

        # Get the TUnfold one for reference, although its errors will be wrong
        if getattr(self, 'folded_unfolded_tunfold', None) is None:
            self.folded_unfolded_tunfold = self.GetFoldedOutput("folded_unfolded_tunf")

        if getattr(self, 'folded_unfolded', None) is None:
            oflow = False

            # Get unfolded results as array
            unfolded_vector, _ = cu.th1_to_ndarray(self.get_output(), oflow)

            # Multiply
            # Note that we need to transpose from row vec to column vec
            folded_vec = self.probability_ndarray.dot(unfolded_vector.T)

            # Convert vector to TH1
            self.folded_unfolded = cu.ndarray_to_th1(folded_vec.T, has_oflow_x=oflow)

            # Error propagation: if y = Ax, with covariance matrices Vyy and Vxx,
            # respectively, then Vyy = (A*Vxx)*A^T
            unfolded_covariance_matrix, _ = cu.th2_to_ndarray((self.get_ematrix_total()), oflow_x=oflow, oflow_y=oflow)
            result = self.probability_ndarray.dot(unfolded_covariance_matrix)
            folded_covariance = result.dot(self.probability_ndarray.T)
            folded_errors = self.make_hist_from_diagonal_errors(folded_covariance)
            self.update_hist_bin_error(h_orig=folded_errors, h_to_be_updated=self.folded_unfolded)

        return self.folded_unfolded

    def fold_generator_level(self, hist_truth):
        oflow = False
        # Convert hist to vector
        gen_vec, gen_vec_err = cu.th1_to_ndarray(hist_truth, oflow_x=oflow)

        # Multiply
        # Note that we need to transpose from row vec to column vec
        folded_vec = self.probability_ndarray.dot(gen_vec.T)

        # Convert vector to TH1
        folded_mc_truth = cu.ndarray_to_th1(folded_vec.T, has_oflow_x=oflow)

        # Error propagation: if y = Ax, with covariance matrices Vyy and Vxx,
        # respectively, then Vyy = (A*Vxx)*A^T
        vxx, _ = cu.th2_to_ndarray(self.make_diag_cov_hist_from_errors(hist_truth, inverse=False), oflow)
        result = self.probability_ndarray.dot(vxx)
        folded_covariance = result.dot(self.probability_ndarray.T)
        folded_errors = self.make_hist_from_diagonal_errors(folded_covariance)
        self.update_hist_bin_error(h_orig=folded_errors, h_to_be_updated=folded_mc_truth)
        return folded_mc_truth

    def get_folded_mc_truth(self):
        """Get response_matrix * MC truth"""
        if getattr(self, 'folded_mc_truth', None) is None:
            self.folded_mc_truth = self.fold_generator_level(self.hist_truth)
        return self.folded_mc_truth

    def do_numpy_comparison(self, output_dir):
        """Do unregularised unfolding and simple matrix inversion with numpy to compare"""
        vyyinv_ndarray = self.get_vyy_inv_ndarray()

        rsp_inv_ndarray = np.linalg.pinv(self.probability_ndarray)

        y_ndarray = np.zeros(shape=(self.GetNy(),1))
        y = self.GetY()
        for i in range(self.GetNy()):
            y_ndarray[i,0] = y(i, 0)

        # calculate my own result simple inversion
        # nb this is really dodgy as non-square matrices don't have an inverse technically...
        result = rsp_inv_ndarray @ y_ndarray
        # print(result.shape)
        result_hist = cu.ndarray_to_th1(result.T)

        # calculate full non-regularised result
        Einv = self.probability_ndarray.T @ vyyinv_ndarray @ self.probability_ndarray
        E = np.linalg.pinv(Einv, rcond=1E-160)  # can't use .inv as singular matrix
        rhs = self.probability_ndarray.T @ vyyinv_ndarray @ y_ndarray
        proper_x = E @ rhs

        E_hist = cu.ndarray_to_th2(E)
        E_hist.SetTitle("E = (A^{T}V^{-1}_{yy}A)^{-1}")
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        E_hist.Draw("COLZ")
        cu.set_french_flag_palette()
        canv.SetRightMargin(0.2)
        canv.SaveAs(os.path.join(output_dir, "E.pdf"))
        canv.Clear()

        E_inv_hist = cu.ndarray_to_th2(Einv)
        E_inv_hist.SetTitle("E^{-1} = A^{T}V^{-1}_{yy}A")
        E_inv_hist.Draw("COLZ")
        ROOT.gStyle.SetPalette(ROOT.kBird)
        canv.SetRightMargin(0.2)
        canv.SetLogz()
        E_inv_hist.SetMinimum(1E-10)
        canv.SaveAs(os.path.join(output_dir, "E_inv.pdf"))
        canv.Clear()

        rhs_hist = cu.ndarray_to_th1(rhs.T)
        canv.SetRightMargin(0.12)
        rhs_hist.SetTitle("A^{T}V^{-1}_{yy}y;Generator bin;N")
        rhs_hist.Draw("HIST")

        canv.SaveAs(os.path.join(output_dir, "rhs.pdf"))

        proper_x_hist = cu.ndarray_to_th1(proper_x.T)

        unfolded_hist = self.get_output().Clone("bah")
        cu.remove_th1_errors(unfolded_hist)
        cu.remove_th1_errors(result_hist)
        cu.remove_th1_errors(proper_x_hist)

        conts = [
            Contribution(unfolded_hist, line_color=ROOT.kBlue, marker_color=ROOT.kBlue, label='TUnfold'),
            # Contribution(result_hist, line_color=ROOT.kGreen, marker_color=ROOT.kGreen, label='Simple numpy inversion', subplot=unfolded_hist),
            Contribution(proper_x_hist, line_color=ROOT.kRed, marker_color=ROOT.kRed, label='Full numpy inversion', subplot=unfolded_hist, line_style=2),
        ]
        # title = "%s\n%s region" % (region['jet_algo'], region['label'])
        plot = Plot(conts, what='hist',
                    xtitle='Generator bin',
                    ytitle='N',
                    # title=title,
                    # ylim=(1E-3, 1E8),
                    ylim=(-1E3, 1E4),  # focus on small values, +ve and -ve
                    subplot_type='ratio',
                    subplot_limits=(0.5, 1.5),
                    subplot_title="* / TUnfold")
        plot.default_canvas_size = (800, 600)
        plot.plot("HISTE NOSTACK")
        plot.legend.SetY1NDC(0.8)
        plot.legend.SetX1NDC(0.65)
        plot.legend.SetX2NDC(0.9)
        # plot.set_logy(do_more_labels=False)
        plot.save(os.path.join(output_dir, 'tunfold_vs_numpy.pdf'))

    def calculate_chi2(self, one_hist, other_hist, cov_inv_matrix, detector_space=True, ignore_underflow_bins=True, debugging_dir=None):
        one_vec, _ = cu.th1_to_ndarray(one_hist, False)
        other_vec, _ = cu.th1_to_ndarray(other_hist, False)
        delta = one_vec - other_vec

        first_signal_bin = 1
        if ignore_underflow_bins:
            first_signal_bin = self.detector_distribution.GetStartBin() if detector_space else self.generator_distribution.GetStartBin()
            delta[0][:first_signal_bin-1] = 0. # subtract 1 as numpy indices start at 0, hists start at 1

        if isinstance(cov_inv_matrix, ROOT.TH2):
            v_inv, _ = cu.th2_to_ndarray(cov_inv_matrix)
        else:
            v_inv = cov_inv_matrix

        inter = v_inv.dot(delta.T)
        chi2 = delta.dot(inter)[0][0]
        ndof = len(delta[0][first_signal_bin-1:])
        p = 1-scipy.stats.chi2.cdf(chi2, int(ndof))

        if debugging_dir:
            # print some debugging plots
            debug_plotter = MyUnfolderPlotter(self, False)
            # 1D inputs
            entries = [
                Contribution(one_hist, label='one_hist', line_color=ROOT.kBlack),
                Contribution(other_hist, label='other_hist', line_color=ROOT.kRed, line_style=2)
            ]
            plot = Plot(entries, what='hist', has_data=False,)
            plot.default_canvas_size = (800, 600)
            plot.plot("NOSTACK HISTE")
            l,t = debug_plotter.draw_pt_binning_lines(plot, which='reco' if detector_space else 'gen', axis='x')
            plot.save(os.path.join(debugging_dir, 'one_other_hists.pdf'))

            # Delta, with missing bins if necessary
            delta_hist = cu.ndarray_to_th1(delta)
            entries = [
                Contribution(delta_hist)
            ]
            plot = Plot(entries,
                        what='hist',
                        xtitle='%s bin' % ('Detector' if detector_space else 'Generator'),
                        ytitle='one_hist - other_hist',
                        has_data=False,
                        )
            plot.default_canvas_size = (800, 600)
            plot.plot("NOSTACK HISTE")
            l,t = debug_plotter.draw_pt_binning_lines(plot, which='reco' if detector_space else 'gen', axis='x')
            plot.save(os.path.join(debugging_dir, 'delta.pdf'))

            # Covariance matrix
            canv = ROOT.TCanvas("c", "Inverse covariance matrix V^{-1}", 800, 600)
            obj = cov_inv_matrix
            if not isinstance(cov_inv_matrix, ROOT.TH2):
                obj = cu.ndarray_to_th2(v_inv)
            obj.Draw("COLZ")
            canv.SetLeftMargin(0.15)
            canv.SetRightMargin(0.18)
            canv.SetLogz(1)
            cov_min = obj.GetMinimum(1E-20) / 10
            obj.SetMinimum(cov_min)
            cov_max = obj.GetMaximum() * 5
            obj.SetMaximum(cov_max)

            l,t = debug_plotter.draw_pt_binning_lines(obj, which='reco' if detector_space else 'gen', axis='x', do_underflow=True, do_labels_inside=False, do_labels_outside=True)
            l2,t2 = debug_plotter.draw_pt_binning_lines(obj, which='reco' if detector_space else 'gen', axis='y', do_underflow=True, do_labels_inside=False, do_labels_outside=True)
            canv.SaveAs(os.path.join(debugging_dir, 'cov_inv_matrix.pdf'))


            # Components of delta * V_inv * delta before summing
            components = delta * inter.T
            components_hist = cu.ndarray_to_th1(components)
            y_max = components_hist.GetMaximum() * 1.2
            entries = [
                Contribution(components_hist)
            ]
            plot = Plot(entries,
                        what='hist',
                        xtitle='%s bin' % ('Detector' if detector_space else 'Generator'),
                        ytitle='Component of chi2 (#Delta V^{-1} #Delta)',
                        ylim=(0, y_max),
                        has_data=False,
                        )
            plot.default_canvas_size = (800, 600)
            plot.plot("NOSTACK HIST")
            l,t = debug_plotter.draw_pt_binning_lines(plot, which='reco' if detector_space else 'gen', axis='x')
            plot.save(os.path.join(debugging_dir, 'components.pdf'))

            pt_bin_edges = self.pt_bin_edges_reco if detector_space else self.pt_bin_edges_gen

            for ibin_pt, (pt_low, pt_high) in enumerate(zip(pt_bin_edges[:-1], pt_bin_edges[1:])):
                # plot component (delta * V_inv * delta) for this pt bin
                this_h = self.hist_bin_chopper.get_var_hist_pt_binned(components_hist, ibin_pt, binning_scheme='detector' if detector_space else 'generator')
                entries = [
                    Contribution(this_h, label=None)
                ]
                plot = Plot(entries,
                            what='hist',
                            title='%g < p_{T} %g GeV' % (pt_low, pt_high),
                            xtitle='lambda variable',
                            ytitle='Component of chi2 (#Delta V^{-1} #Delta)',
                            ylim=(0, y_max),
                            has_data=False,
                            )
                plot.default_canvas_size = (800, 600)
                plot.plot("NOSTACK TEXT HIST")
                plot.save(os.path.join(debugging_dir, 'components_pt_bin_%d.pdf' % (ibin_pt)))

                # plot delta for this pt bin
                this_delta = self.hist_bin_chopper.get_var_hist_pt_binned(delta_hist, ibin_pt, binning_scheme='detector' if detector_space else 'generator')
                entries = [
                    Contribution(this_delta, label=None)
                ]
                plot = Plot(entries,
                            what='hist',
                            title='%g < p_{T} %g GeV' % (pt_low, pt_high),
                            xtitle='lambda variable',
                            ytitle='#Delta',
                            # ylim=(0, y_max),
                            has_data=False,
                            )
                plot.default_canvas_size = (800, 600)
                plot.plot("NOSTACK HIST TEXT")
                plot.save(os.path.join(debugging_dir, 'delta_pt_bin_%d.pdf' % (ibin_pt)))

                # plot cov matrix for this bin
                binning = self.detector_binning.FindNode("detectordistribution") if detector_space else self.generator_binning.FindNode("generatordistribution")
                var_bins = np.array(binning.GetDistributionBinning(0))
                pt_bins = np.array(binning.GetDistributionBinning(1))
                # -1 since ndarray is 0 index, th2 are 1-indexed
                start = binning.GetGlobalBinNumber(var_bins[0] * 1.001, pt_bins[ibin_pt]*1.001) - 1
                end = binning.GetGlobalBinNumber(var_bins[-2] * 1.001, pt_bins[ibin_pt]*1.001) - 1
                this_cov = cu.ndarray_to_th2(v_inv[start:end+1,start:end+1])
                canv = ROOT.TCanvas(cu.get_unique_str(), "Inverse covariance matrix V^{-1}", 800, 600)
                this_cov.SetTitle("Inverse covariance matrix V^{-1} for %g < p_{T} < %g GeV" % (pt_low, pt_high))
                this_cov.Draw("COLZ TEXT")
                canv.SetLeftMargin(0.15)
                canv.SetRightMargin(0.18)
                # canv.SetLogz(1)
                # this_cov.SetMinimum(cov_min)
                # this_cov.SetMaximum(cov_max)
                canv.SaveAs(os.path.join(debugging_dir, 'cov_inv_matrix_pt_bin_%d.pdf' % (ibin_pt)))

            # histogram the non-zero components
            h = ROOT.TH1D("h_components", "", 25 if detector_space else 10, -5, 5)
            for i, x in enumerate(components[0]):
                if x > 0:
                    h.Fill(math.sqrt(x) * np.sign(delta[0][i]))
            h.Fit('gaus')
            entries = [
                Contribution(h)
            ]
            plot = Plot(entries,
                        what='hist',
                        xtitle='sign(#Delta) #times #sqrt{#Delta V^{-1} #Delta}',
                        ytitle='N',
                        has_data=False,
                        )
            plot.default_canvas_size = (600, 600)
            plot.plot("NOSTACK")
            # have to plot TH1 not THStack to get fit box, this is ugly hack
            canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
            entries[0].obj.Draw()
            canv.SaveAs(os.path.join(debugging_dir, 'components_pull.pdf'))

        return chi2, ndof, p

    # METHODS FOR JACKKNIFED UNCERTAINTIES
    # --------------------------------------------------------------------------
    def update_stat_response_from_jackknife(self, jackknife_variations):
        """Use jackknife results to update absolute response matrix stat uncert"""
        num_vars = len(jackknife_variations)
        # Factor to scale up to full uncert
        scale_factor = len(jackknife_variations) / (len(jackknife_variations) - 1)
        # Go through each bin of the output distribution,
        # get the RMS of the variations, apply scale factor
        # This then becomes the response uncertainty in that bin
        unfolded_rsp_err = self.get_output().Clone("unfolded_rsp_err")
        all_values = []
        for ix in range(1, self.get_output().GetNbinsX()+1):
            values = [jk['unfolder'].get_output().GetBinContent(ix)
                      for jk in jackknife_variations]
            all_values.append(values)
            scaled_rms = np.std(values, ddof=0) * scale_factor
            unfolded_rsp_err.SetBinError(ix, scaled_rms)

        # Calculate covariance matrix using all "observations"
        # Be careful about ddof - MUST specify it since default is diff to std()
        all_values = np.array(all_values, dtype='float64')
        cov_matrix = np.cov(all_values, rowvar=True, ddof=0)
        cov_matrix *= (scale_factor**2)

        self.unfolded_rsp_err = unfolded_rsp_err

        # Setup response uncert error matrix
        if getattr(self, 'ematrix_stat_response', None) is None:
            # Create a new error matrix if one doesn't exist
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            self.ematrix_stat_response = this_binning.CreateErrorMatrixHistogram("ematrix_stat_response_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")

        # NB this includes overflow bins
        for ix in range(len(cov_matrix)):
            for iy in range(len(cov_matrix)):
                self.ematrix_stat_response.SetBinContent(ix+1, iy+1, cov_matrix[ix][iy])
                self.ematrix_stat_response.SetBinError(ix+1, iy+1, 0)

        # update HistBinChopper objects & cache
        if self.rsp_uncert_name in self.hist_bin_chopper.objects:
            self.hist_bin_chopper.objects[self.rsp_uncert_name] = unfolded_rsp_err
            for k, v in self.hist_bin_chopper._cache.items():
                if self.rsp_uncert_name in k:
                    del self.hist_bin_chopper._cache[k]  # reset HBC cache


    def update_input_stat_uncert_from_jackknife(self, jackknife_variations):
        """Use jackknife results to update absolute input stat uncert

        TODO what about uncert from backgrounds?
        """
        num_vars = len(jackknife_variations)
        # Factor to scale up to full uncert
        scale_factor = len(jackknife_variations) / (len(jackknife_variations) - 1)
        # Go through each bin of the output distribution,
        # get the RMS of the variations, apply scale factor
        # This then becomes the response uncertainty in that bin
        unfolded_stat_err = self.get_output().Clone("unfolded_stat_err")
        all_values = []
        for ix in range(1, self.get_output().GetNbinsX()+1):
            values = [jk['unfolder'].get_output().GetBinContent(ix)
                      for jk in jackknife_variations]
            all_values.append(values)
            scaled_rms = np.std(values, ddof=0) * scale_factor
            unfolded_stat_err.SetBinError(ix, scaled_rms)

        # Calculate covariance matrix using all "observations"
        # Be careful about ddof - MUST specify it since default is diff to std()
        all_values = np.array(all_values, dtype='float64')
        cov_matrix = np.cov(all_values, rowvar=True, ddof=0)
        cov_matrix *= (scale_factor**2)

        self.unfolded_stat_err = unfolded_stat_err

        # Setup input uncert error matrix
        # FIXME this is technically only input stat err, not including background stat error
        # Need to add in the latter somehow
        if getattr(self, 'ematrix_stat', None) is None:
            # Create a new error matrix if one doesn't exist
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            self.ematrix_stat = this_binning.CreateErrorMatrixHistogram("ematrix_stat_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")

        # NB this includes overflow bins
        for ix in range(len(cov_matrix)):
            for iy in range(len(cov_matrix)):
                self.ematrix_stat.SetBinContent(ix+1, iy+1, cov_matrix[ix][iy])
                self.ematrix_stat.SetBinError(ix+1, iy+1, 0)

        # update HistBinChopper objects & cache
        if self.stat_uncert_name in self.hist_bin_chopper.objects:
            self.hist_bin_chopper.objects[self.stat_uncert_name] = unfolded_stat_err
            for k, v in self.hist_bin_chopper._cache.items():
                if self.stat_uncert_name in k:
                    del self.hist_bin_chopper._cache[k]  # reset HBC cache

    # METHODS FOR NORMALISED RESULTS
    # --------------------------------------------------------------------------
    def create_normalised_jackknife_response_uncertainty_per_pt_bin(self, jackknife_variations):
        """Create response uncertainty & error matrices for each gen pt bin for use later.

        Done like the absolute jackknife uncertainty, but only on the normalised
        plots for each pt bin.
        """
        # Add just incase user hasn't done so already
        # We wont use this object - we'll overwrite the cache ourselves
        self.hist_bin_chopper.add_obj(self.rsp_uncert_name, self.get_unfolded_with_ematrix_response())
        self.hist_bin_chopper.add_obj(self.rsp_ematrix_name, self.get_ematrix_stat())

        num_vars = len(jackknife_variations)
        # Factor to scale up to full uncert
        scale_factor = len(jackknife_variations) / (len(jackknife_variations) - 1)
        for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
            # Go through each bin of the output distribution,
            # get the RMS of the variations, apply scale factor
            # This then becomes the response uncertainty in that bin
            # TODO: how to create ematrix - usual x.x^T?
            hbc_args = dict(ind=ibin_pt, binning_scheme='generator')
            bin_variations = []
            for jvar in jackknife_variations:
                self.hist_bin_chopper.add_obj(jvar['label'], jvar['unfolder'].get_output())
                bin_variations.append(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(jvar['label'], **hbc_args))

            this_bin_unfolded_rsp_err = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.rsp_uncert_name, **hbc_args).Clone()
            all_values = []
            for ix in range(1, bin_variations[0].GetNbinsX()+1):
                values = [h.GetBinContent(ix) for h in bin_variations]
                all_values.append(values)
                scaled_rms = np.std(values, ddof=0) * scale_factor
                this_bin_unfolded_rsp_err.SetBinError(ix, scaled_rms)

            # Calculate covariance matrix using all "observations"
            # Be careful about ddof - MUST specify it since default is diff to std()
            all_values = np.array(all_values, dtype='float64')
            cov_matrix = np.cov(all_values, rowvar=True, ddof=0)
            cov_matrix *= (scale_factor**2)

            bins = self.variable_bin_edges_gen
            nbins = len(bins) - 1
            # FIXME which binning to use? index or physical?
            this_bin_unfolded_rsp_ematrix = ROOT.TH2D("ematrix_rsp_bin_%d_%s" % (ibin_pt, cu.get_unique_str()), "", nbins, 0, nbins, nbins,0, nbins)
            for ix in range(nbins):
                for iy in range(nbins):
                    this_bin_unfolded_rsp_ematrix.SetBinContent(ix+1, iy+1, cov_matrix[ix][iy])
                    this_bin_unfolded_rsp_ematrix.SetBinError(ix+1, iy+1, 0)

            if ibin_pt == 0:
                # Check dimensions
                if len(cov_matrix) != this_bin_unfolded_rsp_err.GetNbinsX():
                    raise ValueError("len(cov_matrix) != this_bin_unfolded_rsp_err.GetNbinsX()")
                # Check values
                if not cu.same_floats(this_bin_unfolded_rsp_ematrix.GetBinContent(2, 2), this_bin_unfolded_rsp_err.GetBinError(2)**2):
                    raise ValueError("Mismatch in this_bin_unfolded_rsp_ematrix, this_bin_unfolded_rsp_err")

            # Store in the HistBinChopper
            key = self.hist_bin_chopper._generate_key(self.rsp_uncert_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = this_bin_unfolded_rsp_err

            key = self.hist_bin_chopper._generate_key(self.rsp_ematrix_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = this_bin_unfolded_rsp_ematrix

            # Cleanup all the jackknife variations
            keys_to_del = []
            for k in self.hist_bin_chopper._cache.keys():
                for jvar in jackknife_variations:
                    if jvar['label'] in k:
                        keys_to_del.append(k)
            for k in keys_to_del:
                self.hist_bin_chopper._cache[k]

    def create_normalised_jackknife_input_uncertainty_per_pt_bin(self, jackknife_variations):
        """Create input uncertainty & error matrices for each gen pt bin for use later.

        Done like the absolute jackknife uncertainty, but only on the normalised
        plots for each pt bin.
        """
        # Add just incase user hasn't done so already
        # We wont use this object - we'll overwrite the cache ourselves
        self.hist_bin_chopper.add_obj(self.stat_uncert_name, self.get_unfolded_with_ematrix_stat())
        self.hist_bin_chopper.add_obj(self.stat_ematrix_name, self.get_ematrix_stat())

        num_vars = len(jackknife_variations)
        # Factor to scale up to full uncert
        scale_factor = len(jackknife_variations) / (len(jackknife_variations) - 1)
        for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
            # Go through each bin of the output distribution,
            # get the RMS of the variations, apply scale factor
            # This then becomes the response uncertainty in that bin
            # TODO: how to create ematrix - usual x.x^T?
            hbc_args = dict(ind=ibin_pt, binning_scheme='generator')
            bin_variations = []
            for jvar in jackknife_variations:
                self.hist_bin_chopper.add_obj(jvar['label'], jvar['unfolder'].get_output())
                bin_variations.append(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(jvar['label'], **hbc_args))

            this_bin_unfolded_stat_err = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.stat_uncert_name, **hbc_args).Clone()
            all_values = []
            for ix in range(1, bin_variations[0].GetNbinsX()+1):
                values = [h.GetBinContent(ix) for h in bin_variations]
                all_values.append(values)
                scaled_rms = np.std(values, ddof=0) * scale_factor
                this_bin_unfolded_stat_err.SetBinError(ix, scaled_rms)

            # Calculate covariance matrix using all "observations"
            # Be careful about ddof - MUST specify it since default is diff to std()
            all_values = np.array(all_values, dtype='float64')
            cov_matrix = np.cov(all_values, rowvar=True, ddof=0)
            cov_matrix *= (scale_factor**2)

            bins = self.variable_bin_edges_gen
            nbins = len(bins) - 1
            # FIXME which binning to use? index or physical?
            this_bin_unfolded_stat_ematrix = ROOT.TH2D("ematrix_stat_bin_%d_%s" % (ibin_pt, cu.get_unique_str()), "", nbins, 0, nbins, nbins,0, nbins)
            for ix in range(nbins):
                for iy in range(nbins):
                    this_bin_unfolded_stat_ematrix.SetBinContent(ix+1, iy+1, cov_matrix[ix][iy])
                    this_bin_unfolded_stat_ematrix.SetBinError(ix+1, iy+1, 0)

            if ibin_pt == 0:
                # Check dimensions
                if len(cov_matrix) != this_bin_unfolded_stat_err.GetNbinsX():
                    raise ValueError("len(cov_matrix) != this_bin_unfolded_stat_err.GetNbinsX()")
                # Check values
                if not cu.same_floats(this_bin_unfolded_stat_ematrix.GetBinContent(2, 2), this_bin_unfolded_stat_err.GetBinError(2)**2):
                    raise ValueError("Mismatch in this_bin_unfolded_stat_ematrix, this_bin_unfolded_stat_err")

            # Store in the HistBinChopper
            key = self.hist_bin_chopper._generate_key(self.stat_uncert_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = this_bin_unfolded_stat_err

            key = self.hist_bin_chopper._generate_key(self.stat_ematrix_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = this_bin_unfolded_stat_ematrix

            # Cleanup all the jackknife variations
            keys_to_del = []
            for k in self.hist_bin_chopper._cache.keys():
                for jvar in jackknife_variations:
                    if jvar['label'] in k:
                        keys_to_del.append(k)
            for k in keys_to_del:
                self.hist_bin_chopper._cache[k]


    def create_normalised_scale_syst_uncertainty_per_pt_bin(self, scale_systs):
        """Create scale uncertainty from unfolding with scale variation response matrices.
        Stores hist where error bar is envelope of variations of unfolded result.

        This is done by taking in all the unfolded scale systematics results,
        then getting the normalised result for each pT bin.
        We can then figure out the envelope of the scale variations.
        We can then store this as an extra uncertianty, to be added in quadrature later.

        scale_systs is a list of dicts, the ones produced in unfolding.py
        Each has the form:
        {
            "label": "muR up, muF nominal",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNom', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure,
            "unfolder": MyUnfolder,
        }
        """
        for syst in scale_systs:
            syst['hbc_key_unfolded'] = 'scale_syst_%s_unfolded' % cu.no_space_str(syst['label'])
            self.hist_bin_chopper.add_obj(syst['hbc_key_unfolded'], syst['unfolder'].get_unfolded_with_ematrix_stat())

        # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
        self.hist_bin_chopper.add_obj(self.scale_uncert_name, self.get_unfolded_with_ematrix_stat())

        self.hist_bin_chopper.add_obj('unfolded_stat_err', self.get_unfolded_with_ematrix_stat())

        # print("Doing scale variation")
        for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
            hbc_args = dict(ind=ibin_pt, binning_scheme='generator')
            variations = [
                self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(syst['hbc_key_unfolded'], **hbc_args)
                for syst in scale_systs
            ]

            # Calculate envelope error bar from max variation in each bin
            nominal = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            variations_envelope = nominal.Clone("scale_envelope_pt_bin%d" % ibin_pt)

            # print("pt bin", ibin_pt)
            for ix in range(1, variations_envelope.GetNbinsX()+1):
                max_variation = max([abs(v.GetBinContent(ix) - nominal.GetBinContent(ix))
                                     for v in variations])
                variations_envelope.SetBinError(ix, max_variation)

            # Store in hist_bin_chopper for later
            key = self.hist_bin_chopper._generate_key(self.scale_uncert_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = variations_envelope


    def create_normalised_scale_syst_ematrices_per_pt_bin(self):
        """Create ematrix corresponding to scale uncertainty for each pt bin

        Required create_normalised_scale_syst_uncertainty_per_pt_bin() to be run first

        Calculated as x.x^T, where x is the difference between the scale-uncert
        normalised hist, and the nominal normalised hist

        Stored in self.hist_bin_chopper with key self.scale_uncert_ematrix_name
        """
        # add dummy object to HistBinChopper.objects so check doesn't fail
        self.hist_bin_chopper.add_obj(self.scale_uncert_ematrix_name, self.get_unfolded_with_ematrix_stat())

        for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
            nominal = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err',
                                                                             ibin_pt,
                                                                             binning_scheme='generator')
            scale_hist = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.scale_uncert_name,
                                                                               ibin_pt,
                                                                               binning_scheme='generator')
            scale_shift = self.convert_error_bars_to_error_shift(scale_hist)
            scale_ematrix = cu.shift_to_covariance(scale_shift)
            key = self.hist_bin_chopper._generate_key(self.scale_uncert_ematrix_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = scale_ematrix


    def create_normalised_pdf_syst_uncertainty_per_pt_bin(self, pdf_systs):
        """Create PDF uncertainty from unfolded PDF variations

        This is done by taking in all the unfolded pdf systematics results,
        then getting the normalised result for each pT bin.
        Then we figure out the RMS of these variations.
        We can then store this as an extra uncertianty, to be added in quadrature later.

        pdf_systs is a list of dicts, the ones produced in unfolding.py
        Each has the form:
        {
            "label": "PDF",  # this is a template entry, used for future
            "tfile": os.path.join(source_dir_systs, 'PDFvariationsTrue', qgc.QCD_FILENAME),
            "colour": ROOT.kCyan+2,
            "unfolder": None,
        }
        """
        # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
        self.hist_bin_chopper.add_obj(self.pdf_uncert_name, self.get_unfolded_with_ematrix_stat())

        for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
            # Calculate error by using RMS of variations in each bin of the histogram
            variations = [syst['unfolder'].hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin_pt, binning_scheme='generator')
                          for syst in pdf_systs]

            variations_envelope = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin_pt, binning_scheme='generator').Clone("pdf_%d" % (ibin_pt))

            for ix in range(1, variations_envelope.GetNbinsX()+1):
                # np.std does sqrt((abs(x - x.mean())**2) / (len(x) - ddof)),
                # and the PDF4LHC recommendation is N-1 in the denominator
                rms_ratio = np.std([v.GetBinContent(ix) for v in variations], ddof=1)
                variations_envelope.SetBinError(ix, rms_ratio)

            # Store in hist_bin_chopper for later
            key = self.hist_bin_chopper._generate_key(self.pdf_uncert_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = variations_envelope

    def create_normalised_pdf_syst_ematrices_per_pt_bin(self):
        """Create ematrix corresponding to pdf uncertainty for each pt bin

        Requires create_normalised_pdf_syst_uncertainty_per_pt_bin() to be run first

        Calculated as x.x^T, where x is the difference between the pdf-uncert
        normalised hist, and the nominal normalised hist

        Stored in self.hist_bin_chopper with key self.pdf_uncert_ematrix_name
        """
        # add dummy object to HistBinChopper.objects so check doesn't fail
        self.hist_bin_chopper.add_obj(self.pdf_uncert_ematrix_name, self.get_unfolded_with_ematrix_stat())

        for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
            nominal = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err',
                                                                             ibin_pt,
                                                                             binning_scheme='generator')
            pdf_hist = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.pdf_uncert_name,
                                                                               ibin_pt,
                                                                               binning_scheme='generator')
            pdf_shift = self.convert_error_bars_to_error_shift(pdf_hist)
            pdf_ematrix = cu.shift_to_covariance(pdf_shift)
            key = self.hist_bin_chopper._generate_key(self.pdf_uncert_ematrix_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = pdf_ematrix


    @staticmethod
    def remove_error_bars(h):
        for i in range(1, h.GetNbinsX()+1):
            h.SetBinError(i, 0)

    @staticmethod
    def convert_error_bars_to_error_shift(h):
        """Create histogram with bin contents equal to error bar on h,
        and 0 error bars"""
        h_new = h.Clone(cu.get_unique_str())
        for i in range(1, h.GetNbinsX()+1):
            h_new.SetBinContent(i, h.GetBinError(i))
            h_new.SetBinError(i, 0)
        return h_new

    @staticmethod
    def convert_error_shift_to_error_bars(h_unshifted, h_shifted):
        """Create histogram with bin contents from h_unshifted,
        and error bars from bin values of h_shifted"""
        h = h_unshifted.Clone(cu.get_unique_str())
        for i in range(1, h_unshifted.GetNbinsX()+1):
            h.SetBinError(i, h_shifted.GetBinContent(i))
        return h

    def setup_normalised_experimental_systs_per_pt_bin(self):
        """Setup normalised experimental uncertainties

        In particular recalculates uncertainties, since the systematic one are non-trivial.
        We must re-calculate the systematic shifts on the _normalised_ result,
        by comparing the normalised nominal and shifted hists for each pt bin,
        then add those in quadrature per bin of each histogram.

        From the difference in normalised hists, we can then recalculate the
        error (covariance) matrix, like as is done in TUnfold
        """
        self.hist_bin_chopper.add_obj('unfolded_stat_err', self.get_unfolded_with_ematrix_stat())

        for exp_syst in self.exp_systs:
            if exp_syst.syst_shifted is None:
                raise ValueError("Need to run _post_process() or get_syst_shifted_hist() first to fill this")
            self.hist_bin_chopper.add_obj(exp_syst.syst_shifted_label, exp_syst.syst_shifted)

            # add these dummy obj to to HistBinChopper for later, but it isn't used
            # Just to bypass internal checks that it exists in its cached objects
            # when e.g. get_pt_bin_normed_div_bin_width() called
            self.hist_bin_chopper.add_obj(exp_syst.syst_shift_label, exp_syst.syst_shift)
            self.hist_bin_chopper.add_obj(exp_syst.syst_ematrix_label, exp_syst.syst_shifted)

            # Now manually recalculate the syst shifts and store them
            # And from this calculate the error matrix for this pt bin and store
            for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
                syst_hist = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(exp_syst.syst_shifted_label, ibin_pt, binning_scheme='generator').Clone()
                nominal_hist = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin_pt, binning_scheme='generator')
                syst_hist.Add(nominal_hist, -1)
                self.remove_error_bars(syst_hist)
                key = self.hist_bin_chopper._generate_key(exp_syst.syst_shift_label,
                                                          ind=ibin_pt,
                                                          axis='pt',
                                                          do_norm=True,
                                                          do_div_bin_width=True,
                                                          binning_scheme='generator')
                self.hist_bin_chopper._cache[key] = syst_hist

                syst_ematrix = cu.shift_to_covariance(syst_hist)
                key = self.hist_bin_chopper._generate_key(exp_syst.syst_ematrix_label,
                                                          ind=ibin_pt,
                                                          axis='pt',
                                                          do_norm=True,
                                                          do_div_bin_width=True,
                                                          binning_scheme='generator')
                self.hist_bin_chopper._cache[key] = syst_ematrix


    @staticmethod
    def get_sub_th2(h2d, start_bin, end_bin):
        """Create square TH2D from sub-matrix of h2d, from start_bin to end_bin (inclusive)"""
        nbins = end_bin - start_bin + 1
        h2d_new = ROOT.TH2D("h2d_"+cu.get_unique_str(), "", nbins, 0, nbins, nbins, 0, nbins)
        for ix_new, ix in enumerate(range(start_bin, end_bin+1), 1):
            for iy_new, iy in enumerate(range(start_bin, end_bin+1), 1):
                value = h2d.GetBinContent(ix, iy)
                err = h2d.GetBinError(ix, iy)
                h2d_new.SetBinContent(ix_new, iy_new, value)
                h2d_new.SetBinError(ix_new, iy_new, err)
        return h2d_new

    @staticmethod
    def scale_th2_bin_widths(h2d, bins):
        """Scale bins of a square TH2 by bin widths

        bins is a list of bin edges, must have 1 more value than the number of bins in h2d
        """
        if len(bins) != h2d.GetNbinsX()+1:
            print(bins)
            print(h2d.GetNbinsX())
            raise ValueError("Wrong number of bins to scale x axis")
        if len(bins) != h2d.GetNbinsY()+1:
            raise ValueError("Wrong number of bins to scale y axis")
        for ix, (binx_low, binx_high) in enumerate(zip(bins[:-1], bins[1:]), 1):
            for iy, (biny_low, biny_high) in enumerate(zip(bins[:-1], bins[1:]), 1):
                width_x = binx_high - binx_low
                width_y = biny_high - biny_low
                scale = width_x * width_y
                value = h2d.GetBinContent(ix, iy)
                err = h2d.GetBinError(ix, iy)
                h2d.SetBinContent(ix, iy, value / scale)
                h2d.SetBinError(ix, iy, err / scale)

    def setup_normalised_results_per_pt_bin(self):
        """Setup final normalised results per pt bin with all uncertainties.

        Experimental, model & PDF normalised systs should have already been setup.
        """
        self.hist_bin_chopper.add_obj('hist_truth', self.hist_truth)
        self.hist_bin_chopper.add_obj('unfolded', self.get_output())
        self.hist_bin_chopper.add_obj('unfolded_stat_err', self.get_unfolded_with_ematrix_stat())
        self.hist_bin_chopper.add_obj('unfolded_rsp_err', self.get_unfolded_with_ematrix_response())

        # add dummy objects to fool check
        self.hist_bin_chopper.add_obj(self.stat_ematrix_name, self.get_ematrix_stat())
        self.hist_bin_chopper.add_obj(self.rsp_ematrix_name, self.get_ematrix_stat())
        self.hist_bin_chopper.add_obj(self.total_ematrix_name, self.get_ematrix_stat())

        # For each pt bin, recalculate total error in quadrature and store in unfolded hist
        for ibin_pt, pt in enumerate(self.pt_bin_edges_gen[:-1]):
            first_bin = ibin_pt == 0
            hbc_args = dict(ind=ibin_pt, binning_scheme='generator')

            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_rsp_err', **hbc_args)

            unfolded_hist_bin_abs = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', **hbc_args)
            norm = unfolded_hist_bin_abs.Integral("width") / unfolded_hist_bin_stat_errors.Integral("width")

            # create stat & rsp err covariance matrices for this pt bin,
            # if they haven't been calculated by jackknife methods,
            # scaling by overall normalisation and bin widths
            binning = self.generator_binning.FindNode("generatordistribution")
            var_bins = self.variable_bin_edges_gen
            # FIXME what to do if non-sequential bin numbers?!
            # the 1.0001 is to ensure we're def inside this bin
            start_bin = binning.GetGlobalBinNumber(var_bins[0]*1.0001, pt*1.0001)
            end_bin = binning.GetGlobalBinNumber(var_bins[-2]*1.0001, pt*1.0001)  # -2 since the last one is the upper edge of the last bin
            stat_key = self.hist_bin_chopper._generate_key(self.stat_ematrix_name,
                                                           ind=ibin_pt,
                                                           axis='pt',
                                                           do_norm=True,
                                                           do_div_bin_width=True,
                                                           binning_scheme='generator')
            if stat_key not in self.hist_bin_chopper._cache:
                # Get the stat error matrix from TUnfold, then select the sub-matrix
                # for this pt bin, then scale by normalisation and bin widths
                stat_ematrix = self.get_sub_th2(self.get_ematrix_stat(), start_bin, end_bin)
                stat_ematrix.Scale(1./(norm*norm))
                self.scale_th2_bin_widths(stat_ematrix, var_bins)
                self.hist_bin_chopper._cache[stat_key] = stat_ematrix

            rsp_key = self.hist_bin_chopper._generate_key(self.rsp_ematrix_name,
                                                          ind=ibin_pt,
                                                          axis='pt',
                                                          do_norm=True,
                                                          do_div_bin_width=True,
                                                          binning_scheme='generator')
            if rsp_key not in self.hist_bin_chopper._cache:
                # if it has been setup already, it was from jackknife
                # otherwise we use the one from TUnfold
                rsp_ematrix = self.get_sub_th2(self.get_ematrix_stat_response(), start_bin, end_bin)
                rsp_ematrix.Scale(1./(norm*norm))
                self.scale_th2_bin_widths(rsp_ematrix, var_bins)
                self.hist_bin_chopper._cache[rsp_key] = rsp_ematrix

            # Calculate total ematrix
            total_ematrix = self.hist_bin_chopper._cache[stat_key].Clone()
            total_ematrix.Add(self.hist_bin_chopper._cache[rsp_key])

            for exp_syst in self.exp_systs:
                if first_bin:
                    print("Adding", exp_syst.label, "ematrix to total normalised ematrix...")
                total_ematrix.Add(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(exp_syst.syst_ematrix_label, **hbc_args))

            if self.scale_uncert_ematrix_name in self.hist_bin_chopper.objects:
                if first_bin:
                    print("Adding scale ematrix to total normalised ematrix")
                total_ematrix.Add(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.scale_uncert_ematrix_name, **hbc_args))

            if self.pdf_uncert_ematrix_name in self.hist_bin_chopper.objects:
                if first_bin:
                    print("Adding pdf ematrix to total normalised ematrix")
                total_ematrix.Add(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.pdf_uncert_ematrix_name, **hbc_args))

            key = self.hist_bin_chopper._generate_key(self.total_ematrix_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = total_ematrix

            error_bar_hists = [unfolded_hist_bin_stat_errors, unfolded_hist_bin_rsp_errors]

            # convert all shifts to error bars
            for exp_syst in self.exp_systs:
                if first_bin:
                    print("Adding", exp_syst.label, "uncertainty to normalised result...")
                # Here we access the things we just manually put in the cache - must match up with key!
                # Don't worry about it being normed etc - that is just so keys agree, and it matches
                # up with the nominal result (which we do want normed_div_bin_width)
                syst_shift = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(exp_syst.syst_shift_label, **hbc_args)
                error_bar_hists.append(self.convert_error_shift_to_error_bars(unfolded_hist_bin_stat_errors, syst_shift))

            # Add in scale syst
            if self.scale_uncert_name in self.hist_bin_chopper.objects:
                if first_bin:
                    print("Adding scale uncertainty to normalised result...")
                error_bar_hists.append(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.scale_uncert_name, **hbc_args))

            # Add in PDF syst
            if self.pdf_uncert_name in self.hist_bin_chopper.objects:
                if first_bin:
                    print("Adding PDF uncertainty to normalised result...")
                error_bar_hists.append(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.pdf_uncert_name, **hbc_args))

            # Get normalised hist with nominal unfolded value, and change error bars
            # to be quadrature sum of those we want (stat+rsp+systs)
            h_total = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)
            for i in range(1, h_total.GetNbinsX()+1):
                err2 = sum([pow(h.GetBinError(i), 2) for h in error_bar_hists])
                h_total.SetBinError(i, math.sqrt(err2))
            # if first_bin:
                # print("total ematrix diags:", [h_total.GetBinError(i) for i in range(1, nbins+1)])

            # Sanity check
            if not cu.same_floats(h_total.GetBinError(3)**2, total_ematrix.GetBinContent(3, 3)):
                print("h_total:", h_total.GetBinError(3)**2)
                print("total_ematrix:", total_ematrix.GetBinContent(3, 3))
                raise ValueError("Disagreement between h_total and total_ematrix: you screwed it up somewhere")

            # Update cache
            key = self.hist_bin_chopper._generate_key('unfolded',
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=True,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = h_total

    # METHODS FOR ABSOLUTE PER LAMBDA BIN RESULTS
    # --------------------------------------------------------------------------
    def create_scale_syst_uncertainty_per_lambda_bin(self, scale_systs):
        """Create scale uncertainty from unfolding with scale variation response matrices.
        Stores hist where error bar is envelope of variations of unfolded result.

        Note that because it's absolute values, we don't have to worry about normalising

        scale_systs is a list of dicts, the ones produced in unfolding.py
        Each has the form:
        {
            "label": "muR up, muF nominal",
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNom', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure,
            "unfolder": MyUnfolder,
        }
        """
        for syst in scale_systs:
            syst['hbc_key_unfolded'] = 'scale_syst_%s_unfolded' % cu.no_space_str(syst['label'])
            self.hist_bin_chopper.add_obj(syst['hbc_key_unfolded'], syst['unfolder'].get_unfolded_with_ematrix_stat())

        # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
        self.hist_bin_chopper.add_obj(self.scale_uncert_name, self.get_unfolded_with_ematrix_stat())

        self.hist_bin_chopper.add_obj('unfolded_stat_err', self.get_unfolded_with_ematrix_stat())

        # print("Doing scale variation")
        for ibin_var in range(len(self.variable_bin_edges_gen[:-1])):
            hbc_args = dict(ind=ibin_var, binning_scheme='generator')
            variations = [
                self.hist_bin_chopper.get_lambda_bin_div_bin_width(syst['hbc_key_unfolded'], **hbc_args)
                for syst in scale_systs
            ]

            # Calculate envelope error bar from max variation in each bin
            nominal = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', **hbc_args)
            variations_envelope = nominal.Clone("scale_envelope_lambda_bin%d" % ibin_var)

            # print("pt bin", ibin_var)
            for ix in range(1, variations_envelope.GetNbinsX()+1):
                max_variation = max([abs(v.GetBinContent(ix) - nominal.GetBinContent(ix))
                                     for v in variations])
                variations_envelope.SetBinError(ix, max_variation)

            # Store in hist_bin_chopper for later
            key = self.hist_bin_chopper._generate_key(self.scale_uncert_name,
                                                      ind=ibin_var,
                                                      axis='lambda',
                                                      do_norm=False,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = variations_envelope

    def create_scale_syst_ematrices_per_lambda_bin(self):
        pass

    def create_pdf_syst_uncertainty_per_lambda_bin(self, pdf_systs):
        """Create PDF uncertainty from unfolded PDF variations

        This is done by taking in all the unfolded pdf systematics results,
        then we figure out the RMS of these variations.
        We can then store this as an extra uncertianty, to be added in quadrature later.

        pdf_systs is a list of dicts, the ones produced in unfolding.py
        Each has the form:
        {
            "label": "PDF",  # this is a template entry, used for future
            "tfile": os.path.join(source_dir_systs, 'PDFvariationsTrue', qgc.QCD_FILENAME),
            "colour": ROOT.kCyan+2,
            "unfolder": None,
        }
        """
        # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
        self.hist_bin_chopper.add_obj(self.pdf_uncert_name, self.get_unfolded_with_ematrix_stat())

        for ibin_var in range(len(self.variable_bin_edges_gen[:-1])):
            hbc_args = dict(ind=ibin_var, binning_scheme="generator")
            # Calculate error by using RMS of variations in each bin of the histogram
            variations = [syst['unfolder'].hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', **hbc_args)
                          for syst in pdf_systs]

            variations_envelope = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', **hbc_args).Clone("pdf_%d" % (ibin_var))

            for ix in range(1, variations_envelope.GetNbinsX()+1):
                # np.std does sqrt((abs(x - x.mean())**2) / (len(x) - ddof)),
                # and the PDF4LHC recommendation is N-1 in the denominator
                rms_ratio = np.std([v.GetBinContent(ix) for v in variations], ddof=1)
                variations_envelope.SetBinError(ix, rms_ratio)

            # Store in hist_bin_chopper for later
            key = self.hist_bin_chopper._generate_key(self.pdf_uncert_name,
                                                      ind=ibin_var,
                                                      axis='lambda',
                                                      do_norm=False,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = variations_envelope

    def create_pdf_syst_ematrices_per_lambda_bin(self):
        pass

    def setup_absolute_results_per_lambda_bin(self):
        """Setup final absolute results per pt bin with all uncertainties.

        Experimental, model & PDF absolute systs should have already been setup.
        """
        self.hist_bin_chopper.add_obj('hist_truth', self.hist_truth)
        self.hist_bin_chopper.add_obj('unfolded', self.get_output())
        self.hist_bin_chopper.add_obj('unfolded_stat_err', self.get_unfolded_with_ematrix_stat())
        self.hist_bin_chopper.add_obj('unfolded_rsp_err', self.get_unfolded_with_ematrix_response())

        # add dummy objects to fool check
        self.hist_bin_chopper.add_obj(self.stat_ematrix_name, self.get_ematrix_stat())
        self.hist_bin_chopper.add_obj(self.rsp_ematrix_name, self.get_ematrix_stat())
        self.hist_bin_chopper.add_obj(self.total_ematrix_name, self.get_ematrix_stat())

        # For each lambda bin, recalculate total error in quadrature and store in unfolded hist
        for ibin_var, var in enumerate(self.variable_bin_edges_gen[:-1]):
            first_bin = ibin_var == 0
            hbc_args = dict(ind=ibin_var, binning_scheme='generator')

            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', **hbc_args)
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_rsp_err', **hbc_args)

            # # create stat & rsp err covariance matrices for this pt bin,
            # # if they haven't been calculated by jackknife methods,
            # # scaling by overall normalisation and bin widths
            # binning = self.generator_binning.FindNode("generatordistribution")
            # pt_bins = self.pt_bin_edges_gen
            # # FIXME what to do if non-sequential bin numbers?!
            # # the 1.0001 is to ensure we're def inside this bin
            # start_bin = binning.GetGlobalBinNumber(var*1.0001, pt_bins[0]*1.0001)
            # end_bin = binning.GetGlobalBinNumber(var*1.0001, pt_bins[-2]*1.0001)  # -2 since the last one is the upper edge of the last bin
            # stat_key = self.hist_bin_chopper._generate_key(self.stat_ematrix_name,
            #                                                ind=ibin_var,
            #                                                axis='lambda',
            #                                                do_norm=False,
            #                                                do_div_bin_width=True,
            #                                                binning_scheme='generator')
            # if stat_key not in self.hist_bin_chopper._cache:
            #     # Get the stat error matrix from TUnfold, then select the sub-matrix
            #     # for this pt bin, then scale by normalisation and bin widths
            #     stat_ematrix = self.get_sub_th2(self.get_ematrix_stat(), start_bin, end_bin)
            #     self.scale_th2_bin_widths(stat_ematrix, var_bins)
            #     self.hist_bin_chopper._cache[stat_key] = stat_ematrix

            # rsp_key = self.hist_bin_chopper._generate_key(self.rsp_ematrix_name,
            #                                               ind=ibin_var,
            #                                               axis='lambda',
            #                                               do_norm=False,
            #                                               do_div_bin_width=True,
            #                                               binning_scheme='generator')
            # if rsp_key not in self.hist_bin_chopper._cache:
            #     # if it has been setup already, it was from jackknife
            #     # otherwise we use the one from TUnfold
            #     rsp_ematrix = self.get_sub_th2(self.get_ematrix_stat_response(), start_bin, end_bin)
            #     self.scale_th2_bin_widths(rsp_ematrix, var_bins)
            #     self.hist_bin_chopper._cache[rsp_key] = rsp_ematrix

            # # Calculate total ematrix
            # total_ematrix = self.hist_bin_chopper._cache[stat_key].Clone()
            # total_ematrix.Add(self.hist_bin_chopper._cache[rsp_key])

            # for exp_syst in self.exp_systs:
            #     if first_bin:
            #         print("Adding", exp_syst.label, "ematrix to total absolute ematrix...")
            #     total_ematrix.Add(self.hist_bin_chopper.get_lambda_bin_div_bin_width(exp_syst.syst_ematrix_label, **hbc_args))

            # if self.scale_uncert_ematrix_name in self.hist_bin_chopper.objects:
            #     if first_bin:
            #         print("Adding scale ematrix to total absolute ematrix")
            #     total_ematrix.Add(self.hist_bin_chopper.get_lambda_bin_div_bin_width(self.scale_uncert_ematrix_name, **hbc_args))

            # if self.pdf_uncert_ematrix_name in self.hist_bin_chopper.objects:
            #     if first_bin:
            #         print("Adding pdf ematrix to total absolute ematrix")
            #     total_ematrix.Add(self.hist_bin_chopper.get_lambda_bin_div_bin_width(self.pdf_uncert_ematrix_name, **hbc_args))

            # key = self.hist_bin_chopper._generate_key(self.total_ematrix_name,
            #                                           ind=ibin_var,
            #                                           axis='lambda',
            #                                           do_norm=False,
            #                                           do_div_bin_width=True,
            #                                           binning_scheme='generator')
            # self.hist_bin_chopper._cache[key] = total_ematrix

            error_bar_hists = [unfolded_hist_bin_stat_errors, unfolded_hist_bin_rsp_errors]

            # convert all shifts to error bars
            for exp_syst in self.exp_systs:
                if first_bin:
                    print("Adding", exp_syst.label, "uncertainty to absolute result...")
                # Here we access the things we just manually put in the cache - must match up with key!
                # Don't worry about it being normed etc - that is just so keys agree, and it matches
                # up with the nominal result (which we do want normed_div_bin_width)
                # Create shift due to this syst
                syst_shift = self.hist_bin_chopper.get_lambda_bin_div_bin_width(exp_syst.syst_shifted_label, **hbc_args).Clone()
                syst_shift.Add(unfolded_hist_bin_stat_errors, -1)
                error_bar_hists.append(self.convert_error_shift_to_error_bars(unfolded_hist_bin_stat_errors, syst_shift))

            # Add in scale syst
            if self.scale_uncert_name in self.hist_bin_chopper.objects:
                if first_bin:
                    print("Adding scale uncertainty to absolute result...")
                error_bar_hists.append(self.hist_bin_chopper.get_lambda_bin_div_bin_width(self.scale_uncert_name, **hbc_args))

            # Add in PDF syst
            if self.pdf_uncert_name in self.hist_bin_chopper.objects:
                if first_bin:
                    print("Adding PDF uncertainty to absolute result...")
                error_bar_hists.append(self.hist_bin_chopper.get_lambda_bin_div_bin_width(self.pdf_uncert_name, **hbc_args))

            # Get absolute hist with nominal unfolded value, and change error bars
            # to be quadrature sum of those we want (stat+rsp+systs)
            h_total = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', **hbc_args)
            for i in range(1, h_total.GetNbinsX()+1):
                err2 = sum([pow(h.GetBinError(i), 2) for h in error_bar_hists])
                h_total.SetBinError(i, math.sqrt(err2))
            # if first_bin:
                # print("total ematrix diags:", [h_total.GetBinError(i) for i in range(1, nbins+1)])

            # Sanity check
            # if not cu.same_floats(h_total.GetBinError(3)**2, total_ematrix.GetBinContent(3, 3)):
            #     print("h_total:", h_total.GetBinError(3)**2)
            #     print("total_ematrix:", total_ematrix.GetBinContent(3, 3))
            #     raise ValueError("Disagreement between h_total and total_ematrix: you screwed it up somewhere")

            # Update cache
            key = self.hist_bin_chopper._generate_key('unfolded',
                                                      ind=ibin_var,
                                                      axis='lambda',
                                                      do_norm=False,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = h_total


def pickle_region(region, output_filename, infos=True, convert_tfile_to_str=True):
    """pickle a region dict

    You should use this + unpickle_region() to ensure correct compression algo used

    Parameters
    ----------
    region : dict
        region dict, with unfolder + other infos, systematics, etc
    output_filename : str
        pickle filename
    infos : bool, optional
        Print out sizes of components
    convert_tfile_to_str : bool, optional
        Convert TFile to their filepath names before pickling
    """
    if convert_tfile_to_str:
        # recursively change TFile objects back to filenames
        def _convert_tfile_to_str(d):
            for k in d.keys():
                # TODO: write pickler class?
                if isinstance(d[k], ROOT.TFile):
                    filename = d[k].GetName()
                    if infos:
                        print(" - closing", k, filename)
                    d[k].Close()
                    d[k] = filename
                elif isinstance(d[k], dict):
                    _convert_tfile_to_str(d[k])
                elif isinstance(d[k], list):
                    for x in d[k]:
                        _convert_tfile_to_str(x)
        _convert_tfile_to_str(region)
    if infos:
        print("")
        print("region sizes:")
        print("-"*80)
        cu.print_dict_item_sizes(region, recursive=True)
        print("-"*80)
    # LZMA for huge space saving, but very slow when unpickling lots of big dicts
    # i.e. with PDF + exp + model systs
    # with lzma.open(output_filename, "wb") as f:
    # gzip for speed and some space saving
    with gzip.open(output_filename, "wb") as f:
        pickle.dump(region, f, protocol=2) # protocol 2 means very compatible across python versions


def unpickle_region(pickle_filename):
    """Retreive region dict from pickle file"""
    if not os.path.isfile(pickle_filename):
        print("! Warning ! cannot find unfolding pickle file", pickle_filename, ' - skipping')
        return None
    # with lzma.open(pickle_filename, 'r') as f:
    with gzip.open(pickle_filename, 'r') as f:
        unpickled_region = pickle.load(f)
    return unpickled_region


def unpack_slim_unfolding_root_file(input_tfile, region_name, angle_name, pt_bins):
    tdir = "%s/%s" % (region_name, angle_name)
    unfolding_stat_err_hists = []
    unfolding_stat_err_ematrices = []
    unfolding_total_err_hists = []
    unfolding_total_err_ematrices = []
    truth_hists = []
    alt_truth_hists = []
    for ibin_pt, _ in enumerate(pt_bins[:-1]):
        unfolding_stat_err_hist = input_tfile.Get("%s/unfolded_stat_err_norm_divBinWidth_%d" % (tdir, ibin_pt))
        unfolding_stat_err_hist.SetDirectory(0)
        unfolding_stat_err_hists.append(unfolding_stat_err_hist)

        unfolding_stat_err_ematrix = input_tfile.Get("%s/unfolded_stat_ematrix_norm_divBinWidth_%d" % (tdir, ibin_pt))
        unfolding_stat_err_ematrix.SetDirectory(0)
        unfolding_stat_err_ematrices.append(unfolding_stat_err_ematrix)

        unfolding_total_err_hist = input_tfile.Get("%s/unfolded_norm_divBinWidth_%d" % (tdir, ibin_pt))
        unfolding_total_err_hist.SetDirectory(0)
        unfolding_total_err_hists.append(unfolding_total_err_hist)

        unfolding_total_err_ematrix = input_tfile.Get("%s/unfolded_total_ematrix_norm_divBinWidth_%d" % (tdir, ibin_pt))
        unfolding_total_err_ematrix.SetDirectory(0)
        unfolding_total_err_ematrices.append(unfolding_total_err_ematrix)

        truth_hist = input_tfile.Get("%s/hist_truth_norm_divBinWidth_%d" % (tdir, ibin_pt))
        truth_hist.SetDirectory(0)
        truth_hists.append(truth_hist)

        alt_truth_hist = input_tfile.Get("%s/alt_hist_truth_norm_divBinWidth_%d" % (tdir, ibin_pt))
        alt_truth_hist.SetDirectory(0)
        alt_truth_hists.append(alt_truth_hist)

    return dict(
        unfolding_stat_err_hists=unfolding_stat_err_hists,
        unfolding_stat_ematrics=unfolding_stat_err_ematrices,
        unfolding_total_err_hists=unfolding_total_err_hists,
        unfolding_total_ematrices=unfolding_total_err_ematrices,
        truth_hists=truth_hists,
        alt_truth_hists=alt_truth_hists,
    )


# def unfolder_from_tdir(tdir):
#     """Recover MyUnfolder from tdirectory

#     Massive pain, prefer pickling instead
#     """
#     warnings.warn("Favour pickling/unpickling instead of ROOT unpacking", DeprecationWarning)

#     unfolder = MyUnfolder(response_map=cu.get_from_tfile(tdir, "response_map"),
#                           variable_bin_edges_reco=np.array(cu.get_from_tfile(tdir, "variable_bin_edges_reco")),
#                           variable_bin_edges_gen=np.array(cu.get_from_tfile(tdir, "variable_bin_edges_gen")),
#                           variable_name=str(cu.get_from_tfile(tdir, "variable_name").GetTitle()),
#                           pt_bin_edges_reco=np.array(cu.get_from_tfile(tdir, "pt_bin_edges_reco")),
#                           pt_bin_edges_gen=np.array(cu.get_from_tfile(tdir, "pt_bin_edges_gen")),
#                           pt_bin_edges_underflow_reco=np.array(cu.get_from_tfile(tdir, "pt_bin_edges_underflow_reco")),
#                           pt_bin_edges_underflow_gen=np.array(cu.get_from_tfile(tdir, "pt_bin_edges_underflow_gen")),
#                           orientation=int(cu.get_from_tfile(tdir, "orientation")[0]),
#                           constraintMode=int(cu.get_from_tfile(tdir, "constraintMode")[0]),
#                           regMode=int(cu.get_from_tfile(tdir, "regMode")[0]),
#                           densityFlags=int(cu.get_from_tfile(tdir, "densityFlags")[0]),
#                           distribution=str(cu.get_from_tfile(tdir, "distribution").GetTitle()),
#                           axisSteering=str(cu.get_from_tfile(tdir, "axisSteering").GetTitle()))

#     obj_names = cu.get_list_of_element_names(tdir)
#     for name in obj_names:
#         obj = tdir.Get(name)
#         if name.startswith("background_reco_binning_"):
#             bg_name = cu.str_restore_space(name.replace("background_reco_binning_", ""))
#             unfolder.backgrounds[bg_name] = obj

#         elif name.startswith("background_gen_binning_"):
#             bg_name = cu.str_restore_space(name.replace("background_gen_binning_", ""))
#             unfolder.backgrounds_gen_binning[bg_name] = obj

#         elif name.startswith("syst_map_"):
#             syst_name = cu.str_restore_space(name.replace("syst_map_", ""))
#             unfolder.syst_maps[syst_name] = obj

#         elif name.startswith("syst_ematrix_"):
#             syst_name = cu.str_restore_space(name.replace("syst_ematrix_", ""))
#             unfolder.syst_ematrices[syst_name] = obj

#         elif name.startswith("syst_shift_"):
#             syst_name = cu.str_restore_space(name.replace("syst_shift_", ""))
#             unfolder.syst_shifts[syst_name] = obj

#         elif name.startswith("syst_shifted_unfolded_"):
#             syst_name = cu.str_restore_space(name.replace("syst_shifted_unfolded_", ""))
#             unfolder.systs_shifted[syst_name] = obj

#     for attr_name in MyUnfolder._simple_attr:
#         obj = tdir.Get(attr_name)
#         setattr(unfolder, attr_name, obj)

#     unfolder.setup_normalised_results_per_pt_bin()

#     return unfolder


# def unpack_unfolding_root_file(input_tfile, region, angle, do_alt_response=True, do_model_systs=True, do_pdf_systs=True):
#     """Unpack Unfolders, systematics, etc from ROOT file

#     Prefer unpickling instead!
#     """
#     warnings.warn("Favour pickling/unpickling instead of ROOT unpacking", DeprecationWarning)

#     input_tdir_name = "%s/%s" % (region['name'], angle.var)
#     input_tdir = input_tfile.Get(input_tdir_name)
#     cu.check_root_obj(input_tdir)
#     unfolder = unfolder_from_tdir(input_tdir)
#     print("...Loaded main unfolder")

#     list_of_obj = cu.get_list_of_element_names(input_tdir)

#     # Get unregularised unfolder, if available
#     unreg_tdir = [x for x in list_of_obj if x.startswith("unreg_unfolder")]
#     unreg_unfolder = None
#     if len(unreg_tdir) == 1:
#         unreg_unfolder = unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, unreg_tdir[0])))
#         print("...Loaded comparison unregularised unfolder")

#     # Update if experimental systs
#     region['experimental_systematics'] = [k for k in region['experimental_systematics']
#                                           if k['label'] in unfolder.systs_shifted]

#     # Get alternate response object, if it exists
#     alt_unfolder = None
#     alt_hist_truth = None
#     alt_hist_reco = None
#     alt_hist_reco_bg_subtracted = None
#     alt_hist_reco_bg_subtracted_gen_binning = None
#     if do_alt_response:
#         alt_tdir_names = [x for x in list_of_obj if x.startswith("alt_response_")]
#         if len(alt_tdir_names)  == 1:
#             alt_tdir_name = alt_tdir_names[0]
#             alt_hist_truth = input_tfile.Get(os.path.join(input_tdir_name, alt_tdir_name, "alt_hist_mc_gen"))
#             alt_hist_reco = input_tfile.Get(os.path.join(input_tdir_name, alt_tdir_name, "alt_hist_mc_reco"))
#             alt_hist_reco_bg_subtracted = input_tfile.Get(os.path.join(input_tdir_name, alt_tdir_name, "alt_hist_mc_reco_bg_subtracted"))
#             alt_hist_reco_bg_subtracted_gen_binning = input_tfile.Get(os.path.join(input_tdir_name, alt_tdir_name, "alt_hist_mc_reco_bg_subtracted_gen_binning"))

#             # Need to check actually unfolder stored, and not just the parts above
#             alt_tdir = input_tfile.Get(os.path.join(input_tdir_name, alt_tdir_name))
#             alt_unf_obj = [x for x in cu.get_list_of_element_names(alt_tdir)]
#             if 'response_map' in alt_unf_obj:
#                 alt_unfolder = unfolder_from_tdir(alt_tdir)
#                 region['alt_unfolder'] = alt_unfolder
#                 alt_unfolder_name = alt_tdir_name.replace("alt_response_", "")
#                 if cu.no_space_str(region['alt_mc_label']) != alt_unfolder_name:
#                     raise RuntimeError("Bad unpacking of alt response unfolder: expected %s, got %s" % (region['alt_mc_label'], alt_unfolder_name))
#                 print("...Loaded alt unfolder", alt_tdir_name)

#         if len(alt_tdir_names) > 1:
#             raise RuntimeError(">1 alt_response?! %s" % (alt_tdir_names))

#     # Get model systs
#     # print(list_of_obj)
#     if do_model_systs:
#         model_tdirs = [x for x in list_of_obj if x.startswith("modelSyst_")]
#         if len(model_tdirs) > 0:
#             for model_tdir_name in model_tdirs:
#                 syst_name = model_tdir_name.replace("modelSyst_", "")
#                 this_one = [x for x in region['model_systematics'] if cu.no_space_str(x['label']) == syst_name]
#                 if len(this_one) == 0:
#                     print("No entry for model systematic", syst_name, "- skipping")
#                     continue
#                 # TODO: check it agrees with region dict?
#                 this_one[0]['unfolder'] = unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, model_tdir_name)))
#             print("...Loaded", len(model_tdirs), "model systematic unfolders")
#     # remove entries without an unfolder
#     region['model_systematics'] = [k for k in region['model_systematics']
#                                    if k.get('unfolder', None) is not None]
#     # setup normalised errors
#     if len(region['model_systematics']) > 0:
#         print("Doing normalised scale uncertainties")
#         unfolder.create_normalised_scale_syst_uncertainty_per_pt_bin(region['model_systematics'])

#     # Get PDF systs
#     # For some reason, this is done as a list instead of dict
#     if do_pdf_systs:
#         pdf_tdirs = [x for x in list_of_obj if x.startswith("pdfSyst_")]
#         if len(pdf_tdirs) > 0:
#             # Remove original, construct all other
#             region['pdf_systematics'] = []

#             for pdf_tdir_name in pdf_tdirs:
#                 pdf_name = pdf_tdir_name.replace("pdfSyst_", "")
#                 region['pdf_systematics'].append({
#                     'label': pdf_name,
#                     'unfolder': unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, pdf_tdir_name))),
#                     'colour': ROOT.kCyan+2,
#                 })
#             print("...Loaded", len(pdf_tdirs), "PDF systematic unfolders")
#     # remove entries without an unfolder
#     region['pdf_systematics'] = [k for k in region['pdf_systematics']
#                                  if k.get('unfolder', None) is not None]
#     # setup normalised errors
#     if len(region['pdf_systematics']) > 0:
#         print("Doing normalised PDF uncertainties")
#         unfolder.create_normalised_pdf_syst_uncertainty_per_pt_bin(region['pdf_systematics'])

#     unfolder.setup_normalised_results_per_pt_bin()

#     return dict(
#         unfolder=unfolder,
#         unreg_unfolder=unreg_unfolder,
#         alt_unfolder=alt_unfolder,
#         alt_hist_truth=alt_hist_truth,
#         alt_hist_reco=alt_hist_reco,
#         alt_hist_reco_bg_subtracted=alt_hist_reco_bg_subtracted,
#         alt_hist_reco_bg_subtracted_gen_binning=alt_hist_reco_bg_subtracted_gen_binning,
#     )


class HistBinChopper(object):
    """Get histogram for pt or variable bin, and cache it in dict, so can be used later"""

    def __init__(self, generator_binning, detector_binning):
        self.generator_binning = generator_binning
        if self.generator_binning is not None:
            # do this once as expensive
            self.generator_binning_var_bins = np.array(self.generator_binning.GetDistributionBinning(0))
            self.generator_binning_pt_bins = np.array(self.generator_binning.GetDistributionBinning(1))

        self.detector_binning = detector_binning
        if self.detector_binning is not None:
            self.detector_binning_var_bins = np.array(self.detector_binning.GetDistributionBinning(0))
            self.detector_binning_pt_bins = np.array(self.detector_binning.GetDistributionBinning(1))

        self.objects = {}
        self._cache = {}
        self._cache_integral = {}

    def get_binning(self, binning_scheme):
        """Get TUnfoldBinning, lambda var bins, pt bins for binning_scheme = 'generator' or 'detector'"""
        if binning_scheme not in ['generator', 'detector']:
            raise ArgumentError('binning_scheme must be "generator" or "detector"')
        thing = self.generator_binning if binning_scheme == "generator" else self.detector_binning
        if thing is None:
            raise RuntimeError("No valid TUnfoldBinning object for binning scheme '%s'" % binning_scheme)
        var_bins = self.generator_binning_var_bins if binning_scheme == 'generator' else self.detector_binning_var_bins
        pt_bins = self.generator_binning_pt_bins if binning_scheme == 'generator' else self.detector_binning_pt_bins
        return thing, var_bins, pt_bins

    def add_obj(self, name, obj):
        # TODO: allow overwrite?
        if name not in self.objects and obj is not None:
            self.objects[name] = obj

    def update(self, other):
        """Update this object's cached things with those from another HistBinChopper"""
        self.objects.update(other.objects)
        self._cache.update(other._cache)
        self._cache_integral.update(other._cache_integral)

    def get_var_hist_pt_binned(self, hist1d, ibin_pt, binning_scheme='generator'):
        """Get hist of variable for given pt bin from massive 1D hist that TUnfold makes"""
        # FIXME: assume no underflow?!
        binning, var_bins, pt_bins = self.get_binning(binning_scheme)
        h = ROOT.TH1D("h_%d_%s" % (ibin_pt, cu.get_unique_str()), "", len(var_bins)-1, var_bins)
        for var_ind, var_value in enumerate(var_bins[:-1], 1):
            this_val = var_value * 1.001  # ensure its inside
            bin_num = binning.GetGlobalBinNumber(this_val, pt_bins[ibin_pt]*1.001)
            h.SetBinContent(var_ind, hist1d.GetBinContent(bin_num))
            h.SetBinError(var_ind, hist1d.GetBinError(bin_num))
        return h

    def get_var_2d_hist_pt_binned(self, hist2d, ibin_pt, binning_scheme='generator'):
        """Get 2d hist for given pt bin from massive 2D hist"""
        # FIXME: assume no underflow?!
        binning, var_bins, pt_bins = self.get_binning(binning_scheme)
        h = ROOT.TH2D("h2d_%d_%s" % (ibin_pt, cu.get_unique_str()), "", len(var_bins)-1, var_bins, len(var_bins)-1, var_bins)
        for var_ind, var_value in enumerate(var_bins[:-1], 1):
            this_val = var_value * 1.001  # ensure its inside
            bin_num = binning.GetGlobalBinNumber(this_val, pt_bins[ibin_pt]*1.001)
            for var_ind2, var_value2 in enumerate(var_bins[:-1], 1):
                this_val2 = var_value2 * 1.001  # ensure its inside
                bin_num2 = binning.GetGlobalBinNumber(this_val2, pt_bins[ibin_pt]*1.001)
                h.SetBinContent(var_ind, var_ind2, hist2d.GetBinContent(bin_num, bin_num2))
                h.SetBinError(var_ind, var_ind2, hist2d.GetBinError(bin_num, bin_num2))
        return h

    def get_pt_hist_var_binned(self, hist1d, ibin_var, binning_scheme='generator'):
        """Get hist of pt for given variable bin from massive 1D hist that TUnfold makes"""
        # FIXME: assume no underflow?!
        binning, var_bins, pt_bins = self.get_binning(binning_scheme)
        h = ROOT.TH1D("h_%d_%s" % (ibin_var, cu.get_unique_str()), "", len(pt_bins)-1, pt_bins)
        for pt_ind, pt_value in enumerate(pt_bins[:-1], 1):
            this_val = pt_value * 1.001  # ensure its inside
            bin_num = binning.GetGlobalBinNumber(var_bins[ibin_var]*1.001, this_val)
            h.SetBinContent(pt_ind, hist1d.GetBinContent(bin_num))
            h.SetBinError(pt_ind, hist1d.GetBinError(bin_num))
        return h

    def get_pt_2d_hist_var_binned(self, hist2d, ibin_var, binning_scheme='generator'):
        """Get 2d hist for given variable bin from massive 2D hist"""
        # FIXME: assume no underflow?!
        binning, var_bins, pt_bins = self.get_binning(binning_scheme)
        h = ROOT.TH1D("h2d_%d_%s" % (ibin_var, cu.get_unique_str()), "", len(pt_bins)-1, pt_bins)
        for pt_ind, pt_value in enumerate(pt_bins[:-1], 1):
            this_val = pt_value * 1.001  # ensure its inside
            bin_num = binning.GetGlobalBinNumber(var_bins[ibin_var]*1.001, this_val)
            for pt_ind2, pt_value2 in enumerate(pt_bins[:-1], 1):
                this_val2 = pt_value * 1.001  # ensure its inside
                bin_num2 = binning.GetGlobalBinNumber(var_bins[ibin_var]*1.001, this_val2)
                h.SetBinContent(pt_ind, pt_ind2, hist2d.GetBinContent(bin_num, bin_num2))
                h.SetBinError(pt_ind, pt_ind2, hist2d.GetBinError(bin_num, bin_num2))
        return h

    def get_bin_plot(self, name, ind, axis, do_norm=False, do_div_bin_width=False, binning_scheme='generator'):
        """Get plot for given bin (index=ind) of specified axis.

        Note, only considers signal region

        Parameters
        ----------
        name : str
            Name of object to use
        ind : int
            Bin index (0-indexed, 0 = 1st signal region bin)
        axis : str
            'pt' or 'lambda', i.e. axis to get bin ind of
        do_norm : bool, optional
            Normalise to unity
        do_div_bin_width : bool, optional
            Divide by bin width
        binning_scheme : str, optional
            'generator' or 'detector'

        Returns
        -------
        TYPE
            Description

        Raises
        ------
        KeyError
            Description
        """
        if name not in self.objects:
            raise KeyError("No '%s' in HistBinChopper.objects, only: %s" % (name, list(self.objects.keys())))
        if self.objects[name] is None:
            raise RuntimeError("HistBinChopper.objects[%s] is None" % name)

        key = self._generate_key(name, ind, axis, do_norm, do_div_bin_width, binning_scheme)
        if key not in self._cache:
            if axis == 'lambda':
                self._cache[key] = self.get_pt_hist_var_binned(self.objects[name], ind, binning_scheme)
            else:
                self._cache[key] = self.get_var_hist_pt_binned(self.objects[name], ind, binning_scheme)

            # havent done div bin width or normalising yet
            self._cache_integral[key] = self._cache[key].Integral()

            if do_div_bin_width and do_norm:
                self._cache[key] = qgp.normalise_hist_divide_bin_width(self._cache[key])
            elif do_div_bin_width and not do_norm:
                self._cache[key] = qgp.hist_divide_bin_width(self._cache[key])
            elif not do_div_bin_width and do_norm:
                qgp.normalise_hist(self._cache[key])

        return self._cache[key]

    def get_bin_integral(self, name, ind, axis, binning_scheme='generator'):
        """Get integral for bin `ind` of object `name`"""
        if name not in self.objects:
            raise KeyError("No '%s' in HistBinChopper.objects, only: %s" % (name, list(self.objects.keys())))
        key = self._generate_key(name, ind, axis, False, False, binning_scheme)
        if key not in self._cache_integral:
            self.get_bin_plot(name, ind, axis, False, False, binning_scheme)
        return self._cache_integral[key]

    @staticmethod
    def _generate_key(name, ind, axis, do_norm, do_div_bin_width, binning_scheme):
        """Generate consistent name for these args, options as in get_bin_plot()"""
        if axis not in ['pt', 'lambda']:
            raise ArgumentError('_generate_key(): axis must be "pt" or "lambda"')
        key = name + "_%s_bin_%d_%s" % (axis, ind, binning_scheme)
        if do_norm:
            key += "_norm"
        if do_div_bin_width:
            key += "_divBinWidth"
        return key

    # TODO: remove these? just use get_bin_plot instead?
    def get_pt_bin(self, name, ind, binning_scheme='generator'):
        return self.get_bin_plot(name, ind, axis='pt', do_norm=False, do_div_bin_width=False, binning_scheme=binning_scheme)

    def get_pt_bin_div_bin_width(self, name, ind, binning_scheme='generator'):
        return self.get_bin_plot(name, ind, axis='pt', do_norm=False, do_div_bin_width=True, binning_scheme=binning_scheme)

    def get_pt_bin_normed_div_bin_width(self, name, ind, binning_scheme='generator'):
        return self.get_bin_plot(name, ind, axis='pt', do_norm=True, do_div_bin_width=True, binning_scheme=binning_scheme)

    def get_lambda_bin(self, name, ind, binning_scheme='generator'):
        return self.get_bin_plot(name, ind, axis='lambda', do_norm=False, do_div_bin_width=False, binning_scheme=binning_scheme)

    def get_lambda_bin_div_bin_width(self, name, ind, binning_scheme='generator'):
        return self.get_bin_plot(name, ind, axis='lambda', do_norm=False, do_div_bin_width=True, binning_scheme=binning_scheme)

    def get_lambda_bin_normed_div_bin_width(self, name, ind, binning_scheme='generator'):
        return self.get_bin_plot(name, ind, axis='lambda', do_norm=True, do_div_bin_width=True, binning_scheme=binning_scheme)


class ExpSystematic(object):
    """Class to hold info about an experimental systematic"""

    def __init__(self, label, syst_map=None, syst_shift=None, syst_shifted=None, syst_ematrix=None):
        self.label = label
        self.label_no_spaces = cu.no_space_str(label)
        self.syst_map = syst_map

        # shift to norminal unfolded result(from TUnfold)
        self.syst_shift = syst_shift
        # these labels are for HistBinChopper
        self.syst_shift_label = 'syst_shift_%s' % (self.label_no_spaces)

        # norminal + shift
        self.syst_shifted = syst_shifted
        self.syst_shifted_label = 'syst_shifted_%s' % (self.label_no_spaces)

        # error matrix (= shift * shift^T)
        self.syst_ematrix = syst_ematrix
        self.syst_ematrix_label = 'syst_ematrix_%s' % (self.label_no_spaces)

