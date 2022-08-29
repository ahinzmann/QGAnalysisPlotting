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
import pickle
import gzip
import warnings
from functools import partial
from collections import namedtuple
from bisect import bisect_right, bisect_left

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot

import common_utils as cu
import qg_general_plots as qgp
from my_unfolder_plotter import MyUnfolderPlotter


# monkey-patch warning formatter
warnings.formatwarning = cu._formatwarning


# This doesn't seem to work...sigh
np.set_printoptions(edgeitems=3, infstr='Infinity',
                    linewidth=220, nanstr='nan', precision=6,
                    suppress=False, threshold=sys.maxsize, formatter=None)

My_Style.cd()
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

# hold (pt, variable) pair
PtVar = namedtuple("PtVar", ["pt", "var"])


class PtVarBinning(object):
    """Hold a set of pt, variable binning. Assumes separate underflow pT region."""

    def __init__(self,
                 variable_bin_edges,  # 'variable' refers to e.g. ptD, LHA
                 variable_name,
                 pt_bin_edges_signal,
                 pt_bin_edges_underflow,
                 binning_name,
                 binning_underflow_name,
                 binning_signal_name,
                 var_uf=False,  # _uf = underflow
                 var_of=True,  # _of = overflow
                 pt_uf=False,
                 pt_of=True):
        self.variable_name = variable_name
        self.variable_bin_edges = np.array(variable_bin_edges, dtype='float')
        self.nbins_variable = len(variable_bin_edges)-1 if variable_bin_edges is not None else 0

        self.pt_bin_edges_signal = np.array(pt_bin_edges_signal, dtype='float')
        self.nbins_pt = len(pt_bin_edges_signal) - 1 if pt_bin_edges_signal is not None else 0

        self.pt_bin_edges_underflow = np.array(pt_bin_edges_underflow, dtype='float')
        self.nbins_pt_underflow = len(pt_bin_edges_underflow)-1 if pt_bin_edges_underflow is not None else 0

        self.var_uf = var_uf
        self.var_of = var_of
        self.pt_uf = pt_uf
        self.pt_of = pt_of

        self.pt_name = "pt"

        self.binning_name = binning_name
        self.binning = ROOT.TUnfoldBinning(self.binning_name)

        # Do pt underflow ourselves as separate region
        self.binning_underflow_name = binning_underflow_name
        self.distribution_underflow = self.binning.AddBinning(self.binning_underflow_name)
        if self.variable_bin_edges is not None:
            self.distribution_underflow.AddAxis(self.variable_name, self.nbins_variable, self.variable_bin_edges,
                                                self.var_uf, self.var_of)
        self.distribution_underflow.AddAxis(self.pt_name, self.nbins_pt_underflow, self.pt_bin_edges_underflow,
                                            self.pt_uf, False)

        # Signal pt region
        self.binning_signal_name = binning_signal_name
        self.distribution = self.binning.AddBinning(self.binning_signal_name)
        if self.variable_bin_edges is not None:
            self.distribution.AddAxis(self.variable_name, self.nbins_variable, self.variable_bin_edges,
                                      self.var_uf, self.var_of)
        self.distribution.AddAxis(self.pt_name, self.nbins_pt, self.pt_bin_edges_signal,
                                  False, self.pt_of)

        # Hold maps of global bin number < > physical bin edges
        self.global_bin_to_physical_val_map = dict()
        self.physical_val_to_global_bin_map = dict()

        self._cache_global_bin_mapping()
        # print("bin 1:", self.global_bin_to_physical_val_map[1])

    def __str__(self):
        s = """PtVarBinning(
    variable_bin_edges={variable_bin_edges},
    variable_name={variable_name},
    pt_bin_edges_signal={pt_bin_edges_signal},
    pt_bin_edges_underflow={pt_bin_edges_underflow},
    binning_name={binning_name},
    binning_underflow_name={binning_underflow_name},
    binning_signal_name={binning_signal_name},
    var_uf={var_uf},
    var_of={var_of},
    pt_uf={pt_uf},
    pt_of={pt_of})""".format(**self.__dict__)
        return s

    def is_signal_region(self, pt):
        return pt >= self.pt_bin_edges_signal[0]

    def get_distribution(self, pt):
        return self.distribution if self.is_signal_region(pt) else self.distribution_underflow

    def get_pt_bins(self, is_signal_region):
        return self.pt_bin_edges_signal if is_signal_region else self.pt_bin_edges_underflow

    def get_variable_bins(self, pt):
        # pt in args to duck type with PtVarPerPtBinning method
        return self.variable_bin_edges

    def get_start_bin(self, is_signal_region):
        if is_signal_region:
            return self.physical_bin_to_global_bin(pt=self.pt_bin_edges_signal[0], var=self.variable_bin_edges[0])
        else:
            return self.physical_bin_to_global_bin(pt=self.pt_bin_edges_underflow[0], var=self.variable_bin_edges[0])

    def _cache_global_bin_mapping(self):
        """Create maps of global bin <> physical bin values,
        by iterating through all the physical bins (inc oflow)

        Have to do this as the TUnfoldBinning one isnt good.
        """
        all_pt_bins = list(chain(self.pt_bin_edges_underflow[:-1], self.pt_bin_edges_signal))
        for ibin_pt, pt in enumerate(all_pt_bins[:-1]):
            this_pt = pt+1E-6
            pt_bin = (pt, all_pt_bins[ibin_pt+1])
            binning_obj = self.get_distribution(this_pt)
            if self.variable_bin_edges is not None:
                for ibin_var, variable in enumerate(self.variable_bin_edges[:-1]):
                    global_bin = binning_obj.GetGlobalBinNumber(variable+1E-6, this_pt)
                    self.global_bin_to_physical_val_map[global_bin] = PtVar(pt=pt_bin,
                                                                            var=(variable, self.variable_bin_edges[ibin_var+1]))
                if self.var_of:
                    variable = self.variable_bin_edges[-1]
                    global_bin = binning_obj.GetGlobalBinNumber(variable+1E-6, this_pt)
                    self.global_bin_to_physical_val_map[global_bin] = PtVar(pt=pt_bin,
                                                                            var=(variable, np.inf))
            else:
                global_bin = binning_obj.GetGlobalBinNumber(this_pt)
                self.global_bin_to_physical_val_map[global_bin] = PtVar(pt=pt_bin,
                                                                        var=None)

        if self.pt_of:
            pt = self.pt_bin_edges_signal[-1]
            this_pt = pt+1E-6
            pt_bin = (pt, np.inf)
            binning_obj = self.get_distribution(this_pt)
            if self.variable_bin_edges is not None:
                for ibin_var, variable in enumerate(self.variable_bin_edges[:-1]):
                    global_bin = binning_obj.GetGlobalBinNumber(variable+1E-6, this_pt)
                    self.global_bin_to_physical_val_map[global_bin] = PtVar(pt=pt_bin,
                                                                            var=(variable, self.variable_bin_edges[ibin_var+1]))
                if self.var_of:
                    variable = self.variable_bin_edges[-1]
                    global_bin = binning_obj.GetGlobalBinNumber(variable+1E-6, this_pt)
                    self.global_bin_to_physical_val_map[global_bin] = PtVar(pt=pt_bin,
                                                                            var=(variable, np.inf))
            else:
                global_bin = binning_obj.GetGlobalBinNumber(this_pt)
                self.global_bin_to_physical_val_map[global_bin] = PtVar(pt=pt_bin,
                                                                        var=None)

        # now invert
        self.physical_val_to_global_bin_map = {v: k for k, v in self.global_bin_to_physical_val_map.items()}

    def physical_bin_to_global_bin(self, pt, var):
        dist = self.get_distribution(pt)
        return dist.GetGlobalBinNumber(var, pt)

    def global_bin_to_physical_bin(self, global_bin_number):
        """returns PtVar obj"""
        if global_bin_number not in self.global_bin_to_physical_val_map:
            raise KeyError("No global bin %d" % global_bin_number)
        return self.global_bin_to_physical_val_map[global_bin_number]

    def get_first_pt_overflow_global_bin(self):
        """Get global bin corresponding to first pt overflow bin,
        or the last bin+1 if no pt overflow configured"""
        if not self.pt_of:
            return max(self.global_bin_to_physical_val_map.keys())+1

        pt = self.pt_bin_edges_signal[-1] + 1E-6
        binning_obj = self.get_distribution(pt)
        if self.variable_bin_edges is not None:
            global_bin = binning_obj.GetGlobalBinNumber(self.variable_bin_edges[0]+1E-6, pt)
        else:
            global_bin = binning_obj.GetGlobalBinNumber(pt)
        return global_bin

    def get_total_number_of_bins(self):
        return len(self.global_bin_to_physical_val_map)


def flatten_unique(list_of_lists):
    # uses itertool recipe
    return sorted(list(set(chain.from_iterable(list_of_lists))))


class PtVarPerPtBinning(object):
    """Hold a set of pt, variable binning. Separate variable binning per pT bin.
    Assumes separate underflow pT region.

    NB no pT overflow, since we might want to do it with a different lambda binning,
    so assume user has added it themselves

    Structure of pt_*_bin_config:

    [
        ([50, 65], [0, 0.1, 0.2, 0.5, 1.0]), # [0] is pt bins edge, [1] is the variable bins
        ([65, 88], [0, 0.1, 0.25, 0.4, 1.0]),
        ([88, 200], [0, 0.1, 0.3, 0.75, 1.0]),
    ]
    """

    def __init__(self,
                 variable_name,
                 pt_underflow_bin_config,
                 pt_signal_bin_config,
                 binning_name,
                 binning_underflow_name,
                 binning_signal_name,
                 var_uf=False,  # _uf = underflow
                 var_of=True):
        self.variable_name = variable_name
        self.pt_underflow_bin_config = pt_underflow_bin_config
        self.pt_signal_bin_config = pt_signal_bin_config

        # convert pairs of bin edges into ordered list
        self.pt_bin_edges_underflow = np.array(flatten_unique(b[0] for b in pt_underflow_bin_config))
        self.nbins_pt_underflow = len(self.pt_bin_edges_underflow)-1

        self.pt_bin_edges_signal = np.array(flatten_unique(b[0] for b in pt_signal_bin_config))
        self.nbins_pt = len(self.pt_bin_edges_signal)-1

        self.var_uf = var_uf
        self.var_of = var_of
        self.pt_uf = False
        self.pt_of = False  # since we handle it ourselves

        self.pt_name = "pt"

        self.binning_name = binning_name
        self.binning_underflow_name = binning_underflow_name
        self.binning_signal_name = binning_signal_name

        self.binning = ROOT.TUnfoldBinning(self.binning_name)

        # Do pt underflow ourselves as separate region
        # TODO merge into one?
        self.underflow_distributions = []
        for pt_ind, pt in enumerate(self.pt_bin_edges_underflow[:-1]):
            distribution_underflow = self.binning.AddBinning(self.binning_underflow_name + str(pt_ind))
            variable_bin_edges = self.pt_underflow_bin_config[pt_ind][1]
            distribution_underflow.AddAxis(self.variable_name,
                                           len(variable_bin_edges)-1,
                                           np.array(variable_bin_edges, 'd'),
                                           self.var_uf, self.var_of)
            pt_bin_edges_signal = np.array([pt, self.pt_bin_edges_underflow[pt_ind+1]], 'd')
            distribution_underflow.AddAxis(self.pt_name,
                                           1,
                                           pt_bin_edges_signal,
                                           False, False)
            # pt underflow must always be false, because we have to specify
            # the lambda binning for it explicitly, so it will be part of the pt_signal_bin_config
            self.underflow_distributions.append(distribution_underflow)

        # Signal pt region
        self.signal_distributions = []
        for pt_ind, pt in enumerate(self.pt_bin_edges_signal[:-1]):
            distribution_signal = self.binning.AddBinning(self.binning_signal_name + str(pt_ind))
            variable_bin_edges = self.pt_signal_bin_config[pt_ind][1]
            distribution_signal.AddAxis(self.variable_name,
                                        len(variable_bin_edges)-1,
                                        np.array(variable_bin_edges, 'd'),
                                        self.var_uf, self.var_of)
            pt_bin_edges_signal = np.array([pt, self.pt_bin_edges_signal[pt_ind+1]], 'd')
            distribution_signal.AddAxis(self.pt_name,
                                        1,
                                        pt_bin_edges_signal,
                                        False, False)
            # pt overflow must always be false, because we have to specify
            # the lambda binning for it explicitly, so it will be part of the pt_signal_bin_config
            self.signal_distributions.append(distribution_signal)

        # Hold maps of global bin number < > physical bin edges
        self.global_bin_to_physical_val_map = dict()
        self.physical_val_to_global_bin_map = dict()

        self._cache_global_bin_mapping()
        # print("bin 1:", self.global_bin_to_physical_val_map[1])


    def __str__(self):
        s = """PtVarPerPtBinning(
    variable_name={variable_name},
    pt_underflow_bin_config={pt_underflow_bin_config},
    pt_signal_bin_config={pt_signal_bin_config},
    binning_name={binning_name},
    binning_underflow_name={binning_underflow_name},
    binning_signal_name={binning_signal_name},
    var_uf={var_uf},
    var_of={var_of})""".format(**self.__dict__)
        return s

    def is_signal_region(self, pt):
        return pt >= self.pt_bin_edges_signal[0]

    def get_distribution(self, pt):
        pt_bins = self.pt_bin_edges_signal if self.is_signal_region(pt) else self.pt_bin_edges_underflow
        distributions = self.signal_distributions if self.is_signal_region(pt) else self.underflow_distributions
        pos = bisect_right(pt_bins, pt) - 1
        if pos > (len(distributions)-1):
            raise IndexError("Can't find distribution for pt %f in bins %s with pos %d (len = %s)" % (pt, pt_bins, pos, len(distributions)))
        return distributions[pos]

    def get_pt_bins(self, is_signal_region):
        return self.pt_bin_edges_signal if is_signal_region else self.pt_bin_edges_underflow

    def get_variable_bins(self, pt):
        pt_bins = self.pt_bin_edges_signal if self.is_signal_region(pt) else self.pt_bin_edges_underflow
        config = self.pt_signal_bin_config if self.is_signal_region(pt) else self.pt_underflow_bin_config
        pos = bisect_right(pt_bins, pt)-1
        if pos > (len(config)-1):
            raise IndexError("Can't find variable bins for pt %f" % pt)
        thing = config[pos][1]
        if not isinstance(thing, np.ndarray):
            return np.array(thing)
        else:
            return thing

    def get_start_bin(self, is_signal_region):
        pt = self.pt_bin_edges_underflow[0]
        if is_signal_region:
            pt = self.pt_bin_edges_signal[0]
        return self.physical_bin_to_global_bin(pt=pt, var=self.get_variable_bins(pt)[0])

    def _cache_global_bin_mapping(self):
        """Create maps of global bin <> physical bin values,
        by iterating through all the physical bins (inc oflow)

        Have to do this as the TUnfoldBinning one isnt good.
        """
        for ibin_pt, config in enumerate(chain(self.pt_underflow_bin_config, self.pt_signal_bin_config)):
            this_pt = config[0][0] + 1E-6
            pt_bin = tuple(config[0])
            binning_obj = self.get_distribution(this_pt)
            var_bin_edges = config[1]
            for var_low, var_high in zip(var_bin_edges[:-1], var_bin_edges[1:]):
                global_bin = binning_obj.GetGlobalBinNumber(var_low+1E-6, this_pt)
                self.global_bin_to_physical_val_map[global_bin] = PtVar(pt=pt_bin,
                                                                        var=(var_low, var_high))
            if self.var_of:
                variable = var_bin_edges[-1]
                global_bin = binning_obj.GetGlobalBinNumber(variable+1E-6, this_pt)
                self.global_bin_to_physical_val_map[global_bin] = PtVar(pt=pt_bin,
                                                                        var=(variable, np.inf))

        # now invert
        self.physical_val_to_global_bin_map = {v: k for k, v in self.global_bin_to_physical_val_map.items()}

    def physical_bin_to_global_bin(self, pt, var):
        dist = self.get_distribution(pt)
        return dist.GetGlobalBinNumber(var, pt)

    def global_bin_to_physical_bin(self, global_bin_number):
        """Returns PtVar object"""
        if global_bin_number not in self.global_bin_to_physical_val_map:
            raise KeyError("No global bin %d" % global_bin_number)
        return self.global_bin_to_physical_val_map[global_bin_number]

    def get_first_pt_overflow_global_bin(self):
        """Get global bin corresponding to first pt overflow bin"""
        return self.physical_bin_to_global_bin(pt=9999999., var=100000)

    def get_total_number_of_bins(self):
        return len(self.global_bin_to_physical_val_map)


class BinningHandler(object):
    """Class to handle both detector and generator binnings"""

    def __init__(self,
                 detector_ptvar_binning,
                 generator_ptvar_binning):

        self.detector_binning_name = "detector"
        self.detector_ptvar_binning = detector_ptvar_binning

        self.generator_binning_name = "generator"
        self.generator_ptvar_binning = generator_ptvar_binning

        self.binning_mapping = {
            self.generator_binning_name: self.generator_ptvar_binning,
            self.detector_binning_name: self.detector_ptvar_binning,
        }

    def __str__(self):
        s = """BinningHandler(
    detector_ptvar_binning={detector_ptvar_binning},
    generator_ptvar_binning={generator_ptvar_binning})""".format(**self.__dict__)
        return s

    def get_binning_scheme(self, binning_scheme):
        valid_args = list(self.binning_mapping.keys())
        if binning_scheme not in valid_args:
            raise ValueError("binning_scheme should be one of: %s, instead got %s" % (", ".join(valid_args), binning_scheme))
        return self.binning_mapping[binning_scheme]

    def global_bin_to_physical_bin(self, global_bin_number, binning_scheme):
        return self.get_binning_scheme(binning_scheme).global_bin_to_physical_bin(global_bin_number)

    def physical_bin_to_global_bin(self, pt, var, binning_scheme):
        return self.get_binning_scheme(binning_scheme).physical_bin_to_global_bin(pt=pt, var=var)

    def get_first_pt_overflow_global_bin(self, binning_scheme):
        return self.get_binning_scheme(binning_scheme).get_first_pt_overflow_global_bin()

    def get_pt_bins(self, binning_scheme, is_signal_region):
        return self.get_binning_scheme(binning_scheme).get_pt_bins(is_signal_region)

    def get_variable_bins(self, pt, binning_scheme):
        return self.get_binning_scheme(binning_scheme).get_variable_bins(pt)

    def get_total_number_of_bins(self, binning_scheme):
        return self.get_binning_scheme(binning_scheme).get_total_number_of_bins()

    @property
    def variable_name(self):
        return self.detector_ptvar_binning.variable_name

    # FIXME
    # def save_binning(self, print_xml=True, txt_filename=None):
    #     """Save binning scheme to txt and/or print XML to screen"""
    #     if txt_filename:
    #         with open(txt_filename, 'w') as of:
    #             of.write("GEN BINNING\n")
    #             of.write("--------------------\n")
    #             for name, region in [
    #                 ("generator_distribution", self.generator_distribution),
    #                 ("generator_distribution_underflow", self.generator_distribution_underflow)]:
    #                 of.write("%s: bins %d - %d\n" % (name, region.GetStartBin(), region.GetEndBin()))
    #             of.write("\nDETECTOR BINNING\n")
    #             of.write("--------------------\n")
    #             for name, region in [
    #                 ("detector_distribution", self.detector_distribution),
    #                 ("detector_distribution_underflow", self.detector_distribution_underflow)]:
    #                 of.write("%s: bins %d - %d\n" % (name, region.GetStartBin(), region.GetEndBin()))
    #
    #     if print_xml:
    #         # don't know how to create a ofstream in python :( best we can do is ROOT.cout
    #         ROOT.TUnfoldBinningXML.ExportXML(self.detector_binning, ROOT.cout, True, False)
    #         ROOT.TUnfoldBinningXML.ExportXML(self.generator_binning, ROOT.cout, False, True)


class InputHandler(object):
    """Class to handle inputs for MyUnfolder, deal with fakes, etc

    Good for if you need to rebin fakes, etc
    """
    def __init__(self,
                 input_hist,
                 hist_truth=None,
                 hist_mc_reco=None,
                 hist_mc_fakes=None,
                 hist_fake_fraction=None):

        self.input_hist = input_hist.Clone() if input_hist else None
        self.hist_truth = hist_truth.Clone() if hist_truth else None  # TODO: should this be here?
        self.hist_mc_reco = hist_mc_reco.Clone() if hist_mc_reco else None
        self.hist_mc_fakes = hist_mc_fakes.Clone() if hist_mc_fakes else None
        self.hist_fake_fraction = hist_fake_fraction.Clone() if hist_fake_fraction else None  # incase the user wants to specify noe instead
        # these get modified by MyUnfolder.subtract_background()
        self.input_hist_bg_subtracted = input_hist.Clone() if input_hist else None
        self.hist_mc_reco_bg_subtracted = hist_mc_reco.Clone() if hist_mc_reco else None

    def setup_fake_fraction(self):
        """Calculate fake fraction template using hist_mc_fakes and hist_mc_reco"""
        if not self.hist_mc_fakes:
            raise ValueError("Need hist_mc_fakes for setup_fake_fraction()")
        if not self.hist_mc_reco:
            raise ValueError("Need hist_mc_reco for setup_fake_fraction()")
        self.hist_fake_fraction = self.hist_mc_fakes.Clone("fake_fraction_"+cu.get_unique_str())
        self.hist_fake_fraction.Divide(self.hist_mc_reco)

    def calc_fake_hist(self, hist):
        """Calculate fakes background from hist using internal fake_fraction"""
        if not self.hist_fake_fraction:
            self.setup_fake_fraction()
        hist_n_bins = hist.GetNbinsX()
        fake_n_bins = self.hist_fake_fraction.GetNbinsX()
        if hist_n_bins != fake_n_bins:
            raise ValueError("Mimsmatch in number of bins in calc_fake_hist(): %d vs %d" % (hist_n_bins, fake_n_bins))
        fake_hist = hist.Clone(cu.get_unique_str())
        fake_hist.Multiply(self.hist_fake_fraction)
        return fake_hist

    def subtract_fake_hist(self, hist):
        """Calculate and subtract fake background from hist

        Returns boths the fake-subtracted hist, and the calculated fakes hist
        """
        fake_hist = self.calc_fake_hist(hist)
        new_hist = hist.Clone(cu.get_unique_str())
        new_hist.Add(fake_hist, -1)
        return new_hist, fake_hist


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
        "vyy_inv_no_bg_subtraction_th2",
        # cov matr
        "vyy_no_bg_subtraction_th2",
    ]

    def __init__(self,
                 response_map,
                 binning_handler,
                 orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                 constraintMode=ROOT.TUnfold.kEConstraintArea,
                 regMode=ROOT.TUnfold.kRegModeCurvature,
                 densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth,
                 distribution='generatordistribution',
                 axisSteering='*[b]'):
        """MyUnfolder constructor.

        Parameters
        ----------
        response_map : ROOT.TH2
            Response map
        binning_handler : BinningHandler
            BinningHandler object to contain gen & reco binnings
        orientation : int, optional
            Orientation of response map.
            ROOT.TUnfold.kHistMapOutputHoriz if gen on x axis.
        constraintMode : int, optional
            Optional area constraint term to enforce when unfolding
        regMode : int, optional
            Regularisation mode when setting up L matrix
        densityFlags : int, optional
            How to handel bin widths when setting up regularisation L matrix
        distribution : str, optional
            Name of TUnfoldBinning distribution to use when regularising
        axisSteering : str, optional
            Axis options for specifying which axes are used in regularising

        Raises
        ------
        RuntimeError
            If response_map has entries in the 0th gen bin
        """
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

        self.binning_handler = binning_handler
        self.variable_name = binning_handler.variable_name
        self.variable_name_safe = cu.no_space_str(self.variable_name)

        self.orientation = orientation
        self.constraintMode = constraintMode
        self.regMode = regMode
        self.densityFlags = densityFlags
        self.distribution = distribution
        self.axisSteering = axisSteering

        self.check_binning_consistency()

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

        # store inputs
        self.input_handler = None
        self.input_handler_gen_binning = None

        self.tau = 0  # to be set by user later, via TauScanner or LCurveScanner
        self.backgrounds = {}  # gets filled with subtract_background()
        self.backgrounds_gen_binning = {}  # gets filled with subtract_background_gen_binning()

        self.unfolded = None  # set in get_output(), total error
        self.unfolded_stat_err = None  # set in get_unfolded_with_ematrix_stat()
        self.unfolded_rsp_err = None  # set in get_unfolded_with_ematrix_response()

        self.exp_systs = []  # list of ExpSystematic objects to hold systematics
        self.exp_systs_normed = []  # list of ExpSystematic objects to hold normalised systematics (created in normalise_all_systs())

        # distribution name for getting various error, rho, etc matrices
        # use "generator" for signal + underflow region
        # "generatordistribution" only for ???
        self.output_distribution_name = self.binning_handler.get_binning_scheme('generator').binning_name

        self.folded_unfolded = None  # set in get_folded_unfolded()
        self.folded_unfolded_tunfold = None  # set in get_folded_unfolded()
        self.folded_mc_truth = None  # set in get_folded_mc_truth()

        self.L_matrix_entries = []  # set in setup_L_matrix_curvature()

        self.truth_template = None  # set by TruthTemplateMaker, for regularisation

        # For producing normalised distributions
        self.hist_bin_chopper = HistBinChopper(binning_handler=self.binning_handler)

        # For setting/getting various uncerts from HistBinChopper
        self.no_uncert_name = "unfolded_no_err"

        self.stat_uncert_name = 'unfolded_stat_err'
        self.stat_ematrix_name = "stat_ematrix"

        self.rsp_uncert_name = 'unfolded_rsp_err'
        self.rsp_ematrix_name = "rsp_ematrix"

        self.scale_uncert_name = "scale_uncert"
        self.scale_uncert_ematrix_name = "scale_uncert_ematrix"

        self.pdf_uncert_name = "pdf_uncert"
        self.pdf_uncert_ematrix_name = "pdf_uncert_ematrix"

        self.total_ematrix_name = "total_ematrix"

    def check_binning_consistency(self):
        """"Check binning scheme matches with binning in response matrix

        Copied from TUnfoldDensity code
        """
        if (getattr(self, "response_map", None) is not None
                and getattr(self, "generator_binning", None) is not None):
            if self.orientation == ROOT.TUnfold.kHistMapOutputHoriz:
                n_map_output_bins = self.response_map.GetXaxis().GetNbins()
                n_map_input_bins = self.response_map.GetYaxis().GetNbins()
            else:
                n_map_output_bins = self.response_map.GetYaxis().GetNbins()
                n_map_input_bins = self.response_map.GetXaxis().GetNbins()

            n_output_bins_t = abs(self.generator_binning.GetTH1xNumberOfBins(True))
            # I don't know why it needs both?
            n_output_bins_f = abs(self.generator_binning.GetTH1xNumberOfBins(False))
            if (n_output_bins_t != n_map_output_bins) and (n_output_bins_f != n_map_output_bins):
                raise ValueError("Output binning incompatible number of bins: "
                                 "axis %d binning scheme %d (%d)" % (n_map_output_bins, n_output_bins_t, n_output_bins_f))

            n_input_bins_t = abs(self.detector_binning.GetTH1xNumberOfBins(True))
            n_input_bins_f = abs(self.detector_binning.GetTH1xNumberOfBins(False))
            if (n_input_bins_t != n_map_input_bins) and (n_input_bins_f != n_map_input_bins):
                raise ValueError("Input binning incompatible number of bins: "
                                 "axis %d binning scheme %d (%d)" % (n_map_input_bins, n_input_bins_t, n_input_bins_f))

    # rewire these to point to BinningHandler instead
    # one day I'll remove references to these
    @property
    def variable_bin_edges_reco(self):
        return self.binning_handler.detector_ptvar_binning.variable_bin_edges

    @property
    def nbins_variable_reco(self):
        return self.binning_handler.detector_ptvar_binning.nbins_variable

    @property
    def pt_bin_edges_reco(self):
        return self.binning_handler.detector_ptvar_binning.pt_bin_edges_signal

    @property
    def nbins_pt_reco(self):
        return self.binning_handler.detector_ptvar_binning.nbins_pt

    @property
    def pt_bin_edges_underflow_reco(self):
        return self.binning_handler.detector_ptvar_binning.pt_bin_edges_underflow

    @property
    def nbins_pt_underflow_reco(self):
        return self.binning_handler.detector_ptvar_binning.nbins_pt_underflow

    @property
    def detector_binning(self):
        return self.binning_handler.detector_ptvar_binning.binning

    @property
    def detector_distribution_underflow(self):
        return self.binning_handler.detector_ptvar_binning.distribution_underflow

    @property
    def detector_distribution(self):
        return self.binning_handler.detector_ptvar_binning.distribution


    @property
    def variable_bin_edges_gen(self):
        return self.binning_handler.generator_ptvar_binning.variable_bin_edges

    @property
    def nbins_variable_gen(self):
        return self.binning_handler.generator_ptvar_binning.nbins_variable

    @property
    def pt_bin_edges_gen(self):
        return self.binning_handler.generator_ptvar_binning.pt_bin_edges_signal

    @property
    def nbins_pt_gen(self):
        return self.binning_handler.generator_ptvar_binning.nbins_pt

    @property
    def pt_bin_edges_underflow_gen(self):
        return self.binning_handler.generator_ptvar_binning.pt_bin_edges_underflow

    @property
    def nbins_pt_underflow_gen(self):
        return self.binning_handler.generator_ptvar_binning.nbins_pt_underflow

    @property
    def generator_binning(self):
        return self.binning_handler.generator_ptvar_binning.binning

    @property
    def generator_distribution_underflow(self):
        return self.binning_handler.generator_ptvar_binning.distribution_underflow

    @property
    def generator_distribution(self):
        return self.binning_handler.generator_ptvar_binning.distribution

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
    #     self._check_save_to_tfile(tfile, self.response_map, "response_map")
    #     self._check_save_to_tfile(tfile, self.orientation, "orientation")
    #     self._check_save_to_tfile(tfile, self.constraintMode, "constraintMode")
    #     self._check_save_to_tfile(tfile, self.regMode, "regMode")
    #     self._check_save_to_tfile(tfile, self.densityFlags, "densityFlags")
    #     self._check_save_to_tfile(tfile, self.distribution, "distribution")
    #     self._check_save_to_tfile(tfile, self.axisSteering, "axisSteering")

    #     self._check_save_to_tfile(tfile, self.binning_handler, "binning_handler")
    #     self._check_save_to_tfile(tfile, self.variable_name, "variable_name")

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

        def _write(name, ibin, binning, root_name):
            try:
                tfile.WriteTObject(self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(name, ibin, binning), root_name)
            except KeyError:
                print("Missing", name, ibin)
                pass

        for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
            _write(self.stat_uncert_name, ibin_pt, 'generator', "unfolded_stat_err_norm_divBinWidth_%d" % ibin_pt)
            _write("unfolded", ibin_pt, 'generator', "unfolded_norm_divBinWidth_%d" % ibin_pt)
            _write("hist_truth", ibin_pt, 'generator', "hist_truth_norm_divBinWidth_%d" % ibin_pt)
            _write("alt_hist_truth", ibin_pt, 'generator', "alt_hist_truth_norm_divBinWidth_%d" % ibin_pt)
            _write(self.stat_ematrix_name, ibin_pt, 'generator', 'unfolded_stat_ematrix_norm_divBinWidth_%d' % ibin_pt)
            _write(self.rsp_ematrix_name, ibin_pt, 'generator', 'unfolded_rsp_ematrix_norm_divBinWidth_%d' % ibin_pt)
            _write(self.total_ematrix_name, ibin_pt, 'generator', 'unfolded_total_ematrix_norm_divBinWidth_%d' % ibin_pt)

    def __getstate__(self):
        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()

        def _del_state(key):
            if key in state and state[key] is not None:
                del state[key]
        # Remove the large entries from being pickled
        # Normally they can be reconstructed from ROOT objects instead
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

    def slim_down(self,
                  keep_response_map=False,
                  keep_cov_matrices=False,
                  keep_1d_hists=False,
                  keep_backgrounds=False,
                  keep_folded_hists=False):
        """Delete various large member objects to reduce its size"""

        def _del_attr(name, owner=self):
            if hasattr(owner, name) and getattr(owner, name) is not None:
                try:
                    delattr(owner, name)
                except AttributeError as e:
                    print(cu.pcolors.RED + "Error using _del_attr(%s):" % name + cu.pcolors.ENDC)
                    raise e

        if not keep_response_map:
            for x in ['response_map',
                      'probability_matrix',
                      '_probability_ndarray']:
                _del_attr(x)

        if not keep_cov_matrices:
            for x in ['vxx_th2',
                      'vxx_inv_th2',
                      'ematrix_stat',
                      'ematrix_tau',
                      'ematrix_stat_response',
                      'ematrix_tunfold_total',
                      'ematrix_input',
                      'rhoij_total']:
                _del_attr(x)

        if not keep_1d_hists:
            for handler_name in ['input_handler', 'input_handler_gen_binning']:
                handler = getattr(self, handler_name)
                for x in ['input_hist',
                          'input_hist_bg_subtracted',
                          'hist_mc_reco',
                          'hist_mc_reco_bg_subtracted',
                          'hist_mc_fakes',
                          'hist_fake_fraction']:
                    _del_attr(x, owner=handler)

        if not keep_backgrounds:
            self.backgrounds.clear()
            self.backgrounds_gen_binning.clear()

        if not keep_folded_hists:
            for x in ['folded_unfolded_tunfold',
                      'folded_unfolded',
                      'folded_mc_truth']:
                _del_attr(x)

    # SETUP INPUT, BACKGROUNDS
    # --------------------------------------------------------------------------
    def set_input(self,
                  input_handler,
                  input_handler_gen_binning=None,
                  bias_factor=0,
                  error_unconstrained_truth_bins=True):
        """Set hist to be unfolded, plus other basic hists

        Also allow other args to be passed to TUnfoldSys::SetInput
        """
        self.input_handler = input_handler
        self.input_handler_gen_binning = input_handler_gen_binning

        # Check dimensions are consistent between 1D hists and response_map
        n_reco_bins = self.response_map.GetNbinsY() if self.orientation == ROOT.TUnfold.kHistMapOutputHoriz else self.response_map.GetNbinsX()
        if self.input_handler.input_hist.GetNbinsX() != n_reco_bins:
            raise ValueError("Mistmatch in input_handler.input_hist bins: %d vs %d" % (n_reco_bins, self.input_handler.input_hist.GetNbinsX()))
        if self.input_handler.hist_mc_reco and self.input_handler.hist_mc_reco.GetNbinsX() != n_reco_bins:
            raise ValueError("Mistmatch in input_handler.hist_mc_reco bins: %d vs %d" % (n_reco_bins, self.input_handler.hist_mc_reco.GetNbinsX()))
        if self.input_handler.hist_mc_fakes and self.input_handler.hist_mc_fakes.GetNbinsX() != n_reco_bins:
            raise ValueError("Mistmatch in input_handler.hist_mc_fakes bins: %d vs %d" % (n_reco_bins, self.input_handler.hist_mc_fakes.GetNbinsX()))

        n_gen_bins = self.response_map.GetNbinsX() if self.orientation == ROOT.TUnfold.kHistMapOutputHoriz else self.response_map.GetNbinsY()
        if self.input_handler.hist_truth.GetNbinsX() != n_gen_bins:
            raise ValueError("Mistmatch in input_handler.hist_truth bins: %d vs %d" % (n_gen_bins, self.input_handler.hist_truth.GetNbinsX()))
        if self.input_handler_gen_binning.input_hist.GetNbinsX() != n_gen_bins:
            raise ValueError("Mistmatch in input_handler_gen_binning.input_hist bins: %d vs %d" % (n_gen_bins, self.input_handler_gen_binning.input_hist.GetNbinsX()))
        if self.input_handler_gen_binning.hist_mc_reco and self.input_handler_gen_binning.hist_mc_reco.GetNbinsX() != n_gen_bins:
            raise ValueError("Mistmatch in input_handler_gen_binning.hist_mc_reco bins: %d vs %d" % (n_gen_bins, self.input_handler_gen_binning.hist_mc_reco.GetNbinsX()))
        if self.input_handler_gen_binning.hist_mc_fakes and self.input_handler_gen_binning.hist_mc_fakes.GetNbinsX() != n_gen_bins:
            raise ValueError("Mistmatch in input_handler_gen_binning.hist_mc_fakes bins: %d vs %d" % (n_gen_bins, self.input_handler_gen_binning.hist_mc_fakes.GetNbinsX()))

        input_status = self.SetInput(self.input_hist, bias_factor)

        # input_status = nError1+10000*nError2
        # nError1: number of bins where the uncertainty is zero.
        # these bins either are not used for the unfolding (if
        # oneOverZeroError==0) or 1/uncertainty is set to oneOverZeroError.</li>

        # nError2: return values>10000 are fatal errors, because the
        # unfolding can not be done. The number nError2 corresponds to the
        # number of truth bins which are not constrained by data points.
        n_error_1 = input_status % 10000
        n_error_2 = math.floor(input_status / 10000.)
        if n_error_2 > 0:
            if error_unconstrained_truth_bins:
                raise ValueError("%d unconstrained truth bins - cannot unfold" % n_error_2)
            else:
                warnings.warn("%d unconstrained truth bins - cannot unfold" % n_error_2)

    # bodge up old interface - FIXME
    @property
    def input_hist_original(self):
        return self.input_handler.input_hist

    @property
    def input_hist(self):
        return self.input_handler.input_hist

    @property
    def input_hist_bg_subtracted(self):
        return self.input_handler.input_hist_bg_subtracted

    @property
    def input_hist_gen_binning(self):
        return self.input_handler_gen_binning.input_hist

    @property
    def input_hist_gen_binning_bg_subtracted(self):
        return self.input_handler_gen_binning.input_hist_bg_subtracted

    @property
    def hist_truth(self):
        return self.input_handler.hist_truth

    @property
    def hist_mc_reco(self):
        return self.input_handler.hist_mc_reco

    @property
    def hist_mc_reco_bg_subtracted(self):
        return self.input_handler.hist_mc_reco_bg_subtracted

    @property
    def hist_mc_reco_gen_binning(self):
        return self.input_handler_gen_binning.hist_mc_reco

    @property
    def hist_mc_reco_gen_binning_bg_subtracted(self):
        return self.input_handler_gen_binning.hist_mc_reco_bg_subtracted

    @property
    def hist_mc_fakes(self):
        return self.input_handler.hist_mc_fakes

    @property
    def hist_mc_fakes_gen_binning(self):
        return self.input_handler_gen_binning.hist_mc_fakes


    def check_input(self, hist_y=None, raise_error=False):
        """Check the input for unconstrained output bins. Returns unconstrained output gen bins,
        along with their dependent reco bins

        Uses A^t * V_yy. Look for bins where the output bin is 0,
        because the sum(response_i * V_i) has 0s, and look for entries
        where response_i != 0 but V_i = 0

        Copied from TUnfold::SetInput(), but tidied up a bit

        Parameters
        ----------
        hist_y : TH1D, optional
            Hist to check, otherwise uses input Vyy from set_input()
        raise_error : bool, optional
            Raise Error if fails check (i.e. unconstrained output bins).
            Otherwise prints warning

        Returns
        -------
        list[int, list[int]]
            List of unconstrained [gen bins, [dependent reco bins]]

        Raises
        ------
        RuntimeError
            If hist_y is None, and set_input() wasn't already called
        ValueError
            If raise_error=True and we have unconstrained bins
        """
        if hist_y is None and self.GetVyy() is None:
            raise RuntimeError("hist_y and fVyy are empty, either specify hist_y or call SetInput() first")

        rowVyy1 = ROOT.std.vector('int')(self.GetNy())
        colVyy1 = ROOT.std.vector('int')(self.GetNy())
        dataVyy1 = ROOT.std.vector('double')(self.GetNy())
        dataVyyDiag = ROOT.std.vector('double')(self.GetNy())

        nVarianceZero = 0
        nVyy1 = 0

        for iy in range(self.GetNy()):
            # diagonals
            if hist_y is not None:
                dy = hist_y.GetBinError(iy + 1)
            else:
                dy = self.GetVyy()(iy, iy)
            dy2 = dy*dy
            if dy2 <= 0.0:
                nVarianceZero += 1

            rowVyy1[nVyy1] = iy
            colVyy1[nVyy1] = 0
            dataVyyDiag[iy] = dy2
            if dy2 > 0.0:
                dataVyy1[nVyy1] = dy2
                nVyy1 += 1

        vecV = self.CreateSparseMatrix(self.GetNy(), 1, nVyy1, rowVyy1, colVyy1, dataVyy1)

        # mAtV is a nx1 matrix
        mAtV = self.MultiplyMSparseTranspMSparse(self.GetA(), vecV)

        nError2 = 0
        for i in range(mAtV.GetNrows()):
            if mAtV.GetRowIndexArray()[i] == mAtV.GetRowIndexArray()[i + 1]:
                nError2 += 1

        ignored_input_bins = []
        if nError2 > 0:
            # check whether data points with zero error are responsible
            # a_rows[i] has the start index (in a_cols, a_data)
            # for row i, a_rows[i+1] is the end index + 1
            # a_cols[a_rows[i]]...a_cols[a_rows[i+1]-1]
            # is the set of column indices that have !=0 data
            # a_data[a_rows[i]]...a_data[a_rows[i+1]-1]
            # is the set of corresponding data
            fA = self.GetA()
            aT = ROOT.TMatrixDSparse(fA.GetNcols(), fA.GetNrows())
            aT.Transpose(fA)

            a_rows = aT.GetRowIndexArray()
            a_cols = aT.GetColIndexArray()
            a_data = aT.GetMatrixArray()

            # for each row in the resultant mAtV, look for 0 entries.
            # for each of those, look at all the A^T and dataVyyDiag pairs
            # (i.e. iterate over columns) and see which are 0
            for row in range(mAtV.GetNrows()):
                # printout every row's sum entries
                # printf("++++ ROW %d\n", row)
                # const sIndex = a_rows[row]
                # const eIndex = a_rows[row+1]
                # for (i=sIndex i<eIndex i++) {
                #     col = a_cols[i]
                #     cout << "dataVyyDiag[" << col << "]: " << dataVyyDiag[col] << endl
                #     printf("a^T(%d,%d) = %.4e\n", row, col, a_data[i])
                # }

                if mAtV.GetRowIndexArray()[row] == mAtV.GetRowIndexArray()[row + 1]:
                    binlist = ROOT.TString("no data to constrain output bin ")
                    binlist += self.GetOutputBinName(self.GetBinFromRow(row))
                    binlist += " depends on ignored input bins "
                    this_ignored_bin = [self.GetBinFromRow(row), []]

                    # Find all dependent entries
                    sIndex = a_rows[row]
                    eIndex = a_rows[row+1]
                    for i in range(sIndex, eIndex):
                        col = a_cols[i]
                        if dataVyyDiag[col] > 0.:
                            continue
                        if a_data[i] == 0:
                            continue
                        # cout << "dataVyyDiag[" << col << "]: " << dataVyyDiag[col] << endl
                        # printf("a^T(%d,%d) = %.4e\n", row, col, a_data[i])
                        binlist += " "
                        binlist += str(col+1) # +1 for ROOT binning
                        this_ignored_bin[1].append(col+1)

                    ignored_input_bins.append(this_ignored_bin)
                    warnings.warn("check_input: %s" % binlist)

            msg = "check_input: One output bins is not constrained by any data."
            if nError2 > 1:
                msg = "check_input: %d/%d output bins are not constrained by any data." % (nError2, mAtV.GetNrows())
            if raise_error:
                raise ValueError(msg)
            else:
                warnings.warn(msg)

        return ignored_input_bins

    def find_null_gen_bins(self, include_underflow=False, include_overflow=False):
        """Look for columns (gen bins) in response matrix with 0 entries, for removal
        (because y does not contain any information on that bin in x)
        """
        # quick way: get X projection (including y = 0 bin) and iterate throught that
        response_x_proj = self.response_map.ProjectionX("rsp_projX", 0, -1)
        null_gen_bins = []
        first_bin = 0 if include_underflow else 1
        last_bin = response_x_proj.GetNbinsX()+(2 if include_overflow else 1)
        for ix in range(first_bin, last_bin):
            val = response_x_proj.GetBinContent(ix)
            if val == 0:
                null_gen_bins.append(ix)
        return null_gen_bins

    def find_non_null_gen_bins(self, include_underflow=False, include_overflow=False):
        """Look for columns (gen bins) in response matrix with != 0 entries
        """
        # quick way: get X projection (including y = 0 bin) and iterate throught that
        response_x_proj = self.response_map.ProjectionX("rsp_projX", 0, -1)
        non_null_gen_bins = []
        first_bin = 0 if include_underflow else 1
        last_bin = response_x_proj.GetNbinsX()+(2 if include_overflow else 1)
        for ix in range(first_bin, last_bin):
            val = response_x_proj.GetBinContent(ix)
            if val != 0:
                non_null_gen_bins.append(ix)
        return non_null_gen_bins

    @staticmethod
    def make_bin_mapping(orig_bins, current_binning_scheme, new_binning_scheme):
        """Make map of global bin numbers, corresponding to same physical bin in new_binning_scheme

        Parameters
        ----------
        orig_bins : list[int]
            Description
        new_binning_scheme : PtVarBinnins, PtVarPerPtBinning
            Description
        """
        new_bin_mapping = {}
        for ob in orig_bins:
            # ob is a global bin in the current binning scheme
            physical_bin = current_binning_scheme.global_bin_to_physical_bin(ob)
            # physical bin is PtVar obj
            pt = physical_bin.pt[0]+1E-6
            var = physical_bin.var[0]+1E-6
            # map to overflow bin for this pt bin
            new_bin = new_binning_scheme.physical_bin_to_global_bin(pt=pt, var=var)
            # print('found', ob, '=', physical_bin, '->', new_bin, "=", new_binning_scheme.global_bin_to_physical_bin(new_bin))
            new_bin_mapping[ob] = new_bin
        return new_bin_mapping

    def make_null_bin_mapping(self, null_bins, new_binning_scheme=None):
        """Make mapping of {null global bin : overflow global bin}, for the correct pt bin

        Use new_binning_scheme to convert from physical value to global bin
        for the overflow bin
        """
        # figure out which pt bin we're in,
        # then get the index of the variable overflow bin in the pt bin
        # this way we can move any entries into that one
        new_binning_scheme = new_binning_scheme or self.binning_handler.get_binning_scheme('generator')
        if not new_binning_scheme.var_of:
            warnings.warn("new_binning_scheme doesn't have var_of=True, no idea what to do here")
            return {}
        null_bin_mapping = {}
        gen_binning = self.binning_handler.get_binning_scheme('generator')
        for nb in null_bins:
            # nb is a global bin in the current binning scheme
            physical_bin = gen_binning.global_bin_to_physical_bin(nb)
            # physical bin is PtVar obj
            pt = physical_bin.pt[0]+1E-6
            var = physical_bin.var[0]+1E-6
            # map to overflow bin for this pt bin
            new_bin = new_binning_scheme.physical_bin_to_global_bin(pt=pt, var=999999999)
            # print('found', nb, '=', physical_bin, '->', new_bin, "=", new_binning_scheme.global_bin_to_physical_bin(new_bin))
            null_bin_mapping[nb] = new_bin
        return null_bin_mapping


    def subtract_background(self, hist, name, scale=1.0, scale_err=0.0):
        """Subtract background source from input hist"""
        # Save into dict of components - needed? since TUnfoldDensity does this as well
        self.backgrounds[name] = hist.Clone()
        self.backgrounds[name].Scale(scale)
        # Also save total input subtracted
        self.input_handler.input_hist_bg_subtracted.Add(self.backgrounds[name], -1)
        # sanity check that none end up < 0
        for ix in range(1, self.input_hist_bg_subtracted.GetNbinsX()+1):
            val = self.input_hist_bg_subtracted.GetBinContent(ix)
            if val < 0:
                raise ValueError("self.input_hist_bg_subtracted bin %d has contents <0: %g" % (ix, val))
        self.SubtractBackground(hist, name, scale, scale_err) # NB TUnfold doesn't modify input_hist, nor hist

    def subtract_background_gen_binning(self, hist, name, scale=1.0):
        """Subtract background source with gen binning from input hist

        NB doesn't affect TUnfold, only for own book keeping
        """
        # Save into dict of components
        self.backgrounds_gen_binning[name] = hist.Clone()
        self.backgrounds_gen_binning[name].Scale(scale)
        # Also save total input subtracted
        self.input_handler_gen_binning.input_hist_bg_subtracted.Add(hist, -1*scale)

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

    def convert_reco_binned_hist_to_gen_binned(self, hist):
        """Convert a hist with detector level binning to gen level binning"""
        new_hist = self.generator_binning.CreateHistogram(cu.get_unique_str())
        # print('convert_reco_binned_hist_to_gen_binned: bin(1) low edge:', new_hist.GetBinLowEdge(1))
        # Iterate through all the finer reco bins, and for each determine
        # the corresponding gen bin, and add it to it from hist
        #
        # TODO: account for overflow in each axis?
        for pt in chain(self.pt_bin_edges_underflow_reco[:-1], self.pt_bin_edges_reco[:-1]):
            variable_bins = self.binning_handler.get_variable_bins(pt, 'detector')
            for var in variable_bins[:-1]:
                gen_bin = self.binning_handler.physical_bin_to_global_bin(pt=pt, var=var, binning_scheme='generator')
                det_bin = self.binning_handler.physical_bin_to_global_bin(pt=pt, var=var, binning_scheme='detector')
                # print("Converting bin", pt, var, ":", det_bin, "->", gen_bin)
                new_hist.SetBinContent(gen_bin, new_hist.GetBinContent(gen_bin) +  hist.GetBinContent(det_bin))
                new_hist.SetBinError(gen_bin, np.hypot(new_hist.GetBinError(gen_bin), hist.GetBinError(det_bin)))

        return new_hist

    def get_failed_reco(self, as_fraction=True):
        """Get hist of passGen & !passReco, ie 0th bin

        as_fraction = True, then divide by all passGen
        """
        origin = self.response_map
        # origin = self.get_probability_matrix()
        nbins = origin.GetNbinsX()
        ax = origin.GetXaxis()
        first_bin = ax.GetBinLowEdge(1)
        last_bin = ax.GetBinLowEdge(ax.GetNbins()+1)
        h_failed = ROOT.TH1D("h_failed_"+cu.get_unique_str(), "", nbins, first_bin, last_bin)
        for i in range(1, nbins+1):
            h_failed.SetBinContent(i, origin.GetBinContent(i, 0))
            h_failed.SetBinError(i, origin.GetBinError(i, 0))

        if as_fraction:
            h_all = origin.ProjectionX()
            h_failed.Divide(h_all)

        return h_failed

    # SETUP REGULARISATION STUFF
    # --------------------------------------------------------------------------
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
            raise ValueError("calculate_pt_bin_factors: 'which' arg should be 'unfolded' or 'gen'")

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
        pt_val = self.pt_bin_edges_gen[0]
        gen_binning = self.binning_handler.get_binning_scheme('generator')
        first_var = gen_binning.get_variable_bins(pt_val)[0]
        last_var = gen_binning.get_variable_bins(pt_val)[-1]
        start_bin = gen_binning.physical_bin_to_global_bin(var=first_var, pt=pt_val)
        end_bin = gen_binning.physical_bin_to_global_bin(var=last_var - 1E-6, pt=pt_val)  # need to be inside last bin - or include overflow?
        first_bin_integral = hist.Integral(start_bin, end_bin)  # ROOTs integral is inclusive of last bin

        bin_factors = {}
        # Add 1s for the 1st pt bin
        for var in gen_binning.get_variable_bins(pt_val)[:-1]:
            this_bin = gen_binning.physical_bin_to_global_bin(pt=pt_val, var=var)
            bin_factors[this_bin] = 1

        # Iterate through pt bins, figure out integral, scale according to first bin
        for pt_val in self.pt_bin_edges_gen[1:-1]:
            variable_bins = gen_binning.get_variable_bins(pt_val)
            start_bin = gen_binning.physical_bin_to_global_bin(var=variable_bins[0], pt=pt_val)
            end_bin = gen_binning.physical_bin_to_global_bin(var=variable_bins[-1] - 1E-6, pt=pt_val)
            integral = hist.Integral(start_bin, end_bin)
            sf = first_bin_integral / integral

            # Store bin factor for each lambda bin
            for var in variable_bins[:-1]:
                this_bin = gen_binning.physical_bin_to_global_bin(pt=pt_val, var=var)
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
                global_bin = self.generator_distribution_underflow.GetGlobalBinNumber(lambda_var+0.001, pt_val+0.001)
                lambda_width = self.variable_bin_edges_gen[lambda_ind+1] - self.variable_bin_edges_gen[lambda_ind]
                pt_width = self.pt_bin_edges_underflow_gen[pt_ind+1] - self.pt_bin_edges_underflow_gen[pt_ind]
                results[global_bin] = (lambda_width, pt_width)
            # signal region
            for pt_ind, pt_val in enumerate(self.pt_bin_edges_gen[:-1]):
                global_bin = self.generator_distribution.GetGlobalBinNumber(lambda_var+0.001, pt_val+0.001)
                lambda_width = self.variable_bin_edges_gen[lambda_ind+1] - self.variable_bin_edges_gen[lambda_ind]
                pt_width = self.pt_bin_edges_gen[pt_ind+1] - self.pt_bin_edges_gen[pt_ind]
                results[global_bin] = (lambda_width, pt_width)
        return results

    def setup_L_matrix_curvature(self, ref_hist, axis="both"):
        """Setup custom L matrix for curvature regularisation.

        ref_hist is the reference hist to determine factors for differential

        axis should be one of "both", "pt", or "angle",
        to determine along which axis/es the curvature is calculated.
        """
        gen_node = self.generator_distribution
        nr_counter = 0

        # unfolded_max = ref_hist.GetMaximum()

        L_matrix_entries = []

        # Add regularisation across pt bins, per lambda bin
        if axis in ["both", "pt"]:
            if isinstance(self.binning_handler.get_binning_scheme('generator'), PtVarPerPtBinning):
                warnings.warn(cu.pcolors.RED + "Cannot do setup_L_matrix_curvature() across pt bins with PtVarPerPtBinning" + cu.pcolors.ENDC)
            else:
                for ilambda in range(len(self.variable_bin_edges_gen[:-1])):
                    for ipt in range(len(self.pt_bin_edges_gen[:-3])):
                        pt_cen = self.pt_bin_edges_gen[ipt+1]+0.001  # add a tiny bit to make sure we're in the bin properly
                        lambda_cen = self.variable_bin_edges_gen[ilambda]+0.001

                        bin_ind_pt_down = gen_node.GetGlobalBinNumber(lambda_cen, self.pt_bin_edges_gen[ipt]+0.001)
                        bin_ind_pt_up = gen_node.GetGlobalBinNumber(lambda_cen, self.pt_bin_edges_gen[ipt+2]+0.001)

                        bin_ind_cen = gen_node.GetGlobalBinNumber(lambda_cen, pt_cen)

                        # bin_ind_var_down = gen_node.GetGlobalBinNumber(self.variable_bin_edges_gen[ilambda], pt_cen)
                        # bin_ind_var_up = gen_node.GetGlobalBinNumber(self.variable_bin_edges_gen[ilambda+2], pt_cen)

                        # print("Adding L matrix entry", nr_counter)
                        # print('lambda:', self.variable_bin_edges_gen[ilambda], 'pt:', (self.pt_bin_edges_gen[ipt], self.pt_bin_edges_gen[ipt+1], self.pt_bin_edges_gen[ipt+2]))

                        # pt_bin_width_down = pt_bin_edges_gen[ipt+1] - pt_bin_edges_gen[ipt]
                        # pt_bin_width_up = pt_bin_edges_gen[ipt+2] - pt_bin_edges_gen[ipt+1]
                        # factor = (pt_bin_width_down + pt_bin_width_up)
                        # value_pt_down = bin_factors[bin_ind_pt_down]
                        # value_pt_up = bin_factors[bin_ind_pt_up]
                        # ref_hist = unreg_self.unfolded

                        val_down = ref_hist.GetBinContent(bin_ind_pt_down)
                        value_pt_down = 1./val_down if val_down != 0 else 0

                        val_up = ref_hist.GetBinContent(bin_ind_pt_up)
                        value_pt_up = 1./val_up if val_up != 0 else 0

                        # value_pt_down = bin_factors[bin_ind_pt_down]
                        # value_pt_up = bin_factors[bin_ind_pt_up]
                        value_pt_cen = - (value_pt_down + value_pt_up)

                        val_cen = ref_hist.GetBinContent(bin_ind_cen)
                        value_pt_cen = -2. / val_cen if val_cen != 0 else 0

                        # print(bin_ind_pt_down, value_pt_down, bin_ind_cen, value_pt_cen, bin_ind_pt_up, value_pt_up)
                        L_args = [bin_ind_pt_down, value_pt_down, bin_ind_cen, value_pt_cen, bin_ind_pt_up, value_pt_up]
                        L_matrix_entries.append(L_args)
                        self.AddRegularisationCondition(*L_args)
                        nr_counter += 1
                        print(L_args)

                        # value_pt_down = unfolded_max/ref_hist.GetBinContent(bin_ind_pt_down)
                        # value_pt_up = unfolded_max/ref_hist.GetBinContent(bin_ind_pt_up)
                        # value_var_down = unfolded_max/ref_hist.GetBinContent(bin_ind_var_down)
                        # value_var_up = unfolded_max/ref_hist.GetBinContent(bin_ind_var_up)
                        # value_cen = - (value_pt_down + value_pt_up + value_var_down + value_var_up)
                        # print(bin_ind_pt_down, value_pt_down, bin_ind_cen, value_cen, bin_ind_pt_up, value_pt_up)
                        # indices = [bin_ind_pt_down, bin_ind_var_down, bin_ind_cen, bin_ind_pt_up, bin_ind_var_up]
                        # row_data = [value_pt_down, value_var_down, value_cen, value_pt_up, value_var_up]
                        # self.AddRegularisationCondition(5, array('i', indices), array('d', row_data))
                        # print(indices, row_data)

        # Add regularisation across lambda bins, per pt bin
        if axis in ["both", "angle"]:
            for ipt in range(len(self.pt_bin_edges_gen[:-1])):
                for ilambda in range(len(self.variable_bin_edges_gen[:-3])):
                    pt_cen = self.pt_bin_edges_gen[ipt]+0.001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)
                    lambda_cen = self.variable_bin_edges_gen[ilambda+1]+0.001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)

                    bin_ind_lambda_down = gen_node.GetGlobalBinNumber(self.variable_bin_edges_gen[ilambda]+0.001, pt_cen)
                    bin_ind_lambda_up = gen_node.GetGlobalBinNumber(self.variable_bin_edges_gen[ilambda+2]+0.001, pt_cen)

                    bin_ind_cen = gen_node.GetGlobalBinNumber(lambda_cen, pt_cen)

                    # bin_ind_var_down = gen_node.GetGlobalBinNumber(self.variable_bin_edges_gen[ilambda], pt_cen)
                    # bin_ind_var_up = gen_node.GetGlobalBinNumber(self.variable_bin_edges_gen[ilambda+2], pt_cen)

                    # print("Adding L matrix entry", nr_counter)
                    # print('pt:', self.pt_bin_edges_gen[ipt], 'lambda:', (self.variable_bin_edges_gen[ilambda], self.variable_bin_edges_gen[ilambda+1], self.variable_bin_edges_gen[ilambda+2]))

                    # pt_bin_width_down = variable_bin_edges_gen[ilambda+1] - variable_bin_edges_gen[ilambda]
                    # pt_bin_width_up = variable_bin_edges_gen[ilambda+2] - variable_bin_edges_gen[ilambda+1]
                    # factor = (pt_bin_width_down + pt_bin_width_up)
                    # value_lambda_down = bin_factors[bin_ind_lambda_down]
                    # value_lambda_up = bin_factors[bin_ind_lambda_up]
                    # ref_hist = unreg_self.unfolded

                    val_down = ref_hist.GetBinContent(bin_ind_lambda_down)
                    value_lambda_down = 1./val_down if val_down != 0 else 0

                    val_up = ref_hist.GetBinContent(bin_ind_lambda_up)
                    value_lambda_up = 1./val_up if val_up != 0 else 0

                    # value_lambda_down = bin_factors[bin_ind_lambda_down]
                    # value_lambda_up = bin_factors[bin_ind_lambda_up]
                    # value_lambda_cen = - (value_lambda_down + value_lambda_up)

                    val_cen = ref_hist.GetBinContent(bin_ind_cen)
                    value_lambda_cen = -2. / val_cen if val_cen != 0 else 0

                    # print(bin_ind_lambda_down, value_lambda_down, bin_ind_cen, value_lambda_cen, bin_ind_lambda_up, value_lambda_up)
                    L_args = [bin_ind_lambda_down, value_lambda_down, bin_ind_cen, value_lambda_cen, bin_ind_lambda_up, value_lambda_up]
                    L_matrix_entries.append(L_args)
                    self.AddRegularisationCondition(*L_args)
                    nr_counter += 1
                    print(L_args)

        self.L_matrix_entries = L_matrix_entries

    def setup_L_matrix_derivative(self, ref_hist, axis="both"):
        """Setup custom L matrix for derivative regularisation.

        ref_hist is the reference hist to determine factors for differential

        axis should be one of "both", "pt", or "angle",
        to determine along which axis/es the curvature is calculated.
        """
        gen_node = self.generator_distribution
        nr_counter = 0

        # unfolded_max = ref_hist.GetMaximum()

        L_matrix_entries = []

        # Add regularisation across pt bins, per lambda bin
        if axis in ["both", "pt"]:
            if isinstance(self.binning_handler.get_binning_scheme('generator'), PtVarPerPtBinning):
                warnings.warn(cu.pcolors.RED + "Cannot do setup_L_matrix_derivative() across pt bins with PtVarPerPtBinning" + cu.pcolors.ENDC)
            else:
                for ilambda in range(len(self.variable_bin_edges_gen[:-1])):
                    for ipt in range(len(self.pt_bin_edges_gen[:-2])):
                        lambda_cen = self.variable_bin_edges_gen[ilambda]+0.001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)

                        bin_ind_pt_down = gen_node.GetGlobalBinNumber(lambda_cen, self.pt_bin_edges_gen[ipt]+0.001)
                        bin_ind_pt_up = gen_node.GetGlobalBinNumber(lambda_cen, self.pt_bin_edges_gen[ipt+1]+0.001)

                        # print("Adding L matrix entry", nr_counter)
                        # print('lambda:', self.variable_bin_edges_gen[ilambda], 'pt:', (self.pt_bin_edges_gen[ipt], self.pt_bin_edges_gen[ipt+1]))

                        val_down = ref_hist.GetBinContent(bin_ind_pt_down)
                        value_pt_down = 1./val_down if val_down != 0 else 0

                        val_up = ref_hist.GetBinContent(bin_ind_pt_up)
                        value_pt_up = 1./val_up if val_up != 0 else 0

                        L_args = [bin_ind_pt_down, -value_pt_down, bin_ind_pt_up, value_pt_up]
                        L_matrix_entries.append(L_args)
                        self.AddRegularisationCondition(*L_args)
                        nr_counter += 1
                        print(L_args)

        # Add regularisation across lambda bins, per pt bin
        if axis in ["both", "angle"]:
            for ipt in range(len(self.pt_bin_edges_gen[:-1])):
                for ilambda in range(len(self.variable_bin_edges_gen[:-2])):
                    pt_cen = self.pt_bin_edges_gen[ipt]+0.001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)

                    bin_ind_lambda_down = gen_node.GetGlobalBinNumber(self.variable_bin_edges_gen[ilambda]+0.001, pt_cen)
                    bin_ind_lambda_up = gen_node.GetGlobalBinNumber(self.variable_bin_edges_gen[ilambda+1]+0.001, pt_cen)

                    # print("Adding L matrix entry", nr_counter)
                    # print('pt:', self.pt_bin_edges_gen[ipt], 'lambda:', (self.variable_bin_edges_gen[ilambda], self.variable_bin_edges_gen[ilambda+1]))

                    val_down = ref_hist.GetBinContent(bin_ind_lambda_down)
                    value_lambda_down = 1./val_down if val_down != 0 else 0

                    val_up = ref_hist.GetBinContent(bin_ind_lambda_up)
                    value_lambda_up = 1./val_up if val_up != 0 else 0

                    L_args = [bin_ind_lambda_down, -value_lambda_down, bin_ind_lambda_up, value_lambda_up]
                    L_matrix_entries.append(L_args)
                    self.AddRegularisationCondition(*L_args)
                    nr_counter += 1
                    print(L_args)

        self.L_matrix_entries = L_matrix_entries

    # HANDLE SYSTEMATIC UNCERTAINTIES
    # --------------------------------------------------------------------------
    def get_all_exp_syst_labels(self):
        return [x.label for x in self.exp_systs]

    def has_exp_syst(self, label):
        items = [x for x in self.exp_systs if x.label == label]
        if len(items) == 0:
            return False
        if len(items) > 1:
            raise ValueError("Found >1 exp systematic with label '%s': %s" % (label, [x.label for x in items]))
        return True

    def get_exp_syst(self, label):
        items = [x for x in self.exp_systs if x.label == label]
        if len(items) == 0:
            raise ValueError("Found no exp systematic with label '%s', only have: %s" % (label, [x.label for x in self.exp_systs]))
        if len(items) > 1:
            raise ValueError("Found >1 exp systematic with label '%s': %s" % (label, [x.label for x in items]))
        return items[0]

    def add_sys_error(self, map_syst, label, mode=ROOT.TUnfoldDensity.kSysErrModeMatrix):
        """Add systematic error via response map, arguments as per AddSysError()"""
        # don't store map, since we don't actually use it, and it wastes memory
        this_syst = ExpSystematic(label=label, syst_map=None)
        self.AddSysError(map_syst, label, self.orientation, mode)
        self.exp_systs.append(this_syst)

    def get_syst_shift(self, syst_label):
        """Get shift in result due to a particular systeamtic

        Label must be same as used to add it in add_sys_error()
        """
        # TODO: syst_label -> ExpSystematic obj
        print("get_syst_shift('%s')" % syst_label)
        this_syst = self.get_exp_syst(syst_label)
        if this_syst.syst_shift is None:
            hist = self.GetDeltaSysSource(this_syst.label,
                                          this_syst.syst_shift_label,  # name given to hist obj
                                          "",
                                          self.output_distribution_name,  # must be the same as what's used in get_output
                                          self.axisSteering,
                                          self.use_axis_binning)
            this_syst.syst_shift = hist  # cache shifts
        return this_syst.syst_shift

    def get_syst_shifted_hist(self, syst_label, unfolded=None):
        """Get histogram with systematic shift applied to bin contents

        Can specify starting hist, otherwise assumes unfolded w/no error
        """
        # TODO: syst_label -> ExpSystematic obj
        this_syst = self.get_exp_syst(syst_label)
        if this_syst.syst_shifted is None:
            hist_shift = self.get_syst_shift(syst_label).Clone(this_syst.syst_shifted_label)
            unfolded = unfolded or self.get_unfolded_with_no_errors()
            hist_shift.Add(unfolded)  # TODO what about errors?
            this_syst.syst_shifted = hist_shift
        return this_syst.syst_shifted

    def get_syst_error_hist(self, syst_label, unfolded=None):
        """Get histogram with systematic shift as error bars

        Can specify starting hist, otherwise assumes unfolded w/no error
        """
        this_syst = self.get_exp_syst(syst_label)
        if this_syst.syst_error_bar is None:
            ref = unfolded or self.get_unfolded_with_no_errors()
            new_hist = self.convert_error_shift_to_error_bars(ref, self.get_syst_shift(syst_label))
            this_syst.syst_error_bar = new_hist
        return this_syst.syst_error_bar

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
            print("self.unfolded is a", type(self.unfolded), "and has", self.unfolded.GetNbinsX(), "x bins, bin1 low edge:", self.unfolded.GetBinLowEdge(1))

            # Sanity check incase unfolding went really wrong - e.g. if rank
            # of matrix E is way smaller than what is expected
            # May mean SetEpsMatrix() is needed to help inversion
            if 0 < max_chi2_ndf < (self.chi2sys / self.Ndf):
                raise RuntimeError("chi2sys / Ndf > %d - unfolding is rubbish! Maybe use SetEpsMatrix()? "
                                   "Especially if you get 'rank of matrix E X expect Y' warning" % max_chi2_ndf)

        return self.unfolded

    def inspect_output(self):
        for ix in range(1, self.get_output().GetNbinsX()+1):
            val = self.get_output().GetBinContent(ix)
            if val < 0:
                err = self.get_output().GetBinError(ix)
                msg = "Unfolded bin %d = %s has value %f +- %f" % (ix, self.binning_handler.global_bin_to_physical_bin(ix, 'generator'), val, err)
                print(cu.pcolors.OKCYAN + msg + cu.pcolors.ENDC)

    def _post_process(self):
        """Do some standard things & store various things that are done after unfolding"""
        self.get_ematrix_input()
        self.get_ematrix_stat_response()
        self.get_ematrix_stat()
        # self.get_ematrix_tau()
        self.get_ematrix_tunfold_total()
        self.get_rhoij_total()
        self.get_probability_matrix()
        self.update_unfolded_with_ematrix_tunfold_total()
        self.get_unfolded_with_ematrix_stat()
        self.get_unfolded_with_ematrix_response()
        self.get_folded_unfolded()
        self.get_folded_mc_truth()
        # these 3 for bottom-line tests
        self.get_vxx_th2()
        self.get_vxx_inv_th2()
        # self.get_vyy_no_bg_subtraction_th2()
        # self.get_vyy_inv_no_bg_subtraction_th2()

        for exp_syst in self.exp_systs:
            # setup all the internal maps
            self.get_syst_shift(exp_syst.label)
            self.get_ematrix_syst(exp_syst.label)
            self.get_syst_shifted_hist(exp_syst.label)
            self.get_syst_error_hist(exp_syst.label)

        self.hist_bin_chopper.add_obj('hist_truth', self.hist_truth)
        self.hist_bin_chopper.add_obj('unfolded', self.get_output())
        self.hist_bin_chopper.add_obj(self.stat_uncert_name, self.get_unfolded_with_ematrix_stat())
        self.hist_bin_chopper.add_obj(self.rsp_uncert_name, self.get_unfolded_with_ematrix_response())

    @staticmethod
    def make_hist_from_diagonal_errors(h2d, do_sqrt=True, set_errors=True, offset=0.):
        """Make 1D hist, with errors or contents set to diagonal elements from h2d

        Can be TH2 or numpy.ndarray, cos we have to use both
        Yes that is majorly wack

        set_errors: True to set TH1 bin errors, otherwise sets bin contents
        offset is on bin edge from 0 (TUnfold is 0.5)
        """
        if isinstance(h2d, ROOT.TH2):
            nbins = h2d.GetNbinsX()
            hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, offset, nbins+offset)
            for i in range(1, nbins+1):
                err = h2d.GetBinContent(i, i)
                if do_sqrt and err > 0:
                    err = math.sqrt(err)
                if set_errors:
                    hnew.SetBinContent(i, 0)
                    hnew.SetBinError(i, err)
                else:
                    hnew.SetBinContent(i, err)
                    hnew.SetBinError(i, 0)
            return hnew
        elif isinstance(h2d, np.ndarray):
            nbins = h2d.shape[0]
            hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, offset, nbins+offset)
            for i in range(1, nbins+1):
                err = h2d[i-1, i-1]
                if do_sqrt and err > 0:
                    err = math.sqrt(err)
                if set_errors:
                    hnew.SetBinContent(i, 0)
                    hnew.SetBinError(i, err)
                else:
                    hnew.SetBinContent(i, err)
                    hnew.SetBinError(i, 0)
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

    def update_unfolded_with_ematrix_tunfold_total(self):
        """Update unfolded hist with total errors from total error matrix"""
        error_total_1d = self.make_hist_from_diagonal_errors(self.get_ematrix_tunfold_total(), do_sqrt=True)  # note that bin contents = 0, only bin errors are non-0
        self.update_hist_bin_error(h_orig=error_total_1d, h_to_be_updated=self.get_output())

    def get_unfolded_with_no_errors(self):
        """Make copy of unfolded, but with no error bars"""
        if getattr(self, 'unfolded_no_err', None) is None:
            self.unfolded_no_err = self.get_output().Clone("unfolded_no_err")
            cu.remove_th1_errors(self.unfolded_no_err)
        return self.unfolded_no_err

    def get_unfolded_with_ematrix_stat(self):
        """Make copy of unfolded, but only stat errors"""
        if getattr(self, 'unfolded_stat_err', None) is None:
            error_1d = self.make_hist_from_diagonal_errors(self.get_ematrix_stat(), do_sqrt=True)  # note that bin contents = 0, only bin errors are non-0
            self.unfolded_stat_err = self.get_output().Clone("unfolded_stat_err")
            self.update_hist_bin_error(h_orig=error_1d, h_to_be_updated=self.unfolded_stat_err)
        return self.unfolded_stat_err

    def get_unfolded_with_ematrix_response(self):
        """Create unfolded with error bars from response matrix uncertainties"""
        if getattr(self, 'unfolded_rsp_err', None) is None:
            error_stat_response = self.make_hist_from_diagonal_errors(self.get_ematrix_stat_response(), do_sqrt=True)  # note that bin contents need to be correct, otherwise won't normalise correctly
            self.unfolded_rsp_err = self.get_output().Clone("unfolded_rsp_err")
            self.update_hist_bin_error(h_orig=error_stat_response, h_to_be_updated=self.unfolded_rsp_err)
        return self.unfolded_rsp_err

    def get_unfolded_with_ematrix_total(self):
        """Create unfolded with error bars from total uncertainties"""
        return self.get_syst_error_hist('Total')
        if getattr(self, 'unfolded_total_err', None) is None:
            error_stat_response = self.make_hist_from_diagonal_errors(self.get_ematrix_total_absolute(), do_sqrt=True)  # note that bin contents need to be correct, otherwise won't normalise correctly
            self.unfolded_total_err = self.get_output().Clone("unfolded_total_err")
            self.update_hist_bin_error(h_orig=error_stat_response, h_to_be_updated=self.unfolded_total_err)
        return self.unfolded_total_err

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

    def get_ndarray_signal_region_no_overflow(self, arr, xbinning='generator', ybinning='generator'):
        """Generic method to create matrix that removes pt underflow bins,
        and only signal region pt/lambda bins (i.e no over/underflow)

        very slow as we have to manually iterate over all the bins
        I'm sure there's a quicker way

        Parameters
        ----------
        arr : ROOT.TH2, TH1, or numpy.ndarray

        Returns
        -------
        numpy.ndarray
        """
        if not isinstance(arr, (ROOT.TH2, ROOT.TH1, np.ndarray)):
            raise TypeError("arr should be TH2/1 or ndarray")

        is_1d = ((isinstance(arr, ROOT.TH1) and not isinstance(arr, ROOT.TH2)) or  # important otherwise TH2 will eval True, as subclass of TH1
                 (isinstance(arr, np.ndarray) and arr.shape[0] == 1) or
                 (ybinning is None))

        if ybinning is None:
            ybinning = xbinning # dummy value

        def get_num_bins(binning_scheme):
            n_bins = 0
            binning_obj = self.binning_handler.get_binning_scheme(binning_scheme)
            for pt_bin in binning_obj.get_pt_bins(is_signal_region=True)[:-1]:
                var_bins = binning_obj.get_variable_bins(pt_bin)
                # print("For pt_bin", pt_bin, "I have var bins", var_bins, "adding", len(var_bins)-1, "to the sum")
                n_bins += len(var_bins)-1
            return n_bins

        if is_1d:
            nrows = 1
        else:
            # nrows = self.nbins_pt_gen*self.nbins_variable_gen if ybinning == 'generator' else self.nbins_pt_reco*self.nbins_variable_reco
            nrows = get_num_bins(ybinning)
        # ncols = self.nbins_pt_gen*self.nbins_variable_gen if xbinning == 'generator' else self.nbins_pt_reco*self.nbins_variable_reco
        ncols = get_num_bins(xbinning)
        output_arr = np.zeros(shape=(nrows, ncols))

        x_index = -1
        # x_distribution = self.generator_distribution if xbinning == 'generator' else self.detector_distribution
        # pt_bin_edges_x = self.pt_bin_edges_gen if xbinning == 'generator' else self.pt_bin_edges_reco
        pt_bin_edges_x = self.binning_handler.get_pt_bins(binning_scheme=xbinning, is_signal_region=True)

        # y_distribution = self.generator_distribution if ybinning == 'generator' else self.detector_distribution
        # pt_bin_edges_y = self.pt_bin_edges_gen if ybinning == 'generator' else self.pt_bin_edges_reco
        pt_bin_edges_y = self.binning_handler.get_pt_bins(binning_scheme=ybinning, is_signal_region=True)

        for pt_ind_x, pt_x in enumerate(pt_bin_edges_x[:-1]):
            variable_bin_edges_x = self.binning_handler.get_variable_bins(pt_x, xbinning)
            for var_ind_x, var_x in enumerate(variable_bin_edges_x[:-1]):
                global_bin_x = self.binning_handler.physical_bin_to_global_bin(pt=pt_x, var=var_x, binning_scheme=xbinning)
                # global_bin_x = x_distribution.GetGlobalBinNumber(var_x+0.001, pt_x+0.001)
                x_index += 1
                if is_1d:
                    if isinstance(arr, ROOT.TH1):
                        output_arr[0, x_index] = arr.GetBinContent(global_bin_x)
                    else:
                        output_arr[0, x_index] = arr[0, global_bin_x-1]  # -1 as ndarray 0-index, global bin 1-indexed
                else:
                    y_index = -1
                    for pt_ind_y, pt_y in enumerate(pt_bin_edges_y[:-1]):
                        variable_bin_edges_y = self.binning_handler.get_variable_bins(pt_y, ybinning)
                        for var_ind_y, var_y in enumerate(variable_bin_edges_y[:-1]):
                            global_bin_y = self.binning_handler.physical_bin_to_global_bin(pt=pt_y, var=var_y, binning_scheme=ybinning)
                            # global_bin_y = y_distribution.GetGlobalBinNumber(var_y+0.001, pt_y+0.001)
                            y_index += 1
                            if isinstance(arr, ROOT.TH2):
                                output_arr[y_index, x_index] = arr.GetBinContent(global_bin_x, global_bin_y)
                            else:
                                output_arr[y_index, x_index] = arr[global_bin_y-1, global_bin_x-1]  # -1 as ndarray 0-index, global bin 1-indexed
        return output_arr

    # LOTS OF COVARIANCE MATRIX FUNCTIONS
    # --------------------------------------------------------------------------
    def get_ematrix_input(self):
        """Get error matrix due to statistics from thing being unfolded"""
        if getattr(self, "ematrix_input", None) is None:
            self.ematrix_input = self.GetEmatrixInput("ematrix_input_"+cu.get_unique_str(), "", self.output_distribution_name, "*[]", self.use_axis_binning)
            # print("ematrix_input1st bin low edge:", self.ematrix_input.GetXaxis().GetBinLowEdge(1))
            # print("output bin low edge:", self.get_output().GetXaxis().GetBinLowEdge(1))
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

    def get_ematrix_tunfold_total(self):
        """Total error matrix from TUnfold, from stat+systs you gave it
        - doesn't include extras like scale or PDF"""
        if getattr(self, "ematrix_tunfold_total", None) is None:
            self.ematrix_tunfold_total = self.GetEmatrixTotal("ematrix_tunfold_total_"+cu.get_unique_str(), "", self.output_distribution_name, "*[]", self.use_axis_binning)
            print("ematrix_tunfold_total is:", type(self.ematrix_tunfold_total), "with #xbins:", self.ematrix_tunfold_total.GetNbinsX(), "bin1 low edge:", self.ematrix_tunfold_total.GetXaxis().GetBinLowEdge(1))
        return self.ematrix_tunfold_total

    @property
    def ematrix_tunfold_total_ndarray(self):
        cached_attr_name = '_ematrix_total_ndarray'
        if not hasattr(self, cached_attr_name):
            arr, _ = cu.th2_to_ndarray(self.get_ematrix_tunfold_total())
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
            print("ematrix stat has shape", self.ematrix_stat.GetNbinsX(), self.ematrix_stat.GetNbinsY())
        return self.ematrix_stat

    def get_ematrix_total_absolute(self):
        return self.get_ematrix_syst("Total")

    def get_ematrix_total_normalised(self):
        return self.get_ematrix_syst("Total_Norm")

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
                return None
        return this_syst.syst_ematrix

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

    def get_ematrix_tunfold_total_inv(self):
        """Total error matrix inverted, from stat+systs"""
        return self.InvertMSparseSymmPos(self.GetSummedErrorMatrixXX(), False)  # TMatrixDSparse

    def get_vyy_inv_ndarray(self):
        """Get inverse of V_yy (ie input). Note done after BG-subtraction"""
        if getattr(self, 'vyy_inv_ndarray', None) is None:
            if getattr(self, 'vyy_inv_tmatrix', None) is None:
                self.vyy_inv_tmatrix = self.GetVyyInv()
            self.vyy_inv_ndarray = cu.tmatrixdsparse_to_ndarray(self.vyy_inv_tmatrix)
        return self.vyy_inv_ndarray

    def get_vyy_no_bg_subtraction_th2(self):
        """Same as get_vyy_inv_ndarray() but before BG-subtraction"""
        if getattr(self, 'vyy_no_bg_subtraction_th2', None) is None:
            self.vyy_no_bg_subtraction_th2 = self.make_diag_cov_hist_from_errors(self.input_hist, inverse=False)
        return self.vyy_no_bg_subtraction_th2

    def get_vyy_no_bg_subtraction_ndarray(self):
        # before bg-subtraction
        if getattr(self, 'vyy_no_bg_subtraction_ndarray', None) is None:
            self.vyy_no_bg_subtraction_ndarray, _ = cu.th2_to_ndarray(self.get_vyy_no_bg_subtraction_th2())
        return self.vyy_no_bg_subtraction_ndarray

    def get_vyy_inv_no_bg_subtraction_th2(self):
        # before bg-subtraction
        if getattr(self, 'vyy_inv_no_bg_subtraction_th2', None) is None:
            self.vyy_inv_no_bg_subtraction_th2 = self.make_diag_cov_hist_from_errors(self.input_hist, inverse=True)
            # this_binning = self.detector_binning.FindNode('detector')
            # self.vyy_inv_no_bg_subtraction_th2 = this_binning.CreateErrorMatrixHistogram("ematrix_vyyinv_no_bg_subtraction_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
            # for i in range(1, self.input_hist.GetNbinsX()+1):
            #     bin_err = self.input_hist.GetBinError(i)
            #     new_err = 1./(bin_err*bin_err) if bin_err != 0 else 0
            #     self.vyy_inv_no_bg_subtraction_th2.SetBinContent(i, i, new_err)
        return self.vyy_inv_no_bg_subtraction_th2

    def get_vyy_inv_no_bg_subtraction_ndarray(self):
        if getattr(self, 'vyy_inv_no_bg_subtraction_ndarray', None) is None:
            self.vyy_inv_no_bg_subtraction_ndarray, _ = cu.th2_to_ndarray(self.get_vyy_inv_no_bg_subtraction_th2())
        return self.vyy_inv_no_bg_subtraction_ndarray

    def get_vxx_th2(self):
        if getattr(self, 'vxx_th2', None) is None:
            print("Getting vxx_th2")
            # Have to manually create hist first, awkward
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            # Gnahhhhhhh
            self.vxx_th2 = this_binning.CreateErrorMatrixHistogram("ematrix_vxx_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
            self.ErrorMatrixToHist(self.vxx_th2, self.GetVxx())
            print('vxx_th2.GetXaxis().GetBinLowEdge(1):', self.vxx_th2.GetXaxis().GetBinLowEdge(1))
        return self.vxx_th2

    def get_vxx_ndarray(self):
        if getattr(self, 'vxx_ndarray', None) is None:
            self.vxx_ndarray, _ = cu.th2_to_ndarray(self.get_vxx_th2())
        return self.vxx_ndarray

    def get_vxx_inv_th2(self):
        if getattr(self, 'vxx_inv_th2', None) is None:
            print("Getting vxx_inv_th2")
            # Have to manually create hist first, awkward
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            # Gnahhhhhhh
            self.vxx_inv_th2 = this_binning.CreateErrorMatrixHistogram("ematrix_vxxinv_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
            self.ErrorMatrixToHist(self.vxx_inv_th2, self.GetVxxInv())
            print('vxx_inv_th2.GetXaxis().GetBinLowEdge(1):', self.vxx_inv_th2.GetXaxis().GetBinLowEdge(1))
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
    def calculate_singular_max_min(matrix, non_zero=False):
        """Calculate max & min singular values condition number as per StatsComm guidelines

        These are found using TDecompSVD.
        (we ignore the builtin condition() method as it calculates it differently)
        """
        print("Calculating singular values on matrix size (nrows, ncols): (%d, %d)" % (matrix.GetNrows(), matrix.GetNcols()))
        if matrix.GetNcols() > matrix.GetNrows():
            raise RuntimeError("Condition number only for matrix where # rows >= # cols")

        svd = ROOT.TDecompSVD(matrix)
        sig = svd.GetSig()  # by construction, ordered descending
        sigma_max = sig[0]
        sigma_min = max(0, sig[sig.GetNrows()-1])
        print("sigma_max:", sigma_max, "sigma_min:", sigma_min)
        if sigma_min == 0:
            print("sigma_min > 0:", min([x for x in sig if x > 0]))
            if non_zero:
                sigma_min = min([x for x in sig if x > 0])
        return sigma_max, sigma_min

    def print_condition_number(self, remove_underflow_bins=False, remove_overflow_bins=False):
        """Store & print response matrix condition number and some advice

        Defined as sigma_max / max(0, sigma_min), where sigma_{max/min} are the
        largest/smallest singular values.
        These are also stored for later usage if needed (since expensive to calc)
        """
        if getattr(self, 'condition_number', None) is None:
            sigma_max, sigma_min = self.calculate_singular_max_min(cu.th2_to_tmatrixd(self.get_probability_matrix()), non_zero=True)
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

        # try a version cutting out underflow bins
        if remove_underflow_bins:
            start_ind_gen = self.generator_distribution.GetStartBin() - 1  # -1 as numpy starts at 0, TUnfold at 1
            start_ind_det = self.detector_distribution.GetStartBin() - 1  # -1 as numpy starts at 0, TUnfold at 1
            print("Signal Start ind gen:", start_ind_gen)
            print("Signal Start ind det:", start_ind_det)
            print("prob ndarray.shape:", self.probability_ndarray.shape)
            probability_ndarray_signal_pt = np.copy(self.probability_ndarray)[start_ind_det:, start_ind_gen:]
            sigma_max, sigma_min = self.calculate_singular_max_min(cu.th2_to_tmatrixd(cu.ndarray_to_th2(probability_ndarray_signal_pt)), non_zero=True)
            if sigma_min == 0:
                # avoid DivisionError
                print("Minmum singular value = 0, condition number = Infinity")
                num = np.inf
            else:
                num = sigma_max / sigma_min
            print("Condition number without pt < 50 bins:", num)

        # now also do a version removing the overflow bins in each variable & pT:
        if remove_overflow_bins:
            probability_ndarray_no_oflow = self.get_ndarray_signal_region_no_overflow(self.get_probability_matrix(), xbinning='generator', ybinning='detector')
            sigma_max, sigma_min = self.calculate_singular_max_min(cu.th2_to_tmatrixd(cu.ndarray_to_th2(probability_ndarray_no_oflow)), non_zero=True)
            if sigma_min == 0:
                # avoid DivisionError
                print("Minmum singular value = 0, condition number = Infinity")
                num = np.inf
            else:
                num = sigma_max / sigma_min
            print("Condition number without pt < 50 bins and without overflow bins:", num)

    def get_response_normed_by_detector_pt(self):
        if getattr(self, 'response_map_normed_by_detector_pt', None) is None:
            normed_response_map = self.response_map.Clone("response_map_normed_by_detector_pt")
            sums = []
            pt_bins = list(self.pt_bin_edges_underflow_reco) + list(self.pt_bin_edges_reco)[1:]
            binning_scheme = self.binning_handler.get_binning_scheme("detector")
            if not binning_scheme.pt_of:
                # if no pt overflow, no extra bins
                pt_bins = pt_bins[:-1]

            # print(bins)
            # First get the sum over all bins corresponding to this reco pT bin
            # Means summing over a number of reco lambda bins, and all generator bins
            for ibin_pt, (pt, pt_next) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
                # print(ibin_pt, pt, pt_next)
                this_sum = 0
                var = binning_scheme.get_variable_bins(pt)[0]+0.001
                global_bin_reco = binning_scheme.physical_bin_to_global_bin(var=var, pt=pt+0.001)
                global_bin_reco_next = binning_scheme.physical_bin_to_global_bin(var=var, pt=pt_next+0.001)

                for i_reco in range(global_bin_reco, global_bin_reco_next):
                    for global_bin_gen in range(self.response_map.GetNbinsX()+1):
                        this_sum += self.response_map.GetBinContent(global_bin_gen, i_reco)
                sums.append(this_sum)
            # print(sums)

            # Now set the new bin contents and error bars by scaling using these sums
            for ibin_pt, (pt, pt_next) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
                var = binning_scheme.get_variable_bins(pt)[0]+0.001
                global_bin_reco = binning_scheme.physical_bin_to_global_bin(var=var, pt=pt+0.001)
                global_bin_reco_next = binning_scheme.physical_bin_to_global_bin(var=var, pt=pt_next+0.001)
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
            self.folded_unfolded = cu.ndarray_to_th1(folded_vec.T, has_oflow_x=oflow, offset=0.5)

            # Error propagation: if y = Ax, with covariance matrices Vyy and Vxx,
            # respectively, then Vyy = (A*Vxx)*A^T
            unfolded_covariance_matrix, _ = cu.th2_to_ndarray((self.get_ematrix_tunfold_total()), oflow_x=oflow, oflow_y=oflow)
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
        folded_mc_truth = cu.ndarray_to_th1(folded_vec.T, has_oflow_x=oflow, offset=0.5)
        print('folded_mc_truth bin(1) low edge:', folded_mc_truth.GetBinLowEdge(1))

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

        y_ndarray = np.zeros(shape=(self.GetNy(), 1))
        y = self.GetY()
        for i in range(self.GetNy()):
            y_ndarray[i, 0] = y(i, 0)

        # calculate my own result simple inversion
        # nb this is really dodgy as non-square matrices don't have an inverse technically...
        result = rsp_inv_ndarray @ y_ndarray
        # print(result.shape)
        result_hist = cu.ndarray_to_th1(result.T, offset=0.5)

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

        rhs_hist = cu.ndarray_to_th1(rhs.T, offset=0.5)
        canv.SetRightMargin(0.12)
        rhs_hist.SetTitle("A^{T}V^{-1}_{yy}y;Generator bin;N")
        rhs_hist.Draw("HIST")

        canv.SaveAs(os.path.join(output_dir, "rhs.pdf"))

        proper_x_hist = cu.ndarray_to_th1(proper_x.T, offset=0.5)

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

    def calculate_chi2(self, one_hist, other_hist, cov_inv_matrix,
                       cov_matrix=None, detector_space=True, ignore_underflow_bins=True,
                       has_underflow=True, has_overflow=True, debugging_dir=None):
        """Calculate chi2 stats between one_hist and other_hist,
        using cov_inv_matrix: chi2 = (delta)^T * cov_inv_matrix * delta,
        where delta = one_hist-other_hist

        Can also do debugging plots, to plot individual parts of calculation

        Parameters
        ----------
        one_hist : ROOT.TH1, numpy.ndarray
            1D input hist
        other_hist : ROOT.TH1, numpy.ndarray
            Other 1D input hist
        cov_inv_matrix : ROOT.TH2, numpy.ndarray
            Inverse covariance matrix
        cov_matrix : None, optional
            Covariance matrix. Only used for plotting V^-1 * V
        detector_space : bool, optional
            True if working in detector binning, otherwise working in generator binning
        ignore_underflow_bins : bool, optional
            If True, remove pt underflow region from calculation (input hists & matrix)
        has_underflow : bool, optional
            If underflow bins are included in the input hists/cov matrix
            Only used for plotting purposes
        has_overflow : bool, optional
            If over bins are included in the input hists/cov matrix
            Only used for plotting purposes
        debugging_dir : None, optional str
            Output dir for debugging plot.
            If None or "", plots will not be made

        Returns
        -------
        float, int, float
            chi2, nbins, p-value
        """
        one_vec = one_hist
        if isinstance(one_hist, ROOT.TH1):
            one_vec, _ = cu.th1_to_ndarray(one_hist, False)
        other_vec = other_hist
        if isinstance(other_hist, ROOT.TH1):
            other_vec, _ = cu.th1_to_ndarray(other_hist, False)
        delta = one_vec - other_vec

        bs = 'detector' if detector_space else 'generator'
        first_signal_bin = 1
        if ignore_underflow_bins:
            # set to 0 all the pt underflow bins
            first_signal_bin = self.binning_handler.get_binning_scheme(bs).get_start_bin(is_signal_region=True)
            delta[0][:first_signal_bin-1] = 0. # subtract 1 as numpy indices start at 0, hists start at 1

        if isinstance(cov_inv_matrix, ROOT.TH2):
            v_inv, _ = cu.th2_to_ndarray(cov_inv_matrix)
        else:
            v_inv = cov_inv_matrix

        inter = v_inv.dot(delta.T)
        chi2 = delta.dot(inter)[0][0]
        ndof = len((one_vec[0][first_signal_bin-1:] != 0) | (other_vec[0][first_signal_bin-1:] != 0))  # only consider n bins where at least one has data - if both 0, don't count it
        p = 1-scipy.stats.chi2.cdf(chi2, int(ndof))

        # Print some debugging plots
        # --------------------------
        if debugging_dir:
            debug_plotter = MyUnfolderPlotter(self, is_data=False)

            # 1D inputs
            h_one = one_hist
            h_other = other_hist
            if isinstance(one_hist, np.ndarray):
                h_one = cu.ndarray_to_th1(one_hist, offset=0.5)
            if isinstance(other_hist, np.ndarray):
                h_other = cu.ndarray_to_th1(other_hist, offset=0.5)
            entries = [
                Contribution(h_one, label='one_hist', line_color=ROOT.kBlack),
                Contribution(h_other, label='other_hist', line_color=ROOT.kRed, line_style=2)
            ]
            plot = Plot(entries, what='hist', has_data=False,)
            plot.default_canvas_size = (800, 600)
            plot.plot("NOSTACK HIST")

            offset = 0
            if not has_underflow:
                first_signal_bin = self.binning_handler.get_binning_scheme(bs).get_start_bin(is_signal_region=True)
                offset = -first_signal_bin + 2
                # print("Setting offset to", offset)
            l, t = debug_plotter.draw_pt_binning_lines(plot,
                                                       which='reco' if detector_space else 'gen',
                                                       axis='x',
                                                       do_underflow=has_underflow,
                                                       do_overflow=has_overflow,
                                                       offset=offset)

            plot.save(os.path.join(debugging_dir, 'one_other_hists.pdf'))

            plot.set_logy(override_check=True)
            plot.save(os.path.join(debugging_dir, 'one_other_hists_logY.pdf'))

            # Delta, with missing bins if necessary
            delta_hist = cu.ndarray_to_th1(delta, offset=0.5)
            entries = [
                Contribution(delta_hist)
            ]
            plot = Plot(entries,
                        what='hist',
                        xtitle='%s bin' % ('Detector' if detector_space else 'Generator'),
                        ytitle='#Delta = one_hist - other_hist',
                        has_data=False,
                        )
            plot.default_canvas_size = (800, 600)
            plot.plot("NOSTACK HIST")
            l,t = debug_plotter.draw_pt_binning_lines(plot,
                                                      which='reco' if detector_space else 'gen',
                                                      axis='x',
                                                      do_underflow=has_underflow,
                                                      do_overflow=has_overflow,
                                                      offset=offset)
            plot.save(os.path.join(debugging_dir, 'delta.pdf'))

            plot.set_logy(override_check=True)
            plot.save(os.path.join(debugging_dir, 'delta_logY.pdf'))

            # Inverse Covariance matrix
            canv = ROOT.TCanvas("c", "Inverse covariance matrix V^{-1}", 800, 600)
            obj = cov_inv_matrix
            if not isinstance(cov_inv_matrix, ROOT.TH2):
                obj = cu.ndarray_to_th2(v_inv, offset=0.5)
            obj.SetTitle("Inverse covariance matrix V^{-1}")
            obj.Draw("COLZ")
            canv.SetLeftMargin(0.15)
            canv.SetRightMargin(0.18)

            l,t = debug_plotter.draw_pt_binning_lines(obj,
                                                      which='reco' if detector_space else 'gen',
                                                      axis='x',
                                                      do_underflow=has_underflow,
                                                      do_overflow=has_overflow,
                                                      do_labels_inside=False,
                                                      do_labels_outside=True,
                                                      offset=offset+1)
            l2,t2 = debug_plotter.draw_pt_binning_lines(obj,
                                                        which='reco' if detector_space else 'gen',
                                                        axis='y',
                                                        do_underflow=has_underflow,
                                                        do_overflow=has_overflow,
                                                        do_labels_inside=False,
                                                        do_labels_outside=True,
                                                        offset=offset+1)
            canv.SaveAs(os.path.join(debugging_dir, 'cov_inv_matrix_linZ.pdf'))

            canv.SetLogz(1)
            cov_min = obj.GetMinimum(1E-20) / 10
            obj.SetMinimum(cov_min)
            cov_max = obj.GetMaximum() * 5
            obj.SetMaximum(cov_max)
            canv.SaveAs(os.path.join(debugging_dir, 'cov_inv_matrix_logZ.pdf'))

            # Covariance matrix
            if cov_matrix is not None:
                canv = ROOT.TCanvas("cc", "Covariance matrix V", 800, 600)
                obj = cov_matrix
                if not isinstance(cov_matrix, ROOT.TH2):
                    obj = cu.ndarray_to_th2(cov_matrix, offset=0.5)
                obj.SetTitle("Covariance matrix V")
                obj.Draw("COLZ")
                canv.SetLeftMargin(0.15)
                canv.SetRightMargin(0.18)

                l,t = debug_plotter.draw_pt_binning_lines(obj,
                                                          which='reco' if detector_space else 'gen',
                                                          axis='x',
                                                          do_underflow=has_underflow,
                                                          do_overflow=has_overflow,
                                                          do_labels_inside=False,
                                                          do_labels_outside=True,
                                                          offset=offset+1)
                l2,t2 = debug_plotter.draw_pt_binning_lines(obj,
                                                            which='reco' if detector_space else 'gen',
                                                            axis='y',
                                                            do_underflow=has_underflow,
                                                            do_overflow=has_overflow,
                                                            do_labels_inside=False,
                                                            do_labels_outside=True,
                                                            offset=offset+1)
                canv.SaveAs(os.path.join(debugging_dir, 'cov_matrix_linZ.pdf'))

                canv.SetLogz(1)
                cov_min = obj.GetMinimum(1E-20) / 10
                obj.SetMinimum(cov_min)
                cov_max = obj.GetMaximum() * 5
                obj.SetMaximum(cov_max)
                canv.SaveAs(os.path.join(debugging_dir, 'cov_matrix_logZ.pdf'))

                # Plot product V V^-1 to test inverse quality
                v = cov_matrix
                if not isinstance(cov_matrix, np.ndarray):
                    v = cu.th2_to_ndarray(cov_matrix)[0]

                v_inv = cov_inv_matrix
                if not isinstance(cov_inv_matrix, np.ndarray):
                    v_inv = cu.th2_to_ndarray(cov_inv_matrix)[0]

                prod = v @ v_inv
                # print("product shape:", prod.shape)
                canv.Clear()
                prod_th2 = cu.ndarray_to_th2(prod, offset=0.5)
                prod_th2.SetTitle("V V^{-1}")
                prod_th2.Draw("COLZ")
                l,t = debug_plotter.draw_pt_binning_lines(prod_th2,
                                                          which='reco' if detector_space else 'gen',
                                                          axis='x',
                                                          do_underflow=has_underflow,
                                                          do_overflow=has_overflow,
                                                          do_labels_inside=False,
                                                          do_labels_outside=True,
                                                          offset=offset+1)
                l2,t2 = debug_plotter.draw_pt_binning_lines(prod_th2,
                                                            which='reco' if detector_space else 'gen',
                                                            axis='y',
                                                            do_underflow=has_underflow,
                                                            do_overflow=has_overflow,
                                                            do_labels_inside=False,
                                                            do_labels_outside=True,
                                                            offset=offset+1)
                canv.SetLogz(False)
                canv.SaveAs(os.path.join(debugging_dir, "cov_cov_inv_product_linZ.pdf"))
                canv.SetLogz(True)
                canv.SaveAs(os.path.join(debugging_dir, "cov_cov_inv_product_logZ.pdf"))

            # Components of delta * V_inv * delta before summing
            components = delta * inter.T
            components_hist = cu.ndarray_to_th1(components, offset=0.5)
            # y_max = components_hist.GetMaximum() * 1.2
            entries = [
                Contribution(components_hist)
            ]
            plot = Plot(entries,
                        what='hist',
                        xtitle='%s bin' % ('Detector' if detector_space else 'Generator'),
                        ytitle='Component of chi2 (#Delta V^{-1} #Delta)',
                        ylim=None,
                        has_data=False,
                        )
            plot.default_canvas_size = (800, 600)
            plot.plot("NOSTACK HIST")
            l,t = debug_plotter.draw_pt_binning_lines(plot,
                                                      which='reco' if detector_space else 'gen',
                                                      axis='x',
                                                      do_underflow=has_underflow,
                                                      do_overflow=has_overflow,
                                                      offset=offset)
            plot.save(os.path.join(debugging_dir, 'components.pdf'))

            # pt_bin_edges = self.pt_bin_edges_reco if detector_space else self.pt_bin_edges_gen

            # for ibin_pt, (pt_low, pt_high) in enumerate(zip(pt_bin_edges[:-1], pt_bin_edges[1:])):
            #     # plot component (delta * V_inv * delta) for this pt bin
            #     this_h = self.hist_bin_chopper.get_var_hist_pt_binned(components_hist, ibin_pt, binning_scheme='detector' if detector_space else 'generator')
            #     entries = [
            #         Contribution(this_h, label=None)
            #     ]
            #     plot = Plot(entries,
            #                 what='hist',
            #                 title='%g < p_{T} %g GeV' % (pt_low, pt_high),
            #                 xtitle='lambda variable',
            #                 ytitle='Component of chi2 (#Delta V^{-1} #Delta)',
            #                 ylim=(0, y_max),
            #                 has_data=False,
            #                 )
            #     plot.default_canvas_size = (800, 600)
            #     plot.plot("NOSTACK TEXT HIST")
            #     plot.save(os.path.join(debugging_dir, 'components_pt_bin_%d.pdf' % (ibin_pt)))

            #     # plot delta for this pt bin
            #     this_delta = self.hist_bin_chopper.get_var_hist_pt_binned(delta_hist, ibin_pt, binning_scheme='detector' if detector_space else 'generator')
            #     entries = [
            #         Contribution(this_delta, label=None)
            #     ]
            #     plot = Plot(entries,
            #                 what='hist',
            #                 title='%g < p_{T} %g GeV' % (pt_low, pt_high),
            #                 xtitle='lambda variable',
            #                 ytitle='#Delta',
            #                 # ylim=(0, y_max),
            #                 has_data=False,
            #                 )
            #     plot.default_canvas_size = (800, 600)
            #     plot.plot("NOSTACK HIST TEXT")
            #     plot.save(os.path.join(debugging_dir, 'delta_pt_bin_%d.pdf' % (ibin_pt)))

            #     # plot cov matrix for this bin
            #     binning = self.detector_binning.FindNode("detectordistribution") if detector_space else self.generator_binning.FindNode("generatordistribution")
            #     var_bins = np.array(binning.GetDistributionBinning(0))
            #     pt_bins = np.array(binning.GetDistributionBinning(1))
            #     # -1 since ndarray is 0 index, th2 are 1-indexed
            #     start = binning.GetGlobalBinNumber(var_bins[0]+0.001, pt_bins[ibin_pt]+0.1) - 1
            #     end = binning.GetGlobalBinNumber(var_bins[-2]+0.001, pt_bins[ibin_pt]+0.1) - 1
            #     this_cov = cu.ndarray_to_th2(v_inv[start:end+1,start:end+1])
            #     canv = ROOT.TCanvas(cu.get_unique_str(), "Inverse covariance matrix V^{-1}", 800, 600)
            #     this_cov.SetTitle("Inverse covariance matrix V^{-1} for %g < p_{T} < %g GeV" % (pt_low, pt_high))
            #     this_cov.Draw("COLZ TEXT")
            #     canv.SetLeftMargin(0.15)
            #     canv.SetRightMargin(0.18)
            #     # canv.SetLogz(1)
            #     # this_cov.SetMinimum(cov_min)
            #     # this_cov.SetMaximum(cov_max)
            #     canv.SaveAs(os.path.join(debugging_dir, 'cov_inv_matrix_pt_bin_%d.pdf' % (ibin_pt)))

            # histogram the non-zero components, but with their sign
            # print('components:', components)
            # new_components = np.sqrt(components[0]) * np.sign(delta[0])
            # print('new_components:', new_components)
            # x_low_lim = new_components.min()
            # x_high_lim = new_components.max()
            # print(x_low_lim, x_high_lim)
            # x_range = x_high_lim - x_low_lim
            # x_high_lim += (0.05 * x_range)
            # x_low_lim -= (0.05 * x_range)
            # h = ROOT.TH1D("h_components", "", 25 if detector_space else 10, x_low_lim, x_high_lim)
            # for i, x in enumerate(components[0]):
            #     # if x > 0:
            #     # print(x)
            #     h.Fill(math.sqrt(x) * np.sign(delta[0][i]))
            hist, bins = np.histogram(components, bins=10)
            nbins = len(bins)-1
            h = ROOT.TH1D("h_components", "", nbins, bins)
            for i, x in enumerate(hist, 1):
                h.SetBinContent(i, x)
            h.Fit('gaus', 'q')
            entries = [
                Contribution(h)
            ]
            plot = Plot(entries,
                        what='hist',
                        # xtitle='sign(#Delta) #times #sqrt{#Delta V^{-1} #Delta}',
                        xtitle='#Delta #times V^{-1} #times #Delta',
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
        # num_vars = len(jackknife_variations)
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

            bins = self.binning_handler.get_variable_bins(pt=self.pt_bin_edges_gen[ibin_pt], binning_scheme='generator')
            nbins = len(bins) - 1
            # FIXME which binning to use? index or physical?
            this_bin_unfolded_rsp_ematrix = ROOT.TH2D("ematrix_rsp_bin_%d_%s" % (ibin_pt, cu.get_unique_str()), "", nbins, 0, nbins, nbins, 0, nbins)
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
                del self.hist_bin_chopper._cache[k]

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

            bins = self.binning_handler.get_variable_bins(pt=self.pt_bin_edges_gen[ibin_pt], binning_scheme='generator')
            nbins = len(bins) - 1
            # FIXME which binning to use? index or physical?
            this_bin_unfolded_stat_ematrix = ROOT.TH2D("ematrix_stat_bin_%d_%s" % (ibin_pt, cu.get_unique_str()), "", nbins, 0, nbins, nbins, 0, nbins)
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
                del self.hist_bin_chopper._cache[k]

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
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
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
            # nominal = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err',
            #                                                                 ibin_pt,
            #                                                                 binning_scheme='generator')
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

            variations_envelope = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin_pt, binning_scheme='generator').Clone("pdf_%d" % ibin_pt)

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
    def convert_error_bars_to_error_shift(h):
        """Create histogram with bin contents equal to error bar on h,
        and 0 error bars"""
        h_new = h.Clone(cu.get_unique_str())
        for i in range(1, h.GetNbinsX()+1):
            h_new.SetBinContent(i, h.GetBinError(i))
            h_new.SetBinError(i, 0)
        return h_new

    @staticmethod
    def convert_error_shift_to_error_bars(h_nominal, h_shift):
        """Create histogram with bin contents from h_nominal,
        and error bars from bin values of h_shift"""
        h = h_nominal.Clone(cu.get_unique_str())
        for i in range(1, h_nominal.GetNbinsX()+1):
            h.SetBinError(i, h_shift.GetBinContent(i))
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
                cu.remove_th1_errors(syst_hist)
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
        pt_bins = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=True)
        for ibin_pt, pt in enumerate(pt_bins[:-1]):
            first_bin = ibin_pt == 0
            hbc_args = dict(ind=ibin_pt, binning_scheme='generator', is_signal_region=True)

            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_rsp_err', **hbc_args)

            unfolded_hist_bin_abs = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', **hbc_args)

            norm = 1
            if unfolded_hist_bin_abs.Integral("width") == 0:
                warnings.warn("setup_normalised_results_per_pt_bin: unfolded_hist_bin_abs.Integral = 0 in bin {0} for pt {1}".format(ibin_pt, pt))
            if unfolded_hist_bin_stat_errors.Integral("width") == 0:
                warnings.warn("setup_normalised_results_per_pt_bin: unfolded_hist_bin_stat_errors.Integral = 0 in bin {0} for pt {1}".format(ibin_pt, pt))
            else:
                norm = unfolded_hist_bin_abs.Integral("width") / unfolded_hist_bin_stat_errors.Integral("width")

            # create stat & rsp err covariance matrices for this pt bin,
            # if they haven't been calculated by jackknife methods,
            # scaling by overall normalisation and bin widths
            var_bins = self.binning_handler.get_variable_bins(pt, binning_scheme='generator')
            # FIXME what to do if non-sequential bin numbers?!
            # the 0.001 is to ensure we're def inside this bin
            start_bin = self.binning_handler.physical_bin_to_global_bin(var=var_bins[0] + 1E-6, pt=pt + 1E-6, binning_scheme='generator')
            end_bin = self.binning_handler.physical_bin_to_global_bin(var=var_bins[-2] + 1E-6, pt=pt + 1E-6, binning_scheme='generator')  # -2 since the last one is the upper edge of the last bin
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
                # print(i, "old err:", h_total.GetBinError(i), "new err:", math.sqrt(err2))
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


    # METHODS FOR JACOBIAN ERRORS
    # --------------------------------------------------------------------------
    def create_normalisation_jacobian_np(self):
        """Create Jacobian matrix for normalisation

        For transformation of variable x -> y = f(x), then J(i,j) = df_i / dx_j.

        For 1D normalisation, where y_i = x_i / sum(x_j), with N = sum(x_j)
        that means
        df_i / dx_j = (N * delta(i,j) - x_i) / N^2  (delta is kroenecker delta)

        For 2D normalisation, it's a block matrix per pT bin.
        """
        if getattr(self, "jacobian", None) is None:
            h = self.get_output()
            nbins = h.GetNbinsX()

            # to keep track of the possible under/overflow bins
            nbins_uflow = self.generator_distribution_underflow.GetDistributionNumberOfBins()
            nbins_signal = self.generator_distribution.GetDistributionNumberOfBins()
            assert(nbins == nbins_uflow+nbins_signal)

            J = np.zeros(shape=(nbins, nbins))

            pt_axis_ind = 2 # 0 is lambda axis, 1 is pt?
            nbins_pt_uflow = self.nbins_pt_underflow_gen + int(self.generator_distribution_underflow.HasUnderflow(pt_axis_ind)) + int(self.generator_distribution_underflow.HasOverflow(pt_axis_ind))
            nbins_pt_signal = self.nbins_pt_gen + int(self.generator_distribution.HasUnderflow(pt_axis_ind)) + int(self.generator_distribution.HasOverflow(pt_axis_ind))
            nbins_pt = nbins_pt_uflow+nbins_pt_signal

            variable_axis_ind = 0 # 0 is lambda axis, 1 is pt?
            nbins_variable = self.nbins_variable_gen + int(self.generator_distribution.HasUnderflow(variable_axis_ind)) + int(self.generator_distribution.HasOverflow(variable_axis_ind))

            # print(nbins_pt_uflow, nbins_pt_signal, nbins_variable)

            for pt_ind in range(0, nbins_pt):
                # normalising by pt bin, so get the integral for this pt bin
                # +1 cos it's ROOT
                start_ind = (pt_ind*nbins_variable)+1
                end_ind = start_ind+nbins_variable-1
                N = h.Integral(start_ind, end_ind)
                # print("N", N, "for", start_ind, end_ind)
                # for i in range(start_ind, end_ind+1):
                #     print("   i:", i, h.GetBinContent(i))
                # for each block, iterate over all bins
                start_x = start_y = nbins_variable * pt_ind

                for var_ind_x in range(0, nbins_variable):
                    for var_ind_y in range(0, nbins_variable):
                        bin_ind = (nbins_variable * pt_ind) + var_ind_y
                        # print(bin_ind)
                        bin_content = h.GetBinContent(bin_ind + 1)
                        if bin_content < 0:
                            # protection incase of -ve bins
                            print("Warning: jacobian had -ve bin content (%f) for index %d" % (bin_content, bin_ind+1))
                            bin_content = 0
                        if var_ind_x == var_ind_y:
                            val = (N - bin_content) / N**2
                        else:
                            val = -bin_content / N**2
                        J[start_y + var_ind_y][start_x + var_ind_x] = val
                        # print("J[%d][%d] =" % (start_y+var_ind_y, start_x+var_ind_x), val)

            self.jacobian = J
        return self.jacobian

    def get_jacobian_th2(self):
        if getattr(self, "jacobian_th2", None) is None:
            self.jacobian_th2 = cu.ndarray_to_th2(self.create_normalisation_jacobian_np(), offset=-0.5)
        return self.jacobian_th2

    def get_jacobian_np(self):
        if getattr(self, "jacobian", None) is None:
            self.jacobian_ = cu.th2_to_ndarray(self.get_jacobian_th2())
        return self.jacobian

    def create_absolute_scale_uncertainty(self, scale_systs):
        """Create scale shift, shifted results, & cov matrix, on absolute result"""
        # TODO: better to not store indiv hists in class?
        # if getattr(self, 'scale_shift', None) is None:
        if not self.has_exp_syst("Scale"):
            print("Setting up absolute scale ucnert")
            nominal = self.get_output()
            shift = nominal.Clone("abs_scale_shift")
            shifted_up = nominal.Clone("abs_scale_shifted_up")
            shifted_down = nominal.Clone("abs_scale_shifted_down")
            for ix in range(1, shift.GetNbinsX()+1):
                # for each bin in the output distribution,
                # get the largest deviation across all scale variations
                nom = nominal.GetBinContent(ix)
                this_shift = max([abs(syst['unfolder'].get_output().GetBinContent(ix) - nom)
                                  for syst in scale_systs])
                shift.SetBinContent(ix, this_shift)
                shift.SetBinError(ix, 0)
                shifted_up.SetBinContent(ix, nom+this_shift)
                shifted_up.SetBinError(ix, 0)
                shifted_down.SetBinContent(ix, nom-this_shift)
                shifted_down.SetBinError(ix, 0)
            self.scale_shift = shift
            self.scale_shifted_up = shifted_up
            self.scale_shifted_down = shifted_down

            scale_syst = ExpSystematic(label="Scale")
            scale_syst.syst_shift = shift
            scale_syst.syst_shifted = shifted_up

            # now construct error matrix usinv x * X^T
            shift_np, _ = cu.th1_to_ndarray(shift)
            cov_np = shift_np.T @ shift_np
            assert(cov_np.shape != (1, 1))  # get the dimensions right
            cov_hist = cu.ndarray_to_th2(cov_np, offset=-0.5)
            scale_syst.syst_ematrix = cov_hist
            self.exp_systs.append(scale_syst)
            self.ematrix_scale = cov_hist

        return self.scale_shift, self.scale_shifted_up, self.scale_shifted_down

    def create_absolute_pdf_uncertianty(self, pdf_systs):
        """Create pdf shift, shifted results,  & cov matrix on absolute result"""
        # TODO: better to not store indiv hists in class?

        # if getattr(self, 'pdf_shift', None) is None:
        if not self.has_exp_syst("PDF"):
            print("Setting up absolute PDF ucnert")
            nominal = self.get_output()
            shift = nominal.Clone("abs_pdf_shift")
            shifted_up = nominal.Clone("abs_pdf_shifted_up")
            shifted_down = nominal.Clone("abs_pdf_shifted_down")
            for ix in range(1, shift.GetNbinsX()+1):
                nom = nominal.GetBinContent(ix)
                # error is the RMS of the variations
                # np.std does sqrt((abs(x - x.mean())**2) / (len(x) - ddof)),
                # and the PDF4LHC recommendation is N-1 in the denominator
                rms = np.std([syst['unfolder'].get_output().GetBinContent(ix) for syst in pdf_systs], ddof=1)
                shift.SetBinContent(ix, rms)
                shift.SetBinError(ix, 0)
                shifted_up.SetBinContent(ix, nom+rms)
                shifted_up.SetBinError(ix, 0)
                shifted_down.SetBinContent(ix, nom-rms)
                shifted_down.SetBinError(ix, 0)
            self.pdf_shift = shift
            self.pdf_shifted_up = shifted_up
            self.pdf_shifted_down = shifted_down

            pdf_syst = ExpSystematic(label="PDF")
            pdf_syst.syst_shift = shift
            pdf_syst.syst_shifted = shifted_up

            # now construct error matrix usinv x * X^T
            shift_np, _ = cu.th1_to_ndarray(shift)
            cov_np = shift_np.T @ shift_np
            assert(cov_np.shape != (1, 1))  # get the dimensions right
            cov_hist = cu.ndarray_to_th2(cov_np, offset=-0.5)
            pdf_syst.syst_ematrix = cov_hist
            self.exp_systs.append(pdf_syst)
            self.ematrix_pdf = cov_hist

        # TODO use ExpSystematic obj
        return self.pdf_shift, self.pdf_shifted_up, self.pdf_shifted_down

    def create_absolute_total_uncertainty(self):
        """Get total absolute error matrix from:

        - input stats
        - response matrix stats
        - background stats
        - experimental systs
        (above from GetEmatrixTotal())
        - combined scale syst
        - combined pdf syst

        then create the corresponding shift & shifted distributions
        """
        # stat + rsp stat + exp systs
        print("create_absolute_total_uncertainty()")
        ematrix_total = self.get_ematrix_tunfold_total().Clone("ematrix_total")
        for n in ["Scale", "PDF"]:
            # scale systs if exists
            if self.has_exp_syst(n):
                ematrix_total.Add(self.get_exp_syst(n).syst_ematrix)
        # store
        total_syst = ExpSystematic(label="Total")
        total_syst.syst_ematrix = ematrix_total
        total_syst.syst_shift = self.make_hist_from_diagonal_errors(ematrix_total, do_sqrt=True, set_errors=False, offset=0.5)
        # get_syst_shifted_hist() will handle syst_shifted
        self.exp_systs.append(total_syst)

    def normalise_ematrix(self, ematrix):
        """Generic method to normalise any aboslute cov matrix using Jacobian

        Can handle numpy array, or TH2
        """
        J = self.create_normalisation_jacobian_np()
        if isinstance(ematrix, ROOT.TH2):
            cov_abs, _ = cu.th2_to_ndarray(ematrix)
        else:
            cov_abs = ematrix
        cov_norm = J @ cov_abs @ J.T
        ematrix_norm = cu.ndarray_to_th2(cov_norm, offset=-0.5)
        return ematrix_norm

    def normalise_all_systs(self):
        """Create normalised versions of all ExpSystematic stored in this instance"""
        print("normalise_all_systs()")
        new_systs = []
        for exp_syst in self.exp_systs:
            norm_exp_syst = ExpSystematic(label=exp_syst.label + "_Norm")
            # normalise the ematrix - everything comes from that
            norm_exp_syst.syst_ematrix = self.normalise_ematrix(exp_syst.syst_ematrix)
            norm_exp_syst.syst_shift = self.make_hist_from_diagonal_errors(norm_exp_syst.syst_ematrix, do_sqrt=True, set_errors=False, offset=0.5)
            # the shifted hists and error bar hists don't make sense here
            # - the normalised errors only apply to normalised bins
            new_systs.append(norm_exp_syst)
            self.exp_systs_normed.append(norm_exp_syst)
        self.exp_systs.extend(new_systs)
        print([e.label for e in self.exp_systs])

    def setup_absolute_results(self):
        """Load up HistBinChopper to with absolute results"""
        print("setup_absolute_results()")
        self.hist_bin_chopper.add_obj(self.no_uncert_name, self.get_unfolded_with_no_errors())

        # stat unc
        self.hist_bin_chopper.add_obj(self.stat_uncert_name, self.get_unfolded_with_ematrix_stat())

        # rsp unc
        self.hist_bin_chopper.add_obj(self.rsp_uncert_name, self.get_unfolded_with_ematrix_response())

        # all the exp syst ones, including scale, pdf, total
        for exp_syst in self.exp_systs:
            label = exp_syst.label
            print("... exp_syst:", label)
            if "_Norm" in label:
                continue
            # with shift as error bar
            self.hist_bin_chopper.add_obj(exp_syst.syst_error_bar_label, self.get_syst_error_hist(label))
            # shifted
            self.hist_bin_chopper.add_obj(exp_syst.syst_shifted_label, self.get_syst_shifted_hist(label, self.get_unfolded_with_no_errors()))
            if label == "Total":
                print("get_pt_bin_normed_div_bin_width(%s):" % exp_syst.syst_error_bar_label)
                h = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(exp_syst.syst_error_bar_label, ind=0, binning_scheme='generator')
                # cu.print_th1_bins(h)
                h_abs = self.hist_bin_chopper.get_pt_bin(exp_syst.syst_error_bar_label, ind=0, binning_scheme='generator')
                # cu.print_th1_bins(h_abs)
                N = h_abs.Integral()
                # print("N:", N)
                # print("Jacobian diagonals vs 1/N:")
                # for ix in range(1, h_abs.GetNbinsX()+1):
                #     print(ix, (N-h_abs.GetBinContent(ix))/(N**2), 1./N)
                # FIXME finish this

    def setup_normalised_results(self):
        # For each pt bin, we get the nominal result
        # Then we add / set error bar using the already-calculated normalised shifts
        print("setup_normalised_results()")

        # need to create normalised stat and rsp uncert, and add to HBC
        ematrix_stat_norm = self.normalise_ematrix(self.get_ematrix_stat())
        hist_stat_norm_shift = self.make_hist_from_diagonal_errors(ematrix_stat_norm, do_sqrt=True, set_errors=False)
        hist_stat_norm_err = self.make_hist_from_diagonal_errors(ematrix_stat_norm, do_sqrt=True, set_errors=True)
        self.hist_bin_chopper.add_obj('norm_stat_shift', hist_stat_norm_shift) # already normed!!
        self.hist_bin_chopper.add_obj('norm_stat_err', hist_stat_norm_err)

        ematrix_rsp_norm = self.normalise_ematrix(self.get_ematrix_stat_response())
        hist_rsp_norm_shift = self.make_hist_from_diagonal_errors(ematrix_rsp_norm, do_sqrt=True, set_errors=False)
        hist_rsp_norm_err = self.make_hist_from_diagonal_errors(ematrix_rsp_norm, do_sqrt=True, set_errors=True)
        self.hist_bin_chopper.add_obj('norm_rsp_shift', hist_rsp_norm_shift)
        self.hist_bin_chopper.add_obj('norm_rsp_err', hist_rsp_norm_err)

        # Add objects, will replace cache ourselves
        for exp_syst in self.exp_systs_normed:
            label = exp_syst.label
            unnorm = lambda x : x.replace("_Norm", "")
            unnorm = lambda x : x  # whyy
            print("setup_normalised_results(): Adding", unnorm(exp_syst.syst_error_bar_label))
            self.hist_bin_chopper.add_obj(unnorm(exp_syst.syst_error_bar_label), self.get_syst_error_hist(label))
            print("setup_normalised_results(): Adding", unnorm(exp_syst.syst_shifted_label))
            self.hist_bin_chopper.add_obj(unnorm(exp_syst.syst_shifted_label), self.get_syst_shifted_hist(label, self.get_unfolded_with_no_errors()))
            # add shift?

        for ibin in range(0, self.nbins_pt_gen):
            # iterate over each generator pt bin, and for each create normalised
            # hists and fill them
            # Also store the error matrix for each pt bin

            hbc_args = dict(ind=ibin, binning_scheme='generator')

            # functions to create key for storing the normed results
            key_gen = partial(self.hist_bin_chopper._generate_key,
                              ind=ibin,
                              axis='pt',
                              do_norm=True,
                              do_div_bin_width=False,
                              binning_scheme='generator')
            key_gen_div_width = partial(self.hist_bin_chopper._generate_key,
                                        ind=ibin,
                                        axis='pt',
                                        do_norm=True,
                                        do_div_bin_width=True,
                                        binning_scheme='generator')

            # get normalised unfolded hist for this bin
            nominal = self.hist_bin_chopper.get_pt_bin_normed(self.no_uncert_name, **hbc_args)

            # unfolded bin contents + stat unc error bars
            stat_err = self.hist_bin_chopper.get_pt_bin("norm_stat_err", **hbc_args)  # NOT get_normed - already normed!
            stat_err_hist = nominal.Clone(cu.get_unique_str())
            self.update_hist_bin_error(h_to_be_updated=stat_err_hist, h_orig=stat_err)
            self.hist_bin_chopper._cache[key_gen(self.stat_uncert_name)] = stat_err_hist
            # same but div bin width
            stat_err_hist_div_bin_width = qgp.hist_divide_bin_width(stat_err_hist)
            print("Updating", key_gen_div_width(self.stat_uncert_name))
            self.hist_bin_chopper._cache[key_gen_div_width(self.stat_uncert_name)] = stat_err_hist_div_bin_width
            # cu.print_th1_bins(stat_err_hist_div_bin_width)

            # unfolded bin contents + rsp unc error bars
            rsp_err = self.hist_bin_chopper.get_pt_bin("norm_rsp_err", **hbc_args)  # NOT get_normed - already normed!
            rsp_err_hist = nominal.Clone(cu.get_unique_str())
            self.update_hist_bin_error(h_to_be_updated=rsp_err_hist, h_orig=rsp_err)
            self.hist_bin_chopper._cache[key_gen(self.rsp_uncert_name)] = rsp_err_hist
            # same but div bin width
            rsp_err_hist_div_bin_width = qgp.hist_divide_bin_width(rsp_err_hist)
            print("Updating", key_gen_div_width(self.rsp_uncert_name))
            self.hist_bin_chopper._cache[key_gen_div_width(self.rsp_uncert_name)] = rsp_err_hist_div_bin_width

            # all the exp syst ones, including scale, pdf, total
            for exp_syst in self.exp_systs:
                label = exp_syst.label
                if "_Norm" not in label:
                    continue

                # only deal with Normed systs, but for the key generation
                # we need the un-normed name:
                unnorm = lambda x : x.replace("_Norm", "")

                # combine normalised unfolded hist contents with already normalised errors
                this_err_hist = nominal.Clone(cu.get_unique_str())
                this_err = self.hist_bin_chopper.get_pt_bin(exp_syst.syst_error_bar_label, **hbc_args)
                self.update_hist_bin_error(h_to_be_updated=this_err_hist, h_orig=this_err)
                self.hist_bin_chopper._cache[key_gen(unnorm(exp_syst.syst_error_bar_label))] = this_err_hist

                # same but divide by bin width
                this_err_hist_div_bin_width = qgp.hist_divide_bin_width(this_err_hist)
                self.hist_bin_chopper._cache[key_gen_div_width(unnorm(exp_syst.syst_error_bar_label))] = this_err_hist_div_bin_width

                # if "Total" in label:
                    # print("updating", key_gen_div_width(unnorm(exp_syst.syst_error_bar_label)))
                    # cu.print_th1_bins(this_err_hist_div_bin_width)

                # do shifted version
                this_shifted_hist = nominal.Clone(cu.get_unique_str())
                this_shifted = self.hist_bin_chopper.get_pt_bin(exp_syst.syst_error_bar_label, **hbc_args)
                self.update_hist_bin_error(h_to_be_updated=this_shifted_hist, h_orig=this_shifted)
                self.hist_bin_chopper._cache[key_gen(unnorm(exp_syst.syst_shifted_label))] = this_shifted_hist

                # same but divide by bin width
                this_shifted_hist_div_bin_width = qgp.hist_divide_bin_width(this_shifted_hist)
                self.hist_bin_chopper._cache[key_gen_div_width(unnorm(exp_syst.syst_shifted_label))] = this_shifted_hist_div_bin_width


    # METHODS FOR ABSOLUTE PER PT BIN RESULTS
    # --------------------------------------------------------------------------

    def create_scale_syst_uncertainty_per_pt_bin(self, scale_systs):
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
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
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
                self.hist_bin_chopper.get_pt_bin_div_bin_width(syst['hbc_key_unfolded'], **hbc_args)
                for syst in scale_systs
            ]

            # Calculate envelope error bar from max variation in each bin
            nominal = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', **hbc_args)
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
                                                      do_norm=False,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = variations_envelope

    def create_scale_syst_ematrices_per_pt_bin(self):
        """Create ematrix corresponding to scale uncertainty for each pt bin

        Required create_scale_syst_uncertainty_per_pt_bin() to be run first

        Calculated as x.x^T, where x is the difference between the scale-uncert
        hist, and the nominal hist

        Stored in self.hist_bin_chopper with key self.scale_uncert_ematrix_name
        """
        # add dummy object to HistBinChopper.objects so check doesn't fail
        self.hist_bin_chopper.add_obj(self.scale_uncert_ematrix_name, self.get_unfolded_with_ematrix_stat())

        for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
            # nominal = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err',
            #                                                                  ibin_pt,
            #                                                                  binning_scheme='generator')
            scale_hist = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.scale_uncert_name,
                                                                               ibin_pt,
                                                                               binning_scheme='generator')
            scale_shift = self.convert_error_bars_to_error_shift(scale_hist)
            scale_ematrix = cu.shift_to_covariance(scale_shift)
            key = self.hist_bin_chopper._generate_key(self.scale_uncert_ematrix_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=False,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = scale_ematrix

    def create_pdf_syst_uncertainty_per_pt_bin(self, pdf_systs):
        """Create PDF uncertainty from unfolded PDF variations

        This is done by taking in all the unfolded pdf systematics results,
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
            variations = [syst['unfolder'].hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', ibin_pt, binning_scheme='generator')
                          for syst in pdf_systs]

            variations_envelope = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', ibin_pt, binning_scheme='generator').Clone("pdf_%d" % ibin_pt)

            for ix in range(1, variations_envelope.GetNbinsX()+1):
                # np.std does sqrt((abs(x - x.mean())**2) / (len(x) - ddof)),
                # and the PDF4LHC recommendation is N-1 in the denominator
                rms_ratio = np.std([v.GetBinContent(ix) for v in variations], ddof=1)
                variations_envelope.SetBinError(ix, rms_ratio)

            # Store in hist_bin_chopper for later
            key = self.hist_bin_chopper._generate_key(self.pdf_uncert_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=False,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = variations_envelope

    def create_pdf_syst_ematrices_per_pt_bin(self):
        """Create ematrix corresponding to pdf uncertainty for each pt bin

        Requires create_pdf_syst_uncertainty_per_pt_bin() to be run first

        Calculated as x.x^T, where x is the difference between the pdf-uncert
        normalised hist, and the nominal normalised hist

        Stored in self.hist_bin_chopper with key self.pdf_uncert_ematrix_name
        """
        # add dummy object to HistBinChopper.objects so check doesn't fail
        self.hist_bin_chopper.add_obj(self.pdf_uncert_ematrix_name, self.get_unfolded_with_ematrix_stat())

        for ibin_pt in range(len(self.pt_bin_edges_gen[:-1])):
            # nominal = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err',
            #                                                                  ibin_pt,
            #                                                                  binning_scheme='generator')
            pdf_hist = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.pdf_uncert_name,
                                                                               ibin_pt,
                                                                               binning_scheme='generator')
            pdf_shift = self.convert_error_bars_to_error_shift(pdf_hist)
            pdf_ematrix = cu.shift_to_covariance(pdf_shift)
            key = self.hist_bin_chopper._generate_key(self.pdf_uncert_ematrix_name,
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=False,
                                                      do_div_bin_width=True,
                                                      binning_scheme='generator')
            self.hist_bin_chopper._cache[key] = pdf_ematrix

    def setup_absolute_results_per_pt_bin(self):
        """Setup final absolute results per pt bin with all uncertainties.

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

            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', **hbc_args)
            unfolded_hist_bin_rsp_errors = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_rsp_err', **hbc_args)

            # unfolded_hist_bin_abs = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', **hbc_args)
            # norm = unfolded_hist_bin_abs.Integral("width") / unfolded_hist_bin_stat_errors.Integral("width")

            # # create stat & rsp err covariance matrices for this pt bin,
            # # if they haven't been calculated by jackknife methods,
            # # scaling by overall normalisation and bin widths
            # binning = self.generator_binning.FindNode("generatordistribution")
            # var_bins = self.variable_bin_edges_gen
            # # FIXME what to do if non-sequential bin numbers?!
            # # the 0.001 is to ensure we're def inside this bin
            # start_bin = binning.GetGlobalBinNumber(var_bins[0]+0.001, pt+0.001)
            # end_bin = binning.GetGlobalBinNumber(var_bins[-2]+0.001, pt+0.001)  # -2 since the last one is the upper edge of the last bin
            # stat_key = self.hist_bin_chopper._generate_key(self.stat_ematrix_name,
            #                                                ind=ibin_pt,
            #                                                axis='pt',
            #                                                do_norm=False,
            #                                                do_div_bin_width=True,
            #                                                binning_scheme='generator')
            # if stat_key not in self.hist_bin_chopper._cache:
            #     # Get the stat error matrix from TUnfold, then select the sub-matrix
            #     # for this pt bin, then scale by normalisation and bin widths
            #     stat_ematrix = self.get_sub_th2(self.get_ematrix_stat(), start_bin, end_bin)
            #     stat_ematrix.Scale(1./(norm*norm))
            #     self.scale_th2_bin_widths(stat_ematrix, var_bins)
            #     self.hist_bin_chopper._cache[stat_key] = stat_ematrix

            # rsp_key = self.hist_bin_chopper._generate_key(self.rsp_ematrix_name,
            #                                               ind=ibin_pt,
            #                                               axis='pt',
            #                                               do_norm=False,
            #                                               do_div_bin_width=True,
            #                                               binning_scheme='generator')
            # if rsp_key not in self.hist_bin_chopper._cache:
            #     # if it has been setup already, it was from jackknife
            #     # otherwise we use the one from TUnfold
            #     rsp_ematrix = self.get_sub_th2(self.get_ematrix_stat_response(), start_bin, end_bin)
            #     rsp_ematrix.Scale(1./(norm*norm))
            #     self.scale_th2_bin_widths(rsp_ematrix, var_bins)
            #     self.hist_bin_chopper._cache[rsp_key] = rsp_ematrix

            # # Calculate total ematrix
            # total_ematrix = self.hist_bin_chopper._cache[stat_key].Clone()
            # total_ematrix.Add(self.hist_bin_chopper._cache[rsp_key])

            # for exp_syst in self.exp_systs:
            #     if first_bin:
            #         print("Adding", exp_syst.label, "ematrix to total normalised ematrix...")
            #     total_ematrix.Add(self.hist_bin_chopper.get_pt_bin_div_bin_width(exp_syst.syst_ematrix_label, **hbc_args))

            # if self.scale_uncert_ematrix_name in self.hist_bin_chopper.objects:
            #     if first_bin:
            #         print("Adding scale ematrix to total normalised ematrix")
            #     total_ematrix.Add(self.hist_bin_chopper.get_pt_bin_div_bin_width(self.scale_uncert_ematrix_name, **hbc_args))

            # if self.pdf_uncert_ematrix_name in self.hist_bin_chopper.objects:
            #     if first_bin:
            #         print("Adding pdf ematrix to total normalised ematrix")
            #     total_ematrix.Add(self.hist_bin_chopper.get_pt_bin_div_bin_width(self.pdf_uncert_ematrix_name, **hbc_args))

            # key = self.hist_bin_chopper._generate_key(self.total_ematrix_name,
            #                                           ind=ibin_pt,
            #                                           axis='pt',
            #                                           do_norm=False,
            #                                           do_div_bin_width=True,
            #                                           binning_scheme='generator')
            # self.hist_bin_chopper._cache[key] = total_ematrix

            error_bar_hists = [unfolded_hist_bin_stat_errors, unfolded_hist_bin_rsp_errors]

            # convert all shifts to error bars
            for exp_syst in self.exp_systs:
                if first_bin:
                    print("Adding", exp_syst.label, "uncertainty to normalised result...")
                # Here we access the things we just manually put in the cache - must match up with key!
                # Don't worry about it being normed etc - that is just so keys agree, and it matches
                # up with the nominal result (which we do want normed_div_bin_width)
                # Create shift due to this syst
                syst_shift = self.hist_bin_chopper.get_pt_bin_div_bin_width(exp_syst.syst_shifted_label, **hbc_args).Clone()
                syst_shift.Add(unfolded_hist_bin_stat_errors, -1)
                error_bar_hists.append(self.convert_error_shift_to_error_bars(unfolded_hist_bin_stat_errors, syst_shift))
            # Add in scale syst
            if self.scale_uncert_name in self.hist_bin_chopper.objects:
                if first_bin:
                    print("Adding scale uncertainty to normalised result...")
                error_bar_hists.append(self.hist_bin_chopper.get_pt_bin_div_bin_width(self.scale_uncert_name, **hbc_args))

            # Add in PDF syst
            if self.pdf_uncert_name in self.hist_bin_chopper.objects:
                if first_bin:
                    print("Adding PDF uncertainty to normalised result...")
                error_bar_hists.append(self.hist_bin_chopper.get_pt_bin_div_bin_width(self.pdf_uncert_name, **hbc_args))

            # Get normalised hist with nominal unfolded value, and change error bars
            # to be quadrature sum of those we want (stat+rsp+systs)
            h_total = self.hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', **hbc_args)
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
                                                      ind=ibin_pt,
                                                      axis='pt',
                                                      do_norm=False,
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
            "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
            "colour": ROOT.kAzure,
            "unfolder": MyUnfolder,
        }
        """
        if isinstance(self.binning_handler.get_binning_scheme('generator'), PtVarPerPtBinning):
            warnings.warn(cu.pcolors.RED + "Cannot do create_scale_syst_uncertainty_per_lambda_bin() with PtVarPerPtBinning" + cu.pcolors.ENDC)
            return

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
        if isinstance(self.binning_handler.get_binning_scheme('generator'), PtVarPerPtBinning):
            warnings.warn(cu.pcolors.RED + "Cannot do create_pdf_syst_uncertainty_per_lambda_bin() with PtVarPerPtBinning" + cu.pcolors.ENDC)
            return

        # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
        self.hist_bin_chopper.add_obj(self.pdf_uncert_name, self.get_unfolded_with_ematrix_stat())

        for ibin_var in range(len(self.variable_bin_edges_gen[:-1])):
            hbc_args = dict(ind=ibin_var, binning_scheme="generator")
            # Calculate error by using RMS of variations in each bin of the histogram
            variations = [syst['unfolder'].hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded', **hbc_args)
                          for syst in pdf_systs]

            variations_envelope = self.hist_bin_chopper.get_lambda_bin_div_bin_width('unfolded_stat_err', **hbc_args).Clone("pdf_%d" % ibin_var)

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
        if isinstance(self.binning_handler.get_binning_scheme('generator'), PtVarPerPtBinning):
            warnings.warn(cu.pcolors.RED + "Cannot do setup_absolute_results_per_lambda_bin() with PtVarPerPtBinning" + cu.pcolors.ENDC)
            return

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
            # # the 0.001 is to ensure we're def inside this bin
            # start_bin = binning.GetGlobalBinNumber(var+0.001, pt_bins[0]+0.001)
            # end_bin = binning.GetGlobalBinNumber(var+0.001, pt_bins[-2]+0.001)  # -2 since the last one is the upper edge of the last bin
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


class HistBinChopper(object):
    """Get histogram for pt or variable bin, and cache it in dict, so can be used later"""

    def __init__(self, binning_handler=None):
        self.binning_handler = binning_handler

        self.objects = {}
        self._cache = {}
        self._cache_integral = {}

    def add_obj(self, name, obj):
        # TODO: allow overwrite?
        if name not in self.objects and obj is not None:
            self.objects[name] = obj

    def update(self, other):
        """Update this object's cached things with those from another HistBinChopper"""
        self.objects.update(other.objects)
        self._cache.update(other._cache)
        self._cache_integral.update(other._cache_integral)

    def get_var_hist_pt_binned(self, hist1d, ibin_pt, binning_scheme='generator', is_signal_region=True):
        """Get hist of variable for given pt bin from massive 1D hist that TUnfold makes

        Parameters
        ----------
        hist1d : ROOT.TH1D
            Big 1D histogram to chop up
        ibin_pt : int
            Pt bin number: 0-indexed for each of signal region, underflow region
            # TODO: make this a physical value instead?
        binning_scheme : str, optional
            "generator" or "detector" for respective binning scheme
        is_signal_region : bool, optional
            True if signal region, false for pt underflow region

        Returns
        -------
        ROOT.TH1D
        """
        binning_obj = self.binning_handler.get_binning_scheme(binning_scheme)
        pt = binning_obj.get_pt_bins(is_signal_region)[ibin_pt]
        var_bins = binning_obj.get_variable_bins(pt)
        h = ROOT.TH1D("h_%d_%s" % (ibin_pt, cu.get_unique_str()), "", len(var_bins)-1, np.array(var_bins))
        for var_ind, var_value in enumerate(var_bins[:-1], 1):
            bin_num = binning_obj.physical_bin_to_global_bin(var=var_value+1E-6, pt=pt+1E-6)  # +1E-6 to ensure its inside
            h.SetBinContent(var_ind, hist1d.GetBinContent(bin_num))
            h.SetBinError(var_ind, hist1d.GetBinError(bin_num))
        return h

    def get_var_2d_hist_pt_binned(self, hist2d, ibin_pt, binning_scheme='generator', is_signal_region=True):
        """Get 2d hist for given pt bin from massive 2D hist

        Same options as get_var_hist_pt_binned()
        """
        binning_obj = self.binning_handler.get_binning_scheme(binning_scheme)
        pt = binning_obj.get_pt_bins(is_signal_region)[ibin_pt]
        var_bins = binning_obj.get_variable_bins(pt)
        h = ROOT.TH2D("h2d_%d_%s" % (ibin_pt, cu.get_unique_str()), "", len(var_bins)-1, var_bins, len(var_bins)-1, var_bins)
        for var_ind_x, var_value_x in enumerate(var_bins[:-1], 1):
            bin_num_x = binning_obj.physical_bin_to_global_bin(var=var_value_x+1E-6, pt=pt+1E-6)  # +1E-6 to ensure its inside
            for var_ind_y, var_value_y in enumerate(var_bins[:-1], 1):
                bin_num_y = binning_obj.physical_bin_to_global_bin(var=var_value_y+1E-6, pt=pt+1E-6)
                h.SetBinContent(var_ind_x, var_ind_y, hist2d.GetBinContent(bin_num_x, bin_num_y))
                h.SetBinError(var_ind_x, var_ind_y, hist2d.GetBinError(bin_num_x, bin_num_y))
        return h

    def get_pt_hist_var_binned(self, hist1d, ibin_var, binning_scheme='generator', is_signal_region=True):
        """Get hist of pt for given variable bin from massive 1D hist that TUnfold makes

        NB only makes sense if this variable bin exists across all pt bins, otherwise impossible to do
        """
        binning_obj = self.binning_handler.get_binning_scheme(binning_scheme)
        pt_bins = binning_obj.get_pt_bins(is_signal_region)
        if isinstance(binning_obj, PtVarPerPtBinning):
            # check if this variable bin exists for all pt bins
            var_value = [binning_obj.get_variable_bins(pt)[ibin_var] for pt in pt_bins]
            if len(set(var_value)) != 1:
                raise ValueError("Cannot use get_pt_hist_var_binned: binning scheme has different bin edge for ibin_var=%d" % ibin_var)

        var_value = binning_obj.get_variable_bins(pt_bins[0])[ibin_var]
        h = ROOT.TH1D("h_%d_%s" % (ibin_var, cu.get_unique_str()), "", len(pt_bins)-1, pt_bins)
        for pt_ind, pt_value in enumerate(pt_bins[:-1], 1):
            bin_num = binning_obj.physical_bin_to_global_bin(var=var_value+1E-6, pt=pt_value+1E-6)
            h.SetBinContent(pt_ind, hist1d.GetBinContent(bin_num))
            h.SetBinError(pt_ind, hist1d.GetBinError(bin_num))
        return h

    def get_pt_2d_hist_var_binned(self, hist2d, ibin_var, binning_scheme='generator', is_signal_region=True):
        """Get 2d hist for given variable bin from massive 2D hist"""
        binning_obj = self.binning_handler.get_binning_scheme(binning_scheme)
        pt_bins = binning_obj.get_pt_bins(is_signal_region)
        if isinstance(binning_obj, PtVarPerPtBinning):
            # check if this variable bin exists for all pt bins
            var_value = [binning_obj.get_variable_bins(pt)[ibin_var] for pt in pt_bins]
            if len(set(var_value)) != 1:
                raise ValueError("Cannot use get_pt_hist_var_binned: binning scheme has different bin edge for ibin_var=%d" % ibin_var)

        var_value = binning_obj.get_variable_bins(pt_bins[0])[ibin_var]
        h = ROOT.TH1D("h2d_%d_%s" % (ibin_var, cu.get_unique_str()), "", len(pt_bins)-1, pt_bins)
        for pt_ind_x, pt_value_x in enumerate(pt_bins[:-1], 1):
            bin_num_x = binning_obj.physical_bin_to_global_bin(var=var_value+1E-6, pt=pt_value_x+1E-6)
            for pt_ind_y, pt_value_y in enumerate(pt_bins[:-1], 1):
                bin_num_y = binning_obj.physical_bin_to_global_bin(var=var_value+1E-6, pt=pt_value_y+1E-6)
                h.SetBinContent(pt_ind_x, pt_ind_y, hist2d.GetBinContent(bin_num_x, bin_num_y))
                h.SetBinError(pt_ind_x, pt_ind_y, hist2d.GetBinError(bin_num_x, bin_num_y))
        return h

    def get_bin_plot(self, name, ind, axis, do_norm=False, do_div_bin_width=False, binning_scheme='generator', is_signal_region=True):
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

        key = self._generate_key(name, ind, axis, do_norm, do_div_bin_width, binning_scheme, is_signal_region)
        if key not in self._cache:
            if axis == 'lambda':
                self._cache[key] = self.get_pt_hist_var_binned(hist1d=self.objects[name],
                                                               ibin_var=ind,
                                                               binning_scheme=binning_scheme,
                                                               is_signal_region=is_signal_region)
            else:
                self._cache[key] = self.get_var_hist_pt_binned(hist1d=self.objects[name],
                                                               ibin_pt=ind,
                                                               binning_scheme=binning_scheme,
                                                               is_signal_region=is_signal_region)

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
    def _generate_key(name, ind, axis, do_norm, do_div_bin_width, binning_scheme, is_signal_region=True):
        """Generate consistent name for these args, options as in get_bin_plot()"""
        if axis not in ['pt', 'lambda']:
            raise ValueError('_generate_key(): axis must be "pt" or "lambda"')
        key = name + "_%s_bin_%d_%s" % (axis, ind, binning_scheme)
        if do_norm:
            key += "_norm"
        if do_div_bin_width:
            key += "_divBinWidth"
        if not is_signal_region:
            key += "_ptUnderflow"
        return key

    # TODO: remove these? just use get_bin_plot instead?
    def get_pt_bin(self, name, ind, binning_scheme='generator', is_signal_region=True):
        return self.get_bin_plot(name=name, ind=ind, axis='pt', do_norm=False, do_div_bin_width=False, binning_scheme=binning_scheme, is_signal_region=is_signal_region)

    def get_pt_bin_div_bin_width(self, name, ind, binning_scheme='generator', is_signal_region=True):
        return self.get_bin_plot(name=name, ind=ind, axis='pt', do_norm=False, do_div_bin_width=True, binning_scheme=binning_scheme, is_signal_region=is_signal_region)

    def get_pt_bin_normed(self, name, ind, binning_scheme='generator', is_signal_region=True):
        return self.get_bin_plot(name=name, ind=ind, axis='pt', do_norm=True, do_div_bin_width=False, binning_scheme=binning_scheme, is_signal_region=is_signal_region)

    def get_pt_bin_normed_div_bin_width(self, name, ind, binning_scheme='generator', is_signal_region=True):
        return self.get_bin_plot(name=name, ind=ind, axis='pt', do_norm=True, do_div_bin_width=True, binning_scheme=binning_scheme, is_signal_region=is_signal_region)

    def get_lambda_bin(self, name, ind, binning_scheme='generator', is_signal_region=True):
        return self.get_bin_plot(name=name, ind=ind, axis='lambda', do_norm=False, do_div_bin_width=False, binning_scheme=binning_scheme, is_signal_region=is_signal_region)

    def get_lambda_bin_div_bin_width(self, name, ind, binning_scheme='generator', is_signal_region=True):
        return self.get_bin_plot(name=name, ind=ind, axis='lambda', do_norm=False, do_div_bin_width=True, binning_scheme=binning_scheme, is_signal_region=is_signal_region)

    def get_lambda_bin_normed(self, name, ind, binning_scheme='generator', is_signal_region=True):
        return self.get_bin_plot(name=name, ind=ind, axis='lambda', do_norm=True, do_div_bin_width=False, binning_scheme=binning_scheme, is_signal_region=is_signal_region)

    def get_lambda_bin_normed_div_bin_width(self, name, ind, binning_scheme='generator', is_signal_region=True):
        return self.get_bin_plot(name=name, ind=ind, axis='lambda', do_norm=True, do_div_bin_width=True, binning_scheme=binning_scheme, is_signal_region=is_signal_region)


class ExpSystematic(object):
    """Class to hold info about an experimental systematic"""

    def __init__(self, label, syst_map=None, syst_shift=None, syst_shifted=None, syst_ematrix=None, syst_error_bar=None):
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

        # nominal hist with shift as error bar
        self.syst_error_bar = syst_error_bar
        self.syst_error_bar_label = 'syst_err_%s' % (self.label_no_spaces)


class TruthTemplateMaker(object):
    """Create truth-level data template from fitting MCs templates at detector level,
    fitting per gen pT bin.

    Yes, "template" is used for both the MC shapes being fitted, and the final
    distribution. I'm lazy.
    """

    def __init__(self, binning_handler, output_dir):
        self.binning_handler = binning_handler
        self.variable_name = binning_handler.variable_name
        self.hist_bin_chopper = HistBinChopper(binning_handler=self.binning_handler)

        self.templates = []  # store MC template to be fitted
        self.data_label_reco = "hist_data_reco"
        self.data_label_gen = "hist_data_gen"
        # for fitting to reco level
        self.fits = []
        self.fit_results = []
        # for fitting to gen level as cross-check
        self.fits_gen = []
        self.fit_results_gen = []

        self.components = []  # scaled MC components for each pT bin
        self.output_dir = output_dir  # for plots
        self.truth_template = None  # final 1D template

    def set_input(self, data_hist_reco):
        """Set thing that gets fitted. Should be at detector-level, but with gen binning"""
        self.hist_bin_chopper.add_obj(self.data_label_reco, data_hist_reco)

    def set_input_gen(self, data_hist_gen):
        """Set thing to be fitted at gen level - for cross-checking with MC"""
        self.hist_bin_chopper.add_obj(self.data_label_gen, data_hist_gen)

    def add_mc_template(self, name, hist_reco, hist_gen, colour=None):
        """Add MC template to be part of the fit.
        Both hists should have gen binning.
        """
        reco_label = "mc_reco_%s" % (name)
        gen_label = "mc_gen_%s" % (name)
        self.templates.append(dict(
            name=name,
            reco_label=reco_label,
            hist_reco=hist_reco,
            gen_label=gen_label,
            hist_gen=hist_gen,
            colour=colour
        ))
        self.hist_bin_chopper.add_obj(reco_label, hist_reco.Clone())
        self.hist_bin_chopper.add_obj(gen_label, hist_gen.Clone())

    # Create fit function from templates
    # Should be used with functools.partial
    @staticmethod
    def data_distribution_fn(x, pars, hists):
        xx = x[0]
        # get content of the histograms for this point in all hists,
        # do linear sum
        return sum([pars[i]*h.GetBinContent(h.GetXaxis().FindFixBin(xx))
                    for i, h in enumerate(hists)])

    def create_template(self):
        """Create gen level template using MCs, by fitting to data at detector level

        Note that all hists should use gen binning, since that's what we need
        to produce the gen-level template

        We do this per pT bin invidivdually, then create one big distribution

        """
        # do the fits per pT bin, including plotting the fit
        self.do_fits()
        # plot scale factors vs pt bin
        first_bin = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=False)[0]
        last_bin = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=True)[-1]
        self.plot_fit_results_vs_pt_bin("Fit to reco data, %g < p_{T} < %G GeV" % (first_bin, last_bin),
                                        os.path.join(self.output_dir, "reco_gen_bin_fit_factors.pdf"))
        self.construct_truth_template()
        return self.truth_template

    def plot_fit(self, hist_data, mc_hists, labels, title, filename, other_contributions=None):
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        this_hist_data = hist_data.Clone("Data")
        this_hist_data.SetMarkerStyle(cu.Marker.get('triangleUp'))
        # draw first to get stats box as not drawn with THStack
        this_hist_data.Draw()
        canv.Update()
        all_func_hists = []
        stats_box = this_hist_data.GetListOfFunctions().FindObject("stats")
        func = this_hist_data.GetListOfFunctions().At(0)
        hist_fit = hist_data.Clone("Fit")
        for i in range(1, hist_data.GetNbinsX()+1):
            hist_fit.SetBinContent(i, func.Eval(hist_data.GetBinCenter(i)))
            hist_fit.SetBinError(i, 0)

        lw = 2

        entries = [
            Contribution(h,
                         label=labels[ind],
                         line_color=self.templates[ind]['colour'],
                         line_width=lw,
                         marker_color=self.templates[ind]['colour'])
            for ind, h in enumerate(mc_hists)
        ]
        entries.extend([
            Contribution(hist_data,
                         label="Data",
                         line_color=ROOT.kBlack,
                         line_width=lw,
                         marker_color=ROOT.kBlack,
                         marker_style=cu.Marker.get('triangleUp')),
            Contribution(hist_fit,
                         label=hist_fit.GetName(),
                         line_color=ROOT.kAzure+1,
                         line_width=lw,
                         line_style=2,
                         marker_color=ROOT.kAzure+1,
                         subplot=hist_data)
        ])

        if other_contributions:
            entries.extend(other_contributions)

        subplot_lim = [0.5, 1.5]
        # set subplot limits to cover all contributions
        y_min, y_max = 999, -999
        for ent in entries:
            if ent.subplot is None:
                continue
            ratio = ent.obj.Clone()
            ratio.Divide(ent.subplot)
            y_max = max(y_max, ratio.GetMaximum())
            y_min = max(y_min, ratio.GetMinimum(1E-20))

        subplot_lim = [min(0.5, y_min), max(1.5, y_max)]
        # subplot_lim = [y_min, y_max]

        plot = Plot(entries, what='hist',
                    title=title,
                    xtitle=self.variable_name,
                    ytitle="N / bin width",
                    subplot_type="ratio",
                    subplot_title="Fit / data",
                    subplot_limits=subplot_lim)
        plot.plot("NOSTACK HIST E")
        plot.main_pad.cd()
        stats_box.SetBorderSize(0)
        stats_box.SetFillStyle(0)
        stats_box.SetY1NDC(0.65)
        stats_box.SetY2NDC(0.83)
        stats_box.SetX2NDC(0.9)
        stats_box.Draw()
        plot.legend.SetY2NDC(0.75)
        plot.legend.SetY1NDC(0.5)
        plot.legend.SetX2NDC(0.9)
        plot.save(filename)

    def do_fits(self):
        """Do the MC fits to data one per gen pT bin"""

        # reset result holders
        self.fits = []
        self.fit_results = []
        self.components = []

        n_components = len(self.templates)

        if n_components < 2:
            raise RuntimeError("Cannot do git with < 2 MC components")

        labels = [t['name'] for t in self.templates]
        # Fit underflow pT bins
        pt_bins_uflow = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=False)
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(pt_bins_uflow[:-1], pt_bins_uflow[1:])):
            print("Fitting uflow", bin_edge_low, bin_edge_high)
            hbc_args = dict(ind=ibin, binning_scheme='generator', is_signal_region=False)
            hist_data = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.data_label_reco, **hbc_args).Clone()
            if hist_data.Integral() == 0:
                self.fits.append(None)
                self.fit_results.append(None)
                self.components.append(None)
                print("... empty data, skipping fit")
                continue
            mc_hists = [self.hist_bin_chopper.get_pt_bin_div_bin_width(t['reco_label'], **hbc_args).Clone()
                        for t in self.templates]
            xmin, xmax = 0, hist_data.GetNbinsX()
            # create TF1 using MC histograms
            f = ROOT.TF1("reco_fit_gen_ubin_%d" % ibin, partial(TruthTemplateMaker.data_distribution_fn, hists=mc_hists), xmin, xmax, n_components)
            f.SetNpx(10000)
            # Actually do the fit to data
            fit_result = hist_data.Fit(f, "QEMS")
            self.fits.append(f)
            self.fit_results.append(fit_result)
            # scaled versions of component hists
            for i in range(n_components):
                mc_hists[i].Scale(f.GetParameter(i))
            self.plot_fit(hist_data,
                          mc_hists,
                          labels=labels,
                          title="Fit to reco data\n%g < p_{T} < %g GeV" % (bin_edge_low, bin_edge_high),
                          filename=os.path.join(self.output_dir, "reco_fit_gen_bin_uflow_%d.pdf" % ibin))

        # Fit signal pT bins
        pt_bins_signal = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=True)
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(pt_bins_signal[:-1], pt_bins_signal[1:])):
            print("Fitting", bin_edge_low, bin_edge_high)
            hbc_args = dict(ind=ibin, binning_scheme='generator', is_signal_region=True)
            hist_data = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.data_label_reco, **hbc_args).Clone()
            mc_hists = [self.hist_bin_chopper.get_pt_bin_div_bin_width(template['reco_label'], **hbc_args).Clone()
                        for template in self.templates]
            xmin, xmax = 0, hist_data.GetNbinsX()
            # create TF1 using MC histograms
            f = ROOT.TF1("reco_fit_gen_bin_%d" % ibin, partial(TruthTemplateMaker.data_distribution_fn, hists=mc_hists), xmin, xmax, len(mc_hists))
            f.SetNpx(10000)
            # Actually do the fit to data
            fit_result = hist_data.Fit(f, "QEMS")
            self.fits.append(f)
            self.fit_results.append(fit_result)
            # scaled versions of component hists
            for i in range(n_components):
                mc_hists[i].Scale(f.GetParameter(i))
            self.plot_fit(hist_data,
                          mc_hists,
                          labels=labels,
                          title="Fit to reco data\n%g < p_{T} < %g GeV" % (bin_edge_low, bin_edge_high),
                          filename=os.path.join(self.output_dir, "reco_fit_gen_bin_%d.pdf" % ibin))

    def plot_fit_results_vs_pt_bin(self, title, output_filename):
        n = len(self.fits)
        x = array('d', list(range(n)))
        ex = array('d', [0 for f in self.fits])
        multi_gr = ROOT.TMultiGraph(cu.get_unique_str(), "%s;pt bin index;Fit parameter" % title)
        marker = cu.Marker()
        for ind_t, (template, mark) in enumerate(zip(self.templates, marker.cycle())):
            y = array('d', [f.GetParameter(ind_t) if f else 0 for f in self.fits])
            ey = array('d', [f.GetParError(ind_t) if f else 0 for f in self.fits])
            gr = ROOT.TGraphErrors(n, x, y, ex, ey)
            gr.SetTitle(template['name'])
            gr.SetMarkerStyle(mark)
            col = template['colour']
            if col:
                gr.SetMarkerColor(col)
                gr.SetLineColor(col)
            multi_gr.Add(gr)

        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        canv.SetTicks(1, 1)
        multi_gr.Draw("ALP")
        multi_leg = ROOT.gPad.BuildLegend()
        multi_leg.SetFillStyle(0)
        canv.SaveAs(output_filename)

    def construct_truth_template(self):
        # Construct a new truth-level distribution using the fit factors
        self.truth_template = None

        new_truth_hists = []

        # underflow pt bins
        pt_bins_uflow = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=False)
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(pt_bins_uflow[:-1], pt_bins_uflow[1:])):
            print("Creating template", ibin, bin_edge_low, bin_edge_high)

            hbc_args = dict(ind=ibin, binning_scheme='generator', is_signal_region=False)
            sum_hist = None

            if ibin == 0:
                # add the lowest pt bin manually since it isn't filled at reco level
                # extrapolate from lowest pt bins to get fit factors for this bin
                last_ind = 3
                degree = 1
                pt_params = [0.5*(pt_bins_uflow[i]+pt_bins_uflow[i+1]) for i in range(1, last_ind)]
                center_bin_0 = 0.5*(pt_bins_uflow[0] + pt_bins_uflow[1])

                for ind_t, template in enumerate(self.templates):
                    params = [f.GetParameter(ind_t) for f in self.fits[1:last_ind]]
                    fit_coeff = np.polyfit(pt_params, params, deg=degree)
                    w = np.poly1d(fit_coeff)(center_bin_0)
                    print("Extrapolated w%d for lowest pt bin:" % ind_t, w)
                    # note no div bin width, as that's what TUnfold uses
                    hist_gen = self.hist_bin_chopper.get_pt_bin(template['gen_label'], **hbc_args).Clone()
                    hist_gen.Scale(w)
                    if sum_hist is None:
                        sum_hist = hist_gen
                    else:
                        sum_hist.Add(hist_gen)

            else:
                # Get the scaled truth-level distributions, using fit factors
                # from doing the fit
                # Sum together into one hist for this pt bin
                for ind_t, template in enumerate(self.templates):
                    # note no div bin width, as that's what TUnfold uses
                    hist_gen = self.hist_bin_chopper.get_pt_bin(template['gen_label'], **hbc_args).Clone()
                    print("Scaling template", ind_t, "by", self.fits[ibin].GetParameter(ind_t))
                    hist_gen.Scale(self.fits[ibin].GetParameter(ind_t))
                    if sum_hist is None:
                        sum_hist = hist_gen
                    else:
                        sum_hist.Add(hist_gen)

            sum_hist.SetName("template_truth_%d" % ibin)
            new_truth_hists.append(sum_hist)

        # signal pt bins
        global_ibin = len(pt_bins_uflow) - 1
        pt_bins_signal = self.binning_handler.get_pt_bins(binning_scheme="generator", is_signal_region=True)
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(pt_bins_signal[:-1], pt_bins_signal[1:])):
            print("Creating template", ibin, global_ibin, bin_edge_low, bin_edge_high)

            hbc_args = dict(ind=ibin, binning_scheme='generator', is_signal_region=True)
            sum_hist = None

            # Get the scaled truth-level distributions, using fit factors
            # from doing the fit
            # Sum together into one hist for this pt bin
            for ind_t, template in enumerate(self.templates):
                # note no div bin width, as that's what TUnfold uses
                hist_gen = self.hist_bin_chopper.get_pt_bin(template['gen_label'], **hbc_args).Clone()
                print("Scaling template", ind_t, "by", self.fits[global_ibin].GetParameter(ind_t))
                hist_gen.Scale(self.fits[global_ibin].GetParameter(ind_t))
                if sum_hist is None:
                    sum_hist = hist_gen
                else:
                    sum_hist.Add(hist_gen)
            sum_hist.SetName("template_truth_%d" % global_ibin)
            new_truth_hists.append(sum_hist)
            global_ibin += 1

        # Create 1 big absolute distribution at gen level by combining all these
        # individual binned hists
        # It's a bit convoluted to get the correct binning
        truth_template = self.templates[0]['hist_gen'].Clone("truth_template")
        truth_template.Reset()

        print(type(pt_bins_uflow), pt_bins_uflow)
        print(type(pt_bins_signal), pt_bins_signal)
        all_pt_bins = list(chain(pt_bins_uflow[:-1], pt_bins_signal))
        print(all_pt_bins)
        for pt_ind, (pt_low, pt_high) in enumerate(zip(all_pt_bins[:-1], all_pt_bins[1:])):
            var_bins = self.binning_handler.get_variable_bins(pt=pt_low+1E-6, binning_scheme='generator')
            for bin_ind, var_value in enumerate(var_bins, 1):
                # bin_ind refers to bin in the template TH1 (hence start at 1), global_bin refers to global bin number
                global_bin = self.binning_handler.physical_bin_to_global_bin(pt=pt_low + 1E-6, var=var_value + 1E-6, binning_scheme='generator')
                truth_template.SetBinContent(global_bin, new_truth_hists[pt_ind].GetBinContent(bin_ind))
                truth_template.SetBinError(global_bin, new_truth_hists[pt_ind].GetBinError(bin_ind))
            # TODO: deal with overflow?

        self.truth_template = truth_template
        self.hist_bin_chopper.add_obj("truth_template", self.truth_template)

    def check_template_at_gen(self):
        """Do fits at gen level to input from set_input_gen(), and compare with
        fit factors from reco fit"""
        if self.data_label_gen not in self.hist_bin_chopper.objects:
            raise RuntimeError("No %s: you need to call set_input_gen() first" % self.data_label_gen)

        if self.truth_template is None:
            raise RuntimeError("No truth template: need to call create_template() first")

        self.do_gen_fits()
        # plot scale factors vs pt bin, along with reco ones
        first_bin = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=False)[0]
        last_bin = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=True)[-1]
        self.plot_gen_vs_reco_fit_results_vs_pt_bin("Comparing fit to reco and gen, %g < p_{T} < %G GeV" % (first_bin, last_bin),
                                                    os.path.join(self.output_dir, "reco_vs_gen_bin_fit_factors.pdf"))

    def do_gen_fits(self):
        """Do the MC fits to gen-level data one per gen pT bin"""

        # reset result holders
        self.fits_gen = []
        self.fit_results_gen = []

        n_components = len(self.templates)

        if n_components < 2:
            raise RuntimeError("Cannot do git with < 2 MC components")

        labels = [t['name'] for t in self.templates]
        # Fit underflow pT bins
        pt_bins_uflow = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=False)
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(pt_bins_uflow[:-1], pt_bins_uflow[1:])):
            print("Fitting uflow", bin_edge_low, bin_edge_high)
            hbc_args = dict(ind=ibin, binning_scheme='generator', is_signal_region=False)
            hist_data = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.data_label_gen, **hbc_args).Clone()
            if hist_data.Integral() == 0:
                self.fits_gen.append(None)
                self.fit_results_gen.append(None)
                print("... empty data, skipping fit")
                continue
            xmin, xmax = 0, hist_data.GetNbinsX()
            mc_hists = [self.hist_bin_chopper.get_pt_bin_div_bin_width(t['gen_label'], **hbc_args).Clone()
                        for t in self.templates]
            # create TF1 using MC histograms
            f = ROOT.TF1("gen_fit_gen_ubin_%d" % ibin, partial(TruthTemplateMaker.data_distribution_fn, hists=mc_hists), xmin, xmax, n_components)
            f.SetNpx(10000)
            # Actually do the fit to data
            fit_result = hist_data.Fit(f, "QEMS")
            self.fits_gen.append(f)
            self.fit_results_gen.append(fit_result)
            # scaled versions of component hists
            for i in range(n_components):
                mc_hists[i].Scale(f.GetParameter(i))

            # Get the fit from reco fit values
            other_contributions = [
                Contribution(self.hist_bin_chopper.get_pt_bin_div_bin_width("truth_template", **hbc_args),
                             label="Fit using reco fit params",
                             line_width=2,
                             line_color=ROOT.kOrange-3,
                             line_style=3,
                             subplot=hist_data)
            ]
            self.plot_fit(hist_data,
                          mc_hists,
                          other_contributions=other_contributions,
                          labels=labels,
                          title="Fit to gen data\n%g < p_{T} < %g GeV" % (bin_edge_low, bin_edge_high),
                          filename=os.path.join(self.output_dir, "gen_fit_gen_bin_uflow_%d.pdf" % ibin))

        # Fit signal pT bins
        pt_bins_signal = self.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=True)
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(pt_bins_signal[:-1], pt_bins_signal[1:])):
            print("Fitting", bin_edge_low, bin_edge_high)
            hbc_args = dict(ind=ibin, binning_scheme='generator', is_signal_region=True)
            hist_data = self.hist_bin_chopper.get_pt_bin_div_bin_width(self.data_label_gen, **hbc_args).Clone()
            mc_hists = [self.hist_bin_chopper.get_pt_bin_div_bin_width(template['gen_label'], **hbc_args).Clone()
                        for template in self.templates]
            xmin, xmax = 0, hist_data.GetNbinsX()
            # create TF1 using MC histograms
            f = ROOT.TF1("gen_fit_gen_bin_%d" % ibin, partial(TruthTemplateMaker.data_distribution_fn, hists=mc_hists), xmin, xmax, len(mc_hists))
            f.SetNpx(10000)
            # Actually do the fit to data
            fit_result = hist_data.Fit(f, "QEMS")
            self.fits_gen.append(f)
            self.fit_results_gen.append(fit_result)
            # scaled versions of component hists
            for i in range(n_components):
                mc_hists[i].Scale(f.GetParameter(i))

            # Get the fit from reco fit values
            other_contributions = [
                Contribution(self.hist_bin_chopper.get_pt_bin_div_bin_width("truth_template", **hbc_args),
                             label="Fit using reco fit params",
                             line_width=2,
                             line_color=ROOT.kOrange-3,
                             line_style=3,
                             subplot=hist_data)
            ]
            self.plot_fit(hist_data,
                          mc_hists,
                          other_contributions=other_contributions,
                          labels=labels,
                          title="Fit to gen data\n%g < p_{T} < %g GeV" % (bin_edge_low, bin_edge_high),
                          filename=os.path.join(self.output_dir, "gen_fit_gen_bin_%d.pdf" % ibin))

    def plot_gen_vs_reco_fit_results_vs_pt_bin(self, title, output_filename):
        n = len(self.fits)
        x = array('d', list(range(n)))
        ex = array('d', [0 for f in self.fits])
        multi_gr = ROOT.TMultiGraph(cu.get_unique_str(), "%s;pt bin index;Fit parameter" % title)
        marker = cu.Marker()
        for ind_t, (template, mark) in enumerate(zip(self.templates, marker.cycle())):
            y = array('d', [f.GetParameter(ind_t) if f else 0 for f in self.fits])
            ey = array('d', [f.GetParError(ind_t) if f else 0 for f in self.fits])
            gr = ROOT.TGraphErrors(n, x, y, ex, ey)
            gr.SetTitle(template['name'] + " [RECO]")
            gr.SetMarkerStyle(mark)
            col = template['colour']
            if col:
                gr.SetMarkerColor(col)
                gr.SetLineColor(col)
            multi_gr.Add(gr)

        for ind_t, (template, mark) in enumerate(zip(self.templates, marker.cycle(filled=False))):
            y = array('d', [f.GetParameter(ind_t) if f else 0 for f in self.fits_gen])
            ey = array('d', [f.GetParError(ind_t) if f else 0 for f in self.fits_gen])
            gr = ROOT.TGraphErrors(n, x, y, ex, ey)
            gr.SetTitle(template['name'] + " [GEN]")
            gr.SetMarkerStyle(mark)
            col = template['colour']
            if col:
                gr.SetMarkerColor(col)
                gr.SetLineColor(col)
                gr.SetLineStyle(2)
            multi_gr.Add(gr)

        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        canv.SetTicks(1, 1)
        multi_gr.Draw("ALP")
        multi_leg = ROOT.gPad.BuildLegend()
        multi_leg.SetFillStyle(0)
        canv.SaveAs(output_filename)
