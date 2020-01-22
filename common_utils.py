"""Set of common functions that are used in loads of scripts."""


import re
import ROOT
import os
from subprocess import call
from sys import platform as _platform
import numpy as np
import math
import argparse
from itertools import chain
from collections import OrderedDict


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.TH1.SetDefaultSumw2(True)


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Make argparse respect space formatting (newlines, etc) AND show defaults"""
    pass


def open_pdf(pdf_filename):
    """Open a PDF file using system's default PDF viewer."""
    if _platform.startswith("linux"):
        # linux
        call(["xdg-open", pdf_filename])
    elif _platform == "darwin":
        # OS X
        call(["open", pdf_filename])
    elif _platform == "win32":
        # Windows
        call(["start", pdf_filename])


#
# Filepath/directory fns
#
def cleanup_filepath(filepath):
    """Resolve any sym links, env vars, ~, etc, and return absolute path."""
    thing = os.path.abspath(os.path.expandvars(os.path.expanduser(filepath)))
    if os.path.islink(thing):
        # because otherwise readlink throws if not a link
        return os.readlink(thing)
    else:
        return thing


def get_full_path(filepath):
    """Return absolute directory of filepath.
    Resolve any environment vars, ~, sym links(?)"""
    return os.path.dirname(cleanup_filepath(filepath))


def check_file_exists(filepath):
    """Check if file exists. Can do absolute or relative file paths."""
    return os.path.isfile(cleanup_filepath(filepath))


def check_dir_exists(filepath):
    """Check if directory exists."""
    return os.path.isdir(cleanup_filepath(filepath))


def check_dir_exists_create(filepath):
    """Check if directory exists. If not, create it."""
    if not check_dir_exists(filepath):
        os.makedirs(cleanup_filepath(filepath))


#
# ROOT specific fns, like opening files safely
#


class TFileCacher(object):
    """Simple class to cache TFile & objects already "gotten" from TFile

    Can either do mytfile.get(x) or just mytfile[x]

    TODO: Should this just inherit from TFile instead?
    """
    def __init__(self, filename):
        self.tfile = open_root_file(filename)
        self.objects = {}

    def get(self, obj_name):
        return self[obj_name]

    def __getitem__(self, obj_name):
        # really hoping all these have unique names
        if obj_name not in self.objects:
            self.objects[obj_name] = get_from_tfile(self.tfile, obj_name)
        return self.objects[obj_name]


def open_root_file(filename, mode="READ"):
    """Safe way to open ROOT file. Could be improved."""
    clean_filename = cleanup_filepath(filename)
    if mode in ["READ", "UPDATE"]:
        if not check_file_exists(clean_filename):
            raise IOError("No such file %s" % clean_filename)
    f = ROOT.TFile(clean_filename, mode)
    if f.IsZombie() or not f:
        raise IOError("Can't open TFile %s" % clean_filename)
    return f


def exists_in_file(tfile, obj_name):
    """Check if object exists in TFile.

    Also handles directory structure, e.g. histograms/central/pt_1
    """
    parts = obj_name.split("/")
    current_obj = tfile
    for p in parts:
        if current_obj.GetListOfKeys().Contains(p):
            current_obj = current_obj.Get(p)
        else:
            return False
    return True


def get_from_tfile(tfile, obj_name, info=False):
    """Get some object from ROOT TFile with checks."""
    if info:
        print("Getting %s" % obj_name)
    obj = tfile.Get(obj_name)
    if obj == None:
        raise IOError("No object named %s in %s" % (obj_name, tfile.GetName()))
    else:
        return obj


def grab_obj_from_file(file_name, obj_name):
    """Get object names obj_name from ROOT file file_name"""
    input_file = open_root_file(file_name)
    obj = get_from_tfile(input_file, obj_name)
    if isinstance(obj, (ROOT.TH1, ROOT.TH2)):
        # THIS ORDER IS VERY IMPORTANT TO AVOID MEMORY LEAKS
        new_obj = obj.Clone(get_unique_str())
        new_obj.SetDirectory(0)  # Ownership kludge
        input_file.Close()
        return new_obj
    else:
        return obj


def check_exp(n):
    """
    Checks if number has stupidly larger exponent

    Can occur is using buffers - it just fills unused bins with crap
    """

    from math import fabs, log10, frexp
    m, e = frexp(n)
    return fabs(log10(pow(2, e))) < 10


def get_xy(graph):
    """
    Return lists of x, y points from a graph, because it's such a PITA
    """
    xarr = list(np.ndarray(graph.GetN(), 'd', graph.GetX()))
    yarr = list(np.ndarray(graph.GetN(), 'd', graph.GetY()))
    return xarr, yarr


def get_exey(graph):
    """
    Return lists of errors on x, y points from a graph, because it's such a PITA
    """
    xarr = list(np.ndarray(graph.GetN(), 'd', graph.GetEX()))
    yarr = list(np.ndarray(graph.GetN(), 'd', graph.GetEY()))
    return xarr, yarr


def th2_to_arr(h, do_errors=False):
    """Convert TH2 to 2D numpy array"""
    arr = np.zeros((h.GetNbinsX(), h.GetNbinsY()), dtype=float)
    err_arr = np.zeros((h.GetNbinsX(), h.GetNbinsY()), dtype=float)
    for x_ind in range(1, h.GetNbinsX() + 1):
        for y_ind in range(1, h.GetNbinsY() + 1):
            arr[x_ind-1][y_ind-1] = h.GetBinContent(x_ind, y_ind)
            if do_errors:
                err_arr[x_ind-1][y_ind-1] = h.GetBinError(x_ind, y_ind)
    return arr, err_arr


def make_normalised_TH2(hist, norm_axis, recolour=True, do_errors=False):
    norm_axis = norm_axis.upper()
    possible_opts = ['X', 'Y']
    if norm_axis not in possible_opts:
        raise RuntimeError("norm_axis must be one of %s" % possible_opts)
    norm_axis = norm_axis.upper()

    # easiest way to cope with x or y is to just get a 2D matrix of values,
    # can then do transpose if necessary
    arr, err_arr = th2_to_arr(hist, do_errors)

    if norm_axis == 'Y':
        arr = arr.T
        err_arr = err_arr.T
    if recolour:
        # can set so the maximum in each bin is the same,
        # scale other bins accordingly
        # this retain the colour scheme for each set of bins
        maxes = arr.max(axis=1)
        maxesT = maxes[:, None]
        arr = np.divide(arr, maxesT, where=maxesT!=0, out=np.zeros_like(arr))
        err_arr = np.divide(err_arr, maxesT, where=maxesT!=0, out=np.zeros_like(arr))
        # for ind, xbin in enumerate(arr):
        #     if xbin.max() > 0:
        #         factor = xbin.max()
        #         arr[ind] = xbin / factor
        #         err_arr[ind] = err_arr[ind] / factor
    else:
        # alternatively, can rescale so sum over bins = 1
        for ind, xbin in enumerate(arr):
            summed = arr.sum(axis=1)
            summedT = summed[:, None]
            arr = np.divide(arr, summedT, where=summedT!=0, out=np.zeros_like(arr))
            err_arr = np.divide(err_arr, summedT, where=summedT!=0, out=np.zeros_like(arr))
            # if xbin.sum() != 0:
            #     factor = xbin.sum()
            #     arr[ind] = xbin / factor
            #     err_arr[ind] /= factor

    if norm_axis == 'Y':
        arr = arr.T
        err_arr = err_arr.T

    # Create new hist object - MUST do it this way to get Z range correct
    new_histname = hist.GetName() + "_norm" + norm_axis
    # hnew = ROOT.TH2F(hist)  # TODO: determine if TH2F, TH2D...
    if type(hist) == ROOT.TH2F:
        hnew = ROOT.TH2F(hist)  # TODO: determine if TH2F, TH2D...
    elif type(hist) == ROOT.TH2D:
        hnew = ROOT.TH2D(hist)  # TODO: determine if TH2F, TH2D...
    else:
        raise RuntimeError("Unknown 2D hist type")
    hnew.SetName(new_histname)
    for x_ind, x_arr in enumerate(arr, 1):
        for y_ind, val in enumerate(x_arr, 1):
            hnew.SetBinContent(x_ind, y_ind, val)
            if do_errors:
                hnew.SetBinError(x_ind, y_ind, err_arr[x_ind-1][y_ind-1])
#     hnew.SetAxisRange(0.5, 1., 'Z')
    return hnew


def th1_to_arr(hist):
    return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX()+1)]), np.array([hist.GetBinError(i) for i in range(1, hist.GetNbinsX()+1)])


def get_list_of_element_names(thing):
    """Get list of key names for given thing in a ROOT file"""
    return [x.GetName() for x in thing.GetListOfKeys()]


def get_list_of_objects_in_dir(filename, dirname):
    """Get list of elements in a TDirectory"""
    f = open_root_file(filename)  # need to keep this in memory otherwise the get_list... segfaults
    d = get_from_tfile(f, dirname)
    return get_list_of_element_names(d)


def get_unique_str():
    return ROOT.TUUID().AsString()


def thstack_to_th1(hst):
    """Add all hists in a THStack into one TH1"""
    htotal = None
    hists = hst.GetHists()
    for h in hists:
        if htotal is None:
            htotal = h.Clone(get_unique_str())
        else:
            htotal.Add(h)
    return htotal


def get_hist_mean_rel_error(hist):
    """Get average relative error from histogram bins (i.e. error / contents)

    Parameters
    ----------
    hist : ROOT.TH1 (or descendents)
        Description

    Returns
    -------
    float
        average relative error
    """
    nbins = hist.GetNbinsX()
    contents = np.array([hist.GetBinContent(i) for i in range(1, nbins+1)])
    errors = np.array([hist.GetBinError(i) for i in range(1, nbins+1)])
    # only care about bins with something in them
    mask = contents > 0
    # if no bins have any contents, just return 0
    if not mask.any():
        return 0.
    rel_err = np.divide(errors[mask], contents[mask])
    rel_err[np.isinf(rel_err)] = 0.
    rel_err[np.isnan(rel_err)] = 0.
    return np.average(rel_err)


def get_jet_config_from_dirname(dirname):
    """Get jet configuration from string

    Parameters
    ----------
    dirname : str
        Description

    Returns
    -------
    (str, str)
        jet radius and PUS

    Raises
    ------
    RuntimeError
        If it can't find it
    """

    m = re.search(r'ak([\d]+)(chs|puppi)', dirname.lower())
    if m and len(m.groups()) == 2:
        return m.groups()
    else:
        print("dirname:", dirname)
        if m:
            print(m)
        raise RuntimeError("cannot determine jet config from dirname %s" % dirname)


def get_bin_edges(hist, axis):
    """Get array of bin edges from hist. Must specify which axis to use."""
    axis = axis.lower()
    if axis not in ['x', 'y']:
        raise RuntimeError("get_bin_edges axis must be x or y")
    ax = hist.GetXaxis() if axis == "x" else hist.GetYaxis()
    bins = [ax.GetBinLowEdge(i) for i in range(1, ax.GetNbins()+2)]
    return bins


class Marker(object):
    
    shape_dict = OrderedDict()
    shape_dict['circle'] = {'filled': 20, 'open': 24}
    shape_dict['square'] = {'filled': 21, 'open': 25}
    shape_dict['triangleUp'] = {'filled': 22, 'open': 26}
    shape_dict['triangleDown'] = {'filled': 23, 'open': 32}
    shape_dict['star'] = {'filled': 29, 'open': 30}
    shape_dict['diamond'] = {'filled': 33, 'open': 27}
    shape_dict['cross'] = {'filled': 34, 'open': 28}
    shape_dict['crossX'] = {'filled': 47, 'open': 46}
    shape_dict['doubleDiamond'] = {'filled': 43, 'open': 42}
    shape_dict['3star'] = {'filled': 39, 'open': 37}
    

    def __init__(self, shape='circle', filled=True):
        self.fill_state = filled
        self.shape_state = shape

    @staticmethod
    def get(shape, filled=True):
        if shape not in Marker.shape_dict.keys():
            raise RuntimeError("Unknown marker shape %s" % (shape))        
        return Marker.shape_dict[shape]['filled' if filled else 'open']

    def cycle(self, filled=True, cycle_filling=False, only_cycle_filling=False):
        if only_cycle_filling:
            self.fill_state = filled
            while True:
                yield self.get(self.shape_state, self.fill_state)
                self.fill_state = not self.fill_state
        else:
            for k in chain(Marker.shape_dict.keys()):
                self.fill_state = filled
                yield self.get(k, self.fill_state)
                if cycle_filling:
                    self.fill_state = not self.fill_state
                    yield self.get(k, self.fill_state)





