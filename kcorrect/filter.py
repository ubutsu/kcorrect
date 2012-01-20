#!/usr/bin/env python2.6
"""
Filter utilities.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
from kcorrect.globaldefs import DEFAULT_FILTER_DIR, FTYPE, FILTEREXT
from kcorrect.utils.io import path_to_file
from kcorrect.utils.yanny import YannyFileRead
import numpy as np
import os

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def load_filter(fname):
    """
    Load a filter transmission curve stored in a Yanny .par file.

    Parameters
    ----------
    fname : str
        A full path to a .par file.

    Notes
    -----
    See IDL procedure utils/k_load_filters.pro (kcorrect v4_2).
    """
    o = YannyFileRead(fname)

    # appears dangerous, but assume for now there is only one
    # structure per .par file corresponding to a response curve; the
    # IDL version uses yanny_readone, which assumes there is only one
    # structure anyways
    struct = o.structs[o.structs.keys()[0]]

    # use lambda and pass for response curve
    return (np.asarray(struct['lambda'], dtype=FTYPE),
            np.asarray(struct['pass'], dtype=FTYPE))


def get_filters(filter_dirs=None, full_path=False):
    """
    Get the list of currently available filters.

    This function searches for all files with their name ending .par
    in ${KCORRECT_DIR}/data/filters.  If filter_dirs is given, the
    files are searched in those directories.  Note that this function
    does not check whether a .par file contains a filter transmission
    information.

    Parameters
    ----------
    filter_dirs : list
        List of str paths to directories with filter .par files.
    full_path : bool
        True to return the full paths to the filters.

    Returns
    -------
    out : list
        List of filters available with FilterList.
    """
    filter_dirs = [DEFAULT_FILTER_DIR] if filter_dirs is None else filter_dirs
    filters = []
    for fdir in filter_dirs:
        for each in os.listdir(fdir):
            if not each.endswith(FILTEREXT):
                continue
            each = each if full_path else each[:-len(FILTEREXT)]
            filters.append('/'.join([fdir, each]) if full_path else each)
    filters.sort()
    return filters


class FilterList(object):
    """
    A set of filters is always maintained as a FilterList object.

    Parameters
    ----------
    filters : list
        List of filter file names (including .par extension).
    filter_dirs : list
        List of additional directories to search for filter files.
    """

    def __init__(self, filters, filter_dirs=None):
        filter_dirs = [] if filter_dirs is None else filter_dirs
        self.filter_dirs = filter_dirs + [DEFAULT_FILTER_DIR]

        # find full paths
        full_paths = []
        for each in filters:
            each = ''.join([each, FILTEREXT])
            fpath = path_to_file(each, self.filter_dirs)
            if not fpath:
                dirs = ', '.join(self.filter_dirs)
                raise ValueError('Filter %s not found in %s.' % (each, dirs))
            full_paths.append(fpath)

        self.filters = filters
        self.full_paths = full_paths

    def __eq__(self, other):
        if self.filter_dirs != other.filter_dirs:
            return False
        if self.full_paths != other.full_paths:
            return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return '('+ ', '.join(self.filters) +')'

    def __len__(self):
        return len(self.filters)

    def load_filters(self):
        """
        Load filter response curves.

        Returns
        -------
        The response curve (lambda, pass) for each filter will be
        returned as a tuple in order defined at instantiation.
        """
        ret = []
        for path in self.full_paths:
            x, y = load_filter(path)
            ret.append((x, y))
        return ret

    def plot(self):
        """
        Plot filter response curves.
        """
        for i, (x, y) in enumerate(self.load_filters()):
            plt.plot(x, y, '.-', label=self.filters[i])
        plt.xlabel('Wavelength [Angstrom]')
        plt.ylabel('Throughput')
        plt.legend(loc='best')
        plt.show()

    def fwhm(self):
        """
        Compute FWHM (min and max wavelengths).
        """
        fwhms = []
        for i, (x, y) in enumerate(self.load_filters()):
            hm = np.max(y) / 2.
            m = y > hm
            x = x[m]
            fwhms.append((x[0], x[-1]))
        return fwhms
