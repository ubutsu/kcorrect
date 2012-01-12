#!/usr/bin/env python2.6
"""
Utility functions for spectroscopy.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.globaldefs import FTYPE
from kcorrect.projectiontable import make_projection_table


def lamb_to_edges(lamb_c):
    """
    Convert pixel centers to pixels edges.

    Parameters
    ----------
    lamb_c : array_like
        Array of wavelength pixel centers.

    Notes
    -----
    It is assumed that pixel widths are of the same size.

    See IDL procedure utils/k_lambda_to_edges.pro (kcorrect v4_2).
    """
    lamb_c = np.asarray(lamb_c, dtype=FTYPE)
    n = lamb_c.size
    lamb_e = np.empty(n + 1, dtype=FTYPE)
    lamb_e[1:n] = 0.5 * ( lamb_c[0:n-1] + lamb_c[1:n] )
    lamb_e[0] = lamb_c[0] - (lamb_e[1] - lamb_c[0])
    lamb_e[n] = lamb_c[n-1] + (lamb_c[n-1] - lamb_e[n-1])
    return lamb_e


def project_filters(lambe, flux, filter_list):
    """
    Project a spectrum onto a set of filters.

    Parameters
    ----------
    lambe : array_like
        Wavelengths in angstroms at pixel edges [len(lambe) = len(flux) + 1]
    flux : array_like
        Flux in units of ergs/cm^2/s/A at pixel centers.
    filter_list : A FilterList object
        Filters through which the spectrum is projected.

    Notes
    -----
    See IDL procedure fit/k_project_filters.pro (kcorrect v4_2).
    """
    flux = np.asarray(flux, dtype=FTYPE)
    ndim = 2
    if len(flux.shape) == 1:
        flux = np.array([flux])
        ndim = 1
    elif len(flux.shape) == 0:
        raise ValueError('Too few data points.')
    if len(flux.shape) != 2:
        raise ValueError('Array must be two-dimensional.')
    r = make_projection_table((lambe, flux), filter_list, zrange=(0., 0., 1))
    return np.transpose(r[1])[0]
