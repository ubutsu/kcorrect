#!/usr/bin/env python2.6
"""
Utility functions for statistics.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.globaldefs import FTYPE


def sigma2ivar(sigma):
    """
    Convert sigma to inverse variance.

    Converts sigma (i.e., 1 sigma uncertainty or standard deviation)
    into inverse variance thru ivar = 1 / sigma**2.  If sigma <= 0,
    ivar is set to zero for that entry.

    Parameters
    ----------
    sigma : array_like
        Input sigma array.

    Returns
    -------
    ivar : array_like
        Output inverse variance array.
    """
    sigma = np.ma.array(sigma, copy=False, dtype=FTYPE)
    sigma = np.ma.masked_less_equal(sigma, 0.)
    ivar = 1. / sigma / sigma
    return ivar.filled(0.)


def ivar2sigma(ivar):
    """
    Convert inverse variance into sigma.

    Converts inverse variance into sigma thru ivar = 1 / sigma**2.
    This routine sets sigma to zero if ivar is zero.

    Parameters
    ----------
    ivar : array_like
        Input inverse variance array.

    Returns
    -------
    ivar : array_like
        Output sigma array.
    """
    ivar = np.ma.array(ivar, copy=False, dtype=FTYPE)
    ivar = np.ma.masked_equal(ivar, 0.)
    sigma = 1. / np.ma.sqrt(ivar)
    return sigma.filled(0.)
