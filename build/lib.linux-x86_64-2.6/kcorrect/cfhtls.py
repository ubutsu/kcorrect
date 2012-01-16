#!/usr/bin/env python2.6
"""
KCorrect utility for CFHTLS.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.photoz import PhotoZ
from kcorrect.kcorrect import KCorrect
from kcorrect.globaldefs import FTYPE
from kcorrect.filter import FilterList


def cfhtlsmag_to_maggie(mag, mag_sigma=None):
    """
    Converts CFHTLS photometry to AB maggies.

    This function converts CFHTLS photometry to AB maggies.  CFHTLS
    photometry are assumed to be on the AB magnitude scale.

    If the uncertainty is not provided, the approximate uncertainties
    are computed based on the input photometry and the minimum
    uncertainties of 0.02 AB mag will be added in quadrature.

    Parameters
    ----------
    mag : array_like
        List of CFHTLS magnitudes.
    mag_sigma : array_like, optional
        List of uncertainties in DEEP2 magnitudes; defaults to a null
        list.

    Returns
    -------
    out : ndarray
        AB maggies and corresponding ivars are returned.

    Notes
    -----
    See deep_to_maggies.pro of kcorrect v4_2.
    """
    mag = np.asarray(mag, dtype=FTYPE)
    maggie = 10**(-0.4 * mag)

    if mag_sigma is None:
        # uncertainties not provided; estimate them.
        mbase, dmbase = 24., 0.05
        min_error = 0.02
        maggie_ivar = (1. / (min_error**2
                             + (dmbase * 10**(-0.4 * (mbase - mag)))**2)
                       / (0.4 * np.log(10.))**2 / maggie**2)
    else:
        mag_sigma = np.asarray(mag_sigma, dtype=FTYPE)
        if mag_sigma.shape != mag.shape:
            raise ValueError("Array shape mismatch for mag_sigma.")
        maggie_ivar = 1. / (0.4 * np.log(10.) * maggie * mag_sigma)**2

    return maggie, maggie_ivar


class CFHTLSFilterList(FilterList):
    
    def __init__(self):
        filters = ['sedfit_megacam_u',
                   'sedfit_megacam_g',
                   'sedfit_megacam_r',
                   'sedfit_megacam_i',
                   'sedfit_megacam_z',
                   'sedfit2_wircam_j',
                   'sedfit2_wircam_h',
                   'sedfit2_wircam_ks']
        super(CFHTLSFilterList, self).__init__(filters)


class CFHTLSPhotoZ(PhotoZ):

    def set_input(self, cfhtlsmag, cfhtlsmag_sigma=None, *args, **kwargs):
        filter_list = CFHTLSFilterList()
        maggie, maggie_ivar = cfhtlsmag_to_maggie(cfhtlsmag, cfhtlsmag_sigma)
        super(CFHTLSPhotoZ, self).set_input(filter_list, maggie,
                                            maggie_ivar, *args, **kwargs)

