#!/usr/bin/env python2.6
"""
KCorrect utility routines for DEEP2
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


def deepmag_to_maggie(mag, mag_sigma=None):
    """
    Converts DEEP2 photometry to AB maggies

    This function converts DEEP2 photometry to AB maggies.  DEEP2
    photometry are assumed to be on the AB magnitude scale.

    If the uncertainty is not provided, the approximate uncertainties
    are computed based on the input photometry and the minimum
    uncertainties of 0.02 AB mag will be added in quadrature.  AB
    maggies and corresponding ivars are returned.

    INPUT

    mag -- List of DEEP2 magnitudes
    mag_sigma -- List of uncertainties in DEEP2 magnitudes; defaults to
                 a null list

    REFERENCE

    deep_to_maggies.pro of kcorrect v4_1_4
    """
    mag = np.asarray(mag, dtype=FTYPE)
    maggie = 10**(-0.4 * mag)

    if mag_sigma is None:
        # uncertainties not provided; estimate them

        mbase, dmbase = 24., 0.05
        min_error = 0.02
        maggie_ivar = (1. / (min_error**2
                             + (dmbase * 10**(-0.4 * (mbase - mag)))**2)
                       / (0.4 * np.log(10.))**2 / maggie**2)
        return maggie, maggie_ivar

    mag_sigma = np.asarray(mag_sigma, dtype=FTYPE)

    if mag_sigma.shape == mag.shape:
        # uncertainties are given, so use them
        maggie_ivar = 1. / (0.4 * np.log(10.) * maggie * mag_sigma)**2
    else:
        raise ValueError("Invalid mag_sigma array shape.")

    return maggie, maggie_ivar


class DEEPFilterList(FilterList):
    
    def __init__(self):
        filters = ['deep_B', 'deep_R', 'deep_I']
        super(DEEPFilterList, self).__init__(filters)


class DEEPKCorrect(KCorrect):
    """
    k-correction for DEEP BRI photometry.
    """

    def set_input(self, redshift, deepmag, deepmag_sigma=None, *args, **kwargs):
        """
        Set input data

        The input DEEP BRI photometry will be converted into AB
        maggies via DEEP.deepmag_to_maggie; see the documentation of
        the function for detail.

        INPUT

        redshift -- List of redshifts
        deepmag -- List of DEEP BRI magnitudes
        deepmag_sigma -- List of 1-sigma uncertainties in DEEP BRI
                         magnitudes; defaults to an empty list
        """
        filter_list = DEEPFilterList()
        maggie, maggie_ivar = deepmag_to_maggie(deepmag, deepmag_sigma)
        super(DEEPKCorrect, self).set_input(redshift, filter_list, maggie,
                                            maggie_ivar, *args, **kwargs)


class DEEPPhotoZ(PhotoZ):

    def set_input(self, deepmag, deepmag_sigma=None, *args, **kwargs):
        """
        Set input data

        The input DEEP BRI photometry will be converted into AB
        maggies via DEEP.deepmag_to_maggie; see the documentation of
        the function for detail.

        INPUT

        redshift -- List of redshifts
        deepmag -- List of DEEP BRI magnitudes
        deepmag_sigma -- List of 1-sigma uncertainties in DEEP BRI
                         magnitudes; defaults to an empty list
        """
        filter_list = DEEPFilterList()
        maggie, maggie_ivar = deepmag_to_maggie(deepmag, deepmag_sigma)
        super(DEEPPhotoZ, self).set_input(filter_list, maggie,
                                          maggie_ivar, *args, **kwargs)

