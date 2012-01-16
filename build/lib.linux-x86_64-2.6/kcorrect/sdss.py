#!/usr/bin/env python2.6
"""
KCorrect utility routines for SDSS
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.kcorrect import KCorrect
from kcorrect.photoz import PhotoZ
from kcorrect.globaldefs import FTYPE
from kcorrect.filter import FilterList
from kcorrect.utils.photo import add_absigma_to_maggie_ivar
from kcorrect.utils.stats import ivar2sigma, sigma2ivar


class SDSSFilterList(FilterList):
    
    def __init__(self):
        filters = ['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']
        super(SDSSFilterList, self).__init__(filters)


class SDSS(object):
    """
    Utility class for SDSS 

    Following members define the AB conversions:

    names -- Filter names
  
    abfix   -- [du,dg,dr,di,dz] where in x band x(AB) = x(SDSS) + dx
    abfix02 -- Similar to abfix but obsolte year 2002 version

    photo_err -- Uncertainties in photometry

    bvalues -- Softening 'b' parameters used in the definition of Luptitude
    """
    names = ['u', 'g', 'r', 'i', 'z']

    abfix = [-0.036, 0.012, 0.010, 0.028,  0.040]
    abfix02 = [-0.042, 0.036, 0.015, 0.013, -0.002]

    photo_err = [0.05, 0.02, 0.02, 0.02, 0.03]

    bvalues = [1.4e-10, 0.9e-10, 1.2e-10, 1.8e-10, 7.4e-10]


def lup_to_maggie(lup, lup_sig=[], sdss_filters='ugriz', bvalues=SDSS.bvalues):
    """
    Convert luptitude (asinh mag) to AB maggie

    Conversion from luptitude to maggie is:

    maggie = 2 b sinh( - ln(b) - 0.4 ln(10) lup)

    This is linear when maggie ~ b, identical to maggie for maggie >> b.

    Conversion of the errors is:

    maggie_err = 2 b cosh( -ln(b) - 0.4 ln(10) lup) (0.4 ln 10) lup_err
    
    This function returns AB maggies for the given set of SDSS
    filters.  If lup_sig is also given, the function returns the
    inverse variances of AB maggies.  The size of each entry in the
    input luptitude list must be the same as the number of SDSS
    filters specified in sdss_filters as a string (e.g., 'ugriz',
    'urz', etc.).

    REFERENCE

    k_lups2maggies.pro of kcorrect v4_1_4

    INPUT

    lup -- List of SDSS luptitudes
    lup_sig -- 1 sigma uncertanties in luptitudes
    sdss_filters -- SDSS filters as a string; defaults to 'ugriz'
    bvalues -- Softening 'b' parameters; defaults to SDSS.bvalues
    """
    # check input data
    lup, lup_sig = np.atleast_2d(lup), np.atleast_2d(lup_sig)
    if len(sdss_filters) != lup.shape[1]:
        raise ValueError('Input lup size incompatible with sdss_filters.')

    # get bvalues
    bvalues = dict(zip(SDSS.names, bvalues))
    b = np.zeros(len(sdss_filters), dtype=FTYPE)
    for i in range(len(sdss_filters)):
        b[i] = bvalues[sdss_filters[i]]

    # compute maggie
    maggie = 2. * b * np.sinh(-np.log(b) - 0.4 * np.log(10.) * lup)
    if lup_sig.size != 0:
        # uncertainties given, so compute the ivar of maggies
        if lup.shape != lup_sig.shape:
            raise ValueError('Array size different for lup and lup_sig.')
        maggie_sigma = (2. * b * np.cosh(-np.log(b) - 0.4 * np.log(10.) * lup)
                        * 0.4 * np.log(10.) * lup_sig)
        return maggie, sigma2ivar(maggie_sigma)

    return maggie


def fix_lup_to_maggie(lup, lup_sig, photo_err=SDSS.photo_err, abfix=SDSS.abfix):
    """
    Convert SDSS database ugriz luptitude (asinh mag) to AB maggie

    This function is different from lup_to_maggie in that it corrects
    erratic entries.  It also adds a photometric calibration
    uncertainties and applies a correction to adjust to AB photometry.

    REFERENCE
  
    k_sdssfix.pro of kcorrect v4_1_4
    """
    # the following line should set ivar to zero if sigma <= 0, which
    # is basically what k_sdss_err2ivar.pro does.
    lup_ivar = sigma2ivar(lup_sig)
    mask = np.equal(lup, -9999.) * np.not_equal(lup_ivar, 0.)
    lup_ivar = np.where(mask, 0., lup_ivar)
    lup_sig = ivar2sigma(lup_ivar)
    maggie, maggie_ivar = lup_to_maggie(lup, lup_sig)
    mask = np.equal(lup_ivar, 0.)
    maggie = np.where(mask, 0., maggie)
    maggie_ivar = np.where(mask, 0., maggie_ivar)
    # add photometric calibration errors
    maggie_ivar = add_absigma_to_maggie_ivar(maggie, maggie_ivar, photo_err)
    # correct for AB offsets
    abfix = np.asarray(abfix, dtype=FTYPE)
    maggie = maggie * 10**(-0.4 * abfix)
    maggie_ivar = maggie_ivar * 10**(0.8 * abfix)
    return maggie, maggie_ivar


class SDSSKCorrect(KCorrect):

    def set_input(self, redshift, lup, lup_sig, red, *args, **kwargs):
        """
        Set input data

        The input luptitudes are reddening corrected and also
        processed via SDSS.fix_lup_to_maggie before the best fit SEDs
        are computed.

        INPUT

        redshift -- List of redshifts
        lup -- List of (uncorrected for reddening) luptitudes in ugriz
        lup_sig -- List of 1-sigma luptitude uncertainties in ugriz
        red -- List of reddening corrections in ugriz in magnitudes
        """
        filter_list = SDSSFilterList()
        lup = np.asarray(lup, dtype=FTYPE)
        red = np.asarray(red, dtype=FTYPE)

        # deredden the luptitudes; apparently this is okay in the
        # regime where the luptitudes and traditional log magnitudes
        # agree closely
        lup -= red

        # keep the reddening
        self.red = red

        # KCorrect.set_data expects maggies
        maggie, maggie_ivar = fix_lup_to_maggie(lup, lup_sig)
        
        super(SDSSKCorrect, self).set_input(redshift, filter_list,
                                            maggie, maggie_ivar,
                                            *args, **kwargs)


class SDSSPhotoZ(PhotoZ):

    def set_input(self, lup, lup_sig, red, *args, **kwargs):
        filter_list = SDSSFilterList()
        lup = np.asarray(lup, dtype=FTYPE)
        red = np.asarray(red, dtype=FTYPE)

        lup -= red

        self.red = red

        maggie, maggie_ivar = fix_lup_to_maggie(lup, lup_sig)
        
        super(SDSSPhotoZ, self).set_input(filter_list,
                                          maggie, maggie_ivar,
                                          *args, **kwargs)
  
