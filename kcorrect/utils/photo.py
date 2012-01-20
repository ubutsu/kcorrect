#!/usr/bin/env python2.6
"""
Utility functions for photometry.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.clib import k_fit_nonneg, k_fit_photoz, k_reconstruct_maggies
from kcorrect.globaldefs import FTYPE, KCORRECT_DIR
from kcorrect.utils.spec import project_filters, lamb_to_edges
from kcorrect.utils.io import read_basel


def add_absigma_to_maggie_ivar(maggie, maggie_ivar, absigma):
    """
    Add AB magnitude uncertainties to inverse variances in maggies.

    This function adds in quadrature the input AB magnitude
    uncertainties to inverse variances in maggies.  Returns the
    resulting inverse variances in maggies.

    Parameters
    ----------
    maggie, maggie_ivar : array_like
        AB maggies and their inverse variances.
    absigma : array_like
        1-sigma uncertainties in AB magnitudes.

    Notes
    -----
    See utils/k_minerror.pro of kcorrect v4_2.
    """
    min_err = np.ma.array(absigma, copy=False, dtype=FTYPE)
    maggie = np.ma.array(maggie, copy=False, dtype=FTYPE)
    maggie_ivar = np.ma.array(maggie_ivar, copy=False, dtype=FTYPE)
    maggie_ivar = np.ma.masked_equal(maggie_ivar, 0.)
    factor = 2.5 / np.ma.log(10.)
    err = factor / np.ma.sqrt(maggie_ivar) / maggie
    err2 = err**2 + min_err**2
    maggie_ivar = factor * factor / (maggie * maggie * err2)
    return maggie_ivar.filled(0.)


def vega2ab(filter_list, ref='kurucz'):
    """
    Return the offsets between Vega and AB magnitudes.

    This function returns a conversion from Vega manitudes to AB
    magnitudes such that

    m(AB) = m(Vega) + vega2ab

    where vega2ab is the offsets returned by this function.

    Parameters
    ----------
    filter_list : a FilterList object
        List of filters for photometry.
    ref : str
        Either 'kurucz' or 'hayes' to change stellar atmosphere model
        of reference star; defaults to 'kurucz'.

    Notes
    -----
    See IDL procedure utils/k_vega2ab.pro (kcorrect v4_2).
    """
    if ref == 'kurucz':
        veganame = 'lcbvega.ori'
        fname = ''.join([KCORRECT_DIR, '/data/basel/', veganame])
        lamb, flux = read_basel(fname)
        lamb = lamb * 10.
        cspeed = 2.99792e+18
        flux = np.pi * 4. * flux * cspeed / lamb**2
        flux = flux * 6.5043898e-17
    elif ref == 'hayes':
        fname = ''.join([KCORRECT_DIR, '/data/filters/',
                         'hoggraw/hayes/hayes.txt'])
        f = open(fname)
        lamb, lflux = [], []
        for each in f:
            if each.startswith('#'):
                continue
            x, y = each.split()
            lamb.append(float(x))
            lflux.append(float(y))
        lamb = np.asarray(lamb, dtype=FTYPE)
        flux = np.asarray([10**(-0.4 * np.asarray(lflux)) * 4.65e-9])
    else:
        raise ValueError('Invalid ref: must be "kurucz" (default) or "hayes".')
    lamb_e = lamb_to_edges(lamb)
    maggies = project_filters(lamb_e, flux, filter_list)
    maggies = maggies[0]
    mag = -2.5 * np.log10(maggies)
    return mag


def ab2maggie(mag, dmag=None):
    """
    Convert AB magnitudes into AB maggies.

    If the uncertainties are given as dmagnitude, the uncertainties in
    maggie are also computed.
    
    Parameters
    ----------
    mag : array_like
        Array of AB magnitudes.
    dmag : array_like, optional
        Array of uncertainties in AB magnitudes.
    """
    mag = np.asarray(mag, dtype=FTYPE)
    maggie = 10**(-0.4 * mag)

    if dmag is None:
        return maggie

    dmag = np.asarray(dmag, dtype=FTYPE)

    if dmag.shape == mag.shape:
        dmaggie = 0.4 * np.log(10.) * np.fabs(maggie * dmag)
        return maggie, dmaggie
    else:
        raise ValueError('Invalid dmag array shape.')


def maggie2ab(maggie, dmaggie=None):
    """
    Convert AB maggies into AB magnitudes.

    If the uncertainties are given as the second argument, the
    uncertainties in magnitude are also returned.

    Parameters
    ----------
    maggie : array_like
        Array of AB maggies.
    dmaggie : array_like, optional
        Array of uncertainties in AB maggies.
    """
    maggie = np.asarray(maggie).astype(FTYPE)
    mag = -2.5 * np.log10(maggie)

    if dmaggie is None:
        return mag

    dmaggie = np.asarray(dmaggie, dtype=FTYPE)

    if dmaggie.shape == maggie.shape:
        dmag = 2.5 * np.fabs(dmaggie / maggie) / np.log(10.)
        return mag, dmag
    else:
        raise ValueError('Invalid dmaggie array shape.')


def fit_nonneg(maggies, maggies_ivar, redshift, ptable,
               maxiter=50000, tolerance=1e-6):
    """
    Notes
    -----
    See IDL procedure fit/k_fit_nonneg.pro in kcorrect v4_2.
    """
    return k_fit_nonneg(maggies, maggies_ivar, redshift,
                        ptable[0], ptable[1], maxiter, tolerance, 0)


def fit_photoz(maggies, maggies_ivar, ptable, zpriors=None, lpriors=None, 
               maxiter=50000, tolerance=1e-6):
    lpriors = np.zeros(2, dtype=FTYPE) if lpriors is None else lpriors
    zpriors = np.array([0., 1000.], dtype=FTYPE) if zpriors is None else zpriors
    return k_fit_photoz(maggies,
                        maggies_ivar,
                        zpriors,
                        lpriors,
                        ptable[0], ptable[1],
                        maxiter, tolerance, 0)


def reconstruct_maggie(coeffs, redshift, ptable):
    """
    Compute maggies given the fit coefficients and redshift.

    Parameters
    ----------
    coeffs -- 
    redshift --
    ptable --

    Notes
    -----
    The definition of an AB maggie is

    m = -2.5 log10( maggie )     (1)

    where m is in AB magnitudes.  From the definition of AB magnitude,
    this leads to

              \int dlamb lamb f_lamb(lamb) R(lamb)
    maggie = --------------------------------------
             \int dlamb lamb g_lamb(R,lamb) R(lamb)

    where lamb is the wavelength, R the filter transmission curve, and
    g is the reference spectrum, i.e., for the AB scale, the reference
    spectrum is flat in fnu space.

    What this function returns, however, is

                      \int dlamb (1+z) lamb L_lamb[lamb/(1+z)] R(lamb)
    reconst. maggie = --------------------------------------------------
                           \int dlamb lamb g_lamb(R,lamb) R(lamb)

    (see Eq.(13) of Hogg et al. 2002, astro-ph/0210394v1) where z is
    the redshift, and L_lamb is the rest-frame spectrum whose absolute
    scale is determined to match the observed flux at the redshift for
    which the given coeffs were obtained via non-negative fitting.
    This definition makes the expression for the k correction (not
    band-shifted) sensible, i.e.,

                               /  reconst. maggie(z = z) \
    k correction = -2.5 log10 | ------------------------- | .
                               \  reconst. maggie(z = 0) /

    This means, however, the reconstructed maggies cannot be converted
    into *observed* AB magnitude simply via (1), as the normalization
    of the term (1+z) lamb L_lamb[lamb/(1+z)] is only correct for the
    redshift at which coeffs are obtained via the nonnegative fitting.

    See IDL procedure fit/k_reconstruct_maggies.pro in kcorrect v4_2.
    """
    return k_reconstruct_maggies(coeffs, redshift, ptable[0], ptable[1])


def abmag2flam(mag, lam, magerr=None):
    """
    Convert AB magnitude into flux in erg / s / cm^2 / angstrom.

    Parameters
    ----------
    mag : array_like
        Photometry in AB magnitude.
    lam : array_like
        Wavelength in angstrom at which the flux is measured.
    magerr : array_like, optional
        AB magnitude uncertainty.
    """
    flam = 10**( (mag + 48.60) / -2.5 ) * (+3e18 / lam**2)
    if magerr is None:
        return flam
    dflam = flam * (1. - 10**(-0.4 * magerr))
    return flam, dflam
