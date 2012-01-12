#!/usr/bin/env python2.6
"""
Photometric conversion utilities
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
from kcorrect.globaldefs import FTYPE, KCORRECT_DIR
from kcorrect.projectiontable import project_filters
from kcorrect.utils.spec import lamb_to_edges
from kcorrect.utils.io import read_basel
import numpy as np


__all__ = ['vega2ab', 'ab2maggie', 'maggie2ab']


def vega2ab(filter_list, ref='kurucz'):
    """
    Return the offsets between Vega and AB magnitudes

    This function returns a conversion from Vega manitudes to AB
    magnitudes such that

    m(AB) = m(Vega) + vega2ab

    where vega2ab is the offsets returned by this function.

    INPUT

    filter_list -- A FilterList object
    ref -- Either 'kurucz' or 'hayes' to change stellar atmosphere model
           of reference star; defaults to 'kurucz'

    REFERENCE

    IDL procedure k_vega2ab.pro (kcorrect v4_1_4)
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
        raise ValueError('ref is either "kurucz" (default) or "hayes"')
    lamb_e = lamb_to_edges(lamb)
    maggies = project_filters(lamb_e, flux, filter_list)
    maggies = maggies[0]
    mag = -2.5 * np.log10(maggies)
    return mag


def ab2maggie(mag, dmag=None):
    """
    Convert AB magnitudes into AB maggies

    If the uncertainties are given as dmagnitude, the uncertainties in
    maggie are also computed.

    INPUT

    mag -- AB magnitudes
    dmag -- (Optional) uncertainties in AB magnitudes
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
    Convert AB maggies into AB magnitudes

    If the uncertainties are given as the second argument, the
    uncertainties in magnitude are also returned.

    INPUT

    maggie -- AB maggies
    dmaggie -- (Optional) uncertainties in AB maggies
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
