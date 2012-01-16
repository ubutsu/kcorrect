#!/usr/bin/env python2.6
"""
Utility functions for cosmology.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
from collections import Iterable
import numpy as np
import kcorrect.clib as CL
from kcorrect.globaldefs import COSMO_DEFAULT, FTYPE


# Deal with an input which can be a scalar float or array of floats.
# Return value will be the same as input.
def flexinput(func):
    def f(x, *args, **kwargs):
        isscalar = not isinstance(x, Iterable)
        x = np.asarray([x] if isscalar else x, dtype=FTYPE)
        y = func(x, *args, **kwargs)
        return float(y[0]) if isscalar else y.reshape(x.shape)
    f.__name__ = func.__name__
    f.__doc__ = func.__doc__
    return f


@flexinput
def ztodm(z, cosmo=COSMO_DEFAULT):
    """
    ztodm(z, cosmo=COSMO_DEFAULT)

    Compute distance modulus to the redshift.

    Parameters
    ----------
    z : array_like
        Input redshift array.
    cosmo : (Omega_matter, Omega_lambda, h)
        Input cosmology, where h = 1 is normalized to 100 km/s/Mpc.
    """
    Om, Ol, h = cosmo
    dm = np.empty(z.size, dtype=FTYPE)
    for i, val in enumerate(z.flat):
        dm[i] = CL.ztodm(val, (Om, Ol)) - 5. * np.log10(h)
    return dm


@flexinput
def ztov(z, cosmo=COSMO_DEFAULT):
    """
    ztov(z, cosmo=COSMO_DEFAULT)

    Compute total comoving volume in Mpc^3 out to the redshift.

    Parameters
    ----------
    z : array_like
        Input redshift
    cosmo : (Omega_matter, Omega_lambda, h)
        Input cosmology, where h = 1 is normalized to 100 km/s/Mpc.
    """
    Om, Ol, h = cosmo
    vol = np.empty(z.size, dtype=FTYPE)
    for i, val in enumerate(z.flat):
        vol[i] = CL.ztoV(float(val), (Om, Ol)) / h / h / h
    return vol

