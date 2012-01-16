#!/usr/bin/env python2.6
"""
Utility functions for file IO.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import os
import numpy as np
#import kcorrect.clib as CL
from kcorrect.globaldefs import FTYPE


def path_to_file(fname, dirs):
    """
    Get the path to a file.

    The file is searched in given directoris in order.  The path
    (directory plus file name) for the first file found in search
    directories will be returned; otherwise, None will be returned.

    Parameters
    ----------
    fname : str
        File name to search for.
    dirs : list of directories
        List of directories to search in.
    """
    for each in dirs:
        path = '/'.join([each, fname])
        if os.path.exists(path):
            return path
    return None


## DEPRECATED: USE clib.k_read_ascii_table
## def read_ascii_table(fname):
##     """Read an ASCII file in the Blanton's standard format

## DESCRIPTION

##   Read an ASCII file in the following format, which is:
 
##     <ndim> <size_{0}> ... <size_{ndim-1}>
##     <entry_0>
##     <entry_1>
##     ...
##     <entry_n>
##     ...
##     <entry_{size_0*size_1*..*size_{ndim-1}-1>
 
##   where the table element [k,j,i] (for ndim==3) would be the entry
##   element n, where n=i*size_2*size1+j*size2+k.

## REFERENCE

##   IDL procedure k_read_ascii_table.pro (kcorrect v4_1_4)
##     """
##     f = open(fname)
##     s = f.readline()
##     ts = s.split()
##     ndim = int(ts[0])
##     sizes = []
##     nelements = 1
##     for n in ts[1:]:
##         sizes.append(int(n))
##         nelements *= int(n)
##     table = np.zeros(nelements).astype(FTYPE)
##     i = 0
##     for each in f:
##         table[i] = float(each)
##         i += 1
##     f.close()
##     return np.reshape(table, sizes)


def read_basel(fname):
    """
    Read a spectrum from a Basel spectrum file.

    Parameters
    ----------
    fname : str
        Path to the Basel spectrum file.

    Notes
    -----
    See IDL procedure seds/k_read_basel.pro (kcorrect v4_2).
    """
    npoint = 1221  # hard-coded!!!
    
    f = open(fname)
    elems = f.read().split()
    lamb = np.empty(npoint, dtype=FTYPE)
    for i in range(npoint):
        lamb[i] = float(elems[i])

    flux, modelno, teff, logg, mh, vturb, xh = [], [], [], [], [], [], []
    ite = iter(elems[npoint:])
    try:
        while True:
            modelno.append( float(ite.next()) )
            teff.append( float(ite.next()) )
            logg.append( float(ite.next()) )
            mh.append( float(ite.next()) )
            vturb.append( float(ite.next()) )
            xh.append( float(ite.next()) )
            tmpflux = np.zeros(npoint, dtype=FTYPE)
            for i in range(npoint):
                tmpflux[i] = float(ite.next())
            flux.append(tmpflux)
    except StopIteration:
        pass
    return np.asarray(lamb, dtype=FTYPE), np.asarray(flux, dtype=FTYPE)
