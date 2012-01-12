#!/usr/bin/env python2.6
"""
Utility functions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import os
import numpy as np
import kcorrect.clib as CL
from kcorrect.globaldefs import COSMO_DEFAULT, FTYPE


##############################################################################
# file io utilities


def path_to_file(fname, dirs):
    """
    Get the path to a file 

    The file is searched in given directoris in order.  If the file is
    not found, None will be returned; otherwise, the path (i.e.,
    directory plus file name) will be returned.

    INPUT

    fname -- File name
    dirs -- List of directories
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
    Read a spectrum from a Basel spectrum file

    INPUT

    fname -- Path to the Basel spectrum file

    REFERENCE

    IDL procedure k_read_basel.pro (kcorrect v4_1_4)
    """
    f = open(fname)
    elems = f.read().split()
    lamb = np.zeros(1221, dtype=FTYPE)
    for i in xrange(0, 1221):
        lamb[i] = float(elems[i])

    flux, modelno, teff, logg, mh, vturb, xh = [], [], [], [], [], [], []
    ite = iter(elems[1221:])
    try:
        while True:
            modelno.append( float(ite.next()) )
            teff.append( float(ite.next()) )
            logg.append( float(ite.next()) )
            mh.append( float(ite.next()) )
            vturb.append( float(ite.next()) )
            xh.append( float(ite.next()) )
            tmpflux = np.zeros(1221, dtype=FTYPE)
            for i in xrange(1221):
                tmpflux[i] = float(ite.next())
            flux.append(tmpflux)
    except StopIteration:
        pass
    return (np.asarray(lamb, dtype=FTYPE), np.asarray(flux, dtype=FTYPE))


if __name__ == '__main__':
    area = 0.85 # square degrees
    area = area * (np.pi/180.)**2

    print (comvol(1.0) - comvol(0.8)) * (area / (4.*np.pi))

    print( distmod(0.4))
    print( distmod([[0.2,0.4,0.6],[0.1,0.3,0.5]]))
    print( distmod(0.0))
