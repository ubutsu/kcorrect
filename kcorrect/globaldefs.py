#!/usr/bin/env python2.6
"""
Global definitions
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import os
from numpy import float32


# make sure KCORRECT_DIR is set
KCORRECT_DIR = os.getenv('KCORRECT_DIR')
if KCORRECT_DIR is None:
    print('KCORRECT_DIR not defined.')
    exit()
elif not os.path.exists(KCORRECT_DIR):
    print('Problem accessing %s.' % KCORRECT_DIR)
    exit()


# default directories for some data sets
DEFAULT_FILTER_DIR = KCORRECT_DIR + '/data/filters'
DEFAULT_TEMPLATE_DIR = KCORRECT_DIR + '/data/templates'

# default cosmology (Omega_matter, Omega_lambda, h_100)
COSMO_DEFAULT = (0.3, 0.7, 0.7)

# default redshift range (zmin, zmax, # of points)
#ZRANGE_DEFAULT = (0.0, 2.0, 1000)
ZRANGE_DEFAULT = (0.0, 6.0, 3000)

# template name
VNAME = 'default'

# within the kcorrect C library, float C type is consistently used.
FTYPE = float32

# filter definition file extension
FILTEREXT = '.par'
