#!/usr/bin/env python2.6
"""
Photometric redshift routine
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.photo import Photo
from kcorrect.utils.photo import fit_photoz


MAXITER = 3000


class PhotoZ(Photo):

    #def __init__(self, *args, **kwargs):
    #    super(PhotoZ, self).__init__(*args, **kwargs)

    def set_input(self, filter_list, maggie, maggie_ivar, *args, **kwargs):
        self.maxiter = kwargs['maxiter'] if 'maxiter' in kwargs else MAXITER
        self.lpriors = kwargs['lpriors'] if 'lpriors' in kwargs else None
        self.zpriors = kwargs['zpriors'] if 'zpriors' in kwargs else None

        super(PhotoZ, self).set_input_photo(filter_list, maggie, maggie_ivar)

        self.fit_photoz(self.maxiter)

    def fit_photoz(self, maxiter=MAXITER):
        r = fit_photoz(self.input_maggie, self.input_maggie_ivar,
                       self.ptable[self.filter_list], 
                       lpriors=self.lpriors, zpriors=self.zpriors,
                       maxiter=maxiter)
        self.redshift, self.coeffs, self.chi2, self.iter = r
