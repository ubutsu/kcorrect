#!/usr/bin/env python2.6
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.filter import FilterList
from kcorrect.globaldefs import COSMO_DEFAULT, FTYPE
from kcorrect.projectiontable import PTABLE_MASTER


class Photo(object):
    """
    Base class for basic photometry

    AB maggies are handy to deal with internally, but they are not a
    conventional photometric system observers are used to.  One wants
    to subclass Photo to implement a basic set of photometric
    conversion from the input photometry to AB maggies by overriding
    set_data and convert2maggie methods.

    This class also implements other essential photometry routines.
    """

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        ptable : ProjectionTableDB object
            Defaults to ProjectionTableDB().  If a non-standard
            projection table (i.e., with different redshift intervals
            or template spectra, etc.) is to be used in reconstructing
            AB maggies, then a custom ProjectionTableDB object can be
            specified here.
        cosmo : (Omega_matter, Omega_Lambda, h_100)
        """
        self.ptable = kwargs['ptable'] if 'ptable' in kwargs else PTABLE_MASTER
        self.cosmo = kwargs['cosmo'] if 'cosmo' in kwargs else COSMO_DEFAULT
        self.set_input(*args, **kwargs)

    def __len__(self):
        return self.input_maggie.shape[0]

    def set_input(self, filter_list, maggie, maggie_ivar, *args, **kwargs):
        self.set_input_photo(filter_list, maggie, maggie_ivar, *args, **kwargs)

    def set_input_photo(self, filter_list, maggie, maggie_ivar,
                        *args, **kwargs):
        """
        Set input data

        Parameters
        ----------
        filter_list : A FilterList object
        maggie : List of maggies
        maggie_ivar : List of maggies inverse variances
        """
        if not isinstance(filter_list, FilterList):
            raise TypeError('Expecting a FilterList object for filter_list.')

        # ensure the correct shape for 2D photometry data arrays
        maggie = np.atleast_2d(np.asarray(maggie, dtype=FTYPE))
        maggie_ivar = np.atleast_2d(np.asarray(maggie_ivar, dtype=FTYPE))
        if maggie.shape != maggie_ivar.shape:
            s = ('maggie and maggie_ivar must have the same array shape.')
            raise ValueError(s)
        if maggie.shape[1] != len(filter_list):
            s = ('maggie must be a 2d array of size n * m,'
                 ' where n = # of objects and m = # of filters.')
            raise ValueError(s)

        # input photometry are already in units of maggies
        self.filter_list = filter_list
        self.input_maggie = maggie
        self.input_maggie_ivar = maggie_ivar

    # def convert_to_maggie(self, photo, photo_ivar):
    #     """
    #     Convert input photometry to AB maggies and their inverse
    #     variances

    #     This function converts input photometry, which are often in
    #     magnitudes, into AB maggies internally used.  Override this
    #     method within a subclass if input photometry data are not
    #     already in AB maggies.

    #     The method returns AB maggies and their inverse variances.

    #     Parameters
    #     ----------
    #     photo, photo_ivar -- Input photometry and inverse variances.
    #     """
    #     # simply return maggies in this base class
    #     return photo, photo_ivar
