#!/usr/bin/env python2.6
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.filter import FilterList
from kcorrect.globaldefs import COSMO_DEFAULT, FTYPE
from kcorrect.projectiontable import PTABLE_MASTER
from kcorrect.utils.cosmology import ztodm
from kcorrect.utils.photo import maggie2ab, reconstruct_maggie


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

    def appmag(self, filter_list=None, redshift=None):
        """
        Compute apparent AB magnitudes from best fit SED models.

        Computes the apparent magnitudes through a set of filters
        convolving the best fit SED models.  Optionally, the objects
        to be observed can be moved to different redshifts than the
        input redshifts; i.e., the best fit SEDs are redshifted
        arbitrarily and *observed* with the set of filters.

        Parameters
        ----------
        filter_list : A FilterList object, optional
            A set of filters.
        redshift : array_like, optional
            A redshift or a 1D array of redshift (with the number of
            elements equal to the number of input objects) to observe
            galaxies at arbitrary redshifts
        """
        filter_list = self.filter_list if filter_list is None else filter_list
        
        if redshift is None:
            # simply observe the input objects with given filters
            maggie = reconstruct_maggie(self.coeffs, self.redshift,
                                        self.ptable[filter_list])
            return maggie2ab(maggie)

        # when redshifts are given, use them
        redshift = np.atleast_1d(redshift).astype(FTYPE)
        if redshift.size == 1:
            redshift = np.ones(self.redshift.shape, dtype=FTYPE) * redshift[0]
        else:
            redshift = np.reshape(redshift, self.redshift.shape)

        maggie_S = reconstruct_maggie(self.coeffs, redshift,
                                      self.ptable[filter_list])
        m_S = maggie2ab(maggie_S)
        dm1 = ztodm(self.redshift, self.cosmo)
        dm2 = ztodm(redshift, self.cosmo)

        # because of what reconstruct_maggie returns, the flux needs
        # to be rescaled by the ratio of distances
        return m_S - dm1 + dm2

    def model_spectrum(self, index=None, lmin=None, lmax=None):
        """
        Returns the best fit SED(s) in the rest frame.

        The output SED is in units of erg s^-1 cm^-2 A^-1.

        Parameters
        ----------
        index : array of indices
            The index(ices) of the object for which spectrum is obtained.
        lmin : float
            Minimum wavelength in angstrom.
        lmax : float
            Maximum wavelength in angstrom.
        """
        lamb, vmatrix = self.ptable[self.filter_list].templates
        lamb = 0.5 * (lamb[:-1] + lamb[1:])
        vmatrix = np.transpose(vmatrix)

        if (lmin is not None) or (lmax is not None):
            mask = np.ones(lamb.shape, dtype=bool)
            mask = mask if lmin is None else mask * np.greater(lamb, lmin)
            mask = mask if lmax is None else mask * np.less(lamb, lmax)
            lamb = lamb[mask]
            vmatrix = np.compress(mask, vmatrix, axis=0)

        index = (np.arange(0, len(self)) if index is None
                 else np.atleast_1d(np.asarray(index, dtype=int)))

        data = []
        for coeffs in self.coeffs[index]:
            flux = np.add.reduce(np.transpose(vmatrix * coeffs))
            data.append(flux)

        return (lamb, np.asarray(data) if len(data) > 1 else data[0])

