#!/usr/bin/env python2.6
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.photo import Photo
from kcorrect.utils.cosmology import ztodm
from kcorrect.utils.photo import maggie2ab, fit_nonneg, reconstruct_maggie
from kcorrect.globaldefs import FTYPE


__all__ = ['KCorrect', 'KCorrectAB']


class KCorrect(Photo):
    """
    k correction class for photometry in maggies

    This is the base class for k correction which implements the most
    basic functionalities.  AB maggies are handy to deal with
    internally, but they are not really the conventional photometric
    system observers are used to.  Therefore one usually wants to
    subclass KCorrect to implement a conversion from the input
    photometry to AB maggies by overriding set_data and convert2maggie
    methods.
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
        maxiter : int
            Maximum number of iterations for non-negative fitting.
        verbosity : int
            Set the level of verbosity for output; currently not
            really used.
        """
        self.maxiter = kwargs['maxiter'] if 'maxiter' in kwargs else 3000
        super(KCorrect, self).__init__(*args, **kwargs)
        self._model_maggie = None
        self._distmod = None

    def set_input(self, redshift, filter_list, maggie, maggie_ivar,
                  *args, **kwargs):
        """
        Set input data

        INPUT

        redshift -- List of redshifts
        filter_list -- A FilterList object
        maggie -- List of AB maggies
        maggie_ivar -- List of inverse variances of AB maggies
        """
        super(KCorrect, self).set_input_photo(filter_list, maggie, maggie_ivar)

        redshift = np.atleast_1d(np.asarray(redshift, dtype=FTYPE))
        self.redshift = redshift.reshape((len(self), 1))

        # zero redshifts for convenience
        self.redshift0 = np.zeros(self.redshift.shape, dtype=FTYPE)

        # get fit coefficients once and for all
        self.do_fit(maxiter=self.maxiter)

    @property
    def model_maggie(self):
        """
        Maggies of the model SED through the input passbands.
        
        model_maggie only differs from input_maggie in that they are
        computed from the best fit SEDs through the input passbands;
        these two maggies should be very close to each other if the
        fits are good.
        """
        if self._model_maggie is not None:
            return self._model_maggie
        self._model_maggie = reconstruct_maggie(self.coeffs, self.redshift,
                                                self.ptable[self.filter_list])
        return self._model_maggie

    @property
    def distmod(self):
        """
        Distance moduli to input redshifts.
        """
        if self._distmod is not None:
            return self._distmod
        dm = ztodm(self.redshift, self.cosmo)
        self._distmod = np.asarray(dm, dtype=FTYPE)
        return self._distmod

    def do_fit(self, maxiter=3000):
        """
        Do non-negative fitting.
        """
        self._model_maggie = None
        r = fit_nonneg(self.input_maggie, self.input_maggie_ivar,
                       self.redshift, self.ptable[self.filter_list],
                       maxiter)
        self.coeffs, self.chi2, self.iter = r

    def kcorrect(self, filter_list_q=None, filter_list_r=None, band_shift=0.):
        """
        Compute the k-corrections between two sets of filters

        This method computes k corrections bewteen two sets of filters
        such that

        m_r = M_q + DM(z) + K_qr(z)

        where m_r is the apparent magnutude through filter r, M_q the
        absolute magnitude through filter q, DM the distance modulus,
        and K_qr is the k correction between filters q and r.

        The number of filters in sets q and r must match.

        INPUT

        filter_list_q, filter_list_r
            -- FilterList objects; if not provided explicitly, q and r
               both defaults to self.filter_list

        band_shift -- Amount of redshift used for computing the
                      absolute magnitudes through band-shifted
                      filters; defaults to zero.
        """
        if not filter_list_q: filter_list_q = self.filter_list
        if not filter_list_r: filter_list_r = self.filter_list

        if len(filter_list_q) != len(filter_list_r):
            raise ValueError('Numbers of filters must match.')

        if filter_list_r != self.filter_list:
            rm2 = reconstruct_maggie(self.coeffs, self.redshift,
                                     self.ptable[filter_list_r])
        else:
            rm2 = self.model_maggie

        if band_shift > 0.:
            bs = np.ones(self.redshift.shape, dtype=FTYPE) * band_shift
            rm1 = reconstruct_maggie(self.coeffs, bs,
                                     self.ptable[filter_list_q])
            rm1 = rm1 / (1. + bs)
        else:
            rm1 = reconstruct_maggie(self.coeffs, self.redshift0,
                                     self.ptable[filter_list_q])
        return 2.5 * np.log10(rm1 / rm2)

    def absmag(self, filter_list=None, band_shift=0.):
        """
        Compute the absolute AB magnitude through a set of filters

        This method returns the absolute AB magnitudes through a set
        of (band-shifted) filters, which defaults to the input
        filters.  The output magnitudes are computed from the input
        photometries where available; otherwise the output magnitudes
        are computed from the model photometries reconstructed from
        the best fit SED.

        INPUT

        filter_list -- A FilterList object; defaults to the input filters
        band_shift -- Amount of redshift used for computing the
                      absolute magnitudes through band-shifted
                      filters; defaults to zero.
        """
        if not filter_list:
            filter_list = self.filter_list

        # get model maggies for abs mag
        bs = np.ones(self.redshift.shape, dtype=FTYPE) * band_shift
        rm = reconstruct_maggie(self.coeffs, bs, self.ptable[filter_list])
        rm = rm / (1. + bs)

        # fill all entries with model
        mag = maggie2ab(rm) - self.distmod
        if len(filter_list) == len(self.filter_list):
            # now replace those entries with good input and well
            # defined ivar
            mask = (np.greater(self.input_maggie,0.) *
                    np.greater(self.input_maggie_ivar, 0.)).astype(int)
            kcorr = self.kcorrect(filter_list_q=filter_list,
                                  band_shift=band_shift)
            mag = ((1 - mask) * mag
                   + mask * (maggie2ab(self.input_maggie)
                             - self.distmod - kcorr))
        return mag

    def appmag(self, filter_list, redshift=[]):
        """
        Compute apparent AB magnitudes through a set of filters

        Computes the apparent magnitudes through a set of filters.
        Optionally, the objects to be observed can be moved to
        different redshifts than the input redshifts; i.e., the best
        fit SEDs are redshifted arbitrarily and *observed* with the
        set of filters.

        INPUT

        filter_list -- A FilterList object
        redshift -- A redshift or a 1D array of redshift (with the
                    number of elements equal to the number of input
                    objects) to observe galaxies at arbitrary
                    redshifts
        """
        redshift = np.atleast_1d(redshift).astype(FTYPE)
        if redshift.size == 0:
            # simply observe the input objects with given filters
            maggie = reconstruct_maggie(self.coeffs, self.redshift,
                                        self.ptable[filter_list])
            return maggie2ab(maggie)

        # when redshifts are given, use them
        if redshift.size == 1:
            redshift = np.ones(self.redshift.shape, dtype=FTYPE) * redshift[0]
        elif redshift.size == self.redshift.size:
            redshift = np.reshape(redshift, self.redshift.shape)
        else:
            raise ValueError('Inconsistent redshift array')
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
            Minimum wavelength.
        lmax : float
            Maximum wavelength.
        """
        index = (np.arange(0, len(self)) if index is None
                 else np.atleast_1d(np.asarray(index, int)))

        lamb, vmatrix = self.ptable[self.filter_list].templates
        lamb = 0.5 * (lamb[:-1] + lamb[1:])

        mask = np.ones(lamb.shape, dtype=bool)
        if lmin is not None:
            mask = mask * np.greater(lamb, lmin)
        if lmax is not None:
            mask = mask * np.less(lamb, lmax)
        lamb = lamb[mask]

        vmatrix = np.compress(mask, np.transpose(vmatrix), axis=0)

        data = []
        for coeffs in self.coeffs[index]:
            flux = np.add.reduce(np.transpose(vmatrix * coeffs))
            data.append(flux)
        return (lamb, data if len(data) > 1 else data[0])


class KCorrectAB(KCorrect):
    """
    K-correction for AB photometry

    Expect AB magnitudes and 1-sigma uncertainties
    """

    def set_input(self, filter_list, abmag, abmag_sigma, redshift,
                  *args, **kwargs):
        """
        Set input data

        If abmag_sigma is not given, the minimum uncertainties of 0.02
        mag will be assumed.

        INPUT

        redshift -- List of redshifts
        filter_list -- FilterList object
        abmag -- List of AB magnitudes
        abmag_sigma -- List of 1-sigma uncertainties in AB magnitudes;
                       defaults to an empty list
        """

        # if uncertainties not given, assume 2% errors
        abmag_sigma = (np.ones(abmag.shape, dtype=FTYPE) * 0.02
                       if abmag_sigma is None else abmag_sigma)

        # KCorrect.set_data expects maggies
        maggie, maggie_sigma = ab2maggie(abmag, abmag_sigma)
        maggie, maggie_ivar = maggie, sigma2ivar(maggie_sigma)

        super(KCorrect, self).set_input(filter_list, maggie, maggie_ivar,
                                        redshift, *args, **kwargs)
