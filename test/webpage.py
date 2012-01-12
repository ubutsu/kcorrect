#!/usr/bin/env python2.6
"""
Testing output from the kcorrect (v4_2) web site examples.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import kcorrect as KC


with_pyplot = True
try:
    import matplotlib.pyplot as plt
except ImportError:
    with_pyplot = False
    print('# need matplotlib.pyplot for plotting a spectrum.')


def main():
    ptable = KC.ProjectionTableDB()

    mg = [1., 4.78, 10.96, 14.45, 19.05]
    mg_ivar = [1100., 28., 7.7, 4.4, 2.5]
    redshift = 0.03

    filters = KC.FilterList(['sdss_u0', 'sdss_g0', 'sdss_r0',
                             'sdss_i0', 'sdss_z0'])

    o1 = KC.KCorrect(redshift, filters, mg, mg_ivar, ptable=ptable)

    print(o1.kcorrect(band_shift=0.1))

    mg = [1., 4.73, 11.26, 14.25, 18.85]
    mg_ivar = [1100., 28., 7.7, 4.4, 2.5]

    o2 = KC.KCorrect(redshift, filters, mg, mg_ivar, ptable=ptable)

    print(o2.kcorrect(band_shift=0.1))

    x, y = o1.model_spectrum(None, 2000., 12000.)
    if with_pyplot:
        plt.plot(x, y)
        plt.show()

    print(o1.model_maggie)

    bessel = KC.FilterList(['bessell_B', 'bessell_V'])
    mags = o1.appmag(bessel)[0]
    vega2ab = KC.vega2ab(bessel)
    Bmag, Vmag = mags - vega2ab
    print(Bmag - Vmag)

    mg = [1., 4.78, 10.96, 14.45, 19.05]
    mg_ivar = [1100., 28., 7.7, 4.4, 2.5]
    filters = KC.FilterList(['sdss_u0', 'sdss_g0', 'sdss_r0',
                             'sdss_i0', 'sdss_z0'])
    d = KC.fit_photoz([mg], [mg_ivar], ptable[filters])
    print(d)


if __name__ == '__main__':
    main()
