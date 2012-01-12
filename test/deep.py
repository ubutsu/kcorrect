#!/usr/bin/env python2.6
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
import kcorrect as KC
from kcorrect.deep import DEEPKCorrect, DEEPPhotoZ


def main():
    d = np.loadtxt('test/deep.dat',
                   dtype=[('z', float), ('dz', float),
                          ('B', float), ('R', float), ('I', float)])
    z = d['z']
    mag = np.transpose([d['B'], d['R'], d['I']])

    # testing photo-z
    pz = DEEPPhotoZ(mag)
    print(pz.redshift)
    print(z)
    exit()
    # photo-z end


    kc = DEEPKCorrect(z, mag, cosmo=(0.3, 0.7, 1.0))

    print('### k-corrections: DEEP BRI ###')
    kcorr = kc.kcorrect()
    for each in kcorr:
        print('%+.6f %+.6f %+.6f' % tuple(each))
    print()

    print('### absolute magnitudes: DEEP BRI ###')
    absmag = kc.absmag()
    for each in absmag:
        print('%+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### k-corrections: DEEP BRI band-shifted to z=1 ###')
    band_shift = 1.0
    kcorr = kc.kcorrect(band_shift=band_shift)
    for each in kcorr:
        print('%+.6f %+.6f %+.6f' % tuple(each))
    print()

    print('### absolute magnitudes: DEEP BRI band-shifted to z=1 ###')
    absmag = kc.absmag(band_shift=band_shift)
    for each in absmag:
        print('%+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### k-corrections: NUB (IDL default) ###')
    fl = KC.FilterList(['galex_NUV', 'bessell_U', 'bessell_B'])
    kcorr = kc.kcorrect(fl)
    for each in kcorr:
        print('%+.6f %+.6f %+.6f' % tuple(each))
    print()

    print('### absolute magnitudes: NUB (IDL default) ###')
    absmag = kc.absmag(fl)
    for each in absmag:
        print('%+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### absolute magnitudes: Bessell UBVRI on Vega scale ###')
    fl = KC.FilterList(['bessell_U', 'bessell_B', 'bessell_V',
                        'bessell_R', 'bessell_I'])
    vega2ab = KC.vega2ab(fl)
    absmag = kc.absmag(fl) - vega2ab
    for each in absmag:
        print('%+.4f %+.4f %+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### apparent magnitude SDSS ugriz if at z = 0.1 ###')
    fl = KC.FilterList(['sdss_u0', 'sdss_g0', 'sdss_r0',
                        'sdss_i0', 'sdss_z0'])
    mag = kc.appmag(fl,redshift=0.1)
    for each in mag:
        print('%+.4f %+.4f %+.4f %+.4f %+.4f' % tuple(each))
    print()


if __name__ == '__main__':
    main()
