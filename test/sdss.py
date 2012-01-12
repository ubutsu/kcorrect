#!/usr/bin/env python2.6
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
import matplotlib.pyplot as plt
import kcorrect as KC
from kcorrect.sdss import SDSSKCorrect, SDSSPhotoZ


def main():
    d = np.loadtxt('test/sdss.dat', dtype=float)
    z = d[:, 0]
    lup = d[:, 2:7]
    dlup = d[:, 7:12]
    red = d[:, 12:17]

    ####
    kc = SDSSPhotoZ(lup, dlup, red)
    print(kc.redshift)
    print(z)
    plt.plot(z, kc.redshift, 'o')
    plt.show()
    ####

    kc = SDSSKCorrect(z, lup, dlup, red, cosmo=(0.3, 0.7, 1.0))

    print('### k-corrections: SDSS ugriz ###')
    kcorr = kc.kcorrect()
    for each in kcorr: print('%+.6f %+.6f %+.6f %+.6f %+.6f' % tuple(each))
    print()

    print('### absolute magnitudes: SDSS ugriz ###')
    absmag = kc.absmag()
    for each in absmag: print('%+.4f %+.4f %+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### k-corrections: SDSS ugriz band-shifted to z=0.1 ###')
    band_shift = 0.1
    kcorr = kc.kcorrect(band_shift=band_shift)
    for each in kcorr: print('%+.6f %+.6f %+.6f %+.6f %+.6f' % tuple(each))
    print()

    print('### absolute magnitudes: SDSS ugriz band-shifted to z=0.1 ###')
    absmag = kc.absmag(band_shift=band_shift)
    for each in absmag: print('%+.4f %+.4f %+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### absolute magnitudes: Bessell UBVRI ###')
    fl = KC.FilterList(['bessell_U','bessell_B','bessell_V',
                        'bessell_R','bessell_I'])
    absmag = kc.absmag(fl)
    for each in absmag: print('%+.4f %+.4f %+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### absolute magnitudes: Bessell UBVRI band-shifted to z=0.1 ###')
    absmag = kc.absmag(fl,band_shift)
    for each in absmag: print('%+.4f %+.4f %+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### apparent magnitudes: DEEP BRI ###')
    fl = KC.FilterList(['deep_B','deep_R','deep_I'])
    fl.plot()
    mag = kc.appmag(fl)
    for each in mag: print('%+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### apparent magnitudes: DEEP BRI on Vega scale ###')
    vega2ab = KC.vega2ab(fl)
    mag = kc.appmag(fl) - vega2ab
    for each in mag: print('%+.4f %+.4f %+.4f' % tuple(each))
    print()

    print('### apparent magnitudes: DEEP BRI if at z = 1 ###')
    mag = kc.appmag(fl,redshift=1.0)
    for each in mag: print('%+.4f %+.4f %+.4f' % tuple(each))
    print()


if __name__ == '__main__':
    main()
