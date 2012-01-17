#!/usr/bin/env python2.6
import glob
import os
import sys


# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

from distutils.core import setup, Extension


version = '20111222'
__version__ = version


if (not hasattr(sys, 'version_info')
    or sys.version_info < (2, 6, 2, 'final', 0)):
    raise SystemExit('Python 2.6.2 or later required to build kcorrect.')


# def add_user_options():
#     """Add user options.

#     Add a command line option --local=<install-dir> which is an
#     abbreviation for 'put all of pyfits in <install-dir>/pyfits'.
#     """
#     if "--help" in sys.argv:
#         print 
#         print ' options:'
#         print '--local=<install-dir> same as --install-lib=<install-dir>'
#     for a in sys.argv:
#         if a.startswith('--local='):
#             instdir = a.split('=')[1]
#             sys.argv.extend(['--install-lib='+instdir,])
#             sys.argv.remove(a)


def main():
    #add_user_options()

    # NumPy path.
    try:
        import numpy
    except ImportError:
        raise SystemExit('numpy needs to be installed.')
    numpy_include = numpy.__path__[0] +'/core/include'

    # C library source codes.
    srcs = ['kcorrect/kcorrect.i',
            'src/gaussj.c',
            'src/iterate_lf.c',
            'src/k_binspec.c',
            'src/k_brent.c',
            'src/k_choldc.c',
            'src/k_cholsl.c',
            'src/k_evolve.c',
            'src/k_fileopen.c',
            'src/k_filter_struct.c',
            'src/k_fit_nonneg.c',
            'src/k_fit_photoz.c',
            'src/k_fit_spec.c',
            'src/k_fit_spec_linear.c',
            'src/k_interpolate.c',
            'src/k_load_filters.c',
            'src/k_locate.c',
            'src/k_midpnt.c',
            'src/k_nonneg_solve.c',
            'src/k_polint.c',
            'src/k_projection_table.c',
            'src/k_qromo.c',
            'src/k_read_ascii_table.c',
            'src/k_reconstruct_maggies.c',
            'src/k_strparse.c',
            'src/k_utils.c',
            'src/k_yanny_readone.c',
            'src/k_zbrent.c',
            'src/lf_WH_interp.c',
            'src/lf_calc_vmax.c',
            'src/lf_eep.c',
            'src/lf_eepfit.c',
            'src/lf_select_eep.c',
            'src/lf_set_AB.c',
            'src/lf_sum_AB.c',
            'src/phierrors_lf.c',
            'src/philike.c',
            'src/ztransform.c']

    # Data files.
    filters = glob.glob('data/filters/*.par')
    templates = glob.glob('data/templates/*.dat')

    # Install data files under the package directory.
    path_data = ('lib/python%d.%d/site-packages/kcorrect/data'
                 % sys.version_info[0:2])

    setup(name = 'kcorrect',
          version = version,
          description = 'Python implemenation of kcorrect by Michael Blanton',
          author = 'Taro Sato',
          author_email = 'ubutsu@gmail.com',
          maintainer = 'Taro Sato',
          maintainer_email = 'ubutsu@gmail.com',
          url = 'http://www.physics.ucsb.edu/~taro/kcorrect/',
          download_url = 'http://www.physics.ucsb.edu/~taro/kcorrect/download',
          license = 'http://cosmo.nyu.edu/mb144/kcorrect/',
          platforms = ['Linux'],
          packages = ['kcorrect',
                      'kcorrect.utils'],
          package_dir = {'kcorrect': 'kcorrect'},
          data_files = [(path_data + '/filters', filters),
                        (path_data + '/templates', templates),
                        (path_data + '/filters/hoggraw/hayes',
                         ['data/filters/hoggraw/hayes/hayes.txt']),
                        (path_data + '/basel',
                         ['data/basel/lcbvega.ori'])],
          ext_modules = [Extension('kcorrect._clib', srcs,
                                   include_dirs=[numpy_include,
                                                 'include'])])
          #include_package_data=True)


if __name__ == "__main__":
    main()
