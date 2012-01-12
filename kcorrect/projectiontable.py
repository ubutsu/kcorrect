#!/usr/bin/env python2.6
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
import numpy as np
from kcorrect.globaldefs import DEFAULT_TEMPLATE_DIR, FTYPE, ZRANGE_DEFAULT
from kcorrect.utils.io import path_to_file
from kcorrect.filter import FilterList
from kcorrect.clib import k_read_ascii_table, k_projection_table


def load_vmatrix(vname='default', file_paths=None):
    """
    Load spectral template information.

    Notes
    -----
    This function loads so-called vmatrix in the IDL version.
    
    See IDL procedure utils/k_load_vmatrix.pro (kcorrect v4_2).
    """
    vfile = 'vmatrix.%s.dat' % vname
    lfile = 'lambda.%s.dat' % vname
    file_paths = [] if file_paths is None else file_paths
    file_paths = file_paths + [DEFAULT_TEMPLATE_DIR]
    vmatrix = k_read_ascii_table(path_to_file(vfile, file_paths))
    lamb = k_read_ascii_table(path_to_file(lfile, file_paths))
    return lamb, vmatrix


def make_projection_table(templates, filter_list, zrange=ZRANGE_DEFAULT):
    """
    Make a table of redshift-dependent projection.

    Returns a tuple of redshift values and corresponding projection
    table.

    Parameters
    ----------
    templates : (lamb_edges, flux)
        A spectrum to be convolved with filters.
    filter_list : a FilterList object
        Filters to convolve with.
    zrange : (z_min, z_max, nz)
        Tuple of z_min, z_max, and a number of redshifts between them.

    Notes
    -----
    This function generates so-called rmatrix in the IDL version.

    See IDL procedure fit/k_projection_table.pro (kcorrect v4_2).
    """
    zmin, zmax, nz = zrange
    zvals = zmin + (zmax - zmin) * (np.arange(nz, dtype=FTYPE)
                                    + 0.5) / (1. * nz)
    filter_curves = filter_list.load_filters()
    band_shift = 0.
    # k_projection_table.pro returns (zvals, rmatrix)
    return zvals, k_projection_table(filter_curves,
                                     templates[0], templates[1],
                                     zvals, band_shift)


class ProjectionTable(tuple):
    """
    An instance of this class keeps the projection table for a set of
    filters and redshifts.
    """

    def __new__(cls, filter_list, zrange=ZRANGE_DEFAULT, vname='default',
                file_paths=None):
        templates = load_vmatrix(vname, file_paths)
        pt = make_projection_table(templates, filter_list, zrange)
        return tuple.__new__(cls, pt)

    def __init__(self, filter_list, zrange=ZRANGE_DEFAULT,
                 vname='default', file_paths=None):
        self.filter_list = filter_list
        self.zrange = zrange
        self.vname = vname
        self.file_paths = file_paths
        self.templates = load_vmatrix(vname, file_paths)


class ProjectionTableDB(dict):
    """
    A collection of ProjectionTable instances.
    """
    
    def __init__(self, zrange=ZRANGE_DEFAULT, vname='default',
                 file_paths=None):
        self.zrange = zrange
        self.vname = vname
        self.file_paths = file_paths
        dict.__init__(self, {})

    def __getitem__(self, key):
        # check the key format
        filter_list = key
        if not isinstance(filter_list, FilterList):
            raise KeyError('Key must be an FilterList instance.')
        # return if the projection table exists for the FilterList
        for fl in self.keys():
            if fl == filter_list:
                return dict.__getitem__(self, fl)
        # projection table does not exist for the FilterList, so make
        # a new one
        self[filter_list] = ProjectionTable(filter_list,
                                            self.zrange,
                                            self.vname,
                                            self.file_paths)
        return dict.__getitem__(self, filter_list)


# global master projection table database
PTABLE_MASTER = ProjectionTableDB()


def test():
    filter_list = FilterList(['sdss_u0.par','sdss_g0.par','sdss_r0.par',
                              'sdss_i0.par','sdss_z0.par'])
    filter_list2 = FilterList(['sdss_u0.par','sdss_g0.par','sdss_r0.par',
                               'sdss_i0.par'])
    zrange = (0.,2.,2000)
    pt = ProjectionTable(filter_list,zrange)
    print(pt.vname)
    print(pt)

    ptdb = ProjectionTableDB()
    a = ptdb[filter_list]
    print('done 1')
    b = ptdb[filter_list]
    print('done 2')
    c = ptdb[filter_list2]
    print('done 3')


if __name__ == '__main__':
    test()
