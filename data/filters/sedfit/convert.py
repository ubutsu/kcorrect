#!/usr/bin/python
"""
Add artificially generated object(s) to an input image.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
from optparse import OptionParser
import numpy as np


HEADER = """# Units:
#  "lambda" is in Angstroms
#  "pass" is the contribution to the detector signal per photon
#         entering the atmosphere of Earth
#
# Bandpass Name(s): 
#
# Instrument: 
#
# Determined by: 
#
# Date of determination: 
#
# Meaning of/Reason for default column: Only column available
#
# Notes:
#
#"""


def main(filename, norm):
    d = np.loadtxt(filename, dtype=[('lambda', float), ('pass', float)])

    if np.any(d['pass'] > 1):
        d['pass'] /= 100.

    if norm is not None:
        d['pass'] = norm * d['pass'] / d['pass'].max()


    print(HEADER)
    print()
    print("typedef struct {")
    print(" double lambda;")
    print(" double pass;")
    print("} KFILTER;")
    print()

    for o in d:
        print('KFILTER    %(lambda).4f    %(pass).7f'
              % o)


if __name__ == '__main__':
    p = OptionParser()

    opts, args = p.parse_args()

    filename = args[0]
    norm = float(args[1]) if len(args) == 2 else None

    main(filename, norm)
