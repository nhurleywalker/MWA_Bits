#!/usr/bin/env python

import decimal
import numpy as np

# Helper functions for formatting tables

# https://stackoverflow.com/questions/38847690/convert-float-to-string-without-scientific-notation-and-false-precision
def float_to_str(f, p):
    ctx = decimal.Context()
    ctx.prec = p
    """
    Convert the given float to a string,
    without resorting to scientific notation
    """
    d1 = ctx.create_decimal(repr(f))
    return format(d1, 'f')

def combine_cols(data, errors, uniform=True):
    '''Look at a data column and an error column, and return a LaTeX-formatted
       single column of data +/- error, to a useful number of decimal places.'''

# In the case where I want all numbers formatted in the same way: looks very nice for tables
    if uniform is True:
# First loop identifies the largest number of decimal places I need to show
        mine = 99
        for j in errors:
            if j is not None:
                e = decimal.Decimal("{0:2.1g}".format(j)).as_tuple().exponent
                if e < mine:
                    mine = e
# Second loop does the formating
        rs = []
        for i, j in zip(data, errors):
            if i is not None and j is not None:
                if mine <= 0:
                    r = ("${0:2.{1}f}\pm{2:2.{1}f}$".format(i, abs(mine), j))
                else:
                    n = decimal.Decimal("{0:2.1g}".format(i)).as_tuple().exponent
                    p = n - mine + 1
                    r = ("${0}\pm{1}$".format(float_to_str(i, p), float_to_str(j, mine)))
            else:
                r = ("--")
            rs.append(r)

# Want a different number of decimal places depending on the errors -- doesn't look nice but is more
# sensible, especially for small datasets
    else:
        rs = []
        for i, j in zip(data, errors):
            if i is not None and j is not None:
                e = decimal.Decimal("{0:2.1g}".format(j)).as_tuple().exponent
                if e <= 0:
                    r = ("${0:2.{1}f}\pm{2:2.{1}f}$".format(i, abs(e), j))
                else:
# precision on errors should always be 1
# precision on the data should be the be the exponent of the error place
# As an input to float_to_str, that is the exponent of the data minus the exponent of the error
                    n = decimal.Decimal("{0:2.1g}".format(i)).as_tuple().exponent
                    p = n - e + 1
                    r = ("${0}\pm{1}$".format(float_to_str(i, p), float_to_str(j, 1)))
            else:
                r = ("--")
            rs.append(r)

    return rs
