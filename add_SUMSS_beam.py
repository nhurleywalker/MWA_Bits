#!/usr/bin/env

import sys
from astropy.io import fits
from numpy import sin, radians

hdu = fits.open(sys.argv[1])

Dec = hdu[0].header["CRVAL2"]

bpa = 0.0
bmin = 45./3600.

bmaj = bmin / sin(abs(radians(Dec)))

hdu[0].header["BMAJ"] = bmaj
hdu[0].header["BMIN"] = bmin
hdu[0].header["BPA"] = bpa

hdu.writeto(sys.argv[2])

