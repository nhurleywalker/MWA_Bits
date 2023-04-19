#!/usr/bin/env python

# Subtract one fits file from another

import sys
import numpy as np

# fits files
try:
    import astropy.io.fits as fits
except ImportError:
    import pyfits as fits

file1=sys.argv[1]
file2=sys.argv[2]
output=sys.argv[3]

print("Subtracting "+file2+" from "+file1)

hdu1=fits.open(file1, naxis=2)
hdu2=fits.open(file2, naxis=2)
data1 = np.squeeze(hdu1[0].data)
data2 = np.squeeze(hdu2[0].data)
hdu1[0].data = data1 - data2
hdu1.writeto(output,overwrite=True)
