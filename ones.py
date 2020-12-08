#!/usr/bin/env python

# Zero a fitsfile

import sys
import numpy as np
import astropy.io.fits as fits

file1=sys.argv[1]
output=sys.argv[2]
val=sys.argv[3]

print("Blanking "+file1)

hdu1=fits.open(file1)
hdu1[0].data = val*np.ones(hdu1[0].data.shape)
hdu1.writeto(output,overwrite=True)
