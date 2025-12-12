#!/usr/bin/env python

# Put NaNs around the edges of a FITS file

import sys
import numpy as np
import astropy.io.fits as fits

file1 = sys.argv[1]

b = int(sys.argv[2])

output = sys.argv[3]

print("NaN-ing edges of "+file1)

hdu1 = fits.open(file1)
dat = hdu1[0].data
print(dat.shape)
if dat.shape[0] == 1:
    l = dat.shape[3]
    dat[:,:,0:b,:] = np.nan
    dat[:,:,:,0:b] = np.nan
    dat[:,:,l-b:,:] = np.nan
    dat[:,:,:,l-b:] = np.nan
else:
    l = dat.shape[0]
    dat[0:b,:] = np.nan
    dat[:,0:b] = np.nan
    dat[l-b:,:] = np.nan
    dat[:,l-b:] = np.nan
    
hdu1.writeto(output,overwrite=True)
