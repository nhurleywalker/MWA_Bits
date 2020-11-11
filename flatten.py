#!/usr/bin/env python

# Script to flatten 4D images

from astropy.io import fits
import sys
import numpy as np

__author__ = "Natasha Hurley-Walker"

if __name__ == '__main__':
    infile = sys.argv[-1:][0]
    image4D = fits.open(infile)
    if len(image4D[0].data.shape)==2:
        print(infile+" already has just two dimensions. Leaving unchanged.")
    else:
        image4D[0].data = np.squeeze(image4D[0].data)
        print("Replacing "+infile+" with flattened version.")
        try:
            image4D[0].header.remove('DATAMIN')
            image4D[0].header.remove('DATAMAX')
            image4D[0].header['NAXIS'] = 2
        except:
            pass
        image4D.writeto(infile,overwrite=True)

    
