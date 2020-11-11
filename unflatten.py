#!/usr/bin/env python

# Script to unflatten flattened images

from astropy.io import fits
import sys
import numpy as np

__author__ = "Natasha Hurley-Walker"

if __name__ == '__main__':
    infile = sys.argv[-1:][0]
    image2D = fits.open(infile)
    if len(image2D[0].data.shape)>2:
        print(infile+" already has more than two dimensions. Leaving unchanged.")
    else:
        image2D[0].data = np.reshape(image2D[0].data,(1,1,image2D[0].data.shape[0],image2D[0].data.shape[1]))
        print("Replacing "+infile+" with unflattened version.")
        try:
            image2D[0].header.remove('DATAMIN')
            image2D[0].header.remove('DATAMAX')
        except:
            pass
        image2D.writeto(infile,overwrite=True)

    
