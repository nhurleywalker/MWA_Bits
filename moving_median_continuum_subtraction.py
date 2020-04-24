#!/usr/bin/env python

# Perform a moving median continuum subtraction on MWA spectral line data

from astropy.io import fits
import argparse
import sys
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Input/output files")
    group1.add_argument("--fitsfile", dest='fitsfile', type=str, default=None,
                        help="Name of the input cube to be continuum-subtracted.")
    group1.add_argument("--binwidth", dest='binwidth', type=float, default=20,
                        help="Number of channels in each bin (default = 20).")
    options = parser.parse_args()

    if len(sys.argv) <= 1 or options.fitsfile is None:
        parser.print_help()
        sys.exit()

    infile = options.fitsfile
    outfile = infile.replace(".fits", "_subtracted.fits")
    contfile = infile.replace(".fits", "_moving_median.fits")
    binwidth = options.binwidth

    hdu = fits.open(infile)
    data = hdu[0].data
    nchans = data.shape[1] #0 pol, 2 Dec, 3 RA
# Replace Slice 55 with a copy of Slice 56 to avoid having to use nanmedian
    data[:,54,:,:] = data[:,55,:,:]

    med = np.zeros(hdu[0].data.shape)

    # Calculate the moving average
    # TODO: thread
    # TODO: implement step sizes >1
    # TODO: re-implement nanmedian with speedups
    for i in range(binwidth/2, nchans - binwidth/2):
        med[:,i,:,:] = np.median(data[:,i-binwidth/2:i+binwidth/2,:,:], axis=1)

    data[:,54,:,:] = np.nan*np.ones(data[:,54,:,:].shape)
    hdu[0].data = data - med
    hdu.writeto(outfile)
    hdu[0].data = med 
    hdu.writeto(contfile)
