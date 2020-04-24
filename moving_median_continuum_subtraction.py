#!/usr/bin/env python

# Perform a moving median continuum subtraction on MWA spectral line data

from astropy.io import fits
import argparse
import sys
import numpy as np

__author__ = "Natasha Hurley-Walker"
__date__ = "2020-04-24"

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
    med = np.zeros(hdu[0].data.shape)

    # Calculate the moving average
    for i in range(binwidth/2, nchans - binwidth/2):
        med[:,i,:,:] = np.nanmedian(data[:,i-binwidth/2:i+binwidth/2,:,:], axis=1)

    # Handle the edges with a single median across (half the binwidth) channels
    # TODO: Use the right numpy function (tile?) to make this easier to read
    # TODO: Or come up with a nicer way of handling this, based on feedback from the group
    med_start = np.nanmedian(data[:,0:binwidth/2,:,:], axis=1)
    med_end = np.nanmedian(data[:,nchans-binwidth/2:nchans,:,:], axis=1)
    for i in range(0, binwidth/2):
        med[:,i,:,:] = med_start
    for i in range(nchans-binwidth/2, nchans):
        med[:,i,:,:] = med_end

    hdu[0].data = data - med
    hdu.writeto(outfile, overwrite=True)
    hdu[0].data = med 
    hdu.writeto(contfile, overwrite=True)
