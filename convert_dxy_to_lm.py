#!/usr/bin/env python

''' Create files in the format required by https://sourceforge.net/p/wsclean/wiki/ATermCorrection/ '''

from astropy.io import fits
import numpy as np

import glob
import sys

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Warping input/output files")
    group1.add_argument("--prefix", dest='prefix', type=str, default=None,
                        help="The prefix for the input fits files; will append _delx.fits and _dely.fits when it looks for the relevant files. Will also be used as the output with _dldm.fits attached.")

    results = parser.parse_args()

    if len(sys.argv) <= 1 or results.prefix is None:
        parser.print_help()
        sys.exit()

    fits_dx = results.prefix+"_delx.fits"
    fits_dy = results.prefix+"_dely.fits"

    out = results.prefix+"_dldm.fits"

    hdu_dx = fits.open(fits_dx)
    hdu_dy = fits.open(fits_dy)

# Convert from pixels to radians
    dx = np.squeeze(np.squeeze(np.radians(hdu_dx[0].data * hdu_dx[0].header["CDELT2"])))
    dy = np.squeeze(np.squeeze(np.radians(hdu_dy[0].data * hdu_dx[0].header["CDELT2"])))

# Move the frequency header information to the 4th axis
    hdu_dx[0].header["CTYPE4"] = hdu_dx[0].header["CTYPE3"]
    hdu_dx[0].header["CRPIX4"] = hdu_dx[0].header["CRPIX3"]
    hdu_dx[0].header["CRVAL4"] = hdu_dx[0].header["CRVAL3"]
    hdu_dx[0].header["CDELT4"] = hdu_dx[0].header["CDELT3"]
    hdu_dx[0].header["CUNIT4"] = hdu_dx[0].header["CUNIT3"]
    hdu_dx[0].header["CTYPE3"] = "MATRIX"
    hdu_dx[0].header["CRPIX3"] = 1.0
    hdu_dx[0].header["CRVAL3"] = 0.0
    hdu_dx[0].header["CDELT3"] = 1.0
    hdu_dx[0].header["CUNIT3"] = ("", "") # Remove comment about frequency

# Insert the new data
    hdu_dx[0].data = np.array([[-dx, dy]])
    hdu_dx.writeto(out, overwrite=True)
