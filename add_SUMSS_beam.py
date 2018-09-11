#!/usr/bin/env python

import sys
from astropy.io import fits
from numpy import sin, radians
import argparse

def add_sumss_beam(fitsfile, outname):
    hdu = fits.open(fitsfile)
    Dec = hdu[0].header["CRVAL2"]

    bpa = 0.0
    bmin = 45./3600.
    bmaj = bmin / sin(abs(radians(Dec)))

    hdu[0].header["BMAJ"] = bmaj
    hdu[0].header["BMIN"] = bmin
    hdu[0].header["BPA"] = bpa

    if outname is None:
        outname = fitfile.replace(".fits", "_wpsf.fits")
    hdu.writeto(outname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Input/output files")
    group1.add_argument("--infits", dest='infits', type=str, default=None,
                        help="The fits image to add a beam to.")
    group1.add_argument("--outfits", dest='outfits', type=str, default=None,
                        help="The output filename (default = _wpsf).")
    results = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()

    add_sumss_beam(results.infits, results.outfits)

