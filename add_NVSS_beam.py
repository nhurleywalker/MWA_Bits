#!/usr/bin/env python

import sys
from astropy.io import fits
from numpy import cos, radians
import argparse
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

def add_nvss_beam(fitsfile, outname):
    hdu = fits.open(fitsfile)

    bpa = 0.0
    bmaj = 45./3600.
    bmin = 45./3600.

    hdu[0].header["BMAJ"] = bmaj
    hdu[0].header["BMIN"] = bmin
    hdu[0].header["BPA"] = bpa

    if outname is None:
        outname = fitsfile.replace(".fits", "_wpsf.fits")
    hdu.writeto(outname, overwrite=True)

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

    add_nvss_beam(results.infits, results.outfits)

