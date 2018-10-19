#!/usr/bin/env python

import sys
from astropy.io import fits
from numpy import cos, radians
import argparse
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

def add_tgss_beam(fitsfile, outname):
    hdu = fits.open(fitsfile)
    w = wcs.WCS(hdu[0].header)
    lon, lat = w.wcs_pix2world([[hdu[0].header["NAXIS1"]/2,hdu[0].header["NAXIS2"]/2]],0).transpose()
    if hdu[0].header["CTYPE1"] == "GLON":
        coords = SkyCoord(lon, lat, frame="galactic", unit = (u.deg, u.deg))
    else:
        coords = SkyCoord(lon, lat, frame="fk5", unit = (u.deg, u.deg))
        
    Dec = coords.fk5.dec.value[0]

#    Dec = hdu[0].header["CRVAL2"]

    bpa = 0.0
    bmin = 25./3600.
    if Dec > 19.:
        bmaj = bmin
    else:
        bmaj = bmin / cos(abs(radians(Dec-19.)))

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

    add_tgss_beam(results.infits, results.outfits)

