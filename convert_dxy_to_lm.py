#!/usr/bin/env python

from astropy.io import fits
import astropy.units as u
from astropy import wcs
import numpy as np

from datetime import datetime
import glob
import sys

import argparse

def RADecToLM(RA, Dec, RA_pc, Dec_pc):
    deltaAlpha = RA - RA_pc
    sinDeltaAlpha = np.sin(deltaAlpha)
    cosDeltaAlpha = np.cos(deltaAlpha)
    sinDec = np.sin(Dec)
    cosDec = np.cos(Dec)
    sinDec0 = np.sin(Dec_pc)
    cosDec0 = np.cos(Dec_pc)
    l = cosDec * sinDeltaAlpha
    m = sinDec*cosDec0 - cosDec*sinDec0*cosDeltaAlpha
    return l, m

def unwrap(RA):
    if RA < 0.:
        RA = 360. + RA
    return RA

vunwrap = np.vectorize(unwrap)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Warping input/output files")
    group1.add_argument("--prefix", dest='prefix', type=str, default=None,
                        help="The prefix for the input fits files; will append _delx.fits and _dely.fits when it looks for the relevant files. Will also be used as the output with _dl.fits and _dm.fits attached.")
    group1.add_argument("--ra", dest='RA_pc', type=float, default=None,
                        help="RA of the phase centre, default is to obtain from centre pixel")
    group1.add_argument("--dec", dest='Dec_pc', type=float, default=None,
                        help="Dec of the phase centre, default is to obtain from centre pixel")

    results = parser.parse_args()

    if len(sys.argv) <= 1 or results.prefix is None:
        parser.print_help()
        sys.exit()

    fits_dx = results.prefix+"_delx.fits"
    fits_dy = results.prefix+"_dely.fits"

    out = results.prefix+"_dldm.fits"
#    out_dm = results.prefix+"_dm.fits"

    if results.RA_pc is None or results.Dec_pc is None:
        hdu = fits.open(fits_dx)
        cent = hdu[0].data.shape[0] / 2
        w = wcs.WCS(hdu[0].header, naxis=2)
        RA_pc , Dec_pc = w.wcs_pix2world([[cent, cent]], 0).transpose()
        hdu.close()
    else:
        RA_pc , Dec_pc = results.RA_pc, results.Dec_pc

    hdu_dx = fits.open(fits_dx)
    hdu_dy = fits.open(fits_dy)

    w = wcs.WCS(hdu_dx[0].header, naxis=2)

    print("creating indices at {0}".format(datetime.now()))
    #create an array but don't set the values (they are random)
    indexes = np.empty( (hdu_dx[0].data.shape[0]*hdu_dx[0].data.shape[1],2),dtype=int)
    #since I know exactly what the index array needs to look like I can construct
    # it faster than list comprehension would allow
    #we do this only once and then recycle it
    idx = np.array([ (j,0) for j in xrange(hdu_dx[0].data.shape[1])])
    j=hdu_dx[0].data.shape[1]
    for i in xrange(hdu_dx[0].data.shape[0]):
        idx[:,1]=i
        indexes[i*j:(i+1)*j] = idx
    #put ALL the pixles into our vectorized functions and minimise our overheads
    print("Finding RAs and Decs at {0}".format(datetime.now()))
    RA, Dec = w.wcs_pix2world(indexes,0).transpose()
    print("Unwrapping RA at {0}".format(datetime.now()))
    RA = vunwrap(RA)

    print("Converting to l,m at {0}".format(datetime.now()))
    l, m = RADecToLM(np.radians(RA), np.radians(Dec), np.radians(RA_pc), np.radians(Dec_pc))

    print("Making copies of dx,dy at {0}".format(datetime.now()))
    dx = np.copy(hdu_dx[0].data).reshape(hdu_dx[0].data.shape[0]*hdu_dx[0].data.shape[1])
    dy = np.copy(hdu_dy[0].data).reshape(hdu_dy[0].data.shape[0]*hdu_dy[0].data.shape[1])

    tind = indexes.T.astype(np.float32)

    tind[0] += dx
    tind[1] += dy

    indexes = tind.T
    print("Finding RAs and Decs of tips at {0}".format(datetime.now()))
    RA, Dec = w.wcs_pix2world(indexes,0).transpose()
    print("Unwrapping RA at {0}".format(datetime.now()))
    RA = vunwrap(RA)
    print("Converting to l,m at {0}".format(datetime.now()))
    l_tips, m_tips = RADecToLM(np.radians(RA), np.radians(Dec), np.radians(RA_pc), np.radians(Dec_pc))

    print("Finding d_l, d_m at {0}".format(datetime.now()))
    d_l = l_tips - l
    d_m = m_tips - m

    print("Reshaping d_l, d_m at {0}".format(datetime.now()))
    d_l = d_l.reshape((hdu_dx[0].data.shape[0],hdu_dx[0].data.shape[1]))
    d_m = d_m.reshape((hdu_dx[0].data.shape[0],hdu_dx[0].data.shape[1]))
    print("Writing d_l, d_m at {0}".format(datetime.now()))
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
    hdu_dx[0].data = np.array([[d_l.astype(np.float32), d_m.astype(np.float32)]])
    hdu_dx.writeto(out, overwrite=True)
    
#    hdu_dx[0].data = d_l.astype(np.float32)
#    hdu_dx.writeto(out_dl, overwrite=True)
#    hdu_dx[0].data = d_m.astype(np.float32)
#    hdu_dx.writeto(out_dm, overwrite=True)

