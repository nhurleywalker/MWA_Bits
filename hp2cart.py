#!/usr/bin/python

# Reproject a healpix map into a Cartesian fits file
# NHW 17/01/16

import sys
import healpy as hp
import matplotlib as ml
ml.use('Agg') # So does not use display
import matplotlib.image as mpimg
import matplotlib.pyplot as plot
import os
import random
import numpy as np

from astropy.table import Table, Column
from astropy.io import fits

from optparse import OptionParser

usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('-i','--input',dest="input",default=None,
                  help="Healpix image <FILE>",metavar="FILE")
parser.add_option('-x','--xpix',dest="xpix",default=2000,
                  help="Number of pixels wide to make the image",type="float")
parser.add_option('-o','--output',dest="output",default="test.fits",
                  help="Output file <FILE>; default = test.fits",metavar="FILE")
(options, args) = parser.parse_args()

if not options.input:
   print "Must specify an input file."
   sys.exit(1)

# Number of pixels L-R, covering 360 deg of RA
xpix=options.xpix
# Resolution in deg
cdelt=360.0/xpix

data=hp.read_map(options.input)
#data=hp.read_map("red_map_hp.fits")
# Dummy figure because cartview always activates matplotlib
dummy_figure=plot.figure(1,figsize=(4,3))
axd = dummy_figure.add_subplot(111)
temp_arr=hp.cartview(data, fig=1, xsize=xpix, cbar=False, title="", coord=["C"], notext=True, return_projected_map=True)
# Unmask and flip
temp_arr=temp_arr.filled(np.nan)

hdulist=fits.open(options.input)
header=hdulist[0].header

# CRVAL2 needs to be zero in order for a plate caree (CAR) projection to resemble a Cartesian grid
header['CRVAL2']=0.0
header['CRPIX2']=xpix/4
header['CDELT2'] = cdelt
header['CTYPE2']='DEC--CAR'

# CRVAL1 should be zero, in the middle of the image
header['CRVAL1']=0.0
header['CRPIX1']=xpix/2
header['CDELT1'] = -cdelt
header['CTYPE1']='RA---CAR'

hdulist[0].data = temp_arr
hdulist[0].header = header
hdulist.writeto(options.output, clobber=True)
