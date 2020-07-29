#!/usr/bin/python

# Regrid fits files to HEALPix
# Tested on MWA images
# Modified from Herschel regridding code (C. North, Cardiff)
# NHW 02/12/2015

import healpy as hp
import astropy.wcs as wcs
from astropy.table import Table,Column
import os
import numpy as np
import sys
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
try:
    from progress.bar import Bar
    progress=True
except ImportError:
    print "No progress bar module found."
    progress=False
from optparse import OptionParser

usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('-f','--filename',dest="filename",default=None,
                  help="Input file to read <FILE>",metavar="FILE")
parser.add_option('-o','--output',dest="output",default="test.fits",
                  help="Output file <FILE>; default = test.fits",metavar="FILE")
parser.add_option('-n','--nside',dest="nside",default=1024,
                  help="Number of sides for healpix sphere (must be a power of 2); default = 1024",type=int)
parser.add_option('-u','--unseen',action="store_true",dest="unseen",default=False,
                  help="Use -1.67E30 as the unseen pixel value instead of 0.0: good for 'normal maps', bad for real-positive maps like RMS maps (default=False)")
(options, args) = parser.parse_args()

d2r=np.pi/180.
r2d=1./d2r

#Author: A.Polino
def is_power2(num):
    'states if a number is a power of two'
    return num != 0 and ((num & (num - 1)) == 0)

def mag(x):
    return int(np.log10(x))

# Modified from
# http://stackoverflow.com/questions/29702010/healpy-pix2ang-convert-from-healpix-index-to-ra-dec-or-glong-glat
def IndexToRADec(index,nside):
    theta,phi=hp.pixelfunc.pix2ang(nside,index)
#    return [np.degrees(np.pi*2.-phi),-np.degrees(theta-np.pi/2.)]
    return [np.degrees(phi),-np.degrees(theta-np.pi/2.)]
def RADecToIndex(RA,decl,nside):
    return hp.pixelfunc.ang2pix(nside,np.radians(-decl+90.),np.radians(360.-RA))

nside=options.nside
if not is_power2(nside):
    print "nside needs to be a power of 2"
    sys.exit(1)

# Open file
if os.path.exists(options.filename):
    hdu_in = pyfits.open(options.filename)
else:
    print options.filename+" doesn't exist!"
    sys.exit(1)

# Is it an MWA snapshot or a flat mosaic?
if len(hdu_in[0].data.shape)>2:
    #format is e.g. RA, Dec, stokes, spectral -- usual for MWA snapshots
    sfxy=False
else:
    # format is e.g. RA, Dec -- usual for mosaics
    sfxy=True

# Find centre of map
header=hdu_in[0].header
w = wcs.WCS(header,naxis=2)
naxis1=header['NAXIS1']
naxis2=header['NAXIS2']
raCtr=header['CRVAL1']
if sfxy:
    decCtr=header['CRVAL2']
# NB: Cannot use decCtr=header['CRVAL2'] for w-snapshot zenith maps
# Instead, what is the co-ordinate of the central pixel?
else:
    decCtr=float(w.wcs_pix2world(int(naxis1/2),int(naxis2/2),1)[1])

# Swarp uses CD instead of CDELT
if 'CDELT1' in header:
    cdelt1 = header['CDELT1']
elif 'CD1_1' in header:
    cdelt1 = header['CD1_1']
else:
    print "Error: Can't find CDELT1 or CD1_1"
if 'CDELT2' in header:
    cdelt2 = header['CDELT2']
elif 'CD2_2' in header:
    cdelt2 = header['CD2_2']
else:
    print "Error: Can't find CDELT2 or CD2_2"

# work out size of map
dx=np.abs(cdelt1*naxis1)
dy=cdelt2*naxis2
print '  map dimensions: (%.3f x %.3f degree)'%(dx,dy)

# Get co-latitude and longitude
lonCtr=raCtr*d2r
colatCtr=(90.-decCtr)*d2r
##Convert to vector
vecCtr=hp.ang2vec(colatCtr,lonCtr)

##Use 10% extra radius (in radians)
discRad=0.55*max([dx,dy])*np.pi/180.
#discRad=2*np.pi

##Get HP pixel numbers around centre of map
pixDisc=hp.query_disc(nside,vecCtr,discRad,True)
print '  found %d HEALPix pixels (%.3f degree radius)'%(len(pixDisc),discRad*r2d)

# Report progress in approximately 1% increments
granularity=10**(mag(len(pixDisc))-2)

validPix=[]
pixMap=[]
if progress:
    bar = Bar('Processing', max=int(len(pixDisc)/granularity))
for p in pixDisc:
    pixH=w.wcs_world2pix([IndexToRADec(p,nside)],1)[0]
# pixH is in x,y (usually, roughly, -ra, dec)
    if pixH[0]>=0 and pixH[0]<naxis1-1 and pixH[1]>=0 and pixH[1]<naxis2-1:
# fits files are in dec, ra
        if sfxy:
            pixVal=hdu_in[0].data[pixH[1],pixH[0]]
        else:
            pixVal=hdu_in[0].data[0,0,pixH[1],pixH[0]]
#        pixVal=map.getIntensity(pixH[0],pixH[1])
        if pixVal == pixVal:
            validPix.append(p)
            pixMap.append(pixVal)
#            sys.exit(0)
    if p%granularity==0:
        if progress:
            bar.next()
        else:
            print str(p)+"/"+str(len(pixDisc))
if progress:
    bar.finish()
print '  written %d valid HEALPix pixels'%(len(validPix))

##Construct table of pixel numbers and values for that map
#Users astropy syntax
pixTable=Table()
pixTable.add_column(Column(validPix,name='index'))
pixTable.add_column(Column(pixMap,name='value'))

##Write that map's pixel values to file
if os.path.exists(options.output):
   os.remove(options.output)
pixTable.write(options.output)

#####THAT'S MADE A LIST OF HEALPIX PIXEL INDICES AND VALUES FOR ONE MAP

###Make full-sky map
npix=hp.nside2npix(nside)
hpMapFullsky=np.zeros(npix)
if options.unseen:
    # Sets them to -1.67x10^30 : not good for the RMS maps!
    hpMapFullsky[:]=hp.UNSEEN

#Populate map with the pixels above
hpMapFullsky[pixTable['index'].data] = pixTable['value'].data

#Write to file (assuming Celestial coords; use coord='G' for Galactic)
hp.write_map(options.output,hpMapFullsky,coord='C')

##plot figure as a test
#plot.figure()
#hp.cartview(hpMapFullsky,norm='log',coord=['C','G'])
