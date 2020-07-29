#!/usr/bin/python

# Make an RGB image out of three input healpix maps
# NHW 16/01/16

import sys
import healpy as hp
import matplotlib as ml
ml.use('Agg') # So does not use display
import matplotlib.image as mpimg
import matplotlib.pyplot as plot
import os
import random

try:
    from numpy import stack
except ImportError:
    print "This script requires at least numpy 1.10 for the array stack functionality."
    print "You are using version "+str(numpy.__version__)
    sys.exit(1)

import numpy as np

from optparse import OptionParser

usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('-r','--red',dest="red",default=None,
                  help="Red image <FILE>",metavar="FILE")
parser.add_option('-g','--green',dest="green",default=None,
                  help="Green image <FILE>",metavar="FILE")
parser.add_option('-b','--blue',dest="blue",default=None,
                  help="Blue image <FILE>",metavar="FILE")
parser.add_option('--rmin',dest="rmin",default=-0.1,
                  help="Minimum value for red color range",type="float")
parser.add_option('--rmax',dest="rmax",default=1,
                  help="Maximum value for red color range",type="float")
parser.add_option('--gmin',dest="gmin",default=-0.1,
                  help="Minimum value for green color range",type="float")
parser.add_option('--gmax',dest="gmax",default=1,
                  help="Maximum value for green color range",type="float")
parser.add_option('--bmin',dest="bmin",default=-0.1,
                  help="Minimum value for blue color range",type="float")
parser.add_option('--bmax',dest="bmax",default=1,
                  help="Maximum value for blue color range",type="float")
parser.add_option('-t','--transfer',dest="tfunc",default="linear",
                  help="Transfer function to use, of linear, asinh, or log",type="string")
parser.add_option('-o','--output',dest="output",default="test_rgb.png",
                  help="Output file <FILE>; default = test_rgb.png",metavar="FILE")
(options, args) = parser.parse_args()

if (not options.red) or (not options.green) or (not options.blue):
   print "Must supply all three of R, G and B images."
   sys.exit(1)

tfuncs=["linear","asinh","log"]

if options.tfunc not in tfuncs:
   print "Transfer function must be linear, asinh, or log."

output=options.output

data=[[hp.read_map(options.red),options.rmin,options.rmax],[hp.read_map(options.green),options.gmin,options.gmax],[hp.read_map(options.blue),options.bmin,options.bmax]]
# Dummy figure because cartview always activates matplotlib
dummy_figure=plot.figure(1,figsize=(4,3))
axd = dummy_figure.add_subplot(111)

png_data=[]
for line in data:
    temp_arr=hp.cartview(line[0], fig=1, xsize=26000, cbar=False, title="", coord=["C","G"], notext=True, return_projected_map=True)
# Unmask and flip
    temp_arr=np.flipud(temp_arr.filled(0.0))
# Saturate any data below minimum or above maximum
    temp_arr[np.where(temp_arr<line[1])]=line[1]
    temp_arr[np.where(temp_arr>line[2])]=line[2]
# Apply transfer function
    if options.tfunc=="asinh":
       temp_arr=np.arcsinh(temp_arr)
    elif options.tfunc=="log":
       temp_arr=np.log(temp_arr)
# Normalise the array to range of 0 to 1 for PNG output
    temp_arr=temp_arr-line[1]
    temp_arr=temp_arr/line[2]
    png_data.append(temp_arr)
# Set alpha to 1 = fully opaque
alpha=np.ones(png_data[0].shape)
# Stack RGB and Alpha
png_array=np.stack([png_data[0],png_data[1],png_data[2],alpha],axis=2)

# Then do a normal matplotlib plot
real_figure=plot.figure(2,figsize=(40,30))
axr = real_figure.add_subplot(111)
axr.set_axis_off()
axr.imshow(png_array)
real_figure.savefig(output,pad_inches=0.0,bbox_inches='tight',dpi=500)

# Crop
tempfile="temp_image"+str(random.randrange(1000000,9999999))+".png"
os.system("convert "+output+" -trim +repage "+tempfile)
os.remove(output)
os.rename(tempfile,output)


