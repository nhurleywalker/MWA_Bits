#!/usr/bin/python

# Plot a healpix image to PNG
# NHW 02/12/2015

import healpy as hp
import numpy as np
import matplotlib as ml
ml.use('Agg') # So does not use display
import matplotlib.pyplot as plot
import os
import random

try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from optparse import OptionParser

usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('-f','--filename',dest="filename",default=None,
                  help="Input map to read <FILE>",metavar="FILE")
parser.add_option('-c','--cmap',dest="cmap",default="Reds_r",
                  help="Matplotlib color map to use; default = Reds reversed (Reds_r)",type="string")
parser.add_option('-o','--output',dest="output",default="test.png",
                  help="Output file <FILE>; default = test.png",metavar="FILE")
(options, args) = parser.parse_args()

# Padding
sbplt_pad_left  = 0  # the left side of the subplots of the figure
sbplt_pad_right = 0    # the right side of the subplots of the figure
sbplt_pad_bottom = 0   # the bottom of the subplots of the figure
sbplt_pad_top = 0      # the top of the subplots of the figure
sbplt_pad_wspace = 0   # the amount of width reserved for blank space between subplots
sbplt_pad_hspace = 0   # the amount of height reserved for white space between subplots

image=hp.read_map(options.filename)

if options.output:
   output=options.output
else:
   output="test.png"

# Other option: return as a numpy array
#dummy_figure=plot.figure(1,figsize=(40,30))
#axd = dummy_figure.add_subplot(111)
## Return as an array
#map_array=hp.cartview(image,fig=1,xsize=21600,latra=[-90,30],min=-0.1,max=1,return_projected_map=True)
#new_map_array=np.flipud(map_array)
##map_array=hp.cartview(image, fig=1, xsize=200,cmap=options.cmap, cbar=False,title="",latra=[-90,30],min=-0.1,max=1,return_projected_map=True)
# Then do a normal matplotlib plot
#real_figure=plot.figure(2,figsize=(16,12))
#axr = real_figure.add_subplot(111)
#axr.set_axis_off()
#axr.imshow(new_map_array,cmap=options.cmap,vmin=-0.1,vmax=1)
#real_figure.savefig(output,pad_inches=0.0,bbox_inches='tight',dpi=500))

# Simpler option (seems to work): just use hp.cartview with various options
# Have to make a massive figure in order to not be limited by pixel number
# At some point I should make this dynamic depending on the size of the healpix image!
plot.figure(1,figsize=(40,30))
hp.cartview(image, fig=1, xsize=26000,cmap=options.cmap,cbar=False,title="",min=-0.05,max=1,coord=["C","G"],notext=True)

if os.path.exists(output):
    os.remove(output)
plot.savefig(output,pad_inches=0.0,bbox_inches='tight',dpi=500) 

tempfile="temp_image"+str(random.randrange(1000000,9999999))+".png"
os.system("convert "+output+" -trim +repage "+tempfile)
os.remove(output)
os.rename(tempfile,output)
