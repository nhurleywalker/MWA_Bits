#!/usr/bin/env python
import os, logging
from optparse import OptionParser #NB zeus does not have argparse!

import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
#import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
#import pylab
from astropy.io import fits
from mwapy import aocal

def get_tile_info(metafits):
    hdus = fits.open(metafits)
    inputs = hdus[1].data
    tiles = inputs[inputs['pol'] == 'X']
    Names = tiles["TileName"]
    North = tiles["North"]
    East = tiles["East"]
    return Names, North, East

def diff(ao, metafits, refant):
    diffs = []
    t_start = 0
    t_end = ao.n_int - 1
# Divide through by refant
    ao = ao / ao[:, refant, :, :][:, np.newaxis, :, :]
    ant_iter = xrange(ao.n_ant)
    for a, antenna in enumerate(ant_iter):
        temp = []
# Only XX and YY
        for pol in 0, 3:
# Difference the complex gains, then convert to angles
           temp.append(np.angle(ao[t_end, antenna, :, pol] / ao[t_start, antenna, :, pol], deg=True))
        diffs.append(temp)
    return diffs

def histo(diffs, outname):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, patches = ax.hist(diffs, bins = 60, range=[-180, 180])
    peak = bins[np.where(n == n.max())][0]
    ax.axvline(x=np.median(diffs), color="red")
    ax.axvline(x=peak, color="orange")
    ax.set_xlabel("Phase change / degrees")
    at = AnchoredText("Median: {0:3.0f}deg\nPeak: {1:3.0f}deg\nStdev: {2:3.0f}deg".format(np.median(diffs), peak, np.std(diffs)),
                  prop=dict(size=8), frameon=True,
                  loc=1,
                  )
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    print np.median(diffs), peak, np.std(diffs)
    fig.savefig(outname)

def phase_map(diffs, metafits, names, outname):
    fig = plt.figure(figsize = (10,8))
    Names, North, East = get_tile_info(metafits)
    ax = fig.add_axes([0.15, 0.1, 0.65, 0.75])
    ax.axis("equal")
    sc = ax.scatter(North, East, marker='o', s=150, linewidths=4, c=diffs, cmap=plt.cm.hsv, vmin = -180., vmax = 180.)
    ax.set_xlabel("East / m")
    ax.set_ylabel("North / m")
    if options.names is True:
        for i, txt in enumerate(Names):
            ax.annotate(txt, (North[i], East[i]))

    cbaxes = fig.add_axes([0.82, 0.1, 0.02, 0.75])
    cb = plt.colorbar(sc, cax = cbaxes, orientation="vertical")
    cb.set_label('Phase change / degrees')
    fig.savefig(outname)

if __name__ == '__main__':
    parser = OptionParser(usage = "usage: %prog binfile" +
    """
    Difference time-based calibration solutions to determine ionspheric variation
    """)
    parser.add_option("--refant", default=127, dest="refant", type="int", help="Default = 127")
    parser.add_option("-m", "--metafits", default=None, dest="metafits", help="metafits file (must be supplied to generate phase map")
    parser.add_option("-v", "--verbose", action="count", dest="verbose", help="-v info, -vv debug")
    parser.add_option("--outdir", default=None, dest="outdir", help="output directory [default: same as binfile]")
    parser.add_option("--names", default=False, dest="names", help="Plot tile names on phase map")
#    parser.add_option("--output", default=None, dest="output", help="output names [default: OBSID_histogram.png and OBSID_phasemap.png")
#    parser.add_option("--marker", default=',', dest="marker", type="string", help="matplotlib marker [default: %default]")
#    parser.add_option("--markersize", default=2, dest="markersize", type="int", help="matplotlib markersize [default: %default]")
    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("incorrect number of arguments")

    filename = args[0]
    ao = aocal.fromfile(filename)
    obsid = filename[0:10]
    histname = obsid+"_histogram.png"
    mapname = obsid+"_phasemap.png"

    diffs = np.array(diff(ao, options.metafits, options.refant))
    print diffs.shape
# Flatten array and delete NaNs for histogram
    histo(diffs[np.logical_not(np.isnan(diffs))].flatten(), histname)
    if options.metafits is not None:
        if os.path.exists(options.metafits):
# Plotting on a single frequency, single pol on map because it's impossible otherwise
            diffs = diffs[:, 0, 15]
# Could also take the average of a few frequencies, but it doesn't change anything
#            diffs = np.average(diffs[:, 0, 12:20], axis=1)
            phase_map(diffs, options.metafits, options.names, mapname)
           
#    plot(ao, os.path.splitext(args[0])[0]+opts.suffix, opts.refant, plot_title = opts.plot_title, outdir=opts.outdir, format=opts.format, amp_max=opts.amp_max, marker=opts.marker, markersize=opts.markersize, verbose=opts.verbose, metafits=opts.metafits)
