#!/usr/bin/env python

# Plot specific source from GLEAM catalogue

import os, sys

# Need a least-squares estimator that gives a useable error estimate
from scipy.optimize import leastsq

import numpy as np
#tables and votables
import astropy.io.fits as fits
from astropy.io.votable import parse_single_table
from astropy.io.votable import writeto as writetoVO
from astropy.table import Table, Column
#import aplpy

import matplotlib.pyplot as plt

from optparse import OptionParser

usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('--catalogue',type="string", dest="catalogue",
                    help="The filename of the catalogue you want to read in.", default=None)
#parser.add_option('--outpng',type="string", dest="outpng",
#                    help="The filename of the output png.", default=None)
#parser.add_option('--plot',action="store_true",dest="make_plots",default=False,
#                  help="Make fit plots? (default = False)")
#parser.add_option('--order',dest="poly_order",default=1,type=int,
#                  help="Set the order of the polynomial fit. (default = 1)")
parser.add_option('--source',dest="source",default=None,type="string",
                  help="Name of the source to plot (use quotes)")
parser.add_option('--list',dest="sourcelist",default=None,type="string",
                  help="List of sources in VO format to plot; will check \"Name\" column; overrides \"--source\" option.")
parser.add_option('--curve',dest="fitcurve", action="store_true", default=False,
                  help="Plot a curved spectrum using beta column (default False)")
#parser.add_option('--mosaic',dest="mosaic",default=None,type="string",
#                  help="Mosaic from which to get a cutout image (default = don't make image plot)")
(options, args) = parser.parse_args()

# http://scipy-cookbook.readthedocs.org/items/FittingData.html
# Define function for calculating a power law
powerlaw = lambda x, amp, index: amp * (x**index)
curvepowlaw = lambda x, amp, index, q: amp*(x**(index))*np.exp(q*np.power(np.log(x),2))

# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
pinit = [-2.0, -0.7]

# define the curve fitting function
fitcurv = lambda p, x: p[0] + p[1] * x + p[2] *x**2
errcurv = lambda p, x, y, err: (y - fitcurv(p, x)) / err
pinit_curv = [-2.0, -0.7, 0.0]

if options.sourcelist:
    if os.path.exists(options.sourcelist):
        temp = parse_single_table(options.sourcelist)
        sources = temp.array["Name"]
        if sources is None:
            print("No sources found in the \"Name\" column of "+options.sourcelist)
            sys.exit(1)
    else:
        print(options.sourcelist+" does not exist.")
        sys.exit(1)
else:
    if options.source is None:
        print("Must select a source to plot!")
        sys.exit(1)
    else:
        sources = [options.source]

if options.catalogue is None:
    print("must specify input catalogue using the --catalogue option")
    sys.exit(1)
else:
    filename, file_extension = os.path.splitext(options.catalogue)
    if file_extension == ".fits":
        temp = fits.open(options.catalogue)
        data = temp[1].data
    elif file_extension == ".vot":
        temp = parse_single_table(options.catalogue)
        data = temp.array

# Frequencies to write out -- use the full band since it shows the range I did the fitting over
freq1=72
freq2=231

# Can't figure out how to read the bloody column names! hardcode
freqs=["076", "084", "092", "099", "107",  "115", "122", "130", "143", "151", "158", "166", "174",  "181", "189","197", "204", "212", "220", "227"]
#freqs=["151", "158", "166", "174",  "181", "189","197", "204", "212", "220", "227"]
#freqs=["204", "212", "220", "227"]
#freqs=["76", "84", "92", "99", "107",  "115", "122", "130", "143", "151", "158", "166", "174",  "181", "189","197", "204", "212", "220", "227"]

for source in sources:
    # Replace any underscores for searching catalogue
    source=source.replace("_"," ")
    # But put them back in for naming the output
    outpng=(source.replace(" ","_"))+".png"

    index = np.where(data["Name"]==source)
    # select fluxes
    flux_list = []
    for x in freqs:
        flux_list.append(np.squeeze(np.ma.log(data["int_flux_"+x][index])))
# HACK
#    flux_list.append(np.ma.log(0.267))
#    flux_list.append(np.ma.log(0.383))
#    flux_list.append(np.ma.log(0.869))
    flux_array = np.ma.masked_array(flux_list)
#    flux_array = np.ravel(np.ma.masked_array(flux_list))
    #flux_array = np.ravel(np.asarray(flux_list))

    # Representative calibration error: 2% at good decs, 3% elsewhere
    if (data["DEJ2000"][index] > 18.5) or (data["DEJ2000"][index] < -72.0):
        calibration_error = 0.03
    else:
        calibration_error = 0.02

    err_list = []
    for x in freqs:
        fitting_error = data["err_int_flux_"+x][index]/data["int_flux_"+x][index]
        err_list.append(np.squeeze(np.sqrt(fitting_error**2 + calibration_error**2)))
# Quick HACK
#    err_list.append(0.1)
#    err_list.append(0.1)
#    err_list.append(0.1)
    flux_errors = np.ravel(np.asarray(err_list))


# Quick HACK
#    freqs.append(1400)
#    freqs.append(843)
#    freqs.append(150)
     
    weights = 1/(flux_errors*flux_errors)
    freq_array = np.log(np.ma.array([float(int(x)) for x in freqs]),dtype="float64")
    #print freq_array
    #print flux_array
    #print flux_errors


    if options.fitcurve:
        fit = leastsq(errcurv, pinit_curv, args=(freq_array, flux_array, flux_errors), full_output=1)
        covar = fit[1]
        if covar is not None:
            P = fit[0]
            #print P
            alpha = P[1]
            beta = P[2]
            amp = np.exp(P[0])
            flux1 = curvepowlaw(freq1,amp,alpha,beta)
            flux2 = curvepowlaw(freq2,amp,alpha,beta)
            #print flux1, flux2
        # Errors
            err_beta = np.sqrt(covar[2][2])
            err_alpha = np.sqrt(covar[1][1])
            err_flux1 = np.sqrt(covar[0][0])*flux1
            err_flux2 = np.sqrt(covar[0][0])*flux2
            residual = errcurv(P,freq_array, flux_array, flux_errors)
            chi2red = sum(np.power(residual,2))/(len(freqs)-len(pinit))
        else:
            beta=None
            alpha=None
            amp=None
            flux1=None
            flux2=None
            err_beta=None
            err_alpha=None
            err_flux1=None
            err_flux2=None

    #indices = np.where(np.bitwise_not(np.isnan(alpha)))
    else:
        fit = leastsq(errfunc, pinit, args=(freq_array, flux_array, flux_errors), full_output=1)
        covar = fit[1]
        if covar is not None:
            P = fit[0]
            alpha=P[1]
            amp = np.exp(P[0])
            flux1=powerlaw(freq1,amp,alpha)
            flux2=powerlaw(freq2,amp,alpha)
        # Errors
            err_alpha = np.sqrt(covar[1][1])
            err_flux1 = np.sqrt(covar[0][0])*flux1
            err_flux2 = np.sqrt(covar[0][0])*flux2
            residual = errfunc(P,freq_array, flux_array, flux_errors)
            chi2red = sum(np.power(residual,2))/(len(freqs)-len(pinit))
        else:
            alpha=None
            amp=None
            flux1=None
            flux2=None
            err_alpha=None
            err_flux1=None
            err_flux2=None

    #indices = np.where(np.bitwise_not(np.isnan(alpha)))
#    if options.mosaic is not None:
#        nplots = 3
#    else:
    nplots = 2

    # Plot
    example=plt.figure(figsize=(10,5))
    ax2=example.add_subplot(1,nplots,2)
    ax1=example.add_subplot(1,nplots,1)
    if alpha is None:
        plt.title("{0:s}".format(data["Name"][index][0]))
        ax2.set_xscale("log")
        ax2.set_yscale("log")
    else:
#        plt.title("{0:s}: alpha={1:3.2f}+/-{2:3.2f} ; reduced chi2={3:4.2f}".format(data["Name"][index][0],alpha,err_alpha,chi2red))
        if options.fitcurve:
            #print np.exp(freq_array), curvepowlaw(np.exp(freq_array), amp, alpha, beta)
#            ax1.plot(np.exp(freq_array), curvepowlaw(np.exp(freq_array), amp, alpha, beta))     # Fit
            #print amp, alpha, beta
            ax1.set_title("{0:s}: beta = {3:3.2f}+/-{4:3.2f} alpha = {1:3.2f}+/-{2:3.2f}".format(data["Name"][index][0],alpha,err_alpha,beta,err_beta))
#            ax2.loglog(np.exp(freq_array), curvepowlaw(np.exp(freq_array), amp, alpha, beta))     # Fit
            ax2.set_title("reduced chi2 = {0:4.2f}".format(chi2red))
        else:
#            ax1.plot(np.exp(freq_array), powerlaw(np.exp(freq_array), amp, alpha))     # Fit
            ax1.set_title("{0:s}: alpha = {1:3.2f}+/-{2:3.2f}".format(data["Name"][index][0],alpha,err_alpha))
#            ax2.loglog(np.exp(freq_array), powerlaw(np.exp(freq_array), amp, alpha))     # Fit
            ax2.set_title("reduced chi2 = {0:4.2f}".format(chi2red))
    ax1.errorbar(np.exp(freq_array), np.exp(flux_array), yerr=flux_errors*np.exp(flux_array), fmt='k.')  # Data
    ax1.set_xlabel("Frequency (MHz)")
    ax1.set_ylabel("Flux density (Jy)")
    ax2.errorbar(np.exp(freq_array), np.exp(flux_array), yerr=flux_errors*np.exp(flux_array), fmt='k.')  # Data
    ax2.set_xlim([0.9*min(np.exp(freq_array)),1.1*max(np.exp(freq_array))])
    ax2.set_ylim([0.9*min(np.exp(flux_array)),1.1*max(np.exp(flux_array))])
    ax2.set_xlabel("Log Frequency (MHz)")
    ax2.set_ylabel("Log Flux density (Jy)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    example.savefig(outpng)

# Too large for aplpy
    #if options.mosaic is not None:
    #    fig=aplpy.FITSFigure(options.mosaic,figure=example,subplot=(1,nplots,3))
    #    RA = data["RAJ2000"][index]
    #    Dec = data["DEJ2000"][index]
    #    fig.recenter(RA,Dec,width=1.0,height=1.0)
    #    fig.show_colorscale(vmin=-0.1*min(np.exp(flux_array)),vmax=max(np.exp(flux_array)),cmap="cubehelix")

