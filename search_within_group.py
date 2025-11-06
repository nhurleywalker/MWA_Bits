#!/usr/bin/env python

from glob import glob
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from astropy.nddata import Cutout2D
import numpy as np
import os

# First we need to read all the tables that we've been fed
# For now, hardcode; later we figure out a way of specifying inputs

hdus = sorted(glob("*warp-corrected_wbeam.fits"))

minra = 0.0
maxra = 360
mindec = -90
maxdec = 90
for h in hdus:
    tab = Table(fits.open(h)[1].data)
    minra = np.nanmax([minra, np.nanmin(tab["ra"])])
    maxra = np.nanmin([maxra, np.nanmax(tab["ra"])])
    mindec = np.nanmax([mindec, np.nanmin(tab["dec"])])
    maxdec = np.nanmin([maxdec, np.nanmax(tab["dec"])])

# Shave off a small margin -- this avoids really terrible parts of the images
minra+=2
maxra-=2
mindec+=2
maxdec-=2

# Create a matched table with STILTS
nin = len(hdus)
xrad = 25 #arcsec
outfile = "test_auto_join.fits"

# Downselect to sources within the region we care about
# Note the slightly frustrating quotes
icmd = "".join([f"icmd{i+1}=\'select \"ra > {minra} && ra < {maxra} && dec < {maxdec} && dec > {mindec}\"\' " for i in range(0, nin)])
# Always include an entry -- a transient could be in a single frame
joincmd = "".join([f"join{i+1}='always' " for i in range(0, nin)])
# Sky crossmatch of 25" was found to be necessary -- this is not the ionosphere, this is sources being decomposed in different ways
valuescmd = "".join([f"values{i+1}='ra dec' " for i in range(0, nin)])
incmd = "".join([f"in{i+1}={hdus[i]} " for i in range(0, nin)])
stiltscmd = f"nin={nin} matcher='sky' params={xrad} multimode='group' out={outfile} {incmd} {icmd} {valuescmd} {joincmd}"
os.system(f"stilts tmatchn {stiltscmd}")

jointab = Table(fits.open(outfile)[1].data)

# Have to select compactness after join or we will occasionally mark slightly extended sources as transient
# We have to do a NOT(NOT EXTENDED) because that also includes the 'NaN' sources, i.e. non-detections
masks = np.empty((nin, len(jointab)), dtype='bool')
for i in range(0, nin):
    masks[i, :] = ~(jointab[f"int_flux_{i+1}"]/jointab[f"peak_flux_{i+1}"] > 1.5)

# The final mask is the sources that are nan or compact
compact_mask = np.nanmin(masks, axis=0)

# We want to find out the average RAs of all the sources and the average Decs of all the sources so that we can look up where they are in the maps if we need to
jointab['mean_ra'] = np.nanmean([jointab[f'ra_{i+1}'] for i in range(0, nin)], axis=0)
jointab['mean_dec'] = np.nanmean([jointab[f'dec_{i+1}'] for i in range(0, nin)], axis=0)
coords = SkyCoord(jointab['mean_ra'], jointab['mean_dec'], frame='fk5', unit=(u.deg, u.deg))
# We also want to know the rough centroid of the observations so we can exclude the edges
cent = SkyCoord(np.median(jointab['mean_ra']), np.median(jointab['mean_dec']), frame='fk5', unit=(u.deg, u.deg))
spatial_mask = coords.separation(cent) < 10*u.deg
# We will later want the average flux density of sources BEFORE we have populated them with zeros
jointab['mean_on_peak_flux'] = np.nanmean([jointab[f'peak_flux_{i+1}'] for i in range(0, nin)], axis=0)
jointab['mean_on_local_rms'] = np.nanmean([jointab[f'local_rms_{i+1}'] for i in range(0, nin)], axis=0)

# Now we have compact transient sources, for every missing entry, go and look up the RMS in the associated RMS map (just take the value at that pixel).
for i in range(0, nin):
    mask = np.isnan(jointab[f'int_flux_{i+1}'])

    rmsmap = hdus[i].replace('_comp_warp-corrected_wbeam.fits', '_warp_rms.fits')
    rmshdu = fits.open(rmsmap)
    w = WCS(rmshdu[0].header, naxis=2)
    index = w.world_to_array_index(coords[mask])
    rms = rmshdu[0].data[index]
    jointab[f'local_rms_{i+1}'][mask] = rms

    img = hdus[i].replace('_comp_warp-corrected_wbeam.fits', '_warp.fits')
    imghdu = fits.open(img)
    w = WCS(imghdu[0].header, naxis=2)
    index = w.world_to_array_index(coords[mask])
    flux = imghdu[0].data[index]
    jointab[f'int_flux_{i+1}'][mask] = flux
    jointab[f'peak_flux_{i+1}'][mask] = flux

# We need some summary stats for the next bit to work
fluxes = np.empty((nin, len(jointab)), dtype='float32')
sigmas = np.empty((nin, len(jointab)), dtype='float32')
for i in range(0, nin):
    fluxes[i, :] = jointab[f"peak_flux_{i+1}"]
    sigmas[i, :] = jointab[f"local_rms_{i+1}"]
wfluxes = fluxes * sigmas**-2
weights = sigmas**-2
jointab['weighted_avg_peak_flux'] = np.nanmean(wfluxes, axis=0) / np.nanmean(weights, axis=0)
jointab['mean_peak_flux'] = np.nanmean(fluxes, axis=0)

innerterm = np.empty((nin, len(jointab)), dtype='float32')
fluxsq = np.empty((nin, len(jointab)), dtype='float32')
# Compute the inner terms for eta and var
for i in range(0, nin):
    innerterm[i, :] = ((jointab[f"peak_flux_{i+1}"] - jointab['weighted_avg_peak_flux'])**2) / (jointab[f"local_rms_{i+1}"]**2)
    fluxsq[i, :]  = jointab[f"peak_flux_{i+1}"]**2

# With our completed table, we can run an eta/V analysis.
jointab['eta'] = (1 / (nin - 1)) * np.nansum(innerterm, axis=0)
jointab['var'] = (1 / jointab['mean_peak_flux']) * np.sqrt((nin/(nin-1)) * (np.nanmean(fluxsq, axis=0) - jointab['mean_peak_flux']**2))

# Identify all internal matches within 2' using STILTS
jointab.write('tmp.fits', format='fits', overwrite=True)

os.system("stilts tmatch1 matcher='sky' values='mean_ra mean_dec' params=120 action=identify in=tmp.fits out=sparse.fits")

# Read that back in, and then write out the final table which should consist only of compact, isolated sources within the useful region of the images
jointab = Table(fits.open('sparse.fits')[1].data)
iso_mask = ~(jointab['GroupSize'] > 1)

final_mask = np.logical_and(np.logical_and(compact_mask, iso_mask), spatial_mask)
jointab[final_mask].write('isocompact_join_table.fits', format='fits', overwrite=True)
# Write everything while we're debugging
#jointab.write('modified_join_table.fits', format='fits', overwrite=True)

# Select interesting sources and make useful plots
eta_mask = jointab['eta'] > 5
var_mask = jointab['var'] > 0.05
transients_mask = np.logical_and(np.logical_and(final_mask, eta_mask), var_mask)

boxwidth = 200
for src in jointab[transients_mask]:
    # Get some source stats
    pos = SkyCoord(src['mean_ra'], src['mean_dec'], unit=(u.deg, u.deg), frame='fk5')
    radec = pos.to_string('hmsdms', sep=":", precision=1)
    ra, dec = pos.to_string('hmsdms', sep="", precision=0).split(" ")
    srcname = f"GPM_J{ra}{dec}"
    eta, var = src['eta'], src['var']
    fluxes = np.array([src[f'peak_flux_{i+1}'] for i in range(0, nin)])
    errs = np.array([src[f'local_rms_{i+1}'] for i in range(0, nin)])
    snr = src['mean_on_peak_flux']/src['mean_on_local_rms']
    # Make a light curve
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    ax.set_title(f"Pos={radec}   S/N={snr:3.1f}   $\eta =${eta:3.1f}   $V =${var:3.1f}")
    if snr > 10:
        color = 'red'
    else:
        color = 'grey'
    ax.errorbar(x=np.arange(0, nin), y=1000*fluxes, yerr=1000*errs, color=color)
    ax.set_ylabel("Flux density / mJy")
    ax.set_xlabel("Time / minutes")
    fig.savefig(f"{srcname}_lightcurve.png", bbox_inches="tight")
    # Make an animated gif
    for i in range(0, nin):
        img = hdus[i].replace('_comp_warp-corrected_wbeam.fits', '_warp.fits')
        imghdu = fits.open(img)
        w = WCS(imghdu[0].header, naxis=2)
        cutout = Cutout2D(imghdu[0].data, pos, (boxwidth, boxwidth), wcs = w)
        fig = plt.figure(figsize=(3,3))
        ax = fig.add_subplot(111)
        ax.imshow(cutout.data, vmin = np.nanmin(errs), vmax = np.nanmax(fluxes), origin="lower")
        plt.axis("off")
        plt.margins(x=0)
        plt.margins(y=0)
        fig.savefig("{0}_{1:02d}.png".format(srcname, i), bbox_inches="tight", dpi=100, pad_inches = 0)
        plt.close(fig)
    os.system(f"convert {srcname}_??.png {srcname}_animation.gif")
    madefiles = glob(f"{srcname}_??.png")
    for f in madefiles:
        os.remove(f)
