#!/usr/bin/env python

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg') # So does not use display
import matplotlib.pylab as plt
from astropy.modeling import models, fitting
from astropy.modeling.fitting import LevMarLSQFitter

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

sig2fwhm = (2.*np.sqrt(2.*np.log(2.)))

N = 20 # of observations
try:
    sigma = float(sys.argv[1])
except IndexError:
    sigma = 120. # of ionospheric variations
fwhm_psf = 3.*60. # of PSF
sig_psf = fwhm_psf / sig2fwhm


rangex = np.arange(-8*sig_psf,8*sig_psf,sig_psf/5)

# Randomised offsets
offsets = sigma*np.random.randn(N)

gaussians=np.empty((N,len(rangex)))
i = 0
for offset in offsets:
    gaussians[i] = gaussian(rangex,offset,sig_psf)
    i+=1

total_gaussian = np.sum(gaussians,axis=0)/N

# Fit a gaussian to the final total to see what the new sigma is
g_init = models.Gaussian1D(amplitude=1, mean=0., stddev=sig_psf)
fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, rangex, total_gaussian)

print "Offset of new PSF = {0:3.2f} arcsec ; FWHM of new PSF = {1:3.2f} arcsec (old was {2:3.2f} arcsec; blur factor = {3:3.2f})".format(g.mean.value,g.stddev.value*sig2fwhm,fwhm_psf,g.stddev.value*sig2fwhm/fwhm_psf)

# Make a nice plot

fig=plt.figure(figsize=(5, 5))
ax = plt.gca()
#ax.set_yscale('log')
ax.plot(rangex,total_gaussian,label="Blurred PSF")
ax.plot(rangex,gaussian(rangex,0,sig_psf),label="Original PSF")
ax.plot(rangex,g(rangex),label="Gaussian fit")
ax.plot([g.mean.value,g.mean.value],[0,1])
ax.set_xlim(-6.*sig_psf,6.*sig_psf)
ax.set_ylim(0,1.3)
ax.legend()
plt.savefig("randomised_gaussian.png")

