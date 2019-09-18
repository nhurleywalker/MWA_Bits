#!/usr/bin/env python

# Simulate some large and varying chromatic background with a narrow-band signal
# Use a moving median to remove it

import matplotlib.pyplot as plt

import numpy as np

nchans = 100

x = np.linspace(0, nchans, nchans)

# Large background: a couple of different sin waves with different wavelengths

sin1 = np.random.rand()*np.sin(x/20.)
sin2 = np.random.rand()*np.sin((x-5)/10.)
sin3 = np.random.rand()*np.sin((x+3)/5.)
sin4 = np.random.rand()*np.sin((x-10)/30.)

# Add some noise
noise = np.random.rand(nchans)

bkg = sin1+sin2+sin3+sin4+noise

# Add a large signal
bkg[int(nchans*np.random.rand())] += 5.

# Calculate the moving average
binwidth = 20

med = np.zeros(nchans)

for i in range(binwidth/2, nchans - binwidth/2):
    med[i] = np.median(bkg[i-binwidth/2:i+binwidth/2])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x, bkg, label="signal + background")
ax.plot(x, med, label="moving median")
ax.plot(x, bkg-med, label="background-subtracted signal")
ax.legend()
fig.savefig("test.png")
