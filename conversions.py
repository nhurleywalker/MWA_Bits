#!/usr/bin/env python

import numpy as np
import math

# Useful conversion functions

def km2pc(D):
    return D/(3.086e13)

def pc2km(D):
    return D*(3.086e13)

def sec2yr(T):
    return T/(365.25*24*60*60)

def yr2sec(T):
    return T*(365.25*24*60*60)

def Jy2mJy(S):
    return S*1000

def mJy2Jy(S):
    return S/1000

def deg2AS(angle):
    return angle*3600

def AM2AS(angle):
    return angle*60.

def AM2deg(angle):
    return angle/60.

def WL2freq(WL):
    return 3.e8/WL

def freq2WL(nu):
    return 3.e8/nu

def m2cm(length):
    return length*100

def cm2m(length):
    return length/100

def mK2K(T):
    return T/1000

#https://science.nrao.edu/facilities/vla/proposing/TBconv
def T2PFD(T, nu, bmaj, bmin):
    ''' Convert a brightness temperature to peak flux density. \
        Expecting T in K, nu in Hz, bmaj and bmin in deg; \
        Will convert to mJy, cm, arcsec." '''
    PFD = T * deg2AS(bmaj) * deg2AS(bmin) /     \
       ( 1.36 * (m2cm(freq2WL(nu)))**2 )
    return mJy2Jy(PFD)

def S2SB(S, a, b):
    ''' Convert a flux density in Jy and the radii subtended \
        (in arcmin) to a surface brightness in W m^-2 Hz^-1 sr&-1.'''
    SB = S / (math.pi * np.radians(AM2deg(a)) * np.radians(AM2deg(b)))
    SB *= 1.e-26
    return SB

