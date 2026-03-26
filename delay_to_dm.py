#!/usr/bin/env python

import numpy as np
import sys
import argparse

def delayToDM(delay, nu1, nu2):
    return abs(delay / (4.15*((nu1**-2) - (nu2**-2))))
    
def DMToDelay(DM, nu1, nu2):
    return 1.e-3*abs(4.15 * DM * ((nu1**-2) - (nu2**-2)))

def t_scat(DM, nu, a=-6.46, b=0.154, c=1.07, alpha=3.86):
    ''' implements t_scat / DM relationship of Bhat et al. 2004 '''
    ltscat = a + b*np.log10(DM) + c*(np.log10(DM)**2) - alpha*np.log10(nu/1000)
    return(10**ltscat)
           

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--delay', dest='delay', default=0.0, type=float,
                        help="Delay in seconds (set to convert delay to DM)")
    parser.add_argument('--DM', dest='DM', default=0.0, type=float,
                        help="DM in pc cm^-3 (set to convert DM to delay)")
    parser.add_argument('--nu1', dest='nu1', default=185, type=float,
                        help="lowest frequency in MHz")
    parser.add_argument('--nu2', dest='nu2', default=215, type=float,
                        help="highest frequency in MHz")
    args = parser.parse_args()

    if args.delay != 0.0:
        DM = delayToDM(1.e3*args.delay, args.nu1/1.e3, args.nu2/1.e3)
        tscat = t_scat(DM, (args.nu2+args.nu1)/2)
        print(f"DM = {DM:4.4f}, scattering timescale={tscat:4.4f} ms")
    elif args.DM != 0.0:
        tscat = t_scat(args.DM, (args.nu2+args.nu1)/2)
        delay = DMToDelay(args.DM, args.nu1/1.e3, args.nu2/1.e3)
        print(f"delay = {delay:4.4f} s, scattering timescale={tscat:4.4f} ms")
    else:
        print("Must set either delay or DM.")
