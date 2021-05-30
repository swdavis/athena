#! /usr/bin/env python

"""
Plot single athena++ spectrum.  Allows specification of multiple inclinations
to be plotted simultaneously.
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Athena++ modules
import athena_mc_spec as mcspec
import athena_read

# Main function
def main(**kwargs):


    # filenames for io
    infile = kwargs.pop('infile')

    # read spectrum as dict from infile
    spectrum = mcspec.read_spectrum(infile)

    h = 6.6262e-27
    everg = 1.6021772e-12
    nuf = spectrum['xfaces']*everg/h
    intens = spectrum['intensity'][0,0,0,:]
    dnu = nuf[1:]-nuf[:-1]
    efin = np.sum(dnu*intens)
    # read hdf5 output
    data = athena_read.athdf("MCTest.out1.00001.athdf",quantities=['Ermc','Cooling'])
    cool = np.average(data['Cooling'])*1.e33
    print cool,efin

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input photon list filename')


    args = parser.parse_args()
    main(**vars(args))
