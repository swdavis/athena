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
import athena_mc_list as mclist
from athena_mc_list import photons
import athena_read

# Main function
def main(**kwargs):


    # filenames for io
    infile = kwargs['infile']

    # read photon list from infile
    phlist = mclist.read_list(infile)
    phots = photons(phlist)
    efin = np.average(phots.weight*phots.energy)
    h = 6.6262e-27
    everg = 1.6021772e-12
    kb = 1.3806580e-16
    #einit = kb*kwargs['temp']*kwargs['x0']
    # read hdf5 output
    data = athena_read.athdf("MCTest.out1.00001.athdf",quantities=['Ermc','Cooling'])
    cool = np.average(data['Cooling'])*1.e33
    #print einit,efin,efin-einit
    #print -cool,(einit-efin)/cool
    print -cool,efin,-cool/efin
# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input photon list filename')
    #parser.add_argument('temp',
    #    type=float,
    #    help='temperature')
    #parser.add_argument('x0',
    #    type=float,
    #    help='x0')

    args = parser.parse_args()
    main(**vars(args))
