#! /usr/bin/env python

"""
Read in photon list and create escape time distribution

"""

# standard python modules
import argparse
import numpy as np

# athena++ modules
import escape_time as etdist
import athena_mc_list as mclist
from athena_mc_list import photons

def radius(phots):
    if (phots.coord == 'cartesian'):
        return np.sqrt(phots.x1**2+phots.x2**2+phots.x3**2)
    if (phots.coord == 'spherical_polar'):
        return phots.x1

# Main function
def main(**kwargs):

    # Filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    # Read photon list
    phlist = mclist.read_list(infile)
    phots = photons(phlist)

    logr = kwargs['log']

    print radius(phots)
    radii = etdist.make_histogram(phots,radius,kwargs['nr'],kwargs['rmin'], \
                                   kwargs['rmax'],logr)

    # Write spectrum to file
    if outfile is None:
        outfile = infile.replace('.list','.rdist')

    etdist.write_histogram(outfile,radii)

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input filename')
    parser.add_argument('nr',
        type=int,
        help='number of radial bins')
    parser.add_argument('rmin',
        type=float,
        help='minimum for radial variable')
    parser.add_argument('rmax',
        type=float,
        help='maximum for radial variable')
    parser.add_argument('--log',
        action='store_true',
        help='bins distributed logarithmically')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for escape times')

    args = parser.parse_args()
    main(**vars(args))
