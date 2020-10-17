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

def esc_time(phots):
    return phots.user[:,0]

# Main function
def main(**kwargs):

    # Filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    # Read photon list
    phlist = mclist.read_list(infile)
    phots = photons(phlist)

    logt = kwargs['log']

    times = etdist.make_histogram(phots,esc_time,kwargs['nt'],\
                                  kwargs['tmin'], kwargs['tmax'],logt)

    # Write spectrum to file
    if outfile is None:
        outfile = infile.replace('.list','.tdist')

    etdist.write_histogram(outfile,times)

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input filename')
    parser.add_argument('nt',
        type=int,
        help='number of time bins')
    parser.add_argument('tmin',
        type=float,
        help='minimum for time variable')
    parser.add_argument('tmax',
        type=float,
        help='maximum for time variable')
    parser.add_argument('--log',
        action='store_true',
        help='bins distributed logarithmically')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for escape times')

    args = parser.parse_args()
    main(**vars(args))
