#! /usr/bin/env python

"""
Read in photon list and create escape time distribution

"""

# standard python modules
import argparse
import numpy as np

# athena++ modules
import athena_mc_spec as mcspec
import athena_mc_list as mclist
from athena_mc_list import photons

def make_histogram(phots,func,n,pmin,pmax,logp=True):
    """
    Makes simple histogram (dict) from photon object
    """

    # Store dist as a dictionary
    data = {}

    # Store total number of photons for refernce
    data['ntot'] = phots.ntot

    # Create bins
    faces = mcspec.build_bins(pmin,pmax,n,logp)
    data['n'] = n
    data['dist'] = np.zeros(n)
    data['faces'] = faces

    # Get x bins
    bins = mcspec.get_bins(func(phots),faces,n,log=logp,uniform=True)
  
    for i,bin in enumerate(bins):
        if (bin >= 0):
            data['dist'][bin] += phots.weight[i]
    data['dist'] /= np.sum(data['dist']*(faces[1:]-faces[:-1]))

    return data

def write_histogram(filename,data):
    """
    Write simple histogram to file as ascii table
    """
    n = data['n']
    out_arr = np.zeros((n,3))
    out_arr[:,0] = data['faces'][:-1]
    out_arr[:,1] = data['faces'][1:]
    out_arr[:,2] = data['dist']
    np.savetxt(filename,out_arr,fmt='%1.6e')


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

    radii = make_histogram(phots,radius,kwargs['nr'],kwargs['rmin'], \
                           kwargs['rmax'],logr)

    # Write spectrum to file
    if outfile is None:
        outfile = infile.replace('.list','.rdist')

    write_histogram(outfile,radii)

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
