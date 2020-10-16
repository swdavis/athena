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

def make_escape_time_dist(phots,nt,tmin,tmax,logt=True):
    """
    Makes escape time distribution (dict) from photon object
    """

    # Store dist as a dictionary
    times = {}

    # Store total number of photons for refernce
    times['ntot'] = phots.ntot

    # Create bins
    tfaces = mcspec.build_bins(tmin,tmax,nt,logt)
    times['nt'] = nt
    times['dist'] = np.zeros(nt)
    times['tfaces'] = tfaces

    # Get x bins
    tbins = mcspec.get_bins(phots.user[:,0],tfaces,nt,log=logt,uniform=True)
  
    for i in range(phots.nphot):
        if (tbins[i] >= 0):
            wght = phots.weight[i]
            times['dist'][tbins[i]] += wght
    #times['dist'] /= (tfaces[1:]-tfaces[:-1])
    times['dist'] /= float(phots.ntot)

    return times

def write_times(filename,times):

    nt = times['nt']
    out_arr = np.zeros((nt,3))
    out_arr[:,0] = times['tfaces'][:-1]
    out_arr[:,1] = times['tfaces'][1:]
    out_arr[:,2] = times['dist']
    np.savetxt(filename,out_arr,fmt='%1.6e')

# Main function
def main(**kwargs):

    # Filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    # Read photon list
    phlist = mclist.read_list(infile)
    phots = photons(phlist)

    logt = kwargs['log']

    times = make_escape_time_dist(phots,kwargs['nt'],kwargs['tmin'], \
                                  kwargs['tmax'],logt)

    # Write spectrum to file
    if outfile is None:
        outfile = infile.replace('.list','.tdist')

    write_times(outfile,times)

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
