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


def make_histogram_2d(phots,nx,np,path0):
    """
    Makes histogram (dict) of average escape time from photon object based
    on path lenth and final energy
    """

    # Store dist as a dictionary
    data = {}

    # Store total number of photons for refernce
    data['ntot'] = phots.ntot

    # set physical parameters
    rad0 = 1.e10
    tau = 60.

    # Create bins
    mfp = rad0/tau
    path0 = mfp*tau**2*3./np.pi**2
    xfaces = mcspec.build_bins(0.01,100.,nx,True)
    pfaces = mcspec.build_bins(0.01,100.,np,True)
    data['nx'] = nx
    data['np'] = np
    data['hist'] = np.zeros((nx,np))
    data['xfaces'] = xfaces
    data['pfaces'] = pfaces

    x = energy_kt(phots)
    binsx = mcspec.get_bins(x,xfaces,nx,log=True,uniform=True)
    binsp = mcspec.get_bins(phots.user[:,0]/path0,pfaces,np,log=True,uniform=True)

    for i,bin in enumerate(bins):
        if (bin >= 0):
            tauabs = -np.log(phots.weight[i])
            data['dist'][bin] += (tauabs/phots.user[i,0])
            number[bin] += 1.

    return data


def write_histogram_2d(filename,data):
    """
    Write simple histogram to file as ascii table
    """
    n = data['n']
    out_arr = np.zeros((n,3))
    out_arr[:,0] = data['faces'][:-1]
    out_arr[:,1] = data['faces'][1:]
    out_arr[:,2] = data['dist']
    np.savetxt(filename,out_arr,fmt='%1.6e')

def energy_kt(phots,temp):
    kb = 1.3807e-16
    return phots.energy/(kb*temp)


# Main function
def main(**kwargs):

    # Filenames for io
    infile = kwargs.pop('infile')
    #outfile = kwargs.pop('outfile')

    # Read photon list
    phlist = mclist.read_list(infile)
    phots = photons(phlist)

    logx = not kwargs['linear']
    temp = kwargs.pop('temp')
    x = energy_kt(phots,temp)

    # set physical parameters
    rad0 = 1.e10
    tau = 60.

    # Create bins
    mfp = rad0/tau
    path0 = mfp*tau**2*3./np.pi**2

    # Make histogram of the average absoprtion opacity and write to file
    hist_2d = make_histogram_2d(phots,x,kwargs['nx'],kwargs['np'], path0)
    outfile = infile.replace('.list','.2d.dist')
    write_histogram_2d(outfile,hist_2d)


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input filename')
    parser.add_argument('temp',
        type=float,
        help='temperature')
    parser.add_argument('--nx',
        type=int,
        default = 100,
        help='number of energy bins')
    parser.add_argument('--np',
        type=int,
        default = 100,
        help='number of path bins')
    parser.add_argument('--linear',
        action='store_true',
        help='bins distributed linearly')
    #parser.add_argument('--outfile',
    #    default=None,
    #    help='output filename for escape spectrum')

    args = parser.parse_args()
    main(**vars(args))
