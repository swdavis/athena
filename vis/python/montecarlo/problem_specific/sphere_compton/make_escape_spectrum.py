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


def make_histogram_path(phots,x,n,pmin,pmax,logp=True):
    """
    Makes histogram (dict) of average escape time from photon object
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
    bins = mcspec.get_bins(x,faces,n,log=logp,uniform=True)
    number = np.zeros(n)
    for i,bin in enumerate(bins):
        if (bin >= 0):
            data['dist'][bin] += phots.weight[i] * phots.user[i,0]
            number[bin] += phots.weight[i]
   
    data['dist'] /= number
    print data['dist']
    return data

def make_histogram_abs(phots,x,n,pmin,pmax,logp=True):
    """
    Makes histogram (dict) of average escape time from photon object
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
    bins = mcspec.get_bins(x,faces,n,log=logp,uniform=True)
    number = np.zeros(n)
    for i,bin in enumerate(bins):
        if (bin >= 0):
            tauabs = -np.log(phots.weight[i])
            #data['dist'][bin] += phots.weight[i] * (tauabs/phots.user[i,0])
            #number[bin] += phots.weight[i]
            data['dist'][bin] += (tauabs/phots.user[i,0])
            number[bin] += 1.
   
    data['dist'] /= number
    print data['dist']
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

def energy_kt(phots,temp):
    kb = 1.3807e-16
    return phots.energy/(kb*temp)

def energy_kev(phots):
    keverg = 1.6021772e-9
    return phots.energy/keverg

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

    # Make histogram of the average absoprtion opacity and write to file
    hist_abs = make_histogram_abs(phots,x,kwargs['nx'],\
                                  kwargs['xmin'], kwargs['xmax'],logx)
    outfile = infile.replace('.list','.opac.dist')
    write_histogram(outfile,hist_abs)

    # Make histogram of the average path length and write to file
    hist_path = make_histogram_path(phots,x,kwargs['nx'],\
                                    kwargs['xmin'], kwargs['xmax'],logx)
    outfile = infile.replace('.list','.path.dist')
    write_histogram(outfile,hist_path)

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
    parser.add_argument('--xmin',
        type=float,
        default=0.01,
        help='minimum for time variable')
    parser.add_argument('--xmax',
        type=float,
        default = 100.,
        help='maximum for time variable')
    parser.add_argument('--linear',
        action='store_true',
        help='bins distributed linearly')
    #parser.add_argument('--outfile',
    #    default=None,
    #    help='output filename for escape spectrum')

    args = parser.parse_args()
    main(**vars(args))
