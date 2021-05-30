"""
Support for creating escape time and radial distributions
"""

# standard python modules
import numpy as np

# athena++ modules
import athena_mc_spec as mcspec

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
    print data['dist']
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




