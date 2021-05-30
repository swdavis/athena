#! /usr/bin/env python

"""
Plot radius distribution

"""

# standard python modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# athena++ modules
import athena_mc_spec as mcspec

def rad_dist(x,t,nmax=1000):

    y = np.exp(-t)
    prob = 0.
    for n in range(1,nmax):
        prob = prob + x*np.sin(n*np.pi*x)*n*y**(n*n)
    return prob
 

# Main function
def main(**kwargs):

   # Use latex labels
    plt.rc('text',usetex=True)
    plt.rc('font', **{'family' :"serif"})

    # Filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    path = kwargs.pop('path')
    radius = kwargs.pop('radius')
    tau = kwargs.pop('tau')
    mfp = radius/tau
  

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Read radial distribution
    rdist = np.loadtxt(infile)
    rmid = 0.5*(rdist[:,0]+rdist[:,1])/radius
    pr = rdist[:,2]*radius

    path0 = mfp*tau**2*3./np.pi**2


    # Compute comparison function
    pr_comp = np.zeros(len(rmid))
    for i,r in enumerate(rmid):
        pr_comp[i] = rad_dist(r,path/path0)
    pr_comp /= np.sum(pr_comp*(rdist[:,1]-rdist[:,0])/radius)
    print pr_comp
    ax.set_ylabel(r"$P(r)$")
    ax.set_xlabel(r"$r$")
    # Plot radius hisotgram
    if (kwargs['notnorm']): 
        c = 2.99792e10
        #tmid *= time0/c
        ax.set_xlabel(r"$t \; \rm(s)$")
    ax.plot(rmid,pr,'.')
    ax.plot(rmid,pr_comp,':')

    # Set axis scales
    ax.set_xscale('linear')
    ax.set_yscale('log')


    # Write distribution to file
    if outfile is None:
        outfile = infile.replace('.rdist','.pdf')
        
    plt.savefig(outfile)
    plt.close()


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input filename')
    parser.add_argument('tau',
        type=float,
        help='optical depth of sphere')
    parser.add_argument('--radius',
        type=float,
        default = 1.e10,
        help='radius of sphere')
    parser.add_argument('--path',
        type=float,
        default = 2.e11,
        help='maximum path length')
    parser.add_argument('--notnorm',
        action='store_true',
        help='Sets r axis to physical units')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for escape time plot')

    args = parser.parse_args()
    main(**vars(args))
