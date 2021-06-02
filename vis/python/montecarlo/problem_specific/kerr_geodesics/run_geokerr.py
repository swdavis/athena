#! /usr/bin/env python

"""
Read in photon list, then initialize and run geokerr
"""

# standard python modules
import argparse
import numpy as np
import os

# athena++ modules
import athena_mc_list as mclist
from athena_mc_list import photons

# Main function
def main(**kwargs):

    # Filenames for io
    listfile = kwargs.pop('listfile')
    spin = float(kwargs.pop('spin'))
    nsteps = int(kwargs.pop('nsteps'))
    # Read photon list
    phlist = mclist.read_list(listfile)
    ph = photons(phlist)

    # Write mccomp.in file input for geokerr with initial parameters
    outfile = "mccomp.in"
    out = open(outfile, 'w')
    out.write("0\n0\n{:d}\n{:e}\n{:d}\n".format(ph.nphot,spin,nsteps))
    for i in range(ph.ntot):
        out.write("{:e} {:e} {:e} {:e} {:e} {:d} {:d}\n"
                  .format(ph.user[i,0],ph.user[i,1],ph.user[i,2],ph.user[i,5],
                          1./ph.x1[i],int(ph.user[i,3]),int(ph.user[i,4])))
    out.close()
    
    # Run geokerr with paramers specified in mccomp.in
    os.system("./geokerr < mccomp.in > mccomp.out");

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('listfile',
        help='input list filename')
    parser.add_argument('spin',
        help='black hole spin')
    parser.add_argument('--nsteps',
        type=int,
        default=500,
        help='number of steps')
    args = parser.parse_args()
    main(**vars(args))
