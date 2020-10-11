#! /usr/bin/env python

"""
Read in athena results to judge convergence towards blackbody emission
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
from os import system

# Athena++ modules
import athena_mc_spec as mcspec

def write_athinput(iseed,nphot,mcranks,file='athinput_params'):
    """
    Write the remainder of the athinput file for convergence test
    """

    outfile = open(file,'w')
    outfile.write("<montecarlo>\n")
    outfile.write("mcranks   = {:d}\n".format(mcranks))
    outfile.write("nphot     = {:d}\n".format(nphot))
    outfile.write("iseed     = {:d}\n".format(iseed))
    outfile.write("emin      = 1.0\n")
    outfile.write("emax      = 200.0\n")
    outfile.write("scattering = none\n")
    outfile.write("emission   = freefree\n")
    outfile.write("absorption = freefree\n")
    outfile.write("polarized = false\n")
    outfile.write("\n<problem>\n")
    outfile.write("temp     = 1.0e5\n")
    outfile.write("taumin   = 1.0e-3\n")
    outfile.write("taumax   = 1.0e4\n")
    outfile.close()

def get_blackbody(temp,nu):
    c = 2.9979246e+10
    kb = 1.3806580e-16
    h = 6.6262e-27
    return 2*h/c**2*nu**3/(np.exp(h*nu/(kb*temp)) - 1.0)

# Main function
def main(**kwargs):

    athena_path="/home/swd8g/athena-swdavis/bin"
    area = 1.e22


    nphots = [kwargs['nmin']]
    step = kwargs['step']
    iseed = []
    for i in range(kwargs['nstep']-1):
        nphots.append(nphots[i]*step)
    iseed = kwargs['iseed']
    lnorm = np.zeros((2,len(nphots)))
    for i,nphot in enumerate(nphots):            
        mcranks = kwargs['mcranks']
        write_athinput(iseed+99*i,nphot,mcranks)
        system("cat athinput_base athinput_params > athinput.mctest")
        com="mpirun -np {:d} ".format(mcranks+1)+athena_path+"/athena -i athinput.mctest"
        print com
        system(com)
        # read spectrum as dict from infile
        spectrum = mcspec.read_spectrum("MCTest.out1.00000.spec")
        xfaces = spectrum['xfaces']
        h = 6.6262e-27
        everg = 1.6021772e-12
        nu = 0.5*(xfaces[1:]+xfaces[:-1])*everg/h
        ibb = get_blackbody(1.e5,nu)*area
        # average spectrum over phi
        intensity = spectrum['intensity']
        intensity = np.sum(intensity[0,:,:,:],axis=0)/float(spectrum['nphi'])
        norm = 0.
        n = 0
        for j in range(spectrum['nmu']):
            for k in range(spectrum['nx']):
                norm += abs(intensity[j,k]-ibb[k])/ibb[k]
                n += 1
        plt.plot(nu,intensity[0,:])
        plt.plot(nu,ibb)
        lnorm[0,i] = float(nphot)
        lnorm[1,i] = norm/float(n)
        

    # save plot to outfile
    np.savetxt(kwargs['outfile'],lnorm)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    print lnorm


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('iseed',
        type = int,
        default = None,
        help='random seed')
    parser.add_argument('nmin',
        type = int,
        default = None,
        help='minimum photon number')
    parser.add_argument('nstep',
        type = int,
        default = None,
        help='number of steps')
    parser.add_argument('step',
        type = int,
        default = 2,
        help='factor by which to increase nphot')
    parser.add_argument('--mcranks',
        type = int,
        default = 10,
        help='mpi ranks to use')
    parser.add_argument('--outfile',
        default="conv.out",
        help='output filename for storing convergence rate')

    args = parser.parse_args()
    main(**vars(args))
