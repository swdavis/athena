#! /usr/bin/env python

"""
Read in athena results to judge convergence towards blackbody emission
as photon number increases
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
from os import system

# Athena++ modules
import athena_read

def write_athinput(iseed,nphot,vel,file='athinput.mctest'):
    """
    Write the remainder of the athinput file for convergence test
    """
    if (vel is None):
        vel = 0.
        boosts = False
    else:
        boosts = True        
    outfile = open(file,'w')
    outfile.write("<comment>\n")
    outfile.write("problem   =  Uniform periodic box\n")
    outfile.write("reference =\n")
    outfile.write("configure = --prob=mctest -mc\n")
    outfile.write("\n")
    outfile.write("<job>\n")
    outfile.write("problem_id = MCTest   # problem ID: basename of output filenames\n")
    outfile.write("\n")
    outfile.write("<output1>\n")
    outfile.write("file_type  = hdf5\n")
    outfile.write("dt         = 1.0e-11\n")
    outfile.write("variable   = mcmom\n")
    outfile.write("\n")
    outfile.write("<time>\n")
    outfile.write("cfl_number = 0.1      # The CFL Number\n")
    outfile.write("nlim       = 1        # cycle limit\n")
    outfile.write("tlim       = 1.0      # time limit\n")
    outfile.write("\n")
    outfile.write("<mesh>\n")
    outfile.write("nx1        = 16       # Number of zones in X1-direction\n")
    outfile.write("x1min      = -5.0e10       # minimum value of X1\n")
    outfile.write("x1max      = 5.0e10      # maximum value of X1\n")
    outfile.write("ix1_bc     = periodic  # inner-X1 boundary flag\n")
    outfile.write("ox1_bc     = periodic  # outer-X1 boundary flag\n")
    outfile.write("ix1_mc_bc     = periodic  # inner-X1 boundary flag\n")
    outfile.write("ox1_mc_bc     = periodic  # outer-X1 boundary flag\n")
    outfile.write("\n")
    outfile.write("nx2        = 16         # Number of zones in X2-direction\n")
    outfile.write("x2min      = -5.0e10      # minimum value of X2\n")
    outfile.write("x2max      = 5.0e10       # maximum value of X2\n")
    outfile.write("ix2_bc     = periodic  # inner-X2 boundary flag\n")
    outfile.write("ox2_bc     = periodic  # outer-X2 boundary flag\n")
    outfile.write("ix2_mc_bc  = periodic  # inner-X2 boundary flag\n")
    outfile.write("ox2_mc_bc  = periodic  # outer-X2 boundary flag\n")
    outfile.write("\n")
    outfile.write("nx3        = 16         # Number of zones in X3-direction\n")
    outfile.write("x3min      = -5.0e10      # minimum value of X2\n")
    outfile.write("x3max      = 5.0e10       # maximum value of X2\n")
    outfile.write("ix3_bc     = periodic  # inner-X3 boundary flag\n")
    outfile.write("ox3_bc     = periodic  # outer-X3 boundary flag\n")
    outfile.write("ix3_mc_bc  = periodic  # inner-X3 boundary flag\n")
    outfile.write("ox3_mc_bc  = periodic  # outer-X3 boundary flag\n")
    outfile.write("\n")
    outfile.write("<hydro>\n")
    outfile.write("gamma = 1.666666666666667 # gamma = C_p/C_v\n")
    outfile.write("iso_sound_speed = 1.0     # isothermal sound speed\n")
    outfile.write("\n")
    outfile.write("<montecarlo>\n")
    outfile.write("nphot     = {:d}\n".format(nphot))
    outfile.write("iseed     = {:d}\n".format(iseed))
    outfile.write("emin      = 0.01\n")
    outfile.write("emax      = 300.\n")
    outfile.write("scattering = none\n")
    outfile.write("emission   = freefree\n")
    outfile.write("absorption = freefree\n")
    outfile.write("polarized = false\n")
    if (boosts):
        outfile.write("boosts     = true\n")
    else:
        outfile.write("boosts     = false\n")
    outfile.write("\n")
    outfile.write("<problem>\n")
    outfile.write("temp     = 1.0e5\n")
    outfile.write("dens     = 3.0e-8\n")
    outfile.write("velocity = {:e}\n".format(vel))
    outfile.write("constdens = true\n")
    outfile.close()

# Main function
def main(**kwargs):

    athena_path="/home/swd8g/athena-swdavis/bin"
    infile = "athinput.mctest"
    er0 = 7.5646e5
    fr0 = 0.
    if (kwargs['vel'] is not None):
        c=2.99e10
        beta = kwargs['vel']
        gamma=1./(1.-beta**2)**0.5
        er = er0*gamma**2*(1+beta**2/3.)
        fr = er0*beta*2.99e10*4./3.*gamma**2
    else:
        er = er0
        fr = fr0

    mcranks = kwargs['mcranks']
    
    nstep = kwargs['nstep']
    nphots = [kwargs['nmin']]
    step = kwargs['step']
    iseed = []
    for i in range(nstep-1):
        nphots.append(nphots[i]*step)
    iseed = kwargs['iseed']
    # Set up array to store norm for convergence evaluation
    lnorm = np.zeros((nstep,5))
    for i,nphot in enumerate(nphots):            
        write_athinput(iseed+99*i,nphot,kwargs['vel'],file=infile)
        com="mpirun -np {:d} ".format(mcranks+1)+athena_path+"/athena -i "+infile
        print com
        system(com)
        # read hdf5 output
        data = athena_read.athdf("MCTest.out1.00001.athdf",quantities=['Ermc','Frmc1','Frmc2','Frmc3'])
        lnorm[i,0] = float(nphot)
        lnorm[i,1] = np.average(abs(data['Ermc']-er))
        lnorm[i,2] = np.average(abs(data['Frmc1']-fr0))
        lnorm[i,3] = np.average(abs(data['Frmc2']-fr0))
        lnorm[i,4] = np.average(abs(data['Frmc3']-fr))

    # save plot to outfile
    np.savetxt(kwargs['outfile'],lnorm)
    
    plt.plot(lnorm[:,0],lnorm[:,1],'r+')
    plt.plot(lnorm[:,0],lnorm[-1,1]*(lnorm[-1,0]/lnorm[:,0])**0.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("conv.pdf")
    print(er/er0,fr/er0)
    print(lnorm)


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
    parser.add_argument('--vel',
        type = float,
        default = None,
        help='boost velocity')
    parser.add_argument('--outfile',
        default="conv.out",
        help='output filename for storing convergence rate')

    args = parser.parse_args()
    main(**vars(args))
