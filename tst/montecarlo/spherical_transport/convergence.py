#! /usr/bin/env python

"""
Read in athena results to judge convergence of general pusher through
spherical grid as stepsize decreases
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
from os import system

# Athena++ modules
import athena_mc_list as mclist
from athena_mc_list import photons

def write_athinput(iseed,nphot,step,file='athinput.sphtran',generalpusher=True):
    """
    Write the remainder of the athinput file for convergence test
    """

    outfile = open(file,'w')
    outfile.write("<comment>\n")
    outfile.write("problem   =  Test movement through spherical grid using general pusher or spherical pusher\n")
    outfile.write("reference =\n")
    outfile.write("configure = --prob=mc_sph_tran -mc --coord=spherical_polar\n")
    outfile.write("\n")
    outfile.write("<job>\n")
    outfile.write("problem_id = sphtran   # problem ID: basename of output filenames\n")
    outfile.write("\n")
    outfile.write("<output1>\n")
    outfile.write("file_type  = phlist\n")
    outfile.write("relativistic = true\n")
    outfile.write("nuser     = 1\n")
    outfile.write("\n")
    outfile.write("<mesh>\n")
    outfile.write("nx1        = 128       # Number of zones in X1-direction\n")
    outfile.write("x1min      = 0.01       # minimum value of X1\n")
    outfile.write("x1max      = 1.0e11     # maximum value of X1\n")
    outfile.write("ix1_mc_bc  = absorb  # inner-X1 boundary flag\n")
    outfile.write("ox1_mc_bc  = escape  # outer-X1 boundary flag\n")
    outfile.write("\n")
    outfile.write("nx2        = 64         # Number of zones in X2-direction\n")
    outfile.write("x2min      = 0.      # minimum value of X2\n")
    outfile.write("x2max      = 3.141592654   # maximum value of X2-\n")
    outfile.write("ix2_mc_bc  = polar  # inner-X2 boundary flag\n")
    outfile.write("ox2_mc_bc  = polar  # outer-X2 boundary flag\n")
    outfile.write("\n")
    outfile.write("nx3        = 64        # Number of zones in X3-direction\n")
    outfile.write("x3min      = 0.0       # minimum value of X3\n")
    outfile.write("x3max      = 6.2831853071 # maximum value of X3\n")
    outfile.write("ix3_mc_bc  = periodic  # inner-X3 boundary flag\n")
    outfile.write("ox3_mc_bc  = periodic  # outer-X3 boundary flag\n")
    outfile.write("\n")
    outfile.write("<hydro>\n")
    outfile.write("gamma = 1.666666666666667 # gamma = C_p/C_v\n")
    outfile.write("\n")
    outfile.write("<montecarlo>\n")
    outfile.write("nphot     = {:d}\n".format(nphot))
    outfile.write("iseed     = {:d}\n".format(iseed))
    outfile.write("scattering = none\n")
    outfile.write("emission   = freefree\n")
    outfile.write("absorption = none\n")
    outfile.write("polarized = false\n")
    if (generalpusher):
        outfile.write("stepsize = {:e}\n".format(step))
        outfile.write("general_pusher = true\n")
        outfile.write("varystep = true\n")
        outfile.write("checkmove = 100000000\n")
    elif (sphpol_pusher):
        outfile.write("stepsize = {:e}\n".format(step))
        outfile.write("sphpol_alt = true\n")
        outfile.write("varystep = true\n")
        outfile.write("checkmove = 100000000\n")
    outfile.write("\n<problem>\n")


def get_blackbody(temp,nu):
    c = 2.99792458e10
    kb = 1.380649e-16
    h = 6.62607015e-27
    return 2*h/c**2*nu**3/(np.exp(h*nu/(kb*temp)) - 1.0)

# Main function
def main(**kwargs):

    path =  kwargs.pop("path")
    athena_path = path+"/bin"

    mcranks = kwargs['mcranks']
    
    nstep = kwargs['nstep']
    nphot = kwargs['nphot']
    iseed = []
    if (nstep > 0):
        steps = [kwargs['step0']]
        ratio = kwargs['ratio']
        for i in range(nstep-1):
            steps.append(steps[i]/ratio)
    else:
        steps =[]
    print("steps: ",steps)
    iseed = kwargs['iseed']
    # Set up array to store norm for convergence evaluation
    error = np.zeros((nstep+1,2))

    # First do a run with standard cell-by-cell spherical pusher integration
    write_athinput(iseed,nphot,0.,sphpol_pusher=True, generalpusher=False)
    com="mpirun -np {:d} ".format(mcranks+1)+athena_path+"/athena -i athinput.sphtran"
    system(com)
    com="python ~/Documents/athena/vis/python/montecarlo/joinlists.py sphtran.out1 {:d} 0 0 -rm".format(mcranks+1)
    system(com)
    # read list as dict from infile
    phots = photons(mclist.read_list("sphtran.out1.list"))
    error[0,0] = 0.
    error[0,1] = np.average(phots.user[:,0])
    
    # Next loop over step size with general pusher prescription
    for i,step in enumerate(steps):            
        i += 1
        write_athinput(iseed+99*i,nphot,step)
        com="mpirun -np {:d} ".format(mcranks+1)+athena_path+"/athena -i athinput.sphtran"
        system(com)
        com="python ~/Documents/athena/vis/python/montecarlo/joinlists.py sphtran.out1 {:d} 0 0 -rm".format(mcranks+1)
        system(com)

        # read list as dict from infile
        phots = photons(mclist.read_list("sphtran.out1.list"))
        error[i,0] = step
        error[i,1] = np.average(phots.user[:,0])

    # save plot to outfile
    np.savetxt(kwargs['outfile'],error)
    
    if (nstep > 0):
        plt.plot(error[1:,0],error[1:,1],'+')
        plt.plot(error[1:,0],error[-1,1]*(error[-1,0]/error[1:,0])**(-1.))
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig("conv.pdf")

    print(error)

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('iseed',
        type = int,
        default = None,
        help='random seed')
    parser.add_argument('nphot',
        type = int,
        default = None,
        help='photon number')
    parser.add_argument('step0',
        type = float,
        default = None,
        help='largest step size')
    parser.add_argument('nstep',
        type = int,
        default = None,
        help='number of step sizez to evaluate')
    parser.add_argument('--ratio',
        type = float,
        default = 10.,
        help='ration by which to decrease step size')
    parser.add_argument('--mcranks',
        type = int,
        default = 5,
        help='mpi ranks to use')
    parser.add_argument('--path',
        default = "/home/bcm2vn/Documents/athena",
        help='path to Athena++ distribution')
    parser.add_argument('--outfile',
        default="conv.out",
        help='output filename for storing convergence rate')

    args = parser.parse_args()
    main(**vars(args))
