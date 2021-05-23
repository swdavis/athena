#! /usr/bin/env python

"""
Read in athena results to judge convergence towards blackbody emission
as photon number increases
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpmath import polylog
from os import system

# Athena++ modules
import athena_read

def bnu_int(en,tgas):
    kb = 1.380649e-16
    kt = kb * tgas
    x = np.exp(-en/(kt))
    if (en == 0.):
        # avoids evaluation of 0 * -inf in log
        return (-3.*kt*en**2*polylog(2,x)-6.*en*kt**2*polylog(3,x)
                -6*kt**3*polylog(4,x))
    else:
        return (en**3*np.log(1.-x)-3.*kt*en**2*polylog(2,x)
                -6.*en*kt**2*polylog(3,x)-6*kt**3*polylog(4,x))

def eta_kappap_ff(rho,tgas,emin,emax):
    """
    Computes the frequency integrated emmisivity for a gas of completely ionized
    H and He over an energy range extending from emin to emax, where theose energies
    are measured in eV
    """
    kb = 1.380649e-16
    mp = 1.67e-24
    me = 9.109e-28
    ec = 4.8032e-10
    c = 2.99792458e10
    h = 6.62607015e-27
    everg = 1.602176634e-12
    
    # efac accoutns for integration limits that are not 0, infinity
    efac = np.exp(-emin*everg/(kb*tgas))-np.exp(-emax*everg/(kb*tgas))
    ef0 = 32.*np.pi*ec**6/(3.*me*h*c**3)*(2.*np.pi*kb/(3.*me))**0.5*efac

    heabund = 0.09
    nh = rho/mp/(1.+4.*heabund)
    nhe = nh*heabund
    nel = nh+2.*nhe

    eta = ef0*nel*(nh+4.*nhe)*tgas**0.5

    # Integrates plank function yields sigma*T^4/pi in limit that
    # emin -> 0 and emax -> infinity
    planck_int = 2.*kb*tgas/h**3/c**2*(bnu_int(emax*everg,tgas)-
                                       bnu_int(emin*everg,tgas))
    kap = eta/(planck_int*4.*np.pi*rho)
    return eta, kap

def write_athinput(iseed,nphot,vel,frame,dens,tgas,emin,emax,file='athinput.mctest'):
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
    outfile.write("frame      = "+frame+"\n")
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
    outfile.write("emin      = {:e}\n".format(emin))
    outfile.write("emax      = {:e}\n".format(emax))
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
    outfile.write("temp     = {:e}\n".format(tgas))
    outfile.write("dens     = {:e}\n".format(dens))
    outfile.write("velocity = {:e}\n".format(vel))
    outfile.write("constdens = true\n")
    outfile.close()

# Main function
def main(**kwargs):

    athena_path="/home/swd8g/athena-swdavis/bin"
    infile = "athinput.mctest"

    mcranks = kwargs['mcranks']    
    nstep = kwargs['nstep']
    nphots = [kwargs['nmin']]
    step = kwargs['step']

    # input energy ranges and physical parameters
    tgas = kwargs.pop("tgas")
    dens = kwargs.pop("dens")

    # set limits
    kb = 1.380649e-16
    emin = 0.001*kb*tgas
    emax = 50.*kb*tgas

    c = 2.99792458e10
    h = 6.62607015e-27
    planck_int = 2.*kb*tgas/h**3/c**2*(bnu_int(emax,tgas)-
                                      bnu_int(emin,tgas))
    # convert to ev for athinput file
    everg = 1.602176634e-12
    emin /= everg
    emax /= everg

    er0 = 4.*np.pi*planck_int/c
    fr0 = 0.

    if ((kwargs['vel'] is not None) and (kwargs['frame'] == 'eulerian')):
        beta = kwargs['vel']
        gamma=1./(1.-beta**2)**0.5
        er = er0*gamma**2*(1+beta**2/3.)
        fr = er0*beta*c*4./3.*gamma**2
    else:
        er = er0
        fr = fr0

    iseed = []
    for i in range(nstep-1):
        nphots.append(nphots[i]*step)
    iseed = kwargs['iseed']
    # Set up array to store norm for convergence evaluation
    output = np.zeros((nstep,10))
    for i,nphot in enumerate(nphots):            
        write_athinput(iseed+99*i,nphot,kwargs['vel'],kwargs['frame'],dens,tgas,emin,emax,file=infile)
        com="mpirun -np {:d} ".format(mcranks+1)+athena_path+"/athena -i "+infile
        print com
        system(com)
        # read hdf5 output
        data = athena_read.athdf("MCTest.out1.00001.athdf",quantities=['Ermc','Frmc1','Frmc2','Frmc3','Cooling','kapjmc'])
        output[i,0] = float(nphot)
        output[i,1] = er0
        output[i,2] = er
        output[i,3] = fr
        output[i,4] = np.average(data['Ermc'])
        output[i,5] = np.average(data['Frmc1'])
        output[i,6] = np.average(data['Frmc2'])
        output[i,7] = np.average(data['Frmc3'])
        output[i,8] = np.average(data['Cooling'])
        output[i,9] = np.average(data['kapjmc'])/dens

    # save plot to outfile
    np.savetxt(kwargs['outfile'],output)
    
    plt.plot(output[:,0],abs(output[:,4]-output[:,1]),'r+')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("conv.pdf")
    print(output)
    eta, kappap = eta_kappap_ff(dens,tgas,emin,emax)
    print kappap,output[-1,9]/kappap

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
    parser.add_argument('--dens',
        type = float,
        default = 3.e-8,
        help='density')
    parser.add_argument('--tgas',
        type = float,
        default = 1.e5,
        help='gas temperature')
    parser.add_argument('--frame',
        default = "eulerian",
        help='boost velocity')
    parser.add_argument('--outfile',
        default="conv.out",
        help='output filename for storing convergence rate')

    args = parser.parse_args()
    main(**vars(args))
