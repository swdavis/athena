#! /usr/bin/env python

"""
Script for running severarla tests of the Athena++ monte carlo module
that can be performed by compiling with
> python configure.py -mc --prob=mctest -mpi -hdf5
> make clean
> make all
This script require mpi libraries are installed and in the user's path
This is a temporary method until a more formal regression test suite
is written
"""

# python standard modules
import argparse
import numpy as np
from os import system,chdir,getcwd

# Main function
def main(**kwargs):

    iseed = 121500
    nphot = 10000
    nstep = 3

    # set path to athena distribution tests
    path = kwargs.pop('path')
    tstpath = path+"/tst/montecarlo"

    # change to rundir for running convergence tests
    curdir = getcwd()
    rundir = curdir+'/rundir'
    system("mkdir -p "+rundir)
    chdir("rundir")

    # Run convergence test towards blackbody spectrum
    print("Running blackbody spectral convergence test.")
    system("python "+tstpath+"/absorption_spectrum/convergence.py {:d} {:d} {:d} 10".format(iseed,nphot,nstep)+" --path "+path)
    conv_abs = np.loadtxt("conv.out")
    system("rm conv.out")

    # Run convergence test of polarized Thomson scattering towards a feautrier
    # solution -- note that convergence will stall at higher number due to
    # slight mismatch between calculation methods
    print("Running polarized thomson scattering converence test.")
    system("cp "+tstpath+"/thomson_polarized_spectrum/feautrier.out .")
    system("python "+tstpath+"/thomson_polarized_spectrum/convergence.py {:d} {:d} {:d} 10 --ffile=feautrier.out".format(iseed,nphot*10,nstep)+" --path "+path)
    conv_scat = np.loadtxt("conv.out")
    system("rm conv.out feautrier.out")

    # Run convergence test for estimate of radiation field without boosts
    print("Running radiation moments test without boosts.")
    system("cp "+path+"/vis/python/athena_read.py .")
    system("python "+tstpath+"/boosts/convergence.py {:d} {:d} {:d} 10".format(iseed,nphot,nstep)+" --path "+path) 
    conv_boost_off = np.loadtxt("conv.out")
    system("rm conv.out")

    # Run convergence test for estimate of radiation field in Eulerian frame
    # with boosts
    print("Running radiation moments test with boosts in Eulerian frame.")
    system("python "+tstpath+"/boosts/convergence.py {:d} {:d} {:d} 10 --vel=0.9".format(iseed,nphot,nstep)+" --path "+path) 
    conv_boost_eul = np.loadtxt("conv.out")
    system("rm conv.out")

    # Run convergence test for estimate of radiation field in comoving frame
    # with boosts
    print("Running radiation moments test with boosts in comoving frame.")
    system("python "+tstpath+"/boosts/convergence.py {:d} {:d} {:d} 10 --vel=0.9 --frame=comoving".format(iseed,nphot,nstep)) 
    conv_boost_com = np.loadtxt("conv.out")
    system("rm conv.out")

    # Run convergence test for cooling functions without compton scattering
    print("Running cooling estimate test without Compton scattering.")
    system("python "+tstpath+"/cartesian_cooling/convergence.py {:d} {:d} {:d} 10 --noscat".format(iseed,nphot,nstep+2)+" --path "+path) 
    conv_cool_abs = np.loadtxt("conv.out")
    system("rm conv.out")

    # Run convergence test for cooling functions with compton scattering
    print("Running cooling estimate test with Compton scattering.")
    system("cp "+tstpath+"/cartesian_cooling/comptontable.out .")
    system("python "+tstpath+"/cartesian_cooling/convergence.py {:d} {:d} {:d} 10".format(iseed,nphot,nstep)+" --path "+path) 
    conv_cool_sct = np.loadtxt("conv.out")
    system("rm conv.out comptontable.out")

    chdir(curdir)
    system("rm -rf rundir")

    np.set_printoptions(linewidth=100) 
    print("Results from the absorption spectrum test: ")
    for i in range(nstep):
        print(conv_abs[i,:])
    print("Results from the scattering spectrum test: ")
    for i in range(nstep):
        print(conv_scat[i,:])
    print("Results for Er, Fr without boosts: ")
    for i in range(nstep):
        print(conv_boost_off[i,:])
    print("Results for Er, Fr with boosts in Eulerian frame: ")
    for i in range(nstep):
        print(conv_boost_eul[i,:])
    print("Results for Er, Fr with boosts in comoving frame: ")
    for i in range(nstep):
        print(conv_boost_eul[i,:])
    print("Results for cooling without scattering: ")
    for i in range(nstep+2):
        print(conv_cool_abs[i,:])
    print("Results for cooling with scattering: ")
    for i in range(nstep):
        print(conv_cool_sct[i,:])

    #print("Results from the cooling: ")
    #for i in range(nstep+1):
    #    print(conv_sphtrans[i,:])

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path',
        default = None,
        help='path to Athena++ distribution')
    parser.add_argument('--mcranks',
        type = int,
        default = 10,
        help='mpi ranks to use')

    args = parser.parse_args()
    main(**vars(args))
