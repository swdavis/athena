#! /usr/bin/env python

"""
Plot escape time distribution

"""

# standard python modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# athena++ modules
import athena_mc_spec as mcspec

def freefree(rho,temp,energy):

    ffnrm = 3.692146e8
    heabund = 0.09
    mp = 1.6726e-24 
    h =  6.6262e-27
    kb = 1.3807e-16

    nhii = rho / (mp*(1.+4.*heabund))
    ne = (1. + 2.*heabund) * nhii

    nu = energy / h
    ehnu = np.exp(-energy / (kb * temp) )

    aff = ffnrm/np.sqrt(temp)/nu**3
    return ne * nhii * aff * (1. - ehnu);

# Main function
def main(**kwargs):

   # Use latex labels
    plt.rc('text',usetex=True)
    plt.rc('font', **{'family' :"serif"})

    # Filenames for io
    infiles = kwargs.pop('infiles')
    print infiles
    outfile = kwargs.pop('outfile')
    type = kwargs.pop('type')

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)



    kb = 1.3807e-16
    keverg = 1.6021772e-9

    # Read escape time distributions from infiles
    for infile in infiles:
        dist = np.loadtxt(infile)
        xmid = 0.5*(dist[:,0]+dist[:,1])
        var = dist[:,2]
        ax.plot(xmid,var,'.')

    if (type == 'opac'):
        ax.set_ylabel(r"$\rm Opacity$")
    if (type == 'path'):
        ax.set_ylabel(r"$\rm Path$")
    ax.set_xlabel(r"$E/(kT)$")

    # Set axis scales
    ax.set_xscale('log')
    ax.set_yscale('log')

    if (type == 'opac'):
        x0 = kwargs.pop('x0')
        tau = kwargs.pop('tau')
        temp = kwargs.pop('temp')
        rad = kwargs.pop('rad')
        heabund = 0.09
        mp = 1.6726e-24 
        sigmat = 6.65248e-25
        kappaes = sigmat * (1. + 2.*heabund) / (mp * (1.+4.*heabund) )
        rho = tau / (kappaes * rad)
        eerg = xmid*kb*temp
        opacff = freefree(rho,temp,eerg)
        einit = x0*kb*temp
        fac = kwargs.pop('fac')
        einit *= fac
        opac0 = freefree(rho,temp,einit)
        #opac = (opacff*opac0)**0.5
        opac = (opac0-opacff)/(2.*np.log(eerg/einit))
        ax.plot(xmid,opac,'xb')
        ax.plot(xmid,opacff,'+r')

    # Write distribution to file
    if outfile is None:
        outfile = 'out.pdf'
        #outfile = infile.replace('.dist','.pdf')
        
    plt.savefig(outfile)
    plt.close()


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infiles',
        nargs = '*',
        help='input filenames')
    parser.add_argument('--type',
        help='type of spectral data')
    parser.add_argument('--x0',
        type=float,
        help='initial normalized energy')
    parser.add_argument('--tau',
        type=float,
        help='optical depth')
    parser.add_argument('--temp',
        default = 3.e6,
        type=float,
        help='temperature of sphere')
    parser.add_argument('--rad',
        default = 1.e10,
        type=float,
        help='radius of sphere')
    parser.add_argument('--fac',
        default = 1.,
        type=float,
        help='e factor')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for plot')

    args = parser.parse_args()
    main(**vars(args))
