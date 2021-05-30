#! /usr/bin/env python

"""
Plot single athena++ spectrum and compare with Kompaneets solution
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import kv
from scipy.special import gamma
from scipy.special import hyp1f1
from scipy.special import hyperu

# Athena++ modules
import athena_mc_spec as mcspec

# whitteker function
def whitteker(alpha,x):
    return np.exp(-x/2)*x**(alpha+2.)*hyperu(alpha,2*alpha+4,x)

# unnormalized Sunyaev Titarchuk solution
def st1980(x,x0,alpha):

    if (x <= x0):
        return (alpha+3.)*gamma(alpha+1.)/gamma(2.*alpha+4.)*np.exp(x0/2)/x0**3*hyp1f1(alpha,2.*alpha+4.,x)*x**(3.+alpha)*np.exp(-x)*whitteker(alpha,x0)
    else:
        return (alpha+3.)*gamma(alpha+1.)/gamma(2.*alpha+4.)*hyp1f1(alpha,2*alpha+4,x0)*x0**(alpha-1.)*x*np.exp(-x/2.)*whitteker(alpha,x)

# normalized Sunyaev Titarchuk solution
def st1980_norm(nu,x0,tcomp,tau):
    """
    Returns the ST1980 Kompaneets solution for a photon initially at
    x0 in an isothermal sphere of temperature tcomp and optical
    depth tau
    """

    h = 6.6262e-27
    c = 2.9979246e+10
    kb = 1.3806580e-16
    me = 9.1093897e-28

    normcomp = h*x0/np.pi
    y = kb*tcomp/(me*c**2)*tau*tau/3.
    gamma = (me*c**2)/kb/tcomp*np.pi**2/3./(tau+2./3.)**2
    alpha = (9./4.+gamma)**0.5-1.5
    #print tcomp,tau,y,alpha,normcomp
    x = h*nu/kb/tcomp
    nnu = len(nu)
    icomp = np.zeros(nnu)
    for i,xi in enumerate(x):
        icomp[i] = normcomp*st1980(xi,x0,alpha)
        
    return icomp

# Main function
def main(**kwargs):


    # Use latex labels
    plt.rc('text',usetex=True)
    plt.rc('font', **{'family' :"serif"})

    # filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')
    normt = kwargs.pop('normt')
    x0 = kwargs.pop('x0')
    temp = kwargs.pop('temp')
    tau = kwargs.pop('tau')


    # read spectrum as dict from infile
    spectrum = mcspec.read_spectrum(infile)

    # plot spectrum
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    x, nu, ax = mcspec.plot_spectrum(spectrum,0,ax,**kwargs)
 
    print("lumin: ",mcspec.get_luminosity(spectrum))

    icomp = normt*st1980_norm(nu,x0,temp,tau)
 
    plt.plot(x,nu*icomp)

    # save plot to outfile
    if outfile is None:
        outfile = infile.replace('.spec','.pdf')
    plt.savefig(outfile)
    plt.close()

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input photon list filename')
    parser.add_argument('temp',
        type = float,
        help='temperature of sphere')
    parser.add_argument('tau',
        type = float,
        help='optical depth of sphere')
    parser.add_argument('x0',
        type = float,
        help='initial normalized photon energy')
    parser.add_argument('--xscale',
        default = 'log',
        help='x-axis scale')
    parser.add_argument('--xmin',
        default = None,
        help='x-axis mimimum')
    parser.add_argument('--xmax',
        default = None,
        help='x-axis maximum')
    parser.add_argument('--yscale',
        default = 'log',
        help='y-axis scale')
    parser.add_argument('--ymin',
        default = None,
        help='y-axis mimimum')
    parser.add_argument('--ymax',
        default = None,
        help='y-axis maximum')
    parser.add_argument('--normt',
        type = float,
        default = 1.,
        help='normalization factor for kompaneets')
    parser.add_argument('--xunit',
        default='kev',
        help='variable to be used for x axis: ev, kev, nu, lambda')
    parser.add_argument('--yunit',
        default='nulnu',
        help='variable to be used for y axis: nulnu, lnu, counts')
    parser.add_argument('--ploterr',
        action='store_true',
        help='plot intensity with error bar')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for spectrum')
    parser.add_argument('--istokes',
        type=int,
        default=0,
        help='component of stokes vector for plotting')

    args = parser.parse_args()
    main(**vars(args))
