#! /usr/bin/env python

"""
Plot single athena++ spectrum.  Allows specification of multiple inclinations
to be plotted simultaneously.
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Athena++ modules
import athena_mc_spec as mcspec

# Main function
def main(**kwargs):

    # Get blackbody parameters
    bbtemp = kwargs.pop('bbtemp')
    bbnorm = kwargs.pop('bbnorm')

    # Use latex labels
    plt.rc('text',usetex=True)
    plt.rc('font', **{'family' :"serif"})

    # filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    # read spectrum as dict from infile
    spectrum = mcspec.read_spectrum(infile)

    def imu_handler(imu):
        if imu == 'sum':
            return [imu]
        if (len(imu) > 1):
            # loop over all imu in the array
            slist = imu.strip(('[]')).split(",")
            ilist = [int(i) for i in slist]
        else:
            ilist = [imu]
        return ilist

    # plot spectrum
    ilist = imu_handler(kwargs.pop('imu'))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for imu in ilist:
        x, nu, ax = mcspec.plot_spectrum(spectrum,imu,ax,**kwargs)
 
    print("lumin: ",mcspec.get_luminosity(spectrum))
    if bbtemp is not None:
        if bbnorm is None:
            bbnorm = 1.
        else:
            bbtemp = float(bbtemp)
            bbnorm = float(bbnorm)
            c = 2.99792458e10
            kb = 1.380649e-16
            h = 6.62607015e-27
            ybb = bbnorm*2*h/c**2*nu**3/(np.exp(h*nu/(kb*bbtemp)) - 1.0)
            if kwargs['yunit'] == 'nulnu':
                plt.plot(x,ybb*nu,linestyle=':')
            if kwargs['yunit'] == 'lnu':
                plt.plot(x,ybb,linestyle=':')
            if kwargs['yunit'] == 'counts':
                plt.plot(x,ybb/(h*nu),linestyle=':')

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
    parser.add_argument('--imu',
        default = None,
        help='index of angle bin to plot')
    parser.add_argument('--iphi',
        default = 'ave',
        help='controls phi bin for plot')
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
    parser.add_argument('--bbtemp',
        default = None,
        help='blackbody temperature')
    parser.add_argument('--bbnorm',
        default = None,
        help='blackbody normalization')

    args = parser.parse_args()
    main(**vars(args))
