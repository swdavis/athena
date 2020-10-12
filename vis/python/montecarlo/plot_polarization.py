#! /usr/bin/env python

"""
Plot athena++ spectrum
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Athena++ modules
import athena_mc_spec as mcspec
#import athena_mc_io as mcio

# Main function
def main(**kwargs):

    # Use latex labels
    plt.rc('text',usetex=True)
    plt.rc('font', **{'family' :"serif"})

    # filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    # read spectrum as dict from infile
    spectrum = mcspec.read_spectrum(infile)

    def imu_handler(imu):
        if (len(imu) > 1):
            # loop over all imu in the array
            slist = imu.strip(('[]')).split(",")
            ilist = [int(i) for i in slist]
        else:
            ilist = [int(imu)]
        return ilist

    # Get imu or imus for plotting different polar angles
    ilist = imu_handler(kwargs.pop('imu'))

    # plot spectrum
    if kwargs['yunit'] == 'specfrac':
        kwargs.pop('yunit')
        fig = plt.figure()
        gs = gridspec.GridSpec(5,1)
        ax1 = fig.add_subplot(gs[0:3,0])
        kwargs1 = dict(kwargs)
        kwargs1['yscale'] = 'log'
        kwargs1['yunit'] = 'nulnu'
        for imu in ilist:
            x, nu, ax1 = mcspec.plot_spectrum(spectrum,imu,ax1,**kwargs1)
            ax1.tick_params(labelbottom=False)
            ax1.set_xlabel("")
        ax2 = fig.add_subplot(gs[3:5,0])
        kwargs2 = dict(kwargs)
        kwargs2['yunit'] = 'frac'
        kwargs2.pop('ymin','ymax')
        kwargs2['ymin'] = 0
        for imu in ilist:
            x, nu, ax2 = mcspec.plot_polarization(spectrum,imu,ax2,**kwargs2)
            plt.tight_layout()
    else:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for imu in ilist:
            x, nu, ax = mcspec.plot_polarization(spectrum,imu,ax,**kwargs)
    
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
        default = 'linear',
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
        default='frac',
        help='variable to be used for y axis: frac, angle')
    parser.add_argument('--ploterr',
        action='store_true',
        help='plot intensity with error bar')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for spectrum')

    args = parser.parse_args()
    main(**vars(args))
