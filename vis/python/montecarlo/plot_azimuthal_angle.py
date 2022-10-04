#! /usr/bin/env python

"""
Plot single athena++ spectrum as function of azimuthal angle.  Allows
specification of multiple frequencies to be plotted simultaneously.
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Athena++ modules
import athena_mc as athenamc


def ix_handler(ix):
    if ix == None:
        return [0]
    if ix == 'sum':
        return [ix]
    if (len(ix) > 1):
        # loop over all imu in the array
        slist = ix.strip(('[]')).split(",")
        ilist = [int(i) for i in slist]
    else:
        ilist = [ix]
    return ilist

def file_handler(infile):
    """
    Parse infile to deterimine file list
    """
    if (len(infile) > 1):
        # loop over all imu in the array
        flist = infile.strip(('[]')).split(",")
    else:
        flist = [infile]
    return flist

def plot_one(spectrum,ax,xunit,yunit,ix,imu,plterr,**kwargs):
    """
    Plot curves corresponding to single spectrum
    """
    #Convert xaxis, if needed
    if (xunit != spectrum['units']):
        athenamc.convert_xaxis(xunit,spectrum)

    # plot spectrum as function mu
    ilist = ix_handler(ix)
    for ix in ilist:
        x,y,yerr,xlabel,ylabel = athenamc.plot_phi(spectrum,ix,imu=imu,
                                 plterr=plterr,xunit=xunit,yunit=yunit)
        athenamc.make_plot(x,y,yerr=yerr,xlabel=xlabel,ylabel=ylabel,ax=ax,**kwargs)


# Main function
def main(**kwargs):

    # Use latex labels
    #plt.rc('text',usetex=True)
    #plt.rc('font', **{'family' :"serif"})

    # filenames for io
    infile = kwargs.pop('infile')
    files = file_handler(infile)
    outfile = kwargs.pop('outfile')
    if outfile is None:
        outfile = files[0].replace('.spec','.png')

    #  Set plot parameters
    plterr = kwargs.pop("ploterr")
    xunit = kwargs.pop("xunit")
    yunit = kwargs.pop("yunit")
    ix = kwargs.pop("ix")
    imu = kwargs.pop("imu")

    # Set axis to be reused
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # plot spectra from all infiles
    for file in files:
        # read spectrum as dict from infile
        spectrum = athenamc.read_spectrum(file)
        print("luminosity: ("+file+")",athenamc.get_luminosity(spectrum))

        # plot curves corresponding to this spectrum
        plot_one(spectrum,ax,xunit,yunit,ix,imu,plterr,**kwargs)

    # save plot to outfile
    plt.savefig(outfile)
    plt.close()

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input photon list filename')
    parser.add_argument('--ix',
        default = None,
        help='index of frequency bin to plot')
    parser.add_argument('--imu',
        default = 'sum',
        help='controls phi bin for plot')
    parser.add_argument('--xscale',
        default = 'linear',
        help='x-axis scale')
    parser.add_argument('--xmin',
        default = None,
        type = float,
        help='x-axis mimimum')
    parser.add_argument('--xmax',
        default = None,
        type = float,
        help='x-axis maximum')
    parser.add_argument('--yscale',
        default = 'linear',
        help='y-axis scale')
    parser.add_argument('--ymin',
        default = None,
        type = float,
        help='y-axis mimimum')
    parser.add_argument('--ymax',
        default = None,
        type = float,
        help='y-axis maximum')
    parser.add_argument('--xunit',
        default='mu',
        help='variable to be used for x axis: mu')
    parser.add_argument('--yunit',
        default='lnu',
                        help='variable to be used for y axis: lnu')
    parser.add_argument('--ploterr',
        action='store_true',
        help='plot intensity with error bar')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for spectrum')

    args = parser.parse_args()
    main(**vars(args))
