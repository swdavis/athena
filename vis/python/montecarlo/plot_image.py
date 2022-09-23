#! /usr/bin/env python

"""
Plot athena++ image
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Athena++ modules
import athena_mc as athenamc

# Main function
def main(**kwargs):

    # Use latex labels
    #plt.rc('text',usetex=True)
    #plt.rc('font', **{'family' :"serif"})

    # filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')
    if outfile is None:
        outfile = infile.replace('.img','.png')

    # read image
    image = athenamc.read_image(infile)

    # Set plot parameters
    iinc = kwargs.pop("iinc")

    # Set axis to be reused
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    athenamc.plot_image(image,iinc,ax=ax,**kwargs)
    #athenamc.plot_image_old(image,iinc,ax=ax)

    # save plot to outfile
    plt.savefig(outfile)
    plt.close()

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
                        help='input photon spectrum filename(s)')
    parser.add_argument('--type',
                        default = 'intensity',
                        help='variable to plot')
    parser.add_argument('--iinc',
                        default = 0,
                        type = int,
                        help='index of angle bin to plot')
    parser.add_argument('--xmin',
                        default = None,
                        type = float,
                        help='x-axis mimimum')
    parser.add_argument('--xmax',
                        default = None,
                        type = float,
                        help='x-axis maximum')
    parser.add_argument('--ymin',
                        default = None,
                        type = float,
                        help='y-axis mimimum')
    parser.add_argument('--ymax',
                        default = None,
                        type = float,
                        help='y-axis maximum')
    parser.add_argument('-c', '--colormap',
                        default='hot',
                        help='name of Matplotlib colormap to use instead of default; \
                              hot is default, twilight is good for cyclic variables \
                              like polarization angle')
    parser.add_argument('--vmin',
                        type=float,
                        default=None,
                        help='data value to correspond to colormap minimum; use \
                              --vmin=<val> if <val> has negative sign')
    parser.add_argument('--vmax',
                        type=float,
                        default=None,
                        help='data value to correspond to colormap maximum; use \
                              --vmax=<val> if <val> has negative sign')
    parser.add_argument('--vnorm',
                        action='store_true',
                        help='flag indicating that intensity should be normalized \
                              to maximum')
    parser.add_argument('--logc',
                        action='store_true',
                        help='flag indicating data should be colormapped logarithmically')
    parser.add_argument('-p', '--pvec',
                        action='store_true',
                        default=False,
                        help='flag indicating that polarization should be plotted')
    parser.add_argument('-a', '--average',
                        action='store_true',
                        default=False,
                        help='flag indicating that polarization should be averaged \
                              over steps')
    parser.add_argument('--step',
                        type=int,
                        default=4,
                        help='flag indicating that polarization should be normalized \
                              to maximum')
    parser.add_argument('--outfile',
                        default=None,
                        help='output filename for spectrum')

    args = parser.parse_args()
    main(**vars(args))
