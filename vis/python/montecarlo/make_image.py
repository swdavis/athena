#! /usr/bin/env python

"""
Read in photon list and create an image from the list.
"""

# python standard modules
import argparse
import numpy as np

# Athena++ modules
import athena_mc as athenamc
from athena_mc import Photons

try:
    import screen
except ModuleNotFoundError:
    pass

# Main function
def main(**kwargs):

    # Filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    # Read photon list
    phlist = athenamc.read_list(infile)
    phots = Photons(phlist)

    # Set parameters
    rcam = kwargs.pop('rcam')
    ninc = kwargs.pop('ninc')
    imin = kwargs.pop('imin')
    imax = kwargs.pop('imax')
    nx = kwargs.pop('nx')
    xmax = kwargs.pop('xmax')
    xmin = kwargs.pop('xmin')
    if (xmin is None):
        xmin = -xmax
    ny = kwargs.pop('ny')
    ymax = kwargs.pop('ymax')
    ymin = kwargs.pop('ymin')
    if (ymin is None):
        ymin = -ymax
    unit = kwargs.pop('unit')
    nen = kwargs.pop('nen')
    emin = kwargs.pop('emin')
    emax = kwargs.pop('emax')

    # check for screening function
    screen_name = kwargs.pop('screen')
    if screen_name != 'no_screen':
        screen_function = getattr(screen, screen_name)
        mask = screen_function(phots)
    else:
        mask = None

    # Create image
    print(f"Generating image from list.")
    image = athenamc.make_image_mc(phots,rcam,ninc,imin,imax,
                                   nen,emin,emax,
                                   nx,xmin,xmax,
                                   ny,ymin,ymax,
                                   unit=unit,mask=mask,**kwargs)

    # Write image to file
    if outfile is None:
        outfile = infile.replace('.list','.img')
    athenamc.write_image(outfile,image)

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input photon list filename')
    parser.add_argument('rcam',
        type=float,
        help='camera radius')
    parser.add_argument('xmax',
        type=float,
        help='maximum x')
    parser.add_argument('ymax',
        type=float,
        help='maximum y')
    parser.add_argument('--ninc',
        default=16,
        type=int,
        help='number of camera inclination bins')
    parser.add_argument('--imin',
        default=-1.,
        type=float,
        help='minimum for inclination variable')
    parser.add_argument('--imax',
        default=1.,
        type=float,
        help='maximum for inclination variable')
    parser.add_argument('--nx',
        type=int,
        default=16,
        help='number of x pixels per row')
    parser.add_argument('--ny',
        type=int,
        default=16,
        help='number of y pixels per column')
    parser.add_argument('--unit',
        default='cm',
        help='unit for x, y arrays')
    parser.add_argument('--xmin',
        type=float,
        default=None,
        help='minimum x')
    parser.add_argument('--ymin',
        type=float,
        default=None,
        help='minimum y')
    parser.add_argument('--nen',
        type=int,
        default=1,
        help='number of camera energy bins')
    parser.add_argument('--emin',
        type=float,
        default=1.e-300,
        help='minimum for energy variable')
    parser.add_argument('--emax',
        type=float,
        default=1.e300,
        help='maximum for energy variable')
    parser.add_argument('--screen',
        default = 'no_screen',
        help = 'name of screen function in screen.py file')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for spectrum')

    args = parser.parse_args()
    main(**vars(args))
