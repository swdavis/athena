#! /usr/bin/env python

"""
Read in photon list and create a spectrum from the list.
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
    """
    Wrapper for running the make_spectrum() function in athena_mc.py. Parameters
    of the spectrum are specified at the command line (with argparse) and passed
    via kwargs.
    """
    # Filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    # spectrum parameters
    nx = kwargs.pop('nx')
    xmin = kwargs.pop('xmin')
    xmax = kwargs.pop('xmax')
    logx = not kwargs.pop('linearx')

    # check for screening function
    screen_name = kwargs.pop('screen')
    if screen_name != 'no_screen':
        screen_function = getattr(screen, screen_name)

    # Read photon list
    reader = athenamc.read_list_generator(infile)
    result = next(reader)  # Get header
    header = result['header']

    spectrum = {}
    nchunk = 0
    list_lum = 0.
    for result in reader:

        phlist = header.copy()
        phlist['list'] = result['chunk']
        phlist['length'] = result['length']

        if kwargs['calclum']:
            list_lum += athenamc.get_luminosity_list(phlist)

        if (nchunk % 20) == 0:
            print(f"Generating spectrum: {result['remaining']} samples remain.")
        nchunk += 1

        # Creat photon object for current chunk
        phots = Photons(phlist)

        if screen_name != 'no_screen':
            mask = screen_function(phots)
        else:
            mask = None

        # Make spectrum from photon phots
        spec = athenamc.make_spectrum(phots, nx,xmin, xmax,logx=logx,
                                      mask=mask, **kwargs)

        spectrum = athenamc.add_spectra(spectrum, spec)

        if result['done']:
            break

    if (kwargs['calclum']):
        print("List luminosity: ", list_lum)
    # Write spectrum to file
    if outfile is None:
        outfile = infile.replace('.list','.spec')
    athenamc.write_spectrum(outfile,spectrum)

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help = 'input photon list filename')
    parser.add_argument('nx',
        type = int,
        help = 'number of x bins')
    parser.add_argument('xmin',
        type = float,
        help = 'minimum for x variable')
    parser.add_argument('xmax',
        type = float,
        help = 'maximum for x variable')
    parser.add_argument('--nmu',
        type = int,
        default = 1,
        help = 'number of cos(theta) bins')
    parser.add_argument('--mumin',
        type = float,
        default = 0.,
        help = 'minimum cos polar angle')
    parser.add_argument('--mumax',
        type = float,
        default = 1.,
        help = 'maximum cos polar angle')
    parser.add_argument('--nphi',
        type = int,
        default = 1,
        help = 'number of phi bins')
    parser.add_argument('--phimin',
        type = float,
        default = 0.,
        help = 'minimum phi')
    parser.add_argument('--phimax',
        type = float,
        default = 2.*np.pi,
        help = 'maximum phi')
    parser.add_argument('--anglebin',
        default = 'cartesian',
        help = 'controls binning: cartesian or spherical')
    parser.add_argument('--xaxis',
        default = 'ev',
        help = 'variable to be used for x axis: ev, kev, nu, lambda')
    parser.add_argument('-linearx',
        action = 'store_true',
        help = 'bins energies distributed logarithmically')
    parser.add_argument('-calclum',
        action = 'store_true',
        help = 'calculate luminosity directrly from list')
    parser.add_argument('--screen',
        default = 'no_screen',
        help = 'name of screen function in screen.py file')
    parser.add_argument('--outfile',
        default = None,
        help = 'output filename for spectrum')
    parser.add_argument('-yerror',
        action = 'store_true',
        help = 'compute intensity errors')

    args = parser.parse_args()
    main(**vars(args))