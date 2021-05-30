#! /usr/bin/env python

"""
Read in photon list and create a spectrum from the list.
"""

# python standard modules
import argparse
import numpy as np

# Athena++ modules
import athena_mc_spec as mcspec
import athena_mc_list as mclist
from athena_mc_list import photons

# Main function
def main(**kwargs):

    # Filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    # Read photon list
    phlist = mclist.read_list(infile)
    phots = photons(phlist)
    
    # Generate spectrum from phots
    nx = kwargs.pop('nx')
    xmax = kwargs.pop('xmax')
    xmin = kwargs.pop('xmin')
    if (xmin is None):
        xmin = -xmax
    temp = kwargs.pop('temp')
    kb = 1.3806504e-16;
    mass = 1.660538782e-24;
    vth = np.sqrt( 2 * kb * temp / mass );
    nu0 = 2.468e15;
    c = 2.997924589e10;
    dopw = nu0 * vth / c;
    xmin = xmin * dopw + nu0
    xmax = xmax * dopw + nu0
    print kwargs
    spectrum = mcspec.make_spectrum(phots,nx,xmin,xmax,'nu',logx=False,**kwargs)
    spectrum['xfaces'] = (spectrum['xfaces'] - nu0) / dopw
 
    # Write spectrum to file
    if outfile is None:
        outfile = infile.replace('.list','.spec')
    mcspec.write_spectrum(outfile,spectrum)

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input photon list filename')
    parser.add_argument('nx',
        type=int,
        help='number of x bins')
    parser.add_argument('temp',
        type=float,
        help='temperature for computing Doppler width')
    parser.add_argument('xmax',
        type=float,
        help='maximum for x variable')
    parser.add_argument('--xmin',
        default = None,
        type=float,
        help='minimum for x variable')
    parser.add_argument('--nmu',
        type=int,
        default=1,
        help='number of cos(theta) bins')
    parser.add_argument('--mumin',
        type=float,
        default=0.,
        help='minimum cos(theta)')
    parser.add_argument('--mumax',
        type=float,
        default=1.,
        help='maximum cos(theta)')
    parser.add_argument('--nphi',
        type=int,
        default=1,
        help='number of phi bins')
    parser.add_argument('--phimin',
        type=float,
        default=0.,
        help='minimum phi')
    parser.add_argument('--phimax',
        type=float,
        default=2.*np.pi,
        help='maximum phi')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for spectrum')
    parser.add_argument('--yerror',
        action='store_true',
        help='compute intensity errors')
    
    args = parser.parse_args()
    main(**vars(args))

