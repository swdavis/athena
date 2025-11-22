#! /usr/bin/env python
"""
Read in photon list(s) and create a spectrum from the list(s).
"""
# python standard modules
import argparse
import numpy as np
import glob
import re
from collections import defaultdict
# Athena++ modules
import athena_mc as athenamc
from athena_mc import Photons
try:
    import screen
except ModuleNotFoundError:
    pass

def process_single_file(infile, screen_function, nx, xmin, xmax, logx, **kwargs):
    """
    Process a single input file and return the spectrum.
    """
    # Read photon list
    reader = athenamc.read_list_generator(infile)
    result = next(reader)  # Get header
    header = result['header']
    spectrum = {}
    nchunk = 0

    print(f"Processing list file: {infile}")
    for result in reader:
        phlist = header.copy()
        phlist['list'] = result['chunk']
        phlist['length'] = result['length']

        if (nchunk % 20) == 0:
            print(f"  Generating spectrum: {result['remaining']} samples remain.")
        nchunk += 1

        # Create photon object for current chunk
        phots = Photons(phlist)
        if screen_function is not None:
            mask = screen_function(phots)
        else:
            mask = None

        # Make spectrum from photon phots
        spec = athenamc.make_spectrum(phots, nx, xmin, xmax, logx=logx,
                                      mask=mask, **kwargs)
        spectrum = athenamc.add_spectra(spectrum, spec)

        if result['done']:
            break

    return spectrum


def main(**kwargs):
    """
    Wrapper for running the make_spectrum() function in athena_mc.py. Parameters
    of the spectrum are specified at the command line (with argparse) and passed
    via kwargs.
    """
    # Filenames for io
    infiles = kwargs.pop('infiles')
    outfile = kwargs.pop('outfile')
    combine = kwargs.pop('combine')
    nouts = kwargs.pop('nouts')

    # spectrum parameters
    nx = kwargs.pop('nx')
    xmin = kwargs.pop('xmin')
    xmax = kwargs.pop('xmax')
    logx = not kwargs.pop('linearx')
    print(nx, xmin, xmax)
    # check for screening function
    screen_name = kwargs.pop('screen')
    if screen_name != 'no_screen':
        screen_function = getattr(screen, screen_name)
    else:
        screen_function = None

    all_files = []
    by_timestep = defaultdict(list)
    for pattern in infiles:
        matched = glob.glob(pattern)
        if matched:
            all_files.extend(matched)
        else:
            all_files.append(pattern)  # Keep as-is if no match

    # search for different proc numbers and output numbers
    proc_numbers = set()
    file_numbers = set()
    pattern = r'proc(\d+)\.(\d+)\.list'
    base_name = None
    for filename in all_files:
        match = re.search(pattern, filename)
        if match:
            proc_num = int(match.group(1))  # First captured group
            file_num = int(match.group(2))  # Second captured group
            proc_numbers.add(proc_num)
            file_numbers.add(file_num)
            by_timestep[file_num].append(filename)
        match = re.match(r'(.+?)\.proc\d+\.\d+\.list', filename)
        if match:
            base_name = match.group(1)

    nouts_inf = len(file_numbers)
    nprocs = len(proc_numbers)

    if not all_files:
        raise ValueError("No input files found!")

    print(f"Found {len(all_files)} input file(s) to process")
    print(f"Inferring {nprocs} processes and {nouts_inf} outputs")
    
    if nouts is not None:
        if nouts_inf != nouts:
            print(f"Error: number of outputs specified ({nouts}) "
                  f"does not match inferred number ({nouts_inf})")
            return
        if combine:
            combined_spectrum = {}
        for i in range(nouts):
            # Combine each output into a single spectrum
            out_spectrum = {}
            for infile in by_timestep[i]:
                spectrum = process_single_file(infile, screen_function, nx, xmin,
                                               xmax, logx, **kwargs)
                out_spectrum = athenamc.add_spectra(out_spectrum, spectrum)
            if combine:
                combined_spectrum = athenamc.add_spectra(combined_spectrum,
                                                        out_spectrum, 'time')
            else:
                # Write output spectrum for this timestep
                if outfile is None:
                    outname = basename+f'.{i:04d}.spec'
                else:
                    outname = outfile.replace('.spec', f'_{i:04d}.spec')
                athenamc.write_spectrum(outname, out_spectrum)
                print(f"Spectrum for output {i} written to: {outname}")
        if combine:
            # Write combined spectrum to file
            if outfile is None:
                outfile = base_name + '_combined.spec'
            athenamc.write_spectrum(outfile, combined_spectrum)
            print(f"Combined spectrum written to: {outfile}")
    elif combine:
        # Combine all files into a single spectrum
        combined_spectrum = {}
        for infile in all_files:
            spectrum = process_single_file(infile, screen_function, nx, xmin,
                                          xmax, logx, **kwargs)
            combined_spectrum = athenamc.add_spectra(combined_spectrum, spectrum)

        # Write combined spectrum to file
        if outfile is None:
            outfile = all_files[0].replace('.list', '_combined.spec')
        athenamc.write_spectrum(outfile, combined_spectrum)
        print(f"Combined spectrum written to: {outfile}")
    else:
        # Process each file separately
        for infile in all_files:
            spectrum = process_single_file(infile, screen_function, nx, xmin,
                                          xmax, logx, **kwargs)

            # Write spectrum to file
            if outfile is None:
                out = infile.replace('.list', '.spec')
            else:
                out = outfile

            athenamc.write_spectrum(out, spectrum)
            print(f"Spectrum written to: {out}")


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('nx',
        type=int,
        help='number of x bins')
    parser.add_argument('xmin',
        type=float,
        help='minimum for x variable')
    parser.add_argument('xmax',
        type=float,
        help='maximum for x variable')
    parser.add_argument('infiles',
        nargs='+',  # Accept one or more files
        help='input photon list filename(s) - supports wildcards like "*.list"')
    parser.add_argument('--nouts',
        type=int,
        default=None,
        help='number of outputs timesteps to process')   
    parser.add_argument('--nmu',
        type=int,
        default=1,
        help='number of cos(theta) bins')
    parser.add_argument('--mumin',
        type=float,
        default=0.,
        help='minimum cos polar angle')
    parser.add_argument('--mumax',
        type=float,
        default=1.,
        help='maximum cos polar angle')
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
    parser.add_argument('--xaxis',
        default='ev',
        help='variable to be used for x axis: ev, kev, nu, lambda')
    parser.add_argument('-linearx',
        action='store_true',
        help='bins energies distributed logarithmically')
    parser.add_argument('-calclum',
        action='store_true',
        help='calculate luminosity directly from list')
    parser.add_argument('--screen',
        default='no_screen',
        help='name of screen function in screen.py file')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for spectrum (used as base for multiple files)')
    parser.add_argument('-yerror',
        action='store_true',
        help='compute intensity errors')
    parser.add_argument('-combine',
                        action='store_true',
        help='combine all input files into a single output spectrum')
    args = parser.parse_args()
    main(**vars(args))
