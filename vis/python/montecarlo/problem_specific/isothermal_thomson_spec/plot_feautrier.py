#! /usr/bin/env python

"""
Plot athena++ spectrum
"""

# python standard modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from scipy import interpolate

# Athena++ modules
import athena_mc as mcspec
import feautrier as feaut

def interp_feaut(mu0,mu,varin):
    """
    Interpolate feautrier solution in angle
    """
    nnu = len(varin[:,0])
    varout  = np.zeros(nnu)
    for i in range(nnu):
        varout[i] = np.interp(mu0,mu,varin[i,:])
    return varout

# Main function
def main(**kwargs):

    # Use latex labels
    #plt.rc('text',usetex=True)
    #plt.rc('font', **{'family' :"serif"})

    # filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    # read spectrum as dict from infile
    spectrum = mcspec.read_spectrum(infile)
    # compute nu, assuming spectrum stored as ev
    xfaces = spectrum['xfaces']
    h = 6.6262e-27
    everg = 1.6021772e-12
    nu = 0.5*(xfaces[1:]+xfaces[:-1])*everg/h
    mumid = 0.5*(spectrum['mufaces'][1:]+spectrum['mufaces'][:-1])

    fnorm = kwargs.pop('fnorm')
    ffile = kwargs.pop('ffile')
    if ffile is None:
        # Compute feautrier solution
        feaut.transfer(nd=128,na=32,nu=nu)
        ffile = "feautrier.out"
    # read in feautier solution
    nuf,muf,intensf,polf = feaut.read_feautrier(ffile)
    evf = nuf*h/everg

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

    # plot bottom spectra
    fig = plt.figure()
    gs = gridspec.GridSpec(5,1)
    ax1 = fig.add_subplot(gs[0:3,0])
    kwargs['yscale'] = 'log'
    yunitt = kwargs.pop('yunittop')
    yunitb = kwargs.pop('yunitbot')
    iphi = kwargs.pop("iphi")
    plterr = kwargs.pop("ploterr")
    xunit = kwargs.pop("xunit")

    for imu in ilist:
        x,y,yerr,xlabel,ylabel = \
        mcspec.plot_frequency(spectrum,imu,iphi=iphi,plterr=plterr,xunit=xunit,
                              yunit=yunitt)
        mcspec.make_plot(x,y,yerr=yerr,xlabel=xlabel,ylabel=ylabel,ax=ax1,**kwargs)
        ax1.tick_params(labelbottom=False)
        ax1.set_xlabel("")
        iinterp = interp_feaut(mumid[imu],muf,intensf)
        ax1.plot(evf,nuf*iinterp*fnorm)

    # plot top spectra
    ax2 = fig.add_subplot(gs[3:5,0])
    kwargs['yscale'] = 'linear'
    for imu in ilist:
        x,y,yerr,xlabel,ylabel = \
        mcspec.plot_frequency(spectrum,imu,iphi=iphi,plterr=plterr,xunit=xunit,
                              yunit=yunitb)
        mcspec.make_plot(x,y,yerr=yerr,xlabel=xlabel,ylabel=ylabel,ax=ax2,**kwargs)
        pinterp = interp_feaut(mumid[imu],muf,polf)
        ax2.plot(evf,pinterp)

    plt.tight_layout()

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
    parser.add_argument('--yunittop',
        default='nulnu',
        help='variable to be used for top y axis: frac, angle')
    parser.add_argument('--yunitbot',
        default='q',
        help='variable to be used for botom y axis: frac, angle')
    parser.add_argument('--ploterr',
        action='store_true',
        help='plot intensity with error bar')
    parser.add_argument('--fnorm',
        type = float,
        default = 1.,
        help='feautrier area normalization')
    parser.add_argument('--ffile',
        default = None,
        help='feautrier input file')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for spectrum')

    args = parser.parse_args()
    main(**vars(args))
