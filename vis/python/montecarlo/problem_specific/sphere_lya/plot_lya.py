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

def plot_spec_lya(spectrum,imu,ax=None,iphi='ave',
                  ploterr=True,yscale='log',xmin=None,
                  xmax=None,ymin=None,ymax=None,**kwargs):
    """
    Plot spectrum. Assumes saving, etc. are performed by the calling function
    """

    if (ax is None):
        # Create figure, axis and assume a single plot window
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    # Set up x axis as bin midpoints
    xfaces = spectrum['xfaces']
    x = 0.5*(xfaces[1:]+xfaces[:-1])

    # Initialize x label
    xlabel = r"$x$"
    ax.set_xlabel(xlabel)

    # Check if error requested and stored
    if ploterr:
        if spectrum['yerror'] != "true":
            print("Warning: error requested but not computed in spectrum.\n")
            ploterr = False


    intensity = spectrum['intensity'][0,:,:,:]
    print intensity
    if ploterr:
        errors = spectrum['errors'][0,:,:,:]

    # Selection for azimuthal angle
    if ((iphi == 'ave') or (iphi == 'sum')):
        norm = 1./float(spectrum['nphi'])
        if iphi == 'sum':
            norm *= 2.*np.pi
        intensity = np.sum(intensity,axis=0)*norm
        if ploterr:
            errors = np.sqrt(np.sum((errors)**2,axis=0))*norm
    else:
        iphi = int(iphi)
        intensity = intensity[iphi,:,:]
        if ploterr:
            errors = errors[iphi,:,:]
    # Selection for polar angle
    if imu == 'sum':
        nmu = spectrum['nmu']
        mumid = 0.5*(spectrum['mufaces'][1:]+spectrum['mufaces'][:-1])
        intensity = np.dot(mumid,intensity)/nmu
        if ploterr:
            errors = np.sqrt(np.dot((mumid)**2,(errors)**2))/nmu
    else:
        imu = int(imu)
        intensity = intensity[imu,:]
        if ploterr:
            errors = errors[imu,:]

    # Set y, yerr, and ylabel according to input units
    ylabel = r"$P(x)$"
    h = 6.6262e-27
    y = intensity
    if ploterr:
        yerr = errors
    ax.set_ylabel(ylabel)
    dnu = (spectrum['xfaces'][1:]-spectrum['xfaces'][:-1])
    yint = np.sum(dnu*y)
    y /= yint
    if (ploterr):
        ax.errorbar(x,y,yerr=yerr,fmt='.',**kwargs)
    else:
        ax.plot(x,y,'.',**kwargs)
    nx = spectrum['nx']
    out = np.zeros((nx,2))
    out[:,0] = x
    out[:,1] = y
    np.savetxt("spec.out",out)
    # Set axis scales
    ax.set_xscale('linear')
    ax.set_yscale(yscale)

    # (re)Set plot ranges
    left,right = ax.get_xlim()
    if xmin is not None:
        left=float(xmin)
    if xmax is not None:
        right=float(xmax)
    ax.set_xlim([left,right])
    left,right = ax.get_ylim()
    if ymin is not None:
        left=float(ymin)
    if ymax is not None:
        right=float(ymax)
    ax.set_ylim([left,right])
    
    # Return x and nu to facilitate evaluation of comparison functions
    # that may plotted by calling function.  Return ax to enable
    # further call to plot on the same axis
    return x, ax


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
        if imu == None:
            return [0]
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
        x, ax = plot_spec_lya(spectrum,imu,ax,**kwargs)
 
    
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
    #parser.add_argument('--xscale',
    #    default = 'linear',
    #    help='x-axis scale')
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
    parser.add_argument('--ploterr',
        action='store_true',
        help='plot intensity with error bar')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for spectrum')

    args = parser.parse_args()
    main(**vars(args))
