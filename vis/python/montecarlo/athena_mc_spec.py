"""
Support for manipulating and plotting Monte Carlo spectra
"""

# standard python modules
import numpy as np
import struct
import matplotlib.pyplot as plt

# athena++ modules
#import athena_mc_list as mclist
from athena_mc_list import photons

def write_spectrum(filename,spectrum):
    """
    Writes spectrum to output file
    """

    # Open outfile
    outfile = open(filename, 'w')

    nx = spectrum['nx']
    nmu = spectrum['nmu']
    nphi = spectrum['nphi']

    # Write header information
    outfile.write("nx={:d}\n".format(nx))
    outfile.write("nmu={:d}\n".format(nmu))
    outfile.write("nphi={:d}\n".format(nphi))
    outfile.write("ntot={:d}\n".format(spectrum['ntot']))
    outfile.write("nintens={:d}\n".format(spectrum['nintens']))
    outfile.write("units="+spectrum['xaxis']+"\n")
    outfile.write("polarized="+spectrum['polarized']+"\n")
    outfile.write("yerror="+spectrum['yerror']+"\n")
    outfile.close()

    # Write binfaces
    outfile = open(filename, 'ab')
    myfmt='>'+'d'*(nx+1)
    bin=struct.pack(myfmt,*(spectrum['xfaces']))
    outfile.write(bin)
    myfmt='>'+'d'*(nmu+1)
    bin=struct.pack(myfmt,*(spectrum['mufaces']))
    outfile.write(bin)
    myfmt='>'+'d'*(nphi+1)
    bin=struct.pack(myfmt,*(spectrum['phifaces']))
    outfile.write(bin)
    # Write data
    nelements = (spectrum['nintens']*nx*nmu*nphi)
    myfmt='>'+'d'*nelements
    bin=struct.pack(myfmt,*(spectrum['intensity'].reshape(nelements)))
    outfile.write(bin)
    if (spectrum['yerror'] == 'true'):
        bin=struct.pack(myfmt,*(spectrum['errors'].reshape(nelements)))
        outfile.write(bin)
    outfile.close()

def read_spectrum(filename):
    """
    Read spectrum and return as a dictionary
    """

    # Read raw data
    with open(filename, 'rb') as data_file:
        raw_data = data_file.read()
    raw_data_ascii = raw_data.decode('ascii', 'replace')

    spectrum = {}
    current_index = 0

    # Function for skipping though the file
    def skip_string(expected_string):
        expected_string_len = len(expected_string)
        if raw_data_ascii[current_index:current_index+expected_string_len] != expected_string:
            raise RuntimeError('File not formatted as expected')
        return current_index+expected_string_len

    current_index = skip_string("nx=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nx'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("nmu=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nmu'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("nphi=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nphi'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("ntot=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['ntot'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("nintens=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['nintens'] = list(map(int,raw_data_ascii[current_index:end_of_line_index].split(' ')))[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("units=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['units'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("polarized=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['polarized'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1
    current_index = skip_string("yerror=")
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    spectrum['yerror'] = raw_data_ascii[current_index:end_of_line_index].split(' ')[0]
    current_index = end_of_line_index + 1

    # Read in faces
    nx = spectrum['nx']
    format_string = '>' + 'd'*(nx+1)
    begin_index = current_index
    end_index = begin_index + 8*(nx+1)
    spectrum['xfaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    nmu = spectrum['nmu']
    format_string = '>' + 'd'*(nmu+1)
    begin_index = end_index
    end_index = begin_index + 8*(nmu+1)
    spectrum['mufaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    nphi = spectrum['nphi']
    format_string = '>' + 'd'*(nphi+1)
    begin_index = end_index
    end_index = begin_index + 8*(nphi+1)
    spectrum['phifaces'] = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))

    # Read intensities
    nintens = spectrum['nintens']
    nelements = nintens*nx*nmu*nphi
    format_string = '>' + 'd'*nelements
    begin_index = end_index
    end_index = begin_index + 8*nelements
    vals = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    spectrum['intensity'] = vals.reshape((nintens,nphi,nmu,nx))
    if (spectrum['yerror'] == 'true'):
        begin_index = end_index
        end_index = begin_index + 8*nelements
        vals = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
        spectrum['errors'] = vals.reshape((nintens,nphi,nmu,nx))
    return spectrum

# Retrun x for desired units
def convert_xaxis(newunit,spectrum):

    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10

    baseunit = spectrum['units']
    xfaces = spectrum['xfaces']
    if (baseunit == 'kev'):
        nu = xfaces*1000.*everg/h
    if (baseunit == 'ev'):
        nu = xfaces*everg/h
    if (baseunit == 'nu'):
        nu = xfaces
    if (baseunit == 'lambda'):
        nu = c/1.e8/xfaces
    if (newunit == 'kev'):
        spectrum['xfaces'] = nu*h/(everg*1000.)
    if (newunit == 'ev'):
        spectrum['xfaces'] = nu*h/everg
    if (newunit == 'nu'):
        spectrum['xfaces'] = nu
    if (newunit == 'lambda'):
        spectrum['xfaces'] = c/nu*1.e8
    spectrum['units'] = newunit

# Return nu
def get_frequency(xunit,xfaces):

    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10

    if (xunit == 'kev'):
        nu = 0.5*(xfaces[1:]+xfaces[:-1])*1000.*everg/h
    if (xunit == 'ev'):
        nu = 0.5*(xfaces[1:]+xfaces[:-1])*everg/h
    if (xunit == 'nu'):
        nu = 0.5*(xfaces[1:]+xfaces[:-1])
    if (xunit == 'lambda'):
        nu = 0.5*(1./xfaces[:-1]+1./xfaces[1:])*c/1.e8
    return nu

def plot_spectrum(spectrum,imu,ax=None,iphi='ave',xunit='kev',yunit='nulnu',
                  ploterr=True,xscale='log',yscale='log',istokes=0,xmin=None,
                  xmax=None,ymin=None,ymax=None,xlabel=None,nu=None,**kwargs):
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
    if (nu is None):
        nu = get_frequency(spectrum['units'],xfaces)

    # Initialize x labels
    if (xlabel is None):
        if (xunit == 'kev'):
            xlabel = r"$E {\rm (keV)}$"
        if (xunit == 'ev'):
            xlabel = r"$E {\rm (eV)}$"
        if (xunit == 'nu'):
            xlabel = r"$\nu {\rm (Hz)}$"
        if (xunit == 'lambda'):
            xlabel = r"$\lambda {\rm (\AA)$"
    ax.set_xlabel(xlabel)


    # Check if error requested and stored
    if ploterr:
        if spectrum['yerror'] != "true":
            print("Warning: error requested but not computed in spectrum.\n")
            ploterr = False

    nintens = spectrum['nintens']
    # Compute intensity spectrum
    if ((istokes > 0) and (spectrum['polarized'] != 'true')):
        print("Warning: polarization output requested for unpolarized spectrum. \
               Plotting unpolarized intensity.")
        istokes = 0
    intensity = spectrum['intensity'][istokes,:,:,:]
    if ploterr:
        errors = spectrum['errors'][istokes,:,:,:]

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
    if (yunit == 'nulnu'):
        ylabel = r"$\nu L_\nu {\rm (erg/s)}$"
        y = intensity*nu
        if ploterr:
            yerr = errors*nu
    if (yunit == 'lnu'):
        ylabel = r"$L_\nu {\rm (erg/s/Hz)}$"
        y = intensity
        if ploterr:
            yerr = errors
    if (yunit == 'counts'):
        ylabel = r"$N_\nu {\rm (counts/s/Hz)}$"
        h = 6.62607015e-27
        y = intensity/(h*nu)
        if ploterr:
            yerr = errors/(h*nu)
    ax.set_ylabel(ylabel)

    if (ploterr):
        ax.errorbar(x,y,yerr=yerr,fmt='.',**kwargs)
    else:
        ax.plot(x,y,'.',**kwargs)

    # Set axis scales
    ax.set_xscale(xscale)
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
    return x, nu, ax

def compute_fraction_error(intensity,errors=None):
    """
    Compute fractional polarization and error (%)
    """
    i = intensity[0,:]
    q = intensity[1,:]
    u = intensity[2,:]
    frac = np.sqrt(q**2+u**2)/i
    if errors is not None:
        ei = errors[0,:]
        eq = errors[1,:]
        eu = errors[2,:]
        err = np.sqrt(((q*eq)**2+(u*eu)**2)/(q**2+u**2) + (q**2+u**2)*(ei/i)**2)/i
        return frac*100, err*100.
    else:
        return frac, None

def compute_angle_error(intensity,errors=None):
    """
    Compute polarization angle and error if requested (degrees)
    """
    q = intensity[1,:]
    u = intensity[2,:]
    angle = 90./np.pi*np.arctan2(u,q)
    if errors is not None:
        eq = errors[1,:]
        eu = errors[2,:]
        err = 90./np.pi*np.sqrt(((u*eu)**2+(q*eq)**2))/(q**2+u**2)
        return angle, err
    else:
        return angle, None

def compute_q_error(intensity,errors=None):
    """
    Compute q=Q/I and error if requested
    """
    i = intensity[0,:]
    q = intensity[1,:]
    frac = -q/i
    if errors is not None:
        ei = errors[0,:]
        eq = errors[1,:]
        err = np.sqrt(eq**2 + (q**2)*(ei/i)**2)/i
        return frac, err
    else:
        return frac, None

def plot_polarization(spectrum,imu,ax=None,iphi='ave',xunit='kev',yunit='frac',
                      ploterr=True,xscale='log',yscale='linear',xmin=None,xmax=None,
                      ymin=None,ymax=None,nu=None,**kwargs):
    """
    Plot polarization fraction or angle. Assumes saving, etc. are performed by the
    calling function
    """

    if (ax is None):
        # Create figure, axis and assume a single plot window
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    if spectrum['polarized'] != 'true':
         raise RuntimeError('Error: Polarized fraction or angle requested for \
                             unpolarized spectrum.')
    # Assumes saving, etc. are performed by the calling function

    # Set up x axis as bin midpoints
    xfaces = spectrum['xfaces']
    x = 0.5*(xfaces[1:]+xfaces[:-1])
    if (nu is None):
        nu = get_frequency(spectrum['units'],xfaces)

    # Initialize x labels
    xlabel = ""
    if (xunit == 'kev'):
        xlabel = r"$E {\rm (keV)}$"
    if (xunit == 'ev'):
        xlabel = r"$E {\rm (eV)}$"
    if (xunit == 'nu'):
        xlabel = r"$\nu {\rm (Hz)}$"
    if (xunit == 'lambda'):
        xlabel = r"$\lambda {\rm (\AA)$"
    ax.set_xlabel(xlabel)

    # Check if error requested and stored
    if ploterr:
        if spectrum['yerror'] != "true":
            print("Warning: error requested but not computed in spectrum.\n")
            ploterr = False

    nintens = spectrum['nintens']
    intensity = spectrum['intensity']
    if ploterr:
        errors = spectrum['errors']
    else:
        errors = None

    # Selection for azimuthal angle
    if ((iphi == 'ave') or (iphi == 'sum')):
        norm = 1./float(spectrum['nphi'])
        if iphi == 'sum':
            norm *= 2.*np.pi
        intensity = np.sum(intensity,axis=1)*norm
        if ploterr:
            errors = np.sqrt(np.sum((errors)**2,axis=1))*norm
    else:
        iphi = int(iphi)
        intensity = intensity[:,iphi,:,:]
        if ploterr:
            errors = errors[:,iphi,:,:]

    # Selection for polar angle
    imu = int(imu)
    intensity = intensity[:,imu,:]
    if ploterr:
        errors = errors[:,imu,:]

    # Set y, yerr, and ylabel according to input units
    if (yunit == 'frac'):
        ylabel = r"$\rm Pol.\; Fraction$"
        y, yerr = compute_fraction_error(intensity,errors)
    if (yunit == 'angle'):
        ylabel = r"$\rm Pol.\; Angle$"
        y, yerr = compute_angle_error(intensity,errors)
    if (yunit == 'q'):
        ylabel = r"$Q$"
        y, yerr = compute_q_error(intensity,errors)
    ax.set_ylabel(ylabel)

    if (ploterr):
        ax.errorbar(x,y,yerr=yerr,fmt='.',**kwargs)
    else:
        ax.plot(x,y,'.',**kwargs)


    # Set axis scales
    ax.set_xscale(xscale)
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
    return x, nu, ax

def plot_polar_angle(spectrum,ix,ax=None,iphi='ave',xunit='mu',yunit='lnu',
                     ploterr=True,xscale='log',yscale='log',istokes=0,xmin=None,
                     xmax=None,ymin=None,ymax=None,xlabel=None,nu=None,**kwargs):
    """
    Plot polar angle. Assumes saving, etc. are performed by the calling function
    """

    if (ax is None):
        # Create figure, axis and assume a single plot window
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    # Set up x axis as bin midpoints
    xfaces = spectrum['mufaces']
    x = 0.5*(xfaces[1:]+xfaces[:-1])

    # Initialize x labels
    if (xlabel is None):
        if (xunit == 'mu'):
            xlabel = r"$\cos \theta$"
    ax.set_xlabel(xlabel)


    # Check if error requested and stored
    if ploterr:
        if spectrum['yerror'] != "true":
            print("Warning: error requested but not computed in spectrum.\n")
            ploterr = False

    nintens = spectrum['nintens']
    # Compute intensity spectrum
    if ((istokes > 0) and (spectrum['polarized'] != 'true')):
        print("Warning: polarization output requested for unpolarized spectrum. \
               Plotting unpolarized intensity.")
        istokes = 0
    intensity = spectrum['intensity'][istokes,:,:,:]
    if ploterr:
        errors = spectrum['errors'][istokes,:,:,:]

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

    # Selection for frequency
    if ix == 'sum':
        nx = spectrum['nx']
        dx = (spectrum['xfaces'][1:]-spectrum['xfaces'][:-1]).reshape(nx)
        #intensity = np.dot(dx,intensity,axis=1)
        intensity = np.average(intensity,axis=1)
        #if ploterr:
        #    errors = np.sqrt(np.dot((xmid)**2,(errors)**2))/nx
    else:
        ix = int(ix)
        intensity = intensity[:,ix]
        if ploterr:
            errors = errors[:,ix]

    # Set y, yerr, and ylabel according to input units
    if (yunit == 'lnu'):
        ylabel = r"$L_\nu {\rm (erg/s/Hz)}$"
        y = intensity
        if ploterr:
            yerr = errors
    ax.set_ylabel(ylabel)

    if (ploterr):
        ax.errorbar(x,y,yerr=yerr,fmt='.',**kwargs)
    else:
        ax.plot(x,y,'.',**kwargs)

    # Set axis scales
    ax.set_xscale(xscale)
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
    return x, nu, ax


def get_luminosity(spec):
    """
    Computes the integrated luminosity corresponding to a spectrum
    """

    nx = spec['nx']
    nmu = spec['nmu']
    nphi = spec['nphi']
    mumid = 0.5*(spec['mufaces'][1:]+spec['mufaces'][:-1])

    # Compute frequency width and mean energy (in ergs) of bins
    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10
    xaxis = spec['units']
    xfaces = spec['xfaces']
    if (xaxis == 'kev'):
        dnu = (xfaces[1:]-xfaces[:-1])*1000.*everg/h
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*1000.*everg
    if (xaxis == 'ev'):
        dnu = (xfaces[1:]-xfaces[:-1])*everg/h
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*everg
    if (xaxis == 'nu'):
        dnu = (xfaces[1:]-xfaces[:-1])
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*h
    if (xaxis == 'lambda'):
        dnu = (1./xfaces[:-1]-1./xfaces[1:])*c/1.e8
        emid = 0.5*(1./xfaces[:-1]-1./xfaces[1:])*c*h/1.e8

    # compute sum over frequency and solid angle
    lumin = 0.
    for k in range(nx):
        for j in range(nmu):
            for i in range(nphi):
                lumin += dnu[k]*mumid[j]*spec['intensity'][0,i,j,k]
    lumin *= 2.*np.pi/nmu/nphi
    return lumin

def build_bins(xmin,xmax,nx,logx):
    """
    Builds a x-axis grid for binning the photons
    """
    if (logx):
        return np.logspace(np.log10(xmin),np.log10(xmax),nx+1)
    else:
        return np.linspace(xmin,xmax,nx+1)

def get_bins(xphots,xfaces,nx,uniform=True,log=True):

    xbins = np.zeros(len(xphots),dtype=int)-1
    if uniform:
        if log:
            xlfaces = np.log10(xfaces)
            xwidth = xlfaces[nx]-xlfaces[0]
            for i,xphot in enumerate(xphots):
                xbins[i] = int((np.log10(xphot)-xlfaces[0])/xwidth*float(nx))
                if ((xbins[i] < 0) or (xbins[i] >= nx)):
                    xbins[i] = -1
        else:
            xwidth = xfaces[nx]-xfaces[0]
            for i,xphot in enumerate(xphots):
                if np.isinf(xphot):
                    xbins[i] = -1
                else:
                    xbins[i] = int((xphot-xfaces[0])/xwidth*float(nx))
                if ((xbins[i] < 0) or (xbins[i] >= nx)):
                    xbins[i] = -1
        return xbins
    else:
        # use binary search for non uniform data
        return get_bins_binary_search(xphots,xfaces,nx)

def get_bins_binary_search(xphots,xfaces,nx):
    """
    Returns x bin numbers corresponding to xphots for non-uniformly
    binned data via binary search.
    """
    # Exclude values outside of search range
    indsp = (xphots > xfaces[nx]).nonzero()
    indsm = (xphots < xfaces[0]).nonzero()
    xbins = np.searchsorted(xfaces,xphots)-1
    xbins[indsp] = -1
    xbins[indsm] = -1
    return xbins

def get_angle_bins_cartesian(photons,nmu,mufaces,nphi,phifaces):

    if ((nmu == 1) and (mufaces[0] <= 0.) and (mufaces[1] >= 1.0)):
        skipmu = True
    else:
        skipmu = False
    if ((nphi == 1) and (phifaces[0] <= 0.) and (phifaces[1] >= 2.*np.pi)):
        skipphi = True
    else:
        skipphi = False

    if (skipmu and skipphi):
        return np.zeros(photons.nphot,dtype=int),np.zeros(photons.nphot,dtype=int)

    if photons.coord == 'spherical':
        kr = photon.k1
        kth = photon.k2
        kph = photon.k3
        cth = np.cos(photon.x2)
        sth = np.sin(photon.x2)
        cph = np.cos(photon.x3)
        sph = np.sin(photon.x3)
        if (not skipphi):
            kx = kr*sth*cph + kth*cth*cph - kph*sph
            ky = kr*sth*sph + kth*cth*sph - kph*cph
        if (not skipmu):
            kz = kr*cth - kth*sth
    else:
        if (not skipphi):
            kx = photons.k1
            ky = photons.k2
        if (not skipmu):
            kz = photons.k3

    if (skipmu):
        # return 0
        mubins = np.zeros(photons.nphot,dtype=int)
    else:
        # Bin based on k . z
        mu = abs(kz)
        mubins = get_bins(mu,mufaces,nmu,log=False)
    if (skipphi):
        # return 0
        phibins = np.zeros(photons.nphot,dtype=int)
    else:
        if (skipmu):
            mu = abs(kz)
        smu = np.sqrt(1.-mu*mu)

        phi = np.arctan2(ky,kx)
        phi[(phi<0.).nonzero()] += 2.*np.pi
        phibins = get_bins(phi,phifaces,nphi,log=False)


    return mubins,phibins


def make_spectrum(phots,nx,xmin,xmax,xaxis='kev',logx=True,nmu=1,mumin=0,mumax=1.,
                  nphi=1,phimin=0,phimax=2.*np.pi,yerror=True,mask=None,xfunc=None,**kwargs):
    """
    Makes spectrum (dict) from photon object
    """

    # Store spectrum as a dictionary
    spectrum = {}

    # Store total number of photons for refernce
    spectrum['ntot'] = phots.ntot

    # Set x binning variable
    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10
    preset = False
    spectrum['xaxis'] = xaxis
    if (xaxis == 'kev'):
        xphots = phots.energy/everg/1000.
        preset = True
    if (xaxis == 'ev'):
        xphots = phots.energy/everg
        preset = True
    if (xaxis == 'nu'):
        xphots = phots.energy/h
        preset = True
    if (xaxis == 'lambda'):
        xphots = c*h/(phots.energy*1.e8)
        preset = True
    if (not preset):
        if xfunc is None:
            raise RuntimeError('Unrecognized xunit and no xfunc provided')
        else:
            xphots = xfunc(phots.energy,True,**kwargs)

    # Create bins
    xfaces = build_bins(xmin,xmax,nx,logx)
    spectrum['nx'] = nx
    spectrum['xfaces'] = xfaces

    # Get x bins
    xbins = get_bins(xphots,xfaces,nx,log=logx)

    # Get angle bins
    spectrum['nmu'] = nmu
    mufaces = build_bins(mumin,mumax,nmu,False)
    spectrum['mufaces'] = mufaces

    spectrum['nphi'] = nphi
    phifaces = build_bins(phimin,phimax,nphi,False)
    spectrum['phifaces'] = phifaces

    mubins, phibins = get_angle_bins_cartesian(phots,nmu,mufaces,nphi,phifaces)

    # Create intensity grid and loop over photons to add contribution
    nintens = 1
    if phots.polarized:
        spectrum['polarized'] = 'true'
        nintens += 2
    else:
        spectrum['polarized'] = 'false'

    spectrum['nintens'] = nintens
    intensity = np.zeros((nintens,nphi,nmu,nx))
    if yerror:
        errors = np.zeros((nintens,nphi,nmu,nx))

    if (mask is not None):
        xbins[mask] = -1

    for i in range(phots.nphot):
        if ((xbins[i] >= 0) and (mubins[i] >= 0) and (phibins[i] >= 0)):
            wght = phots.weight[i]
            #print phibins[i],mubins[i],xbins[i]
            intensity[0,phibins[i],mubins[i],xbins[i]] += wght
            if phots.polarized:
                intensity[1,phibins[i],mubins[i],xbins[i]] += wght*phots.q[i]
                intensity[2,phibins[i],mubins[i],xbins[i]] += wght*phots.u[i]
            if yerror:
                errors[0,phibins[i],mubins[i],xbins[i]] += (wght)**2
                if phots.polarized:
                    errors[1,phibins[i],mubins[i],xbins[i]] += (wght*phots.q[i])**2
                    errors[2,phibins[i],mubins[i],xbins[i]] += (wght*phots.u[i])**2

    # Compute frequency width and mean energy (in erg) of bins
    h = 6.62607015e-27
    everg = 1.6021772e-12
    c = 2.99792e10
    if (xaxis == 'kev'):
        dnu = (xfaces[1:]-xfaces[:-1])*1000.*everg/h
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*1000.*everg
    if (xaxis == 'ev'):
        dnu = (xfaces[1:]-xfaces[:-1])*everg/h
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*everg
    if (xaxis == 'nu'):
        dnu = (xfaces[1:]-xfaces[:-1])
        emid = 0.5*(xfaces[1:]+xfaces[:-1])*h
    if (xaxis == 'lambda'):
        dnu = (1./xfaces[:-1]-1./xfaces[1:])*c/1.e8
        emid = 0.5*(1./xfaces[:-1]-1./xfaces[1:])*c*h/1.e8
    if (not preset):
        efaces = xfunc(xfaces,False,**kwargs)
        emid = 0.5*(efaces[1:]+efaces[:-1])
        dnu  = (efaces[1:]-efaces[:-1])/h

    # Normalize intensities
    mumid = 0.5*(mufaces[1:]+mufaces[:-1])
    for k in range(nintens):
        for j in range(nphi):
            for i in range(nmu):
                fac = nphi*nmu*emid/(mumid[i]*dnu*2.*np.pi)
                intensity[k,j,i,:] *= fac/phots.ntot
                if yerror:
                    errors[k,j,i,:] *= fac**2/phots.ntot
    spectrum['intensity'] = intensity

    # Finish computing errors on intensities
    if yerror:
        for k in range(nintens):
            errors[k,:,:,:] = 0.675*np.sqrt((errors[k,:,:,:] - (intensity[k,:,:,:])**2)/
                                            phots.ntot)
        spectrum['errors'] = errors

    if yerror:
        spectrum['yerror'] = "true"
    else:
        spectrum['yerror'] = "false"

    return spectrum
