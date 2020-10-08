"""
Plotting support for Monte Carlo spectra
"""

# standard python modules
import numpy as np
import matplotlib.pyplot as plt

# Athena++ modules
import athena_mc_io as mcio
from athena_mc_io import photons

# Retrun xmid, nu for desired units
def convert_xaxis(baseunit,newunit,spectrum):

    h = 6.6262e-27
    everg = 1.6021772e-12
    c = 2.99792e10

    xfaces = spectrum['xfaces']
    if (baseunit == 'kev'):
        nu = 0.5*(xfaces[1:]+xfaces[:-1])*1000.*everg/h
    if (baseunit == 'ev'):
        nu = 0.5*(xfaces[1:]+xfaces[:-1])*everg/h
    if (baseunit == 'nu'):
        nu = 0.5*(xfaces[1:]+xfaces[:-1])
    if (baseunit == 'lambda'):
        nu = 0.5*(1./xfaces[:-1]+1./xfaces[1:])*c/1.e8
    if (newunit == 'kev'):
        xmid = nu*h/(everg*1000.)
    if (newunit == 'ev'):
        xmid = nu*h/everg
    if (newunit == 'nu'):
        xmid = nu
    if (newunit == 'lambda'):
        xmid = c/nu*1.e8
    return xmid,nu

def plotspec(spectrum,imu,iphi='ave',xunit='kev',yunit='nulnu',ploterr=True,
             xscale='log',yscale='log',istokes=0,xmin=None,xmax=None,
             ymin=None,ymax=None,**kwargs):
    """
    Plot spectrum
    """

    # Assumes saving, etc. are performed by the calling function

    # Set up x axis as bin midpoints
    x, nu = convert_xaxis(spectrum['units'],xunit,spectrum)
    
    # Initialize x labels
    if (xunit == 'kev'):
        xlabel = r"$E {\rm (keV)}$"
    if (xunit == 'ev'):
        xlabel = r"$E {\rm (eV)}$"
    if (xunit == 'nu'):
        xlabel = r"$\nu {\rm (Hz)}$"
    if (xunit == 'lambda'):
        xlabel = r"$\lambda {\rm (\AA)$"
    plt.xlabel(xlabel)


    # Check if error requested and stored
    if ploterr:
        if spectrum['yerror'] != "true":
            print("Warning: error requested but not computed in spectrum.\n")
            ploterr = False

    nintens = spectrum['nintens']
    # Compute intensity spectrum
    intensity = spectrum['intensity'][istokes,:,:,:]
    if ploterr:
        errors = spectrum['errors'][istokes,:,:,:]

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
        yerr = errors
    if (yunit == 'counts'):
        ylabel = r"$N_\nu {\rm (counts/s/Hz)}$"
        h = 6.6262e-27
        y = intensity/(h*nu)
        yerr = errors/(h*nu)
    plt.ylabel(ylabel)


    if (ploterr):
        plt.errorbar(x,y,yerr=yerr,fmt='.',**kwargs)
    else:
        plt.plot(x,y,'.',**kwargs)

    # Set axis scales
    plt.xscale(xscale)
    plt.yscale(yscale)

    # (re)Set plot ranges
    left,right = plt.xlim()
    if xmin is not None:
        left=float(xmin)
    if xmax is not None:
        right=float(xmax)
    plt.xlim([left,right])
    left,right = plt.ylim()
    if ymin is not None:
        left=float(ymin)
    if ymax is not None:
        right=float(ymax)
    plt.ylim([left,right])

    return x, nu

def get_luminosity(spec):
    """
    Computes the integrated luminosity corresponding to a spectrum
    """

    nx = spec['nx']
    nmu = spec['nmu']
    nphi = spec['nphi']
    mumid = 0.5*(spec['mufaces'][1:]+spec['mufaces'][:-1])

    # Compute frequency width and mean energy (in ergs) of bins
    h = 6.6262e-27
    everg = 1.6021772e-12
    c = 2.99792e10
    xaxis = spec['xaxis']
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


# function for building bins
def build_bins(xmin,xmax,nx,logx):
    """
    Builds a x-axis grid for binning the photons
    """
    if (logx):
        return np.logspace(np.log10(xmin),np.log10(xmax),nx+1)
    else:
        return np.linspace(xmin,xmax,nx+1)

def get_x_bins(xphots,xfaces,nx):
    """
    Returns x bin numbers corresponding to xphots. 
    """
    xbins = np.zeros(len(xphots),dtype=int)-1
    for i,xphot in enumerate(xphots):
        if ((xphot < xfaces[0]) and (xphot > xfaces[nx])):
            xbins[i] = -1
            break

        # perform binary search
        low = 0
        high = nx+1
        while (low <= high):
            mid = (low+high)/2
            if (xfaces[mid] <= xphot):
                if (xfaces[mid+1] > xphot):
                    xbins[i] = mid
                    break
                else:
                    low = mid+1
            else:
                high = mid-1
        if (xbins[i] != mid):
            if (mid == high):
                xbins[i] = mid
            else:
                print("Warning: binary search failed for x: {:e}.".format(xphot))
                xbins[i] = -1

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
        return np.zeros(photons.nphot),np.zeros(photons.nphot)

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
        mubins = np.zeros(photons.nphot)
    else:
        # Bin based on k . z
        mu = abs(kz)
        mubins = get_x_bins(mu,mufaces,nmu)
    if (skipphi):
        # return 0
        phibins = np.zeros(photons.nphot)
    else:
        if (skipmu):
            mu = abs(kz)
        smu = np.sqrt(1.-mu*mu)

        phi = np.arctan2(ky,kx)
        phi[(phi<0.).nonzero()] += 2.*np.pi
        phibins = get_x_bins(phi,phifaces,nphi)


    return mubins,phibins


def makespec(phots,nx,xmin,xmax,xaxis='kev',logx=True,nmu=1,mumin=0,mumax=1.,nphi=1,
             phimin=0,phimax=2.*np.pi,yerror=True):
    """
    Makes spectrum from photon list
    """

    # Store spectrum as a dictionary
    spectrum = {}

    # Store total number of photons for refernce
    spectrum['ntot'] = phots.ntot

    # Set x binning variable
    h = 6.6262e-27
    everg = 1.6021772e-12
    c = 2.99792e10
    spectrum['xaxis'] = xaxis
    if (xaxis == 'kev'):
        xphots = phots.energy/everg/1000.
    if (xaxis == 'ev'):
        xphots = phots.energy/everg
    if (xaxis == 'nu'):
        xphots = phots.energy/h
    if (xaxis == 'lambda'):
        xphots = c*h/(phots.energy*1.e8)

    # Create bins
    xfaces = build_bins(xmin,xmax,nx,logx)
    spectrum['nx'] = nx
    spectrum['xfaces'] = xfaces

    # Get x bins
    xbins = get_x_bins(xphots,xfaces,nx+1)

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
    spectrum['nintens'] = nintens
    intensity = np.zeros((nintens,nphi,nmu,nx))
    if yerror:
        errors = np.zeros((nintens,nphi,nmu,nx))
    for i in range(phots.nphot):
        intensity[0,phibins[i],mubins[i],xbins[i]] += phots.weight[i]
        if yerror:
            errors[0,phibins[i],mubins[i],xbins[i]] += (phots.weight[i])**2

    # Compute frequency width and mean energy (in ergs) of bins
    h = 6.6262e-27
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

    # Normalize intensities
    mumid = 0.5*(mufaces[1:]+mufaces[:-1])
    for j in range(nphi):
        for i in range(nmu):
            fac = nphi*nmu*emid/(mumid[i]*dnu*2.*np.pi)
            intensity[0,j,i,:] *= fac/phots.ntot
            if yerror:
                errors[0,j,i,:] *= fac**2/phots.ntot

    spectrum['intensity'] = intensity

    # Finish computing errors
    if yerror:
        errors[0,:,:,:] = 0.675*np.sqrt((errors[0,:,:,:] - (intensity[0,:,:,:])**2)/
                                        phots.ntot)
        spectrum['errors'] = errors

    if yerror:
        spectrum['yerror'] = "true"
    else:
        spectrum['yerror'] = "false"

    return spectrum


