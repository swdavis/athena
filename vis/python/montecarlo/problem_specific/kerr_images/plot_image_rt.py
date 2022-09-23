#! /usr/bin/env python

"""
Plot a black hole image
"""

# standard python modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# athena++ modules
import athena_mc as athenamc
from athena_mc import photons

class AthenaError(RuntimeError):
  """General exception class for Athena++ read functions."""
  pass

class AthenaWarning(RuntimeWarning):
  """General warning class for Athena++ read functions."""
  pass

def blackbody(teff, nu):
  """
  Calculates blackbody spectrum
  """
  h = 6.62607015e-27; # planck constant [cgs]
  kb = 1.380649e-16 # bolztmann constant [cgs]
  c = 2.99792458e10 # speed of light [cgs]
  efact = np.exp(h*nu/(kb*teff))
  bnu = 2.*h*nu**3/(c**2*(efact - 1.)) # [cgs] = [erg s^-1 cm^-2 Hz^-1]
  return bnu

def DiskFlux(mass,mdot,radius,abh):
    """
    Returns the flux as a function of radius, mass, accretion rate but without
    relativistic and no-torque correction, which is provided by grcor
    """
    G = 6.67430e-8
    msun = 1.99e33
    kapes = 0.34
    c = 2.99792458e10

    # Minimum radius for last stable circular orbit per unit mass, X0
    abh2 = abh*abh
    signa = 1.0
    z1 = 1.0+(1.0-abh2)**(1.0/3.0)*((1.0+abh)**(1.0/3.0)+(1.0-abh)**(1.0/3.0))
    z2 = np.sqrt(3.0*abh2+z1*z1)
    r0 = 3.0+z2-signa*np.sqrt((3.0-z1)*(3.0+z1+2.0*z2))

    eta = 1.0 - (r0**2 - 2.0*r0 + abh*r0**0.5)/r0/(r0**2 - 3.0*r0 + 2.0*abh*r0**0.5)**0.5

    flux = 1.5*c**5/(G*eta*kapes*msun)/mass/radius**3*mdot
    flux[np.where(radius < r0)] = 0.

    return flux

def grcor(abh,r):
  """
  Pythonized version of dispar.f routine used with kerrtrans
  """

  if(abh < 0.):
    signa=-1.
  else:
    signa=1.0

  abh2 = abh * abh
  r1 = 1.0/r
  r12 = np.sqrt(r1)
  r2 = r1*r1
  a2r2 = abh2*r2
  a4r4 = a2r2*a2r2
  a2r3 = abh2*r2*r1
  ar32 = np.sqrt(a2r3)

  A = 1 + a2r2 + 2.*a2r3
  B = 1 + ar32
  C = 1 - 3.*r1 + 2.*ar32
  D = 1 - 2.*r1 + a2r2
  E = 1 + 4.*a2r2 - 4.*a2r3 + 3*a4r4

  # gravity correction  (see Page & Thorne,'73, eq.35)
  qcor = (1. - 4.*ar32 + 3.*a2r2)/C

  # Minimum radius for last stable circular orbit per unit mass, x0
  z1 = 1.0+(1.0-abh2)**(1.0/3.0)*((1.0+abh)**(1.0/3.0)+(1.0-abh)**(1.0/3.0))
  z2 = np.sqrt(3.0*abh2+z1*z1)
  r0 = 3.0+z2-signa*np.sqrt((3.0-z1)*(3.0+z1+2.0*z2))
  x0 = np.sqrt(r0)

  # roots of x^3 - 3x + 2a = 0
  ca3 = 1.0/3.0 * np.arccos(abh)
  x1 =  2.*np.cos(ca3-np.pi/3.)
  x2 =  2.*np.cos(ca3+np.pi/3.)
  x3 = -2.*np.cos(ca3)

  # FB = '[]' term in eq. (35) of Page&Thorne '73
  x = np.sqrt(r)
  c1 = 3*(x1-abh)*(x1-abh)/(x1*(x1-x2)*(x1-x3))
  c2 = 3*(x2-abh)*(x2-abh)/(x2*(x2-x1)*(x2-x3))
  c3 = 3*(x3-abh)*(x3-abh)/(x3*(x3-x1)*(x3-x2))
  al0 = 1.5*abh * np.log(x/x0)
  al1 = np.log((x-x1)/(x0-x1))
  al2 = np.log((x-x2)/(x0-x2))
  al3 = np.log((x-x3)/(x0-x3))
  fb = (x-x0 - al0 - c1*al1 - c2*al2 - c3*al3)
  Q = fb*(1.0+ar32)*r12/np.sqrt(1.0-3.0*r1+2.0*ar32)

  # temperature correction
  tcor = (Q/B/np.sqrt(C))**0.25

  inds0 = np.where(r < r0)

  #nx = r.shape[0]
  #ny = r.shape[1]

  tcor[inds0] = 0.
  qcor[inds0] = 0.

  return qcor,tcor

def intensityblackholedisk(phots,nx,ny,mbh,abh,mdot,nuobs):
    """
    Generate pixels corresponding to intensity in Novikov-Thorne
    black hole accretion disk
    """

    # set up radius and intensity corrections as nx * ny arrays
    # switch so that x is fastest running index
    nufac = np.transpose(phots.user[:,2].reshape((nx,ny)))
    radius = np.transpose(phots.user[:,3].reshape((nx,ny)))

    tcor = grcor(abh, radius)[1]
    sigma = 5.67e-5 # stefan-boltzmann constant [cgs]
    teff = (DiskFlux(mbh, mdot, radius, abh) / sigma)**(1./4.) * tcor
    indsbad = np.where(teff == 0.)
    teff[indsbad] = 1.e-50
    nufac[indsbad] = 1.e-50
    nuem = nuobs * nufac
    intens_em = blackbody(teff, nuem)
    intens_obs = intens_em / (nufac)**3
    intens_obs[indsbad] = 1.e-50

    return intens_obs

def make_image(phots,nx,ny,abh,thcam):
    """
    Makes image (dict) from photon object
    """

    # Store image as a dictionary
    image = {}

    # Store integration time
    image['dt'] = 0.

    # Ray traced images have only one inclination value
    image['ninc'] = ninc = 1
    cthc = np.cos(thcam/180.*np.pi)
    image['ifaces'] = np.array([cthc,cthc])

    # Currently assume single frequency/energy for simplicity
    image['nen'] = nen = 1
    image['efaces'] = np.array([0.,1.])

    # Store total number of photons for refernce
    image['ntot'] = phots.ntot

    if (nx*ny != phots.ntot):
        raise AthenaError("Error: nx*ny: {:d} not equal to ntot: {:e}."
                          .format(nx*ny,phots.ntot))

    if (phots.user[ny,1] != phots.user[0,1]):
        raise AthenaError("Error: y(ny): {:e} not equal to y(0): {:e}."
                          .format(phots.user[ny,1],phots.user[0,1]))

    # Set up vertical image coordinates
    image['ny'] = ny
    y = phots.user[0:ny,1]
    image['yfaces'] = np.zeros(ny+1)
    image['yfaces'][0] = 0.5*(3.*y[0] - y[1])
    image['yfaces'][ny] = 0.5*(3.*y[-1] - y[-2])
    image['yfaces'][1:ny] = 0.5*(y[1:]+y[:-1])
    #image['y'] = y

    # Set up horizontal image coordinates
    image['nx'] = nx
    x = phots.user[0:nx,1]
    image['xfaces'] = np.zeros(nx+1)
    image['xfaces'][0] = 0.5*(3.*x[0] - x[1])
    image['xfaces'][nx] = 0.5*(3.*x[-1] - x[-2])
    image['xfaces'][1:nx] = 0.5*(x[1:]+x[:-1])
    image['unit'] = 'rg'
    #image['x'] = x

    nintens = 1
    if (phots.polarized):
        image['polarized'] = 'true'
        nintens += 2
    else:
        image['polarized'] = 'false'
    image['nintens'] = nintens

    print("spin:",abh)
    mbh = 1.e9
    mdot = 0.1
    nuobs = 5.e15
    intensity = np.zeros([nintens,ninc,nen,ny,nx])
    intensity[0,0,0,:,:] = intensityblackholedisk(phots,nx,ny,mbh,abh,mdot,nuobs)
    if (phots.polarized):
        image['polarized'] = True
        intensity[1,0,0,:,:] = np.transpose(phots.q.reshape((nx,ny)))
        intensity[2,0,0,:,:] = np.transpose(phots.u.reshape((nx,ny)))
        for i in range(ny):
          inds = np.where(intensity[0,0,0,i,:] <= 1.e-50)
          intensity[1,0,0,i,inds] = 0.
          intensity[2,0,0,i,inds] = 0.

    image['intensity'] = intensity
    return image

def sort_list(phlist):
    """
    Sort list by value of iuser (last user variable)
    """

    # need to sort phlist['list'] which is a numpy array
    phlist['list'] = phlist['list'][np.argsort(phlist['list'][:, -1]),:]

# Main function
def main(**kwargs):

    # Filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')
    if outfile is None:
        outfile = infile.replace(".list",".img")

    sort = kwargs.pop('sort')
    # Read photon list
    phlist = athenamc.read_list(infile)
    if (sort):
      sort_list(phlist)
    phots = photons(phlist)

    # create image
    thcam = kwargs.pop('thcam')
    image = make_image(phots,kwargs['nx'],kwargs['ny'],kwargs['spin'],
                       thcam)

    # write image to file
    athenamc.write_image(outfile,image)

    # plot image and write to file
    athenamc.plot_image(image,0,**kwargs)
    outfile = outfile.replace(".img",".png")
    plt.savefig(outfile, bbox_inches='tight')


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
                        help='base input filename')
    parser.add_argument('nx',
                        type=int,
                        help='number of horizontal pixels')
    parser.add_argument('ny',
                        type=int,
                        help='number of vertical pixels')
    parser.add_argument('spin',
                        type=float,
                        help='black hole spin parameter')
    parser.add_argument('thcam',
                        type=float,
                        help='camera inclination in degrees')
    parser.add_argument('--outfile',
                        default=None,
                        help='output filename')
    parser.add_argument('-c', '--colormap',
                        default='hot',
                        help='name of Matplotlib colormap to use instead of default; \
                              hot is default, twilight is good for cyclic variables \
                              like polarization angle')
    parser.add_argument('-p', '--pvec',
                        action='store_true',
                        default=False,
                        help='flag indicating that polarization should be plotted')
    parser.add_argument('-a', '--average',
                        action='store_true',
                        default=False,
                        help='flag indicating that polarization should be averaged \
                              over steps')
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
    parser.add_argument('--sort',
                        action='store_true',
                        help='sort list by last user variable (should be iphot)')
    parser.add_argument('--step',
                        type=int,
                        default=4,
                        help='flag indicating that polarization should be normalized \
                              to maximum')
    parser.add_argument('--logc',
                        action='store_true',
                        help='flag indicating data should be colormapped logarithmically')

    args = parser.parse_args()
    main(**vars(args))
