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
import athena_mc_list as mclist
from athena_mc_list import photons

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
    c = 2.998e10 # speed of light [cgs]
    efact = np.exp(h*nu/(kb*teff))
    bnu = 2.*h*nu**3/(c**2*(efact - 1.)) # [cgs] = [erg s^-1 cm^-2 Hz^-1]
    return bnu

def DiskFlux(mass,mdot,radius,abh):
    """
    Returns the flux as a function of radius, mass, accretion rate but without relativistic
    and no-torque correction, which is provided by grcor
    """
    G=6.67e-8
    msun=1.99e33
    kapes=0.33
    c=3.0e10

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

#       IDLized version of dispar.f routine used with kerrtrans
#       SWD: Needs to be cleaned up                  
#                                                                                
#       Procedure for computing general-relativistic correction                  
#       factors to gravitational factor (QGRAV) and effective                    
#       temperature (TEFF)                                                       
#       Also calculates all frour quantities in the Riffer-Herlod (RH)           
#       notation - arh, BRH, CRH, DRH                                            
#                                                                                
#       Input:                                                                   
#             abh   - angular momentum (0.98 maximum)                            
#            r   - R/R_g = r/(GM/c^2)                                           
#       Outout:
#             qcor - g-correction  = C/B   in RH notation                        
#             tcor - T-correction  = (D/B)^(1/4)   in RH notation                
#             arh  - A  in RH notation                                           
#             brh  - B  in RH notation                                           
#             crh  - C  in RH notation                                           
#             drh  - D  in RH notation                                           
#                                                                                
#  ----------------                                                              
#  Imput parameters                                                              
#  ----------------                                                              
#                                                                                
#       abh     - specific angular momentum/mass                                 
#                 of the Kerr black hole                                         
#       r       - distance/mass of the Kerr black hole                           
#                                                                                
#  ---------------------------------                                             
#  Set correcion factors A through G  (see Novikov & Thorne,'73, eq.5.4.1a-g)    
#  ---------------------------------                                             
#

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

    A = 1 +   a2r2  + 2.*a2r3
    B = 1 +   ar32
    C = 1 - 3.*r1   + 2.*ar32
    D = 1 - 2.*r1   +    a2r2
    E = 1 + 4.*a2r2 - 4.*a2r3 + 3*a4r4

    #  correction - after Riffert and Harold                                         

    qcor = (1. - 4.*ar32 + 3.*a2r2)/C
#  -----------------------                                                       
#  Set correction factor Q  (see Page & Thorne,'73, eq.35)                       
#  -----------------------                                                       

# Minimum radius for last stable circular orbit per unit mass, X0                
    z1 = 1.0+(1.0-abh2)**(1.0/3.0)*((1.0+abh)**(1.0/3.0)+(1.0-abh)**(1.0/3.0))
    z2 = np.sqrt(3.0*abh2+z1*z1)
    r0 = 3.0+z2-signa*np.sqrt((3.0-z1)*(3.0+z1+2.0*z2))
    x0 = np.sqrt(r0)

    #       Roots of x^3 - 3x + 2a = 0                                               

    ca3 = 1.0/3.0 * np.arccos(abh)
    x1 =  2.*np.cos(ca3-np.pi/3.)
    x2 =  2.*np.cos(ca3+np.pi/3.)
    x3 = -2.*np.cos(ca3)

    #       FB = '[]' term in eq. (35) of Page&Thorne '73                            
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

        
    #  ------------------------------                                                           
    #  Set correction factor for TEFF  (see Novikov & Thorne,'73, eq.5.5.14b)                
    #  ------------------------------                                                           
    
    tcor = (Q/B/np.sqrt(C))**0.25

    inds0 = np.where(r < r0)

    nx = r.shape[0]
    ny = r.shape[1]

    tcor[inds0] = 0.
    qcor[inds0] = 0.

    return qcor,tcor

def intensityblackholedisk(phots,nx,ny,mbh,abh,mdot,nuobs):
    """
    Generate pixels corresponding to intensity in Novikov-Thorne black hole accreion disk
    """
    
    # set up radius and intensity corrections as nx * ny arrays
 
    nufac = phots.user[:,2].reshape((nx,ny))
    radius = phots.user[:,3].reshape((nx,ny))

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

def make_image(phots,nx,ny,polarization=False):

    """
    Makes image (dict) from photon object
    """

    # Store image as a dictionary
    image = {}

    # Store total number of photons for refernce
    image['ntot'] = phots.ntot

    if (nx*ny != phots.ntot):
        raise AthenaError("Error: nx*ny: {:d} not equal to ntot: {:e}.".format(nx*ny,phots.ntot))

    if (phots.user[ny,1] != phots.user[0,1]):
        raise AthenaError("Error: y(ny): {:e} not equal to y(0): {:e}.".format(phots.user[ny,1],phots.user[0,1]))

    # Set up vertical image coordinates
    image['ny'] = ny
    y = phots.user[0:ny,1]
    image['yfaces'] = np.zeros(ny+1)
    image['yfaces'][0] = 2.*y[0] - y[1]
    image['yfaces'][ny] = 2.*y[ny-1] - y[ny-2]
    image['yfaces'][1:ny] = 0.5*(y[1:]+y[:-1])
    image['y'] = y
    
    # Set up horizontal image coordinates
    image['nx'] = nx
    x = phots.user[0:nx,1]
    image['xfaces'] = np.zeros(nx+1)
    image['xfaces'][0] = 2.*x[0] - x[1]
    image['xfaces'][nx] = 2.*x[nx-1] - x[nx-2]
    image['xfaces'][1:nx] = 0.5*(x[1:]+x[:-1])
    image['x'] = x
  
    mbh = 1.e9
    abh = 0.99
    mdot = 0.1
    nuobs = 5.e15
    image['intens'] = intensityblackholedisk(phots,nx,ny,mbh,abh,mdot,nuobs) 
    if (polarization):
        image['q'] = phots.q.reshape((nx,ny))
        image['u'] = phots.u.reshape((nx,ny)) 

    return image

def subsample_polarization(q,u,x,y,step,average):

    nx = len(x)
    ny = len(y)

    if (not average):
        x = x[step/2:nx:step]
        y = y[step/2:ny:step]
        q = q[step/2:nx:step,step/2:ny:step]
        u = u[step/2:nx:step,step/2:ny:step]
        return q,u,x,y
    else:
        # too lazy to work out pythony way of doing this
        xp = np.zeros(nx/step)
        yp = np.zeros(ny/step)
        qp = np.zeros((nx/step,ny/step))
        up = np.zeros((nx/step,ny/step))
        for i in range(nx/step):
            xp[i] = np.average(x[i*step:(i+1)*step])
        for i in range(ny/step):
            yp[i] = np.average(y[i*step:(i+1)*step])
        for i in range(nx/step):
            for j in range(ny/step):
                qp[i,j] = np.average(q[i*step:(i+1)*step,j*step:(j+1)*step])
                up[i,j] = np.average(u[i*step:(i+1)*step,j*step:(j+1)*step])

        return qp,up,xp,yp



def plot_image(image,polarization=False,average=False,step=4,**kwargs):
    """
    Plot an image
    """
    
    vmin = kwargs['vmin']
    vmax = kwargs['vmax']
    cmap = plt.get_cmap(kwargs['colormap'])
    plt.figure()
    if (kwargs['logc']):
        norm = colors.LogNorm()
    else:
        norm = colors.Normalize()
    if kwargs['vnorm']:
        vals = image['intens'] / np.max(image['intens'])
        vmax = 1.
        if vmin is None:
            vmin = 1.e-5
    else:
        vals = image['intens']
        if vmin is None:
            vmin = 1.e-5*np.max(vals)

    x_2d, y_2d = np.meshgrid(image['xfaces'],image['yfaces'],indexing='ij')

    im = plt.pcolormesh(x_2d, y_2d, vals, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
    plt.xlim(image['xfaces'][0],image['xfaces'][image['nx']-1])
    plt.ylim(image['yfaces'][0],image['yfaces'][image['ny']-1])
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.colorbar(im)
    plt.gca().set_aspect('equal')
    if (polarization):
        x = image['x']
        y = image['y']
        q = image['q']
        u = image['u']
        q, u, x, y = subsample_polarization(q,u,x,y,step,average)
    
        x_pol, y_pol = np.meshgrid(x,y,indexing='ij')
        pol_angle = 0.5 * np.arctan2(u,q)
        pol_frac = np.sqrt(q*q+u*u)

        vx = pol_frac*np.cos(pol_angle)
        vy = pol_frac*np.sin(pol_angle)
        plt.quiver(x_pol, y_pol, vx, vy, color='k',headwidth=0, headlength=0, headaxislength=0, scale = None,pivot='middle')

# Main function
def main(**kwargs):
 
    # Filenames for io
    infile = kwargs.pop('infile')
    if kwargs['outfile'] is None:
        kwargs['outfile'] = infile+".png"

    # Read photon list
    phlist = mclist.read_list(infile+'.list')
    phots = photons(phlist)
    print phots.q[np.where(phots.q < 0.9)]
    image = make_image(phots,kwargs['nx'],kwargs['ny'],kwargs['polarization'])

    plot_image(image,**kwargs)
    plt.savefig(kwargs['outfile'], bbox_inches='tight')


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
    parser.add_argument('--outfile',
        default=None,
        help='output filename')
    parser.add_argument('-c', '--colormap',
        default='hot',
        help='name of Matplotlib colormap to use instead of default; highly recommended; \
               try "RdBu_r" or "gist_heat" if looking for suggestions')
    parser.add_argument('-p', '--polarization',
        action='store_true',
        default=False,
        help='flag indicating that polarization should be plotted')
    parser.add_argument('-a', '--average',
        action='store_true',
        default=False,
        help='flag indicating that polarization should be averaged over steps')
    parser.add_argument('--vmin',
        type=float,
        default=None,
        help='data value to correspond to colormap minimum; use --vmin=<val> if <val> has \
          negative sign')
    parser.add_argument('--vmax',
        type=float,
        default=None,
        help='data value to correspond to colormap maximum; use --vmax=<val> if <val> has \
            negative sign')
    parser.add_argument('--vnorm',
        action='store_true',
        help='flag indicating that polarization should be normalized to maximum')
    parser.add_argument('--step',
        type=int,
        default=4,
        help='flag indicating that polarization should be normalized to maximum')
    parser.add_argument('--logc',
        action='store_true',
        help='flag indicating data should be colormapped logarithmically')

    args = parser.parse_args()
    main(**vars(args))

