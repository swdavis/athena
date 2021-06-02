#! /usr/bin/env python

"""
Read in photon trajectories and plot, possibly with geokerr output for comparison
"""

# standard python modules
import argparse
import numpy as np
import matplotlib.pyplot as plt

# athena++ modules
import athena_mc_traj as mctraj

def plot_geokerr(axs,infile='mccomp.out'):
    """
    Read in geokerr output and plot
    """
    
    file = open(infile,"r")
    line = file.readline().split()
    ngeo = int(line[0])
    a = float(line[1])
    for i in range(ngeo):
        line = file.readline().split()
        n = int(line[2])+1
        ri = 1/float(line[3])
        x = np.empty(n)
        y = np.empty(n)
        z = np.empty(n)
        for j in xrange(n-1):
            line = file.readline().split()   
            if (float(line[0]) > 0.0):
                r = 1./float(line[0])
                phi = float(line[3])
                cth = float(line[1])
                sth = np.sqrt(1.-cth*cth)
                x[j] = r * sth*np.cos(phi)
                y[j] = r * sth*np.sin(phi)
                z[j] = r * cth
        x=x[0:n-1]
        y=y[0:n-1]
        z=z[0:n-1]
        axs[0].plot(x,y,'r.', markersize=1.0)
        axs[1].plot(x,z,'r.', markersize=1.0)
        axs[2].plot(y,z,'r.', markersize=1.0)
            

# Main function
def main(**kwargs):


    # Use latex labels
    plt.rc('text',usetex=True)
    plt.rc('font', **{'family' :"serif"})

    # create figure, axis
    fig, axs = plt.subplots(1, 3, figsize=(14,4))

    # Filename for io
    trajfile = kwargs.pop('trajfile')

    # parameters
    a = float(kwargs.pop('spin')) 
    blflag = kwargs.pop('bl')

    # if specified, plot geokerr output for comparison
    geokerrflag = kwargs.pop('geok')
    if (geokerrflag):
        plot_geokerr(axs)

    # Read trajectory list
    traj = mctraj.read_list(trajfile)
    nmark = 500

    # Plot projections
    for i in range(traj['length']):
        ri = traj['list'][i,0,0]
        n = traj['nsteps'][i]
        x = np.zeros(n)
        y = np.zeros(n)
        z = np.zeros(n)
        for j in range(n):
            r = traj['list'][i,j,0]
            th = traj['list'][i,j,1]
            phi = traj['list'][i,j,2]
            if (not blflag):
                a2 = a*a
                phi += -a*0.5/np.sqrt(1.-a2)*(np.log((r-1.-np.sqrt(1.-a2))
                  /(r-1.+np.sqrt(1.-a2)))-np.log((ri-1.-np.sqrt(1.-a2))
                  /(ri-1.+np.sqrt(1.-a2))))
            cth = np.cos(th)
            sth = np.sin(th)
            cph = np.cos(phi)
            sph = np.sin(phi)
            x[j] = r*sth*cph
            y[j] = r*sth*sph
            z[j] = r*cth

        axs[0].plot(x,y,'k.', markersize=2., markevery=nmark)
        axs[1].plot(x,z,'k.', markersize=2., markevery=nmark)
        axs[2].plot(y,z,'k.', markersize=2., markevery=nmark)

    axs[0].set_xlabel(r"$r \, \sin \theta \, \cos \phi$")
    axs[0].set_ylabel(r"$r \, \sin \theta \, \sin \phi$")
    axs[1].set_xlabel(r"$r \, \sin \theta \, \cos \phi$")
    axs[1].set_ylabel(r"$r \, \cos \theta$")
    axs[2].set_xlabel(r"$r \, \sin \theta \, \sin \phi$")
    axs[2].set_ylabel(r"$r \, \cos \theta$")

    # plot outer horizon as sphere
    nh = 500
    rh = 1 + np.sqrt(1.-a*a)
    xh = rh * np.cos(np.linspace(0,2.*np.pi,nh))
    yh = rh * np.sin(np.linspace(0,2.*np.pi,nh))
    for i in range(3):
        axs[i].plot(xh,yh,'k')

    l0 = kwargs.pop('l0')
    if l0 is not None:
        for i in range(3):
            axs[i].set_xlim([-l0,l0])
            axs[i].set_ylim([-l0,l0])

    # write to output file
    plt.savefig(kwargs['outfile'])
    plt.close()

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('trajfile',
        help='input trajectory filename')
    parser.add_argument('spin',
        help='black hole spin')
    parser.add_argument('--bl',
        action='store_true',
        help='assume coordinate boyer lindquist') 
    parser.add_argument('--geok',
        action='store_true',
        help='plot geokerr output') 
    parser.add_argument('--l0',
        type=float,
        default=None,                
        help='window size')
    parser.add_argument('--outfile',
        default='geodesics.pdf',                
        help='outfile pdf')
    args = parser.parse_args()
    main(**vars(args))
