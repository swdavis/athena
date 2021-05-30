#! /usr/bin/env python

"""
Plot escape time distribution

"""

# standard python modules
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# athena++ modules
import athena_mc_spec as mcspec



def tanf(x,tau):
    return np.tan(x) - x/(1.-1.5*tau)

def prob_ct_tau(t,tau):
    """
    Compute (derivative) of cumulative probability
    """
    y = np.exp(-t)
    prob = 0
    n = 1
    sol = optimize.root_scalar(tanf,args=(tau),bracket=[0.51*np.pi,1.49*np.pi])
    lamn = sol.root
    dlamn = lamn
    In = 0.5*tau*(1+(1.5*tau-1)/((1.5*tau-1)**2+lamn*lamn))
    yn2 = (y)**(lamn*lamn/(np.pi*np.pi))*np.cos(lamn)*(1.5*tau**2/(1-1.5*tau)/In)*(lamn*lamn/(np.pi*np.pi))

    # Compute the sum to the limit of double precision        
    while (abs(yn2) > 1.e-17 ): 
        prob += yn2
        n = n+1;
        bracket = [lamn+(1.-0.1/n)*dlamn,lamn+(1.+0.1/n)*dlamn]
        sol = optimize.root_scalar(tanf,args=(tau),bracket=bracket)
        dlamn = sol.root-lamn
        lamn = sol.root
        In = 0.5*tau*(1+(1.5*tau-1)/((1.5*tau-1)**2+lamn*lamn))
        yn2 = (y)**(lamn*lamn/(np.pi*np.pi))*np.cos(lamn)*(1.5*tau**2/(1-1.5*tau)/In)*(lamn*lamn/(np.pi*np.pi))

    return prob
 

# Main function
def main(**kwargs):

   # Use latex labels
    plt.rc('text',usetex=True)
    plt.rc('font', **{'family' :"serif"})

    # Filenames for io
    infile = kwargs.pop('infile')
    outfile = kwargs.pop('outfile')

    radius = kwargs.pop('radius')
    #tau = kwargs.pop('tau')
    #mfp = radius/tau
    #time0 = mfp*tau**2*3./np.pi**2
    time0 = 1.e11

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Read escape time dist
    tdist = np.loadtxt(infile)
    tmid = 0.5*(tdist[:,0]+tdist[:,1])/time0
    pt = tdist[:,2]*time0

    # Compute comparison function
    #pt_comp = np.zeros(len(tmid))
    #for i,t in enumerate(tmid):
    #    pt_comp[i] = prob_ct_tau(t,tau)

    ax.set_ylabel(r"$P(t)$")
    ax.set_xlabel(r"$t$")
    # Plot times
    if (kwargs['notnorm']): 
        c = 2.99792e10
        tmid *= time0/c
        ax.set_xlabel(r"$t \; \rm(s)$")
    ax.plot(tmid,pt,'.')
    #ax.plot(tmid,pt_comp,':')

    # Set axis scales
    ax.set_xscale('linear')
    ax.set_yscale('log')


    # Write distribution to file
    if outfile is None:
        outfile = infile.replace('.tdist','.pdf')
        
    plt.savefig(outfile)
    plt.close()


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input filename')
    parser.add_argument('tau',
        type=float,
        help='optical depth of sphere')
    parser.add_argument('--radius',
        type=float,
        default = 1.e10,
        help='radius of sphere')
    parser.add_argument('--notnorm',
        action='store_true',
        help='Sets t axis to physical units')
    parser.add_argument('--outfile',
        default=None,
        help='output filename for escape time plot')

    args = parser.parse_args()
    main(**vars(args))
