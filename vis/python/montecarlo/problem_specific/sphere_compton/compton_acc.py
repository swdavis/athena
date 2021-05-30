import numpy as np
import matplotlib.pyplot as plt
import random as ran
from bisect import bisect_right
from scipy.special import kv
from scipy.special import gamma
from scipy.special import hyp1f1
from scipy.special import hyperu
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.optimize as opt
import struct
from mpmath import whitw
from mpmath import whitm


def tanf(x,tau):
    """
    Used by prob_ct_tau to evalue time/path distribution
    """
    return np.tan(x) - x/(1.-1.5*tau)

def prob_ct_tau(t,tau):
    """
    Sets up the time/path distributions for a point source at origin in the limit for any tau
    in diffusion limit
    """

    y = np.exp(-t)
    prob = 0.
  
    n = 1
    #sol = sci.optimize.root_scalar(tanf,args=(tau),x0=np.pi*n,x1=(np.pi+0.1)*n)
    
    sol = opt.root_scalar(tanf,args=(tau),bracket=[0.51*np.pi,1.49*np.pi])
    lamn = sol.root
    dlamn = lamn
    In = 0.5*tau*(1+(1.5*tau-1)/((1.5*tau-1)**2+lamn*lamn))
    yn2 = (y)**(lamn*lamn/(np.pi*np.pi))*np.cos(lamn)*(1.5*tau**2/(1-1.5*tau)/In)
    #print n,lamn,lamn/np.pi,yn2
    # Compute the sum to the limit of double precision
    while (abs(yn2) > 1.e-17 ): 
        prob = prob + yn2
        n = n+1;
        #sol=opt.root_scalar(tanf,args=(tau),x0=np.pi*n,x1=(np.pi+0.1)*n)
        bracket=[lamn+(1.-0.1/n)*dlamn,lamn+(1.+0.1/n)*dlamn]
        #print 'l: ',lamn/np.pi,dlamn/np.pi
        sol=opt.root_scalar(tanf,args=(tau),bracket=bracket)
        dlamn = sol.root-lamn
        lamn = sol.root
        In = 0.5*tau*(1+(1.5*tau-1)/((1.5*tau-1)**2+lamn*lamn))
        yn2 = (y)**(lamn*lamn/(np.pi*np.pi))*np.cos(lamn)*(1.5*tau**2/(1-1.5*tau)/In)
        #print n,lamn,lamn/np.pi,yn2,prob

    return prob

def gen_table_ct(ntau,npr,nt,taurange=[10.,100.],trange=[0.05,50],plotonly=False,plow=1.e-4,phigh=1.e-5):
    

    step = (np.log10(taurange[1])-np.log10(taurange[0]))/(ntau-1.)
    tau = 10**(np.arange(ntau)*step+np.log10(taurange[0]))
    prob0 = np.arange(npr)/(npr-1.)
    npr0 = npr/3
    npr1 = npr*2/3
    prange=[prob0[npr0],prob0[npr1]]
    print npr0,npr1,prange
    step = (prange[1]-prange[0])/((npr1-npr0)-1.)
    probm = np.arange(npr1-npr0)*step+prange[0]
    step = (np.log10(prange[0])-np.log10(plow))/(npr0-1.)
    probl = 10**(np.arange(npr0)*step+np.log10(plow))
    step = (np.log10(1.-prange[1])-np.log10(phigh))/((npr-npr1)-1.)
    probh = np.flip(1.-10**(np.arange(npr-npr1)*step+np.log10(phigh)))
    prob = np.zeros(1)
    prob = np.append(prob,probl[0:npr0-1])
    prob = np.append(prob,probm[0:(npr1-npr0-1)])
    prob = np.append(prob,probh)
    prob = np.append(prob,[1.])
    #print prob
    #print len(prob)

    # set up array of times
    step = (np.log10(trange[1])-np.log10(trange[0]))/(nt-1.)
    t = 10**(np.arange(nt)*step+np.log10(trange[0]))
    it = np.arange(nt)
    pt = np.zeros(nt)
    dt =np.zeros(nt)
    dt[0] = t[1]-t[0]
    dt[nt-1] = t[nt-1]-t[nt-2]
    for i in range(1,nt-1):
        dt[i] = 0.5*(t[i+1]-t[i-1])
    out = np.zeros([ntau,npr])

    # main loop
    for i in range(ntau):
        print i,tau[i]
        pt[0] = 0.
        for j in range(1,nt):
            pt[j] = 1.0-prob_ct_tau(t[j],tau[i])
        pt[nt-1] = 1.
        #out[i,0] = t[0]
        ip = np.interp(0.001,pt,it)
        ipm = int(ip)
        a = (0.001-pt[ipm])/(pt[ipm+1]-pt[ipm])
        a1 = 1.-a
        out[i,0] = (a*t[ipm+1]+a1*t[ipm])
        for j in range(1,npr-1):
            ip = np.interp(prob[j],pt,it)
            ipm = int(ip)
            a = (prob[j]-pt[ipm])/(pt[ipm+1]-pt[ipm])
            a1 = 1.-a
            out[i,j] = (a*t[ipm+1]+a1*t[ipm])
            #print i,j,tau[i],prob[j],out[i,j]
        ip = np.interp(0.999,pt,it)
        ipm = int(ip)
        a = (0.999-pt[ipm])/(pt[ipm+1]-pt[ipm])
        a1 = 1.-a
        out[i,npr-1] = (a*t[ipm+1]+a1*t[ipm])
        #out[i,npr-1] = out[i,npr-2]
        #print i,npr-1,tau[i],prob[npr-1],out[i,npr-1]
    if (plotonly):
        plt.plot(prob,out[0,:],'.')
        plt.xlim([0,0.2])
        plt.show()
    else:
        # Write output files
        write_binary("time_table_tau.out",np.log(tau))
        write_binary("time_table_p.out",prob)
        write_binary("time_table_t.out",out)
        np.save("time_table_tau",np.log(tau))
        np.save("time_table_p",prob)
        np.save("time_table_t",out)

def prob_ct_delta():
    """
    Sets up the distributions for a point source at origin in the limit of large tau
    or under the assumption of a perfectly absorbing boundary condition
    """
    nmax = 1000    
    mrwdev = np.zeros(nmax)
    mrwprob = np.zeros(nmax)
    for i in range(nmax-1):
        mrwdev[i] = i/(nmax-1.0)
        n = 1
        sign = 1.0
        yn2 = (mrwdev[i])**(n*n)
    #plt.plot(s,y1,'.')
        # Compute the sum to the limit of double precision
        while (yn2 > 1.e-17 ): 
            mrwprob[i] = mrwprob[i] + sign * yn2
            sign = -sign
            n = n+1;
            yn2 = mrwdev[i]**(n*n)
    
        mrwprob[i] = mrwprob[i] * 2.
        if (mrwprob[i] > 1.0): mrwprob[i] = 1.0;
  
    mrwprob[nmax-1] = 1.0;

    return mrwdev,mrwprob

def get_ct_delta(x,mrwdev,mrwprob):
    """
    Draws random path length for delta function source at origin when provided
    with random number uniformly distributed between 0 and 1 (x)
    """
    bin = bisect_right(mrwprob,x)
    if (bin > len(mrwprob)-1): bin = len(mrwprob)-1
    slope = (x - mrwprob[bin-1]) / (mrwprob[bin] - mrwprob[bin-1]);
    return mrwdev[bin-1]+(mrwdev[bin]-mrwdev[bin-1])*slope

def plot_ct_dist(ndraw=100000,nbin=50,trange=[0,2.],log=False):
    """
    Plot source path length/time distribution for testing purposes
    """
    mdev,mprob = prob_ct_delta()
    path = np.zeros(ndraw)
    for i in range(ndraw):
        x = ran.random()
        y = get_ct_delta(x,mdev,mprob)
        #print x,y
        path[i] = - np.log(y)*3./np.pi**2
        if (log):
            path[i] = np.log10(path[i])
        #print path[i]
    
    hist,binp = np.histogram(path,range=trange,bins=nbin)
    norm = 1./np.amax(hist)
    plt.plot(binp[0:nbin],norm*hist,'.')
    plt.yscale('log')
    plt.show()

def get_x_dist(x10,tau=20.,temp=3.e6,ndraw=1000,urange=[0.01,100.],nmax=100,deltasrc=True):
    """
    Function for using precomputed tabulations of the greens function of the Kompaneets
    equation to return the final energy of a photon emerging from a sphere
    of a given optical depth tau and initial photon energy x10. Developed as a test
    of the method used in the athena++ acceleration implemntation
    """

    step = (np.log10(urange[1])-np.log10(urange[0]))/(nmax-1.)
    x0 = 10**(np.arange(nmax)*step+np.log10(urange[0]))
    iout = np.zeros(nmax)

    dx0 =np.zeros(nmax)
    dx0[0] = x0[1]-x0[0]
    dx0[nmax-1] = x0[nmax-1]-x0[nmax-2]
    for i in range(1,nmax-1):
        dx0[i] = 0.5*(x0[i+1]-x0[i-1])

    c = 2.9979246e+10
    kb = 1.3806580e-16
    me = 9.1093897e-28
    if (deltasrc):
        mdev,mprob = prob_ct_delta()
    t = np.load("comp_table_t.npy")
    #tc = read_binary("compton_table_t.out",len(t))
    x1 = np.load("comp_table_x1.npy")    
    prob = np.load("comp_table_p.npy")
    nt = len(t)
    tmax = t[nt-1]
    tmin = t[0]
    print nt,tmin,tmax
    #nx = len(x)
    nx1 = len(x1)
    npr = len(prob)
    print 't range: ',np.exp(min(t)),np.exp(max(t)),nt
    print 'x1 range: ',np.exp(min(x1)),np.exp(max(x1)),nx1
    print 'p range: ',min(prob),max(prob),npr
    #print np.exp(t)
    itt = np.arange(nt)
    itp = np.arange(npr)
    itx1 = np.arange(nx1)

    table = np.load("comp_table_x.npy")

    # get interpolant for x1
    ix1 = np.interp(np.log(x10),x1,itx1)
    ix1m = int(ix1)
    b = (np.log(x10)-x1[ix1m])/(x1[ix1m+1]-x1[ix1m])
    b1 = 1.-b

    for i in range(ndraw):
        dev = ran.random()
        if (deltasrc):
            y = get_ct_delta(dev,mdev,mprob)
            path = - np.log(y)*3./np.pi**2*tau**2*kb*temp/(me*c**2)
        else:
            path = - np.log(dev)*3./np.pi**2*(tau+2./3.)**2*kb*temp/(me*c**2)
        #print i,path,dev,y,tau
        if (path >= tmax):
            print i,'path > tmax',path,tmax
            path= 0.999*tmax
        if (path < tmin):
            print i,'path < tmin',path,tmin
            path = 1.00001*tmin
        #path = 0.75
        t0 = np.log(path)
        it = np.interp(t0,t,itt)
        itm = int(it)
        a = (t0-t[itm])/(t[itm+1]-t[itm])
        a1 = 1.-a
        dev = ran.random()
        
        #ip = int(dev*(npr-1))
        ip = int(np.interp(dev,prob,itp))
        d = (dev-prob[ip])/(prob[ip+1]-prob[ip])
        d1 = 1.-d
        
        xdev = (b*(d*(a*table[itm+1,ix1m+1,ip+1]+a1*table[itm,ix1m+1,ip+1])
                   +d1*(a*table[itm+1,ix1m+1,ip]  +a1*table[itm,ix1m+1,ip]))
                +b1*(d*(a*table[itm+1,ix1m,ip+1]+a1*table[itm,ix1m,ip+1])
                   +d1*(a*table[itm+1,ix1m,ip]  +a1*table[itm,ix1m,ip])))
        #print itm,a,a1,t[itm],t[itm+1]
        #print ix1m,b,b1,x1[ix1m],x1[ix1m+1]
        #print ip,d,d1,prob[ip],prob[ip+1]
        #print xdev

        if (xdev > x0[nmax-1]):
            print i,'xdev > xmax',xdev,x0[nmax-1]
            xdev = 1.
        if (xdev < x0[0]):
            print i,'xdev < xmin',xdev,x0[0]
            xdev = 0.
        xbin = bisect_right(x0,xdev)
        #print dev,xdev,xbin,ip
        iout[xbin] = iout[xbin] + 1.

    out = np.zeros([nmax,2])
    for j in range(nmax):
        out[j,0]=x0[j]
        out[j,1]=x0[j]*iout[j]/dx0[j]/ndraw

    np.savetxt("x_green.out",out)
    plt.plot(x0,iout)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("intens_acc.pdf")
    #plt.show()
    plt.close()


def comp_x_t(x10,path,temp=3.e6,ndraw=1000,urange=[0.01,100.],nmax=100,fast=False,brisk=True,ylim=None,yscale='log'):
    """
    Function for using precomputed tabulations of the greens function of the Kompaneets
    equation to return the final energy of a photon with path. Developed as a test
    of the method used in the athena++ acceleration implemntation
    """

    plt.rc('text',usetex=True)
    plt.rc('font', **{'family' :"serif"})

    if (brisk):
        init_brisk()
    step = (np.log10(urange[1])-np.log10(urange[0]))/(nmax-1.)
    x0 = 10**(np.arange(nmax)*step+np.log10(urange[0]))
    iout = np.zeros(nmax)
    dx0 =np.zeros(nmax-1)
    xm = np.zeros(nmax-1)
    #dx0[0] = x0[1]-x0[0]
    #dx0[nmax-1] = x0[nmax-1]-x0[nmax-2]
    #for i in range(1,nmax-1):
    #    dx0[i] = 0.5*(x0[i+1]-x0[i-1])
    for i in range(nmax-1):
        dx0[i] = (x0[i+1]-x0[i])
        xm[i] = 0.5*(x0[i+1]+x0[i])

    c = 2.9979246e+10
    kb = 1.3806580e-16
    me = 9.1093897e-28

    t = np.load("comp_table_t.npy")
    x1 = np.load("comp_table_x1.npy")    
    prob = np.load("comp_table_p.npy")
    nt = len(t)
    tmax = 10**t[nt-1]
    tmin = 10**t[0]
    print nt,tmin,tmax
    #nx = len(x)
    nx1 = len(x1)
    npr = len(prob)
    print 't range: ',np.exp(min(t)),np.exp(max(t)),nt
    print 'x1 range: ',np.exp(min(x1)),np.exp(max(x1)),nx1
    print 'p range: ',min(prob),max(prob),npr
    print np.exp(t)
    print np.exp(x1)
    itt = np.arange(nt)
    itp = np.arange(npr)
    itx1 = np.arange(nx1)

    table = np.load("comp_table_x.npy")

    # get interpolant for x1
    ix1 = np.interp(np.log(x10),x1,itx1)
    ix1m = int(ix1)
    b = (np.log(x10)-x1[ix1m])/(x1[ix1m+1]-x1[ix1m])
    b1 = 1.-b

    for i in range(ndraw):
        #print i,path,dev,y,tau
        if (path >= tmax):
            print i,'path > tmax',path,tmax
            path= 0.999*tmax
        if (path < tmin):
            print i,'path < tmin',path,tmin
            path = 1.00001*tmin
        #path = 0.75
        t0 = np.log(path)
        it = np.interp(t0,t,itt)
        itm = int(it)
        a = (t0-t[itm])/(t[itm+1]-t[itm])
        a1 = 1.-a
        dev = ran.random()
        
        #ip = int(dev*(npr-1))
        ip = int(np.interp(dev,prob,itp))
        d = (dev-prob[ip])/(prob[ip+1]-prob[ip])
        d1 = 1.-d
        
        xdev = (b*(d*(a*table[itm+1,ix1m+1,ip+1]+a1*table[itm,ix1m+1,ip+1])
                   +d1*(a*table[itm+1,ix1m+1,ip]  +a1*table[itm,ix1m+1,ip]))
                +b1*(d*(a*table[itm+1,ix1m,ip+1]+a1*table[itm,ix1m,ip+1])
                   +d1*(a*table[itm+1,ix1m,ip]  +a1*table[itm,ix1m,ip])))
        #print 't: ',itm,a,a1,np.exp(t[itm]),np.exp(t[itm+1])
        #print 'x1: ',ix1m,b,b1,np.exp(x1[ix1m]),np.exp(x1[ix1m+1])
        #print ip,d,d1,prob[ip],prob[ip+1]
        #print xdev

        if (xdev > x0[nmax-1]):
            print i,'xdev > xmax',xdev,x0[nmax-1]
            xdev = 1.
        if (xdev < x0[0]):
            print i,'xdev < xmin',xdev,x0[0]
            print itm,a,a1,t[itm],t[itm+1]
            print ix1m,b,b1,x1[ix1m],x1[ix1m+1]
            print ip,d,d1,prob[ip],prob[ip+1],dev
            print xdev
            xdev = 0.
        xbin = bisect_right(x0,xdev)-1
        #print xdev,x0[xbin],x0[xbin+1]
        #print dev,xdev,xbin,ip
        iout[xbin] = iout[xbin] + 1.

    out = np.zeros([nmax-1,2])
    sum = 0.
    for j in range(nmax-1):
        sum = sum + iout[j]
    print sum
    for j in range(nmax-1):
        out[j,0]=xm[j]
        out[j,1]=xm[j]*iout[j]/dx0[j]/ndraw
    #out[:,1] = out[:,1] / sum *x10**2
    ncomp = get_green(path,x10,xm,fast=fast,brisk=brisk)
    sum = 0.
    for j in range(nmax-1):
        sum = sum + xm[j]**2*ncomp[j]*dx0[j]/x10**2
    print sum
    np.savetxt("x_green.out",out)
    plt.plot(xm,out[:,1]*x10**2,'.')
    plt.plot(xm,xm**3*ncomp)
    if (ylim is not None):
        plt.ylim(ylim)
    plt.xscale('log')
    plt.yscale(yscale)
    plt.ylabel(r"$I_\nu {\rm (erg/s/Hz)}$")
    plt.xlabel(r"$x$")
    plt.savefig("intens_acc.pdf")
    #plt.show()
    plt.close()

def plot_table(it,ix1,xlim=None,ylim=None):

    t = np.load("comp_table_t.npy")
    x1 = np.load("comp_table_x1.npy")    
    prob = np.load("comp_table_p.npy")
    nt = len(t)
    nx1 = len(x1)
    npr = len(prob)
    print 't range: ',np.exp(min(t)),np.exp(max(t)),nt
    print 'x1 range: ',np.exp(min(x1)),np.exp(max(x1)),nx1
    print 'p range: ',min(prob),max(prob),npr
    table = np.load("comp_table_x.npy")
   
    for i in it:
        for j in ix1:
            if ((i < nt) and (j < nx1)):
                print 't,x1: ',np.exp(t[i]),np.exp(x1[j])
                plt.plot(prob,table[i,j,:],'.')
    plt.xscale('linear')
    plt.yscale('log')
    if (xlim is not None):
        plt.xlim(xlim)
    if (ylim is not None):
        plt.ylim(ylim)
    plt.show()
            
def get_nx_spec(x10,tau=20.,temp=3.e6,ndraw=1000,urange=[0.01,100.],nmax=100,deltasrc=True):
    """
    Similar to get_x_dist but instead of binning x after drawing a final x and then
    binning, it only draws path/time and then computes the corresponding spectrum 
    directly from the saved greens function solution.  Used for testing
    """

    step = (np.log10(urange[1])-np.log10(urange[0]))/(nmax-1.)
    x0 = 10**(np.arange(nmax)*step+np.log10(urange[0]))
    iout = np.zeros(nmax)
    ioutc = np.zeros(nmax)

    c = 2.9979246e+10
    kb = 1.3806580e-16
    me = 9.1093897e-28
    if (deltasrc):
        mdev,mprob = prob_ct_delta()
    t = np.load("comp_table_t.npy")
    #x1 = np.load("comp_table_x1.npy")    
    x = np.load("comp_table_x.npy")
    nt = len(t)
    nx = len(x)
    tmin = min(t)
    tmax = max(t)
    print 't range: ',min(t),max(t),nt
    print 'x range: ',min(x),max(x),nx
    itt = np.arange(nt)
    itx = np.arange(nx)
    table = np.load("comp_table_old.npy")
    itab = table[:,0,:]

    for i in range(ndraw):
        dev = ran.random()
        if (deltasrc):
            y = get_ct_delta(dev,mdev,mprob)
            path = - np.log(y)*3./np.pi**2*tau**2*kb*temp/(me*c**2)
        else:
            if (i == 0):
                print "using distributed sources" 
            path = - np.log(dev)*3./np.pi**2*(tau+2./3.)**2*kb*temp/(me*c**2)
        #print i,path
        #nxc = get_green(path,x10,x0)
        it = np.interp(path,t,itt)
        if (path > tmax):
            print 't: ',path,tmax
        if (path < tmin):
            print 't: ',path,tmin        
        itm = int(it)
        a = (path-t[itm])/(t[itm+1]-t[itm])
        a1 = 1.-a
        for j in range(nmax):
            ix = np.interp(x0[j],x,itx)
            ixm = int(ix)
            if (ixm >= nx-1):
                ixm = nx-2
            b = (x0[j]-x[ixm])/(x[ixm+1]-x[ixm])
            b1 = 1.-b
            #interpolate n
            n = np.exp( b *(a*table[itm+1,0,ixm+1]+a1*table[itm,0,ixm+1])
                       +b1*(a*table[itm+1,0,ixm]  +a1*table[itm,0,ixm]))

            #print n,nxc[j]
            iout[j] = iout[j]+x0[j]**3*n
            #ioutc[j] = ioutc[j]+x0[j]**3*nxc[j]
        #print path[i]c
                
    out = np.zeros([nmax,2])
    for j in range(nmax):
        out[j,0]=x0[j]
        out[j,1]=iout[j]/ndraw

    # output file to be read in by plot_intens_sums() in mc_analysis.py
    np.savetxt("x_green.out",out)
    plt.plot(x0,iout)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("intens_acc.pdf")
    plt.close()

def integrand_green(u,t,x1,x):
    """
    Integrand for the greens function in terms of whittaker functions directly coputed
    with mpath library
    """
    u5 = 1j*u**0.5
    f = np.exp(-t*u)/u**0.5*np.absolute(gamma(u5-1.5)/gamma(2.*u5))**2*whitw(2,u5,x)*whitw(2,u5,x1)

    return f.real/4./np.pi

def integrand_brisk(u,t,x1,x):
    """
    Same as integrand_green but using the interpolated values for the specific whittaker 
    functions used by the code to speed up the calculation
    """
    u5 = 1j*u**0.5
    f = np.exp(-t*u)/u**0.5*np.absolute(gamma(u5-1.5)/gamma(2.*u5))**2*whitwbrisk(u,x)*whitwbrisk(u,x1)
    return f.real/4./np.pi


def gen_whit_table(nx=200,nu=400,xrange=[1.e-3,1.e3],urange=[1.e-4,1.e3]):
    """
    Generate precomputed table for computing whittaker functions
    """
    step = (np.log10(xrange[1])-np.log10(xrange[0]))/(nx-1.)
    x = 10**(np.arange(nx)*step+np.log10(xrange[0]))
    step = (np.log10(urange[1])-np.log10(urange[0]))/(nu-1.)
    u = 10**(np.arange(nu)*step+np.log10(urange[0]))
    out = np.zeros([nu,nx])
    for i in range(nu):
        u5 = 1j*u[i]**0.5
        for j in range(nx):
            print i,j
            out[i,j] = whitw(2,u5,x[j]).real

    np.save("whit_x",x)
    np.save("whit_u",u)
    np.save("whit_table",out)

ub = 0.
xb = 0.
whb = 0.
iub = 0
ixb = 0
nub = 1
nxb = 1

def init_brisk():
    """
    Initializes table for computing whittaker functions
    """
    global ub 
    ub = np.load("whit_u.npy")
    global nub 
    nub = len(ub)
    global xb 
    xb = np.load("whit_x.npy")
    global nxb 
    nxb = len(xb)
    global whb 
    whb = np.load("whit_table.npy")
    global iub 
    iub = np.arange(nub)
    global ixb 
    ixb = np.arange(nxb)

def whitwbrisk(u,x):
    """
    Impliments faste whittaker function bia interpolation on precomputed table
    """

    iu = np.interp(u,ub,iub)
    ium = int(iu)
    #print u,ium
    a = (u-ub[ium])/(ub[ium+1]-ub[ium])
    a1 = 1.-a
    ix = np.interp(x,xb,ixb)
    ixm = int(ix)
    #print x,ixm
    b = (x-xb[ixm])/(xb[ixm+1]-xb[ixm])
    b1 = 1.-b
    return (b*(a*whb[ium+1,ixm+1]+a1*whb[ium,ixm+1])
            +b1*(a*whb[ium+1,ixm]+a1*whb[ium,ixm]))


def plot_whit(x,urange=[1.e-3,1.e3],xscale='log',xlim=None):
    """
    Plotting for diagnostic purposes
    """
    nmax = 1000
    step = (np.log10(urange[1])-np.log10(urange[0]))/(nmax-1.)
    s = 10**(np.arange(nmax)*step+np.log10(urange[0]))
    y = np.zeros(nmax)
    for j in range(len(x)):
        for i in range(nmax):
            u5 = 1j*s[i]**0.5
            y[i] = whitw(2,u5,x[j]).real
            #print s[i],y[i]]
        plt.plot(s,y,'.')

    if (xlim is not None): plt.xlim(xlim)
    plt.xscale(xscale)
    plt.show()

def plot_integrand(t,x1,x,urange=[1.e-3,1.e3],xscale='log',xlim=None):
    """
    Plotting for diagnostic purposes
    """
    nmax = 1000
    step = (np.log10(urange[1])-np.log10(urange[0]))/(nmax-1.)
    s = 10**(np.arange(nmax)*step+np.log10(urange[0]))
    y = np.zeros(nmax)
    for i in range(nmax):
        y[i] = integrand_i2(s[i],t,x1,x)
        #print s[i],y[i]]
    plt.plot(s,y,'.')
    if (xlim is not None): plt.xlim(xlim)
    plt.xscale(xscale)
    plt.show()


def green_function_slow(t,x1,x,max=1.e3,min=1.e-3):
    """
    Compute greens function with mmap whittaker functions
    """
    integ = integrate.quad(integrand_green,min,max,args=(t,x1,x))
    ret = x1**2/2.*np.exp(-x)*(1.+(1-2./x)*(1.-2./x1)*np.exp(-2.*t))+np.exp(-9./4.*t)/x**2*np.exp((x1-x)/2.)*integ[0]
    return np.amax([ret,0.])


def green_function_brisk(t,x1,x,max=1.e3,min=1.e-3):
    """
    Compute greens function with pretabulated whittaker function
    """
    integ = integrate.quad(integrand_brisk,min,max,args=(t,x1,x))
    ret = x1**2/2.*np.exp(-x)*(1.+(1-2./x)*(1.-2./x1)*np.exp(-2.*t))+np.exp(-9./4.*t)/x**2*np.exp((x1-x)/2.)*integ[0]
    return np.amax([ret,0.])

def integrand_i2(u,t,x1,x):
    """
    Integrand used for computing greens function with Nargirner and Loskutov's approximation
    Only seems to be accurate for relatively low values of t (t < 0.5)
    """
    xp = x*x1
    return np.exp(-t*u**2)*np.sin(2*u*np.log(4*u/np.e/xp**0.5))

def green_function_fast(t,x1,x,max=1.e6):
    """
    Compute the approximate greens function provided by Nargner and Loskutov. Works better
    than explicit form with Whittaker functions for small t, Also faster.  But doesn't
    work at large t.
    """
    i1 = 0.5*(np.pi/t)**0.5*np.exp(-(np.log(x/x1))**2/4./t)
    i20 = integrate.quad(integrand_i2,0.,max,args=(t,x1,x),full_output=1)
    i2 = integrate.quadrature(integrand_i2,0.,max,args=(t,x1,x))
    return (x1/x)**0.5/(np.pi*x)*np.exp((x1-x)/2.-9./4*t)*(i1+i2[0])


def get_green(t,x1,x,max=None,fast=True,verbose=False,brisk=False,dx=None):
    """
    Return green functions compute over an array of final energies.  Function
    returns renormalizes n_x
    """
    nmax = len(x)
    nx = np.zeros(nmax)
    if (dx is None):
        dx = np.zeros(nmax)
        dx[0] = x[1]-x[0]
        dx[nmax-1] = x[nmax-1]-x[nmax-2]
        for i in range(1,nmax-1):
            dx[i] = 0.5*(x[i+1]-x[i-1])
        
    sum =0.
    for i in range(nmax):
        if (fast):
            if (max is None):
                max = 1.e6
            nx[i] = green_function_fast(t,x1,x[i],max=max)
        else:
            if (brisk):
                if (max is None):
                    max = 1.e3
                nx[i] = green_function_brisk(t,x1,x[i],max=max)
            else:
                if (max is None):
                    max = 1.e3
                nx[i] = green_function_slow(t,x1,x[i],max=max)
        sum = sum + dx[i]*nx[i]*x[i]**2
        #print i,sum,x[i],dx[i]
        if (verbose):
            print i,nmax,x[i],nx[i]

    if (sum <= 0.):
        print 'Sum <= 0 in get green: ',sum
        print t,x1
        norm = 1.0
    else:
        norm = x1**2/sum
        #norm = 1.
    #print t,x1,norm
    return nx*norm

def gen_table_old(nt,nx1,nx,trange=[1.e-5,10.],x1range=[0.001,100.],xrange=[0.001,100.],tswitch=0.5,verbose=False):
    """
    Generate input tables used for get_nx_spec().  For that application x1
    should be single valued
    """
    step = (np.log10(trange[1])-np.log10(trange[0]))/(nt-1.)
    t = 10**(np.arange(nt)*step+np.log10(trange[0]))
    if (nx1 == 1):
        x1 = np.zeros(1)
        x1[0] = x1range[0]
    else:
        step = (np.log10(x1range[1])-np.log10(x1range[0]))/(nx1-1.)
        x1 = 10**(np.arange(nx1)*step+np.log10(x1range[0]))
    step = (np.log10(xrange[1])-np.log10(xrange[0]))/(nx-1.)
    x = 10**(np.arange(nx)*step+np.log10(xrange[0]))
    out = np.zeros([nt,nx1,nx])
    print x1
    for i in range(nt):
        if (t[i] < tswitch):
            fast = True
        else:
            fast = False
        for j in range(nx1):
            print i,j
            n = get_green(t[i],x1[j],x,fast=fast,verbose=verbose)
            for k in range(nx):
                #if ((i == 100) and (k == 76)):
                #    print t[i],x[k],n[k]
                out[i,j,k] = np.log(max(n[k],1.e-100))
            #out[i,j,:] = np.log(get_green(t[i],x1[j],x))

    np.save("comp_table_old",out)
    np.save("comp_table_t",t)
    np.save("comp_table_x1",x1)
    np.save("comp_table_x",x)

def gen_table_depricated(nt,nx1,nx,trange=[1.e-4,100.],x1range=[0.001,100.],xrange=[0.001,100.]):
    """
    No longer used.  Implemented for testing
    """
    step = (np.log10(trange[1])-np.log10(trange[0]))/(nt-1.)
    t = 10**(np.arange(nt)*step+np.log10(trange[0]))
    if (nx1 == 1):
        x1 = np.zeros(1)
        x1[0] = x1range[0]
    else:
        step = (np.log10(x1range[1])-np.log10(x1range[0]))/(nx1-1.)
        x1 = 10**(np.arange(nx1)*step+np.log10(x1range[0]))
    
    step = (np.log10(xrange[1])-np.log10(xrange[0]))/(nx-1.)
    x = 10**(np.arange(nx)*step+np.log10(xrange[0]))
    dx =np.zeros(nx)
    for i in range(nx-1):
        dx[i] = x[i+1]-x[i]
    out = np.zeros([nt,nx1,nx])
    print x1
    for i in range(nt):
        for j in range(nx1):
            n = get_green(t[i],x1[j],x)
            sum = 0.
            out[i,j,0] = 0.
            for k in range(1,nx):
                out[i,j,k] = out[i,j,k-1]+n[k-1]*dx[k-1]
            norm = out[i,j,nx-1]
            out[i,j,:] = out[i,j,:]/norm
            print norm,out[i,j,nx-1]

    np.save("comp_table",out)
    np.save("comp_table_t",np.log(t))
    np.save("comp_table_x1",np.log(x1))
    np.save("comp_table_x",x)


def gen_table(nt,nx1,nx,npr,trange=[1.e-4,1.e4],x1range=[0.001,100.],xrange=[0.001,100.],itrange=None,tswitch=0.2,plotonly=False,plim=None):
    """
    Generated tables used by athena++ (or get_x_dist()) to draw final energy based on
    time/path, initial energy, and a uniform probability from 0 to 1.
    """
    init_brisk()
    if (nt == 1):
        t = np.zeros(1)
        t[0] = trange[0]
    else:
        step = (np.log10(trange[1])-np.log10(trange[0]))/(nt-1.)
        t = 10**(np.arange(nt)*step+np.log10(trange[0]))
        if (itrange is not None):
            t=t[itrange[0]:itrange[1]]
            nt = len(t)
            print nt
            print t
    if (nx1 == 1):
        x1 = np.zeros(1)
        x1[0] = x1range[0]
    else:
        step = (np.log10(x1range[1])-np.log10(x1range[0]))/(nx1-1.)
        x1 = 10**(np.arange(nx1)*step+np.log10(x1range[0]))

    plow=1.e-4
    phigh=1.e-4
    prob0 = np.arange(npr)/(npr-1.)
    npr0 = npr/3
    npr1 = npr*2/3
    prange=[prob0[npr0],prob0[npr1]]
    print npr0,npr1,prange
    step = (prange[1]-prange[0])/((npr1-npr0)-1.)
    probm = np.arange(npr1-npr0)*step+prange[0]
    step = (np.log10(prange[0])-np.log10(plow))/(npr0-1.)
    probl = 10**(np.arange(npr0)*step+np.log10(plow))
    step = (np.log10(1.-prange[1])-np.log10(phigh))/((npr-npr1)-1.)
    probh = np.flip(1.-10**(np.arange(npr-npr1)*step+np.log10(phigh)))
    prob = np.zeros(1)
    prob = np.append(prob,probl[0:npr0-1])
    prob = np.append(prob,probm[0:(npr1-npr0-1)])
    prob = np.append(prob,probh)
    prob = np.append(prob,[1.])
    print len(prob),prob
    #prob = np.arange(npr)/(npr-1.)
    step = (np.log10(xrange[1])-np.log10(xrange[0]))/(nx-1.)
    x = 10**(np.arange(nx)*step+np.log10(xrange[0]))
    ix = np.arange(nx)
    px = np.zeros(nx)
    dx =np.zeros(nx-1)
    xm = np.zeros(nx-1)
    #dx[0] = x[1]-x[0]
    #dx[nx-1] = x[nx-1]-x[nx-2]
    #for i in range(1,nx-1):
    #    dx[i] = 0.5*(x[i+1]-x[i-1])
    for i in range(nx-1):
        dx[i] = (x[i+1]-x[i])
        xm[i] = 0.5*(x[i+1]+x[i])

    out = np.zeros([nt,nx1,npr])

    for i in range(nt):
        for j in range(nx1):
            print i,j,t[i],x1[j]
            if (t[i] < tswitch):
                fast = True
                verbose = False
                brisk = False
            else:
                fast = False
                verbose = True
                brisk = True
            nnorm = get_green(t[i],x1[j],xm,fast=fast,brisk=brisk,verbose=verbose,dx=dx)
            px[0] = 0.
            for k in range(1,nx):
                px[k] = px[k-1] + nnorm[k-1] * xm[k-1]**2 * dx[k-1] / x1[j]**2
                #print k,px[k],x[k],dx[k]
            #print px[nx-1]
            #for k in range(nx-1):
            #    print 'p: ',k,x[k],px[k],nnorm[k]
            px[nx-1]=1
            #print 'n-2: ',px[nx-2]
            #out[i,j,0] = x[0]
            plow = 1.e-5
            ip = np.interp(plow,px,ix)
            ipm = int(ip)
            a = (plow-px[ipm])/(px[ipm+1]-px[ipm])
            a1 = 1.-a
            out[i,j,0] = (a*x[ipm+1]+a1*x[ipm])
            for k in range(1,npr-1):
                ip = np.interp(prob[k],px,ix)
                ipm = int(ip)
                a = (prob[k]-px[ipm])/(px[ipm+1]-px[ipm])
                a1 = 1.-a
                out[i,j,k] = (a*x[ipm+1]+a1*x[ipm])
                #if (t[i] > 0.2):
                #    print k,out[i,j,k],x[ipm],x[ipm+1],a,a1
            ip = np.interp(1.-plow,px,ix)
            ipm = int(ip)
            a = ((1.-plow)-px[ipm])/(px[ipm+1]-px[ipm])
            a1 = 1.-a
            out[i,j,npr-1] = (a*x[ipm+1]+a1*x[ipm])
            #out[i,j,0] = out[i,j,1]
            #out[i,j,npr-1] = out[i,j,npr-2]

    if (plotonly):
        plt.subplot(1,2,1)
        plt.plot(prob,out[0,0,:],'.')
        #print out[0,0,:]
        #plt.plot(px,x,':')
        plt.yscale('log')
        if (plim is not None):
            plt.xlim(plim)
        #plt.xscale('log')
        plt.subplot(1,2,2)
        plt.plot(xm,xm**2*nnorm)
        plt.yscale('log')
        plt.xscale('log')
        plt.show()
    else:
        write_binary("compton_table_x.out",out)
        write_binary("compton_table_t.out",np.log(t))
        write_binary("compton_table_x1.out",np.log(x1))
        write_binary("compton_table_p.out",prob)
        np.save("comp_table_x",out)
        np.save("comp_table_t",np.log(t))
        np.save("comp_table_x1",np.log(x1))
        np.save("comp_table_p",prob)

def check_solution(x1,t,nx=200,xrange=[1.e-3,1.e2],fast=True,brisk=True,verbose=False,rat=0.001):
    """
    Check solution of Kompaneets equation
    """
    if (brisk):
        init_brisk()

    if (not fast):
        verbose = True
    # evaluate function
    step = (np.log10(xrange[1])-np.log10(xrange[0]))/(nx-1.)
    x = 10**(np.arange(nx)*step+np.log10(xrange[0]))
    n0 = get_green(t,x1,x,fast=fast,brisk=brisk,verbose=verbose)

    # evaluate time derivative
    dt = rat*t
    nl = get_green(t-dt,x1,x,fast=fast,brisk=brisk,verbose=verbose)
    nu = get_green(t+dt,x1,x,fast=fast,brisk=brisk,verbose=verbose)
    dndt = 0.5*(nu-nl)/dt

    #next evaluate spatial derivatives and rhs
    dndx=np.zeros(nx)
    dndx[0]=(n0[1]-n0[0])/(x[1]-x[0])
    for i in range(1,nx-1):
        dndx[i] = (n0[i+1]-n0[i-1])/(x[i+1]-x[i-1])
    dndx[nx-1] = (n0[nx-1]-n0[nx-2])/(x[nx-1]-x[nx-2])
    paren = x**4*(dndx+n0)
    rhs=np.zeros(nx)
    rhs[0]=(paren[1]-paren[0])/(x[1]-x[0])
    for i in range(1,nx-1):
        rhs[i] = (paren[i+1]-paren[i-1])/(x[i+1]-x[i-1])
    rhs[nx-1]=(paren[nx-1]-paren[nx-2])/(x[nx-1]-x[nx-2])
    rhs = rhs/x**2
    
    out = np.zeros([nx,6])
    out[:,0] = x
    out[:,1] = n0
    out[:,2] = dndt
    out[:,3] = rhs
    out[:,4] = nl
    out[:,5] = nu
    np.save("check_sol",out)

def plot_check_sol(plot=0,infile="check_sol.npy",ylim=None,yscale=None):

    data = np.load(infile)
    x = data[:,0]
    n0 = data[:,1]
    dndt = data[:,2]
    rhs = data[:,3]
    nl = data[:,4]
    nu = data[:,5]
    if (plot == 0):
        plt.plot(x,(rhs-dndt),'.')
        plt.plot(x,n0)
    if (plot == 1):
        plt.plot(x,x**2*n0,'.')
        plt.plot(x,x**2*nl,'.')
        plt.plot(x,x**2*nu,'.')
        if (yscale is not None):
            plt.yscale(yscale)
        else:
            plt.yscale('log')
    if (plot == 2):
        #print rhs
        plt.plot(x,rhs,'.')
        plt.plot(x,dndt,'.')
        plt.plot(x,n0)
    if (ylim is not None):
        plt.ylim(ylim)
    plt.xscale('log')
    nx = len(x)
    dndx=np.zeros(nx)
    dndx[0]=(n0[1]-n0[0])/(x[1]-x[0])
    for i in range(1,nx-1):
        dndx[i] = (n0[i+1]-n0[i-1])/(x[i+1]-x[i-1])
    dndx[nx-1] = (n0[nx-1]-n0[nx-2])/(x[nx-1]-x[nx-2])
    paren = x**4*(dndx+n0)
    print paren
    plt.plot(x,paren)
    rhs=np.zeros(nx)
    rhs[0]=(paren[1]-paren[0])/(x[1]-x[0])
    for i in range(1,nx-1):
        rhs[i] = (paren[i+1]-paren[i-1])/(x[i+1]-x[i-1])
    rhs[nx-1]=(paren[nx-1]-paren[nx-2])/(x[nx-1]-x[nx-2])
    rhs = rhs/x**2
    
    plt.show()

def write_binary(filename,array):
    """
    Write binary files that can be read by Athena++
    """
    file = open(filename, 'wb')
    bout = struct.pack('=%sd' % array.size, *array.flatten('F'))
    file.write(bout);
    file.close()

def read_binary(filename,size):
    """
    Read in files written with write_binary().  Used to check consistancy
    between numpy and binary output file tables
    """
    file = open(filename, 'rb')
    bin = file.read()
    file.close()
    return struct.unpack('=%sd' % size,bin)


def combine_tables(infiles,outfile="comp_table",appendt=True,appendx=False,appendp=False):
    """
    Unfinished utility for combining output tables.  Implemented because
    table computation can be very long on single core.  Currently only
    appends tables as function of time
    """
    nf = len(infiles)
    nt = np.zeros(nf,dtype=int)
    nx1 = np.zeros(nf,dtype=int)
    npr = np.zeros(nf,dtype=int)

    for i in range(nf):
        t = np.load(infiles[i]+"_t.npy")
        x1= np.load(infiles[i]+"_x1.npy")
        prob = np.load(infiles[i]+"_p.npy")
        nt[i] = len(t)
        nx1[i] = len(x1)
        npr[i] = len(prob)
        print nt[i],nx1[i],npr[i]
        if (i == 0):
            tc = t
            x1c = x1
            pc = prob
        else:
            if (appendt):
                tc = np.append(tc,t)

    print np.exp(tc)
    ntc = len(tc)
    nx1c = len(x1c)   
    npc = len(pc)
    print 'final: ',ntc,nx1c,npc
    outc = np.zeros([ntc,nx1c,npc])
    ic = -1
    jc = 0
    kc = 0
    for ifile in range(nf):
        out = np.load(infiles[ifile]+"_x.npy")
        for i in range(nt[ifile]):
            ic = ic + 1
            for j in range(nx1[ifile]):
                for k in range(npr[ifile]):
                    outc[ic,j,k] = out[i,j,k]

    np.save(outfile+"_t",tc)
    np.save(outfile+"_x1",x1c)
    np.save(outfile+"_p",pc)
    np.save(outfile+"_x",outc)


def convert_table(modify=False):
    """
    Utitility conver numpy file tables to binary tables
    """
    t = np.load("comp_table_t.npy")
    x1 = np.load("comp_table_x1.npy")    
    out = np.load("comp_table_x.npy")
    prob = np.load("comp_table_p.npy")
    if (modify):
        nt = len(t)
        nx1 = len(x1)
        npr = len(prob)
        for i in range(nt):
            for j in range(nx1):
                out[i,j,npr-1] = out[i,j,npr-2]

    write_binary("compton_table_x.out",out)
    write_binary("compton_table_t.out",t)
    write_binary("compton_table_x1.out",x1)
    write_binary("compton_table_p.out",prob)

def plot_green(t,x1,urange=[0.01,100.],nmax=50,max=1.e3,linear=False,brisk=False,yscale='linear',ylim=None,normx=True,compare=False,verbose=False):
    """
    Plots greens function for diagnostic purposes.  If compare is true, it will compare mmap method
    with Nargirner method
    """
    if (brisk):
        init_brisk()
    if (linear):
        x = np.arange(nmax)*(urange[1]-urange[0])/(nmax-1)+urange[0]
    else:
        step = (np.log10(urange[1])-np.log10(urange[0]))/(nmax-1.)
        x = 10**(np.arange(nmax)*step+np.log10(urange[0]))
    nx = np.zeros(nmax)
    #nx1 = np.zeros(nmax)
    dx =np.zeros(nmax)
    if (linear):
        dx[0] = x[1]-x[0]
        dx[nmax-1] = x[nmax-1]-x[nmax-2]
        for i in range(1,nmax-1):
            dx[i] = 0.5*(x[i+1]-x[i-1])
    else:
        dx[0] = np.log(x[1]/x[0])
        dx[nmax-1] = np.log(x[nmax-1]/x[nmax-2])
        for i in range(1,nmax-1):
            dx[i] = 0.5*np.log(x[i+1]/x[i-1])

    out = np.zeros(nmax)
    if (compare):
        outc = np.zeros(nmax)
    nt = len(t)
    norm = np.zeros(nt)
    for j in range(nt):
        if (normx):
            out = x**2*get_green(t[j],x1,x,verbose=verbose)
        else:
            out = get_green(t[j],x1,x,verbose=verbose)
        plt.plot(x,out)
        if (compare):
            if (normx):
                outc = x**2*get_green(t[j],x1,x,fast=False,brisk=brisk,verbose=verbose,max=max)
            else:
                outc = get_green(t[j],x1,x,fast=False,brisk=brisk,verbose=verbose,max=max)
            plt.plot(x,outc)

    #plt.plot(t,np.log(norm))
    #pow = np.zeros(nt)
    #for i in range(nt):
        #pow[i] = t[i]*9./4.-0.65
    #    pow[i] = t[i]*9./4.
    #plt.plot(t,pow)
    #plt.plot(x,x**3*nx1)
    if (not linear):
        plt.xscale('log')
    plt.yscale(yscale)
    #plt.yscale('log')
    if (ylim is not None):
        plt.ylim(ylim)
    plt.savefig("green.pdf")
    plt.close()
    plt.show()

def plot_green_brisk(t,x1,urange=[0.01,100.],nmax=50,max=1.e3,linear=False,yscale='linear',ylim=None,normx=True,compare=True,verbose=False):
    """
    Plots greens function for diagnostic purposes.  If compare is true, it will compare 
    the mmap method with the interpolation (brisk) method
    """
    init_brisk()
    print ub
    print xb
    print nxb
    print nub
    if (linear):
        x = np.arange(nmax)*(urange[1]-urange[0])/(nmax-1)+urange[0]
    else:
        step = (np.log10(urange[1])-np.log10(urange[0]))/(nmax-1.)
        x = 10**(np.arange(nmax)*step+np.log10(urange[0]))
    nx = np.zeros(nmax)
    dx =np.zeros(nmax)
    if (linear):
        dx[0] = x[1]-x[0]
        dx[nmax-1] = x[nmax-1]-x[nmax-2]
        for i in range(1,nmax-1):
            dx[i] = 0.5*(x[i+1]-x[i-1])
    else:
        dx[0] = np.log(x[1]/x[0])
        dx[nmax-1] = np.log(x[nmax-1]/x[nmax-2])
        for i in range(1,nmax-1):
            dx[i] = 0.5*np.log(x[i+1]/x[i-1])

    out = np.zeros(nmax)
    outc = np.zeros(nmax)
    nt = len(t)
    norm = np.zeros(nt)
    for j in range(nt):
        if (normx):
            out = x**2*get_green(t[j],x1,x,fast=False,brisk=True,verbose=verbose)
        else:
            out = get_green(t[j],x1,x,fast=False,brisk=True,verbose=verbose)

        plt.plot(x,out)
        print 'out: ',out
        if (compare):
            if (normx):
                outc = x**2*get_green(t[j],x1,x,fast=False,brisk=False,verbose=verbose,max=max)
            else:
                outc = get_green(t[j],x1,x,fast=False,brisk=False,verbose=verbose,max=max)
            plt.plot(x,outc)
            print 'outc: ',outc
    #plt.plot(t,np.log(norm))
    #pow = np.zeros(nt)
    #for i in range(nt):
        #pow[i] = t[i]*9./4.-0.65
    #    pow[i] = t[i]*9./4.
    #plt.plot(t,pow)
    #plt.plot(x,x**3*nx1)
    if (not linear):
        plt.xscale('log')
    plt.yscale(yscale)
    #plt.yscale('log')
    if (ylim is not None):
        plt.ylim(ylim)
    plt.savefig("green.pdf")
    plt.close()
    plt.show()



def plot_green_file(infiles=['x_green.out'],xscale='log',yscale='log',norm=None,ylim=None,xlim=None):
    """
    Plots greens_function from file for diagnostic purposes
    """
    i=0
    for file in infiles:
        inp = np.loadtxt(file)
        kb = 1.3806580e-16
        #energy = 3.e6*kb*inp[:,0]/1.6021772e-12/1000.
        energy = inp[:,0]
        if norm is not None:
            inp[:,1] = inp[:,1]*norm[i]
        plt.plot(energy,energy*inp[:,1],'.')
        i = i + 1

    plt.xscale(xscale)
    plt.yscale(yscale)
    if (xlim is not None):
        plt.xlim(xlim)
    if (ylim is not None):
        plt.ylim(ylim)
    plt.savefig("green.pdf")
    plt.close()
    plt.show()


def rad_dist(x,y,nmax=1000):
    """
    Computes distribution of final radii 
    """
    if (y == 0.):
        nmax = 1000
    else:
        if (y == 1.):
            nmax = 1000
        else:
            nmax = int(np.log(1.e-8)/np.log(y))
    if (nmax < 1000):
        nmax = 1000

    prob = 0.
    for n in range(1,nmax):
        prob = prob + x*np.sin(n*np.pi*x)*n*y**(n*n)
    return prob

def gen_table_time(ny,nx,npr):
    """
    Generates table for determining file radius for diffusion over a given
    time interval
    """
    y = np.arange(ny)/(ny-1.)
    prob = np.arange(npr)/(npr-1.)
    x = np.arange(nx)/(nx-1.)
    ix = np.arange(nx)
    dx = 1./(nx-1.)
    px = np.zeros(nx)
    out = np.zeros([ny,npr])
    for i in range(ny):
        px[0] = 0.
        for j in range(1,nx):
            print i,j,y[i],x[j]
            px[j] = px[j-1] + dx*rad_dist(x[j],y[i])
        norm = 1./px[nx-1]
        px = px*norm
        out[i,0] = x[0]
        for k in range(1,npr-2):
            ip = np.interp(prob[k],px,ix)
            ipm = int(ip)
            a = (prob[k]-px[ipm])/(px[ipm+1]-px[ipm])
            a1 = 1.-a
            out[i,k] = (a*x[ipm+1]+a1*x[ipm])
        out[i,npr-1] = x[nx-1]

    write_binary("radius_table_r.out",out)
    write_binary("radius_table_p.out",prob)
    write_binary("radius_table_t.out",y)

    np.save("radius_table_r",out)
    np.save("radius_table_p",prob)
    np.save("radius_table_t",y)


# The following are now unused functions used to compute and plot greens functions or 
# spectra for certain limits and assumptions.
def whittaker2(alpha,x):
    return (np.pi/x)**0.5*((8*alpha**3+x**3+4*alpha**2*(3.+2.*x)+4.*alpha*(1+x+x**2))*kv(-0.5-alpha,0.5*x)+x*(2.*alpha**2+x**2+2.*alpha*(1+x))*kv(0.5-alpha,0.5*x))/np.sin(alpha*np.pi)/2/gamma(1.-alpha)/gamma(alpha)

def whittaker(a,b,x):
    return np.exp(-x/2)*x**(b+0.5)*hyperu(0.5+b-a,1.+2.*b,x)


def compton_spec(x,x0,alpha):

    if (x <= x0):
        return alpha*(alpha+3.)/(2*alpha+3.)*(x/x0)**(alpha+3.)
    else:
        intfun = np.exp(-x/2)*np.pi**0.5*x**(0.5+alpha)*((8*alpha**3+x**3+4*alpha**2*(3.+2.*x)+4.*alpha*(1+x+x**2))*kv(-0.5-alpha,0.5*x)+x*(2.*alpha**2+x**2+2.*alpha*(1+x))*kv(0.5-alpha,0.5*x))/np.sin(alpha*np.pi)/2/gamma(1.-alpha)
        return intfun*alpha*(alpha+3.)/gamma(2.*alpha+4.)*(x0/x)**alpha

def compton_spec2(x,x0,alpha):

    if (x <= x0):
        return (alpha+3.)*gamma(alpha+1.)/gamma(2.*alpha+4.)*np.exp(x0/2)/x0**3*hyp1f1(alpha,2.*alpha+4.,x)*x**(3.+alpha)*np.exp(-x)*whittaker(2.,alpha+1.5,x0)
    else:
        return (alpha+3.)*gamma(alpha+1.)/gamma(2.*alpha+4.)*hyp1f1(alpha,2*alpha+4,x0)*x0**(alpha-1.)*x*np.exp(-x/2.)*whittaker(2.,alpha+1.5,x)

def compton_spec3(x,x0,alpha):


    low = (alpha+3.)*gamma(alpha+1.)/gamma(2.*alpha+4.)*np.exp(x0/2)/x0**3*hyp1f1(alpha,2.*alpha+4.,x)*x**(3.+alpha)*np.exp(-x)*whittaker(2.,alpha+1.5,x0)

    high = (alpha+3.)*gamma(alpha+1.)/gamma(2.*alpha+4.)*hyp1f1(alpha,2*alpha+4,x0)*x0**(alpha-1.)*x*np.exp(-x/2.)*whittaker(2.,alpha+1.5,x)

    return 1./(1./low**0.5+1./high**0.5)**2
