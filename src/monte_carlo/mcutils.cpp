//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mcutils.cpp
//  \brief implementation of utility functions


// Athena++ headers
#include "mcutils.hpp"
#include <complex>

//------------------------------------------------------------------------------
//! \fn Real BessI0(Real x)
//  \brief Evaluate modified Bessel function In(x) and n=0

Real BessI0(Real x)
{
  Real ax, ans;
  Real y;

  if ((ax=fabs(x)) < 3.75) {
    y = x/3.75,y=y*y;
    ans = 1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
          +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    ans = (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
          +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
          +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
          +y*0.392377e-2))))))));
  }
  return ans;
}

//------------------------------------------------------------------------------
//! \fn Real BessI1(Real x)
//  \brief Evaluate modified Bessel function In(x) and n=1.

Real BessI1(Real x)
{
  Real ax, ans;
  Real y;

  if ((ax=fabs(x)) < 3.75) {
    y = x/3.75,y=y*y;
    ans = ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
          +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  } else {
    y = 3.75/ax;
    ans = 0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
          -y*0.420059e-2));
    ans =  0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
           +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}

//------------------------------------------------------------------------------
//! \fn Real BessK0(Real x)
//  \brief Evaluate modified Bessel function Kn(x) and n=0.

Real BessK0(Real x)
{
  Real y,ans;

  if (x <= 2.0) {
    y = x*x/4.0;
    ans = (-log(x/2.0)*BessI0(x))+(-0.57721566+y*(0.42278420
          +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
          +y*(0.10750e-3+y*0.74e-5))))));
  } else {
    y = 2.0/x;
    ans = (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
          +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
          +y*(-0.251540e-2+y*0.53208e-3))))));
  }
  return ans;
}

//------------------------------------------------------------------------------
//! \fn Real BessK1(Real x)
//  \brief Evaluate modified Bessel function Kn(x) and n=1

Real BessK1(Real x)
{
  Real y,ans;

  if (x <= 2.0) {
    y = x*x/4.0;
    ans = (log(x/2.0)*BessI1(x))+(1.0/x)*(1.0+y*(0.15443144
          +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
          +y*(-0.110404e-2+y*(-0.4686e-4)))))));
  } else {
    y = 2.0/x;
    ans = (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
          +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
          +y*(0.325614e-2+y*(-0.68245e-3)))))));
  }
  return ans;
}

//------------------------------------------------------------------------------
//! \fn Real BessK(int n, Real x)
//  \brief Evaluate modified Bessel function Kn(x) and n >= 0
//
// Note that for x == 0 the functions bessy and bessk are not
// defined and a blank is returned.

Real BessK(int n, Real x) {

  if (n < 0 || x == 0.0) {
    printf("Error in bessel function\n");
  }
  if (n == 0)
    return( BessK0(x) );
  if (n == 1)
    return( BessK1(x) );
 
  Real tox = 2.0/x;
  Real bkm = BessK0(x);
  Real bk = BessK1(x);
  for (int j=1; j<n; ++j) {
    Real bkp = bkm+j*tox*bk;
    bkm = bk;
    bk = bkp;
  }
  return bk;
}

//----------------------------------------------------------------------------------------
//! \fn int mcbisect(Real x, AthenaArray<Real> &array)
// use bisection to search array and return bin

int mcbisect(Real x, AthenaArray<Real> &array) {

  int nmax = array.GetDim1();
  int low = 0, high = nmax-1, mid;

  while(low<=high) {
    mid=(low+high)/2;
    if(array(mid-1) <= x) {
      if(array(mid) > x)
        break;
      else
        low=mid+1;
    }
    else
      high=mid-1;
  }

  if (low >= nmax-1)
    return nmax-2;
  else
    return std::max(mid-1,0);

}


//----------------------------------------------------------------------------------------
//! \fn void ZetaFast(std::complex<double> arg, std::complex<double> &zeta, std::complex<double> &dzeta)
// use bisection to search array and return bin

void ZetaFast(std::complex<double> arg, std::complex<double> &zeta) {

  std::complex<double> b=(0.5,0.80558);
  std::complex<double> a=(0.50556,-0.81462);
  std::complex<double> d2=(0.,7.08981540362206);
  std::complex<double> aux0 = arg;
  if(aux0.imag() < 0.) 
    aux0 = std::conj(aux0);
  zeta = b / (a-aux0) - std::conj(b) / (std::conj(a)+aux0);
  if (aux0.imag() < 0.)
    zeta = std::conj(zeta) + d2 * exp(-arg*arg);

}
