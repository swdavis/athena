//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mcutils.cpp
//! \brief implementation of monte carlo utility functions

// Athena++ headers
#include "mcutils.hpp"
#include <complex>

//------------------------------------------------------------------------------
//! \fn Real BessI0(Real x)
//! \brief Evaluate modified Bessel function In(x) and n=0

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
//! \brief Evaluate modified Bessel function In(x) and n=1.

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
//! \brief Evaluate modified Bessel function Kn(x) and n=0.

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
//! \brief Evaluate modified Bessel function Kn(x) and n >= 0
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
//! \brief  use bisection to search array and return bin

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
//! \fn int mcbisect(Real x, AthenaArray<Real> &array)
//! \brief  use bisection to search array and return bin

int mcbisect(Real x, Real *array, int dim) {

  int nmax = dim;
  int low = 0, high = nmax-1, mid;

  while(low<=high) {
    mid=(low+high)/2;
    if(array[mid-1] <= x) {
      if(array[mid] > x)
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
//! \fn std::complex<Real> ZetaFast(std::complex<Real> arg)
//! \brief fast but approximate evaluation of zeta function

std::complex<Real> ZetaFast(std::complex<Real> arg) {

  std::complex<Real> b=(0.5,0.80558);
  std::complex<Real> a=(0.50556,-0.81462);
  std::complex<Real> d2=(0.,7.08981540362206);
  std::complex<Real> aux0 = arg;
  if(aux0.imag() < 0.)
    aux0 = std::conj(aux0);
  std::complex<Real> zeta = b / (a-aux0) - std::conj(b) / (std::conj(a)+aux0);
  if (aux0.imag() < 0.)
    zeta = std::conj(zeta) + d2 * exp(-arg*arg);
  return zeta;
}

//----------------------------------------------------------------------------------------
//! \fn std::complex<Real> ZetaVoigt(std::complex<Real> arg)
//! \brief fast but approximate evaluation of zeta function

std::complex<Real> ZetaVoigt(std::complex<Real> arg) {
  // Attempted C++ conversion of the code P. Arras obtained from Hammet
  // SWD: not tested with Real = float -- may not work?
  std::complex<Real> d1(0.,1.77245385090551); //   i sqrt(pi)
  std::complex<Real> d2(0.,3.54490770181103); // 2 i sqrt(pi)
  std::complex<Real> d3(0.,7.08981540362206); // 4 i sqrt(pi)
  std::complex<Real> aux0 = arg;
  Real eps = 1.e-7;

  //Authors- Bill Sharp modified for PDP-10 by C. F. F. Karney.
  // (This code has evolved from something Gary Swanson wrote in the
  // mid-60's!  Hence the archaic fortran branch commands.)
  // Modified for the VAX by Greg Hammett (also fixed the asymptotic
  // formula choice of imaginary part and put in an improved power series
  // solution.)  This routine has been checked against the tables in
  // Fried-Conte.
  //
  //     This routine computes the fried-conte plasma dispersion function,
  // f(z), where
  //       f(z) = 1 / sqrt(pi) * def. integral from minus infinity
  //       infinity of exp( -t * t) / (t - z) dt.
  // an equivalent formulation is
  //       f(z) = sqrt(pi) * i * w(z)
  //       where
  //       w(z) = exp(-z * z) * erfc( -i * z) .
  // here erfc (z) = complimentary error function given by
  //       erfc (z) = 2/sqrt(pi) * def. integral from z to infinity of
  //       exp(- t*t) dt.
  // there are three different approximation methods used depending on the
  // value of the argument of the function.  the methods are:
  //       1. continued fraction method: abs(y).ge.1.0
  //       2. asymtotic series method: abs(x).ge.4.0 and abs(y).lt.1.0
  //       3. power series method: abs(x).lt.4.0 and abs(y).lt.1.0.
  // the routine is accurate to 1.e-07.  a new version of the code
  // accurate to 1.e-11 is being developed and is available upon request.
  //     the routine also computes the first derivative of the pdf after
  // the pdf has been computed.  to obtain the first derivative a call
  // must be made to the function dzeta which is an entry point in the
  // subprogram zeta.

  std::complex<Real> z = arg;
  std::complex<Real> zsq = SQR(z);

  if (fabs(z.imag()) >= 1.) {
    // use continued fraction method

    // GWH:  I don't understand the continued fraction method, nor where these
    // formulas came from.  Fried-Conte, Barberio-Coresetti, and Abramowitz and
    // Stegun all have what appear to be different formulas.  But I have
    // checked the results against the tables in Fried-Conte, and against the
    // other two methods along the borders.

    std::complex<Real> z0 = z;
    // the following operate on the conjugate of z
    if (z.imag() < 0.)
      z = std::conj(z);
    std::complex<Real> aux1 = 1.5 - SQR(z);
    std::complex<Real> aux2(0.,0.);
    std::complex<Real> del(1.5,0.);
    std::complex<Real> a1(0.,0.);
    std::complex<Real> a2(-1.,0.);
    std::complex<Real> b1(1.,0.);
    std::complex<Real> b2 = aux1;
    std::complex<Real> c1 = a2 / b2;

    aux1 += 2.;
    aux2 -= del;
    del += 2.;
    std::complex<Real> a3 = aux1 * a2 + aux2 * a1;
    std::complex<Real> b3 = aux1 * b2 + aux2 * b1;
    std::complex<Real> c2 = a3 / b3;
    std::complex<Real> c3 = c2 - c1;

    while (fabs(c3.real())+fabs(c3.imag()) >= eps) {
      a1 = a2;
      a2 = a3;
      b1 = b2;
      b2 = b3;
      c1 = c2;
      aux1 += 2.;
      aux2 -= del;
      del += 2.;
      a3 = aux1 * a2 + aux2 * a1;
      b3 = aux1 * b2 + aux2 * b1;
      c2 = a3 / b3;
      c3 = c2 - c1;
    }
    if (z0.imag() < 0.0) {
      z = z0;
      c2 = std::conj(c2) - d3 * z * exp(-zsq);
    }

    return -(0.5*c2 + 1.) / z;

  } else {
    // ASYMPTOTIC SERIES METHOD: ABS(X).GE.4.0 AND ABS(Y).LT.1.0
    Real xmag = fabs(z.real());
    if (xmag-4. >= 0.) {
      std::complex<Real> term = 1./z;
      std::complex<Real> aux0 = -term;
      std::complex<Real> aux1 = 0.5 * SQR(term);

      std::complex<Real> p(1.0,0.);

      // Use the proper Stokes lines, the Fried and Conte choices are wrong.
      // A derivation of the proper choice of Stokes' lines for the asymptotic
      // formulas can be found in Stix's Plasma Waves book, Miyamoto, or
      // Greg Hammett's notes. in practice this doesn't matter since
      // cexp(-z**2) is so small:
      if (z.imag() < -xmag)
        aux0 += d2 * exp(-SQR(z));
      if ((z.imag() >= -xmag) && (z.imag() < xmag))
        aux0 += d1 * exp(-SQR(z));

      do {
        term = aux1 * term * p;
        aux0 -= term;
        p += 2.;
        //printf("%e %e\n",term.real(),term.imag());
        //printf("%e %e\n",aux0.real(),aux0.imag());
      } while (fabs(term.real())+fabs(term.imag()) >= eps);

      return aux0;

    } else {
      // POWER SERIES METHOD: ABS(X).LT.4.0 AND ABS(Y).LT.1.0

      // instead of the Fried and Conte power series, use the one suggested by
      // Barberio-Corsetti in MATT-773 (1970) based on the other of the two
      // power series representations for the error function.  The Fried and
      // Conte power series has numerical problems for x as large as 4.0.

      Real nterm = 0.;
      std::complex<Real> aux0(1.,0.);
      std::complex<Real> aux1(1.,0.);
      std::complex<Real> term;
      do {
        nterm += 1.;
        aux1 = aux1 * SQR(z) / nterm;
        term = aux1 / (2.*nterm+1.);
        aux0 += term;
      } while (fabs(term.imag())+fabs(term.real()) >= eps);
      return exp(-SQR(z)) * (d1 - 2. * z * aux0);
    }
  }

}
