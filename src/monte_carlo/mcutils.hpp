#ifndef MCUTILS_HPP
#define MCUTILS_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mcutils.hpp
//! \brief prototypes for monte carlo related utility functions
//
// Athena++ classes headers
#include "montecarlo.hpp"
#include <complex>

// function prototypes
Real BessI0(Real x);
Real BessI1(Real x);
Real BessK0(Real x);
Real BessK1(Real x);
Real BessK(int n, Real x);
int mcbisect(Real x, AthenaArray<Real> &array);
int mcbisect(Real x, Real *array, int dim);
std::complex<Real> ZetaFast(std::complex<Real> arg);
std::complex<Real> ZetaVoigt(std::complex<Real> arg);

#endif // MCUTILS_HPP
