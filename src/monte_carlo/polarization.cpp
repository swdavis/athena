//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarization.cpp
//! \brief moving a photon into and out of the polarized scattering basis
//
// See polarization.hpp for the convention these two functions implement.

// C++ headers
#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "polarization.hpp"
#include "montecarlo.hpp"
#include "photon.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ToScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//! \brief rotate the photon into the basis the polarized scattering routines assume
//
// On a spherical grid the wavevector is stored in the local orthonormal spherical basis,
// which rotates along the ray; the scattering routines want cartesian components.  Q and
// U are deliberately untouched: they are already referenced to the global cartesian
// meridian, which does not move.


void ToScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  if (pmcb->topology != MCTOPO_SPHERICAL) return;

  Real cth = std::cos(pphot->x2p[ip]);
  Real sth = std::sin(pphot->x2p[ip]);
  Real cph = std::cos(pphot->x3p[ip]);
  Real sph = std::sin(pphot->x3p[ip]);

  Real k3[3];
  k3[0] = pphot->k1p[ip];
  k3[1] = pphot->k2p[ip];
  k3[2] = pphot->k3p[ip];

  pphot->k1p[ip] = k3[0]*sth*cph + k3[1]*cth*cph - k3[2]*sph;
  pphot->k2p[ip] = k3[0]*sth*sph + k3[1]*cth*sph + k3[2]*cph;
  pphot->k3p[ip] = k3[0]*cth     - k3[1]*sth;

  Real knorm = std::sqrt(SQR(pphot->k1p[ip]) + SQR(pphot->k2p[ip]) + SQR(pphot->k3p[ip]));
  pphot->k1p[ip] /= knorm;
  pphot->k2p[ip] /= knorm;
  pphot->k3p[ip] /= knorm;
}

//----------------------------------------------------------------------------------------
//! \fn void FromScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//! \brief rotate the photon back out of the scattering basis
//
// Exact inverse of ToScatteringBasis, evaluated at the same position, so the round trip
// is the identity to round-off.

void FromScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  if (pmcb->topology != MCTOPO_SPHERICAL) return;

  Real cth = std::cos(pphot->x2p[ip]);
  Real sth = std::sin(pphot->x2p[ip]);
  Real cph = std::cos(pphot->x3p[ip]);
  Real sph = std::sin(pphot->x3p[ip]);

  Real k3[3];
  k3[0] = pphot->k1p[ip];
  k3[1] = pphot->k2p[ip];
  k3[2] = pphot->k3p[ip];

  pphot->k1p[ip] =  k3[0]*sth*cph + k3[1]*sth*sph + k3[2]*cth;
  pphot->k2p[ip] =  k3[0]*cth*cph + k3[1]*cth*sph - k3[2]*sth;
  pphot->k3p[ip] = -k3[0]*sph     + k3[1]*cph;

  Real knorm = std::sqrt(SQR(pphot->k1p[ip]) + SQR(pphot->k2p[ip]) + SQR(pphot->k3p[ip]));
  pphot->k1p[ip] /= knorm;
  pphot->k2p[ip] /= knorm;
  pphot->k3p[ip] /= knorm;
}
