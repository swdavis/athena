//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file atmos.cpp
//! \brief A hydrostatic atmosphere that tracks ionization using a passive scalar.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

void TrackIonization(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Spherical atmosphere in hydrostatic balance
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserExplicitSourceFunction(TrackIonization);
  // Add user mesh data block
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real kb    = 1.380649e-16;
  Real mp    = 1.6726e-24;
  Real clight = 2.997924589e10;

  Real GM   = pin->GetOrAddReal("problem", "GM", 0.);
  Real mmw   = pin->GetOrAddReal("problem", "mmw", 1.); // Mean molecular weight of gas
  Real rin   = pin->GetOrAddReal("problem", "rin", 1.e10); // Radius of base of atmosphere
  Real pbase = pin->GetOrAddReal("problem", "pbase", 10.); // Atmospheric pressure at rin
  Real temp  = pin->GetOrAddReal("problem", "temp", 1.e4); // Isothermal temperature
  Real rhobase = pbase * mmw * mp / kb / temp;
  Real gamma = peos->GetGamma();
  Real gm1   = gamma - 1.0;

  // set up ambient medium at equilibrium for an isothermal atmosphere
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is-NGHOST; i<=ie+NGHOST; i++) {
        Real r = pcoord->x1v(i);
        phydro->u(IDN,k,j,i) = rhobase * exp(GM * mmw * mp / kb / temp * (1./pcoord->x1v(i) - 1./rin));
        Real rho = phydro->u(IDN,k,j,i);
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = phydro->u(IDN,k,j,i) * kb * temp / mmw / mp / gm1; // Fluid internal energy density
        phydro->u(IEN,k,j,i) += 0.5 * SQR(phydro->u(IM1,k,j,i)) / phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5 * SQR(phydro->u(IM2,k,j,i)) / phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5 * SQR(phydro->u(IM3,k,j,i)) / phydro->u(IDN,k,j,i);

        // Calculate photoionization rate
        Real Gamma = 1. / (6. * 60. * 60.); // Hardcoded for now
        Real alpha = 4.18e-13;

        // Set s to be rho times ionization frac (=np * mp)
        pscalars->s(0,k,j,i) = 0.5 * (std::sqrt(Gamma) * std::sqrt(4.*alpha*rho/mp + Gamma) / alpha
                               - Gamma / alpha) * mp;
        pscalars->r(0,k,j,i) = pscalars->s(0,k,j,i) / rho;
      }
    }
  }
}

void TrackIonization(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar) {
  Real mp = 1.6726e-24;
  Real sigma_pi = 6.3e-18;// cm2
  Real Gamma0 = 1. / (6. * 60. * 60.); // s, photoionization rate coefficient
  Real alpha = 4.18e-13; // cm3 s-1 for T=1e4 K, recombination rate coefficient
  Real nR = 1. / alpha / dt;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
        Real rho = prim(IDN,k,j,i);
        Real n_p = cons_scalar(0,k,j,i)/mp;
        Real n_H = rho/mp - n_p;
        Real dr = pmb->pcoord->dx1f(i);
        Real Gamma = Gamma0 / (1. + std::pow(sigma_pi * n_H*dr, 1.5)); // Trammell et al 2011, Fig 9 powerlaw
        printf("%g\n", Gamma);
        Real nC = Gamma / alpha;
        cons_scalar(0,k,j,i) = 0.5 * mp * (-(nC + nR) + std::sqrt((nC + nR)*(nC + nR) + 4. * (nC*rho/mp + nR*n_p)));
      }
    }
  }
  return;
}
