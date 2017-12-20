//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarwedge.cpp
//  \brief implementation of polar wedge exchange BCs in x2 direction

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void PolarWedgeInnerExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &src,
//                         const Real time, const Real dt,
//                         int is, int ie, int js, int je, int ks, int ke)
//  \brief polar wedge boundary conditions, inner x2 boundary

void PolarWedgeInnerExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &src,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Real sign;
  // add coupling variables to active mesh, reflecting v2
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
#pragma simd
        for (int n=0; n<(src.GetDim4()); ++n) {
          sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
          if (src.GetDim4()==6) // pressure tensor
            sign = flip_across_pole_ptensor[n] ? -1.0 : 1.0;
          src(k,js+j-1,i,n) += sign * src(k,js-j,i,n);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void PolarWedgeOuterExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &src,
//                          const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief polar wedge boundary conditions, outer x2 boundary

void PolarWedgeOuterExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &src,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Real sign;
  // add coupling variables to active mesh, reflecting v2
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
#pragma simd
        for (int n=0; n<(src.GetDim4()); ++n) {
          sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
          if (src.GetDim4()==6) // pressure tensor
            sign = flip_across_pole_ptensor[n] ? -1.0 : 1.0;
          src(k,je-j+1,i,n) += sign * src(k,je+j,i,n);
        }
      }
    }
  }

  return;
}
