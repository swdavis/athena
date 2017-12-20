//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outflow_exchange.cpp
//  \brief implementation of outflow exchange BCs in each dimension

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerExchangeX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, inner x1 boundary

void OutflowInnerExchangeX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables to active mesh
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
#pragma simd
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          mcoup(k,j,is,n) += mcoup(k,j,is-i,n);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterExchangeX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                         Real time, Real dt,
//                         int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x1 boundary

void OutflowOuterExchangeX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables to active mesh
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
#pragma simd
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          mcoup(k,j,ie,n) += mcoup(k,j,ie+i,n);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, inner x2 boundary

void OutflowInnerExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables into ghost zones
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
#pragma simd
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          mcoup(k,js,i,n) += mcoup(k,js-j,i,n);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x2 boundary

void OutflowOuterExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables to active mesh
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
#pragma simd
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          mcoup(k,je,i,n) += mcoup(k,je+j,i,n);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerExchangeX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, inner x3 boundary

void OutflowInnerExchangeX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables to active mesh
  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
#pragma simd
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          mcoup(ks,j,i,n) += mcoup(ks-k,j,i,n);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterExchangeX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          Real time, Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x3 boundary

void OutflowOuterExchangeX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables to active mesh
  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
#pragma simd
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          mcoup(ke,j,i,n) += mcoup(ke+k,j,i,n);
        }
      }
    }
  }

  return;
}
