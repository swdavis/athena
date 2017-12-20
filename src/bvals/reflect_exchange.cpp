//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect.cpp
//  \brief implementation of reflecting exchange BCs in each dimension

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"
#include "../hybrid/hybrid.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ReflectInnerExchangeX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x1 boundary

void ReflectInnerExchangeX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables to the active mesh, reflecting v1
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          if (n==(IM1) && mcoup.GetDim4()==NMCOUP) { // only apply this part to phybrid->mcoup
            mcoup(k,j,(is+i-1),IM1) -= mcoup(k,j,is-i,IM1);  // reflect 1-velocity
          } else {
            mcoup(k,j,(is+i-1),n) += mcoup(k,j,is-i,n);
          }
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectOuterExchangeX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x1 boundary

void ReflectOuterExchangeX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables to active mesh, reflecting v1
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          if (n==(IM1) && mcoup.GetDim4()==NMCOUP) { // only apply this part to phybrid->mcoup
            mcoup(k,j,(ie-i+1),IM1) -= mcoup(k,j,ie+i,IM1);  // reflect 1-velocity
          } else {
            mcoup(k,j,(ie-i+1),n) += mcoup(k,j,ie+i,n);
          }
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflecInnerExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                         const Real time, const Real dt,
//                         int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void ReflectInnerExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables to active mesh, reflecting v2
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          if (n==(IM2) && mcoup.GetDim4()==NMCOUP) { // only apply this part to phybrid->mcoup
            mcoup(k,js+j-1,i,IM2) -= mcoup(k,js-j,i,IM2);  // reflect 2-velocity
          } else {
            mcoup(k,js+j-1,i,n) += mcoup(k,js-j,i,n);
          }
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectOuterExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x2 boundary

void ReflectOuterExchangeX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables to active mesh, reflecting v2
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          if (n==(IM2) && mcoup.GetDim4()==NMCOUP) { // only apply this part to phybrid->mcoup
            mcoup(k,je-j+1,i,IM2) -= mcoup(k,je+j,i,IM2);  // reflect 2-velocity
          } else {
            mcoup(k,je-j+1,i,n) += mcoup(k,je+j,i,n);
          }
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectInnerExchangeX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, inner x3 boundary

void ReflectInnerExchangeX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // add coupling variables into ghost zones, reflecting v3
  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          if (n==(IM3) && mcoup.GetDim4()==NMCOUP) { // only apply this part to phybrid->mcoup
            mcoup(ks+k-1,j,i,IM3) -= mcoup(ks-k,j,i,IM3);  // reflect 3-velocity
          } else {
            mcoup(ks+k-1,j,i,n) += mcoup(ks-k,j,i,n);
          }
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReflectOuterExchangeX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
//                          const Real time, const Real dt,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions, outer x3 boundary

void ReflectOuterExchangeX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &mcoup,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones, reflecting v3
  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int n=0; n<(mcoup.GetDim4()); ++n) {
          if (n==(IM3) && mcoup.GetDim4()==NMCOUP) { // only apply this part to phybrid->mcoup
            mcoup(ke-k+1,j,i,IM3) -= mcoup(ke+k,j,i,IM3);  // reflect 3-velocity
          } else {
            mcoup(ke-k+1,j,i,n) += mcoup(ke+k,j,i,n);
          }
        }
      }
    }
  }

  return;
}
