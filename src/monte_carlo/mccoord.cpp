//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//  \brief implementation of functions in class MCCoord

// Athena++ headers
#include "mccoord.hpp"


// constructor
MCCoord::MCCoord(Coordinates *pcoord, MonteCarloBlock *pmcb) {

  x1f.InitWithShallowCopy(pcoord->x1f);
  x2f.InitWithShallowCopy(pcoord->x2f);
  x3f.InitWithShallowCopy(pcoord->x3f);

  // Needed for black hole coordinates
  if (GENERAL_RELATIVITY) {
    bh_mass_ = pcoord->GetMass();
    bh_spin_ = pcoord->GetSpin();
  }

  // Allocate volume array
  int ncells1 = pmcb->nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmcb->nx2 > 1) ncells2 = pmcb->nx2 + 2*(NGHOST);
  if (pmcb->nx3 > 1) ncells3 = pmcb->nx3 + 2*(NGHOST);
  vol.NewAthenaArray(ncells3,ncells2,ncells1);
  // Initialize volume array
  for (int k=pmcb->ks; k<=pmcb->ke; ++k) {
    for (int j=pmcb->js; j<=pmcb->je; ++j) {
      for (int i=pmcb->is; i<=pmcb->ie; ++i) {
        vol(k,j,i) = pcoord->GetCellVolume(k,j,i);
      }}}
  if(pmcb->acceleration) {
    dmin.NewAthenaArray(ncells3,ncells2,ncells1);
    AthenaArray<Real> dw1,dw2,dw3;
    dw1.NewAthenaArray(ncells1);
    dw2.NewAthenaArray(ncells2);
    dw3.NewAthenaArray(ncells3);
    // Initialize dmin array
    for (int k=pmcb->ks; k<=pmcb->ke; ++k) {
      for (int j=pmcb->js; j<=pmcb->je; ++j) {
	pcoord->CenterWidth1(k,j,pmcb->is,pmcb->ie,dw1);
	pcoord->CenterWidth2(k,j,pmcb->is,pmcb->ie,dw2);
	pcoord->CenterWidth3(k,j,pmcb->is,pmcb->ie,dw3);
	for (int i=pmcb->is; i<=pmcb->ie; ++i) {
	  Real dmin0 = std::min(dw1(i),dw2(i));
	  dmin(k,j,i) = std::min(dmin0,dw3(i));
	  //printf("%d %d %d %g\n",k,j,i,dmin(k,j,i));
	}
      }}
    dw1.DeleteAthenaArray();
    dw2.DeleteAthenaArray();
    dw3.DeleteAthenaArray();
  }
}

// constructor
MCCoord::MCCoord(int ncells1, int ncells2, int ncells3, bool acc) {

  x1f.NewAthenaArray(ncells1+1);
  x2f.NewAthenaArray(ncells2+1);
  x3f.NewAthenaArray(ncells3+1);

  vol.NewAthenaArray(ncells3,ncells2,ncells1);
  acceleration = acc;
  if (acceleration) {
    dmin.NewAthenaArray(ncells3,ncells2,ncells1);
  }
}


// destructor
MCCoord::~MCCoord() {

  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();
  vol.DeleteAthenaArray();
  if (acceleration)
    dmin.DeleteAthenaArray();
}
