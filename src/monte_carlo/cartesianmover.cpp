//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cartesianmover.cpp
//  \brief implementation for moving photons through cartesian grid

// Athena++ headers
#include "photon.hpp"
#include "photonmover.hpp"
#include "../mesh/mesh.hpp"
#include "debug.hpp"
#define MAXITER 1000000

// Implementation of Cartesian Photon mover

CartesianMover::CartesianMover(MonteCarloBlock *pmcb) 
  : PhotonMover(pmcb) {

}

CartesianMover::~CartesianMover() {

}

//----------------------------------------------------------------------------------------
//! \fn void CartesianMover::Move(Photon *pphot)
//  \brief Moves photon along straight line specified number of mean free paths or until
//         photon leave monte carlo block

void CartesianMover::Move(Photon *pphot) {

  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  Coordinates *pco = pmy_mcb->pmy_coord;

  // get number of mean free paths photon will travel
  Real TauRemaining = GetOpticalDepth(pran);

  Real xf,yf,zf,dl0;
  Real xi,yi,zi;
  xi = pphot->x[0]; yi = pphot->x[1]; zi = pphot->x[2];
  FinalPositionCartesian(pmcb,pco,pphot,xf,yf,zf,dl0);

  int iter = 0;
  // calculate distances to nearest faces
  while( (TauRemaining > 0.) && (pphot->status == EVOLVING) && (iter < MAXITER)) {
    iter++;
    // Compute distance to all faces
    Real dl, dlx, dly, dlz;
    bool ascend[3];
    if(pphot->k[0] > 0.0) {
      dlx = (pco->x1f(pphot->i1+1) - pphot->x[0]) / pphot->k[0];
      ascend[0] = true;
    } else if(pphot->k[0] < 0.0) {
      dlx = (pco->x1f(pphot->i1) - pphot->x[0]) / pphot->k[0];
      ascend[0] = false;
    } else {
      dlx = HUGE_NUMBER;
      ascend[0] = false;
    }

    if(pphot->k[1] > 0.0) {
      dly = (pco->x2f(pphot->i2+1)  - pphot->x[1]) / pphot->k[1];
      ascend[1] = true;
    } else if(pphot->k[1] < 0.0) {
      dly = (pco->x2f(pphot->i2) - pphot->x[1]) / pphot->k[1];
      ascend[1] = false;
    } else {
      dly = HUGE_NUMBER;
      ascend[1] = false;
    }
    
    if(pphot->k[2] > 0.0) {
      dlz = (pco->x3f(pphot->i3+1) - pphot->x[2])/pphot->k[2];
      ascend[2] = true;
    } else if(pphot->k[2] < 0.0) {
      dlz = (pco->x3f(pphot->i3) - pphot->x[2])/pphot->k[2];
      ascend[2] = false;
    } else {
      dlz = HUGE_NUMBER;
      ascend[2] = false;
    }

    int face;
    NextFace(dlx,dly,dlz,face,dl);

    Real chi = pphot->sct_coef + pphot->abs_coef;
    chi = (chi > TINY_NUMBER) ? chi : TINY_NUMBER;
    if (dl > TauRemaining / chi) { // Photon remains in zone
      UpdatePhotonPositionInZone(pphot,TauRemaining/chi);
      return;
    } else { // Photon moves to next zone and reduce TauRemaining
      MovePhotonToNextZone(pphot,pco,pmcb,dl,face,ascend);
      TauRemaining -= chi * dl;
    }
  }
  if (iter >= MAXITER) {
    std::cout << "Warning: iter exceeded ITERMAX in photon mover." << std::endl;
    pphot->status = DESTROYED;
  }
  Real delta = sqrt(SQR(xf-pphot->x[0])+SQR(yf-pphot->x[1])+SQR(zf-pphot->x[2]));
  if ((delta > 10.)&&(iter < MAXITER)) {
    std::cout << "-----------------------" << std::endl;
    std::cout << delta <<  ' ' << iter << std::endl;
    std::cout << "k: " << pphot->k[0] << ' ' << pphot->k[1] << ' ' << pphot->k[2] << std::endl;
    std::cout << "xi: " << xi << ' ' << yi << ' ' << zi << std::endl;
    std::cout << "xf: " << xf << ' ' << yf << ' ' << zf << ' ' << dl0 << std::endl;
    std::cout << "xp: " << pphot->x[0] << ' ' <<  pphot->x[1] << ' ' <<  pphot->x[2] << std::endl;
  }
}


//----------------------------------------------------------------------------------------
//! \fn void CartesianMover::UpdatePhotonPositionInZone(Photon *pphot)
//  \brief updates photon position without changing zones

void CartesianMover::UpdatePhotonPositionInZone(Photon *pphot, Real dl) {

  for (int i=0; i<3; ++i)
    pphot->x[i] += pphot->k[i] * dl;

}


//----------------------------------------------------------------------------------------
//! \fn void CartesianMover::MovePhotonToNextZone()
//  \brief updates photon position and changes zones

void CartesianMover::MovePhotonToNextZone(Photon *pphot, Coordinates *pco,
  MonteCarloBlock *pmcb, Real dl, int face, bool ascend[3]) {

  // Update positions
  for (int i=0; i<3; ++i)
    pphot->x[i] += pphot->k[i] * dl;

  // Update face(s) and adjust positions to lie exactly on boundary
  if ((face == 0) || (face == 3) || (face == 5) || (face == 6)) { //update x face
    if (ascend[0]) {
      pphot->i1++;
      if(pphot->i1 <= pmcb->ie)
        pphot->x[0] = pco->x1f(pphot->i1);
      else
        pmcb->pbval->BoundaryFunction_[OUTER_X1](pmcb,pco,pphot);
    } else {
      pphot->i1--;
      if(pphot->i1 >= pmcb->is)
        pphot->x[0] = pco->x1f(pphot->i1+1);
      else
        pmcb->pbval->BoundaryFunction_[INNER_X1](pmcb,pco,pphot);
    }
  }
  if ((face == 1) || (face == 3) || (face == 4) || (face == 6)) { //update y face
    if (ascend[1]) {
      pphot->i2++;
      if(pphot->i2 <= pmcb->je)
        pphot->x[1] = pco->x2f(pphot->i2);
      else
        pmcb->pbval->BoundaryFunction_[OUTER_X2](pmcb,pco,pphot);
    } else {
      pphot->i2--;
      if(pphot->i2 >= pmcb->js)
        pphot->x[1] = pco->x2f(pphot->i2+1);
      else
        pmcb->pbval->BoundaryFunction_[INNER_X2](pmcb,pco,pphot);
    } 
  }
  if ((face == 2) || (face == 4) || (face == 5) || (face == 6)) { //update z face
    if (ascend[2]) {
      pphot->i3++;
      if(pphot->i3 <= pmcb->ke)
        pphot->x[2] = pco->x3f(pphot->i3);
      else
        pmcb->pbval->BoundaryFunction_[OUTER_X3](pmcb,pco,pphot);
    } else {
      pphot->i3--;
      if(pphot->i3 >= pmcb->ks)
        pphot->x[2] = pco->x3f(pphot->i3+1);
      else
        pmcb->pbval->BoundaryFunction_[INNER_X3](pmcb,pco,pphot);
    } 
  }

  // Update opacities
  if (pphot->status == EVOLVING) {
    pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot);
    pphot->sct_coef = pmcb->ScatteringOpacity(pmcb,pphot);
  }
}
