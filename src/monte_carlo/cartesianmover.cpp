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
#define MAXITER 10000
//#define DEBUG
#define TAU_TOL_COH 20.

Real DistanceToNearestFace(MCCoord *pco, Photon *pphot);

int ix[MAXITER],iy[MAXITER],iz[MAXITER];
// Implementation of Sphericalpolar Photon mover

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
  MCCoord *pco = pmy_mcb->pcoord;

  // get number of mean free paths photon will travel
  Real tauremaining = GetOpticalDepth(pran);
  Real tau0 = tauremaining;
#ifdef DEBUG
  Real xf,yf,zf,dl0;
  Real xi,yi,zi;
  xi = pphot->x[IMC1]; yi = pphot->x[IMC2]; zi = pphot->x[IMC3];
  FinalPositionCartesian(pmcb,pco,pphot,xf,yf,zf,dl0);
#endif

  CurvalinearToCartesian(pphot);
  Real& kx = pphot->k[IMC1];
  Real& ky = pphot->k[IMC2];
  Real& kz = pphot->k[IMC3];

  int iter = 0;
  // calculate distances to nearest faces
  while( (tauremaining > 0.) && (pphot->status == EVOLVING) && (iter < MAXITER)) {
    ix[iter] = pphot->i1;
    iy[iter] = pphot->i2;
    iz[iter] = pphot->i3;
    iter++;
    // Compute distance to all faces
    Real dlx, dly, dlz;
    bool ascend[3];
    if(kx > 0.0) {
      dlx = (pco->x1f(pphot->i1+1) - pphot->x[IMC1]) / kx;
      ascend[0] = true;
    } else if(kx < 0.0) {
      dlx = (pco->x1f(pphot->i1) - pphot->x[IMC1]) / kx;
      ascend[0] = false;
    } else {
      dlx = HUGE_NUMBER;
      ascend[0] = false;
    }

    if(ky > 0.0) {
      dly = (pco->x2f(pphot->i2+1)  - pphot->x[IMC2]) / ky;
      ascend[1] = true;
    } else if(ky < 0.0) {
      dly = (pco->x2f(pphot->i2) - pphot->x[IMC2]) / ky;
      ascend[1] = false;
    } else {
      dly = HUGE_NUMBER;
      ascend[1] = false;
    }
    
    if(kz > 0.0) {
      dlz = (pco->x3f(pphot->i3+1) - pphot->x[IMC3]) / kz;
      ascend[2] = true;
    } else if(kz < 0.0) {
      dlz = (pco->x3f(pphot->i3) - pphot->x[IMC3]) / kz;
      ascend[2] = false;
    } else {
      dlz = HUGE_NUMBER;
      ascend[2] = false;
    }

    int face;
    NextFace(dlx,dly,dlz,face,dl);

    Real chi = pphot->sct_coef + pphot->abs_coef;
    chi = (chi > TINY_NUMBER) ? chi : TINY_NUMBER;
    if (dl > tauremaining / chi) { // Photon remains in zone
      bool accel_success = false;
      if (acceleration) {
	Real dist;
	dist = DistanceToNearestFace(pco,pphot);
	/*dist = pco->dmin(pphot->i3,pphot->i2,pphot->i1);
	  dist = std::min(dl,dist);*/

	// Try/perform MRW acceleration if optical depth is large enough
	if (pmcb->coherent_scattering) {
	  Real tauacc = 10.;
	  if ((pphot->abs_coef+pphot->sct_coef) * dist > tauacc)
	    accel_success = MRWAcceleration(pphot,pran,dist,tauacc);
	} else {
	  Real tauacc = 10.;
	  if (pmcb->planck_inv_opacity(pphot->i3,pphot->i2,pphot->i1) * dist > tauacc)
	    accel_success = MRWAcceleration(pphot,pran,dist,tauacc);
	}
      }

      // Perform standard displacement if acceleration not atempted or unsuccsessful
      if (!accel_success) {
	if (pphot->status != EVOLVING)
	  return;
	// compute distance remaining in zone
        dl = tauremaining/chi;
	// Update moments
        if (pmcb->moments_flag) {
          pmcb->UpdateMoments(pphot,dl);
        }
        // update position
        for (int i=0; i<4; ++i)
          pphot->x[i] += pphot->k[i] * dl;
      }
      // Perform any user work
      if (UserWorkInMove != NULL) UserWorkInMove(pmcb,pphot,this);
      return;

    } else { // Photon moves to next zone and reduce tauremaining
      if (pmcb->moments_flag) {
	pmcb->UpdateMoments(pphot,dl);
      }      
      // update position
      for (int i=0; i<3; ++i)
	pphot->x[i] += pphot->k[i] * dl;
      tauremaining -= chi * dl;

      // Perform any user work
      if (UserWorkInMove != NULL) UserWorkInMove(pmcb,pphot,this);

      MovePhotonToNextZone(pphot,pco,pmcb,face,ascend);
    }
  }

  if (iter >= MAXITER) {
    std::cout << "Warning: iter exceeded MAXITER " << MAXITER << " in photon mover." 
	      << std::endl;
    std::cout << "tau: " << tau0 << " " << tauremaining << std::endl;
    pphot->PrintPhoton();
    pphot->status = DESTROYED;
  }

#ifdef DEBUG
  Real delta = sqrt(SQR(xf-pphot->x[IMC1])+SQR(yf-pphot->x[IMC2])+SQR(zf-pphot->x[IMC3]));
  Real dmax = 1.e-6*(pco->x3f(pmcb->ke+1)-pco->x3f(pmcb->ks));
  if ((delta > dmax)&&(iter < MAXITER)) {
    std::cout << "-----------------------" << std::endl;
    std::cout << delta <<  ' ' << iter << std::endl;
    std::cout << "k: " << pphot->k[IMC1] << ' ' << pphot->k[IMC2] << ' ' << pphot->k[IMC3] 
              << std::endl;
    std::cout << "xi: " << xi << ' ' << yi << ' ' << zi << std::endl;
    std::cout << "xf: " << xf << ' ' << yf << ' ' << zf << ' ' << dl0 << std::endl;
    std::cout << "xp: " << pphot->x[0] << ' ' <<  pphot->x[1] << ' ' <<  pphot->x[2] << std::endl;
  }
#endif
}

Real DistanceToNearestFace(MCCoord *pco, Photon *pphot) {

  Real dx1p = pco->x1f(pphot->i1+1) - pphot->x[IMC1];
  Real dx1m = pphot->x[IMC1] - pco->x1f(pphot->i1);
  Real dx2p = pco->x2f(pphot->i2+1) - pphot->x[IMC2];
  Real dx2m = pphot->x[IMC2] - pco->x2f(pphot->i2);
  Real dx3p = pco->x3f(pphot->i3+1) - pphot->x[IMC3];
  Real dx3m = pphot->x[IMC3] - pco->x3f(pphot->i3);
  dx1p = (dx1p < dx1m) ? dx1p : dx1m;
  dx2p = (dx2p < dx2m) ? dx2p : dx2m;
  dx3p = (dx3p < dx3m) ? dx3p : dx3m;
  Real dist = (dx1p < dx2p) ? dx1p : dx2p;
  dist = (dist < dx3p) ? dist : dx3p;
  dist = (dist > 0.) ? dist : 0.;
  return dist;
}
