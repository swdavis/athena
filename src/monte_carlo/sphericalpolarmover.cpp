//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sphericalpolarmover.cpp
//  \brief implementation for moving photons through sphericalpolar grid

// Athena++ headers
#include "photon.hpp"
#include "photonmover.hpp"
#include "../mesh/mesh.hpp"
#include "debug.hpp"
#define MAXITER 1000000
#define DEBUG

// Implementation of Sphericalpolar Photon mover

SphericalPolarMover::SphericalPolarMover(MonteCarloBlock *pmcb) 
  : PhotonMover(pmcb) {

}

SphericalPolarMover::~SphericalPolarMover() {

}

//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::Move(Photon *pphot)
//  \brief Moves photon along straight line specified number of mean free paths or until
//         photon leave monte carlo block

void SphericalPolarMover::Move(Photon *pphot) {

  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  Coordinates *pco = pmy_mcb->pmy_coord;

  // get number of mean free paths photon will travel
  Real TauRemaining = GetOpticalDepth(pran);

  // References for momentum vectors
  Real& kx = pphot->kcart[0];
  Real& ky = pphot->kcart[1];
  Real& kz = pphot->kcart[2];
  Real& kr  = pphot->k[0];
  Real& kth = pphot->k[1];
  Real& kph = pphot->k[2];

#ifdef DEBUG
  Real rf,thf,phf,dl0;
  FinalPositionSphericalPolar(pmcb,pco,pphot,rf,thf,phf,dl0);
#endif


  Real cth = cos(pphot->x[1]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[2]);
  Real sph = sin(pphot->x[2]);

  int iter = 0;
  // Move photon until requisite # of mean free paths or escape
  // MAXITER is present to account for (near) infinite trajectories in optically thin,
  // periodic domains.
  while( (TauRemaining > 0.) && (pphot->status == EVOLVING) && (iter < MAXITER)) {
    iter++;

    //Real& kx = pphot->kcart[0];
    //Real& ky = pphot->kcart[1];
    //Real& kz = pphot->kcart[2];
    // Update spherical polar momentum
    kr  = kx * sth * cph + ky * sth * sph + kz * cth;
    kth = kx * cth * cph + ky * cth * sph - kz * sth;
    kph = -kx * sph + ky * cph;
    // Compute cartesian positions
    Real r0 = pphot->x[0];
    Real x0 = r0 * sth * cph;
    Real y0 = r0 * sth * sph;
    Real z0 = r0 * cth;
    
    // Compute distance to all faces
    Real dl, dlr, dlt, dlp;
    bool ascend[3];

    // r face
    if (kr > 0.0) { // Can only intersect outer sphere
      Real ri = pco->x1f(pphot->i1+1) / r0;
      Real det = 1. + (SQR(ri) - 1.) / SQR(kr);
      dlr = kr * r0 * (sqrt(det) - 1.);
      ascend[0] = true;
    } else if (kr < 0.0) { // Can intersect either sphere
      Real ri = pco->x1f(pphot->i1) / r0;
      Real det = 1.0 + (SQR(ri)-1)/SQR(kr);
      if (det < 0.) {
        // ray does not intersect inner sphere so must intersect outer
        ri = pco->x1f(pphot->i1+1) / r0;
        det = 1. + (SQR(ri) - 1.) / SQR(kr);
        dlr = -kr * r0 * (sqrt(det) + 1.); //positive solution
        ascend[0] = true;
      } else {
        // ray intersects inner sphere first
        // one or two solutions -- pick shorter dlr
        dlr = r0 * kr * (sqrt(det) - 1.);
        ascend[0] = false;   
      }
    } else { // kr == 0
      Real ri = pco->x1f(pphot->i1+1) / r0;
      if (ri > 1.) {
        dlr = r0 * sqrt(SQR(ri) - 1.);
	ascend[0] = true;
      } else {
	std::cout << "Warning: kr == 0 and ri < r0, absorbing photon" << std::endl;
        pphot->status = DESTROYED;
        return;
      }
    }
    
    // theta face
    Real thi;
    if (kth > 0) {
      thi = pco->x2f(pphot->i2+1);
      ascend[1] = true;
    } else if (kth < 0) {
      thi = pco->x2f(pphot->i2);
      ascend[1] = false;
    } else {
      if (pphot->x[1] < 0.5*PI) {
        thi = pco->x2f(pphot->i2+1);
	ascend[1] = true;
      } else {
        thi = pco->x2f(pphot->i2);
	ascend[1] = false;
      }
    }
    Real cthi = cos(thi);

    // impact parameters
    Real b2 = SQR(r0) * (1. - SQR(kr));
    Real bz = r0 * (cth - kr * kz);
    Real th_ninf = acos(-kz);
    Real th_pinf = acos(kz);
    Real l_ext = (kz * b2 - r0 * kr * bz) / bz;
    Real cth_ext2 = SQR(kz) + SQR(bz)/b2;
    Real th_ext;
    if (bz < 0) {
      th_ext = acos(-sqrt(cth_ext2));
    } else {
      th_ext = acos(sqrt(cth_ext2));
    }

    // check for turning point
    if ((l_ext > 0) && ((thi-th_ninf) * (thi-th_ext)) > 0.) {
      dlt = HUGE_NUMBER;
    } else if ((l_ext <= 0) && ((thi-th_ext) * (thi-th_pinf)) > 0.) {
      dlt = HUGE_NUMBER;
    } else { //No turning point
      Real det = b2*SQR(cthi)*(cth_ext2-SQR(cthi)); //det should be positive
      Real lm = (-bz*kz - sqrt(det)) / (SQR(kz)-SQR(cthi)) - r0*kr;
      Real lp = (-bz*kz + sqrt(det)) / (SQR(kz)-SQR(cthi)) - r0*kr;
      if (lm <= 0) {
        dlt = lp;
      } else if (lp <= 0) {
        dlt = lm;
      } else {
        if (lp >= lm) 
          dlt = lm;
        else
          dlt = lp;
      }
    }

    // phi face
    Real phii;
    if (kph > 0) {
      phii = pco->x3f(pphot->i3+1);
      ascend[2] = true;
    } else if (kph < 0.0) {
      phii = pco->x3f(pphot->i3);
      ascend[2] = false;
    }
    
    if (kph != 0.) {
      Real tphii = tan(phii);
      if (fabs(kx * tphii - ky) < TINY_NUMBER) {
        dlp = HUGE_NUMBER;
      } else {
        dlp = r0 * sth * (sph-tphii*cth) / (kx*tphii-ky);
        if (dlp < 0.) dlp = HUGE_NUMBER;
      }
    }


    if ((dlr <= dlt) && (dlr <= dlp)) {
      dl = dlr;
    } else if ((dlp <= dlt) && (dlp <= dlr)) {
      dl = dlp;
    } else if ((dlt <= dlr) && (dlt <= dlp)) {
      dl = dlt;
    }

    //checking turning point for theta
    if ((dl > l_ext) && (l_ext > 0.)) {
      if (kth < 0.) {
        cthi = cos(pco->x2f(pphot->i2+1));
      } else if (kth > 0.) {
        cthi = cos(pco->x2f(pphot->i2));
      }
      b2 = SQR(r0) * (1. - SQR(kr));
      bz = r0 * (cth - kr * kz);
      Real det = b2 * SQR(cthi) * (cth_ext2 - SQR(cthi)); //det should be positive
      Real lm = (-bz*kz - sqrt(det)) / (SQR(kz) - SQR(cthi)) - kr*r0;
      Real lp = (-bz*kz + sqrt(det)) / (SQR(kz) - SQR(cthi)) - kr*r0;    
      if (lm <= 0) {
        if (lp > 0)
          dlt = lp;
        else
          dlt = HUGE_NUMBER;
      } else {
        if (lp > 0)
	  dlt = (lm > lp) ? lm : lp;
        else
          dlt = lm;
      }
      if (dlt <= dl) {
	ascend[1] = !(ascend[1]);
      }
    }

    int face;
    NextFace(dlr,dlt,dlp,face,dl);

    Real chi = pphot->sct_coef + pphot->abs_coef;
    chi = (chi > TINY_NUMBER) ? chi : TINY_NUMBER;
    if (dl > TauRemaining / chi) { // Photon remains in zone
      dl = TauRemaining/chi;
      // Update moments
      if (pmcb->moments_flag)
	pmcb->UpdateMoments(pphot,dl);
      // Update postions
      pphot->x[0] = sqrt(SQR(r0) + 2. * dl * pphot->k[0] * r0 + SQR(dl));
      pphot->x[1] = acos((z0 + kz * dl) / pphot->x[0]);
      pphot->x[2] = atan2(y0 + ky * dl,x0 + kx * dl);
      if (pphot->x[2] < 0.)
	pphot->x[2] += 2.*PI;
      return;
    } else { // Photon moves to next zone and reduce TauRemaining
      // Update moments
      if (pmcb->moments_flag)
	pmcb->UpdateMoments(pphot,dl);
      // Update postions
      pphot->x[0] = sqrt(SQR(r0) + 2. * dl * pphot->k[0] * r0 + SQR(dl));
      pphot->x[1] = acos((z0 + kz * dl) / pphot->x[0]);
      pphot->x[2] = atan2(y0 + ky * dl,x0 + kx * dl);
      if (pphot->x[2] < 0.)
	pphot->x[2] += 2.*PI;
    
      TauRemaining -= chi * dl;
      // move photon to next zone and pdate angular positions
      MovePhotonToNextZone(pphot,pco,pmcb,face,ascend);
      cth = cos(pphot->x[1]);
      sth = sqrt(1. - SQR(cth));
      cph = cos(pphot->x[2]);
      sph = sin(pphot->x[2]);
      
    }
  }
  if (iter >= MAXITER) {
    std::cout << "Warning: iter exceeded ITERMAX in photon mover." << std::endl;
    pphot->status = DESTROYED;
  }
#ifdef DEBUG
  /*Real delta = sqrt(SQR(xf-pphot->x[0])+SQR(yf-pphot->x[1])+SQR(zf-pphot->x[2]));
  Real dmax = 1.e-6*(pco->x3f(pmcb->ke+1)-pco->x3f(pmcb->ks));
  if ((delta > dmax)&&(iter < MAXITER)) {
    std::cout << "-----------------------" << std::endl;
    std::cout << delta <<  ' ' << iter << std::endl;
    std::cout << "k: " << pphot->k[0] << ' ' << pphot->k[1] << ' ' << pphot->k[2] << std::endl;
    std::cout << "xi: " << xi << ' ' << yi << ' ' << zi << std::endl; */
    std::cout << "xf: " << rf << ' ' << thf << ' ' << phf << ' ' << dl0 << std::endl;
    std::cout << "xp: " << pphot->x[0] << ' ' <<  pphot->x[1] << ' ' <<  pphot->x[2] << std::endl;
    //}
#endif

}

//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::CartesianToCurvalinear(Photon *pphot)
//  \brief convert k vector from cartesian to curvalinear

void SphericalPolarMover::CartesianToCurvalinear(Photon *pphot) {

  Real cth = cos(pphot->x[1]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[2]);
  Real sph = sin(pphot->x[2]);
  // Compute spherical-polar
  pphot->k[0] = pphot->kcart[0]*sth*cph + pphot->kcart[1]*sth*sph + pphot->kcart[2]*cth;
  pphot->k[1] = pphot->kcart[0]*cth*cph + pphot->kcart[1]*cth*sph - pphot->kcart[2]*sth;
  pphot->k[2] = -pphot->kcart[0]*sph + pphot->kcart[1]*cph;
  
}


//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::CurvalinearToCartesian(Photon *pphot)
//  \brief convert k vector from curvalinear to cartesian

void SphericalPolarMover::CurvalinearToCartesian(Photon *pphot) {

  Real cth = cos(pphot->x[1]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[2]);
  Real sph = sin(pphot->x[2]);
  // Compute cartesian
  pphot->kcart[0] = pphot->k[0]*sth*cph + pphot->k[1]*cth*cph - pphot->k[2]*sph;
  pphot->kcart[1] = pphot->k[0]*sth*sph + pphot->k[1]*cth*sph + pphot->k[2]*cph;
  pphot->kcart[2] = pphot->k[0]*cph - pphot->k[1]*sph;
  
}

