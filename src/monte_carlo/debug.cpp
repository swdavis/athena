//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file debug.cpp
//  \brief temporary functions for debugging and testing during initial development

// Athena++ headers
#include "debug.hpp"

void FinalPositionCartesian(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot,
                            Real &xf, Real &yf, Real &zf, Real &dl) {

  Real dlx, dly, dlz;
  Real xw=-1.,yw=-1.,zw=-1.;

  if (pmcb->mcb_bcs[INNER_X1] != MC_PERIODIC_BNDRY) {
    if(pphot->k[0] > 0.0) {
      dlx = (pco->x1f(pmcb->ie+1) - pphot->x[0]) / pphot->k[0];
    } else if(pphot->k[0] < 0.0) {
      dlx = (pco->x1f(pmcb->is) - pphot->x[0]) / pphot->k[0];
    } else {
      dlx = HUGE_NUMBER;
    }
  } else {
    dlx  = HUGE_NUMBER;
    xw = pco->x1f(pmcb->ie+1)-pco->x1f(pmcb->is);
  }

  if (pmcb->mcb_bcs[INNER_X2] != MC_PERIODIC_BNDRY) {
    if(pphot->k[1] > 0.0) {
      dly = (pco->x2f(pmcb->je+1)  - pphot->x[1]) / pphot->k[1];
    } else if(pphot->k[1] < 0.0) {
      dly = (pco->x2f(pmcb->js) - pphot->x[1]) / pphot->k[1];
    } else {
      dly = HUGE_NUMBER;
    }
  } else {
    dly =  HUGE_NUMBER;
    yw = pco->x2f(pmcb->je+1)-pco->x2f(pmcb->js);
  }

  if (pmcb->mcb_bcs[INNER_X3] != MC_PERIODIC_BNDRY) {
    if(pphot->k[2] > 0.0) {
      dlz = (pco->x3f(pmcb->ke+1) - pphot->x[2])/pphot->k[2];
    } else if(pphot->k[2] < 0.0) {
      dlz = (pco->x3f(pmcb->ks) - pphot->x[2])/pphot->k[2];
    } else {
      dlz = HUGE_NUMBER;
    }
  } else {
    dlz = HUGE_NUMBER;
    zw = pco->x3f(pmcb->ke+1)-pco->x3f(pmcb->ks);
  }

  dl = dlx;
  if (dly < dl)
    dl = dly;
  if (dlz < dl)
    dl = dlz;
  
  // Compute final position
  xf = pphot->x[0] + pphot->k[0] * dl;
  yf = pphot->x[1] + pphot->k[1] * dl;
  zf = pphot->x[2] + pphot->k[2] * dl;

  // Account for periodic boundary condtions
  if (xw > 0.) {
    xf = fmod(xf,xw);
    if (xf >  pco->x1f(pmcb->ie+1)) xf -= xw;
    if (xf <  pco->x1f(pmcb->is)) xf += xw;
  }
  if (yw > 0.) {
    yf = fmod(yf,yw);
    if (yf >  pco->x2f(pmcb->je+1)) yf -= yw;
    if (yf <  pco->x2f(pmcb->js)) yf += yw;
  }
  if (zw > 0.) {
    zf = fmod(zf,zw);
    if (zf >  pco->x3f(pmcb->ke+1)) zf -= zw;
    if (zf <  pco->x3f(pmcb->ks)) zf += zw;
  }
  //std::cout << xf << ' ' << xw << ' ' << yf << ' ' << yw << std::endl;
}

void FinalPositionSphericalPolar(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot,
				 Real &rf, Real &thf, Real &phf, Real &dl) {

  Real cth = cos(pphot->x[1]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[2]);
  Real sph = sin(pphot->x[2]);

  // Convert to cartesian
  //Real kx = pphot->k[0] * sth*cph + pphot->k[1] * cth*cph - pphot->k[2] * sph;
  //Real ky = pphot->k[0] * sth*sph + pphot->k[1] * cth*sph + pphot->k[2] * cph;
  //Real kz = pphot->k[0] * cth - pphot->k[1] * sth;
  Real kx = pphot->kcart[0];
  Real ky = pphot->kcart[1];
  Real kz = pphot->kcart[2];

  // Outer boundary is r = rf -- find dlr to this boundary
  rf = pco->x1f(pmcb->ie+1);
  Real ndr0 = pphot->x[0] * (sth * (kx * cph + ky * sph) + kz * cth);
  Real det = 1.0 + (SQR(rf) - SQR(pphot->x[0])) / SQR(ndr0);
  Real dlr1 = ndr0 * (sqrt(det) - 1.0);
  Real dlr2 = -ndr0 * (sqrt(det) + 1.0);

  if (dlr1 > 0.0) {
    if (dlr2 > 0.0) {
      std::cout << "Warning: both roots positive in FinalPositionSphericalPolar: "
		<< dlr1 << " " << dlr2 << std::endl;
    } else {
      dl = dlr1;
    }
  } else if (dlr2 > 0.0) {
    dl = dlr2;
  }

  // Compute other boundary positions
  //theta
  Real zf = pphot->x[0] * cth + kz * dl;
  thf = acos(zf / rf);

  //phi
  Real xf = pphot->x[0] * sth * cph + kx * dl;
  Real yf = pphot->x[0] * sth * sph + ky * dl;
  phf = atan2(yf,xf);
  if (phf < 0.0)
    phf += 2.*PI;

  //std::cout << "x: " << pphot->x[0] << " " << pphot->x[1] << " " << pphot->x[2] << std::endl;
  //std::cout << "k: " << kx << " " << ky << " " << kz << std::endl;
  //std::cout << "dl: " << dl << " " << dlr1 << " " << dlr2 << " " << det << std::endl;
}
