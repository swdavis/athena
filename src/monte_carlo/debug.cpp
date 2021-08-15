//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file debug.cpp
//! \brief temporary functions for debugging and testing during initial development

// Athena++ headers
#include "debug.hpp"

//----------------------------------------------------------------------------------------
//! \fn void FinalPositionCartesian(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot,
//!                                 Real &xf, Real &yf, Real &zf, Real &dl)
//! \brief set absorption opacity flag

void FinalPositionCartesian(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot,
                            Real &xf, Real &yf, Real &zf, Real &dl) {

  Real dlx, dly, dlz;
  Real xw=-1.,yw=-1.,zw=-1.;

  if (pmcb->mcb_bcs[BoundaryFace::inner_x1] != MC_PERIODIC_BNDRY) {
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

  if (pmcb->mcb_bcs[BoundaryFace::inner_x2] != MC_PERIODIC_BNDRY) {
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

  if (pmcb->mcb_bcs[BoundaryFace::inner_x3] != MC_PERIODIC_BNDRY) {
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
