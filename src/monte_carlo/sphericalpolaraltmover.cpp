//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sphericalpolaraltmover.cpp
//  \implementation for moving photons through sphericalpolar grid

// Athena++ headers
#include "photon.hpp"
#include "photonmover.hpp"
#include "../mesh/mesh.hpp"

// Temporary headers - holdover from GR generalmover
#define tolerance 1.e-5
#define max_iteration 2

//#define DEBUG_SM
//#define NBUFFER 50

//----------------------------------------------------------------------------------------
//! SphericalPolarMover class constructor, derived from PhotonMover base class

SphericalPolarAltMover::SphericalPolarAltMover(MonteCarloBlock *pmcb)
  : PhotonMover(pmcb) {

}

//----------------------------------------------------------------------------------------
//! destructor

SphericalPolarAltMover::~SphericalPolarAltMover() {

}

//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::Move(Photon *pphot, int ips, int ipe)
//! \brief Moves photon using cell-by-cell approach through cartesian grid, then converts to spherical polar

void SphericalPolarAltMover::Move(Photon *pphot, int ips, int ipe) {
  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  PhotonTrajectoryList *ptraj = pmy_mcb->ptraj;

  for (int ip=ips; ip<=ipe; ip++) 
    {
    // get number of mean free paths photon will travel
    Real tauremaining = GetOpticalDepth(pran);

    int count = 0;
    int iter = 0;
    int zone_counter = 0;
    Real chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);
    Real dl = tauremaining / chi;
    Real dmin = pmcb->pcoord->dmin[i3, i2, i1];

    // CM: Set step size to be min of dl and min distance to cell faces
    Real step = (dl < dmin) ? dl : dmin;

    while ( (pphot->statp[ip] == EVOLVING) && (tauremaining > TINY_NUMBER) &&
            (iter < checkmove)) {

      iter++;
      count++;

      SimpleStep(pphot,step,ip);
      tauremaining -= chi * step;

      // SWD: Clean up these checks
      // Check if photon changed zones
      if (UpdateZone(pphot,ip)) {
        UpdateOpacities(pphot,pmcb,ip);
        zone_counter++;
        chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);
      }
      if (pphot->statp[ip] == DESTROYED) {
        pphot->PrintPhoton(ip);
      }

      // Update moments
      if (pmcb->moments_flag) {
        pmcb->UpdateMoments(pphot,step,1.,ip);
      }

      if (pphot->IsNanPhoton(ip)) {
        pphot->PrintPhoton(ip);
        pphot->statp[ip] = DESTROYED;
      }

      // Perform any user work
      if (UserWorkInMove != NULL) UserWorkInMove(pmcb,pphot,this,ip);
      // SWD: put here for now, may need additional flag
      if (ptraj != NULL) ptraj->AddToTrajectory(pphot,ip);

    } // end of photon integration

    /*if (pphot->statp[ip] == ESCAPED) {
      pphot->ep[ip] *= pphot->k0p[ip];
      //pphot->PrintPhoton(ip);
      }*/

    if (iter >= checkmove) {
      std::cout << "Warning: iter exceeded " << checkmove << " in photon mover."
                << std::endl;
      pphot->PrintPhoton(ip);
      std::cout << "tau remaining, chi: " << tauremaining << " " << chi << std::endl;
      pphot->statp[ip] = DESTROYED;
    }
  }

}

// SWD: Deprecated and slated for removal
//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4])
//! \brief convert k vector from curvalinear to cartesian

void SphericalPolarAltMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4]) {

  Real cth = cos(pphot->x[IMC2]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[IMC3]);
  Real sph = sin(pphot->x[IMC3]);
  // Compute cartesian
  kcart[IMC1] = pphot->k[IMC1]*sth*cph + pphot->k[IMC2]*cth*cph - pphot->k[IMC3]*sph;
  kcart[IMC2] = pphot->k[IMC1]*sth*sph + pphot->k[IMC2]*cth*sph + pphot->k[IMC3]*cph;
  kcart[IMC3] = pphot->k[IMC1]*cth - pphot->k[IMC2]*sth;
  Real norm = sqrt(SQR(kcart[IMC1])+SQR(kcart[IMC2])+SQR(kcart[IMC3]));
  kcart[IMC1] /= norm;
  kcart[IMC2] /= norm;
  kcart[IMC3] /= norm;
}

//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarAltMover::UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb, int ip)
//! \brief update opacities after a photon has changed zones

void SphericalPolarAltMover::UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb, int ip) {

  pmy_mcb = pmcb;

    if (pphot->statp[ip] == EVOLVING) {
    // Opacities need to be calculated using comoving frame energy and then transformed
    // back to Eulerian frame when Lorentz Transformations are enabled.
    Real shift;
    int i1 = pphot->i1p[ip];
    int i2 = pphot->i2p[ip];
    int i3 = pphot->i3p[ip];
    if (pmcb->boosts) {
      // Shift photon energy to comoving frame
      shift = pmy_mcb->LorentzTransformFrequencyShift(pphot,ip);
      Real energy = pphot->ep[ip] * shift;
      // compute opacities in comoving frame
      pphot->acp[ip] = pmcb->AbsorptionOpacity(pmcb,i1,i2,i3,energy);
      pphot->scp[ip] = pmcb->ScatteringOpacity(pmcb,i1,i2,i3,energy);
      // Shift opaciteis to Eulerian frame
      pphot->acp[ip] *= shift;
      pphot->scp[ip] *= shift;
    } else {
      // No distinction between comovinng frame and eulerian frame
      pphot->acp[ip] = pmcb->AbsorptionOpacity(pmcb,i1,i2,i3,pphot->ep[ip]);
      pphot->scp[ip] = pmcb->ScatteringOpacity(pmcb,i1,i2,i3,pphot->ep[ip]);
    }

  }
}

//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarAltMover::SimpleStep(Photon *pphot, Real step, int ip)
//! \brief computes a step th

void SphericalPolarAltMover::SimpleStep(Photon *pphot, Real step, int ip) {

  Real x[NCOORD];

  // CM: not sure exactly what x[IMC#] is?
  // CM: add k direction vec * step to the current position
  x[IMC0] = pphot->x0p[ip] += (pphot->k0p[ip])*step;
  x[IMC1] = pphot->x1p[ip] += (pphot->k1p[ip])*step;
  x[IMC2] = pphot->x2p[ip] += (pphot->k2p[ip])*step;
  x[IMC3] = pphot->x3p[ip] += (pphot->k3p[ip])*step;

  // CM: Probably don't need second term for cartesian?
  k_n1[IMC0] = (pphot->k0p[ip]) + 0.5*(pphot->dk0p[ip]) * step;
  k_n1[IMC1] = (pphot->k1p[ip]) + 0.5*(pphot->dk1p[ip]) * step;
  k_n1[IMC2] = (pphot->k2p[ip]) + 0.5*(pphot->dk2p[ip]) * step;
  k_n1[IMC3] = (pphot->k3p[ip]) + 0.5*(pphot->dk3p[ip]) * step;

  // Set photon direction vectors based on current position
  pphot->k0p[ip] = k_n1[IMC0];
  pphot->k1p[ip] = k_n1[IMC1];
  pphot->k2p[ip] = k_n1[IMC2];
  pphot->k3p[ip] = k_n1[IMC3];
}
