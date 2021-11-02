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
//! \brief Moves photon using cell-by-cell approach through spherical polar grid

void SphericalPolarAltMover::Move(Photon *pphot, int ips, int ipe) {
  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  PhotonTrajectoryList *ptraj = pmy_mcb->ptraj;

  for (int ip=ips; ip<=ipe; ip++) {
    // get number of mean free paths photon will travel
    Real tauremaining = GetOpticalDepth(pran);

    Real step = StepSize(pphot,ip);
    int count = 0;
    int iter = 0;
    int zone_counter = 0;
    Real chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);

    while ( (pphot->statp[ip] == EVOLVING) && (tauremaining > TINY_NUMBER) &&
            (iter < checkmove)) {
      //printf("%d %g\n",iter,step);
      //pphot->PrintPhoton(ip);
      iter++;
      count++;

      if (tauremaining > chi * step) {
        VerletStep(pphot,step,ip);
        if (pmy_mcb->pmy_mc->polarized)
          PropogatePolarization(pphot,step,ip);
      } else {
        step = tauremaining / chi;
        VerletStep(pphot,step,ip);
        if (pmy_mcb->pmy_mc->polarized)
          PropogatePolarization(pphot,step,ip);
      }

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
      step = StepSize(pphot,ip);

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
  } // end loop over ip

}

// SWD: Deprecated and slated for removal
//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4])
//! \brief convert k vector from curvalinear to cartesian

void SphericalPolarMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4]) {

  Real cth = cos(pphot->x[IMC2]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[IMC3]);
  Real sph = sin(pphot->x[IMC3]);
  // Compute cartesian
  kcart[IMC1] = pphot->k[IMC1]*sth*cph + pphot->k[IMC2]*cth*cph - pphot->k[IMC3]*sph;
  kcart[IMC2] = pphot->k[IMC1]*sth*sph + pphot->k[IMC2]*cth*sph + pphot->k[IMC3]*cph;
  kcart[IMC3] = pphot->k[IMC1]*cth - pphot->k[IMC2]*sth;
}
