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

// function prototypes
Real DistanceToNearestFace(MCCoord *pco, Photon *pphot, int ip);

//----------------------------------------------------------------------------------------
//! CartesianMover class constructor, derived from PhotonMover base class

CartesianMover::CartesianMover(MonteCarloBlock *pmcb)
  : PhotonMover(pmcb) {

}

//----------------------------------------------------------------------------------------
//! destructor

CartesianMover::~CartesianMover() {

}

//----------------------------------------------------------------------------------------
//! \fn void CartesianMover::Move(Photon *pphot, int ips, int ipe)
//! \brief Moves photon using cell-by-cell approach through spherical polar grid

void CartesianMover::Move(Photon *pphot, int ips, int ipe) {

  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  MCCoord *pco = pmy_mcb->pcoord;

  for (int ip=ips; ip<=ipe; ip++) {

    // get number of mean free paths photon will travel
    Real tauremaining = GetOpticalDepth(pran);
    Real tau0 = tauremaining;

    Real& kx = pphot->k1p[ip];
    Real& ky = pphot->k2p[ip];
    Real& kz = pphot->k3p[ip];

    int iter = 0;
    Real c_cgs = 2.99792458e10;
    // checkmove is needed to account for (near) infinite trajectories that can occur
    // in optically thin, periodic domains.

    while( (tauremaining > 0.) && (pphot->statp[ip] == EVOLVING) && (iter < checkmove) &&
           (pphot->dtp[ip] > 0.) ) {
      iter++;

      // Compute distance to all faces
      Real dlx, dly, dlz;
      bool ascend[3];
      if(kx > 0.0) {
        dlx = (pco->x1f(pphot->i1p[ip]+1) - pphot->x1p[ip]) / kx;
        ascend[0] = true;
      } else if(kx < 0.0) {
        dlx = (pco->x1f(pphot->i1p[ip]) - pphot->x1p[ip]) / kx;
        ascend[0] = false;
      } else {
        dlx = HUGE_NUMBER;
        ascend[0] = false;
      }

      if(ky > 0.0) {
        dly = (pco->x2f(pphot->i2p[ip]+1)  - pphot->x2p[ip]) / ky;
        ascend[1] = true;
      } else if(ky < 0.0) {
        dly = (pco->x2f(pphot->i2p[ip]) - pphot->x2p[ip]) / ky;
        ascend[1] = false;
      } else {
        dly = HUGE_NUMBER;
        ascend[1] = false;
      }

      if(kz > 0.0) {
        dlz = (pco->x3f(pphot->i3p[ip]+1) - pphot->x3p[ip]) / kz;
        ascend[2] = true;
      } else if(kz < 0.0) {
        dlz = (pco->x3f(pphot->i3p[ip]) - pphot->x3p[ip]) / kz;
        ascend[2] = false;
      } else {
        dlz = HUGE_NUMBER;
        ascend[2] = false;
      }

      int face;
      NextFace(dlx,dly,dlz,face,dl);
      Real chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);

      if (dl > tauremaining / chi) { // Photon remains in zone
        bool accel_success = false;
        if (acceleration) {
          Real dist;
          dist = DistanceToNearestFace(pco,pphot,ip);

          // Try/perform MRW acceleration if optical depth is large enough
          if (pmcb->coherent_scattering) {
            Real tauacc = 10.;
            if ((pphot->acp[ip]+pphot->scp[ip]) * dist > tauacc)
              accel_success = MRWAcceleration(pphot,pran,dist,tauacc,ip);
          } else {
            Real tauacc = 10.;
            if (pmcb->planck_inv_opacity(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip])
                * dist > tauacc)
              accel_success = MRWAcceleration(pphot,pran,dist,tauacc,ip);
          }
        }

        // Perform standard displacement if acceleration not atempted or unsuccsessful
        if (!accel_success) {
          if (pphot->statp[ip] != EVOLVING)
            break;
          // compute distance remaining in zone
          dl = tauremaining/chi;
          pphot->dtp[ip] -= dl/c_cgs; // set with k0p instead
          // Account for absorption (if needed) and update moments
          Real etaua = ExpTauAbsorption(pphot->acp[ip],dl);

          if (pmcb->call_moments) {
            pmcb->UpdateMoments(pphot,dl,etaua,ip);
          }
          pphot->wp[ip] *= etaua;
          // update position
          pphot->x0p[ip] += pphot->k0p[ip] * dl;
          pphot->x1p[ip] += pphot->k1p[ip] * dl;
          pphot->x2p[ip] += pphot->k2p[ip] * dl;
          pphot->x3p[ip] += pphot->k3p[ip] * dl;
        }
        // Perform any user work
        if (UserWorkInMove != NULL) UserWorkInMove(pmcb,pphot,this,ip);
        break;

      } else { // Photon moves to next zone and reduce tauremaining
        // Account for absorption (if needed) and update moments
        Real etaua = ExpTauAbsorption(pphot->acp[ip],dl);

        if (pmcb->call_moments) {
          pmcb->UpdateMoments(pphot,dl,etaua,ip);
        }
        pphot->wp[ip] *= etaua;
        // update position
        pphot->x0p[ip] += pphot->k0p[ip] * dl;
        pphot->x1p[ip] += pphot->k1p[ip] * dl;
        pphot->x2p[ip] += pphot->k2p[ip] * dl;
        pphot->x3p[ip] += pphot->k3p[ip] * dl;

        tauremaining -= chi * dl;
        pphot->dtp[ip] -= dl/c_cgs;

        // Perform any user work
        if (UserWorkInMove != NULL) UserWorkInMove(pmcb,pphot,this,ip);
        MovePhotonToNextZone(pphot,pco,pmcb,face,ascend,ip);
      }
    }

    if (iter >= checkmove) {
      std::cout << "Warning: iter exceeded " << checkmove << " in photon mover."
                << std::endl;
      std::cout << "tau: " << tau0 << " " << tauremaining << std::endl;
      pphot->PrintPhoton(ip);
      pphot->statp[ip] = DESTROYED;
    }

  } // loop over photons

}

//----------------------------------------------------------------------------------------
//! \fn Real DistanceToNearestFace(MCCoord *pco, Photon *pphot, int ip)
//! \brief Computes distance to nearest face along current trajectory

Real DistanceToNearestFace(MCCoord *pco, Photon *pphot, int ip) {

  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  Real dx1p = pco->x1f(i1+1) - pphot->x1p[ip];
  Real dx1m = pphot->x1p[ip] - pco->x1f(i1);
  Real dx2p = pco->x2f(i2+1) - pphot->x2p[ip];
  Real dx2m = pphot->x2p[ip] - pco->x2f(i2);
  Real dx3p = pco->x3f(i3+1) - pphot->x3p[ip];
  Real dx3m = pphot->x3p[ip] - pco->x3f(i3);
  dx1p = (dx1p < dx1m) ? dx1p : dx1m;
  dx2p = (dx2p < dx2m) ? dx2p : dx2m;
  dx3p = (dx3p < dx3m) ? dx3p : dx3m;
  Real dist = (dx1p < dx2p) ? dx1p : dx2p;
  dist = (dist < dx3p) ? dist : dx3p;
  dist = (dist > 0.) ? dist : 0.;
  return dist;
}
