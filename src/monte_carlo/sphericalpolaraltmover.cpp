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
    int i3 = pphot->i3p[ip]; // Not a deep copy --- will not update with pphot
    int i2 = pphot->i2p[ip];
    int i1 = pphot->i1p[ip];

    // calculate sines and cosines of theta and phi angles
    Real cth = cos(pphot->x2p[ip]);
    Real sth = sin(pphot->x2p[ip]);
    Real cph = cos(pphot->x3p[ip]);
    Real sph = sin(pphot->x3p[ip]);

    // shorthand for spherical polar coordinates
    Real& kr  = pphot->k1p[ip];
    Real& kth = pphot->k2p[ip];
    Real& kph = pphot->k3p[ip];

    // create kx, ky, kz direction vector
    Real kx = kr*sth*cph + kth*cth*cph - kph*sph;
    Real ky = kr*sth*sph + kth*cth*sph + kph*cph;
    Real kz = kr*cth - kth*sth;

    // create x, y, z position vector
    Real r0 = pphot->x1p[ip];
    Real x0 = r0 * sth * cph;
    Real y0 = r0 * sth * sph;
    Real z0 = r0 * cth;

    // determine step size based on extinction coeff and distance to near face
    Real chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);
    Real dl = tauremaining / chi;
    Real dmin = pmcb->pcoord->dmin(i3, i2, i1);

    // set step size to be min of dl and min distance to cell faces
    Real step = (dl < dmin) ? dl : dmin;

    while ( (pphot->statp[ip] == EVOLVING) && (tauremaining > TINY_NUMBER) &&
            (iter < checkmove)) {

      iter++;
      count++;

      // Move the photon in Cartesian coordinates
      x0 += step * kx;
      y0 += step * ky;
      z0 += step * kz;

      // Deduct distance travelled from remaining optical depth
      tauremaining -= chi * step;

      // Update spherical polar position in pphot
      pphot->x1p[ip] = sqrt(SQR(x0) + SQR(y0) + SQR(z0));
      // SWD: original
      //pphot->x2p[ip] = acos(y0/sqrt(SQR(x0) + SQR(y0) + SQR(z0)));
      // SWD: new, no need to recompute r.  The compiler might be smart enough to not
      // recompute, but don't want to rely on it.
      pphot->x2p[ip] = acos(y0/pphot->x1p[ip]);
      pphot->x3p[ip] = atan2(y0, x0);
      // SWD atan2 will return values less than
      // 0 but code assumes phi: 0 - 2pi
      if (pphot->x3p[ip] < 0.)
        pphot->x3p[ip] += 2.*PI;

      // Check if photon changed zones
      if (UpdateZone(pphot,ip)) {
        // SWD: allowed_move not define. Only used once so better to remove
        //allowed_move = (abs(i3 - pphot->i3p[ip]) <= 1) && (abs(i2 - pphot->i2p[ip]) <= 1) && (abs(i2 - pphot->i2p[ip]) <= 1);
        //if (allowed_move) {

        if ((abs(i3 - pphot->i3p[ip]) <= 1) && (abs(i2 - pphot->i2p[ip]) <= 1)
            && (abs(i2 - pphot->i2p[ip]) <= 1)) {

          // Change of one or zero in each index
          // Update opacities and extinction coefficient
          UpdateOpacities(pphot,pmcb,ip);
          chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);

          // Update step size based on chi and face distances in the new zone
          dl = tauremaining / chi;
          dmin = pmcb->pcoord->dmin(pphot->i3p[ip], pphot->i2p[ip], pphot->i1p[ip]);
          step = (dl < dmin) ? dl : dmin;
        } else {
          // Step back this move, halve the distance, and try again
          x0 -= step * kx;
          y0 -= step * ky;
          z0 -= step * kz;
          tauremaining += chi * step;
          pphot->x1p[ip] = sqrt(SQR(x0) + SQR(y0) + SQR(z0));
          pphot->x2p[ip] = acos(y0/sqrt(SQR(x0) + SQR(y0) + SQR(z0)));
          pphot->x3p[ip] = atan2(y0, x0);
          UpdateZone(pphot,ip); // Updates pphot index
          step /= 2.0;
        }
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
