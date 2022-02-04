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
    // sample number of mean free paths photon will travel
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

    // shorthand for spherical polar direction vector components
    Real& kr  = pphot->k1p[ip];
    Real& kth = pphot->k2p[ip];
    Real& kph = pphot->k3p[ip];

    // create kx, ky, kz direction vector
    Real kx = kr*sth*cph + kth*cth*cph - kph*sph;
    Real ky = kr*sth*sph + kth*cth*sph + kph*cph;
    Real kz = kr*cth - kth*sth;

    // Store original values of r, theta, phi before step
    Real r0 = pphot->x1p[ip];
    Real th0 = pphot->x2p[ip];
    Real ph0 = pphot->x3p[ip];

    // create x, y, z position vector
    Real x0 = r0 * sth * cph;
    Real y0 = r0 * sth * sph;
    Real z0 = r0 * cth;

    // determine step size based on extinction coeff and distance to near face
    Real chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);
    Real chi0 = chi; // extinction coefficient in starting zone
    Real dl = tauremaining / chi; // distance left to travel if chi constant
    Real dmin = pmcb->pcoord->dmin(i3, i2, i1); // distance to nearest face

    // set step size to be min of dl and min distance to cell faces
    Real step = (dl < dmin) ? dl : dmin;

    // WHILE PHOTON IS TAKING STEPS
    // Continue taking steps until tauremaining is zero
    while ( (pphot->statp[ip] == EVOLVING) && (tauremaining > TINY_NUMBER) &&
            (iter < checkmove)) {

      iter++;
      count++;
      printf("\nSTEP #%d \n", iter);
      printf("================================\n");
      printf("i3, i2, i1=%d, %d, %d\n", i3, i2, i1);
      printf("r0, th0, ph0=%e, %e, %e\n", r0, th0, ph0);
      printf("x0, y0, z0=(%e, %e, %e)\n", x0, y0, z0);
      // Take a trial step
      // Move the photon in Cartesian coordinates
      x0 += step * kx;
      y0 += step * ky;
      z0 += step * kz;
      printf("Cartesian step size: (%e, %e, %e)\n", step * kx, step * ky,step * kz);
      // Update spherical polar position in pphot - expensive, do only once
      pphot->x1p[ip] = sqrt(SQR(x0) + SQR(y0) + SQR(z0));
      pphot->x2p[ip] = acos(z0/pphot->x1p[ip]);
      pphot->x3p[ip] = atan2(y0, x0);
      if (pphot->x3p[ip] < 0.)
        pphot->x3p[ip] += 2.*PI;
      printf("r, th, ph=%e, %e, %e\n", pphot->x1p[ip], pphot->x2p[ip],pphot->x3p[ip]);
      // IF PHOTON CHANGES ZONES
      // Check if photon changed zones
      if (UpdateZone(pphot,ip)) {
        printf("Logic: UpdateZone is true\n");
        printf("pphot->i#p[ip]=%d, %d, %d\n", pphot->i3p[ip], pphot->i2p[ip], pphot->i1p[ip]);
        printf("di=%d, %d, %d\n", abs(i1 - pphot->i1p[ip]), abs(i2 - pphot->i2p[ip]), abs(i3 - pphot->i3p[ip]));
        // IF PHOTON CHANGES BY ONE ZONE
        // Check that photon moved by only one zone in each direction
        // Moving through a corner or an edge is allowed by this logic
        if ((abs(i3 - pphot->i3p[ip]) <= 1) && (abs(i2 - pphot->i2p[ip]) <= 1)
            && (abs(i2 - pphot->i2p[ip]) <= 1)) {
          printf("Logic: Photon changed by one zone\n");
          // Use average opacity in starting and ending zone to calculate 
          // amount deducted from tauremaining

          // Update opacity and extinction coefficient in new zone
          UpdateOpacities(pphot,pmcb,ip);
          chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);

          // Find average opacity between original and new zones
          Real chi_avg = 0.5*(chi0 + chi);

          // Deduct distance travelled from remaining optical depth
          tauremaining -= chi_avg * step;
  
          // IF PHOTON MOVES PAST TARGET OPTICAL DEPTH
          if (tauremaining < 0.) {
            // Due to use of chi average, photon has moved too far
            // Take corrective step back using chi in new zone
            step = tauremaining / chi;
            tauremaining -= chi * step;

            // Take corrective step (step is negative)
            x0 += step * kx;
            y0 += step * ky;
            z0 += step * kz;

            // Update spherical polar position since photon loop is terminating
            pphot->x1p[ip] = sqrt(SQR(x0) + SQR(y0) + SQR(z0));
            pphot->x2p[ip] = acos(z0/pphot->x1p[ip]);
            pphot->x3p[ip] = atan2(y0, x0);
            if (pphot->x3p[ip] < 0.)
              pphot->x3p[ip] += 2.*PI;
          }

          // Update all zone variables now that the photon is in a new zone
          // Update zone indices
          i3 = pphot->i3p[ip];
          i2 = pphot->i2p[ip];
          i1 = pphot->i1p[ip];

          // Update values of r, theta, phi before next step
          // (in case next step has to be undone)
          r0 = pphot->x1p[ip];
          th0 = pphot->x2p[ip];
          ph0 = pphot->x3p[ip];

          // Update stepsize for new zone
          dl = tauremaining / chi;
          dmin = pmcb->pcoord->dmin(i3, i2, i1);
          step = (dl < dmin) ? dl : dmin;

          // Update chi now that this is our starting zone for the next step
          chi0 = chi;

        } else {
          printf("Logic: Photon changed by more than one zone\n");
          // PHOTON HAS MOVED THROUGH TOO MANY ZONES
          // Photon has moved by more than one index in a direction
          // Step back this move, halve the distance, and try again
          x0 -= step * kx;
          y0 -= step * ky;
          z0 -= step * kz;
          pphot->x1p[ip] = r0;
          pphot->x2p[ip] = th0;
          pphot->x3p[ip] = ph0;
          pphot->i1p[ip] = i1;
          pphot->i2p[ip] = i2;
          pphot->i3p[ip] = i3;
          step /= 2.0;
        }

      } else {
        printf("Logic: did not cross a zone boundary\n");
        // PHOTON HAS NOT CROSSED A ZONE BOUNDARY ON THIS STEP
        printf("chi=%e, step=%e, prod=%e\n", chi, step, chi*step);
        tauremaining -= chi * step;
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


      printf("tauremaining=%f\n", tauremaining);
      printf("dmin, dl=(%f, %e)\n", dmin, dl);
      printf("step=%f\n", step);
      printf("pphot->i#p[ip]=%d, %d, %d\n", pphot->i3p[ip], pphot->i2p[ip], pphot->i1p[ip]);
      printf("================================\n");


    } // end of photon integration

    // Update k vector
    cth = cos(pphot->x2p[ip]);
    sth = sin(pphot->x2p[ip]);
    cph = cos(pphot->x3p[ip]);
    sph = sin(pphot->x3p[ip]);
    kr  = kx * sth * cph + ky * sth * sph + kz * cth;
    kth = kx * cth * cph + ky * cth * sph - kz * sth;
    kph = -kx * sph + ky * cph;
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
