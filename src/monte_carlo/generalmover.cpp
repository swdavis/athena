//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file grmover.cpp
//! \brief implementation for moving photons via integration with metric and connection

// Athena++ headers
#include "montecarlo.hpp"
#include "photon.hpp"
#include "photonmover.hpp"
#include "../mesh/mesh.hpp"

// SWD: remove all of these
#define NCOORD 4

// GR headers
#define tolerance 1.e-5
#define max_iteration 2

//----------------------------------------------------------------------------------------
//! GeneralMover class constructor, derived from PhotonMover base class

GeneralMover::GeneralMover(MonteCarloBlock *pmcb)
  : PhotonMover(pmcb) {

  step_par = pmy_mcb->stepsize;

}

//----------------------------------------------------------------------------------------
//! destructor

GeneralMover::~GeneralMover() {

}

//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::Move(Photon *pphot, int ips, int ipe)
//! \brief Moves photon using geodesic integration

void GeneralMover::Move(Photon *pphot, int ips, int ipe) {

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
/*
      // Get distance from photon to closest cell face
      Real dmin;
      Real dw3, dw2, dw1;
      Real dx3p = std::min(fabs(pphot->x3p[ip] - pcoord->x3f(pphot->i3p[ip])), 
                           fabs(pphot->x3p[ip] - pcoord->x3f(pphot->i3p[ip] + 1)));
      Real dx2p = std::min(fabs(pphot->x2p[ip] - pcoord->x2f(pphot->i2p[ip])), 
                           fabs(pphot->x2p[ip] - pcoord->x2f(pphot->i2p[ip] + 1)));
      Real dx1p = std::min(fabs(pphot->x1p[ip] - pcoord->x1f(pphot->i1p[ip])), 
                           fabs(pphot->x1p[ip] - pcoord->x1f(pphot->i1p[ip] + 1)));
      dw3 = dx3p * pphot->x1p[ip] * sin(pphot->x2p[ip]);
      dw2 = dx2p * pphot->x1p[ip];
      dw1 = dx1p;
      Real dmin0 = std::min(dw1, dw2);
      dmin = std::min(dmin0, dw3);

      // Calculate accel threshold
      Real tauacc = 10.;
      if (chi * dmin > tauacc) {
//        printf("ACCELERATION: Accel triggered! \n");
        printf("ACCELERATION: taumin: %f    (%f %f %f)\n", chi*dmin, dw1*chi, dw2*chi, dw3*chi);
      } else {
        printf("------------: taumin: %f     (%f %f %f)\n", chi*dmin, dw1*chi, dw2*chi, dw3*chi);
      }     

      // Perform standard displacement if acceleration not attempted or unsuccessful
*/
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

// SWD: The conversion from Curvalinear  to Cartesian should be generalized or removed
//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4])
//! \brief convert k vector from curvalinear to cartesian

void GeneralMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4]) {

  Real r = pphot->x[IMC1];
  Real cth = cos(pphot->x[IMC2]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[IMC3]);
  Real sph = sin(pphot->x[IMC3]);

  Real nth = pphot->k[IMC2]*r;
  Real nph = pphot->k[IMC3]*r*sth;
  // Compute cartesian
  kcart[IMC1] = (pphot->k[IMC1]*sth*cph+nth*cth*cph-nph*sph);
  kcart[IMC2] = (pphot->k[IMC1]*sth*sph+nth*cth*sph+nph*cph);
  kcart[IMC3] = (pphot->k[IMC1]*cth-nth*sth);
  Real norm = sqrt(SQR(kcart[IMC1])+SQR(kcart[IMC2])+SQR(kcart[IMC3]));
  kcart[IMC1] /= norm;
  kcart[IMC2] /= norm;
  kcart[IMC3] /= norm;

}

//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb, int ip)
//! \brief update opacities after a photon has changed zones

void GeneralMover::UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb, int ip) {

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
      pphot->acp[ip] = pmcb->AbsorptionOpacity(pmcb,pphot,ip);
      pphot->scp[ip] = pmcb->ScatteringOpacity(pmcb,pphot,ip);
      // Shift opaciteis to Eulerian frame
      pphot->acp[ip] *= shift;
      pphot->scp[ip] *= shift;
    } else {
      // No distinction between comovinng frame and eulerian frame
      pphot->acp[ip] = pmcb->AbsorptionOpacity(pmcb,pphot,ip);
      pphot->scp[ip] = pmcb->ScatteringOpacity(pmcb,pphot,ip);
    }

  }
}

//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::VerletStep(Photon *pphot, Real step, int ip)
//! \brief performs a single verlet integration step

void GeneralMover::VerletStep(Photon *pphot, Real step, int ip) {

  Real k_n1[NCOORD],k_n1_copy[NCOORD];
  Real dk_n1[NCOORD];
  Real error;
  Real x[NCOORD], k0[NCOORD], dk0[NCOORD];

  // SWD: Need to think about how to do this better without invoking recurssion
  // SWD: Need to rename variables and clean this up with new scheme

  x[IMC0] = pphot->x0p[ip] += (pphot->k0p[ip])*step + 0.5*(pphot->dk0p[ip])*SQR(step);
  x[IMC1] = pphot->x1p[ip] += (pphot->k1p[ip])*step + 0.5*(pphot->dk1p[ip])*SQR(step);
  x[IMC2] = pphot->x2p[ip] += (pphot->k2p[ip])*step + 0.5*(pphot->dk2p[ip])*SQR(step);
  x[IMC3] = pphot->x3p[ip] += (pphot->k3p[ip])*step + 0.5*(pphot->dk3p[ip])*SQR(step);

  k_n1[IMC0] = (pphot->k0p[ip]) + 0.5*(pphot->dk0p[ip]) * step;
  k_n1[IMC1] = (pphot->k1p[ip]) + 0.5*(pphot->dk1p[ip]) * step;
  k_n1[IMC2] = (pphot->k2p[ip]) + 0.5*(pphot->dk2p[ip]) * step;
  k_n1[IMC3] = (pphot->k3p[ip]) + 0.5*(pphot->dk3p[ip]) * step;

  k0[IMC0] = pphot->k0p[ip];
  k0[IMC1] = pphot->k1p[ip];
  k0[IMC2] = pphot->k2p[ip];
  k0[IMC3] = pphot->k3p[ip];

  dk0[IMC0] = pphot->dk0p[ip];
  dk0[IMC1] = pphot->dk1p[ip];
  dk0[IMC2] = pphot->dk2p[ip];
  dk0[IMC3] = pphot->dk3p[ip];


  // Update gamma for current location
  pcoord->Connect(x, gamma);
  int n_iteration = 0;

  // SWD: not clear this while loops is accomplishing anything
  do {
    n_iteration += 1;
    error = 0.;
    for (int i=0;i<NCOORD;i++) {
      k_n1_copy[i] = k_n1[i];
    }

    for (int k=0;k<NCOORD;k++) {
      // off diagonal elements
      dk_n1[k] =
        -2. * (k_n1_copy[IMC0] *
               (gamma[k][IMC0][IMC1] * k_n1_copy[IMC1] +
                gamma[k][IMC0][IMC2] * k_n1_copy[IMC2] +
                gamma[k][IMC0][IMC3] * k_n1_copy[IMC3])
               +
               k_n1_copy[IMC1] * (gamma[k][IMC1][IMC2] * k_n1_copy[IMC2] +
                                  gamma[k][IMC1][IMC3] * k_n1_copy[IMC3]) +
               k_n1_copy[IMC2] * gamma[k][IMC2][IMC3] * k_n1_copy[IMC3]);
      // diagonal elements
      dk_n1[k] -=
        (gamma[k][IMC0][IMC0] * k_n1_copy[IMC0] * k_n1_copy[IMC0] +
         gamma[k][IMC1][IMC1] * k_n1_copy[IMC1] * k_n1_copy[IMC1] +
         gamma[k][IMC2][IMC2] * k_n1_copy[IMC2] * k_n1_copy[IMC2] +
         gamma[k][IMC3][IMC3] * k_n1_copy[IMC3] * k_n1_copy[IMC3]);

      k_n1[k] = k0[k] + 0.5 * (dk0[k] + dk_n1[k]) * step;

      error += fabs(k_n1_copy[k] - k_n1[k]) / (k_n1[k]);
    }
  } while ((error > tolerance) && (n_iteration < max_iteration));

  // SWD probably should not do this here
  // update photon energy due to evolving k_t (coordinate frame)
  //pphot->ep[ip] *= k_n1[IMC0]/(pphot->k[IMC0]);

  pphot->k0p[ip] = k_n1[IMC0];
  pphot->k1p[ip] = k_n1[IMC1];
  pphot->k2p[ip] = k_n1[IMC2];
  pphot->k3p[ip] = k_n1[IMC3];

  pphot->dk0p[ip] = dk_n1[IMC0];
  pphot->dk1p[ip] = dk_n1[IMC1];
  pphot->dk2p[ip] = dk_n1[IMC2];
  pphot->dk3p[ip] = dk_n1[IMC3];

}

//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::PropogatePolarization(Photon *pphot, Real step, int ip)
//! \brief propogates polarization tensor a single step

void GeneralMover::PropogatePolarization(Photon *pphot, Real step, int ip) {

  // SWD: Gamma does not need recomputing
  //Real gamma[NCOORD][NCOORD][NCOORD];
  // Store gamma in Coord to prevent recalculation
  //pcoord->Connect(pphot->x, gamma);


  std::complex<Real> ptcopy[4][4];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      ptcopy[i][j] = pphot->polten[i*4+j][ip];
    }
  }

  Real kp[4];
  kp[IMC0] = pphot->k0p[ip];
  kp[IMC1] = pphot->k1p[ip];
  kp[IMC2] = pphot->k2p[ip];
  kp[IMC3] = pphot->k3p[ip];

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          // eq. 16 of Moscibrodzka & Gammie in vacuum
          pphot->polten[i*4+j][ip] += -(gamma[i][k][l] * ptcopy[k][j] +
                                        gamma[j][k][l] * ptcopy[i][k]) *
                                       kp[l] * step;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real GeneralMover::StepSize(Photon *pphot, int ip)
//! \brief computes stepsize based on size of current zone

// SWD: Requires updates
// return the stepsize based on the current zone and k-vector
// this should be updated with every iteration since k continuously changes
Real GeneralMover::StepSize(Photon *pphot, int ip) {

  if (!pphot->pmy_mcb->varystep_flag) {
    return step_par; // keep step constant
  }

  Real small = 1.e-20;
  Real kx1 = (fabs(pphot->k1p[ip]) > small) ? fabs(pphot->k1p[ip]) : small;
  Real kx2 = (fabs(pphot->k2p[ip]) > small) ? fabs(pphot->k2p[ip]) : small;
  Real kx3 = (fabs(pphot->k3p[ip]) > small) ? fabs(pphot->k3p[ip]) : small;

  // SWD: May want to store as dx1, etc.
  Real stepx1 = ((pcoord->x1f(pphot->i1p[ip]+1)-pcoord->x1f(pphot->i1p[ip]))/kx1);
  Real stepx2 = ((pcoord->x2f(pphot->i2p[ip]+1)-pcoord->x2f(pphot->i2p[ip]))/kx2);
  Real stepx3 = ((pcoord->x3f(pphot->i3p[ip]+1)-pcoord->x3f(pphot->i3p[ip]))/kx3);

  Real step = (stepx1 < stepx2) ? stepx1 : stepx2;
  step = (step < stepx3) ? step : stepx3;

  return step * step_par;
}
