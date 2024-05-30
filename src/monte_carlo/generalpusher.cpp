//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file grpusher.cpp
//! \brief implementation for moving photons via integration with metric and connection

// Athena++ headers
#include "montecarlo.hpp"
#include "photon.hpp"
#include "photonpusher.hpp"
#include "../mesh/mesh.hpp"

// SWD: remove all of these
#define NCOORD 4

// GR headers
#define tolerance 1.e-5
#define max_iteration 2

//----------------------------------------------------------------------------------------
//! GeneralPusher class constructor, derived from PhotonPusher base class

GeneralPusher::GeneralPusher(MonteCarloBlock *pmcb)
  : PhotonPusher(pmcb) {

  step_par = pmy_mcb->stepsize;

}

//----------------------------------------------------------------------------------------
//! destructor

GeneralPusher::~GeneralPusher() {

}

//----------------------------------------------------------------------------------------
//! \fn void GeneralPusher::Move(Photon *pphot, int ips, int ipe)
//! \brief Moves photon using geodesic integration

void GeneralPusher::Move(Photon *pphot, int ips, int ipe) {

  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  PhotonTrajectoryList *ptraj = pmy_mcb->ptraj;

  for (int ip=ips; ip<=ipe; ip++) {
    // get number of mean free paths photon will travel
    Real tauremaining = GetOpticalDepth(pran);
    Real step = StepSize(pphot,ip);
    Real path_length;
    Real k1, k2, k3;
    int iter = 0;
    Real c_cgs = 2.99792458e10;
    Real chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);
    while ( (pphot->statp[ip] == EVOLVING) && (tauremaining > TINY_NUMBER) &&
            (iter < checkmove) && (pphot->dtp[ip] > 0.) ) {
      iter++;

      bool accel_success = false;
      if ((acceleration) && (resonance)) {
        // Get distance from photon to closest cell face
        Real dl;
        Real dw3, dw2, dw1;
        Real dx3f = fabs(pcoord->x3f(pphot->i3p[ip]) - pcoord->x3f(pphot->i3p[ip] + 1));
        Real dx2f = fabs(pcoord->x2f(pphot->i2p[ip]) - pcoord->x2f(pphot->i2p[ip] + 1));
        Real dx1f = fabs(pcoord->x1f(pphot->i1p[ip]) - pcoord->x1f(pphot->i1p[ip] + 1));
        Real x1v = (pcoord->x1f(pphot->i1p[ip]) + pcoord->x1f(pphot->i1p[ip] + 1))/2.;
        Real x2v = (pcoord->x2f(pphot->i2p[ip]) + pcoord->x2f(pphot->i2p[ip] + 1))/2.;
        dw3 = dx3f * x1v * sin(x2v);
        dw2 = dx2f * x1v;
        Real dmin0 = std::min(dx1f, dw2);
        dl = std::min(dmin0, dw3); // Distance to nearest face

        Real tauacc = 1000.; //BCM: make this an input parameter
        // Try to perform MRW acceleration if optical depth is large enough
        if (dl*chi > tauacc) {
          MRWResonanceAcceleration(pphot,pran,dl,tauacc,path_length,k1,k2,k3,ip);
          accel_success = true;
        }
      }

      if (!accel_success) {// Acceleration not triggered - take standard step
        if (tauremaining > chi * step) { // Photon hasn't yet reached tauremaining
          VerletStep(pphot,step,ip);
          if (pmy_mcb->pmy_mc->polarized)
            PropogatePolarization(pphot,step,ip);
          tauremaining -= chi * step;
        } else { // Photon has reached end of tauremaining - step to make it 0
          step = tauremaining / chi;
          VerletStep(pphot,step,ip);
          if (pmy_mcb->pmy_mc->polarized)
            PropogatePolarization(pphot,step,ip);
          tauremaining = 0.;
        }
        pphot->dtp[ip] -= step/c_cgs;
        // Update moments
        if (pmcb->call_moments) {
          pmcb->UpdateMoments(pphot,step,ip);
        }
      } else {
        // Photon has been given a new position on sphere of radius dl
        // Set exit parameters and continue the loop over photons
        step = dl;
        tauremaining = 0.;
        if (pmcb->call_moments) {
          pmcb->UpdateMomentsAcceleration(pphot,step,path_length,k1,k2,k3,1.,ip);
        }
      }

      // Check if photon changed zones
      if (UpdateZone(pphot,ip)) {
        UpdateOpacities(pphot,pmcb,ip);
        chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);
      }

      if (pphot->IsNanPhoton(ip)) {
        pphot->statp[ip] = DESTROYED;
        pphot->PrintPhoton("Photon returned Nan in general pusher",ip);
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
      std::cout << "Warning: iter exceeded " << checkmove << " in photon pusher."
                << std::endl;
      pphot->PrintPhoton(ip);
      std::cout << "tau remaining, chi: " << tauremaining << " " << chi << std::endl;
      pphot->statp[ip] = DESTROYED;
    }
  } // end loop over ip
}

// SWD: The conversion from Curvalinear  to Cartesian should be generalized or removed
//----------------------------------------------------------------------------------------
//! \fn void GeneralPusher::CurvalinearToCartesian(Photon *pphot, Real kcart[4])
//! \brief convert k vector from curvalinear to cartesian

void GeneralPusher::CurvalinearToCartesian(Photon *pphot, Real kcart[4]) {

  /*
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
  */
}

//----------------------------------------------------------------------------------------
//! \fn void GeneralPusher::UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb, int ip)
//! \brief update opacities after a photon has changed zones

void GeneralPusher::UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb, int ip) {

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
//! \fn void GeneralPusher::VerletStep(Photon *pphot, Real step, int ip)
//! \brief performs a single verlet integration step

void GeneralPusher::VerletStep(Photon *pphot, Real step, int ip) {

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
//! \fn void GeneralPusher::PropogatePolarization(Photon *pphot, Real step, int ip)
//! \brief propogates polarization tensor a single step

void GeneralPusher::PropogatePolarization(Photon *pphot, Real step, int ip) {

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
//! \fn Real GeneralPusher::StepSize(Photon *pphot, int ip)
//! \brief computes stepsize based on size of current zone

// SWD: Requires updates
// return the stepsize based on the current zone and k-vector
// this should be updated with every iteration since k continuously changes
Real GeneralPusher::StepSize(Photon *pphot, int ip) {

  if (!pmy_mcb->varystep_flag) {
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
