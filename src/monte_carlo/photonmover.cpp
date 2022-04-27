//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photonmover.cpp
//! \brief implementation for photon moving functions

// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "photon.hpp"
#include "photonmover.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"

//----------------------------------------------------------------------------------------
//! PhotonMover base class constructor, built from  MonteCarloBlock

PhotonMover::PhotonMover(MonteCarloBlock *pmcb) {

  pmy_mcb = pmcb;
  pcoord = NULL;
  UserWorkInMove = pmcb->pmy_mc->UserWorkInMove;
  checkmove = pmcb->pmy_mc->checkmove;

  // MRW acceleration
  acceleration = pmcb->acceleration;
  boosts = pmcb->boosts;
  resonance = (pmcb->scattering_meth == SCATRES);
  compton = (!pmcb->coherent_scattering) && (!resonance);
  time_acc = pmcb->time_acc;

  if (acceleration) {
    InitializeMRWDist();
    if (compton)
      ReadComptonGreensFunction();
    if (time_acc) {
      ReadRadiusDistribution();
      ReadTimeDistribution();
    }
  }
}

//----------------------------------------------------------------------------------------
//! destructor

PhotonMover::~PhotonMover() {

  if (acceleration) {
    mrwprob.DeleteAthenaArray();
    mrwdev.DeleteAthenaArray();
    if (compton) {
      mrwxf.DeleteAthenaArray();
      mrwt.DeleteAthenaArray();
      mrwxi.DeleteAthenaArray();
      mrwp.DeleteAthenaArray();
    }
    if (time_acc) {
      mrwrt.DeleteAthenaArray();
      mrwrp.DeleteAthenaArray();
      mrwrr.DeleteAthenaArray();
    }
    mrwta.DeleteAthenaArray();
    mrwtp.DeleteAthenaArray();
    mrwtt.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonMover::Move(Photon *pphot, int ips, int ipe)
//  \brief base class move does nothing

void PhotonMover::Move(Photon *pphot, int ips, int ipe) {

}

//----------------------------------------------------------------------------------------
//! \fn bool PhotonMover::MRWAcceleration(Photon *pphot, MCRandom *pran, Real dist, 
//!                                       Real tauacc, int ip)
//! \brief Accelerate photon diffusion with modified random walk method

bool PhotonMover::MRWAcceleration(Photon *pphot, MCRandom *pran, Real dist, Real tauacc,
                                  int ip) {

  MonteCarloBlock *pmcb = pmy_mcb;
  bool accel_success = true;

  if (resonance) {
    Real r0 = dist;

    // Line constants
    Real melectron = 9.10938215e-28;
    Real charge = 4.80320427e-10;
    Real osc_strength = 0.4164;
    Real nu0 = 2.468e15;
    Real c = 2.99792458e10;
    Real h = 6.62607015e-27;
    Real kb = 1.380649e-16;
    Real mass = 1.660538782e-24;

    if (COORDINATE_SYSTEM == "spherical_polar") {

      // ********* FREQUENCY REDISTRIBUTION *********
      // Sample outgoing frequency from Dijkstra et al 2006 solution
      int &i1 = pphot->i1p[ip];
      int &i2 = pphot->i2p[ip];
      int &i3 = pphot->i3p[ip];

      // Cell properties
      Real tgas = pmcb->tgas(i3,i2,i1);
      Real rho = pmcb->rho(i3,i2,i1);

      // Derived parameters
      Real vth = sqrt( 2 * kb * tgas / mass);
      Real doppwidth = nu0 * vth / c;
      Real lorwidth = 6.265e8/(4.*PI);
      Real a = lorwidth / doppwidth;
      Real k = (rho/mass) * PI*charge*charge / (melectron*c) * osc_strength;
      Real tau0 = k * r0 / sqrt(PI) / doppwidth;
      //if (a*tau0 > 1.) {
      //  return false;
      //}

      // Sample sigma and convert to frequency
      Real x_s = (pphot->ep[ip] / h - nu0)/doppwidth;
      Real sigma_s = sqrt(2./3.) * PI/a * std::pow(x_s, 3.)/3.;
      Real samp_sigma = tau0 * 2./sqrt(PI) * std::atanh(2.*pran->uniform() - 1.) + sigma_s;
      Real x = std::cbrt(3. * sqrt(3./2.) * a / PI * samp_sigma);
      Real nu = doppwidth * x + nu0;
      pphot->ep[ip] = h * nu;

      // ********* POSITION *********
      // position packet on sphere of radius r0
      Real mu = 2.*pran->uniform()-1.0;

      // Local angles within the sphere of radius r0
      Real lsth = sqrt(1.0-mu*mu);
      Real lphi = 2.*PI*pran->uniform();

      // convert to cartesian
      // Global simulation angles based on photon position
      Real cth = cos(pphot->x2p[ip]);
      Real sth = sqrt(1. - SQR(cth));
      Real cph = cos(pphot->x3p[ip]);
      Real sph = sin(pphot->x3p[ip]);
      Real r = pphot->x1p[ip];
      Real x0 = r * sth * cph;
      Real y0 = r * sth * sph;
      Real z0 = r * cth;
      x0 += lsth*cos(lphi) * r0;
      y0 += lsth*sin(lphi) * r0;
      z0 += mu * r0;
      pphot->x1p[ip] = sqrt(SQR(x0)+SQR(y0)+SQR(z0));
      pphot->x2p[ip] = acos(z0 / pphot->x1p[ip]);
      pphot->x3p[ip] = atan2(y0,x0);
      if (pphot->x3p[ip] < 0.)
        pphot->x3p[ip] += 2.*PI;

      // ********* DIRECTION *********
      // Sample outgoing angles to local normal - zero ingoing flux, so must be outward
      Real sq3 = 2.*sqrt(3.);
      Real xi = pran->uniform();
      Real samp_cth = (2./sq3)*(sqrt(-sq3 * xi - 3.*xi + sq3 + 4.) - 1.);
      Real samp_sth = sqrt(1.0 - samp_cth*samp_cth);
      Real samp_phi = 2.*PI*pran->uniform();

      // Local cos and sin of phi - positions in sphere of radius r0
      Real lcph = cos(lphi);
      Real lsph = sin(lphi);

      // Cartesian direction vectors in local sphere
      // n = er * samp_cth + etheta * samp_sth * cos(samp_phi) + ephi * samp_sth * sin(samp_phi)
      Real nx = lsth * lcph * samp_cth + mu * lcph * samp_sth * cos(samp_phi) - lsph * samp_sth * sin(samp_phi);
      Real ny = lsth * lsph * samp_cth + mu * lsph * samp_sth * cos(samp_phi) + lcph * samp_sth * sin(samp_phi);
      Real nz = mu * samp_cth - lsth * samp_sth * cos(samp_phi);

      // Updated global coordinates after the move onto surf of sphere
      cth = cos(pphot->x2p[ip]);
      sth = sqrt(1. - SQR(cth));
      cph = cos(pphot->x3p[ip]);
      sph = sin(pphot->x3p[ip]);

      // Global sphpol direction vectors from local cartesian on the sphere
      pphot->k1p[ip] = nx * sth * cph + ny * sth * sph + nz * cth;
      pphot->k2p[ip] = nx * cth * cph + ny * cth * sph - nz * sth;
      pphot->k3p[ip] = -nx * sph + ny * cph;

    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in function [PhotonMover::MRWAcceleration]"
            <<std::endl<< "Specified coordinate system not implemented for resonance acceleration" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
    }


  } else {
    // draw from path length distribution and reduce weight accordingly
    //Real mrw = MRWDist(pran);
    //while (mrw <= 0.)
    //  mrw = MRWDist(pran);
    //Real delta = -log(mrw);
    Real delta = InterpPathTime(pphot->scp[ip]*dist,pran->uniform());
    //printf("delta: %g\n",delta);
    Real chi;
    if (!compton)
      chi = 3. * (pphot->acp[ip]+pphot->scp[ip]) / SQR(PI);
    else {
      //chi = 3. * pmcb->planck_inv_opacity(pphot->i3,pphot->i2,pphot->i1) / SQR(PI);
      chi = 3. * pphot->scp[ip] / SQR(PI);
    }
    Real ct,r0;
    Real beta[3], beta2, gamma, gonembdk;
    if (boosts) {
      // tranform relevant quanitites to comoving frame
      beta[0] = pmcb->vel(0,pphot->i3p[ip],pphot->i2,pphot->i1) / 2.99792458e10;
      beta[1] = pmcb->vel(1,pphot->i3p[ip],pphot->i2,pphot->i1) / 2.99792458e10;
      beta[2] = pmcb->vel(2,pphot->i3p[ip],pphot->i2,pphot->i1) / 2.99792458e10;
      beta2 = SQR(beta[0]) + SQR(beta[1]) + SQR(beta[2]);
      gamma = 1./sqrt(1.-beta2);
      Real bdk = (pphot->k[0] * beta[0] + pphot->k[1] * beta[1] + pphot->k[2] * beta[2]);
      gonembdk = gamma * (1. - bdk);
      chi *= (1.+4.*beta2)/gonembdk;
      pphot->abs_coef /= gonembdk;
      pphot->sct_coef /= gonembdk;
      pphot->energy *= gonembdk;
      // multiply beta by gamma
      Real betagamma = sqrt(beta2)*gamma;
      for(int i=0; i<3; ++i)
        beta[i] *= gamma;
      // Compute radius of sphere in comoving frame
      if (betagamma > 0.) {
        r0 = 0.5*(sqrt(1.+4.*chi*dist*delta*betagamma)-1.)/(delta*betagamma*chi);
        Real tau;
        if (!compton)
          tau = (pphot->abs_coef+pphot->sct_coef) * r0;
        else {
          //tau = pmcb->planck_inv_opacity(pphot->i3,pphot->i2,pphot->i1) *r0;
          tau = pphot->sct_coef*r0;
        }
        if (tau > tauacc) {
          ct = delta*SQR(r0)*chi;
        } else {
          accel_success = false;
        }
      } else {
        // zone has zero velocity so use method for static MRW
        for(int i=0; i<3; ++i)
          beta[i] = 0.;
        ct = delta*SQR(dist)*chi;
        //Real tau = pmcb->planck_inv_opacity(pphot->i3,pphot->i2,pphot->i1) * dist;
        //ct = -log(pran->uniform())*3./SQR(PI)*SQR(tau+2./3.);
        r0 = dist;
      }
    } else {
      // use method for static MRW
      for(int i=0; i<3; ++i)
        beta[i] = 0.;
      ct = delta*SQR(dist)*chi;
      /*if (!compton) {
        Real tau = (pphot->abs_coef+pphot->sct_coef) * dist;
        ct = -log(pran->uniform())*3./SQR(PI)*SQR(tau+2./3.) /
              (pphot->abs_coef+pphot->sct_coef);
      } else {
        Real tau =  pphot->sct_coef * dist;
        //printf("%g %g\n",tau,dist);
        ct = -log(pran->uniform())*3./SQR(PI)*SQR(tau+2./3.) /
              pmcb->planck_inv_opacity(pphot->i3,pphot->i2,pphot->i1);
        }*/
      r0 = dist;
    }
    if (accel_success) {
      Real tauabs;
      //pphot->path += ct;
      if (!compton) {
        tauabs = ct*pphot->abs_coef;
        pphot->weight *= exp(-tauabs);
      } else {

        //tauabs = ct*sqrt(pphot->abs_coef*pmcb->planck_opacity(pphot->i3,pphot->i2,pphot->i1));
        //tauabs = ct*pmcb->planck_opacity(pphot->i3,pphot->i2,pphot->i1);
        //tauabs = ct*pphot->abs_coef;
        Real opaci = pphot->abs_coef;
        Real c = 2.99792458e10;
        Real kb = 1.380649e-16;
        Real me = 9.1093897e-28;
        Real temp = pmcb->tgas(pphot->i3,pphot->i2,pphot->i1);
        Real xi = pphot->energy  / (kb *temp);
        Real ypar = pphot->sct_coef*ct*kb*temp/(me*c*c);
        Real xf = InterpComptonEnergy(xi,ypar,pran->uniform());

        pphot->energy = xf * kb * temp;
        //tauabs = ct*pmcb->planck_opacity(pphot->i3,pphot->i2,pphot->i1);
        //pphot->energy = PlanckDist(pmcb->tgas(pphot->i3,pphot->i2,pphot->i1),pran);
        pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot,ip);
        Real opacf = pphot->abs_coef;
        //Real opacf = std::max(pphot->abs_coef,pmcb->planck_opacity(pphot->i3,pphot->i2,
        //pphot->i1));
        //tauabs = ct*sqrt(opaci*pphot->abs_coef);
        //tauabs = ct*sqrt(opaci*(opacf));
        //tauabs = ct*(opaci-opacf)/(2.*log(xf/xi));
        Real ct0 = ypar/ct*log(xf/xi);
        Real ct1 = ct;
        //tauabs = ct0*(opaci-opacf)/(2.*log(xf/xi))+opacf*ct1;
        if (xi < 1.) {
          if (xf > 1.) {
            tauabs = ct1*(opaci-2./3.*opacf)/(2.*log(xf/xi));
          } else  {
            tauabs = ct1*(opaci-opacf)/(2.*log(xf/xi));
          }
        } else {
          if (xf > 1.) {
            tauabs = ct1*(opaci-opacf)/(3.*log(xf/xi));
          } else  {
            tauabs = ct1*(opaci-opacf)/(2.*log(xf/xi));
          }
        }
        //tauabs = opaci/8./ypar*ct;

        Real y = exp(4.*ypar);
        Real tau0 = tauabs;
        //tauabs = (-SQR(4-xi)+4.*y*(xi-4)*xi+SQR(y)*(16.+8.*xi+(4.*ypar-3.)*
        //         SQR(xi)))/32./SQR(y)*opaci/ypar*ct;
        //tauabs = (1.-1./SQR(y))/8.*opaci/ypar*ct;
        //tauabs = fabs(opaci-opacf)/8./ypar*ct;
        Real us = fabs(0.25*log(xf/xi));
        Real u0 = ct/ypar;
        //tauabs = fabs(opaci-opacf)/8.*u0;
        //printf("%g %g %g %g %g %g %g\n",ct,tauabs,tau0,opaci,opacf,xi,xf);
        pphot->weight *= exp(-tauabs);
      }

      // position packet on sphere of radius r0
      Real mu = 2.*pran->uniform()-1.0;
      Real stheta = sqrt(1.0-mu*mu);
      Real phi = 2.*PI*pran->uniform();
      if (COORDINATE_SYSTEM == "cartesian") {
        pphot->x[0] += stheta*cos(phi) * r0 + beta[0] * ct;
        pphot->x[1] += stheta*sin(phi) * r0 + beta[1] * ct;
        pphot->x[2] += mu * r0 + beta[2] * ct;
      } else if (COORDINATE_SYSTEM == "spherical_polar") {
        // convert to carteisan
        Real cth = cos(pphot->x[1]);
        Real sth = sqrt(1. - SQR(cth));
        Real cph = cos(pphot->x[2]);
        Real sph = sin(pphot->x[2]);
        Real r = pphot->x[0];
        Real x0 = r * sth * cph;
        Real y0 = r * sth * sph;
        Real z0 = r * cth;
        Real betax = beta[0] * sth*cph + beta[1] * cth*cph + beta[2] * sph;
        Real betay = beta[0] * sth*sph + beta[1] * cth*sph + beta[2] * cph;
        Real betaz = beta[0] * cth + beta[1] * sth;
        x0 += stheta*cos(phi) * r0 + betax * ct;
        y0 += stheta*sin(phi) * r0 + betay * ct;
        z0 += mu * r0 + betaz * ct;
        pphot->x[0] = sqrt(SQR(x0)+SQR(y0)+SQR(z0));
        pphot->x[1] = acos(z0 / pphot->x[0]);
        pphot->x[2] = atan2(y0,x0);
        if (pphot->x[2] < 0.)
          pphot->x[2] += 2.*PI;
      }

      // Check if photon has left original zone and update
      bool newzone = UpdateZone(pphot,0); //SWDFIX
      if (newzone) {
        // Check if photon is absorbed or escape due to boundary condition
        if (pphot->status != EVOLVING)
          return false;
      }
      if (newzone || compton) {
        // update opacity if zone or energy has changed
        pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot,ip);
        pphot->sct_coef = pmcb->ScatteringOpacity(pmcb,pphot,ip);
      }

      // update direction assuming isotropic random direction in comoving frame
      mu = 2.*pran->uniform()-1.0;
      stheta = sqrt(1.0-mu*mu);
      phi = 2.*PI*pran->uniform();
      pphot->k[0] = stheta*cos(phi);
      pphot->k[1] = stheta*sin(phi);
      pphot->k[2] = mu;
      if (boosts) {
        //transform back to Eulerian frame
        for(int i=0; i<3; ++i) {
          // undo multiply by gamma above and flip sign
          beta[i] /= -gamma;
        }
        if(beta2 > 0.) {
          Real bdk = (pphot->k[0] * beta[0] + pphot->k[1] * beta[1] +
                      pphot->k[2] * beta[2]);
          Real gonembdk = gamma * (1. - bdk);
          Real aber = gamma*(1.-gamma*bdk/(gamma+1.));
          pphot->k[0] = (pphot->k[0] - aber * beta[0]) / gonembdk;
          pphot->k[1] = (pphot->k[1] - aber * beta[1]) / gonembdk;
          pphot->k[2] = (pphot->k[2] - aber * beta[2]) / gonembdk;
          pphot->energy *= gonembdk;
          pphot->abs_coef /= gonembdk;
          pphot->sct_coef /= gonembdk;
        }
      }

      if (pphot->IsNanPhoton())
        pphot->PrintPhoton();

    } else {
      // return properties to eulerian frame if modified
      if (boosts) { // should always be true if !accel_success
        pphot->energy /= gonembdk;
        pphot->abs_coef *= gonembdk;
        pphot->sct_coef *= gonembdk;
      }
    }
  }
  return accel_success;

}

//----------------------------------------------------------------------------------------
//! \fn Real PhotonMover::GetOpticalDepth(MCRandom *pran)
//! \brief return exponentially distributed optical depth variable

Real PhotonMover::GetOpticalDepth(MCRandom *pran) {

  Real dev = pran->uniform();
  while(dev <= 0.)
    dev=pran->uniform();
  //std::cout << dev << std::endl;
  return -log(dev);
}

//----------------------------------------------------------------------------------------
//! \fn Real PhotonMover::GetExtinctionCoefficient(Real ac, Real sc)
//! \brief returns total opacity or scattering opacity depending on method

Real PhotonMover::GetExtinctionCoefficient(Real ac, Real sc) {

  Real chi;
  if (pmy_mcb->absorption_meth == ABSTAU) {
    chi = sc;
  } else {
    chi = sc + ac;
  }
  return (chi > TINY_NUMBER) ? chi : TINY_NUMBER;
}

//----------------------------------------------------------------------------------------
//! \fn Real PhotonMover::ExpTauAbsorption(Real ac, Real dl)
//! \brief Computes e^-tau_abs

Real PhotonMover::ExpTauAbsorption(Real ac, Real dl) {

  if (pmy_mcb->absorption_meth == ABSTAU) {
    return exp(-ac * dl);
  } else {
    return 1.;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonMover::NextFace(Real dx1, Real dx2, Real dx3, int &face, Real &dx)
//! \brief returns flag with next face and distance to next face

void PhotonMover::NextFace(Real dx1, Real dx2, Real dx3, int &face, Real &dx) {

// face tells which cell coordinates need to be updatde
//   x:   0
//   y:   1
//   z:   2
//   xy:  3
//   yz:  4
//   xz:  5
//   xyz: 6

  // check for positiviity
  /*if (dx1 < 0.) {
    dx1 = HUGE_NUMBER;
    printf("Warning: dx1 < 0\n");
  }
  if (dx2 < 0.) {
    dx2 = HUGE_NUMBER;
    printf("Warning: dx2 < 0\n");
  }
  if (dx3 < 0.) {
    dx3 = HUGE_NUMBER;
    printf("Warning: dx3 < 0\n");
    }*/

  dx = dx1;

  if(dx2 < dx) {
    dx = dx2;
    if(dx3 < dx) {
      dx = dx3;
      face = 2;
      return;
    } else if(dx3 > dx) {
      face = 1;
      return;
    } else {
      face = 4;
      return;
    }
  } else if(dx2 > dx) {
    if(dx3 < dx) {
      dx = dx3;
      face = 2;
      return;
    } else if(dx3 > dx) {
      face = 0;
      return;
    } else {
      face = 5;
      return;
    }
  } else {
    if(dx3 < dx) {
      dx = dx3;
      face = 2;
      return;
    } else if(dx3 > dx) {
      face = 3;
      return;
    } else {
      face = 6;
      return;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonMover::MovePhotonToNextZone(Photon *pphot, MCCoord *pco,
//!                           MonteCarloBlock *pmcb, int face, bool ascend[3], int ip))
//! \brief updates photon zone when face is known

void PhotonMover::MovePhotonToNextZone(Photon *pphot, MCCoord *pco, MonteCarloBlock *pmcb,
                                       int face, bool ascend[3], int ip) {

  // Update face(s) and adjust positions to lie exactly on boundary
  if ((face == 0) || (face == 3) || (face == 5) || (face == 6)) {
    //update x1 face
    if (ascend[0]) {
      pphot->i1p[ip]++;
      if(pphot->i1p[ip] <= pmcb->ie)
        pphot->x1p[ip] = pco->x1f(pphot->i1p[ip]);
      else {
        pmcb->pbval->BoundaryFunction_[BoundaryFace::outer_x1](pmcb,pco,pphot,ip);
      }
    } else {
      pphot->i1p[ip]--;
      if(pphot->i1p[ip] >= pmcb->is)
        pphot->x1p[ip] = pco->x1f(pphot->i1p[ip]+1);
      else {
        pmcb->pbval->BoundaryFunction_[BoundaryFace::inner_x1](pmcb,pco,pphot,ip);
      }
    }
  }
  if ((face == 1) || (face == 3) || (face == 4) || (face == 6)) {
    //update x2 face
    if (ascend[1]) {
      pphot->i2p[ip]++;
      if(pphot->i2p[ip] <= pmcb->je)
        pphot->x2p[ip] = pco->x2f(pphot->i2p[ip]);
      else {
        pmcb->pbval->BoundaryFunction_[BoundaryFace::outer_x2](pmcb,pco,pphot,ip);
      }
    } else {
      pphot->i2p[ip]--;
      if(pphot->i2p[ip] >= pmcb->js)
        pphot->x2p[ip] = pco->x2f(pphot->i2p[ip]+1);
      else {
        pmcb->pbval->BoundaryFunction_[BoundaryFace::inner_x2](pmcb,pco,pphot,ip);
      }
    }
  }
  if ((face == 2) || (face == 4) || (face == 5) || (face == 6)) {
    //update x3 face
    if (ascend[2]) {
      pphot->i3p[ip]++;
      if(pphot->i3p[ip] <= pmcb->ke)
        pphot->x3p[ip] = pco->x3f(pphot->i3p[ip]);
      else {
        pmcb->pbval->BoundaryFunction_[BoundaryFace::outer_x3](pmcb,pco,pphot,ip);
      }
    } else {
      pphot->i3p[ip]--;
      if(pphot->i3p[ip] >= pmcb->ks)
        pphot->x3p[ip] = pco->x3f(pphot->i3p[ip]+1);
      else {
        pmcb->pbval->BoundaryFunction_[BoundaryFace::inner_x3](pmcb,pco,pphot,ip);
      }
    }
  }

  // Update opacities
  if (pphot->statp[ip] == EVOLVING) {
    // Opacities need to be calculated using comoving frame energy and then transformed
    // back to Eulerian frame when Lorentz Transformations are enabled.
    int &i1 = pphot->i1p[ip];
    int &i2 = pphot->i2p[ip];
    int &i3 = pphot->i3p[ip];
    if (pmy_mcb->boosts) {
      // Shift photon energy to comoving frame
      Real shift = pmy_mcb->LorentzTransformFrequencyShift(pphot,ip);
      pphot->ep[ip] *= shift;
      // compute opacities in comoving frame
      pphot->acp[ip] = pmcb->AbsorptionOpacity(pmcb,pphot,ip);
      pphot->scp[ip] = pmcb->ScatteringOpacity(pmcb,pphot,ip);
      // Shift energy back to Eulerian frame
      pphot->ep[ip] /= shift;
      // Shift opacities to Eulerian frame
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
//! \fn bool PhotonMover::UpdateZone(photon *pphot, int ip)
//! \brief check/updates photon zone after displacement

bool PhotonMover::UpdateZone(Photon *pphot, int ip) {

  bool change = false;
  MonteCarloBlock *pmcb = pmy_mcb;
  bool update = false;

  if (pphot->x1p[ip] >= pcoord->x1f(pphot->i1p[ip]+1)) {
    update = true;
    while (pphot->x1p[ip] >= pcoord->x1f(pphot->i1p[ip]+1)) {
      pphot->i1p[ip]++;
      if(pphot->i1p[ip] > pmcb->ie)
        pmcb->pbval->BoundaryFunction_[BoundaryFace::outer_x1](pmcb,pcoord,pphot,ip);
      if (pphot->statp[ip] != EVOLVING) {
        break;
      }
    }
  } else if (pphot->x1p[ip] < pcoord->x1f(pphot->i1p[ip])) {
    update = true;
    while (pphot->x1p[ip] < pcoord->x1f(pphot->i1p[ip])) {
      pphot->i1p[ip]--;
      if(pphot->i1p[ip] < pmcb->is)
        pmcb->pbval->BoundaryFunction_[BoundaryFace::inner_x1](pmcb,pcoord,pphot,ip);
      if (pphot->statp[ip] != EVOLVING) {
        break;
      }
    }
  }
  if (pphot->x2p[ip] >= pcoord->x2f(pphot->i2p[ip]+1)) {
    update = true;
    while (pphot->x2p[ip] >= pcoord->x2f(pphot->i2p[ip]+1)) {
      pphot->i2p[ip]++;
      if(pphot->i2p[ip] > pmcb->je)
        pmcb->pbval->BoundaryFunction_[BoundaryFace::outer_x2](pmcb,pcoord,pphot,ip);
      if (pphot->statp[ip] != EVOLVING) {
        break;
      }
    }
  } else if (pphot->x2p[ip] < pcoord->x2f(pphot->i2p[ip])) {
    update = true;
    while (pphot->x2p[ip] < pcoord->x2f(pphot->i2p[ip])) {
      pphot->i2p[ip]--;
      if(pphot->i2p[ip] < pmcb->js)
        pmcb->pbval->BoundaryFunction_[BoundaryFace::inner_x2](pmcb,pcoord,pphot,ip);
      if (pphot->statp[ip] != EVOLVING) {
        break;
      }
    }
  }
  if (pphot->x3p[ip] >= pcoord->x3f(pphot->i3p[ip]+1)) {
    update = true;
    while (pphot->x3p[ip] >= pcoord->x3f(pphot->i3p[ip]+1)) {
      pphot->i3p[ip]++;
      if(pphot->i3p[ip] > pmcb->ke)
        pmcb->pbval->BoundaryFunction_[BoundaryFace::outer_x3](pmcb,pcoord,pphot,ip);
      if (pphot->statp[ip] != EVOLVING) {
        break;
      }
    }
  } else if (pphot->x3p[ip] < pcoord->x3f(pphot->i3p[ip])) {
    update = true;
    while (pphot->x3p[ip] < pcoord->x3f(pphot->i3p[ip])) {
      pphot->i3p[ip]--;
      if(pphot->i3p[ip] < pmcb->ks)
        pmcb->pbval->BoundaryFunction_[BoundaryFace::inner_x3](pmcb,pcoord,pphot,ip);
      if (pphot->statp[ip] != EVOLVING) {
        break;
      }
    }
  }
  // Returns true if zone changes, false otherwise
  return update;

}


//----------------------------------------------------------------------------------------
//! \fn void PhotonMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4])
//! \brief convert k vector from curvalinear to cartesian

void PhotonMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4]) {

  // Default corresponds to Cartesian so just copy
  for (int i=0; i<4; ++i)
    kcart[i] = pphot->k[i];
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonMover::InitializeMWDist(void)
//! \brief initialize modified randon walk path length distribution

void PhotonMover::InitializeMRWDist(void) {

  nmax = 1000;
  mrwprob.NewAthenaArray(nmax);
  mrwdev.NewAthenaArray(nmax);

  for(int i=0; i<nmax; ++i) {
    mrwdev(i) = static_cast<Real>(i)/static_cast<Real>(nmax-1);
    mrwprob(i) = 0.0;
  }
  for(int i=0; i<nmax-1; ++i) {
    int n = 1;
    Real sign = 1.0;
    Real yn2 = pow(mrwdev(i),n*n);
    // Compute the sum to the limit of Real precision
    while (yn2 > 1.e-17 ) {
      mrwprob(i) += sign * yn2;
      sign *= -1.0;
      n += 1;
      yn2 = pow(mrwdev(i),n*n);
    }
    mrwprob(i) *= 2.;
    if (mrwprob(i) > 1.0) mrwprob(i)=1.0;
  }
  mrwprob(nmax-1) = 1.0;


}

//----------------------------------------------------------------------------------------
//! \fn Real MRWDist(MCRandom *pran)
//! \brief get modified randon walk path length

Real PhotonMover::MRWDist(MCRandom *pran) {

  Real x0 = pran->uniform();

  // Perform a binary search
  int low =0, high = nmax-1, mid;
  while(low<=high) {
    mid=(low+high)/2;
    if(mrwprob(mid-1) <= x0) {
      if(mrwprob(mid) > x0)
        break;
      else
        low=mid+1;
    }
    else
      high=mid-1;
  }

  // Replace binary search with initial guess ?
  //int i = static_cast<int>(x0*static_cast<Real>(nmax));

  // use linear interpolation to find location
  if (mid == 0)
    return mrwdev(0);
  else if (low == nmax) {
    return mrwdev(nmax-1);
  } else {
    Real slope = (x0 - mrwprob(mid-1)) / (mrwprob(mid) - mrwprob(mid-1));
    return mrwdev(mid-1)+(mrwdev(mid)-mrwdev(mid-1)) * slope;
  }


}

//----------------------------------------------------------------------------------------
//! \fn void  PhotonMover::ReadComptonGreensFunction()
//! \brief Reads in pre-tabulated binary table used in MRW acceleration with Compton scat.

void PhotonMover::ReadComptonGreensFunction(void) {

  int nt=200,nxi=100,np=100;
  //int nt=50,nxi=2,np=100;

  FILE *pfile;
  double fdat;
  // Read in time array
  if((pfile = fopen("compton_table_t.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadComptonGreensFunction]"
          <<std::endl<< "Input file compton_table_t.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwt.NewAthenaArray(nt);
  for (int i=0; i<nt; i++) {
    fread(&fdat, sizeof(double), 1, pfile);
    mrwt(i) = static_cast<Real>(fdat);
    //printf("t: %15.8e \n",exp(mrwt(i)));
  }
  fclose(pfile);
  // Read in initial dimensionless energy array
  if((pfile = fopen("compton_table_x1.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadComptonGreensFunction]"
          <<std::endl<< "Input file compton_table_x1.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwxi.NewAthenaArray(nxi);
  for (int i=0; i<nxi; i++) {
    fread(&fdat, sizeof(double), 1, pfile);
    mrwxi(i) = static_cast<Real>(fdat);
    //printf("xi: %15.8e \n",exp(mrwxi(i)));
  }
  fclose(pfile);
  // Read in probability array
  if((pfile = fopen("compton_table_p.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadComptonGreensFunction]"
          <<std::endl<< "Input file compton_table_p.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwp.NewAthenaArray(np);
  for (int i=0; i<np; i++) {
    fread(&fdat, sizeof(double), 1, pfile);
    mrwp(i) = static_cast<Real>(fdat);
  }
  fclose(pfile);
  // Read in output dimensionless energy array
  if((pfile = fopen("compton_table_x.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadComptonGreensFunction]"
          <<std::endl<< "Input file compton_table_x.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwxf.NewAthenaArray(np,nxi,nt);
  for (int k=0; k<np; k++) {
    for (int j=0; j<nxi; j++) {
      for (int i=0; i<nt; i++) {
        fread(&fdat, sizeof(double), 1, pfile);
        mrwxf(k,j,i) = static_cast<Real>(fdat);
      }}}
  fclose(pfile);

}


//----------------------------------------------------------------------------------------
//! \fn void  PhotonMover::ReadRadiusDistribution()
//! \brief Reads in pre-tabulated binary table used in MRW acceleration with advection

void PhotonMover::ReadRadiusDistribution(void) {

  int ny=100,np=100;

  FILE *pfile;
  double fdat;
  // Read in y (y=exp(-t)) array
  if((pfile = fopen("radius_table_t.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadRadiusDistribution]"
          <<std::endl<< "Input file radius_table_t.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwrt.NewAthenaArray(ny);
  for (int i=0; i<ny; i++) {
    fread(&fdat, sizeof(double), 1, pfile);
    mrwrt(i) = static_cast<Real>(fdat);
  }
  fclose(pfile);
  // Read in probability array
  if((pfile = fopen("radius_table_p.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadTRadiusDistribution]"
          <<std::endl<< "Input file radius_table_p.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwrp.NewAthenaArray(np);
  for (int i=0; i<np; i++) {
    fread(&fdat, sizeof(double), 1, pfile);
    mrwrp(i) = static_cast<Real>(fdat);
  }
  fclose(pfile);
  // Read in output dimensionless radius
  if((pfile = fopen("radius_table_r.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadRadiusDistribution]"
          <<std::endl<< "Input file radius_table_r.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwrr.NewAthenaArray(np,ny);
  for (int j=0; j<np; j++) {
      for (int i=0; i<nt; i++) {
        fread(&fdat, sizeof(double), 1, pfile);
        mrwrr(j,i) = static_cast<Real>(fdat);
      }}
  fclose(pfile);

}


//----------------------------------------------------------------------------------------
//! \fn void  PhotonMover::ReadTimeDistribution()
//! \brief Reads in pre-tabulated binary table used in MRW acceleration with advection

void PhotonMover::ReadTimeDistribution(void) {

  int ntau=100,np=400;

  FILE *pfile;
  double fdat;

  if((pfile = fopen("time_table_tau.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadTimeDistribution]"
          <<std::endl<< "Input file time_table_y.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwta.NewAthenaArray(ntau);
  for (int i=0; i<ntau; i++) {
    fread(&fdat, sizeof(double), 1, pfile);
    mrwta(i) = static_cast<Real>(fdat);
  }
  fclose(pfile);
  // Read in probability array
  if((pfile = fopen("time_table_p.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadTimeDistribution]"
          <<std::endl<< "Input file time_table_p.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwtp.NewAthenaArray(np);
  for (int i=0; i<np; i++) {
    fread(&fdat, sizeof(double), 1, pfile);
    mrwtp(i) = static_cast<Real>(fdat);
  }
  fclose(pfile);
  // Read in output dimensionless radius
  if((pfile = fopen("time_table_t.out","r")) == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [PhotonMover:ReadTimeDistribution]"
          <<std::endl<< "Input file time_table_t.out could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  mrwtt.NewAthenaArray(np,ntau);
  for (int j=0; j<np; j++) {
      for (int i=0; i<ntau; i++) {
        fread(&fdat, sizeof(double), 1, pfile);
        mrwtt(j,i) = static_cast<Real>(fdat);
      }}
  fclose(pfile);

}

//----------------------------------------------------------------------------------------
//! \fn Real PhotonMover::InterpComptonEnergy(Real x0, Real time, Real prob)
//! \brief return energy for photon after MRW with compton scattering

Real PhotonMover::InterpComptonEnergy(Real xi, Real time, Real prob) {

  // get interpolant for prob
  int ip = mcbisect(prob,mrwp);
  Real a = (prob-mrwp(ip))/(mrwp(ip+1)-mrwp(ip));
  Real a1 = 1.-a;

  if (log(xi) > mrwxi(mrwxi.GetDim1()-1)) {
    printf("xi: %g %g\n",xi, exp(mrwxi(mrwxi.GetDim1()-1)));
  }
  // get interpolant for x0
  int ixi = mcbisect(log(xi),mrwxi);
  Real b = (log(xi)-mrwxi(ixi))/(mrwxi(ixi+1)-mrwxi(ixi));
  Real b1 = 1.-b;

  if (log(time) > mrwt(mrwt.GetDim1()-1)) {
    printf("time: %g %g\n",time, exp(mrwt(mrwt.GetDim1()-1)));
  }
  // get interpolant for t
  int it = mcbisect(log(time),mrwt);
  Real c = (log(time)-mrwt(it))/(mrwt(it+1)-mrwt(it));
  Real c1 = 1.-c;

  /*printf("%d %d %d\n",ip,ixi,it);
  printf("%d %g %g %g %g\n",ip,a,a1,mrwp(ip),mrwp(ip+1));
  printf("%d %g %g %g %g\n",ixi,b,b1,mrwxi(ixi),mrwxi(ixi+1));
  printf("%d %g %g %g %g\n",it,c,c1,mrwt(it),mrwt(it+1));
  printf("final: %g\n",(b*(c*(a*mrwxf(ip+1,ixi+1,it+1)+a1*mrwxf(ip,ixi+1,it+1))
         +c1*(a*mrwxf(ip+1,ixi+1,it)+a1*mrwxf(ip,ixi+1,it)))
         +b1*(c*(a*mrwxf(ip+1,ixi,it+1)+a1*mrwxf(ip,ixi,it+1))
         +c1*(a*mrwxf(ip+1,ixi,it)+a1*mrwxf(ip,ixi,it)))));*/
  return (b*(c*(a*mrwxf(ip+1,ixi+1,it+1)+a1*mrwxf(ip,ixi+1,it+1))
          +c1*(a*mrwxf(ip+1,ixi+1,it)+a1*mrwxf(ip,ixi+1,it)))
          +b1*(c*(a*mrwxf(ip+1,ixi,it+1)+a1*mrwxf(ip,ixi,it+1))
          +c1*(a*mrwxf(ip+1,ixi,it)+a1*mrwxf(ip,ixi,it))));
}

//----------------------------------------------------------------------------------------
//! \fn Real PhotonMover::InterpPathTime(Real tau, Real prob)
//! \brief return time/path for photon undergoing MRW at given tau

Real PhotonMover::InterpPathTime(Real tau, Real prob) {

  // get interpolant for prob
  int ip = mcbisect(prob,mrwtp);
  //int ip = static_cast<int>(prob*static_cast<Real>(200));
  //if (ip == 200) ip = 199;
  Real a = (prob-mrwtp(ip))/(mrwtp(ip+1)-mrwtp(ip));
  Real a1 = 1.-a;

  // get interpolant for tau
  int it = mcbisect(log(tau),mrwta);
  Real b = (log(tau)-mrwta(it))/(mrwta(it+1)-mrwta(it));
  Real b1 = 1.-b;
  //printf("%d %d %g %g %g %g\n",ip,it,mrwtt(ip+1,it+1),mrwtt(ip,it+1),
  //	 mrwtt(ip+1,it),mrwtt(ip,it));
  return b*(a*mrwtt(ip+1,it+1)+a1*mrwtt(ip,it+1))
         +b1*(a*mrwtt(ip+1,it)+a1*mrwtt(ip,it));

}
