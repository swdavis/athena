//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!  \file emission.cpp
//!  \brief implementation of photon emission functions

// Athena++ headers
#include "montecarlo.hpp"
#include "../defs.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../globals.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real InitializeEmissionFreefree(MonteCarloBlock *pmcb)
//! \brief Initialize emission array and return minimum

Real InitializeEmissionFreeFree(MonteCarloBlock *pmcb) {

  Real heabund = 0.09; // Should have more general EOS functions
  Real kb = 1.380649e-16;
  Real mp = 1.67262192369e-24;
  Real eta0 = 1.032521e-11;
  Real g = 1.0; // Gaunt factor

  //eta0 *= 12.;  // Added to match the Athena++ prescription

  Real ncells = static_cast<Real>(pmcb->pmy_mc->ncells);
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  Real emm_min = SQR(HUGE_NUMBER);
  Real emm_max = -HUGE_NUMBER;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        Real temp = pmcb->tgas(k,j,i);
        Real nh = pmcb->rho(k,j,i)/mp/(1.+4.*heabund);
        Real nhe = nh*heabund;
        Real ne = nh + 2.*nhe;
        Real vol = pmcb->pcoord->vol(k,j,i);
        pmcb->emission(k,j,i) = eta0/sqrt(temp)*ne*(nh+4.*nhe)*g*vol*ncells;
        if (pmcb->emission(k,j,i) > emm_max) emm_max = pmcb->emission(k,j,i);
        if (pmcb->emission(k,j,i) < emm_min) emm_min = pmcb->emission(k,j,i);
      }}}
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE,&emm_min,1,MPI_ATHENA_REAL,MPI_MIN,
               MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&emm_max,1,MPI_ATHENA_REAL,MPI_MAX,
               MPI_COMM_WORLD);
#endif
  if (Globals::my_rank == 0) {
    std::cout << "Emission array range (min, max): " << emm_min << " " << emm_max
              << std::endl;
  }
  return emm_min;
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonInitFreeFree(MonteCarloBlock *pmcb, Photon *pphot, Real lemin,
//!                             Real lemax, int ips, int ipe))
//! \brief initialize photon consistent with free-free emission

void PhotonEmitFreeFree(MonteCarloBlock *pmcb, Photon *pphot, Real lemin, Real lemax,
                        int ip)
{
  Real kb = 1.380649e-16;
  MCRandom *pran = pmcb->pran;

  // Scheme in which packets are drawn from a uniform distribution in log E
  // requires weight = exp(-x) note log(10)=2.30258509299
  Real dev = exp((lemax-lemin)*pran->uniform()+lemin);
  pphot->ep[ip] = dev;
  Real x = dev / (kb * pmcb->tgas(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip]));

  // Initialize weight
  pphot->wp[ip] *= exp(-x) * (lemax-lemin);

  if (pmcb->pmy_mc->polarized) {
    // Initialize Stokes vector
    pphot->sip[ip] = 1.0;
    pphot->sup[ip] = 0.0;
    pphot->sqp[ip] = 0.0;
    pphot->svp[ip] = 0.0;
  }

  // Generate initial angle parameters
  Real phi = 2. * PI * pran->uniform();
  Real cphi = cos(phi);
  Real sphi = sin(phi);
  Real cth = 2. * pran->uniform() - 1.;
  Real sth = sqrt(1. - SQR(cth));

  // Initialize wave vector with isotropic distribution
  pphot->k0p[ip] = 1.;
  pphot->k1p[ip] = sth*cphi;
  pphot->k2p[ip] = sth*sphi;
  pphot->k3p[ip] = cth;
}


//----------------------------------------------------------------------------------------
//! \fn void GetZonePositionCartesian(Photon *pphot, MCRandom *pran, MCCoord *pcoord,
//!                                   int ip)
//! \brief choose random position within cartesian gridzone ip

void GetZonePositionCartesian(Photon *pphot, MCRandom *pran, MCCoord *pcoord, int ip) {

  Real xl = pcoord->x1f(pphot->i1p[ip]); Real dx = pcoord->x1f(pphot->i1p[ip]+1)-xl;
  Real yl = pcoord->x2f(pphot->i2p[ip]); Real dy = pcoord->x2f(pphot->i2p[ip]+1)-yl;
  Real zl = pcoord->x3f(pphot->i3p[ip]); Real dz = pcoord->x3f(pphot->i3p[ip]+1)-zl;

  pphot->x0p[ip] = 1.;
  pphot->x1p[ip] = xl+pran->uniform()*dx;
  pphot->x2p[ip] = yl+pran->uniform()*dy;
  pphot->x3p[ip] = zl+pran->uniform()*dz;

}

//----------------------------------------------------------------------------------------
//! \fn void GetZonePositionSphericalPolar(Photon *pphot, MCRandom *pran, MCCoord *pcoord,
//!                                        int ip)
//! \brief choose random position within spherical-polar gridzone ip

void GetZonePositionSphericalPolar(Photon *pphot, MCRandom *pran, MCCoord *pcoord,
                                   int ip) {
  Real rl = pcoord->x1f(pphot->i1p[ip]), rh = pcoord->x1f(pphot->i1p[ip]+1);
  pphot->x1p[ip]  = pow(pran->uniform()*(rh*rh*rh-rl*rl*rl)+rl*rl*rl,1./3.);
  Real cthh = cos(pcoord->x2f(pphot->i2p[ip]));
  Real cthl = cos(pcoord->x2f(pphot->i2p[ip]+1));
  Real cth = cthl + pran->uniform() * (cthh-cthl);
  pphot->x2p[ip] = acos(cth);
  Real pl = pcoord->x3f(pphot->i3p[ip]); Real dp = pcoord->x3f(pphot->i3p[ip]+1)-pl;
  pphot->x3p[ip] = pl+pran->uniform()*dp;

}

//----------------------------------------------------------------------------------------
//! \fn Real PlanckDist(Real temp, MCRandom *pran)
//! \brief returns energy distributed according to Planck function

Real PlanckDist(Real temp, MCRandom *pran)
{
  // Method of choosing the energy of the initial photon which a Planck spectrum
  // distribution. See Pozdnyakov et al. sec 9.4.  Originally, Fleck and Cumming (1971)

  Real x1 = 1.202 * pran->uniform();
  Real x2 = pran->uniform();
  Real x3 = pran->uniform();
  Real x4 = pran->uniform();

  Real sum = 1.0;
  int alpha = 1;
  while (x1 >= sum) {
    alpha++;
    sum = sum + 1.0 / (alpha * alpha * alpha);
  }

  Real kb = 1.380649e-16;
  return -kb * temp * log(x2 * x3 * x4) / ((Real)alpha);

}
