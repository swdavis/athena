//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!  \file emission.cpp

// Athena++ headers
#include "montecarlo.hpp"
#include "../defs.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn void InitializeEmissionFreefree(MonteCarloBlock *pmbc, MeshBlock *pmb)
//  \brief Initialize emission array for static monte carlo calculation

void InitializeEmissionFreeFree(MonteCarloBlock *pmcb) {
  //void InitializeEmissionFreeFree(MonteCarloBlock *pmcb, MeshBlock *pmb) {


  Real heabund = 0.0; // Should have more general EOS functions
  Real mp = 1.6726e-24;
  Real eta0 = 1.032521e-11;
  Real g = 1.0;

  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        Real temp = pmcb->tgas(k,j,i);
        Real nhii = pmcb->rho(k,j,i)/mp/(1.0+4.0*heabund);
        Real ne = (1.0+2.0*heabund) * nhii;
        Real vol = pmcb->pmy_coord->GetCellVolume(k,j,i);
        //std::cout << vol << std::endl;
        pmcb->emission(k,j,i) = eta0/sqrt(temp)*ne*nhii*g*vol;
      }}}
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonInitFreeFree(MonteCarloBlock *pmcb, Photon *pphot)
//  \brief initialize energy, direction, polarization and weight of the photon
//         consistent with free-free emission

void PhotonEmitFreeFree(MonteCarloBlock *pmcb, Photon *pphot)
{
  Real kb = 1.3807e-16;
  MCRandom *pran = pmcb->pran;

  // Scheme in which packets are drawn from a uniform distribution in log E
  // requires weight = exp(-x) note log(10)=2.30258509299
  Real dev = exp(2.30258509299*(pmcb->elog*pran->uniform()+pmcb->eminlog));  

  pphot->energy = dev;
  Real x = dev / (kb * pmcb->tgas(pphot->i3,pphot->i2,pphot->i1));

  // Initialize weight
  pphot->weight *= exp(-x) * 2.30258509299 * pmcb->elog;

  // Initialize Stokes vector
  pphot->stokes[0] = 1.0;
  pphot->stokes[1] = 0.0;
  pphot->stokes[2] = 0.0;

  // Generate initial angle parameters
  Real cphi = cos(2. * PI * pran->uniform());
  Real sphi = sqrt(1. - SQR(cphi));

  Real cth = 2. * pran->uniform() - 1.;
  Real sth = sqrt(1. - SQR(cth));

  // Initialize wave vector with isotropic distribution
  pphot->k[0] = sth*cphi;
  pphot->k[1] = sth*sphi;
  pphot->k[2] = cth;
}


//----------------------------------------------------------------------------------------
//! \fn void GetZonePositionCartesian(Photon *pphot, MCRandom *pran, Coordinates *pcoord)
//  \brief choose random position within cartesian gridzone

void GetZonePositionCartesian(Photon *pphot, MCRandom *pran, Coordinates *pcoord) {

  Real xl = pcoord->x1f(pphot->i1), dx = pcoord->dx1f(pphot->i1);
  Real yl = pcoord->x2f(pphot->i2), dy = pcoord->dx2f(pphot->i2);
  Real zl = pcoord->x3f(pphot->i3), dz = pcoord->dx3f(pphot->i3);

  pphot->x[0] = xl+pran->uniform()*dx;
  pphot->x[1] = yl+pran->uniform()*dy;
  pphot->x[2] = zl+pran->uniform()*dz;

}
