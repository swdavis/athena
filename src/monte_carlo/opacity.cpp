//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!  \file emission.cpp

// Athena++ headers
#include "montecarlo.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real NoOpacity(MonteCarloBlock *pmcb, Photon *pphot)
//  \brief return zero for extinction coeffictent

Real NoOpacity(MonteCarloBlock *pmcb, Photon *pphot) {

  return 0.;
}

//----------------------------------------------------------------------------------------
//! \fn Real FreeFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot)
//  \brief calculation extinction coefficient for free-free absorption

Real FreeFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot) {

  Real ffnrm = 3.692146e8;
  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.6726e-24; 
  Real h =  6.6262e-27;
  Real kb = 1.3807e-16;

  Real nhii = pmcb->rho(pphot->i3,pphot->i2,pphot->i1) / (mp*(1.+4.*heabund));
  Real ne = (1. + 2.*heabund) * nhii;

  Real nu = pphot->energy / h;
  Real temp = pmcb->tgas(pphot->i3,pphot->i2,pphot->i1);
  Real ehnu = exp(-pphot->energy / (kb * temp) );

  Real aff = ffnrm/sqrt(temp)/pow(nu,3);
 
  return ne * nhii * aff * (1. - ehnu);

}

//----------------------------------------------------------------------------------------
//! \fn Real ThomsonOpacity(MonteCarloBlock *pmcb, Photon *pphot)
//  \brief calculation extinction coefficient for Thomson scattering

Real ThomsonOpacity(MonteCarloBlock *pmcb, Photon *pphot) {

  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.6726e-24; 
  Real sigmat = 6.65248e-25;

  Real kappaes = sigmat * (1. + 2.*heabund) / (mp * (1.+4.*heabund) );

  return kappaes * pmcb->rho(pphot->i3,pphot->i2,pphot->i1);
}
