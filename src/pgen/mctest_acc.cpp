//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mctest_acc.cpp
//  \brief Problem generator for simple monte carlo problem
//
//========================================================================================

#include <iostream> // temporary for testing

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../monte_carlo/montecarlo.hpp"
#include "../monte_carlo/photon.hpp"
#include "../monte_carlo/photonmover.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

static Real lmax;
static Real energy0;
static bool first = true;
//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief monte carlo test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real rideal = 8.314e7;
  Real temp = pin->GetReal("problem","temp");
  Real rho = pin->GetReal("problem","rho");
  Real vel = pin->GetOrAddReal("problem","velocity",0.);
  Real gamma = peos->GetGamma();
  Real c = 2.9979e10;
  vel *= c;

  Real xlow = pin->GetReal("mesh","x3min");
  Real xhigh = pin->GetReal("mesh","x3max");
  lmax = xhigh-xlow;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {       
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = rho*vel;
        phydro->u(IEN,k,j,i) = rideal*rho*temp/(gamma-1.0);
      }
    }
  }
  // add kinetic energy
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
      }}}

}

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){

}

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin){
  

  // Set codetocgs here
  //EnrollUserEmissionFunction();
  
}

void MonteCarloBlock::InitializePhoton(Photon *pphot) {


  MCRandom *pran = pphot->pmy_mcb->pran;
  MCCoord *pco = pphot->pmy_mcb->pcoord;

  // Set status flag
  pphot->status = EVOLVING;


  //pphot->energy = energy0;

  pphot->x[0] = 0.0;
  pphot->i1 = -1;
  for(int i=pphot->pmy_mcb->is; i<=pphot->pmy_mcb->ie; i++) {
    if ((pphot->x[0] > pco->x1f(i)) && (pphot->x[0] <= pco->x1f(i+1)))
      pphot->i1 = i;
  }
  if (pphot->i1 < 0) pphot->weight = -1.0;
 
  pphot->x[1] = 0.0;
  pphot->i2 = -1;
  for(int i=pphot->pmy_mcb->js; i<=pphot->pmy_mcb->je; i++) {
    if ((pphot->x[1] > pco->x2f(i)) && (pphot->x[1] <= pco->x2f(i+1)))
      pphot->i2 = i;
  }
  if (pphot->i2 < 0) pphot->weight = -1.0;
 
  pphot->x[2] = 0.0;
  pphot->i3 = -1;
  for(int i=pphot->pmy_mcb->ks; i<=pphot->pmy_mcb->ke; i++) {
    if ((pphot->x[2] > pco->x3f(i)) && (pphot->x[2] <= pco->x3f(i+1)))
      pphot->i3 = i;
  }
  if (pphot->i3 < 0) pphot->weight = -1.0;


  // Utilize free-free emission function in emission.cpp
  //PhotonEmitFreeFree(this,pphot);

  // Initialize Stokes vector
  pphot->stokes[0] = 1.0;
  pphot->stokes[1] = 0.0;
  pphot->stokes[2] = 0.0;

  // Generate initial angle parameters
  Real phi = 2. * PI * pran->uniform();
  Real cphi = cos(phi);
  Real sphi = sin(phi);

  Real cth = 2. * pran->uniform() - 1.;
  Real sth = sqrt(1. - SQR(cth));

  // Initialize wave vector with isotropic distribution
  pphot->kcart[0] = sth*cphi;
  pphot->kcart[1] = sth*sphi;
  pphot->kcart[2] = cth;
  
  // eweight is a constant weighting factor which accounts for the
  // emissivity of the grid zone in which the photon was emitted
  pphot->eweight = 1.0;
  pphot->weight = 1.0;

  Real dev = exp(2.30258509299*(elog*pran->uniform()+eminlog));  
  pphot->energy = dev;

  if (pphot->weight < 0.0) { 
    pphot->status = DESTROYED;
  } else {
    // Initialize the absorption and scattering extinction coefficients
    // to the values appropriate in the emitted zone
    pphot->abs_coef = AbsorptionOpacity(this,pphot);
    pphot->sct_coef = ScatteringOpacity(this,pphot);
    if (first) {
      printf("Tau_es: %g\n",pphot->sct_coef*lmax);
      printf("Tau_abs: %g\n",pphot->abs_coef*lmax);
      first = false;
    }
  }
}

