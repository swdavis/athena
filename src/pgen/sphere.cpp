//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mctest.cpp
//  \brief Problem generator for  monte carlo through uniform sphere
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

static Real rad0,path0;
static Real energy0;
static bool srcdist;
static bool first = true;
//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief monte carlo test problem generator
//========================================================================================
void SphericalOrTimedEscape(MonteCarloBlock *pmcb, Photon *phot);

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real rideal = 8.314e7;
  Real c = 2.9979e10;
  Real temp = pin->GetReal("problem","temp");
  Real tau = pin->GetReal("problem","tau");
  Real rad0 = pin->GetReal("problem","rad0");
  Real vel = pin->GetOrAddReal("problem","velocity",0.);
  Real gamma = peos->GetGamma();
  vel *= c;

  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.6726e-24; 
  Real sigmat = 6.65248e-25;
  Real kappaes = sigmat * (1. + 2.*heabund) / (mp * (1.+4.*heabund) );

  Real rho = tau / (kappaes * rad0);
  printf("rho: %g %g %g %g\n",rho,kappaes,rad0,tau);
  // density is non-zero only in sphere
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

void MonteCarloBlock::InitUserMonteCarloBlockData(ParameterInput *pin){

  // Set variables 
  srcdist =pin->GetOrAddBoolean("problem","srcdist",false);
  rad0 = pin->GetReal("problem","rad0");
  path0 = pin->GetOrAddReal("problem","path0",HUGE_NUMBER);
  Real xi = pin->GetReal("problem","energy0");
  Real temp = pin->GetReal("problem","temp");
  Real kb = 1.3807e-16;
  energy0 = kb*temp*xi;
  printf("Energy initial (keV, x): %g %g\n",energy0/1.6021772e-12/1000.,xi);

  // enroll function to cease photon propogation based on escape radius
  // or total integration time
  EnrollUserWorkInMove(SphericalOrTimedEscape);
}

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin){
  
}

void MonteCarloBlock::InitializePhoton(Photon *pphot) {


  MCCoord *pco = pphot->pmy_mcb->pcoord;
  // Set status flag
  pphot->status = EVOLVING;

  // Choose random intial position, weights, energy, and direction
  // for photon emission.  In this version an equal number of photons
  // is emitted in  each grid zone.  The relative emission from each grid 
  // zone is then accounted for by a weighting factor cweight. 

  // cweight is a constant weighting factor which accounts for the
  // emissivity of the grid zone in which the photon was emitted
  if (zone_weight_flag) {
    pphot->eweight = 1.0;
    pphot->weight = 1.0;
  }

  pphot->energy = energy0;
  pphot->path = 0.;

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
 
  pphot->x[2] = 0.;
  pphot->i3 = -1;
  for(int i=pphot->pmy_mcb->ks; i<=pphot->pmy_mcb->ke; i++) {
    if ((pphot->x[2] > pco->x3f(i)) && (pphot->x[2] <= pco->x3f(i+1)))
      pphot->i3 = i;
  }
  if (pphot->i3 < 0) pphot->weight = -1.0;

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
  pphot->k[0] = sth*cphi;
  pphot->k[1] = sth*sphi;
  pphot->k[2] = cth;

  if (srcdist) {
    Real x = pran->uniform();
    while (x <= 0.)
      x = pran->uniform();
    Real dev = pran->uniform();
    //while (sin(PI*x)*x < dev) {
    // Yields nearly exponential escape time distribution
    while (sin(PI*x)/(x*PI) < dev) {
      x = pran->uniform();
      while (x <= 0.)
	x = pran->uniform();
      dev = pran->uniform();
    }
    //printf("x: %g\n",x);
    Real r0 = x*rad0;
    phi = 2. * PI * pran->uniform();
    cphi = cos(phi);
    sphi = sin(phi);

    cth = 2. * pran->uniform() - 1.;
    sth = sqrt(1. - SQR(cth));
    pphot->x[0] += r0*sth*cphi;
    pphot->x[1] += r0*sth*sphi;
    pphot->x[2] += r0*cth;
  }
  
  if (pphot->weight < 0.0) pphot->status = DESTROYED;
  
  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = AbsorptionOpacity(this,pphot);
  pphot->sct_coef = ScatteringOpacity(this,pphot);
  //pphot->PrintPhoton();
  if (first) {
    printf("%g \n",pphot->abs_coef);
    printf("taus, taua: %g %g\n",pphot->sct_coef*rad0,pphot->abs_coef*rad0);
    first = false;
  }
}

void SphericalOrTimedEscape(MonteCarloBlock *pmcb, Photon *pphot) {

  // First check radius condition
  Real r = sqrt(SQR(pphot->x[0])+SQR(pphot->x[1])+SQR(pphot->x[2]));
  if (r >= rad0) {
    Real dr = r-rad0;
    pphot->path -= dr;
    for (int i=0; i<3; i++) {
      // assume cartesian for now
      pphot->x[i] -= pphot->k[i]*dr;
    }
    pphot->status = ESCAPED;
    pphot->face = FACE_UNDEF;
  } 
  // Then check path condition -- ensures path is not over estimated
  if (pphot->path >= path0) {
    Real dp = pphot->path - path0;
    pphot->path = path0;
    for (int i=0; i<3; i++) {
      // assume cartesian for now
      pphot->x[i] -= pphot->k[i]*dp;
    }
    pphot->status = ESCAPED;
    pphot->face = FACE_UNDEF;
  }
}
