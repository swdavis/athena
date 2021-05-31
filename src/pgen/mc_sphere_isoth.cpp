//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mctest.cpp
//  \brief Problem generator for  monte carlo through uniform sphere
//
//========================================================================================

// C/C++ headers
#include <iostream> 
#include <stdexcept>

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
#include "../globals.hpp"

#if !MONTE_CARLO_STATIC
#error "This problem requires monte carlo"
#endif

namespace {
  // Global variables
  Real rad0,path0;
  Real energy0;
  bool srcdist;
  int i1start,i2start,i3start;
  Real logemin, logemax;

  // function headers
  void SphericalEscape(MonteCarloBlock *pmcb, Photon *phot, PhotonMover *pmover);
  void TimedEscape(MonteCarloBlock *pmcb, Photon *phot, PhotonMover *pmover);
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief monte carlo test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real rideal = 8.314e7;
  Real c = 2.99792458e10;
  Real temp = pin->GetReal("problem","temp");
  Real tau = pin->GetReal("problem","tau");
  Real rad0 = pin->GetReal("problem","radius");
  Real vel = pin->GetOrAddReal("problem","velocity",0.);
  Real gamma = peos->GetGamma();
  vel *= c;

  Real heabund = 0.09; //hardcode for now (should be parameter)
  //Real heabund = 0.;
  Real mp = 1.6726e-24; 
  Real sigmat = 6.65248e-25;
  Real kappaes = sigmat * (1. + 2.*heabund) / (mp * (1.+4.*heabund) );

  Real rho = tau / (kappaes * rad0);
  //printf("rho: %g %g %g %g\n",rho,kappaes,rad0,tau);
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

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){

  // Initialize user variable for time 
  nuser_var = 1;

  // If path is set in problem generator, terminate photon integration based on time
  // but if not terminated based on radius
  Real path = pin->GetOrAddReal("problem","path",-1.);
  if (path > 0.) {
    EnrollUserWorkInMove(TimedEscape);
  } else
    EnrollUserWorkInMove(SphericalEscape);

}

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {

  if (emission_meth == EMISNONE) {
    Real x0 = pin->GetReal("problem","x0");
    Real temp = pin->GetReal("problem","temp");
    Real kb = 1.380649e-16;
    energy0 = kb*temp*x0;
  } else if (emission_meth == EMISFF) {
    // Set the energy boundaries for free-free emission
    Real everg = 1.6021772e-12;
    logemin = log(everg*pin->GetReal("problem", "emin"));
    logemax = log(everg*pin->GetReal("problem", "emax"));
  }
  // Set variables 
  srcdist =pin->GetOrAddBoolean("problem","srcdist",false);
  rad0 = pin->GetReal("problem","radius");
  path0 = pin->GetOrAddReal("problem","path",-1.);

  // Deterime cell of initial photon, which is asssumed to include 
  // if origin if more than one cell is specified for each direction
  i1start = -1;
  for(int i=is; i<=ie; i++) {
    if ((0. > pcoord->x1f(i)) && (0. <= pcoord->x1f(i+1)))
      i1start = i;
  }
  i2start = -1;
  for(int i=js; i<=je; i++) {
    if ((0. > pcoord->x2f(i)) && (0. <= pcoord->x2f(i+1)))
      i2start = i;
  }
  i3start = -1;
  for(int i=ks; i<=ke; i++) {
    if ((0. > pcoord->x3f(i)) && (0. <= pcoord->x3f(i+1)))
      i3start = i;
  }
  if ((i1start < 0) || (i2start < 0) || (i3start < 0)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MonteCarloProblemGenerator" << std::endl
        << "Origin not found within domain." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}


void MonteCarloBlock::InitializePhoton(Photon *pphot) {

  // Set status flag
  pphot->status = EVOLVING;

  pphot->user_var[0] = 0.;

  // Initialize cell number
  pphot->i1 = i1start;
  pphot->i2 = i2start;
  pphot->i3 = i3start;

  // Initialize Photon weights, energy, direction, polarization
  if (pmy_mc->emission_meth == EMISNONE) {
    pphot->weight = 1.0;
    pphot->energy = energy0;
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

    // Get initial position of photon
    if (srcdist) {
      // Model an exponential weight time distribuion
      Real x = pran->uniform();
      while (x <= 0.)
        x = pran->uniform();
      Real dev = pran->uniform();
      while (sin(PI*x)*x < 0.57923*dev) {
      // Yields nearly exponential escape time distribution
      //while (sin(PI*x)/(x*PI) < dev) {
        x = pran->uniform();
        while (x <= 0.)
          x = pran->uniform();
        dev = pran->uniform();
      }
      
      Real r0 = x*rad0;
      Real phi = 2. * PI * pran->uniform();
      Real cphi = cos(phi);
      Real sphi = sin(phi);
      
      Real cth = 2. * pran->uniform() - 1.;
      Real sth = sqrt(1. - SQR(cth));
      pphot->x[0] = r0*sth*cphi;
      pphot->x[1] = r0*sth*sphi;
      pphot->x[2] = r0*cth;
    } else {
      // Initialize photon at the origin
      pphot->x[0] = 0.;
      pphot->x[1] = 0.;
      pphot->x[2] = 0.;
    }

  } else if (pmy_mc->emission_meth == EMISFF) {
    pphot->weight = emission(pphot->i3,pphot->i2,pphot->i1);
    PhotonEmitFreeFree(this,pphot,logemin,logemax);
    Real r0 = pow(pran->uniform()*rad0*rad0*rad0,1./3.);
    Real phi = 2. * PI * pran->uniform();
    Real cphi = cos(phi);
    Real sphi = sin(phi);      
    Real cth = 2. * pran->uniform() - 1.;
    Real sth = sqrt(1. - SQR(cth));
    pphot->x[0] = r0*sth*cphi;
    pphot->x[1] = r0*sth*sphi;
    pphot->x[2] = r0*cth;
  }


  
  for (int i=0; i<pphot->nuser_var; i++)
    pphot->user_var[i] = 0.;

  if (pphot->weight < 0.0) pphot->status = DESTROYED;
  
  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = AbsorptionOpacity(this,pphot);
  pphot->sct_coef = ScatteringOpacity(this,pphot);

}

namespace {

// Used to evalue photons time (path) length distribution as fixed spherical
// escape surface
void SphericalEscape(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover) {

  // Update path length for user output
  pphot->user_var[0] += pmover->dl;
  // First check radius condition
  Real r = sqrt(SQR(pphot->x[0])+SQR(pphot->x[1])+SQR(pphot->x[2]));
 
  if (r >= rad0) {
    Real dr = r-rad0;
    pphot->user_var[0] -= dr;
    for (int i=0; i<3; i++) {
      // assumes cartesian for now
      pphot->x[i] -= pphot->k[i]*dr;
    }
    pphot->status = ESCAPED;
    pphot->face = FACE_UNDEF;
  } 
 
}

// Used to test photons radial distributions after a fixed travel time
void TimedEscape(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover) {

  // Update path length for user output
  pphot->user_var[0] += pmover->dl;

  // First check radius condition
  Real r = sqrt(SQR(pphot->x[0])+SQR(pphot->x[1])+SQR(pphot->x[2]));
  if (r >= rad0) {
    //Real dr = r-rad0;
    //pphot->user_var[0] -= dr;
    //for (int i=0; i<3; i++) {
    //  // assume cartesian for now
    //  pphot->x[i] -= pphot->k[i]*dr;
    //}
    pphot->status = ESCAPED;
    pphot->face = FACE_UNDEF;
  }
  // Then check path condition -- ensures path is not over estimated
  if (pphot->user_var[0] >= path0) {
    //printf("%g %g %g %g\n",pphot->x[0],pphot->x[1],pphot->x[2],pphot->user_var[0]);
    Real dp = pphot->user_var[0] - path0;
    pphot->user_var[0] = path0;
    for (int i=0; i<3; i++) {
      // assume cartesian for now
      pphot->x[i] -= pphot->k[i]*dp;
    }
    pphot->status = ESCAPED;
    pphot->face = FACE_UNDEF;
  }
}

} //namespace
