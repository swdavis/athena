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

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

namespace {
  // Global variables
  static Real rad0,path0;
  static Real energy0;
  static bool first = true;
  static int i1start,i2start,i3start;
}

// function headers
void SphericalEscape(MonteCarloBlock *pmcb, Photon *phot, PhotonMover *pmover);
void TimedEscape(MonteCarloBlock *pmcb, Photon *phot, PhotonMover *pmover);

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief monte carlo test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real rideal = 8.314e7;
  Real c = 2.997924589e10;
  Real temp = pin->GetReal("problem","temp");
  Real tau = pin->GetReal("problem","tau");
  Real rad0 = pin->GetReal("problem","radius");
  Real vel = pin->GetOrAddReal("problem","velocity",0.);
  Real gamma = peos->GetGamma();
  vel *= c;

  Real kb = 1.3806504e-16;
  Real mass = 1.660538782e-24;
  Real vth = sqrt(2.*kb*temp/mass);
  Real nu0 = 2.468e15;
  Real dopw = nu0 * vth / c;
  Real kappa = ResLinePre() / (mass*sqrt(PI)*dopw);
  Real rho = tau / (kappa * rad0);
  printf("kappa: %g %g %g\n",kappa,ResLinePre(),rho/mass);
  printf("voigt: %g\n",XsecVoigt(nu0,temp));
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
  nuser_var = 1;
  if (emission_meth == EMISUSER) {
    Real x0 = pin->GetReal("problem","x0");
    Real temp = pin->GetReal("problem","temp");
    Real kb = 1.3806504e-16;
    Real mass = 1.660538782e-24;
    Real vth = sqrt( 2. * kb * temp / mass);
    Real c = 2.997924589e10;
    Real nu0 = 2.468e15;
    Real dopw = nu0 * vth / c;
    Real lorw = 6.265e8/(4.*PI);
    printf("dop, lor, a: %e %e %e\n",dopw,lorw,lorw/dopw);
    Real h = 6.6260755e-27;
    energy0 = h * (nu0 + dopw * x0);

  }

}

void MonteCarloBlock::InitUserMonteCarloBlockData(ParameterInput *pin){

  // enroll function to cease photon propogation based on escape radius
  // or total integration time
  rad0 = pin->GetReal("problem","radius");
  path0 = pin->GetOrAddReal("problem","path",-1.);
  if (path0 > 0.) {
    EnrollUserWorkInMove(TimedEscape);
  } else
    EnrollUserWorkInMove(SphericalEscape);

}

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin){
  
}

void MonteCarloBlock::InitializePhoton(Photon *pphot) {


  // Set status flag
  pphot->status = EVOLVING;

  if (first) {
    // Deterime zone of initial photon -- asssumed to be zone that includes origin
  
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
      msg << "### FATAL ERROR in InitUserMonteCarloBlockData" << std::endl
          << "Origin not found within domain." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
  pphot->i1 = i1start;
  pphot->i2 = i2start;
  pphot->i3 = i3start;

  // Initialize Photon weights, energy, direction, polarization
  if (pmy_mc->emission_meth == EMISUSER) {
    pphot->eweight = 1.0;
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
 
    // Initialize photon at the origin
    pphot->x[0] = 0.;
    pphot->x[1] = 0.;
    pphot->x[2] = 0.;
    

  }

  
  for (int i=0; i<pphot->nuser_var; i++)
    pphot->user_var[i] = 0.;

  if (pphot->weight < 0.0) pphot->status = DESTROYED;
  
  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = AbsorptionOpacity(this,pphot);
  pphot->sct_coef = ScatteringOpacity(this,pphot);

  if (first) {
    if ((Globals::my_rank == 0) || (Globals::my_rank == 1))
      printf("taus, taua: %g %g\n",pphot->sct_coef*rad0,pphot->abs_coef*rad0);
    first = false;
  }
}

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
      // assume cartesian for now
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
    Real dr = r-rad0;
    pphot->user_var[0] -= dr;
    for (int i=0; i<3; i++) {
      // assume cartesian for now
      pphot->x[i] -= pphot->k[i]*dr;
    }
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
