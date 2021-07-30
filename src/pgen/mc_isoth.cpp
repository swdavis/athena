//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mc_isoth.cpp
//  \brief Problem generator for monte carlo isothermal atmosphere 
//
//========================================================================================

// C++ headers
#include <iostream> // SWD: temporary for testing

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

#if !MONTE_CARLO_STATIC
#error "This problem requires monte carlo"
#endif

namespace {
  // Global variables
  bool tnorm;
  Real logemin, logemax;
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief monte carlo test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Determine density via optical depth or constant density
  bool constdens = pin->GetOrAddBoolean("problem","constdens",false);
  Real rho, tau, taumin, taumax;
  if (constdens) {
    rho = pin->GetOrAddReal("problem","dens",-1.);
    tau = pin->GetOrAddReal("problem","tau",-1.);
  } else {
    taumin = pin->GetReal("problem","taumin");
    taumax = pin->GetReal("problem","taumax");
  }

  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.6726e-24; 
  Real sigmat = 6.65248e-25;
  Real kappaes = sigmat * (1. + 2.*heabund) / (mp * (1.+4.*heabund) );
  if (tau > 0.) {
    Real length;
    if (COORDINATE_SYSTEM == "cartesian") {
      Real xlow = pin->GetReal("mesh","x3min");
      Real xhigh = pin->GetReal("mesh","x3max");
      length = xhigh-xlow;
    } else {
      length = pin->GetReal("mesh","x1max");
    }
    rho = tau / (kappaes * length);
  }

  AthenaArray<Real> tau1d,dens1d;
  if (!constdens) {
    Real xlow, xhigh;
    int nx;
    if (COORDINATE_SYSTEM == "cartesian") {
      xlow = pin->GetReal("mesh","x3min");
      xhigh = pin->GetReal("mesh","x3max");
      nx = pin->GetInteger("mesh","nx3");
    } else {
      bool radial = pin->GetOrAddBoolean("problem","radial","true");
      if (radial) {
        xlow = pin->GetReal("mesh","x1min");
        xhigh = pin->GetReal("mesh","x1max");
        nx = pin->GetInteger("mesh","nx1");
      }
    }
    tau1d.NewAthenaArray(nx); 
    dens1d.NewAthenaArray(nx);
    Real dx = (xhigh-xlow) / static_cast<Real>(nx);
    Real step = log10(taumax/taumin) / static_cast<Real>(nx-1);
    for (int i=0; i<nx; ++i) {
      tau1d(i) = log10(taumin) + step * static_cast<Real>(i);
      tau1d(i) = pow(10.,tau1d(i));
    }
    dens1d(0) = tau1d(0) / dx / kappaes;
    for (int i=1; i<nx; ++i) {
      dens1d(i) = (tau1d(i)-tau1d(i-1) ) / (dx * kappaes);
    }
  }

  // Assume constant velocity provided as fraction of speed of light
  Real vel = pin->GetOrAddReal("problem","velocity",0.);
  Real c = 2.99792458e10;
  vel *= c;

  // Assume constant temperature and ideal gas
  Real gamma = peos->GetGamma();
  Real rideal = 8.314e7;
  Real tgas = pin->GetReal("problem","temp");
  // Set initial conditions
  if (COORDINATE_SYSTEM == "cartesian") {
    // density varies in the z direction
    for (int k=ks; k<=ke; k++) {
      if (!constdens) rho = dens1d(ke-k);
      for (int j=js; j<=je; j++) {
	for (int i=is; i<=ie; i++) {
	  phydro->u(IDN,k,j,i) = rho;
	  phydro->u(IM1,k,j,i) = 0.0;
	  phydro->u(IM2,k,j,i) = 0.0;
	  phydro->u(IM3,k,j,i) = rho*vel;
	  phydro->u(IEN,k,j,i) = rideal*rho*tgas/(gamma-1.0);
	}
      }
    }
  } else if  (COORDINATE_SYSTEM == "spherical_polar") {
    bool radial = pin->GetOrAddBoolean("problem","radial","true");

    if (radial) {
      // density varies in the r direction
      for (int k=ks; k<=ke; k++) {
	for (int j=js; j<=je; j++) {
	  for (int i=is; i<=ie; i++) {
	    if (!constdens) rho = dens1d(ie-i);
	    phydro->u(IDN,k,j,i) = rho;
	    phydro->u(IM1,k,j,i) = rho*vel;
	    phydro->u(IM2,k,j,i) = 0.0;
	    phydro->u(IM3,k,j,i) = 0.0;
	    phydro->u(IEN,k,j,i) = rideal*rho*tgas/(gamma-1.0);
	  }
	}
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


void MonteCarloBlock::InitializePhoton(Photon *pphot) {

  // Set status flag

  pphot->status = EVOLVING;

  // Choose a random cell for emission
  Real nx1 = static_cast<Real>(ie-is+1);
  Real nx2 = static_cast<Real>(je-js+1);
  Real nx3 = static_cast<Real>(ke-ks+1);

  pphot->i1 = static_cast<int>(pran->uniform()*nx1)+is;
  pphot->i2 = static_cast<int>(pran->uniform()*nx2)+js;
  pphot->i3 = static_cast<int>(pran->uniform()*nx3)+ks;

  // Set weight according to the emission array, which is the relative number of photons
  // per unit time emitted in each cell
  pphot->weight = emission(pphot->i3,pphot->i2,pphot->i1);

  // Obtain initial position within zone
  GetZonePosition(pphot,pran,pcoord);

  // Obtain intitial energy, polarization, direction and weight
  // Utilize free-free emission function in emission.cpp
  if(tnorm) {
    Real logtg = log(tgas(pphot->i3,pphot->i2,pphot->i1));
    PhotonEmitFreeFree(this,pphot,logemin+logtg,logemax+logtg);
  } else{
    PhotonEmitFreeFree(this,pphot,logemin,logemax);
  }
  //  PhotonEmitFreeFree(this,pphot,logemin,logemax);
  // Convert k unit vector to k^\alpha
  if (general_mover_flag) {
    pphot->k[IMC0] = 1.;
    pphot->k[IMC2] /= pphot->x[IMC1];
    pphot->k[IMC3] /= (pphot->x[IMC1]*sin(pphot->x[IMC2]));
    for(int i=0; i<4; i++) pphot->dk[i] = 0.;
  }
  if (pphot->weight < 0.0) pphot->status = DESTROYED;

  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = AbsorptionOpacity(this,pphot);
  pphot->sct_coef = ScatteringOpacity(this,pphot);
 
}

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {


  // Set the energy boundaries for free-free emission
  tnorm = pin->GetOrAddBoolean("problem","tnorm",false);
  if (tnorm) {
    // interpret as xmin/xmax with x=E/(kb*T)
    Real kb = 1.380649e-16;
    logemin = log(kb*pin->GetReal("problem", "emin"));
    logemax = log(kb*pin->GetReal("problem", "emax"));
  } else {
    // interpret as emin/emax in eV
    Real everg = 1.6021772e-12;
    logemin = log(everg*pin->GetReal("problem", "emin"));
    logemax = log(everg*pin->GetReal("problem", "emax"));
  }
  //Real everg = 1.6021772e-12;
  //logemin = log(everg*pin->GetReal("problem", "emin"));
  //logemax = log(everg*pin->GetReal("problem", "emax"));

}

