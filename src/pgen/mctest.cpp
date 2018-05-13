//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mctest.cpp
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

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief monte carlo test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real rideal = 8.314e7;
  Real temp = pin->GetReal("problem","temp");
  //Real rho = pin->GetReal("problem","rho");
  Real taumin = pin->GetReal("problem","taumin");
  Real taumax = pin->GetReal("problem","taumax");
  Real gamma = peos->GetGamma();
  Real zlow = pin->GetReal("mesh","x3min");
  Real zhigh = pin->GetReal("mesh","x3max");
  int nz = pin->GetInteger("mesh","nx3");

  Real dz = (zhigh-zlow) / static_cast<Real>(nz);
  Real step = log10(taumax/taumin) / static_cast<Real>(nz-1);
       
  AthenaArray<Real> tau1d,dens1d;
  
  // Chosen to match initialize.py in mcgrid
  tau1d.NewAthenaArray(nz); 
  dens1d.NewAthenaArray(nz);
  for (int i=0; i<nz; ++i) {
    tau1d(i) = log10(taumin) + step * static_cast<Real>(i);
    tau1d(i) = pow(10.,tau1d(i));
  }
  Real kapes = 0.33;
  dens1d(0) = tau1d(0) / dz / kapes;
  for (int i=1; i<nz; ++i) {
    dens1d(i) = (tau1d(i)-tau1d(i-1) ) / (dz * kapes);
  }

  Real rho;
  // Set initial conditions
  for (int k=ks; k<=ke; k++) {
    rho = dens1d(ke-k);
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = rideal*rho*temp/(gamma-1.0);
      }
    }
  }

}

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){

}

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin){

  // Set codetocgs here
  //EnrollUserEmissionFunction();
  
}

void MonteCarloBlock::InitializePhoton(Photon *pphot) {

  // Set status flag

  pphot->status = EVOLVING;

  // Choose random intial position, weights, energy, and direction
  // for photon emission.  In this version an equal number of photons
  // is emitted in  each grid zone.  The relative emission from each grid 
  // zone is then accounted for by a weighting factor cweight. 

  Real nx1 = static_cast<Real>(ie-is);
  Real nx2 = static_cast<Real>(je-js);
  Real nx3 = static_cast<Real>(ke-ks);

  pphot->i1 = static_cast<int>(pran->uniform()*nx1)+is;
  pphot->i2 = static_cast<int>(pran->uniform()*nx2)+js;
  pphot->i3 = static_cast<int>(pran->uniform()*nx3)+ks;
  //pphot->i3 = 32+ks;
  
  // cweight is a constant weighting factor which accounts for the
  // emissivity of the grid zone in which the photon was emitted
  if (zone_weight_flag) {
    pphot->eweight = emission(pphot->i3,pphot->i2,pphot->i1);
    //pphot->weight = pphot->eweight;
    pphot->weight = 1.0;
  }

  //std::cout << "test: " << pphot->weight << ' ' << pphot->i1 << ' ' 
  //          << pphot->i2 << ' ' << pphot->i3 << std::endl;

  // Obtain initial position within zone
  GetZonePosition(pphot,pran,pmy_coord);
  //std::cout << pphot->x[0] << ' ' << pphot->x[1] << ' ' << pphot->x[2]
  //          << std::endl;

  // Obtain intitial energy, polarization, direction and weight
  // Utilize free-free emission function in emission.cpp
  PhotonEmitFreeFree(this,pphot);


  if (pphot->weight < 0.0) pphot->status = DESTROYED;

  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = AbsorptionOpacity(this,pphot);
  pphot->sct_coef = ScatteringOpacity(this,pphot);
  //std::cout << "Init: " << pphot->i3-ks << ' ' << pphot->abs_coef << ' ' <<  pphot->sct_coef << std::endl;

}

