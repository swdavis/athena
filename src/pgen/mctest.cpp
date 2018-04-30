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

  Real temp = pin->GetReal("problem","temp");
  Real rho = pin->GetReal("problem","rho");

  // Set initial conditions
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = rho*temp;
      }
    }
  }

}

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin){

  //EnrollUserEmissionFunction();
  
}

void MonteCarloBlock::InitializePhoton(MeshBlock *pmb, Photon *pphot) {

  // Set status flag

  pphot->status = EVOLVING;

  // Choose random intial position, weights, energy, and direction
  // for photon emission.  In this version an equal number of photons
  // is emitted in  each grid zone.  The relative emission from each grid 
  // zone is then accounted for by a weighting factor cweight. 
  //pmb->block_size.nx1;

  int nx1 = pmb->block_size.nx1;
  int nx2 = pmb->block_size.nx2;
  int nx3 = pmb->block_size.nx3;

  pphot->i1 = static_cast<int>(pran->uniform()*static_cast<Real>(nx1));
  pphot->i2 = static_cast<int>(pran->uniform()*static_cast<Real>(nx2));
  pphot->i3 = static_cast<int>(pran->uniform()*static_cast<Real>(nx3));

  // cweight is a constant weighting factor which accounts for the
  // emissivity of the grid zone in which the photon was emitted
  if (zone_weight_flag)
    pphot->weight = emission(pphot->i3,pphot->i2,pphot->i1);
    //weight = pmb->etat[];
  else
    pphot->weight = 1.0;
  
  std::cout << "test: " << pphot->weight << ' ' << pphot->i2 << ' ' 
            << pphot->i2 << ' ' << pphot->i3 << std::endl;

  // Obtain initial position within zone
  //get_position_uniform(-1,pG,pPack);
   
  // Obtain intitial energy, polarization, direction and weight
  //photon_emit_freefree(ind,pG,pPack);

  // All packet propogations ceases below a certain cutoff weight WMIN
  // We rescale cweight to account for the differing initial weights so that
  // all photons begin with an equivalent weight of unity, and therefore, all
  // packets are treated equivalently */

  if (pphot->weight < 0.0) pphot->status = DESTROYED;

  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  //abs_coef = absopac(pG->temp[ind],pG->dens[ind],pPack->energy);
  //sct_coef = sctopac(pG->temp[ind],pG->dens[ind],pPack->energy);
}

