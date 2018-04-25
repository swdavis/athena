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
#include "../monte_carlo/monte_carlo.hpp"
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

void Photon::InitializePhoton(MeshBlock *pmb) {

  // Set status flag

  status = EVOLVING;

  // Choose random intial position, weights, energy, and direction
  // for photon emission.  In this version an equal number of photons
  // is emitted in  each grid zone.  The relative emission from each grid 
  // zone is then accounted for by a weighting factor cweight. 
  //pmb->block_size.nx1;

  int nx1 = pmb->block_size.nx1;
  int nx2 = pmb->block_size.nx2;
  int nx3 = pmb->block_size.nx3;

  izone[0] = static_cast<int>(pmy_mc->pran->uniform()*static_cast<Real>(nx1));
  izone[1] = static_cast<int>(pmy_mc->pran->uniform()*static_cast<Real>(nx2));
  izone[2] = static_cast<int>(pmy_mc->pran->uniform()*static_cast<Real>(nx3));

  //initialize_zone_prob(pmb,pPack);

  // cweight is a constant weighting factor which accounts for the
  // emissivity of the grid zone in which the photon was emitted
  if (pmy_mc->zone_weight)
    weight = pmy_mc->pran->uniform();
    //weight = pmb->etat[];
  else
    weight = 1.0;
  
  std::cout << "test: " << weight << ' ' << izone[0] << ' ' 
            << izone[1] << ' ' << izone[2] << std::endl;

  // Obtain initial position within zone
  //get_position_uniform(-1,pG,pPack);
   
  // Obtain intitial energy, polarization, direction and weight
  //photon_emit_freefree(ind,pG,pPack);

  // All packet propogations ceases below a certain cutoff weight WMIN
  // We rescale cweight to account for the differing initial weights so that
  // all photons begin with an equivalent weight of unity, and therefore, all
  // packets are treated equivalently */

  if (weight < 0.0) status = DESTROYED;

  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  //abs_coef = absopac(pG->temp[ind],pG->dens[ind],pPack->energy);
  //sct_coef = sctopac(pG->temp[ind],pG->dens[ind],pPack->energy);
}

