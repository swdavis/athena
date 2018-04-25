//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo_tasks.cpp
//  \brief mesh level monte carlo control

// needed for vector of pointers in D
#include <vector>

// Athena++ headers
#include "monte_carlo_tasks.hpp"
#include "monte_carlo.hpp"
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"

// constructor, initializes data structures and parameters

MonteCarloTasks::MonteCarloTasks(Mesh *pm, ParameterInput *pin) {
  
  pmy_mesh_ = pm;

}

// destructor

MonteCarloTasks::~MonteCarloTasks() {

}

void MonteCarloTasks::LaunchPhotons(Mesh *pmesh) {
  MeshBlock *pmb = pmesh->pblock;

  // Looop over all MeshBlocks
  int nmb = pmesh->GetNumMeshBlocksThisRank(Globals::my_rank);

  // Looop over all MeshBlocks
  pmb = pmesh->pblock;
  for (int i=0; i<nmb; ++i) {
    pmb->pmc->TransferPhotons();
    pmb=pmb->next;
  }

  return;
}
