//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//  \brief implementation of functions in class MonteCarlo, MCRandom

#include <gsl/gsl_randist.h>
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "montecarlo.hpp"

#include "../globals.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"


// constructor, initializes data structures and parameters

MonteCarlo::MonteCarlo(ParameterInput *pin, Mesh *pmesh) {

  MonteCarloBlock *pfirst;
 
  pmy_mesh = pmesh;

  emission_meth = GetEmissionFlag(pin->GetOrAddString("montecarlo","emission","error"));
  if (emission_meth ==  EMISFF) {
    InitEmission = InitializeEmissionFreeFree;
  }

  // initialize monte carlo block structure to match mesh
  MeshBlock *pmb = pmesh->pblock;
  pblock = new MonteCarloBlock(pmb, this, pin);
  pfirst = pblock;
  //pblock->pmy_mc = this; // link block to this mc
 
  pmb=pmb->next;
  while (pmb != NULL)  {
    pblock->next = new MonteCarloBlock(pmb, this, pin);
    pblock->next->prev = pblock;
    pblock = pblock->next;
    //pblock->pmy_mc = this;
    pmb=pmb->next;
  }
  pblock = pfirst;

}

// destructor

MonteCarlo::~MonteCarlo() {

  while(pblock->next != NULL)
    delete pblock->next;
  delete pblock;

}

enum EmissionFlag MonteCarlo::GetEmissionFlag(std::string input_string) {
  if (input_string == "user") {
    return EMISUSER;
  } else if (input_string == "freefree") {
    return EMISFF;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetEmissionFlag" << std::endl
        << "Input string=" << input_string << " not valid emission type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}

//----------------------------------------------------------------------------------------
//! \fn oid MonteCarlo::EnrollUserEmissionInitialization(EmisFunc_t emissfunc)
//  \brief Enroll a user-defined function for initializing emission methods

void MonteCarlo::EnrollUserEmissionInitialization(EmisFunc_t emissfunc) {
  InitEmission = emissfunc;
}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::LaunchPhotons()
//  \brief start evolving photons in each monte carlo block

void MonteCarlo::LaunchPhotons() {
 
  MonteCarloBlock *pmcb = pblock;

  // Looop over all MeshBlocks
  //int nmb = pmesh->GetNumMeshBlocksThisRank(Globals::my_rank);
  //pmb = pmesh->pblock;
  //for (int i=0; i<nmb; ++i) {
  //  pmb->pmcb->InitEmission(pmb->pmcb,pmb);
  //  pmb=pmb->next;
  //}

  // Initialize emission over all blocks
  //((pmcb).*(InitEmission))();
  InitEmission(pmcb);
  pmcb = pmcb->next;
  while (pmcb != NULL) {
    //pmcb->InitEmission();
    InitEmission(pmcb);
    pmcb = pmcb->next;
  }

  // transfer photons overall blocks
  pmcb = pblock;
  pmcb->TransferPhotons();
  pmcb = pmcb->next;
  while (pmcb != NULL) {
    pmcb->TransferPhotons();
    pmcb = pmcb->next;
  }

  return;
}


// constructor

MCRandom::MCRandom(int iseed) {
  dev = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(dev, iseed);
}

// destructor

MCRandom::~MCRandom() {

}

Real MCRandom::uniform() {

  return static_cast<Real>(gsl_rng_uniform(dev));
}
