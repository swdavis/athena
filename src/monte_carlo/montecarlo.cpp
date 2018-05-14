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
#include "../hydro/hydro.hpp"

// constructor, initializes data structures and parameters

MonteCarlo::MonteCarlo(ParameterInput *pin, Mesh *pmesh) {

  MonteCarloBlock *pfirst;
 
  pmy_mesh = pmesh;

  pmcout = new MCOutput(this,pin);

  InitEmission=NULL;
  GetTemperature=NULL;

  // Set flags that control emission, absorption and scattering
  emission_meth = GetEmissionFlag(pin->GetOrAddString("montecarlo","emission","error"));
  if (emission_meth ==  EMISFF) {
    InitEmission = InitializeEmissionFreeFree;
  }
  absorption_meth = GetAbsorptionFlag(pin->GetOrAddString("montecarlo","absorption",
                                                          "error"));
  scattering_meth = GetScatteringFlag(pin->GetOrAddString("montecarlo","scattering",
                                                          "error"));
   // read bc flags for each of the 6 boundaries.
  mc_bcs[INNER_X1] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ix1_mc_bc","escape"));
  mc_bcs[OUTER_X1] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ox1_mc_bc","escape"));
  mc_bcs[INNER_X2] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ix2_mc_bc","escape"));
  mc_bcs[OUTER_X2] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ox2_mc_bc","escape"));
  mc_bcs[INNER_X3] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ix3_mc_bc","escape"));
  mc_bcs[OUTER_X3] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ox3_mc_bc","escape"));

  moments_flag = pin->GetOrAddBoolean("montecarlo","moments",true);
  lorentz_trans_flag = pin->GetOrAddBoolean("montecarlo","lorentz_trans",true);

  // Create and intitialize randon number generator
  iseed = pin->GetInteger("montecarlo","iseed");

  // initialize monte carlo block structure to match mesh
  MeshBlock *pmb = pmesh->pblock;
  pblock = new MonteCarloBlock(pmb, this, pin);
  pfirst = pblock;
  pmb=pmb->next;
  while (pmb != NULL)  {
    pblock->next = new MonteCarloBlock(pmb, this, pin);
    pblock->next->prev = pblock;
    pblock = pblock->next;
    pmb=pmb->next;
  }
  pblock = pfirst;

}

// destructor

MonteCarlo::~MonteCarlo() {

  delete pmcout;

  while(pblock->next != NULL)
    delete pblock->next;
  delete pblock;

}

//----------------------------------------------------------------------------------------
//! \fn enum AbsorptionFlag GetAbsorptionFlag(std::string input_string)
//  \brief set absorption flag

enum AbsorptionFlag GetAbsorptionFlag(std::string input_string) {
  if (input_string == "user") {
    return ABSUSER;
  } else if (input_string == "none") {
    return ABSNONE;
  } else if (input_string == "freefree") {
    return ABSFF;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetAbsorptionFlag" << std::endl
        << "Input string=" << input_string << " not valid absorption type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}

//----------------------------------------------------------------------------------------
//! \fn enum ScatteringFlag GetScatteringFlag(std::string input_string)
//  \brief set scatering flag

enum ScatteringFlag GetScatteringFlag(std::string input_string) {
  if (input_string == "user") {
    return SCATUSER;
  } else if (input_string == "none") {
    return SCATNONE;
  } else if (input_string == "isotropic") {
    return SCATISO;
  } else if (input_string == "thomson") {
    return SCATTHOM;
  } else if (input_string == "compton") {
    return SCATCOMP;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetAbsorptionFlag" << std::endl
        << "Input string=" << input_string << " not valid scattering type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}

//----------------------------------------------------------------------------------------
//! \fn enum EmissionFlag GetEmissionFlag(std::string input_string)
//  \brief set emission flag

enum EmissionFlag GetEmissionFlag(std::string input_string) {
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
//! \fn enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string)
//  \brief set boundary flag

enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string) {

  if (input_string == "periodic") {
    return MC_PERIODIC_BNDRY;
  } else if (input_string == "escape") {
    return MC_ESCAPE_BNDRY;
  } else if (input_string == "absorb") {
    return MC_ABSORB_BNDRY;
  } else if (input_string == "polar") {
    return MC_POLAR_BNDRY;
  } else if (input_string == "reflect") {
    return MC_REFLECT_BNDRY;
  } else if (input_string == "user") {
    return MC_USER_BNDRY;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetMCBoundaryFlag" << std::endl
        << "Input string=" << input_string << " not valid boundary type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}
//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetDensity(MonteCarloBlock *pmcb)
//  \brief Make hard copy of density from MeshBlock to MonteCarloBlock. Uses hard copy
//  so that rho is always in cgs units

void MonteCarlo::GetDensity(MonteCarloBlock *pmcb) {

  
  // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        pmcb->rho(k,j,i) = pmcb->codetocgs_rho * pmcb->pmy_block->phydro->u(IDN,k,j,i);
      }}}
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetVelocities(MonteCarloBlock *pmcb)
//  \brief Make hard copy of velocites from MeshBlock to MonteCarloBlock. Uses hard copy
//  so that velocities is always fraction of speed of light

void MonteCarlo::GetVelocity(MonteCarloBlock *pmcb) {

  
  // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        Real rho = pmcb->pmy_block->phydro->u(IDN,k,j,i);
        pmcb->vel(0,k,j,i) = pmcb->codetoc_vel*pmcb->pmy_block->phydro->u(IM1,k,j,i)/rho;
        pmcb->vel(1,k,j,i) = pmcb->codetoc_vel*pmcb->pmy_block->phydro->u(IM2,k,j,i)/rho;
        pmcb->vel(2,k,j,i) = pmcb->codetoc_vel*pmcb->pmy_block->phydro->u(IM3,k,j,i)/rho;     
        // transform to cartesian if not cartesian
      }}}
}

//----------------------------------------------------------------------------------------
//! \fn void DefaultGetTemperature(MonteCarloBlock *pmcb)
//  \brief default function for computing temperature if no user function provided.
//  Assumes that code values correspond to cgs with simple equation of state.

void DefaultGetTemperature(MonteCarloBlock *pmcb) {

  Real rideal = 8.314e7;
  Hydro* phydro = pmcb->pmy_block->phydro;

   // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        pmcb->tgas(k,j,i) = phydro->w(IEN,k,j,i)/phydro->w(IDN,k,j,i)/rideal;

      }}}

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserEmissionInitialization(EmisFunc_t emissfunc)
//  \brief Enroll a user-defined function for initializing emission methods

void MonteCarlo::EnrollUserEmissionInitialization(EmisFunc_t emissfunc) {
  InitEmission = emissfunc;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserGetTemperature(TempFunc_t tempfunc)
//  \brief Enroll a user-defined function for computing temperature

void MonteCarlo::EnrollUserGetTemperature(TempFunc_t tempfunc) {
  GetTemperature = tempfunc;
}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::LaunchPhotons()
//  \brief start evolving photons in each monte carlo block
//
//  Temporary implementation

void MonteCarlo::LaunchPhotons() {
 
  MonteCarloBlock *pmcb = pblock;

  if (InitEmission == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in LaunchPhotons()" << std::endl
        << "InitEmission function pointer not set." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (GetTemperature == NULL)
    GetTemperature = DefaultGetTemperature;

  
  // Initialize variables over all blocks
  GetDensity(pmcb);
  GetTemperature(pmcb);
  //(pmcb->*(pmcb->GetTemperature2))();
  if (lorentz_trans_flag) GetVelocity(pmcb);
  InitEmission(pmcb);
  pmcb = pmcb->next;
  while (pmcb != NULL) {
    GetDensity(pmcb);
    GetTemperature(pmcb);
    if (lorentz_trans_flag) GetVelocity(pmcb);
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
