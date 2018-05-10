//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file montecarloblock.cpp
//  \brief implementation of functions in class MCBoundaryValues

#include <stdexcept>  // runtime_error

// Athena++ headers
#include "mcbvals.hpp"

// constructor, initializes data structures and parameters

MCBoundaryValues::MCBoundaryValues(MonteCarloBlock *pmcb, ParameterInput *pin) {
  
  pmy_mcb = pmcb;

// Set BC functions for each of the 6 boundaries -----------------------------------------
  // Inner x1
  switch(pmcb->mcb_bcs[INNER_X1]) {
    case MC_PERIODIC_BNDRY:
      BoundaryFunction_[INNER_X1] = PeriodicInnerX1;
      break;
    case MC_ESCAPE_BNDRY:
      BoundaryFunction_[INNER_X1] = Escape;
      break;
    case MC_ABSORB_BNDRY:
      BoundaryFunction_[INNER_X1] = Absorb;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MCBoundaryValues constructor" << std::endl
          << "Flag ix1_mc_bc=" << pmcb->mcb_bcs[INNER_X1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  // outer x1
  switch(pmcb->mcb_bcs[OUTER_X1]) {
    case MC_PERIODIC_BNDRY:
      BoundaryFunction_[OUTER_X1] = PeriodicOuterX1;
      break;
    case MC_ESCAPE_BNDRY:
      BoundaryFunction_[OUTER_X1] = Escape;
      break;
    case MC_ABSORB_BNDRY:
      BoundaryFunction_[OUTER_X1] = Absorb;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MCBoundaryValues constructor" << std::endl
          << "Flag ix1_mc_bc=" << pmcb->mcb_bcs[OUTER_X1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  // Inner x2
  switch(pmcb->mcb_bcs[INNER_X2]) {
    case MC_PERIODIC_BNDRY:
      BoundaryFunction_[INNER_X2] = PeriodicInnerX2;
      break;
    case MC_ESCAPE_BNDRY:
      BoundaryFunction_[INNER_X2] = Escape;
      break;
    case MC_ABSORB_BNDRY:
      BoundaryFunction_[INNER_X2] = Absorb;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MCBoundaryValues constructor" << std::endl
          << "Flag ix1_mc_bc=" << pmcb->mcb_bcs[INNER_X2] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  // outer x2
  switch(pmcb->mcb_bcs[OUTER_X2]) {
    case MC_PERIODIC_BNDRY:
      BoundaryFunction_[OUTER_X2] = PeriodicOuterX2;
      break;
    case MC_ESCAPE_BNDRY:
      BoundaryFunction_[OUTER_X2] = Escape;
      break;
    case MC_ABSORB_BNDRY:
      BoundaryFunction_[OUTER_X2] = Absorb;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MCBoundaryValues constructor" << std::endl
          << "Flag ix1_mc_bc=" << pmcb->mcb_bcs[OUTER_X2] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  // Inner x3
  switch(pmcb->mcb_bcs[INNER_X3]) {
    case MC_PERIODIC_BNDRY:
      BoundaryFunction_[INNER_X3] = PeriodicInnerX3;
      break;
    case MC_ESCAPE_BNDRY:
      BoundaryFunction_[INNER_X3] = Escape;
      break;
    case MC_ABSORB_BNDRY:
      BoundaryFunction_[INNER_X3] = Absorb;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MCBoundaryValues constructor" << std::endl
          << "Flag ix1_mc_bc=" << pmcb->mcb_bcs[INNER_X3] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  // outer x3
  switch(pmcb->mcb_bcs[OUTER_X3]) {
    case MC_PERIODIC_BNDRY:
      BoundaryFunction_[OUTER_X3] = PeriodicOuterX3;
      break;
    case MC_ESCAPE_BNDRY:
      BoundaryFunction_[OUTER_X3] = Escape;
      break;
    case MC_ABSORB_BNDRY:
      BoundaryFunction_[OUTER_X3] = Absorb;
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in MCBoundaryValues constructor" << std::endl
          << "Flag ix1_mc_bc=" << pmcb->mcb_bcs[OUTER_X3] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

}


// destructor

MCBoundaryValues::~MCBoundaryValues() {

}


// Currently assumes single block -- ***will need to be changed***

//----------------------------------------------------------------------------------------
//! \fn void PeriodicInnerX1(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot)
//  \brief periodic boundary conditions, inner x1 boundary

void PeriodicInnerX1(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot) {

  pphot->i1 = pmcb->ie;
  pphot->x[0] = pco->x1f(pphot->i1+1);

}

//----------------------------------------------------------------------------------------
//! \fn void PeriodicOuterX1(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot)
//  \brief periodic boundary conditions, outer x1 boundary

void PeriodicOuterX1(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot) {

  pphot->i1 = pmcb->is;
  pphot->x[0] = pco->x1f(pphot->i1);

}

//----------------------------------------------------------------------------------------
//! \fn void PeriodicInnerX2(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot)
//  \brief periodic boundary conditions, inner x2 boundary

void PeriodicInnerX2(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot) {

  pphot->i2 = pmcb->je;
  pphot->x[1] = pco->x2f(pphot->i2+1);

}

//----------------------------------------------------------------------------------------
//! \fn void PeriodicOuterX2(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot)
//  \brief periodic boundary conditions, outer x2 boundary

void PeriodicOuterX2(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot) {

  pphot->i2 = pmcb->js;
  pphot->x[1] = pco->x2f(pphot->i2);

}

//----------------------------------------------------------------------------------------
//! \fn void PeriodicInnerX3(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot)
//  \brief periodic boundary conditions, inner x3 boundary

void PeriodicInnerX3(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot) {

  pphot->i3 = pmcb->ke;
  pphot->x[2] = pco->x3f(pphot->i3+1);

}

//----------------------------------------------------------------------------------------
//! \fn void PeriodicOuterX3(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot)
//  \brief periodic boundary conditions, outer x3 boundary

void PeriodicOuterX3(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot) {

  pphot->i3 = pmcb->ks;
  pphot->x[2] = pco->x3f(pphot->i3);

}


//----------------------------------------------------------------------------------------
//! \fn void Escape(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot)
//  \brief mark photon as escaped

void Escape(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot) {

  pphot->status = ESCAPED;

}

//----------------------------------------------------------------------------------------
//! \fn void Absorb(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot)
//  \brief mark photon as escaped

void Absorb(MonteCarloBlock *pmcb, Coordinates *pco, Photon *pphot) {

  pphot->status = DESTROYED;

}
