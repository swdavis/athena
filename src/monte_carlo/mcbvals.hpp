#ifndef MCBVALS_HPP
#define MCBVALS_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mcbvals.hpp
//  \brief defines MCBoundaryValues prototypes for boundary condtion functions
//
// Current design focusses on implementing static post-processing so this implementation
// will evolve

// Athena++ classes headers
#include "../athena.hpp"
#include "montecarlo.hpp"

class Photon;
class MCCoord;

//----------------------------------------------------------------------------------------
// function pointer prototypes for boundary conditions set at runtime
typedef void (*MCBValFunc_t)(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);

//---------------------- prototypes for boundary functions -------------------------------
void PeriodicInnerX1(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
void PeriodicOuterX1(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
void PeriodicInnerX2(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
void PeriodicOuterX2(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
void PeriodicInnerX3(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
void PeriodicOuterX3(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
//void PeriodicWedgeInnerX3(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
//void PeriodicWedgeOuterX3(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
void Escape(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
void Absorb(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
void Polar(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
//----------------------------------------------------------------------------------------
//! \class MCBoundaryValues
//  \brief BVals data and functions for monte carlo

class MCBoundaryValues {
public:
  MCBoundaryValues(MonteCarloBlock *pmcb, ParameterInput *pin);
  ~MCBoundaryValues();

  MonteCarloBlock *pmy_mcb;

  MCBValFunc_t BoundaryFunction_[6];
  
};

#endif // MCBVALS_HPP
