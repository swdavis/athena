#ifndef MCDEBUG_HPP
#define MCDEBUG_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file debug.hpp
//  \brief temporary functions for debugging and testing during initial development

// Athena++ classes headers
#include "../athena.hpp"
#include "montecarlo.hpp"

void FinalPositionCartesian(MonteCarloBlock *pmb, Coordinates *pco, Photon *pphot,
                            Real &xf, Real &yf, Real &zf, Real &dl);

#endif
