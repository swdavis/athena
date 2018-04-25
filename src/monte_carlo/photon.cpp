//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//  \brief implementation of functions in class Photon

// Athena++ headers
#include "photon.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// constructor, initializes data structures and parameters

Photon::Photon(MonteCarlo *pmc) {

  pmy_mc = pmc;

}

// destructor

Photon::~Photon() {

}
