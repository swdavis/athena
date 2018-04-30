//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photonmover.cpp
//  \brief implementation for photon moving functions

// Athena++ headers
#include "photon.hpp"
#include "photonmover.hpp"
#include "../mesh/mesh.hpp"

// Implementation of base class

PhotonMover::PhotonMover(MonteCarloBlock *pmcb) {

  pmy_mcb = pmcb;
}

PhotonMover::~PhotonMover() {

}

void PhotonMover::Move(MeshBlock *pmb, Photon *pphot) {

}


// Implementation of Cartesian Photon mover

CartesianMover::CartesianMover(MonteCarloBlock *pmcb) 
  : PhotonMover(pmcb) {

}

CartesianMover::~CartesianMover() {

}

void CartesianMover::Move(MeshBlock *pmb, Photon *pphot) {

}
