#ifndef PHOTONMOVER_HPP
#define PHOTONMOVER_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photonmover.hpp
//  \brief defines abstract base derived classes for moving photons

// Athena++ classes headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "photon.hpp"
#include "montecarlo.hpp"

class MeshBlock;
class ParameterInput;
class MonteCarloBlock;
class Photon;

//! \class PhotonMover
//  \brief abstract base class for all derived classes

class PhotonMover {
public:
  PhotonMover(MonteCarloBlock *pmcb);
  ~PhotonMover();
  // data
  MonteCarloBlock *pmy_mcb;

  // functions
  virtual void Move(MeshBlock *pmb, Photon *pphot);
};

//! \class CartesianMover
//  \brief derived class for Cartesian coordinates

class CartesianMover : public PhotonMover {
public:
  CartesianMover(MonteCarloBlock *pmcb);
  ~CartesianMover();

  // functions
  void Move(MeshBlock *pmb, Photon *pphot);
};

#endif // PHOTONMOVER_HPP
