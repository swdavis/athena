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
#include "monte_carlo.hpp"

class MeshBlock;
class ParameterInput;
class MonteCarlo;
class Photon;

//! \class PhotonMover
//  \brief abstract base class for all derived classes

class PhotonMover {
public:
  PhotonMover(MonteCarlo *pmc);
  ~PhotonMover();
  // data
  MonteCarlo *pmy_mc;

  // functions
  virtual void Move(MeshBlock *pmb, Photon *pphot);
};

//! \class CartesianMover
//  \brief derived class for Cartesian coordinates

class CartesianMover : public PhotonMover {
public:
  CartesianMover(MonteCarlo *pmc);
  ~CartesianMover();

  // functions
  void Move(MeshBlock *pmb, Photon *pphot);
};

#endif // PHOTONMOVER_HPP
