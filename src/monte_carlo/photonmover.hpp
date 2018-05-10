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
#include "../coordinates/coordinates.hpp"
#include "photon.hpp"
#include "montecarlo.hpp"

class MeshBlock;
class ParameterInput;
class MonteCarloBlock;
class Photon;
class Coordinate;

// Function prototypes
Real GetOpticalDepth(MCRandom *pran);

//! \class PhotonMover
//  \brief abstract base class for all derived classes

class PhotonMover {
public:
  PhotonMover(MonteCarloBlock *pmcb);
  ~PhotonMover();
  // data
  MonteCarloBlock *pmy_mcb;

  // functions
  virtual void Move(Photon *pphot);
  virtual Real GetOpticalDepth(MCRandom *pran);
  virtual void NextFace(Real dx1, Real dx2, Real dx3, int &face, Real &dx);
  virtual void UpdatePhotonPositionInZone(Photon *pphot, Real dl);
  virtual void MovePhotonToNextZone(Photon *pphot, Coordinates *pco,
                 MonteCarloBlock *pmcb, Real dl, int face, bool ascend[3]);

};

//! \class CartesianMover
//  \brief derived class for Cartesian coordinates

class CartesianMover : public PhotonMover {
public:
  CartesianMover(MonteCarloBlock *pmcb);
  ~CartesianMover();

  // functions
  void Move(Photon *pphot);
  void UpdatePhotonPositionInZone(Photon *pphot, Real dl);
  void MovePhotonToNextZone(Photon *pphot, Coordinates *pco,
         MonteCarloBlock *pmcb, Real dl, int face, bool ascend[3]);
};

#endif // PHOTONMOVER_HPP
