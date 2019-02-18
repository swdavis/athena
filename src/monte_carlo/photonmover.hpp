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

//---------------------- prototypes for photon moving ------------------------------------
Real GetOpticalDepth(MCRandom *pran);
//---------------------- prototypes for accleration via mrw ------------------------------

//! \class PhotonMover
//  \brief abstract base class for all derived classes

class PhotonMover {
public:
  PhotonMover(MonteCarloBlock *pmcb);
  ~PhotonMover();
  // data
  MonteCarloBlock *pmy_mcb;
  // Arrays for MRW acceleration
  AthenaArray<Real> mrwprob;
  AthenaArray<Real> mrwdev;
  bool acceleration;
  bool lorentz_transform;
  // functions
  virtual void Move(Photon *pphot);
  virtual Real GetOpticalDepth(MCRandom *pran);
  virtual void NextFace(Real dx1, Real dx2, Real dx3, int &face, Real &dx);
  virtual void MovePhotonToNextZone(Photon *pphot, MCCoord *pco,
    MonteCarloBlock *pmcb, int face, bool ascend[3]);
  virtual void UpdateZone(Photon *pphot);
  virtual void CartesianToCurvalinear(Photon *pphot);
  virtual void CurvalinearToCartesian(Photon *pphot);
  virtual void InitializeMRWDist(void);
  virtual bool MRWAcceleration(Photon *pphot, MCRandom *pran, Real dist, Real tauacc);
  virtual Real MRWDist(MCRandom *pran);
  
};

//! \class CartesianMover
//  \brief derived class for Cartesian coordinates

class CartesianMover : public PhotonMover {
public:
  CartesianMover(MonteCarloBlock *pmcb);
  ~CartesianMover();

  // functions
  void Move(Photon *pphot);

};

//! \class SphericalPolarMover
//  \brief derived class for spherical-polar coordinates

class SphericalPolarMover : public PhotonMover {
public:
  SphericalPolarMover(MonteCarloBlock *pmcb);
  ~SphericalPolarMover();

  // functions
  void Move(Photon *pphot);
  void CartesianToCurvalinear(Photon *pphot);
  void CurvalinearToCartesian(Photon *pphot);
};


#endif // PHOTONMOVER_HPP
