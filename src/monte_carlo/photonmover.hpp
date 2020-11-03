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
#include "mcutils.hpp"

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

  Real dl; //displacement
  MonteCarloBlock *pmy_mcb;
  MCCoord *pcoord;

  // function pointers
  UserMoveFunc_t UserWorkInMove;

  // Arrays for MRW acceleration
  AthenaArray<Real> mrwprob;
  AthenaArray<Real> mrwdev;
  AthenaArray<Real> mrwxf;
  AthenaArray<Real> mrwt;
  AthenaArray<Real> mrwxi;
  AthenaArray<Real> mrwp;
  AthenaArray<Real> mrwrt;
  AthenaArray<Real> mrwrp;
  AthenaArray<Real> mrwrr;
  AthenaArray<Real> mrwta;
  AthenaArray<Real> mrwtp;
  AthenaArray<Real> mrwtt;

  int nmax,nxi,np,nt; //used for acceleration arrays

  bool acceleration;
  bool boosts;
  bool compton;
  bool time_acc;

  // functions
  virtual void Move(Photon *pphot);
  virtual Real GetOpticalDepth(MCRandom *pran);
  virtual void NextFace(Real dx1, Real dx2, Real dx3, int &face, Real &dx);
  virtual void MovePhotonToNextZone(Photon *pphot, MCCoord *pco,
    MonteCarloBlock *pmcb, int face, bool ascend[3]);
  virtual bool UpdateZone(Photon *pphot);
  // SWD Move Cartesian conversions to coordinates
  virtual void CartesianToCurvalinear(Photon *pphot);
  virtual void CurvalinearToCartesian(Photon *pphot);
  virtual void InitializeMRWDist(void);
  // Acceleration methods
  virtual bool MRWAcceleration(Photon *pphot, MCRandom *pran, Real dist, Real tauacc);
  virtual Real MRWDist(MCRandom *pran);
  virtual void ReadComptonGreensFunction(void);
  virtual Real InterpComptonEnergy(Real x0, Real time, Real prob);
  virtual void ReadRadiusDistribution(void);
  virtual void ReadTimeDistribution(void);
  virtual Real InterpPathTime(Real tau, Real prob);
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

//! \class GeneralMover
//  \brief derived class for general covariant coordinates

class GeneralMover : public PhotonMover {
public:
  GeneralMover(MonteCarloBlock *pmcb);
  ~GeneralMover();

  Real step_par;

  // functions
  void Move(Photon *pphot);
  void CartesianToCurvalinear(Photon *pphot);
  void CurvalinearToCartesian(Photon *pphot);
  void UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb);
  bool UpdateZone(Photon *pphot);
  void VerletStep(Photon *pphot, Real step);
  void PropogatePolarization(Photon *nphot);
  Real StepSize(Photon *pphot);

};

#endif // PHOTONMOVER_HPP
