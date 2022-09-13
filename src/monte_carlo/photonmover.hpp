#ifndef PHOTONMOVER_HPP
#define PHOTONMOVER_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photonmover.hpp
//! \brief defines abstract base derived classes for moving photons

// Athena++ headers
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

//----------------------------------------------------------------------------------------
//! \class PhotonMover
//! \brief abstract base class for all derived classes

class PhotonMover {
public:
  PhotonMover(MonteCarloBlock *pmcb);
  ~PhotonMover();
  // data

  Real dl; // current displacement
  int checkmove; // check/terminate move

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
  bool resonance;
  bool compton;
  bool time_acc;

  // functions
  virtual void Move(Photon *pphot, int ips, int ipe);
  virtual Real GetOpticalDepth(MCRandom *pran);
  virtual Real GetExtinctionCoefficient(Real ac, Real sc);
  virtual Real ExpTauAbsorption(Real ac, Real dl);
  virtual void NextFace(Real dx1, Real dx2, Real dx3, int &face, Real &dx);
  virtual void MovePhotonToNextZone(Photon *pphot, MCCoord *pco,
               MonteCarloBlock *pmcb, int face, bool ascend[3], int ip);
  virtual bool UpdateZone(Photon *pphot, int ip);
//  virtual bool UpdateSingleZone(Photon *pphot, int ip, bool *multizone);
  virtual void CurvalinearToCartesian(Photon *pphot, Real kcart[4]);
  virtual void InitializeMRWDist(void);
  // Acceleration methods
  virtual Real SampleEscapeTime(MCRandom *pran, Real decayRate, Real sphereRadius,
                                Real diffusionTime);
  virtual bool MRWAcceleration(Photon *pphot, MCRandom *pran, Real dist, Real tauacc,
                               int ip);
  virtual Real MRWResonanceAcceleration(Photon *pphot, MCRandom *pran, Real dist, 
                               Real tauacc, int ip);
  virtual Real MRWDist(MCRandom *pran);
  virtual void ReadComptonGreensFunction(void);
  virtual Real InterpComptonEnergy(Real x0, Real time, Real prob);
  virtual void ReadRadiusDistribution(void);
  virtual void ReadTimeDistribution(void);
  virtual Real InterpPathTime(Real tau, Real prob);
};

//----------------------------------------------------------------------------------------
//! \class CartesianMover
//! \brief derived class for moving in Cartesian coordinates

class CartesianMover : public PhotonMover {
public:
  CartesianMover(MonteCarloBlock *pmcb);
  ~CartesianMover();

  // functions
  void Move(Photon *pphot, int ips, int ipe);

};

//----------------------------------------------------------------------------------------
//! \class SphericalPolarMover
//! \brief derived class for moving in spherical-polar coordinates

class SphericalPolarMover : public PhotonMover {
public:
  SphericalPolarMover(MonteCarloBlock *pmcb);
  ~SphericalPolarMover();

  // functions
  void Move(Photon *pphot, int ips, int ipe);
  void CurvalinearToCartesian(Photon *pphot, Real kcart[4]);
};

//----------------------------------------------------------------------------------------
//! \class SphericalPolarAltMover
//! \brief derived class for moving in spherical-polar coordinates

class SphericalPolarAltMover : public PhotonMover {
public:
  SphericalPolarAltMover(MonteCarloBlock *pmcb);
  ~SphericalPolarAltMover();

  Real step_par;
  Real gamma[NCOORD][NCOORD][NCOORD];

  // functions
  bool UpdateSingleZone(Photon *pphot, int ip, bool *multizone);
  void Move(Photon *pphot, int ips, int ipe);
  void UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb, int ip);
};

//----------------------------------------------------------------------------------------
//! \class GeneralMover
//! \brief derived class for moving in general coordinates

class GeneralMover : public PhotonMover {
public:
  GeneralMover(MonteCarloBlock *pmcb);
  ~GeneralMover();

  Real step_par;
  Real gamma[NCOORD][NCOORD][NCOORD];

  // functions
  void Move(Photon *pphot, int ips, int ipe);
  void CurvalinearToCartesian(Photon *pphot, Real kcart[4]);
  void UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb, int ip);
  void VerletStep(Photon *pphot, Real step, int ip);
  void PropogatePolarization(Photon *nphot, Real step, int ip);
  Real StepSize(Photon *pphot, int ip);

};

#endif // PHOTONMOVER_HPP
