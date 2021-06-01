#ifndef PHOTON_HPP
#define PHOTON_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photon.hpp
//  \brief definitions for Photon class

#include <complex>

// Athena++ classes headers
#include "../athena.hpp"
#include "montecarlo.hpp"

class MonteCarloBlock;

// photon status identifiers
enum PhotonStatus {EVOLVING = 0, ESCAPED = 1, ABSORBED = 2, DESTROYED = 3};
enum {IMC1 = 0, IMC2 = 1, IMC3 = 2, IMC0 = 3};

//! \class Photon
//  \brief phton data and functions

class Photon {
public:
  Photon(MonteCarloBlock *pmcb, int nuser);
  ~Photon();

  // data
  MonteCarloBlock* pmy_mcb; // ptr to MonteCarlo containing this Photon

  int i1,i2,i3; // zone indicies currently containing photon
  int status; // photon status (escaped, absorbed, evolving)
  int nuser_var; // number of user variables
  enum BoundaryFace face;

  // SWD: x can always include time, k could always include energy?
  // SWD: kcart maybe should be deprecated
  Real x[4];  // current photon position in spacetime
  Real k[4];  // photon direction (momentum vector) curvalinear
  Real dk[4]; // the change in photon direction used for general mover
  Real kcart[3]; // photon direction in cartesian coordinates
  Real stokes[4];  // stokes vectors
  Real weight; // photon statistical weight
  Real energy;  // photon energy
  Real *user_var; // storage for user variables
  AthenaArray<Real> trajectory; // Store trajectory
  Real sct_coef, abs_coef;  //scattering and absoprtion coefficients
  std::complex<Real> polten[4][4]; // the polarization tensor

  // functions
  void CopyPhoton(Photon *pphot);
  void PrintPhoton();
  bool IsNanPhoton();
  void AllocateUserVariables(int n);

};
#endif // PHOTON_HPP
