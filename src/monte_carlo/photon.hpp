#ifndef PHOTON_HPP
#define PHOTON_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photon.hpp
//  \brief definitions for Photon class

// Athena++ classes headers
#include "../athena.hpp"
#include "montecarlo.hpp"

//

class MonteCarloBlock;

// photon status identifiers
enum PhotonStatus {EVOLVING = 0, ESCAPED = 1, DESTROYED = 2};

//! \class Photon
//  \brief phton data and functions

class Photon {
public:
  Photon(MonteCarloBlock *pmcb);
  ~Photon();

  // data
  MonteCarloBlock* pmy_mcb; // ptr to MonteCarlo containing this Photon

  int i1,i2,i3; // zone indicies currently containing photon
  int status; // photon status (escaped, absorbed, evolving)

  Real x[3];  // current photon position
  Real k[3];  // photon direction (momentum vector) curvalinear
  Real kcart[3]; // photon direction in cartesian coordinates
  Real stokes[3];  // stokes vectors for linear polarization
  Real weight, eweight; // photon weights
  Real energy;  // photon energy

  Real sct_coef, abs_coef;  //scattering and absoprtion coefficients

  // functions
  

};
#endif // PHOTON_HPP
