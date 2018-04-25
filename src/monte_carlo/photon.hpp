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
#include "monte_carlo.hpp"

class MonteCarlo;

// photon status identifiers
enum PhotonStatus {EVOLVING = 0, ESCAPED = 1, DESTROYED = 2};

//! \class Photon
//  \brief phton data and functions

class Photon {
public:
  Photon(MonteCarlo *pmc);
  ~Photon();

  // data
  MonteCarlo* pmy_mc; // ptr to MonteCarlo containing this Photon

  int izone[3]; // Zone currently containing photon
  int status; // photon status (escaped, absorbed, evolving)

  Real x[3];  // current photon position
  Real k[3];  // photon direction (momentum vector)
  Real stokes[3];  // stokes vectors for linear polarization
  Real weight, cweight; // photon weights
  Real energy;  // photon energy

  Real sct_coef, abs_coef;  //scattering and absoprtion coefficients

  // functions
  void InitializePhoton(MeshBlock *pmb); // defined in problem gen.

};
#endif // PHOTON_HPP
