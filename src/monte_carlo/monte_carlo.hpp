#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.hpp
//  \brief definitions for MonteCarlo class

#include <gsl/gsl_randist.h>

// Athena++ classes headers
#include "../athena.hpp"
//#include "../athena_arrays.hpp"
#include "photon.hpp"

class MeshBlock;
class ParameterInput;
class Photon;
class PhotonMover;


// Current design focusses on implementing static post-processing so this class
// implementation will likely evolve.

enum {TOEUL=0, TOCOM=1};

//! \class MCRandom
//  \brief monte carlo random number generator

class MCRandom {
public:
  MCRandom(int iseed);
  ~MCRandom();

  gsl_rng *dev;
  
  Real uniform();
};

//! \class MonteCarlo
//  \brief monte carlo functions and data

class MonteCarlo {
public:
  MonteCarlo(MeshBlock *pmb, ParameterInput *pin);
  ~MonteCarlo();

  // data
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this MonteCarlo
  Photon* pphoton; // ptr to photon packet
  PhotonMover* pmover; // ptr to photon mover

  int ntot;  // total number of photons for this block;
  bool zone_weight; // flag for zone weighting
  MCRandom *pran;

  // functions
  void TransferPhotons();  // Transfer photons on this block

};



#endif // MONTECARLO_HPP
