//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//  \brief implementation of functions in class Photon

// Athena++ headers
#include "photon.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// constructor, initializes data structures and parameters

Photon::Photon(MonteCarloBlock *pmcb) {

  pmy_mcb = pmcb;
  weight = 1.0;
  face = FACE_UNDEF;
}

// destructor

Photon::~Photon() {

}

// rewrite this as a constructor
void Photon::CopyPhoton(Photon *pphot) {

  i1 = pphot->i1;
  i2 = pphot->i2;
  i3 = pphot->i3;
  status = pphot->status;
  for(int i=0; i<3; ++i) {
    x[i] = pphot->x[i];
    k[i] = pphot->k[i];
    stokes[i] = pphot->stokes[i];
  }

  weight = pphot->weight;
  eweight = pphot->eweight;
  energy = pphot->energy;
  sct_coef = pphot->sct_coef;
  abs_coef = pphot->abs_coef;
    
}


// for debugging purposes
void Photon::PrintPhoton() {

  std::cout << "----------------------------" << std::endl
            << "Energy, weights: " << energy << " " << weight
	    << " " << eweight << std::endl
	    << "i: " << i1 << " " << i2 << " " << i3 <<std::endl
	    << "x: " << x[0] << " " << x[1] << " " << x[2] <<std::endl
	    << "k: " << k[0] << " " << k[1] << " " << k[2] <<std::endl
	    << "kcart: " << kcart[0] << " " << kcart[1] << " "
	    << kcart[2] <<std::endl
	    << "stokes: " << stokes[0] << " " << stokes[1] << " "
	    << stokes[2] <<std::endl;
}
