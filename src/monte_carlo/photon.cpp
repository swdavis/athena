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

Photon::Photon(MonteCarloBlock *pmcb, int nuser) {

  pmy_mcb = pmcb;
  weight = 1.0;
  face = FACE_UNDEF;
  nuser_var = nuser;
  if (nuser > 0)
    user_var = new Real[nuser];
  else
    user_var = NULL;

}

// destructor

Photon::~Photon() {
  
  if (user_var != NULL) delete [] user_var;

}

//----------------------------------------------------------------------------------------
//! \fn void Photon::CopyPhoton()
//  \brief Initialize photon from another photon

// rewrite this as a constructor? Not currently used
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

//----------------------------------------------------------------------------------------
//! \fn void Photon::IsNanPhoton()
//  \brief check for Nan in photon properties

bool Photon::IsNanPhoton() {

  if (isnan(weight)) return true;
  if (isnan(eweight)) return true;
  if (isnan(energy)) return true;
  for (int i=0; i<3; ++i) {
    if (isnan(x[i])) return true;
    if (isnan(k[i])) return true;
    if (isnan(stokes[i])) return true;
  }
  if (isnan(sct_coef)) return true;
  if (isnan(abs_coef)) return true;

  return false;
}

//----------------------------------------------------------------------------------------
//! \fn void Photon::PrintPhoton()
//  \brief print key properites

void Photon::PrintPhoton() {
  // Used primarily for debugging
  std::cout << "----------------------------" << std::endl
            << "Energy, weights: " << energy << " " << weight
	    << " " << eweight << std::endl
	    << "i: " << i1 << " " << i2 << " " << i3 <<std::endl
	    << "x: " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl
	    << "k: " << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << std::endl
	    << "kcart: " << kcart[0] << " " << kcart[1] << " "
	    << kcart[2] <<std::endl
	    << "stokes: " << stokes[0] << " " << stokes[1] << " "
	    << stokes[2] << std::endl
	    << "opacity: " << sct_coef << " " << abs_coef << std::endl;
  if (nuser_var > 0) {
    std::cout << "User vars:";
      for (int i=0; i<nuser_var; i++) {
        std::cout << " " << user_var[i];
      }
      std::cout << std::endl;
  }
  if (status == EVOLVING)
    std::cout << "EVOLVING" << std::endl;
  else if (status == ESCAPED)
    std::cout << "ESCAPED" << std::endl;
  else if (status == DESTROYED)
    std::cout << "DESTROYED" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void Photon::AllocateUserVariables(int n)
//  \brief allocate memory for user variables

void Photon::AllocateUserVariables(int n) {
  
  if (n > 0)
    user_var = new Real[n];
  nuser_var = n;
  

}
