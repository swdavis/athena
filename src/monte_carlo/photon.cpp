//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//  \brief implementation of functions in class Photon

// C++ Standard Libraries
#include <vector>
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "photon.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"


int Photon::nint = 5;
int Photon::nreal = 20;
int Photon::naux = 0;
int Photon::nwork = 0;
int Photon::ipid = 0;
int Photon::istatp = 1, Photon::ii1p = 2, Photon::ii2p = 3, Photon::ii3p = 4;
int Photon::ix0p = 0, Photon::ix1p = 1, Photon::ix2p = 2, Photon::ix3p = 3;
int Photon::ik0p = 4, Photon::ik1p = 5, Photon::ik2p = 6, Photon::ik3p = 7;
int Photon::idk0p = 8, Photon::idk1p = 9, Photon::idk2p = 10, Photon::idk3p = 11;
int Photon::iep = 12, Photon::iwp = 13, Photon::iscp = 14, Photon::iacp = 15;
int Photon::isip = 16, Photon::isqp = 17, Photon::isup = 18, Photon::isvp = 19;


// constructor, initializes data structures and parameters

Photon::Photon(MonteCarloBlock *pmcb, int nuser, int len_limit)
  // Allocate space for particle data.
  : intprop(new std::vector<int> [nint]), realprop(new std::vector<Real> [nreal]),
    aux(new std::vector<Real> [naux]), work(new std::vector<Real> [nwork]),
    nphot(npar),pid(intprop[ipid]),
    statp(intprop[istatp]), i1p(intprop[ii1p]), i2p(intprop[ii2p]), i3p(intprop[ii3p]),
    x0p(realprop[ix0p]), x1p(realprop[ix1p]), x2p(realprop[ix2p]), x3p(realprop[ix3p]),
    k0p(realprop[ik0p]), k1p(realprop[ik1p]), k2p(realprop[ik2p]), k3p(realprop[ik3p]),
    dk0p(realprop[idk0p]), dk1p(realprop[idk1p]), dk2p(realprop[idk2p]), dk3p(realprop[idk3p]),
    ep(realprop[iep]), wp(realprop[iwp]), scp(realprop[iscp]), acp(realprop[iacp]),
    sip(realprop[isip]), sqp(realprop[isqp]), sup(realprop[isup]), svp(realprop[isvp]) {

  pmy_mcb = pmcb;
  weight = 1.0;
  face = FACE_UNDEF;
  nphot_limit = len_limit;
  nuser_var = nuser;
  npar = 0;
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
  energy = pphot->energy;
  sct_coef = pphot->sct_coef;
  abs_coef = pphot->abs_coef;
    
}

//----------------------------------------------------------------------------------------
//! \fn void Photon::IsNanPhoton()
//  \brief check for Nan in photon properties

bool Photon::IsNanPhoton() {

  if (isnan(weight)) return true;
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
//! \fn void Photon::PrintPhoton(int ip)
//  \brief print key properites

void Photon::PrintPhoton(int ip) {
  // Used primarily for debugging
  std::cout << "----------------------------" << std::endl
            << "Energy, weight: " << ep[ip] << " " << wp[ip] << std::endl
	    << "i: " << i1p[ip] << " " << i2p[ip] << " " << i3p[ip] <<std::endl
	    << "x: " << x1p[ip] << " " << x2p[ip] << " " << x3p[ip] << " " << x0p[ip] 
            << std::endl
	    << "k: " << k1p[ip] << " " << k2p[ip] << " " << k3p[ip] << " " << k0p[ip] 
            << std::endl
	    << "stokes: " << sip[ip] << " " << sqp[ip] << " " << sup[ip] << std::endl
	    << "opacity: " << scp[ip] << " " << acp[ip] << std::endl;
  if (nuser_var > 0) {
    std::cout << "User vars:";
      for (int i=0; i<nuser_var; i++) {
        std::cout << " " << user_var[i];
      }
      std::cout << std::endl;
  }
  if (statp[ip] == EVOLVING)
    std::cout << "EVOLVING" << std::endl;
  else if (statp[ip] == ESCAPED)
    std::cout << "ESCAPED" << std::endl;
  else if (statp[ip] == DESTROYED)
    std::cout << "DESTROYED" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void Photon::PrintPhoton()
//  \brief print key properites

void Photon::PrintPhoton() {
  // Used primarily for debugging
  std::cout << "----------------------------" << std::endl
            << "Energy, weight: " << energy << " " << weight << std::endl
	    << "i: " << i1 << " " << i2 << " " << i3 <<std::endl
	    << "x: " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl
	    << "k: " << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << std::endl
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


//----------------------------------------------------------------------------------------
//! \fn void Photon::PopulateWorkingArrays(int n)
//  \brief Copies vector element n to arrays used by scattering functions

void Photon::VectorsToWorkingArrays(int n) {
   
  // Copy integer variables
  status = statp[n];
  i1 = i1p[n];
  i2 = i2p[n];
  i3 = i3p[n];

  // Copy real variables
  x[IMC0] = x0p[n];
  x[IMC1] = x1p[n];
  x[IMC2] = x2p[n];
  x[IMC3] = x3p[n];
  k[IMC0] = k0p[n];
  k[IMC1] = k1p[n];
  k[IMC2] = k2p[n];
  k[IMC3] = k3p[n];
  dk[IMC0] = dk0p[n];
  dk[IMC1] = dk1p[n];
  dk[IMC2] = dk2p[n];
  dk[IMC3] = dk3p[n];
  energy = ep[n];
  weight = wp[n];
  sct_coef = scp[n];
  abs_coef = acp[n];
  stokes[0] = sip[n];
  stokes[1] = sqp[n];
  stokes[2] = sup[n];
  stokes[3] = svp[n];

}

//----------------------------------------------------------------------------------------
//! \fn void Photon::WorkingArraysToVectors(int n)
//  \brief Copies vector element n to arrays used by scattering functions

void Photon::WorkingArraysToVectors(int n) {

  // Copy integer variables
  statp[n] = status;
  i1p[n] = i1;
  i2p[n] = i2;
  i3p[n] = i3;

  // Copy real vairables
  x0p[n] = x[IMC0];
  x1p[n] = x[IMC1];
  x2p[n] = x[IMC2];
  x3p[n] = x[IMC3]; 
  k0p[n] = k[IMC0];
  k1p[n] = k[IMC1];
  k2p[n] = k[IMC2];
  k3p[n] = k[IMC3];
  dk0p[n] = dk[IMC0];
  dk1p[n] = dk[IMC1];
  dk2p[n] = dk[IMC2];
  dk3p[n] = dk[IMC3];
  ep[n] = energy;
  wp[n] = weight;
  scp[n] = sct_coef;
  acp[n] = abs_coef;
  sip[n] = stokes[0];
  sqp[n] = stokes[1];
  sup[n] = stokes[2];
  svp[n] = stokes[3];

}

//----------------------------------------------------------------------------------------
//! \fn void Photon::IsNanPhoton(int ip)
//  \brief check for Nan in photon properties

bool Photon::IsNanPhoton(int ip) {

  if (isnan(wp[ip])) return true;
  if (isnan(ep[ip])) return true;
  if (isnan(x0p[ip])) return true;
  if (isnan(x1p[ip])) return true;
  if (isnan(x2p[ip])) return true;
  if (isnan(x3p[ip])) return true;
  if (isnan(k0p[ip])) return true;
  if (isnan(k1p[ip])) return true;
  if (isnan(k2p[ip])) return true;
  if (isnan(k3p[ip])) return true;
  if (isnan(sip[ip])) return true;
  if (isnan(sqp[ip])) return true;
  if (isnan(sup[ip])) return true;
  if (isnan(scp[ip])) return true;
  if (isnan(acp[ip])) return true;

  return false;
}


// Everything that follows was stolen from Particles class and will be replaced by
// Particle Class routines when Photons is transformed to derived class

//--------------------------------------------------------------------------------------
//! \fn void Particles::Resize(int new_npar)
//! \brief changes number of particles.

void Photon::Resize(int new_npar) {
  // Resize the particle arrays.
  for (int i = 0; i < nint; ++i)
    intprop[i].resize(new_npar);
  for (int i = 0; i < nreal; ++i)
    realprop[i].resize(new_npar);
  for (int i = 0; i < naux; ++i)
    aux[i].resize(new_npar);
  for (int i = 0; i < nwork; ++i)
    work[i].resize(new_npar);

  // Flag new particles.
  for (int k = npar; k < new_npar; ++k)
    pid[k] = -1;

  // Update number of particles.
  npar = new_npar;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::RemoveOneParticle(int k)
//! \brief removes particle k in the block.

void Photon::RemoveOneParticle(int k) {
  if (0 <= k && k < npar) {
    if (--npar != k) {
      // Replace the k-th particle by the last particle.
      for (int j = 0; j < nint; ++j)
        intprop[j][k] = intprop[j].back();
      for (int j = 0; j < nreal; ++j)
        realprop[j][k] = realprop[j].back();
      for (int j = 0; j < naux; ++j)
        aux[j][k] = aux[j].back();
      for (int j = 0; j < nwork; ++j)
        work[j][k] = work[j].back();
    }
    // Remove the last particle.
    for (int j = 0; j < nint; ++j)
      intprop[j].pop_back();
    for (int j = 0; j < nreal; ++j)
      realprop[j].pop_back();
    for (int j = 0; j < naux; ++j)
      aux[j].pop_back();
    for (int j = 0; j < nwork; ++j)
      work[j].pop_back();
  } else {
    // Throw error when index k is invalid.
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::RemoveOneParticle]" << std::endl
        << "\tk = " << k << ", npar = " << npar << std::endl
        << "Index k is out of range. " << std::endl;
    throw std::runtime_error(msg.str().c_str());
    //ATHENA_ERROR(msg);
  }
}
