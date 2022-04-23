//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//! \brief implementation of functions in class Photon

// C++ Standard Libraries
#include <vector>
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "photon.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"


int Photon::nint = 7;
int Photon::nreal = 20;
int Photon::naux = 0;
int Photon::nwork = 0;
int Photon::ncplx = 16;
int Photon::ipid = 0;
int Photon::inscp = 1, Photon::istatp = 2, Photon::itrp = 3;
int Photon::ii1p = 4, Photon::ii2p = 5, Photon::ii3p = 6;
int Photon::ix0p = 0, Photon::ix1p = 1, Photon::ix2p = 2, Photon::ix3p = 3;
int Photon::ik0p = 4, Photon::ik1p = 5, Photon::ik2p = 6, Photon::ik3p = 7;
int Photon::idk0p = 8, Photon::idk1p = 9, Photon::idk2p = 10, Photon::idk3p = 11;
int Photon::iep = 12, Photon::iwp = 13, Photon::iscp = 14, Photon::iacp = 15;
int Photon::isip = 16, Photon::isqp = 17, Photon::isup = 18, Photon::isvp = 19;

//----------------------------------------------------------------------------------------
//! Photon constructor

Photon::Photon(MonteCarloBlock *pmcb, int nuser, int len_limit)
  // Allocate space for photon data.
  : intprop(new std::vector<int> [nint]), realprop(new std::vector<Real> [nreal]),
    aux(new std::vector<Real> [naux]), work(new std::vector<Real> [nwork]),
    user(new std::vector<Real> [nuser]),
    polten(new std::vector<std::complex<Real>> [ncplx]),
    nphot(npar),pid(intprop[ipid]),nscp(intprop[inscp]), statp(intprop[istatp]),
    trp(intprop[itrp]), i1p(intprop[ii1p]), i2p(intprop[ii2p]), i3p(intprop[ii3p]),
    x0p(realprop[ix0p]), x1p(realprop[ix1p]), x2p(realprop[ix2p]), x3p(realprop[ix3p]),
    k0p(realprop[ik0p]), k1p(realprop[ik1p]), k2p(realprop[ik2p]), k3p(realprop[ik3p]),
    dk0p(realprop[idk0p]), dk1p(realprop[idk1p]), dk2p(realprop[idk2p]),
    dk3p(realprop[idk3p]),
    ep(realprop[iep]), wp(realprop[iwp]), scp(realprop[iscp]), acp(realprop[iacp]),
    sip(realprop[isip]), sqp(realprop[isqp]), sup(realprop[isup]), svp(realprop[isvp]) {

  pmy_mcb = pmcb;

  nphot_limit = len_limit;
  nuser_var = nuser;
  npar = 0;


}

//----------------------------------------------------------------------------------------
//! destructor

Photon::~Photon() {

}

//----------------------------------------------------------------------------------------
//! \fn void Photon::PrintPhoton(int ip)
//! \brief print key properites

void Photon::PrintPhoton(int ip) {
  // Used primarily for debugging
  std::cout << "----------------------------" << std::endl
            << "Energy, weight: " << ep[ip] << " " << wp[ip] << std::endl
            << "i: " << i1p[ip] << " " << i2p[ip] << " " << i3p[ip] <<std::endl
            << "x: " << x1p[ip] << " " << x2p[ip] << " " << x3p[ip] << " " << x0p[ip]
            << std::endl
            << "k: " << k1p[ip] << " " << k2p[ip] << " " << k3p[ip] << " " << k0p[ip]
            << std::endl
            << "dk: " << dk1p[ip] << " " << dk2p[ip] << " " << dk3p[ip] << " " << dk0p[ip]
            << std::endl
            << "stokes: " << sip[ip] << " " << sqp[ip] << " " << sup[ip] << std::endl
            << "opacity: " << scp[ip] << " " << acp[ip] << std::endl;
  if (nuser_var > 0) {
    std::cout << "User vars:";
      for (int i=0; i<nuser_var; i++) {
        std::cout << " " << user[i][ip];
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
//! \fn void Photon::IsNanPhoton(int ip)
//! \brief check for Nan in photon properties

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
  for (int i = 0; i < nuser_var; ++i)
    user[i].resize(new_npar);
  for (int i = 0; i < ncplx; ++i)
    polten[i].resize(new_npar);

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
      for (int j = 0; j < nuser_var; ++j)
        user[j][k] = user[j].back();
      for (int j = 0; j < ncplx; ++j)
        polten[j][k] = polten[j].back();
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
    for (int j = 0; j < nuser_var; ++j)
      user[j].pop_back();
    for (int j = 0; j < ncplx; ++j)
      polten[j].pop_back();

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

//----------------------------------------------------------------------------------------
//! \fn void Photon::PolarizationToTetrad(std::complex<Real> ttet[4][4], Real ecov[4][4],
//!                                       const int ip)
//!
//! \brief transform complex tensor from coordinate frame to tetrad frame

void Photon::PolarizationToTetrad(std::complex<Real> ttet[4][4], Real ecov[4][4],
                                  const int ip) {

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      ttet[i][j] = std::complex<Real>(0.,0.);

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      for (int k = 0; k < 4; k++)
        for (int l = 0; l < 4; l++) {
          ttet[i][j] += polten[k*4+l][ip] * ecov[i][k] * ecov[j][l];
        }

}

//----------------------------------------------------------------------------------------
//! \fn void PolarizationToCoord(std::complex<Real> ttet[4][4], Real econ[4][4],
//!                              const int ip)
//!
//! \brief transform complex tensor from tetrad frame to coordinate frame

void Photon::PolarizationToCoord(std::complex<Real> ttet[4][4], Real econ[4][4],
                                 const int ip) {

  for(int i = 0; i < NCOORD; i++)
    for(int j = 0; j < NCOORD; j++)
      polten[i*4+j][ip] = std::complex<Real>(0.,0.);

  for(int i = 0; i < NCOORD; i++)
    for(int j = 0; j < NCOORD; j++)
      for(int k = 0; k < NCOORD; k++)
        for(int l = 0; l < NCOORD; l++) {
          polten[i*4+j][ip] += ttet[k][l] * econ[k][i] * econ[l][j];
        }

}
