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

// class variable initialization
bool Photon::initialized = false;
bool Photon::polarized = false;
bool Photon::general_mover_flag = false;

int Photon::inscp = -1, Photon::istatp = -1, Photon::itrp = -1;
int Photon::ii1p = -1, Photon::ii2p = -1, Photon::ii3p = -1;
int Photon::ix0p = -1, Photon::ix1p = -1, Photon::ix2p = -1, Photon::ix3p = -1;
int Photon::ik0p = -1, Photon::ik1p = -1, Photon::ik2p = -1, Photon::ik3p = -1;
int Photon::idk0p = -1, Photon::idk1p = -1, Photon::idk2p = -1, Photon::idk3p = -1;
int Photon::iep = -1, Photon::iwp = -1, Photon::iscp = -1, Photon::iacp = -1;
int Photon::isip = -1, Photon::isqp = -1, Photon::isup = -1, Photon::isvp = -1;
int Photon::iuserp = -1, Photon::ipolp = -1;
//----------------------------------------------------------------------------------------
//! Photon constructor

Photon::Photon(MonteCarloBlock *pmcb, ParameterInput *pin)
  : Particles(pmcb->pmy_block, pin),
  // Allocate space for photon data via initialization list
    //user(new std::vector<Real> [pmcb->pmy_mc->nuser_var]),
    //polten(new std::vector<std::complex<Real>> [ncplx]),
    nphot(npar),nscp(intprop[inscp]), statp(intprop[istatp]),
    trp(intprop[itrp]), i1p(intprop[ii1p]), i2p(intprop[ii2p]), i3p(intprop[ii3p]),
    x0p(rp[ix0p]), x1p(rp[ix1p]), x2p(rp[ix2p]), x3p(rp[ix3p]),
    k0p(rp[ik0p]), k1p(rp[ik1p]), k2p(rp[ik2p]), k3p(rp[ik3p]),
    dk0p(rp[idk0p]), dk1p(rp[idk1p]), dk2p(rp[idk2p]),
    dk3p(rp[idk3p]),
    ep(rp[iep]), wp(rp[iwp]), scp(rp[iscp]), acp(rp[iacp]),
  sip(rp[isip]), sqp(rp[isqp]), sup(rp[isup]), svp(rp[isvp]) {

  pmy_mcb = pmcb;
  nphot_limit = pmcb->pmy_mc->max_phots_init;
  nuser_var = pmcb->pmy_mc->nuser_var;
  user = &(rp[iuserp]);
  polten = &(cplxprop[ipolp]);
  npar = 0;


}

//----------------------------------------------------------------------------------------
//! destructor

Photon::~Photon() {

}

//----------------------------------------------------------------------------------------
//! \fn void Photon::PrintPhoton(int ip)
//! \brief print key photon properites

void Photon::PrintPhoton(int ip) {
  // Used primarily for debugging
  std::cout << "----------------------------" << std::endl
            << "Energy, weight: " << ep[ip] << " " << wp[ip] << std::endl
            << "i: " << i1p[ip] << " " << i2p[ip] << " " << i3p[ip] <<std::endl
            << "x: " << x1p[ip] << " " << x2p[ip] << " " << x3p[ip] << " " << x0p[ip]
            << std::endl
            << "k: " << k1p[ip] << " " << k2p[ip] << " " << k3p[ip] << " " << k0p[ip]
            << std::endl;
  if (general_mover_flag) {
    std::cout << "dk: " << dk1p[ip] << " " << dk2p[ip] << " " << dk3p[ip] << " "
              << dk0p[ip] << std::endl;
  }
  if (polarized) {
    std:: cout << "stokes: " << sip[ip] << " " << sqp[ip] << " " << sup[ip] << std::endl
               << "opacity: " << scp[ip] << " " << acp[ip] << std::endl;
  }
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
  if (polarized) {
    if (isnan(sip[ip])) return true;
    if (isnan(sqp[ip])) return true;
    if (isnan(sup[ip])) return true;
    if (isnan(svp[ip])) return true;
  }
  if (isnan(scp[ip])) return true;
  if (isnan(acp[ip])) return true;

  return false;
}

//--------------------------------------------------------------------------------------
//! \fn void Photon::AllocatePhotons(int nphot)
//! \brief Allocates photons

// SWD: Temporary --> converts protected function to public :(
void Photon::AllocatePhotons(int nphot) {
  // Call Resize function
  Resize(nphot);
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

//--------------------------------------------------------------------------------------
//! \fn Photon::Initialize(MonteCarloBlock *pmcb, ParameterInput *pin)
//! \brief initializes the Photon class.
// SWD: Change name to distinguish with InitializePhoton?
void Photon::Initialize(MonteCarlo *pmc, ParameterInput *pin) {

  // Initialize first the parent class.
  Particles::Initialize(pmc->pmy_mesh, pin);

  if (initialized) return;

  // Add particle ID and status flags, other int parameters.
  inscp = AddIntProperty("nscp");
  istatp = AddIntProperty("statp");
  itrp = AddIntProperty("trp");

  // Add photon position.
  ix0p = AddRealProperty("x0");
  ix1p = AddRealProperty("x1");
  ix2p = AddRealProperty("x2");
  ix3p = AddRealProperty("x3");

  // Add photon momentum.
  ik0p = AddRealProperty("k0");
  ik1p = AddRealProperty("k1");
  ik2p = AddRealProperty("k2");
  ik3p = AddRealProperty("k3");

  if (pmc->general_mover_flag) {
    general_mover_flag = true;
    // Add change in photon momentum.
    idk0p = AddRealProperty("dk0");
    idk1p = AddRealProperty("dk1");
    idk2p = AddRealProperty("dk2");
    idk3p = AddRealProperty("dk3");
  }

  // Add energy, weight, and opacities
  iep = AddRealProperty("ep");
  iwp = AddRealProperty("wp");
  iscp = AddRealProperty("scp");
  iacp = AddRealProperty("acp");

  if (pmc->polarized) {
    polarized = true;
    // Add stokes vectors
    isip = AddRealProperty("sip");
    isqp = AddRealProperty("sqp");
    isup = AddRealProperty("sup");
    isvp = AddRealProperty("svp");
    if (general_mover_flag) {
      // Add complex polarization tensor
      for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
          int idummy = AddComplexProperty("pol"+std::to_string(i)+std::to_string(j));
          if ( (i==0) && (j==0))
            ipolp = idummy;
          }
      }
    }
  }

  // Add particle position indices.
  ii1p = AddIntProperty("i1p");
  ii2p = AddIntProperty("i2p");
  ii3p = AddIntProperty("i3p");

  // Add nuser variables
  for (int i=0; i<pmc->nuser_var; i++) {
    int idummy = AddRealProperty("user"+std::to_string(i));
    if (i==0)
      iuserp = idummy;
  }

#ifdef MPI_PARALLEL
  // Get my MPI communicator.
  MPI_Comm_dup(MPI_COMM_WORLD, &my_comm);
#endif

  initialized = true;
}
