#ifndef PHOTON_HPP
#define PHOTON_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photon.hpp
//! \brief definitions for Photon class

// C++ libraries
#include <complex>
#include <vector>

// Athena++ classes headers
#include "../athena.hpp"
#include "montecarlo.hpp"

class MonteCarloBlock;

// photon status identifiers
enum PhotonStatus {EVOLVING = 0, ESCAPED = 1, ABSORBED = 2, DESTROYED = 3};
enum {IMC1 = 0, IMC2 = 1, IMC3 = 2, IMC0 = 3};

//---------------------------------------------------------------------------------------
//! \class Photon
//! \brief photon data and functions

class Photon {
public:
  Photon(MonteCarloBlock *pmcb, int nuser, int len_limit);
  ~Photon();

  // data
  MonteCarloBlock* pmy_mcb; // ptr to MonteCarlo containing this Photon

  //std::complex<Real> polten[4][4]; // the polarization tensor

  // ---------- New implementation -------------------------------------

  int nuser_var;
  int npar;
  int nphot_limit;
  int &nphot;
  static int nint;
  static int nreal;
  static int naux;
  static int nwork;
  static int ncplx;

  static int ipid;
  static int istatp, inscp, itrp;
  static int ii1p, ii2p, ii3p;
  static int ix0p, ix1p, ix2p, ix3p;
  static int ik0p, ik1p, ik2p, ik3p;
  static int idk0p, idk1p, idk2p, idk3p;
  static int iep, iwp, iscp, iacp;
  static int isip, isqp, isup, isvp;

  std::vector<int> *intprop;   //!>   integer properties
  std::vector<Real> *realprop; //!>   real properties
  std::vector<Real> *aux;      //!>   auxiliary properties (communicated when
                               //!>     particles moving to another meshblock)
  std::vector<Real> *work;     //!>   working arrays (not communicated)
  std::vector<Real> *user;     //!>   user variable arrays
  std::vector<std::complex<Real>> *polten;

  std::vector<int> &pid;                  //!>   particle ID
  std::vector<int> &statp, &nscp, &trp;
  std::vector<int> &i1p, &i2p, &i3p;
  std::vector<Real> &x0p, &x1p, &x2p, &x3p;
  std::vector<Real> &k0p, &k1p, &k2p, &k3p;
  std::vector<Real> &dk0p, &dk1p, &dk2p, &dk3p;
  std::vector<Real> &ep, &wp, &scp, &acp;
  std::vector<Real> &sip, &sqp, &sup, &svp;

  void Resize(int new_npar);
  void RemoveOneParticle(int k);
  void PrintPhoton(int ip);
  bool IsNanPhoton(int ip);
  void PolarizationToTetrad(std::complex<Real> ttet[4][4], Real ecov[4][4], const int ip);
  void PolarizationToCoord(std::complex<Real> ttet[4][4], Real econ[4][4], const int ip);
};
#endif // PHOTON_HPP
