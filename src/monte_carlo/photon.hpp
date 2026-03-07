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
enum PhotonStatus {EVOLVING = 0, ESCAPED = 1, ABSORBED = 2, DESTROYED = 3, BUFFERED = 4};
enum {IMC0 = 0, IMC1 = 1, IMC2 = 2, IMC3 = 3};
//enum {IMC1 = 0, IMC2 = 1, IMC3 = 2, IMC0 = 3};

//---------------------------------------------------------------------------------------
//! \class Photon
//! \brief photon data and functions

class Photon : public Particles {
public:
  Photon(MonteCarloBlock *pmcb, ParameterInput *pin);
  ~Photon();

  // public functions
  //void RemoveOneParticle(int k);
  void PrintPhoton(const std::string &msg, int ip);
  void PrintPhoton(int ip);
  bool IsNanPhoton(int ip);
  void PolarizationToTetrad(std::complex<Real> ttet[4][4], Real ecov[4][4], const int ip);
  void PolarizationToCoord(std::complex<Real> ttet[4][4], Real econ[4][4], const int ip);
  void AllocatePhotons(int nphot);
  void SendToNeighbors();
  void ApplyPeriodicBoundary(Real &x1, Real &x2, Real &x3, int k);
  bool ReceiveFromNeighbors();
  void GetPositionIndices(int ibegin, int iend);
  static void Initialize(MonteCarlo *pmc, ParameterInput *pin);

  // public data
  // SWD: should be reorganized with tighter access control for some variables
  MonteCarloBlock* pmy_mcb; // ptr to MonteCarlo currently containing this Photon

  int nuser_var;
  int nphot_limit;
  int &nphot;

  static int istatp, inscp, ityp;
  static int ii1p, ii2p, ii3p;
  static int ix0p, ix1p, ix2p, ix3p;
  static int ik0p, ik1p, ik2p, ik3p;
  static int idk0p, idk1p, idk2p, idk3p;
  static int iep, iwp, iscp, iacp;
  static int isip, isqp, isup, isvp;
  static int iuserp;
  static int ipolp;
  static int idtp;

  std::vector<int> &statp, &nscp, &type;
  std::vector<int> &i1p, &i2p, &i3p;
  std::vector<Real> &x0p, &x1p, &x2p, &x3p;
  std::vector<Real> &k0p, &k1p, &k2p, &k3p;
  std::vector<Real> &dk0p, &dk1p, &dk2p, &dk3p;
  std::vector<Real> &ep, &wp, &scp, &acp;
  std::vector<Real> &sip, &sqp, &sup, &svp;
  std::vector<Real> &dtp;
  std::vector<Real> *user;     //!>   user variable arrays
  std::vector<std::complex<Real>> *polten; //!> polarization tensor



  static bool initialized;
  static bool polarized;
  static bool general_pusher_flag;


  //#ifdef MPI_PARALLEL
  //static MPI_Comm my_comm;   //!> my MPI communicator
  //ParticleBuffer send_[56];  //!> particle send buffers
  //#endif


};
#endif // PHOTON_HPP
