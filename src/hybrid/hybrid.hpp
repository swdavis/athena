#ifndef HYBRID_HPP
#define HYBRID_HPP
//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//=======================================================================================
//! \file hybrid.hpp
//  \brief defines hybrid class
//  These classes contain data and functions used in Pegasus++
//=======================================================================================

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"

class MeshBlock;
class ParameterInput;

//! \class Hybrid
// \brief data/functions associated with hybrid-PIC algorithm

class Hybrid {
public:
  Hybrid(MeshBlock *pmb, ParameterInput *pin);     // create hybrid class
  ~Hybrid();                                       // destroy hybrid class

  MeshBlock* pmy_block;                            // MeshBlock pointer

  // standard plasma parameters
  Real beta;                                       // plasma beta
  Real ZTeTi;                                      // temperature ratio
  Real gamma;                                      // for electron eos, if adiabatic

  // problem-specific plasma parameters
  Real beta_prp, beta_prl;                         // parallel and perpendicular betas
  Real vinject;                                    // injector speed
  Real tcorr;                                      // driving correlation time
  Real dedt;                                       // energy injection due to driving
  Real eta_Ohm;                                    // Ohmic resistivity
  Real eta_hyper;                                  // hyper-resistivity

  // field-particle coupling arrays
  AthenaArray<Real> fcoup, fcoup_;                 // fields used to move particles
  AthenaArray<Real> mcoup, mcoup_;                 // deposited quantities


  void Initialize();                               // initialize hybrid-PIC scheme
  void NewBlockTimeStep();
  void TransposeMCOUP();
  //void CoumputeEMF();
};


// Definitions for hybrid algorithm
#define NFCOUP 7
#define NMCOUP 4
enum {IE1 = 3, IE2 = 4, IE3 = 5, IEB = 6}; // IB1, IB2, IB3, and mcoup variables
                                           // are in athena.hpp
#endif
