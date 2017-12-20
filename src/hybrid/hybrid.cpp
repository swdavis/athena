//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//=======================================================================================
//! \file hybrid.hpp
//  \brief defines hybrid class
//  These classes contain data and functions used in Pegasus++
//=======================================================================================

// C++ headers
#include <iostream>   // endl
#include <stdexcept>  // runtime_error
#include <sstream>    //stringstream

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "hybrid.hpp"

//! \class Hybrid
// \brief data/functions associated with hybrid-PIC algorithm

Hybrid::Hybrid (MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  int ncells1 = pmb->block_size.nx1 + 2*NGHOST;
  int ncells2 = pmb->block_size.nx2 + 2*NGHOST;
  int ncells3 = pmb->block_size.nx3 + 2*NGHOST;

  fcoup.NewAthenaArray(ncells3, ncells2, ncells1, NFCOUP);
  mcoup.NewAthenaArray(ncells3, ncells2, ncells1, NMCOUP);

  fcoup_.NewAthenaArray(NFCOUP,ncells3, ncells2, ncells1);
  mcoup_.NewAthenaArray(NMCOUP,ncells3, ncells2, ncells1);

  beta = pin->GetReal("problem", "beta");
  ZTeTi = pin->GetOrAddReal("problem", "ZTeTi",1.0);
  gamma = pin->GetOrAddReal("problem", "gamma",5.0/3.0);

  beta_prp = pin->GetOrAddReal("problem", "beta_prp", beta);
  beta_prl = pin->GetOrAddReal("problem", "beta_prl", beta);
  vinject = pin->GetOrAddReal("problem", "vinject", 0.0);
  tcorr = pin->GetOrAddReal("problem", "tcorr", 1.0);
  dedt = pin->GetOrAddReal("problem", "dedt", 0.0);
  eta_Ohm = pin->GetOrAddReal("problem", "eta_Ohm", 0.0);
  eta_hyper = pin->GetOrAddReal("problem", "eta_hyper", 0.0);
}

Hybrid::~Hybrid()
{
  fcoup.DeleteAthenaArray();
  mcoup.DeleteAthenaArray();
}
