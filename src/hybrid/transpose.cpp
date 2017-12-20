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

void Hybrid::TransposeMCOUP()
{
  int ncells1 = pmy_block->block_size.nx1 + 2*NGHOST;
  int ncells2 = pmy_block->block_size.nx2 + 2*NGHOST;
  int ncells3 = pmy_block->block_size.nx3 + 2*NGHOST;

  for (int k=0; k<ncells3; k++) {
    for (int j=0; j<ncells2; j++) { 
      for (int i=0; i<ncells1; i++) {
        mcoup_(IDN,k,j,i) = mcoup(k,j,i,IDN);
        mcoup_(IM1,k,j,i) = mcoup(k,j,i,IM1);
        mcoup_(IM2,k,j,i) = mcoup(k,j,i,IM2);
        mcoup_(IM3,k,j,i) = mcoup(k,j,i,IM3);
  }}}

  return;
}
