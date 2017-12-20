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
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "hybrid.hpp"


//! \class Hybrid
// \brief data/functions associated with hybrid-PIC algorithm

void Hybrid::Initialize()
{
  int ncells1 = pmy_block->block_size.nx1 + 2*NGHOST;
  int ncells2 = pmy_block->block_size.nx2 + 2*NGHOST;
  int ncells3 = pmy_block->block_size.nx3 + 2*NGHOST;

//  for (int k=pmy_block->ks-NGHOST; k<=pmy_block->ke+NGHOST; k++) {
//    for (int j=pmy_block->js-NGHOST; j<=pmy_block->je+NGHOST; j++) {
//      for (int i=pmy_block->is-NGHOST; i<=pmy_block->ie+NGHOST; i++) {

  for (int k=0; k<ncells3; k++) {
    for (int j=0; j<ncells2; j++) { 
      for (int i=0; i<ncells1; i++) {
        mcoup(k,j,i,IDN) = 0;
        mcoup(k,j,i,IM1) = 0;
        mcoup(k,j,i,IM2) = 0;
        mcoup(k,j,i,IM3) = 0;
        fcoup(k,j,i,IB1) = 0;
        fcoup(k,j,i,IB2) = 0;
        fcoup(k,j,i,IB3) = 1;
        fcoup(k,j,i,IE1) = 0;
        fcoup(k,j,i,IE2) = 0;
        fcoup(k,j,i,IE3) = 0;
        fcoup(k,j,i,IEB) = 0;
  }}}
  TransposeMCOUP();
  return;
}
