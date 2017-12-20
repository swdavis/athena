//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//=======================================================================================
//! \file particle.hpp
//  \brief defines particle classes
//  These classes contain data and functions related to particles
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
#include "particle.hpp"
#include "../globals.hpp"

//! \class Particle
// \brief data/functions associated with list of particles

void Particle::Initialize()
{
  return;
}

void Particle::Test()
{ 
  MoveDeposit();
  Deposit();
  Move();
  GetPos();
  PackParticle();
  //ExchangeCounts();
  ReceiveCounts();
  SendCounts();
  WaitCounts();
  //CheckCnt();
  //ExchangeParticle();
  ReceiveParticle();
  SendParticle();
  WaitParticle();
  UnPackParticle();


//  if (nparticle==1){
//    std::cout << pmy_block->pmy_mesh->time << "\t"<< x1(0) << "\t" << x2(0) << std::endl;
//  }

  return;
}
