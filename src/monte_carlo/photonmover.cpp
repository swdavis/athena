//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photonmover.cpp
//  \brief implementation for photon moving functions

// Athena++ headers
#include "photon.hpp"
#include "photonmover.hpp"
#include "../mesh/mesh.hpp"

// Implementation of base class

PhotonMover::PhotonMover(MonteCarloBlock *pmcb) {

  pmy_mcb = pmcb;
}

PhotonMover::~PhotonMover() {

}

void PhotonMover::Move(Photon *pphot) {

}

//----------------------------------------------------------------------------------------
//! \fn Real PhotonMover::GetOpticalDepth(MCRandom *pran)
//  \brief return exponentially distributed optical depth variable

Real PhotonMover::GetOpticalDepth(MCRandom *pran) {

  Real dev = pran->uniform();  
  while(dev <= 0.)
    dev=pran->uniform();
  //std::cout << dev << std::endl;
  return -log(dev);
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonMover::NextFace(Real dx1, Real dx2, Real dx3, int &face, Real &dx)
//  \brief returns flag with next face and distance to next face

void PhotonMover::NextFace(Real dx1, Real dx2, Real dx3, int &face, Real &dx)
{
// face tells which cell coordinates need to be updatde
//   x:   0
//   y:   1
//   z:   2
//   xy:  3
//   yz:  4
//   xz:  5
//   xyz: 6

  dx = dx1;

  if(dx2 < dx) {
    dx = dx2;
    if(dx3 < dx) {
      dx = dx3;
      face = 2;
      return;
    } else if(dx3 > dx) {
      face = 1; 
      return;
    } else {
      face = 4;
      return;
    }
  } else if(dx2 > dx) {
    if(dx3 < dx) {
      dx = dx3;
      face = 2;
      return;
    } else if(dx3 > dx) {
      face = 0;
      return;
    } else {
      face = 5;
      return;
    }
  } else {
    if(dx3 < dx) {
      dx = dx3;
      face = 2;
      return;
    } else if(dx3 > dx) {
      face = 3;
      return;
    } else {
      face = 6;
      return;
    }
  }
}

// virtual functions that must be provied in derived classes

void PhotonMover::UpdatePhotonPositionInZone(Photon *pphot, Real dl) {

}

void PhotonMover::MovePhotonToNextZone(Photon *pphot, Coordinates *pco,
  MonteCarloBlock *pmcb, Real dl, int face, bool higher[3]) {

}
