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

//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::MovePhotonToNextZone()
//  \brief updates photon zone

void PhotonMover::MovePhotonToNextZone(Photon *pphot, Coordinates *pco,
  MonteCarloBlock *pmcb, int face, bool ascend[3]) {
  
  // Update face(s) and adjust positions to lie exactly on boundary
  if ((face == 0) || (face == 3) || (face == 5) || (face == 6)) { //update x1 face
    if (ascend[0]) {
      pphot->i1++;
      if(pphot->i1 <= pmcb->ie)
        pphot->x[0] = pco->x1f(pphot->i1);
      else {
        pmcb->pbval->BoundaryFunction_[OUTER_X1](pmcb,pco,pphot);
        if (pphot->status == ESCAPED) {
          pphot->face = OUTER_X1;
        }
      }
    } else {
      pphot->i1--;
      if(pphot->i1 >= pmcb->is)
        pphot->x[0] = pco->x1f(pphot->i1+1);
      else {
        pmcb->pbval->BoundaryFunction_[INNER_X1](pmcb,pco,pphot);
        if (pphot->status == ESCAPED) {
          pphot->face = INNER_X1;
        }
      }
    }
  }
  if ((face == 1) || (face == 3) || (face == 4) || (face == 6)) { //update x2 face
    if (ascend[1]) {
      pphot->i2++;
      if(pphot->i2 <= pmcb->je)
        pphot->x[1] = pco->x2f(pphot->i2);
      else {
        pmcb->pbval->BoundaryFunction_[OUTER_X2](pmcb,pco,pphot);
        if (pphot->status == ESCAPED) {
          pphot->face = OUTER_X2;
        }
      }
    } else {
      pphot->i2--;
      if(pphot->i2 >= pmcb->js)
        pphot->x[1] = pco->x2f(pphot->i2+1);
      else {
        pmcb->pbval->BoundaryFunction_[INNER_X2](pmcb,pco,pphot);
        if (pphot->status == ESCAPED) {
          pphot->face = INNER_X2;
        }
      }
    } 
  }
  if ((face == 2) || (face == 4) || (face == 5) || (face == 6)) { //update x3 face
    if (ascend[2]) {
      pphot->i3++;
      if(pphot->i3 <= pmcb->ke)
        pphot->x[2] = pco->x3f(pphot->i3);
      else {
        pmcb->pbval->BoundaryFunction_[OUTER_X3](pmcb,pco,pphot);
        if (pphot->status == ESCAPED) {
          pphot->face = OUTER_X2;
        }
      }
    } else {
      pphot->i3--;
      if(pphot->i3 >= pmcb->ks)
        pphot->x[2] = pco->x3f(pphot->i3+1);
      else {
        pmcb->pbval->BoundaryFunction_[INNER_X3](pmcb,pco,pphot);
        if (pphot->status == ESCAPED) {
          pphot->face = INNER_X3;
        }
      }
    } 
  }

  // Update opacities
  if (pphot->status == EVOLVING) {
    pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot);
    pphot->sct_coef = pmcb->ScatteringOpacity(pmcb,pphot);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonMover::CartesianToCurvalinear(Photon *pphot)
//  \brief convert k vector from cartesian to curvalinear

void PhotonMover::CartesianToCurvalinear(Photon *pphot) {

  // Default corresponds to Cartesian so just copy
  for (int i=0; i<3; ++i)
    pphot->k[i] = pphot->kcart[i];
  
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonMover::CurvalinearToCartesian(Photon *pphot)
//  \brief convert k vector from curvalinear to cartesian

void PhotonMover::CurvalinearToCartesian(Photon *pphot) {

  // Default corresponds to Cartesian so just copy
  for (int i=0; i<3; ++i)
    pphot->kcart[i] = pphot->k[i];
  
}
