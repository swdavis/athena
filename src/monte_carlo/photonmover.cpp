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

// Array size for MRW dist computation
int nmax = 1000;

// Implementation of base class

PhotonMover::PhotonMover(MonteCarloBlock *pmcb) {

  pmy_mcb = pmcb;

  // MRW acceleration
  acceleration = pmcb->acceleration;
  lorentz_transform = pmcb->lorentz_transform;
  if (acceleration)
    InitializeMRWDist();

}

PhotonMover::~PhotonMover() {

  if (acceleration) {
    mrwprob.DeleteAthenaArray();
    mrwdev.DeleteAthenaArray();
  }
    
}

void PhotonMover::Move(Photon *pphot) {

}

//----------------------------------------------------------------------------------------
//! \fn bool PhotonMover::MRWAcceleration(Photon *pphot, MCRandom *pran, Real dist, Real)
//  \brief Accelerate photon diffusion with modified random walk method

bool PhotonMover::MRWAcceleration(Photon *pphot, MCRandom *pran, Real dist, Real tauacc) {

  MonteCarloBlock *pmcb = pmy_mcb;
  bool accel_success = true;
  //printf("a");
  // draw from path length distribution and reduce weight accordingly
  Real mrw = MRWDist(pran);          
  while (mrw <= 0.)
    mrw = MRWDist(pran); 
  Real delta = -log(mrw);
  Real chi = 3.*(pphot->abs_coef+pphot->sct_coef)/SQR(PI);
  Real ct,r0;
  Real vx=0.,vy=0.,vz=0.;
  if (lorentz_transform) {
    //printf("%d %d %d \n",pphot->i3,pphot->i2,pphot->i1);
    vx = pmcb->vel(0,pphot->i3,pphot->i2,pphot->i1);
    vy = pmcb->vel(1,pphot->i3,pphot->i2,pphot->i1);
    vz = pmcb->vel(2,pphot->i3,pphot->i2,pphot->i1);
    Real beta = sqrt(SQR(vx)+SQR(vy)+SQR(vz)) / 2.9979e10;
    Real gamma = 1./sqrt(1.-SQR(beta));
    beta *= gamma;
    vx *= gamma;
    vy *= gamma;
    vz *= gamma;
    //chi *= 1.+4.*SQR(beta);
    chi *= SQR(gamma)+SQR(beta)/3.;
    r0 = 0.5*(sqrt(1+4.*chi*dist*delta*beta)-1.)/(delta*beta*chi);
    //printf("%g %g %g %g %g\n",delta,r0,beta,chi,dist);
    if ((pphot->abs_coef+pphot->sct_coef) * r0 > tauacc) {
      ct = delta*SQR(r0)*chi;
      //printf("b ");
    } else {
      //printf("c ");
      //printf("c %g %g %g %g\n",r0,(pphot->abs_coef+pphot->sct_coef) * r0,dist,pco->x1f(3)-pco->x2f(2));
      // Reject delta and use standard update
      accel_success = false;
    }
  } else {
    ct = delta*SQR(dist)*chi;
    r0 = dist;
  }
  if (accel_success) {
    pphot->weight *= exp(-ct*pphot->abs_coef);
  
    // position packet on sphere of radius r0
    Real mu = 2.*pran->uniform()-1.0;
    Real stheta = sqrt(1.0-mu*mu);
    Real phi = 2.*PI*pran->uniform();
    ct /= 2.9979e10;
    pphot->x[0] += stheta*cos(phi) * r0 + vx * ct;
    pphot->x[1] += stheta*sin(phi) * r0 + vy * ct;
    pphot->x[2] += mu * r0 + vz * ct;
    
    UpdateZone(pphot);
    
    // Assume isotropic random direction
    mu = 2.*pran->uniform()-1.0;
    stheta = sqrt(1.0-mu*mu);
    phi = 2.*PI*pran->uniform();
    pphot->k[0] = stheta*cos(phi);
    pphot->k[1] = stheta*sin(phi);
    pphot->k[2] = mu;
  }
  return accel_success;

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

  // check for positiviity
  /*if (dx1 < 0.) {
    dx1 = HUGE_NUMBER;
    printf("Warning: dx1 < 0\n");
  }
  if (dx2 < 0.) {
    dx2 = HUGE_NUMBER;
    printf("Warning: dx2 < 0\n");
  }
  if (dx3 < 0.) {
    dx3 = HUGE_NUMBER;
    printf("Warning: dx3 < 0\n");
    }*/

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
//! \fn void PhotonMover::MovePhotonToNextZone()
//  \brief updates photon zone when face is known

void PhotonMover::MovePhotonToNextZone(Photon *pphot, MCCoord *pco,
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
          pphot->face = OUTER_X3;
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
    // Opacities need to be calculated using comoving frame energy and then transformed
    // back to Eulerian frame when Lorentz Transformations are enabled.
    Real shift;

    if (pmy_mcb->lorentz_transform) {
      // Shift photon energy to comoving frame
      shift = pmy_mcb->LorentzTransformFrequencyShift(pphot);
      pphot->energy *= shift;
      // compute opacities in comoving frame
      pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot);
      pphot->sct_coef = pmcb->ScatteringOpacity(pmcb,pphot);
      // Shift energy back to Eulerian frame
      pphot->energy /= shift;
      // Shift opaciteis to Eulerian frame
      pphot->abs_coef *= shift;
      pphot->sct_coef *= shift;
    } else {
      // No distinction between comovinng frame and eulerian frame
      pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot);
      pphot->sct_coef = pmcb->ScatteringOpacity(pmcb,pphot);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonMover::UpdateZone(photon *pphot)
//  \brief check/updates photon zone after displacement

void PhotonMover::UpdateZone(Photon *pphot) {

  bool change = false;
  MonteCarloBlock *pmcb = pmy_mcb;
  MCCoord *pco = pmcb->pcoord;
  bool update = false;

  if (pphot->x[0] >= pco->x1f(pphot->i1+1)) {
    while (pphot->x[0] >= pco->x1f(pphot->i1+1)) {
      pphot->i1++;
      if(pphot->i1 > pmcb->ie)
	pmcb->pbval->BoundaryFunction_[OUTER_X1](pmcb,pco,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = OUTER_X1;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  } else if (pphot->x[0] < pco->x1f(pphot->i1)) {
    while (pphot->x[0] < pco->x1f(pphot->i1)) {
      pphot->i1--;
      if(pphot->i1 < pmcb->is)
	pmcb->pbval->BoundaryFunction_[INNER_X1](pmcb,pco,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = INNER_X1;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  }
  if (pphot->x[1] >= pco->x2f(pphot->i2+1)) {
    while (pphot->x[1] >= pco->x2f(pphot->i2+1)) {
      pphot->i2++;
      if(pphot->i2 > pmcb->je)
	pmcb->pbval->BoundaryFunction_[OUTER_X2](pmcb,pco,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = OUTER_X2;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  } else if (pphot->x[1] < pco->x2f(pphot->i2)) {
    while (pphot->x[1] < pco->x2f(pphot->i2)) {
      pphot->i2--;
      if(pphot->i2 < pmcb->js)
	pmcb->pbval->BoundaryFunction_[INNER_X2](pmcb,pco,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = INNER_X2;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  }
  if (pphot->x[2] >= pco->x3f(pphot->i3+1)) {
    while (pphot->x[2] >= pco->x3f(pphot->i3+1)) {
      pphot->i3++;
      if(pphot->i3 > pmcb->ke)
	pmcb->pbval->BoundaryFunction_[OUTER_X3](pmcb,pco,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = OUTER_X3;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  } else if (pphot->x[2] < pco->x3f(pphot->i3)) {
    while (pphot->x[2] < pco->x3f(pphot->i3)) {
      pphot->i3--;
      if(pphot->i3 < pmcb->ks)
	pmcb->pbval->BoundaryFunction_[INNER_X3](pmcb,pco,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = INNER_X3;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  }

  // Update opacities
  if (pphot->status == EVOLVING) {
    // Opacities need to be calculated using comoving frame energy and then transformed
    // back to Eulerian frame when Lorentz Transformations are enabled.
    Real shift;

    if (lorentz_transform) {
      // Shift photon energy to comoving frame
      shift = pmcb->LorentzTransformFrequencyShift(pphot);
      pphot->energy *= shift;
      // compute opacities in comoving frame
      pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot);
      pphot->sct_coef = pmcb->ScatteringOpacity(pmcb,pphot);
      // Shift energy back to Eulerian frame
      pphot->energy /= shift;
      // Shift opaciteis to Eulerian frame
      pphot->abs_coef *= shift;
      pphot->sct_coef *= shift;
    } else {
      // No distinction between comovinng frame and eulerian frame
      pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot);
      pphot->sct_coef = pmcb->ScatteringOpacity(pmcb,pphot);
    }
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

//----------------------------------------------------------------------------------------
//! \fn void InitializeMWDist(void)
//  \brief initialize modified randon walk path length distribution

void PhotonMover::InitializeMRWDist(void)
{

  mrwprob.NewAthenaArray(nmax);
  mrwdev.NewAthenaArray(nmax);

  for(int i=0; i<nmax; ++i) {
    mrwdev(i) = static_cast<Real>(i)/static_cast<Real>(nmax-1);
    mrwprob(i) = 0.0;
  }
  for(int i=0; i<nmax-1; ++i) {
    int n = 1;
    Real sign = 1.0;
    Real yn2 = pow(mrwdev(i),n*n);
    // Compute the sum to the limit of Real precision
    while (yn2 > 1.e-17 ) {
      mrwprob(i) += sign * yn2;
      sign *= -1.0;
      n += 1;
      yn2 = pow(mrwdev(i),n*n);
    }
    mrwprob(i) *= 2.;
    if (mrwprob(i) > 1.0) mrwprob(i)=1.0;
  }
  mrwprob(nmax-1) = 1.0;


}

//----------------------------------------------------------------------------------------
//! \fn Real MRWDist(MCRandom *pran)
//  \brief get modified randon walk path length

Real PhotonMover::MRWDist(MCRandom *pran)
{

  Real x0 = pran->uniform();
 
  // Perform a binary search
  int low =0, high = nmax-1, mid;
  while(low<=high) {
    mid=(low+high)/2;
    if(mrwprob(mid-1) <= x0) {
      if(mrwprob(mid) > x0)
        break;
      else
        low=mid+1;
    }
    else
      high=mid-1;
  }

  // Replace binary search with initial guess ?
  //int i = static_cast<int>(x0*static_cast<Real>(nmax));

  // use linear interpolation to find location
  if (mid == 0)
    return mrwdev(0);
  else if (low == nmax) {
    return mrwdev(nmax-1);
  } else {
    Real slope = (x0 - mrwprob(mid-1)) / (mrwprob(mid) - mrwprob(mid-1));
    return mrwdev(mid-1)+(mrwdev(mid)-mrwdev(mid-1)) * slope;
  }


}

