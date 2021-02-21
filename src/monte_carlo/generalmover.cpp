//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file grmover.cpp
//  \brief implementation for moving photons via integration for curvalinear coordinates
//         and general relativity via a metric and connection
// Athena++ headers
#include "montecarlo.hpp"
#include "photon.hpp"
#include "photonmover.hpp"
#include "../mesh/mesh.hpp"
#include "debug.hpp"


// SWD: remove all of these
#define MAXITER 1e8
#define NBUFFER 1000
#define NCOORD 4
//#define VERBOSE

// GR headers
#define tolerance 1.e-5
#define max_iteration 2

// Implementation of general photon mover

GeneralMover::GeneralMover(MonteCarloBlock *pmcb) 
  : PhotonMover(pmcb) {

  step_par = pmy_mcb->stepsize;

}

GeneralMover::~GeneralMover() {

}

//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::Move(Photon *pphot)
//  \brief Moves photon along straight line specified number of mean free paths or until
//         photon leave monte carlo block

void GeneralMover::Move(Photon *pphot) {

  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  PhotonTrajectoryList *ptraj = pmy_mcb->ptraj;

  // get number of mean free paths photon will travel
  Real TauRemaining = GetOpticalDepth(pran);
 
#ifdef DEBUG
  typedef struct {
    Real dl, dlr, dlt, dlp;
    Real cth, sth, cph, sph;
    Real kr, kth, kph;
    Real x,y,z;
    int i,j,k;
    bool ascend[3];
  } debug_t;
  debug_t db[NBUFFER];
#endif

  Real step = StepSize(pphot);
  int count = 0;
  int iter = 0;
  int zone_counter = 0;
  Real chi = pphot->sct_coef + pphot->abs_coef;
  chi = (chi > TINY_NUMBER) ? chi : TINY_NUMBER; // return max

#ifdef VERBOSE
  printf("Tau: %g; chi: %g; chi*dlambda: %g\n", TauRemaining, chi, chi*dlambda);
#endif

  Real cth, sth, cph, sph;
  Real x,y,z,x0,y0,z0;
  cth = cos(pphot->x[1]);
  sth = sqrt(1. - SQR(cth));
  cph = cos(pphot->x[2]);
  sph = sin(pphot->x[2]);
  x0 = pphot->x[0]*sth*cph;
  y0 = pphot->x[0]*sth*sph;
  z0 = pphot->x[0]*cth;

  while ( (pphot->status == EVOLVING) && (iter < MAXITER) && (TauRemaining > 0.) ) {

    iter++;
    count++;

    /*Real cth = cos(pphot->x[1]);
    Real sth = sqrt(1. - SQR(cth));
    Real cph = cos(pphot->x[2]);
    Real sph = sin(pphot->x[2]);
    //pphot->PrintPhoton();
    // Compute cartesian
    pphot->kcart[0] = pphot->k[0]*sth*cph + pphot->k[1]*cth*cph - pphot->k[2]*sph;
    pphot->kcart[1] = pphot->k[0]*sth*sph + pphot->k[1]*cth*sph + pphot->k[2]*cph;
    pphot->kcart[2] = pphot->k[0]*cth - pphot->k[1]*sth;*/
  
   if (TauRemaining > chi * step) {
     VerletStep(pphot,step);
     if (pmy_mcb->pmy_mc->polarized)
       PropogatePolarization(pphot,step);
   } else {
     step = TauRemaining / chi;
     VerletStep(pphot,step);
     if (pmy_mcb->pmy_mc->polarized)
       PropogatePolarization(pphot,step);
   }
   /*cth = cos(pphot->x[1]);
   sth = sqrt(1. - SQR(cth));
   cph = cos(pphot->x[2]);
   sph = sin(pphot->x[2]);*/
   // Compute cartesian
   /*Real kx = pphot->k[0]*sth*cph + pphot->k[1]*cth*cph - pphot->k[2]*sph;
   Real ky  = pphot->k[0]*sth*sph + pphot->k[1]*cth*sph + pphot->k[2]*cph;
   Real kz = pphot->k[0]*cth - pphot->k[1]*sth;
   x = pphot->x[0]*sth*cph;
   y = pphot->x[0]*sth*sph;
   z = pphot->x[0]*cth;*/
   //printf("ks: %g %g %g %g %g %g\n",pphot->kcart[0],kx,pphot->kcart[1],ky,pphot->kcart[2],kz);
   /*x0 += step*pphot->kcart[0];
   y0 += step*pphot->kcart[1];
   z0 += step*pphot->kcart[2];*/
   //printf("xs: %g %g %g %g %g %g\n",x0,x,y0,y,z0,z,step);
   TauRemaining -= chi * step;

   // Check if photon changed zones
   if (UpdateZone(pphot)) {
     UpdateOpacities(pphot, pmcb);
     zone_counter++;
     chi = pphot->sct_coef + pphot->abs_coef;
     chi = (chi > TINY_NUMBER) ? chi : TINY_NUMBER; // return max(chi, TINY_NUMBER)
   } 
   if (pphot->status == DESTROYED) {
     pphot->PrintPhoton();
   }
   // Update moments
   if (pmcb->moments_flag) {
     pmcb->UpdateMoments(pphot,step);
   }
   
   if ((isnan(pphot->k[IMC0])) or (pphot->IsNanPhoton())) {
     pphot->PrintPhoton();
     pphot->status = DESTROYED;
   }

    step = StepSize(pphot);

    // Perform any user work
    if (UserWorkInMove != NULL) UserWorkInMove(pmcb,pphot,this);
    // SWD: put here for now, may need additional flag
    if (ptraj != NULL) ptraj->AddToTrajectory(pphot);
  } // end of photon integration

  // SWD: Try to remove this
  CurvalinearToCartesian(pphot);

  /*if (pphot->status == ESCAPED) {
    pphot->energy *= pphot->k[IMC0];
    //pphot->PrintPhoton();
    }*/

  if (iter >= MAXITER) {
    printf("GeneralMover::Move warning: max iter reached\n");
    pphot->PrintPhoton();
    pphot->status = DESTROYED;
  }



#ifdef VERBOSE
  printf("The photon crossed %d zones, traveling %g after %d iterations.\n", zone_counter, 
  	 dlambda*static_cast<Real>(iter), iter);
  printf("k[IMC0]: %g  k[IMC1]: %g  k[IMC2]: %g  k[IMC3]: %g \n",
  	 pphot->k[IMC0], pphot->k[IMC1], pphot->k[IMC2], pphot->k[IMC3]);
  printf("end GeneralMover::Move\n");
#endif

}

// SWD: The conversion from Cartesian to Curvalinear should be removed entirely
// or handled by MC coordinate class.
//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::CartesianToCurvalinear(Photon *pphot)
//  \brief convert k vector from cartesian to curvalinear

void GeneralMover::CartesianToCurvalinear(Photon *pphot) {

  Real cth = cos(pphot->x[1]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[2]);
  Real sph = sin(pphot->x[2]);
  // Compute spherical-polar
  pphot->k[0] = pphot->kcart[0]*sth*cph + pphot->kcart[1]*sth*sph + pphot->kcart[2]*cth;
  pphot->k[1] = pphot->kcart[0]*cth*cph + pphot->kcart[1]*cth*sph - pphot->kcart[2]*sth;
  pphot->k[2] = -pphot->kcart[0]*sph + pphot->kcart[1]*cph;
  
}


//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::CurvalinearToCartesian(Photon *pphot)
//  \brief convert k vector from curvalinear to cartesian

void GeneralMover::CurvalinearToCartesian(Photon *pphot) {

  Real cth = cos(pphot->x[IMC2]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[IMC3]);
  Real sph = sin(pphot->x[IMC3]);
  // Compute cartesian
  pphot->kcart[0] = (pphot->k[IMC1]*sth*cph + pphot->k[IMC2]*cth*cph - pphot->k[IMC3]*sph) / 
    pphot->k[IMC0];
  pphot->kcart[1] = (pphot->k[IMC1]*sth*sph + pphot->k[IMC2]*cth*sph + pphot->k[IMC3]*cph) / 
    pphot->k[IMC0];
  pphot->kcart[2] = (pphot->k[IMC1]*cth - pphot->k[IMC2]*sth) / pphot->k[IMC0];
 
}

//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb) 
//  \brief update opacities after a photon has changed zones
 
void GeneralMover::UpdateOpacities(Photon *pphot, MonteCarloBlock *pmcb) {

  pmy_mcb = pmcb;

    if (pphot->status == EVOLVING) { 
    // Opacities need to be calculated using comoving frame energy and then transformed 
    // back to Eulerian frame when Lorentz Transformations are enabled.
    Real shift;     
        
    /*if (pmy_mcb->orthotet_flag) {
      /*pmcb->TetradTransform(pphot, 1.0); // to comving frame
      pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot);    
      pphot->sct_coef = pmcb->ScatteringOpacity(pmcb,pphot); 
      pmcb->TetradTransform(pphot, -1.0); // to coordinate frame
      shift = pmy_mcb->TetradTransformFrequencyShift(pphot);
      pphot->energy *= shift;
      pphot->abs_coef = pmcb->AbsorptionOpacity(pmcb, pphot);
      pphot->sct_coef = pmcb->ScatteringOpacity(pmcb, pphot);
      pphot->energy /= shift;
      pphot->abs_coef *= shift;
      pphot->sct_coef *= shift;
    } else if (pmy_mcb->lorentz_transform) {  */  
    if (pmcb->boosts) {
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
//! \fn bool GeneralMover::UpdateZone(Photon *pphot)
//  \brief check if photon has changed zones and update zone indices

bool GeneralMover::UpdateZone(Photon *pphot) {

  bool change = false;
  MonteCarloBlock *pmcb = pmy_mcb;
  bool update = false;

  if (pphot->x[IMC1] >= pcoord->x1f(pphot->i1+1)) {
    update = true;
    while (pphot->x[IMC1] >= pcoord->x1f(pphot->i1+1)) {
      pphot->i1++;
      if(pphot->i1 > pmcb->ie)
	pmcb->pbval->BoundaryFunction_[OUTER_X1](pmcb,pcoord,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = OUTER_X1;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  } else if (pphot->x[IMC1] < pcoord->x1f(pphot->i1)) {
    update = true;
    while (pphot->x[IMC1] < pcoord->x1f(pphot->i1)) {
      pphot->i1--;
      if(pphot->i1 < pmcb->is)
	pmcb->pbval->BoundaryFunction_[INNER_X1](pmcb,pcoord,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = INNER_X1;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  }

  if (pphot->x[IMC2] >= pcoord->x2f(pphot->i2+1)) {
    update = true;
    while (pphot->x[IMC2] >= pcoord->x2f(pphot->i2+1)) {
      pphot->i2++;
      if(pphot->i2 > pmcb->je)
	pmcb->pbval->BoundaryFunction_[OUTER_X2](pmcb,pcoord,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = OUTER_X2;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  } else if (pphot->x[IMC2] < pcoord->x2f(pphot->i2)) {
    update = true;
    while (pphot->x[IMC2] < pcoord->x2f(pphot->i2)) {
      pphot->i2--;
      if(pphot->i2 < pmcb->js)
	pmcb->pbval->BoundaryFunction_[INNER_X2](pmcb,pcoord,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = INNER_X2;
	break;
      }
      if (pphot->status == DESTROYED) {
	break;
      }
    }
  }

  if (pphot->x[IMC3] >= pcoord->x3f(pphot->i3+1)) {
    update = true;
    while (pphot->x[IMC3] >= pcoord->x3f(pphot->i3+1)) {
      pphot->i3++;
      if(pphot->i3 > pmcb->ke)
	pmcb->pbval->BoundaryFunction_[OUTER_X3](pmcb,pcoord,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = OUTER_X3;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  } else if (pphot->x[IMC3] < pcoord->x3f(pphot->i3)) {
    update = true;
    while (pphot->x[IMC3] < pcoord->x3f(pphot->i3)) {
      pphot->i3--;
      if(pphot->i3 < pmcb->ks)
	pmcb->pbval->BoundaryFunction_[INNER_X3](pmcb,pcoord,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = INNER_X3;
	break;
      }
      if (pphot->status == DESTROYED)
	break;
    }
  }
  // Returns true if zone changes, false otherwise
  return update;


}



void GeneralMover::VerletStep(Photon *pphot, Real step) {
   
  Real gamma[NCOORD][NCOORD][NCOORD];
  Real k_n1[NCOORD],k_n1_copy[NCOORD];
  Real dk_n1[NCOORD];
  Real dk, error;
  int n_iteration, i, k, j;

  for (i=0;i<NCOORD;i++) {
    pphot->x[i] += (pphot->k[i])*step + 0.5*(pphot->dk[i])*SQR(step);
    k_n1[i] = (pphot->k[i]) + (pphot->dk[i])*step;
  }
  
  pcoord->Connect(pphot->x, gamma);
  n_iteration = 0;
  
  do {
    
    n_iteration += 1;
    error = 0.;
    
    for (i=0;i<NCOORD;i++) {
      k_n1_copy[i] = k_n1[i];
    }

    for (k=0;k<NCOORD;k++) {  
      // off diagonal elements
      dk_n1[k] = 
	-2. * (k_n1_copy[IMC0] * 
	       (gamma[k][IMC0][IMC1] * k_n1_copy[IMC1] +
		gamma[k][IMC0][IMC2] * k_n1_copy[IMC2] +
		gamma[k][IMC0][IMC3] * k_n1_copy[IMC3])
	       +
	       k_n1_copy[IMC1] * (gamma[k][IMC1][IMC2] * k_n1_copy[IMC2] +
			          gamma[k][IMC1][IMC3] * k_n1_copy[IMC3]) +
	       k_n1_copy[IMC2] * gamma[k][IMC2][IMC3] * k_n1_copy[IMC3]);
      // diagonal elements
      dk_n1[k] -= 
	(gamma[k][IMC0][IMC0] * k_n1_copy[IMC0] * k_n1_copy[IMC0] +
	 gamma[k][IMC1][IMC1] * k_n1_copy[IMC1] * k_n1_copy[IMC1] +
	 gamma[k][IMC2][IMC2] * k_n1_copy[IMC2] * k_n1_copy[IMC2] +
	 gamma[k][IMC3][IMC3] * k_n1_copy[IMC3] * k_n1_copy[IMC3]
	 );
           
      k_n1[k] = pphot->k[k] + (1./2.)*((pphot->dk[k])+dk_n1[k])*step;
      
      error += fabs(k_n1_copy[k] - k_n1[k]) / (k_n1[k]);
    }
  } while ((error > tolerance) && (n_iteration < max_iteration));

  //printf("%g %g %g %g %g %g %g %g %g\n", step, pphot->k[IMC0], pphot->k[IMC1],
  // pphot->k[IMC2], pphot->k[IMC3], k_n1[IMC0], k_n1[IMC1], k_n1[IMC2], k_n1[IMC3]);

  // update photon energy due to evolving k_t (coordinate frame)
  pphot->energy *= k_n1[IMC0]/(pphot->k[IMC0]); 
  
  for (i=0;i<NCOORD;i++) {
    pphot->k[i] = k_n1[i];
    pphot->dk[i] = dk_n1[i];
  }

}


void GeneralMover::PropogatePolarization(Photon *pphot, Real step) {

  Real gamma[NCOORD][NCOORD][NCOORD];
  // Store gamma in Coord to prevent recalculation
  pcoord->Connect(pphot->x, gamma);

  int i, j, k, l;
  std::complex<Real> ptcopy[NCOORD][NCOORD];

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      ptcopy[i][j] = pphot->polten[i][j];

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
	for (int l = 0; l < 4; l++) {
          pphot->polten[i][j] += -(gamma[i][k][l] * ptcopy[k][j] +
          			   gamma[j][k][l] * ptcopy[i][k]) * 
           pphot->k[l] * step;
        }}
      //printf("%d %d %g %g %g %g %g\n",i,j,pphot->polten[i][j].real(),pphot->polten[i][j].real()-ptcopy[i][j].real(),pphot->k[0],step,pphot->x[0]);//,gamma[0][i][j],gamma[1][i][j],,gamma[2][i][j],,gamma[3][i][j]);
    }}
}


// return the stepsize based on the current zone and k-vector
// this should be updated with every iteration since k continuously changes
Real GeneralMover::StepSize(Photon *pphot) {

  if (!pphot->pmy_mcb->varystep_flag) {
    return step_par; // keep step constant
  }

  Real small = 1.e-20;
  
  Real kx1 = (fabs(pphot->k[IMC1]) > small) ? fabs(pphot->k[IMC1]) : small; 
  Real kx2 = (fabs(pphot->k[IMC2]) > small) ? fabs(pphot->k[IMC2]) : small;
  Real kx3 = (fabs(pphot->k[IMC3]) > small) ? fabs(pphot->k[IMC3]) : small;

  // SWD: May want to store as dx1, etc.
  Real stepx1 = ((pcoord->x1f(pphot->i1 + 1) - pcoord->x1f(pphot->i1)) / kx1);
  Real stepx2 = ((pcoord->x2f(pphot->i2 + 1) - pcoord->x2f(pphot->i2)) / kx2);
  Real stepx3 = ((pcoord->x3f(pphot->i3 + 1) - pcoord->x3f(pphot->i3)) / kx3);

  Real step = (stepx1 < stepx2) ? stepx1 : stepx2;
  step = (step < stepx3) ? step : stepx3;

  return step*step_par;
}


