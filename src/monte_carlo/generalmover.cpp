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
//#include "montecarlo.hpp"
#define MAXITER 1e8
//#define DEBUG
//#define OUTTEST_SP
//#define OUTTEST_GK
#define OUTTEST_OG
//#define OUTTEST_TF
#define NBUFFER 1000
#define NCOORD 4

//#define VERBOSE

// GR headers
#define n_array 4
#define tolerance 1.e-5
#define max_iteration 2
#define epsilon 1.e-40
#define tt IMC0
#define slope 1
#define R0 0

static Real dlambda;
static Real a;
static Real ri;

// Implementation of general photon mover

GeneralMover::GeneralMover(MonteCarloBlock *pmcb) 
  : PhotonMover(pmcb) {

  a = pmy_mcb->a;
  dlambda = pmy_mcb->stepsize;

}

GeneralMover::~GeneralMover() {

}

//----------------------------------------------------------------------------------------
//! \fn void GeneralMover::Move(Photon *pphot)
//  \brief Moves photon along straight line specified number of mean free paths or until
//         photon leave monte carlo block

void GeneralMover::Move(Photon *pphot) {

  FILE *file_output = fopen("output_gr.dat", "a");

  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  MCCoord *pco = pmy_mcb->pcoord;

  // get number of mean free paths photon will travel
  Real TauRemaining = GetOpticalDepth(pran);

  // References for momentum vectors
  Real& kx = pphot->kcart[0];
  Real& ky = pphot->kcart[1];
  Real& kz = pphot->kcart[2];
  Real& kr  = pphot->k[IMC1];
  Real& kth = pphot->k[IMC2];
  Real& kph = pphot->k[IMC3];

#ifdef DEBUG
  typedef struct {
    Real dl, dlr, dlt, dlp;
    Real cth, sth, cph, sph;
    Real kr, kth, kph;
    Real kx, ky, kz;
    Real x,y,z;
    int i,j,k;
    bool ascend[3];
  } debug_t;
  debug_t db[NBUFFER];
#endif

  Real r_outer = 1.0 + sqrt(1.0 - SQR(a)) + 1.0e-3;
  
  // Added by Shane
#ifdef OUTTEST_GK

  FILE *outfile1 = fopen("mccomp.in","a");
  FILE *outfile3 = fopen("mccomp.in3", "a");
  Real gcov0[NCOORD][NCOORD];
  Real kphi0, kphif, kt0, ktf, kth0, kthf, xIMC1_old, rf;
  int su,sui,tpr=0;
  Real kphi0_bl, kt0_bl, kth0_bl, delta;
  Real alpha, beta;
  Real e_const, l_const, q_const;

  delta = SQR(pphot->x[IMC1]) - 2. * pphot->x[IMC1] + SQR(a);
  pmcb->Metric(pphot->x,gcov0);

  kt0 = pphot->k[IMC0]*gcov0[IMC0][IMC0] + pphot->k[IMC1]*gcov0[IMC0][IMC1] + 
    pphot->k[IMC2]*gcov0[IMC0][IMC2] + pphot->k[IMC3]*gcov0[IMC0][IMC3];
  kth0 = pphot->k[IMC0]*gcov0[IMC2][IMC0] + pphot->k[IMC1]*gcov0[IMC2][IMC1] + 
    pphot->k[IMC2]*gcov0[IMC2][IMC2] + pphot->k[IMC3]*gcov0[IMC2][IMC3];
  kphi0 = pphot->k[IMC0]*gcov0[IMC3][IMC0] + pphot->k[IMC1]*gcov0[IMC3][IMC1] + 
    pphot->k[IMC2]*gcov0[IMC3][IMC2] + pphot->k[IMC3]*gcov0[IMC3][IMC3];
  
  e_const = -kt0;
  l_const = kphi0;
  q_const = SQR(kth0) + SQR(kphi0 * cos(pphot->x[IMC2])/sin(pphot->x[IMC2])) - 
    SQR(a * kt0 * cos(pphot->x[IMC2]));
  /*printf("Constants of motion before integration (E, l, Q): %g %g %g\n", 
    e_const, l_const, q_const);*/
  
  if (pmcb->kerrschild_flag) {
    Metric_BoyerLindquist(pphot->x, gcov0);
  }

  kt0_bl = (pphot->k[IMC0] - 2. * pphot->x[IMC1] / delta * pphot->k[IMC1])*gcov0[IMC0][IMC0] 
    + pphot->k[IMC1]*gcov0[IMC0][IMC1] + pphot->k[IMC2]*gcov0[IMC0][IMC2] 
    + (pphot->k[IMC3] - a / delta * pphot->k[IMC1])*gcov0[IMC0][IMC3];
  kth0_bl = (pphot->k[IMC0] - 2. * pphot->x[IMC1] / delta * pphot->k[IMC1])*gcov0[IMC2][IMC0]
    + pphot->k[IMC1]*gcov0[IMC2][IMC1] + pphot->k[IMC2]*gcov0[IMC2][IMC2] 
    + (pphot->k[IMC3] - a / delta * pphot->k[IMC1])*gcov0[IMC2][IMC3];
  kphi0_bl = (pphot->k[IMC0] - 2. * pphot->x[IMC1] / delta * pphot->k[IMC1])*gcov0[IMC3][IMC0]
    + pphot->k[IMC1]*gcov0[IMC3][IMC1] + pphot->k[IMC2]*gcov0[IMC3][IMC2] 
    + (pphot->k[IMC3] - a / delta * pphot->k[IMC1])*gcov0[IMC3][IMC3];

  ri = 0.;
  /*if ( (COORDINATE_SYSTEM == "spherical_polar") && (!pmcb->blackhole_flag) )
    ri = pphot->x[IMC1]; // Spherical Polar
    if (pmcb->blackhole_flag) ri = exp(pphot->x[IMC1]); // Kerr-Schild*/
  if (COORDINATE_SYSTEM == "spherical_polar")
    ri = pphot->x[IMC1];

  if (pmy_mcb->kerrschild_flag) {
    Real r = pphot->x[IMC1];
    Real a2 = SQR(a);
    Real phi_bl = pphot->x[IMC3] - a * 0.5 / sqrt(1. - a2) * (log((r - 1. - sqrt(1. - a2)) /
							    (r - 1. + sqrt(1. - a2)))
							- log((ri - 1. - sqrt(1. - a2)) / 
							      (ri - 1. + sqrt(1. - a2))));
    Real delta = SQR(pphot->x[IMC1]) - 2 * pphot->x[IMC1] + SQR(a);
    Real kt_bl = pphot->k[IMC0] - 2. * pphot->x[IMC1] / delta * pphot->k[IMC1];
    Real kph_bl = pphot->k[IMC3] - a / delta * pphot->k[IMC1];

    fprintf(file_output, "0 %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g\n", 
	    pphot->x[IMC0], pphot->x[IMC1], pphot->x[IMC2], phi_bl,
	    kt_bl, pphot->k[IMC1], pphot->k[IMC2], kph_bl);

  } else if (pmy_mcb->boyerlindquist_flag) {
    fprintf(file_output, "0 %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g\n",
	    pphot->x[IMC0], pphot->x[IMC1], pphot->x[IMC2], pphot->x[IMC3],
	    pphot->k[IMC0], pphot->k[IMC1], pphot->k[IMC2], pphot->k[IMC3]);

  } else {
    fprintf(file_output, "0 %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g\n",
	    pphot->x[IMC0], pphot->x[IMC1], pphot->x[IMC2], pphot->x[IMC3],
	    pphot->k[IMC0], pphot->k[IMC1], pphot->k[IMC2], pphot->k[IMC3]);
  }
#endif // #ifdef OUTTEST_GK


#ifdef OUTTEST_TF
  FILE *outtest_tf = fopen("outtest_tf.dat", "a");
  Real psi, coszeta, r0, kt0, kth0, rf, kthf, ktf;
  pphot->k[IMC0] *= -1;
  pphot->k[IMC1] *= -1;
  pphot->k[IMC2] *= -1;
  pphot->k[IMC3] *= -1;
  pmcb->TetradTransform(pphot, 1.0); // 1.0 = to comoving frame
  kt0 = pphot->k[IMC0];
  //printf("kt0: %g\n",kt0);
  pmcb->TetradTransform(pphot, -1.0); // -1.0 = to Eulerian frame
  pphot->k[IMC0] *= -1;
  pphot->k[IMC1] *= -1;
  pphot->k[IMC2] *= -1;
  pphot->k[IMC3] *= -1;

  pmcb->TetradTransform(pphot, 1.0); // 1.0 = to comoving frame

  r0 = pphot->x[IMC1];
  kth0 = pphot->k[IMC2];
  

  pmcb->TetradTransform(pphot, -1.0); // -1.0 = to Eulerian frame

#endif 

  // this should eventually get moved to photon initialization 
  pphot->dk[IMC0] = 0.;
  pphot->dk[IMC1] = 0.;
  pphot->dk[IMC2] = 0.;
  pphot->dk[IMC3] = 0.;

  Stepsize(pphot);

  int count = 0;
  int iter = 0;
  int zone_counter = 0;
  Real chi = pphot->sct_coef + pphot->abs_coef;
  chi = (chi > TINY_NUMBER) ? chi : TINY_NUMBER; // return max
#ifdef VERBOSE
  printf("Tau: %g; chi: %g; chi*dlambda: %g\n", TauRemaining, chi, chi*dlambda);
#endif

#ifdef OUTTEST_SP
  //printf("%g %g\n", pphot->x[IMC1], sin(pphot->x[IMC2]));
  pmcb->TetradTransform(pphot, 1.0);
  /*printf("tetrad frame kt: %g  |n|: %g\n", pphot->k[IMC0], 
    sqrt(SQR(pphot->k[IMC1]) + SQR(pphot->k[IMC2]) + SQR(pphot->k[IMC3])));*/
  CurvalinearToCartesian(pphot);
  Real rf_sp,thf,phf,dl0;
  FinalPositionSphericalPolarGR(pmcb,pco,pphot,rf_sp,thf,phf,dl0);
  pmcb->TetradTransform(pphot, -1.0);
  Real gcov[NCOORD][NCOORD];
  pmcb->Metric(pphot->x, gcov);
/*printf("coordinate frame kt: %g  |n|: %g\n", pphot->k[IMC0], 
	 sqrt(SQR(pphot->k[IMC1]) + SQR(pphot->x[IMC1] * pphot->k[IMC2]) 
	 + SQR(pphot->x[IMC1] * sin(pphot->x[IMC2]) * pphot->k[IMC3])));*/
  FILE *spres = fopen("spres.dat", "a");
  fprintf(spres, "%g %g %g\n", pphot->x[IMC1], pphot->x[IMC2], pphot->x[IMC3]);
  fprintf(spres, "%g %g %g\n\n", rf_sp, thf, phf);
#endif

  while ( (pphot->status == EVOLVING) && (iter < MAXITER) && (TauRemaining > 0.) ) {

    iter++;
    count++;

#ifdef OUTTEST_GK
    fprintf(file_output, "%d ", count);
#endif // #ifdef OUTTEST_GK

   if (TauRemaining > chi * step) {
     VerletStep(pphot);
   } else {
     step = TauRemaining / chi;
     VerletStep(pphot);
   }

   TauRemaining -= chi * step;

#ifdef OUTTEST_GK
   if (iter == 1) {
     if (COORDINATE_SYSTEM == "spherical_polar") {
       if (pphot->x[IMC1] > ri) sui = -1;
       else sui = 1;
       //printf("first step: %g %g %d\n",pphot->x[IMC1],xIMC1_old,sui);
     }
   } else {
     Real sustep;
     if (pphot->x[IMC1] > ri) 
       sustep = -1;
     else 
       sustep = 1;
     if (sui != sustep)
       tpr = 1;
   }
   
#endif

   if (COORDINATE_SYSTEM == "spherical_polar" && (pmcb->kerrschild_flag or pmcb->boyerlindquist_flag)) {
     if (pphot->x[IMC1] <= r_outer) {
       step = -fabs((pphot->x[IMC1] - r_outer) / pphot->k[IMC1]);
       for (int i = 0; i < NCOORD; i++) 
	 pphot->x[i] += pphot->k[i] * step;
       //pphot->x[IMC1] = r_outer + 1.0e-3;
       pphot->status = DESTROYED;
     }
   }

    // Check if photon changed zones
    if (UpdateZone(pphot)) {
      UpdateOpacities(pphot, pmcb);
      zone_counter++;
      chi = pphot->sct_coef + pphot->abs_coef;
      chi = (chi > TINY_NUMBER) ? chi : TINY_NUMBER; // return max(chi, TINY_NUMBER)
      } 

    // Update moments
    if (pmcb->moments_flag) {
      pmcb->UpdateMoments(pphot,step);
    }

    if ((isnan(pphot->k[IMC0])) or (pphot->IsNanPhoton())) {
      pphot->PrintPhoton();
      printf("pphot->x[IMC0]: %g  pphot->x[IMC1]: %g  pphot->x[IMC2]: %g  pphot->x[IMC3]: %g\n",
	     pphot->x[IMC0], pphot->x[IMC1], pphot->x[IMC2], pphot->x[IMC3]);
      printf("pphot->k[IMC0]: %g  pphot->k[IMC1]: %g  pphot->k[IMC2]: %g  pphot->k[IMC3]: %g\n",
	     pphot->k[IMC0], pphot->k[IMC1], pphot->k[IMC2], pphot->k[IMC3]);
      printf("iter: %d\n", iter);
      pphot->status = DESTROYED;
    }

    Stepsize(pphot);

#ifdef OUTTEST_TF

    Real rdisk = 1.0e10; // outer edge for thin optically thick disk
    if ((pphot->x[IMC2] >= (M_PI / 2.0)) and (pphot->x[IMC1] - r_outer > 1.0e-3)) { // photon has crossed plane for the first time
      if (pphot->x[IMC1] < rdisk) {
	step = -(pphot->x[IMC2] - M_PI/2.0) / pphot->k[IMC2];
	for (int i = 0; i < NCOORD; i++) 
	  pphot->x[i] += pphot->k[i] * step;      
	pphot->status = DESTROYED;
      }
    }
    
#endif
    // Perform any user work
    if (UserWorkInMove != NULL) UserWorkInMove(pmcb,pphot,this);

  } // end of photon integration

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


#ifdef OUTTEST_TF

  // transform k^alpha to k^(a) and print relevant quantites
  // if photon terminated because it hit the black hole, then transforming
  // into the tetrad frame will cause errors
  if (pphot->x[IMC1] - r_outer > 1.0e-3) {
    pphot->k[IMC0] *= -1;
    pphot->k[IMC1] *= -1;
    pphot->k[IMC2] *= -1;
    pphot->k[IMC3] *= -1;
    pmcb->TetradTransform(pphot, 1.0); // 1.0 = to comoving frame
    Real ktf = pphot->k[IMC0];
    //printf("ktf: %g\n",ktf);
    pmcb->TetradTransform(pphot, -1.0); // -1.0 = to Eulerian frame
    pphot->k[IMC0] *= -1;
    pphot->k[IMC1] *= -1;
    pphot->k[IMC2] *= -1;
    pphot->k[IMC3] *= -1;

    pmcb->TetradTransform(pphot, 1.0); // 1.0 = to comoving frame
    //fprintf(outtest_tf, "%15.10g %15.10g %15.10g %15.10g\n",
    //	    kth0, pphot->x[IMC1], pphot->k[IMC2], pphot->k[IMC0]);
    fprintf(outtest_tf, "%15.10g %15.10g %15.10g %15.10g\n",
	    kth0, pphot->x[IMC1], pphot->k[IMC2], kt0/ktf);
  } else { // don't tetrad transform if photon is inside the event horizon
    fprintf(outtest_tf, "%15.10g %15.10g %15.10g %15.10g\n",
	    kth0, pphot->x[IMC1], -1.0, -1.0);
  }
  
#endif

#ifdef VERBOSE
  printf("The photon crossed %d zones, traveling %g after %d iterations.\n", zone_counter, 
  	 dlambda*static_cast<Real>(iter), iter);
  printf("k[IMC0]: %g  k[IMC1]: %g  k[IMC2]: %g  k[IMC3]: %g \n",
  	 pphot->k[IMC0], pphot->k[IMC1], pphot->k[IMC2], pphot->k[IMC3]);
  printf("end GeneralMover::Move\n");
#endif

#ifdef OUTTEST_GK
  fprintf(file_output,"\n\n");
  rf = 0.;
  if (COORDINATE_SYSTEM == "spherical_polar") {
    rf = pphot->x[IMC1];
  }

  pmcb->Metric(pphot->x,gcov0);
  kphif = pphot->k[IMC0]*gcov0[IMC3][IMC0] + pphot->k[IMC1]*gcov0[IMC3][IMC1] + 
    pphot->k[IMC2]*gcov0[IMC3][IMC2] + pphot->k[IMC3]*gcov0[IMC3][IMC3];
  ktf = pphot->k[IMC0]*gcov0[IMC0][IMC0] + pphot->k[IMC1]*gcov0[IMC0][IMC1] + 
    pphot->k[IMC2]*gcov0[IMC0][IMC2] + pphot->k[IMC3]*gcov0[IMC0][IMC3];
  kthf = pphot->k[IMC0]*gcov0[IMC2][IMC0] + pphot->k[IMC1]*gcov0[IMC2][IMC1] +
    pphot->k[IMC2]*gcov0[IMC2][IMC2] + pphot->k[IMC3]*gcov0[IMC2][IMC3];
  
  e_const = -ktf;
  l_const = kphif;
  q_const = SQR(kthf) + SQR(kphif * cos(pphot->x[IMC2])/sin(pphot->x[IMC2])) - 
    SQR(a * ktf * cos(pphot->x[IMC2]));
  /*printf("Constants of motion after integration (E, l, Q): %g %g %g\n", 
    e_const, l_const, q_const);*/

  //printf("out: %g %g %g\n",kphi0,kt0,-kphi0/kt0);
  Real ui = 1./ri;
  Real uf = 1./rf;

  // alpha -kphi0/kt0/sin(theta)
  if (pmcb->kerrschild_flag) 
    alpha = -kphi0_bl / kt0_bl;
  else if (pmcb->boyerlindquist_flag)
    alpha = -kphi0/kt0;

  // beta =sqrt(SQR(kth0)/SQR(kt0)+a*cos^2(theta)-\alpha^2*cos^2(theta));
  // Assumes with start at theta = pi/2, then beta^2 = q^2
  if (pmcb->boyerlindquist_flag) {
    beta = sqrt(SQR(kth0)/SQR(kt0));
    if (kth0 > 0)
      beta = -beta;
  } else if (pmcb->kerrschild_flag) {
    beta = sqrt(SQR(kth0_bl)/SQR(kt0_bl));
    if (kth0_bl > 0)
      beta = -beta;
  }

  fprintf(outfile1, "%g %g %g %g %g %d %d\n", alpha, beta, ui, 0., uf, sui, tpr);
  if (su < 0) 
    fprintf(outfile3, "%g %g %g %g %g %d %d\n", -kphif/ktf, 0., uf, 0., ui, -su, tpr); 
  else if (su > 0) 
    fprintf(outfile3, "%g %g %g %g %g %d %d\n", kphif/ktf, 0., uf, 0., ui, -su, tpr); 
  fclose(outfile1);
  fclose(outfile3);
  
  FILE *outfile2 = fopen("angles.in", "a");
  if (rf < r_outer + 1.0e-3) fprintf(outfile2, "%g\n", pphot->x[IMC3]);
  else fprintf(outfile2, "%g\n", pphot->x[IMC3]);
  //fprintf(outfile2, "%g\n", pphot->x[IMC3]);
  fclose(outfile2);
#endif

#ifdef OUTTEST_SP
  if (pphot->status == ESCAPED) {
    Real xf = rf_sp*sin(thf)*cos(phf);
    Real yf = rf_sp*sin(thf)*sin(phf);
    Real zf = rf_sp*cos(thf);
    Real xp =  pphot->x[IMC1]*sin(pphot->x[IMC2])*cos(pphot->x[IMC3]);
    Real yp =  pphot->x[IMC1]*sin(pphot->x[IMC2])*sin(pphot->x[IMC3]);
    Real zp =  pphot->x[IMC1]*cos(pphot->x[IMC2]);
    Real delta = sqrt(SQR(xf-xp)+SQR(yf-yp)+SQR(zf-zp));
    Real dmax = 1.e-8*rf_sp;
    FILE *file_outtest = fopen("output_outtest.dat", "a");
    /*printf("%.3g %.3g %.3g %.3g %.3g %.3g\n", pphot->x[IMC1], pphot->x[IMC2],
	   pphot->x[IMC3], rf_sp, thf, phf);
	   printf("%.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3e\n", xp, yp, zp, xf, yf, zf, delta, dlambda);*/
    fprintf(file_outtest, "%.3g %.3g %.3g %.3g %.3g %.3g %.3g %.3e\n", xp, yp, zp, xf, yf, zf, delta, dlambda);
    fclose(file_outtest);
  }

#endif

#ifdef OUTTEST_TF
  fclose(outtest_tf);
#endif

  fclose(file_output);

}

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
  MCCoord *pco = pmcb->pcoord;
  bool update = false;

  if (pphot->x[IMC1] >= pco->x1f(pphot->i1+1)) {
    update = true;
    while (pphot->x[IMC1] >= pco->x1f(pphot->i1+1)) {
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
  } else if (pphot->x[IMC1] < pco->x1f(pphot->i1)) {
    update = true;
    while (pphot->x[IMC1] < pco->x1f(pphot->i1)) {
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

  if (pphot->x[IMC2] >= pco->x2f(pphot->i2+1)) {
    update = true;
    while (pphot->x[IMC2] >= pco->x2f(pphot->i2+1)) {
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
  } else if (pphot->x[IMC2] < pco->x2f(pphot->i2)) {
    update = true;
    while (pphot->x[IMC2] < pco->x2f(pphot->i2)) {
      pphot->i2--;
      if(pphot->i2 < pmcb->js)
	pmcb->pbval->BoundaryFunction_[INNER_X2](pmcb,pco,pphot);
      if (pphot->status == ESCAPED) {
	pphot->face = INNER_X2;
	break;
      }
      if (pphot->status == DESTROYED) {
	break;
      }
    }
  }

  if (pphot->x[IMC3] >= pco->x3f(pphot->i3+1)) {
    update = true;
    while (pphot->x[IMC3] >= pco->x3f(pphot->i3+1)) {
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
  } else if (pphot->x[IMC3] < pco->x3f(pphot->i3)) {
    update = true;
    while (pphot->x[IMC3] < pco->x3f(pphot->i3)) {
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
  // Returns true if zone changes, false otherwise
  return update;


}

//----------------------------------------------------------------------------------------
// GR functions


void GeneralMover::VerletStep(Photon *pphot) {
   
  Real gamma[NCOORD][NCOORD][NCOORD];
  Real k_n1[NCOORD],k_n1_copy[NCOORD];
  Real dk_n1[NCOORD];
  Real dk, error;
  int n_iteration, i, k, j;

  for (i=0;i<NCOORD;i++) {
    pphot->x[i] += (pphot->k[i])*step + 0.5*(pphot->dk[i])*SQR(step);
    k_n1[i] = (pphot->k[i]) + (pphot->dk[i])*step;
  }

  for (i=0;i<NCOORD;i++) {
    for (j=0;j<NCOORD;j++) {
      for (k=0;k<NCOORD;k++) {
	gamma[i][j][k] = 0.;
      }
    }
  }
  
  pmy_mcb->Connection(pphot->x, gamma);
  n_iteration = 0;
  
  do {
    
    n_iteration += 1;
    error = 0.;
    
    for (i=0;i<NCOORD;i++) {
      k_n1_copy[i] = k_n1[i];
    }
    
    for (k=0;k<NCOORD;k++) {  
      dk_n1[k] = 
	-2. * (k_n1_copy[IMC0] * 
	       (gamma[k][IMC0][IMC1] * k_n1_copy[IMC1] +
		gamma[k][IMC0][IMC2] * k_n1_copy[IMC2] +
		gamma[k][IMC0][IMC3] * k_n1_copy[IMC3])
	       +
	       k_n1_copy[IMC1] * (gamma[k][IMC1][IMC2] * k_n1_copy[IMC2] +
			          gamma[k][IMC1][IMC3] * k_n1_copy[IMC3]) +
	       k_n1_copy[IMC2] * gamma[k][IMC2][IMC3] * k_n1_copy[IMC3]);
     
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

  /*printf("%g %g %g %g %g %g %g %g %g\n", step, pphot->k[IMC0], pphot->k[IMC1],
    pphot->k[IMC2], pphot->k[IMC3], k_n1[IMC0], k_n1[IMC1], k_n1[IMC2], k_n1[IMC3]);*/

  // update photon energy due to evolving k_t (coordinate frame)
  pphot->energy *= k_n1[IMC0]/(pphot->k[IMC0]); 

  for (i=0;i<NCOORD;i++) {
    pphot->k[i] = k_n1[i];
    pphot->dk[i] = dk_n1[i];
  }
  
#ifdef OUTTEST_GK
  if (pmy_mcb->kerrschild_flag) {
    /*Real phi_bl = pphot->x[IMC3] - a * (atan((pphot->x[IMC1] - 1.) / sqrt(1. - SQR(a))) / 
					sqrt(1. - SQR(a)) - 
					atan((ri - 1.) / sqrt(1. - SQR(a))) / 
					sqrt(1. - SQR(a)));*/
    Real r = pphot->x[IMC1];
    Real a2 = SQR(a);
    Real phi_bl = pphot->x[IMC3] - a * 0.5 / sqrt(1. - a2) * (log((r - 1. - sqrt(1. - a2)) /
							    (r - 1. + sqrt(1. - a2)))
							- log((ri - 1. - sqrt(1. - a2)) / 
							      (ri - 1. + sqrt(1. - a2))));
    Real delta = SQR(pphot->x[IMC1]) - 2 * pphot->x[IMC1] + SQR(a);
    Real kt_bl = pphot->k[IMC0] - 2. * pphot->x[IMC1] / delta * pphot->k[IMC1];
    Real kph_bl = pphot->k[IMC3] - a / delta * pphot->k[IMC1];
    fprintf(file_output, "%5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g\n", 
	    pphot->x[IMC0], pphot->x[IMC1], pphot->x[IMC2], phi_bl,
	    kt_bl, pphot->k[IMC1], pphot->k[IMC2], kph_bl); // KS coords

  } else if (pmy_mcb->boyerlindquist_flag) {
    fprintf(file_output, "%5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g\n",
	    pphot->x[IMC0], pphot->x[IMC1], pphot->x[IMC2], pphot->x[IMC3],
	    pphot->k[IMC0], pphot->k[IMC1], pphot->k[IMC2], pphot->k[IMC3]); // BL coords

  } else {
    fprintf(file_output, "%5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g %5.5g\n",
	    pphot->x[IMC0], pphot->x[IMC1], pphot->x[IMC2], pphot->x[IMC3],
	    pphot->k[IMC0], pphot->k[IMC1], pphot->k[IMC2], pphot->k[IMC3]);
  }
#endif // #ifdef OUTTEST_GK

  return;

}

// return the stepsize based on the current zone and k-vector
// this should be updated with every iteration since k continuously changes
void GeneralMover::Stepsize(Photon *pphot) {

  if (!pphot->pmy_mcb->varystep_flag) {
    step = dlambda; // keep pphot->step constant
    return;
  }

  Real stepx1, stepx2, stepx3;
  Real kx1, kx2, kx3;
  Real small = 1.e-20;
  MCCoord *pco = pmy_mcb->pcoord;

  if (pphot->IsNanPhoton()) 
    return;
  
  kx1 = (fabs(pphot->k[IMC1]) > epsilon) ? fabs(pphot->k[IMC1]) : small; // prevents divide by 0
  kx2 = (fabs(pphot->k[IMC2]) > epsilon) ? fabs(pphot->k[IMC2]) : small;
  kx3 = (fabs(pphot->k[IMC3]) > epsilon) ? fabs(pphot->k[IMC3]) : small;
  /*printf("%g %g %g %g %g %g %g %g\n", pphot->k[IMC1], pphot->k[IMC2],
    pphot->k[IMC3], kx1, kx2, kx3, epsilon, small);*/

  stepx1 = ((pco->x1f(pphot->i1 + 1) - pco->x1f(pphot->i1)) / kx1) * dlambda;
  stepx2 = ((pco->x2f(pphot->i2 + 1) - pco->x2f(pphot->i2)) / kx2) * dlambda;
  stepx3 = ((pco->x3f(pphot->i3 + 1) - pco->x3f(pphot->i3)) / kx3) * dlambda;
  
  if (stepx1 < stepx2) {
    if (stepx1 < stepx3)
      step = stepx1;
    else 
      step = stepx3;
  } else if (stepx2 < stepx3) {
    step = stepx2;
  } else
    step = stepx3;
  
  return;

}


void Connect_KerrSchild(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]){

  Real r, r2, cth, sth, cth2, sth2, c2th, s2th;
  Real a2, sigma, sigma2, A, delta;
  void sincos(Real th, Real *sth, Real *cth);

  r = x[IMC1];
  r2 = SQR(r);
  sincos(x[IMC2], &sth, &cth);
  cth2 = SQR(cth);
  sth2 = SQR(sth);
  sincos(2. * x[IMC2], &s2th, &c2th);
  
  a2 = SQR(a);
  sigma = r2 + a2 * cth2;
  sigma2 = SQR(sigma);
  delta = r2 - 2 * r + a2;
  A = SQR(r2 + a2) - a2 * delta * sth2;

  gamma[IMC0][IMC0][IMC0] = -2. * r / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC0][IMC0][IMC1] = -1. / sigma * (1. + 2. * r / sigma) * (1. - 2. * r2 / sigma);
  gamma[IMC0][IMC0][IMC2] = -a2 * r * s2th / sigma2;
  gamma[IMC0][IMC0][IMC3] = 2. * a * r * sth2 / sigma2 * (1. - 2. * r2 / sigma);

  gamma[IMC0][IMC1][IMC0] = gamma[IMC0][IMC0][IMC1];
  gamma[IMC0][IMC1][IMC1] = -2. / sigma * (1. + r / sigma) * (1. - 2. * r2 / sigma);
  gamma[IMC0][IMC1][IMC2] = -a2 * r * s2th / sigma2;
  gamma[IMC0][IMC1][IMC3] = a * sth2 / sigma * (1. + 2. * r / sigma) * 
    (1. - 2. * r2 / sigma);

  gamma[IMC0][IMC2][IMC0] = gamma[IMC0][IMC0][IMC2];
  gamma[IMC0][IMC2][IMC1] = gamma[IMC0][IMC1][IMC2];
  gamma[IMC0][IMC2][IMC2] = -2. * r2 / sigma;
  gamma[IMC0][IMC2][IMC3] = a2 * a * r / sigma2 * sth2 * s2th;

  gamma[IMC0][IMC3][IMC0] = gamma[IMC0][IMC0][IMC3];
  gamma[IMC0][IMC3][IMC1] = gamma[IMC0][IMC1][IMC3];
  gamma[IMC0][IMC3][IMC2] = gamma[IMC0][IMC2][IMC3];
  gamma[IMC0][IMC3][IMC3] = -2. * r * sth2 / sigma * (r + a2 * sth2 / sigma * 
						      (1. - 2. * r2 / sigma));


  gamma[IMC1][IMC0][IMC0] = -delta / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC1][IMC0][IMC1] = 1. / sigma * (1. - 2. * r2 / sigma) * (1. - delta / sigma);
  gamma[IMC1][IMC0][IMC2] = 0.;
  gamma[IMC1][IMC0][IMC3] = a * delta * sth2 / sigma2 * (1. - 2. * r2 / sigma);

  gamma[IMC1][IMC1][IMC0] = gamma[IMC1][IMC0][IMC1];
  gamma[IMC1][IMC1][IMC1] = 1. / sigma * (1. - 2. * r2 / sigma) * (2. - delta / sigma);
  gamma[IMC1][IMC1][IMC2] = -a2 / (2. * sigma) * s2th;
  gamma[IMC1][IMC1][IMC3] = a / sigma * sth2 * (r - (1. - 2. * r2 / sigma) * 
						(1. - delta / sigma));

  gamma[IMC1][IMC2][IMC0] = gamma[IMC1][IMC0][IMC2];
  gamma[IMC1][IMC2][IMC1] = gamma[IMC1][IMC1][IMC2];
  gamma[IMC1][IMC2][IMC2] = -r * delta / sigma;
  gamma[IMC1][IMC2][IMC3] = 0.;
 
  gamma[IMC1][IMC3][IMC0] = gamma[IMC1][IMC0][IMC3];
  gamma[IMC1][IMC3][IMC1] = gamma[IMC1][IMC1][IMC3];
  gamma[IMC1][IMC3][IMC2] = gamma[IMC1][IMC2][IMC3];
  gamma[IMC1][IMC3][IMC3] = -delta / sigma * sth2 * (r + a2 * sth2 / sigma * 
						     (1. - 2. * r2 / sigma));
  

  gamma[IMC2][IMC0][IMC0] = -a2 * r * s2th / (sigma2 * sigma);
  gamma[IMC2][IMC0][IMC1] = -a2 * r * s2th / (sigma2 * sigma);
  gamma[IMC2][IMC0][IMC2] = 0.;
  gamma[IMC2][IMC0][IMC3] = a * r * (r2 + a2) * s2th / (sigma2 * sigma);

  gamma[IMC2][IMC1][IMC0] = gamma[IMC2][IMC0][IMC1];
  gamma[IMC2][IMC1][IMC1] = -a2 * r * s2th / (sigma2 * sigma);
  gamma[IMC2][IMC1][IMC2] = r / sigma;
  /*gamma[IMC2][IMC1][IMC3] = a / (2. * sigma) * (1. + 2. * r * (r2 + a2) / sigma2) * s2th;*/
  gamma[IMC2][IMC1][IMC3] = ((a * cth * sth) / (sigma2 * sigma)) * 
    (r2 * r * (r + 2.) + 2. * a2 * r * (r + 1.) * cth2 + a2 * a2 * cth2 * cth2 
    + 2. * a2 * r * sth2); // from Shane's notebook -- not equal to Takahashi+07 

  gamma[IMC2][IMC2][IMC0] = gamma[IMC2][IMC0][IMC2];
  gamma[IMC2][IMC2][IMC1] = gamma[IMC2][IMC1][IMC2];
  gamma[IMC2][IMC2][IMC2] = -a2 * s2th / (2. * sigma);
  gamma[IMC2][IMC2][IMC3] = 0.;

  gamma[IMC2][IMC3][IMC0] = gamma[IMC2][IMC0][IMC3];
  gamma[IMC2][IMC3][IMC1] = gamma[IMC2][IMC1][IMC3];
  gamma[IMC2][IMC3][IMC2] = gamma[IMC2][IMC2][IMC3];
  /*gamma[IMC2][IMC3][IMC3] = -s2th / (2. * sigma) * (delta + 2. * r * 
    SQR((r2 + a2) / sigma));*/
  gamma[IMC2][IMC3][IMC3] = -(cth * sth / (sigma2 * sigma)) * 
    (a2 * a2 * a2 * cth2 * cth2 * cth2 +
     cth2 * cth2 * (3. * a2 * a2 * r2 + a2 * a2 * a2 * sth2) + 
     cth2 * (3. * a2 * r2 * r2 + 2. * a2 * a2 * r2 * sth2) +
     r * (r2 * r2 * r + a2 * r2 * (r + 4.) * sth2 + 2. * a2 * a2 * sth2 * sth2 +
     a2 * a2 * s2th * s2th)); // frome Shane's notebook


  gamma[IMC3][IMC0][IMC0] = -a / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC3][IMC0][IMC1] = -a / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC3][IMC0][IMC2] = -2. * a * r / sigma2 * cth / sth;
  gamma[IMC3][IMC0][IMC3] = a2 * sth2 / sigma2 * (1. - 2. * r2 / sigma);

  gamma[IMC3][IMC1][IMC0] = gamma[IMC3][IMC0][IMC1];
  gamma[IMC3][IMC1][IMC1] = -a / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC3][IMC1][IMC2] = -a / sigma * (1. + 2. * r / sigma) * cth / sth;
  gamma[IMC3][IMC1][IMC3] = 1. / sigma * (r + a2 * sth2 / sigma * 
					  (1. - 2. * r2 / sigma));

  gamma[IMC3][IMC2][IMC0] = gamma[IMC3][IMC0][IMC2];
  gamma[IMC3][IMC2][IMC1] = gamma[IMC3][IMC1][IMC2];
  gamma[IMC3][IMC2][IMC2] = -a * r / sigma;
  //gamma[IMC3][IMC2][IMC3] = (1. + 2. * r / sigma * ((r2 + a2) / sigma - 1.)) * cth / sth;
  gamma[IMC3][IMC2][IMC3] = ((1. / 4.) * SQR(a2 + 2. * r2 + a2 * c2th) * cth / sth + 
    a2 * r * s2th) / (sigma2); // from Shane's notebook

  gamma[IMC3][IMC3][IMC0] = gamma[IMC3][IMC0][IMC3];
  gamma[IMC3][IMC3][IMC1] = gamma[IMC3][IMC1][IMC3];
  gamma[IMC3][IMC3][IMC2] = gamma[IMC3][IMC2][IMC3];
  gamma[IMC3][IMC3][IMC3] = -a / sigma * sth2 * (r + a2 * sth2 / sigma * 
						 (1. - 2. * r2 / sigma));

					      					      
  /*Real r1, r2, r3, r4, sx, cx;
  Real th, dthdx2, dthdx22, d2thdx22, sth, cth, sth2, cth2, sth4,
    cth4, s2th, c2th;
  Real a2, a3, a4, rho2, irho2, rho22, irho22, rho23, irho23,
    irho23_dthdx2;
  Real fac1, fac1_rho23, fac2, fac3, a2cth2, a2sth2, r1sth2,
    a4cth4;
  // required by broken math.h 
  void sincos(Real th, Real *sth, Real *cth);

  r1 = exp(x[IMC1]);
  r2 = r1 * r1;
  r3 = r2 * r1;
  r4 = r3 * r1;

  sincos(2. * M_PI * x[IMC2], &sx, &cx);

  // HARM-2D MKS 
  th = M_PI * x[IMC2] + 0.5 * (1 - slope) * sx;
  dthdx2 = M_PI * (1. + (1 - slope) * cx);
  d2thdx22 = -2. * M_PI * M_PI * (1 - slope) * sx;

  dthdx22 = dthdx2 * dthdx2;

  sincos(th, &sth, &cth);
  sth2 = sth * sth;
  r1sth2 = r1 * sth2;
  sth4 = sth2 * sth2;
  cth2 = cth * cth;
  cth4 = cth2 * cth2;
  s2th = 2. * sth * cth;
  c2th = 2 * cth2 - 1.;

  a2 = a * a;
  a2sth2 = a2 * sth2;
  a2cth2 = a2 * cth2;
  a3 = a2 * a;
  a4 = a3 * a;
  a4cth4 = a4 * cth4;

  rho2 = r2 + a2cth2;
  rho22 = rho2 * rho2;
  rho23 = rho22 * rho2;
  irho2 = 1. / rho2;
  irho22 = irho2 * irho2;
  irho23 = irho22 * irho2;
  irho23_dthdx2 = irho23 / dthdx2;

  fac1 = r2 - a2cth2;
  fac1_rho23 = fac1 * irho23;
  fac2 = a2 + 2 * r2 + a2 * c2th;
  fac3 = a2 + r1 * (-2. + r1);

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      for (int k = 0; k< NCOORD; k++) {
	gamma[i][j][k] = 0.;
      }
    }
  }

  gamma[IMC0][IMC0][IMC0] = 2. * r1 * fac1_rho23;
  gamma[IMC0][IMC0][IMC1] = r1 * (2. * r1 + rho2) * fac1_rho23;
  gamma[IMC0][IMC0][IMC2] = -a2 * r1 * s2th * dthdx2 * irho22;
  gamma[IMC0][IMC0][IMC3] = -2. * a * r1sth2 * fac1_rho23;

  gamma[IMC0][IMC1][IMC1] = 2. * r2 * (r4 + r1 * fac1 - a4cth4) * irho23;
  gamma[IMC0][IMC1][IMC2] = -a2 * r2 * s2th * dthdx2 * irho22;
  gamma[IMC0][IMC1][IMC3] =
    a * r1 * (-r1 * (r3 + 2 * fac1) + a4cth4) * sth2 * irho23;

  gamma[IMC0][IMC2][IMC2] = -2. * r2 * dthdx22 * irho2;
  gamma[IMC0][IMC2][IMC3] = a3 * r1sth2 * s2th * dthdx2 * irho22;

  gamma[IMC0][IMC3][IMC3] =
    2. * r1sth2 * (-r1 * rho22 + a2sth2 * fac1) * irho23;

  gamma[IMC1][IMC0][IMC0] = fac3 * fac1 / (r1 * rho23);
  gamma[IMC1][IMC0][IMC1] = fac1 * (-2. * r1 + a2sth2) * irho23;
  gamma[IMC1][IMC0][IMC2] = 0.;
  gamma[IMC1][IMC0][IMC3] = -a * sth2 * fac3 * fac1 / (r1 * rho23);

  gamma[IMC1][IMC1][IMC1] =
    (r4 * (-2. + r1) * (1. + r1) +
     a2 * (a2 * r1 * (1. + 3. * r1) * cth4 + a4cth4 * cth2 +
	   r3 * sth2 + r1 * cth2 * (2. * r1 + 3. * r3 -
				    a2sth2))) * irho23;
  gamma[IMC1][IMC1][IMC2] = -a2 * dthdx2 * s2th / fac2;
  gamma[IMC1][IMC1][IMC3] =
    a * sth2 * (a4 * r1 * cth4 + r2 * (2 * r1 + r3 - a2sth2) +
		a2cth2 * (2. * r1 * (-1. + r2) + a2sth2)) * irho23;

  gamma[IMC1][IMC2][IMC2] = -fac3 * dthdx22 * irho2;
  gamma[IMC1][IMC2][IMC3] = 0.;

  gamma[IMC1][IMC3][IMC3] =
    -fac3 * sth2 * (r1 * rho22 - a2 * fac1 * sth2) / (r1 * rho23);

  gamma[IMC2][IMC0][IMC0] = -a2 * r1 * s2th * irho23_dthdx2;
  gamma[IMC2][IMC0][IMC1] = r1 * gamma[IMC2][IMC0][IMC0];
  gamma[IMC2][IMC0][IMC2] = 0.;
  gamma[IMC2][IMC0][IMC3] = a * r1 * (a2 + r2) * s2th * irho23_dthdx2;

  gamma[IMC2][IMC1][IMC1] = r2 * gamma[IMC2][IMC0][IMC0];
  gamma[IMC2][IMC1][IMC2] = r2 * irho2;
  gamma[IMC2][IMC1][IMC3] =
    (a * r1 * cth * sth *
     (r3 * (2. + r1) +
      a2 * (2. * r1 * (1. + r1) * cth2 + a2 * cth4 +
	    2 * r1sth2))) * irho23_dthdx2;

  gamma[IMC2][IMC2][IMC2] =
    -a2 * cth * sth * dthdx2 * irho2 + d2thdx22 / dthdx2;
  gamma[IMC2][IMC2][IMC3] = 0.;

  gamma[IMC2][IMC3][IMC3] =
    -cth * sth * (rho23 +
		  a2sth2 * rho2 * (r1 * (4. + r1) + a2cth2) +
		  2. * r1 * a4 * sth4) * irho23_dthdx2;

  gamma[IMC3][IMC0][IMC0] = a * fac1_rho23;
  gamma[IMC3][IMC0][IMC1] = r1 * gamma[IMC3][IMC0][IMC0];
  gamma[IMC3][IMC0][IMC2] = -2. * a * r1 * cth * dthdx2 / (sth * rho22);
  gamma[IMC3][IMC0][IMC3] = -a2sth2 * fac1_rho23;

  gamma[IMC3][IMC1][IMC1] = a * r2 * fac1_rho23;
  gamma[IMC3][IMC1][IMC2] =
    -2 * a * r1 * (a2 + 2 * r1 * (2. + r1) +
		   a2 * c2th) * cth * dthdx2 / (sth * fac2 * fac2);
  gamma[IMC3][IMC1][IMC3] = r1 * (r1 * rho22 - a2sth2 * fac1) * irho23;

  gamma[IMC3][IMC2][IMC2] = -a * r1 * dthdx22 * irho2;
  gamma[IMC3][IMC2][IMC3] =
    dthdx2 * (0.25 * fac2 * fac2 * cth / sth +
	      a2 * r1 * s2th) * irho22;

  gamma[IMC3][IMC3][IMC3] = (-a * r1sth2 * rho22 + a3 * sth4 * fac1) * irho23;	*/

  return;

}

void Connect_BoyerLindquist(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]) {

  // equations for the connection coefficients come from Frutos-Alfaro et al. (2012)

  Real m, r, j, th, phi;
  Real sth, sth2, cth, cth2;
  Real rho2, delta, rs;
  Real a2, r2, rho4, rho6;
  
  void sincos(Real th, Real *sth, Real *cth);
  
  r = x[IMC1];
  th = x[IMC2];
  phi = x[IMC3];
  m = 1.0;
  j = a*m;

  sincos(th, &sth, &cth);
  cth2 = SQR(cth);
  sth2 = 1. - cth2;

  r2 = SQR(r);
  a2 = SQR(a);
  rho2 = r2 + a2*cth2;
  rs = 2.*m; 
  delta = r2 - rs*r + a2;

  rho4 = SQR(rho2);
  rho6 = rho4*rho2;

  // Real A = SQR(r2+a2) - a2*delta*sth2;

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      for (int k = 0; k < NCOORD; k++) {
	gamma[i][j][k] = 0.;
      }
    }
  }

  gamma[IMC0][IMC0][IMC1] = rs/(2.*rho4*delta)*(r2+a2)*(2.*r2-rho2);
  gamma[IMC0][IMC0][IMC2] = -2.*a*j*r/rho4*sth*cth;
  
  gamma[IMC0][IMC1][IMC3] = -j*sth2/(rho4*delta)*(rho2*(r2-a2)+2.*r2*(r2+a2));

  gamma[IMC0][IMC2][IMC3] = 2.*a2*j*r/rho4*cth*sth2*sth;

  gamma[IMC1][IMC0][IMC0] = rs*delta/(2.*rho6)*(2.*r2-rho2);
  gamma[IMC1][IMC0][IMC3] = -j*delta/rho6*(2.*r2-rho2)*sth2;
  
  gamma[IMC1][IMC1][IMC1] = 1./(rho2*delta)*(rho2*(rs/2.-r)+r*delta);
  gamma[IMC1][IMC1][IMC2] = -a2/rho2*sth*cth;
  
  gamma[IMC1][IMC2][IMC2] = -r*delta/rho2;

  gamma[IMC1][IMC3][IMC3] = -delta*sth2/rho6*(r*rho4-a*j*(2.*r2-rho2)*sth2);
  
  gamma[IMC2][IMC0][IMC0] = -2.*a*j*r/rho6*sth*cth;
  gamma[IMC2][IMC0][IMC3] = 2.*j*r/rho6*(r2+a2)*sth*cth;

  gamma[IMC2][IMC1][IMC1] = a2/(rho2*delta)*sth*cth;
  gamma[IMC2][IMC1][IMC2] = r/rho2;
 
  gamma[IMC2][IMC2][IMC2] = gamma[IMC1][IMC1][IMC2];

  gamma[IMC2][IMC3][IMC3] = -sth*cth/rho6*(rho4*delta+rs*r*SQR(r2+a2));
  //gamma[IMC2][IMC2][IMC3] = -sth*cth/rho6*(A*rho2+(r2+a2)*a2*r*rs*sth2);

  gamma[IMC3][IMC0][IMC1] = j/(rho4*delta)*(2.*r2-rho2);
  gamma[IMC3][IMC0][IMC2] = -2.*j*r*cth/(rho4*sth);

  gamma[IMC3][IMC1][IMC3] = 1./(rho4*delta)*(r*rho2*(rho2-rs*r)-a*j*sth2*(2.*r2-rho2));

  gamma[IMC3][IMC2][IMC3] = cth/(rho4*sth)*(rho4+2.*a*j*r*sth2);

  return;
  
}

void Connect_Cartesian(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]) {

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      for (int k = 0; k < NCOORD; k++) {
	gamma[i][j][k] = 0.;
      }
    }
  }

  return;

}

void Connect_SphericalPolar(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]) {

  void sincos(Real th, Real *sth, Real *cth);
  Real sin,cos;
  sincos(x[IMC2],&sin,&cos);

  for(int i = 0; i < NCOORD; i++) {
    for(int j = 0; j < NCOORD; j++) {
      for(int k = 0; k < NCOORD; k++) {
	gamma[i][j][k]=0;
      }
    }
  }

  gamma[IMC1][IMC2][IMC2] = -x[IMC1];
  gamma[IMC1][IMC3][IMC3] = -x[IMC1]*sin*sin;
  gamma[IMC2][IMC1][IMC2] = 1./x[IMC1];
  gamma[IMC2][IMC2][IMC1] = 1./x[IMC1];
  gamma[IMC2][IMC3][IMC3] = -sin*cos;
  gamma[IMC3][IMC1][IMC3] = 1./x[IMC1];
  gamma[IMC3][IMC2][IMC3] = cos/sin;
  gamma[IMC3][IMC3][IMC1] = 1./x[IMC1];
  gamma[IMC3][IMC3][IMC2] = cos/sin;

  return;

}

void Connect_Cylindrical(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]) {

  for(int i = 0; i < NCOORD; i++) {
    for(int j = 0; j < NCOORD; j++) {
      for(int k = 0; k < NCOORD; k++) {
	gamma[i][j][k]=0;
      }
    }
  }

  gamma[IMC1][IMC2][IMC2] = -x[IMC1];
  gamma[IMC2][IMC1][IMC2] = 1./x[IMC1];
  gamma[IMC2][IMC2][IMC1] = 1./x[IMC1];

  return;

}


void Metric_KerrSchild(Real x[NCOORD], Real gcov[NCOORD][NCOORD]){

  Real sth, cth, s2, rho2, sigma, A, delta, a2;
  Real r, th, sth2, cth2, r2;
  /* required by broken math.h */
  void sincos(Real th, Real *sth, Real *cth);

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      gcov[i][j] = 0.;
    }
  }

  r = x[IMC1];
  r2 = SQR(r);
  th = x[IMC2];
  sincos(th, &sth, &cth);
  cth2 = SQR(cth);
  sth2 = SQR(sth);
  a2 = SQR(a);

  sigma = r2 + a2 * cth2;
  delta = r2 - 2 * r + a2;
  A = SQR(r2 + a2) - a2 * delta * sth2;
  
  gcov[IMC0][IMC0] = -1. * (1. - 2. * r / sigma);
  gcov[IMC0][IMC1] = 2 * r / sigma;
  gcov[IMC0][IMC3] = -2. * a * r * sth2 / sigma;

  gcov[IMC1][IMC0] = gcov[IMC0][IMC1];
  gcov[IMC1][IMC1] = 1. + 2. * r / sigma;
  gcov[IMC1][IMC3] = -a * sth2 * (1. + 2. * r / sigma);

  gcov[IMC2][IMC2] = sigma;

  gcov[IMC3][IMC0] = gcov[IMC0][IMC3];
  gcov[IMC3][IMC1] = gcov[IMC1][IMC3];
  gcov[IMC3][IMC3] = A * sth2 / sigma;

  /*
  Real tfac, rfac, hfac, pfac;
  Real hslope = 1.;
  r = exp(x[IMC1]) + R0;
  th = M_PI * x[IMC2] + ((1. - hslope) / 2.) * sin(2. * M_PI * x[IMC2]);

  sincos(th, &sth, &cth);
  sth = fabs(sth) + 1.e-40;
  s2 = sth * sth;
  rho2 = r * r + a * a * cth * cth;

  // transformation for Kerr-Schild -> modified Kerr-Schild 
  tfac = 1.;
  rfac = r - R0;
  hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * x[IMC2]);
  pfac = 1.;

  gcov[IMC0][IMC0] = (-1. + 2. * r / rho2) * tfac * tfac;
  gcov[IMC0][IMC1] = (2. * r / rho2) * tfac * rfac;
  gcov[IMC0][IMC3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

  gcov[IMC1][IMC0] = gcov[IMC0][IMC1];
  gcov[IMC1][IMC1] = (1. + 2. * r / rho2) * rfac * rfac;
  gcov[IMC1][IMC3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

  gcov[IMC2][IMC2] = rho2 * hfac * hfac;

  gcov[IMC3][IMC0] = gcov[IMC0][IMC3];
  gcov[IMC3][IMC1] = gcov[IMC1][IMC3];
  gcov[IMC3][IMC3] =
    s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;*/

  return;

}

void Metric_KerrSchild_Up(Real x[NCOORD], Real gcon[NCOORD][NCOORD]) {
  
  // equations come from Takahasi (2007) Appendix

  Real sth, cth, s2, rho2, sigma, delta, a2;
  Real r, th, sth2, cth2, r2;
  /* required by broken math.h */
  void sincos(Real th, Real *sth, Real *cth);

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      gcon[i][j] = 0.;
    }
  }

  r = x[IMC1];
  r2 = SQR(r);
  th = x[IMC2];
  sincos(th, &sth, &cth);
  cth2 = SQR(cth);
  sth2 = SQR(sth);
  a2 = SQR(a);

  sigma = r2 + a2 * cth2;
  delta = r2 - 2 * r + a2;

  gcon[IMC0][IMC0] = -(1. + (2. * r / sigma));
  gcon[IMC0][IMC1] = 2. * r / sigma;

  gcon[IMC1][IMC0] = gcon[IMC0][IMC1];
  gcon[IMC1][IMC1] = delta / sigma;
  gcon[IMC1][IMC3] = a / sigma;

  gcon[IMC2][IMC2] = 1. / sigma;
  
  gcon[IMC3][IMC1] = gcon[IMC1][IMC3];
  gcon[IMC3][IMC3] = 1. / (sigma * sth2);

  return;

}

void Metric_BoyerLindquist(Real x[NCOORD], Real gcov[NCOORD][NCOORD]) { 

  // equation for the Metric comes from the inside cover of Hartle
  
  Real m, r, j, th, phi;
  Real sth, sth2, cth, cth2;
  Real r2, a2;
  Real rho2, delta;

  void sincos(Real th, Real *sth, Real *cth);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      gcov[i][j] = 0;
    }
  }

  r = x[IMC1];
  th = x[IMC2];
  phi = x[IMC3];
  m = 1.0;
  j = a*m;

  sincos(th, &sth, &cth);
  cth2 = SQR(cth);
  sth2 = SQR(sth);

  r2 = SQR(r);
  a2 = SQR(a);
  rho2 = r2 + a2*cth2;
  delta = r2 - 2. * m* r + a2;

  gcov[IMC0][IMC0] = (-1. + 2.*m*r/rho2);
  gcov[IMC1][IMC1] = rho2/delta;
  gcov[IMC2][IMC2] = rho2;
  gcov[IMC3][IMC3] = (r2 + a2 + 2.*m*r*a2*sth2/rho2)*sth2;

  gcov[IMC0][IMC3] = -2.*m*a*r*sth2/rho2;
  gcov[IMC3][IMC0] = gcov[IMC0][IMC3];

  return;

}

void Metric_BoyerLindquist_Up(Real x[NCOORD], Real gcon[NCOORD][NCOORD]) {

  // Equation comes from ColinsCosmos.com/wiki/boyer-lindquist-coordinates, which
  // sites Frolov & Novikov Section D.1 (but I don't have access to this book)

  Real r, th;
  Real sth, sth2, cth, cth2;
  Real r2, a2;
  Real rho2, delta;

  void sincos(Real th, Real *sth, Real *cth);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      gcon[i][j] = 0;
    }
  }

  r = x[IMC1];
  th = x[IMC2];

  sincos(th, &sth, &cth);
  cth2 = SQR(cth);
  sth2 = SQR(sth);

  r2 = SQR(r);
  a2 = SQR(a);
  rho2 = r2 + a2 * cth2;
  delta = r2 - 2. * r + a2;

  gcon[IMC0][IMC0] = -1. / delta * (r2 + a2 + 2. * r * a2 * sth2 / rho2);
  gcon[IMC1][IMC1] = delta / rho2;
  gcon[IMC2][IMC2] = 1. / rho2;
  gcon[IMC3][IMC3] = (delta - a2 * sth2) / (rho2 * delta * sth2);

  gcon[IMC0][IMC3] = -2. * r * a / (rho2 * delta);
  gcon[IMC3][IMC0] = gcon[IMC0][IMC3];

  return;

}

void Metric_Cartesian(Real x[NCOORD], Real gcov[NCOORD][NCOORD]) {

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      if (i == j) {
	if (i == IMC0)
	  gcov[i][i] = -1.;
	else
	  gcov[i][i] = 1.;
      } else
	gcov[i][j] = 0.;
    }
  }

  return;

}	  

void Metric_SphericalPolar(Real x[NCOORD], Real gcov[NCOORD][NCOORD]) {

  int m, n;
  Real sin, cos;
  void sincos(Real t, Real *s, Real *c);
  sincos(x[IMC2], &sin, &cos);
  for (m = 0; m < 4; m++) {
    for (n = 0; n < 4; n++) {
      gcov[m][n] = 0;
    }
  }
  gcov[IMC0][IMC0] = -1;
  gcov[IMC1][IMC1] = 1;
  gcov[IMC2][IMC2] = x[IMC1] * x[IMC1];
  gcov[IMC3][IMC3] = x[IMC1] * x[IMC1] * sin * sin;

  return;

}


void Metric_Cylindrical(Real x[NCOORD], Real gcov[NCOORD][NCOORD]) {
  int m, n;
  for (m = 0; m < 4; m++) {
    for (n = 0; n < 4; n++) {
      gcov[m][n] = 0;
    }
  }
  gcov[IMC0][IMC0] = -1;
  gcov[IMC1][IMC1] = 1;
  gcov[IMC2][IMC2] = x[IMC1] * x[IMC1];
  gcov[IMC3][IMC3] = 1;

  return;

}


// Given initial position x^alpha, alpha, beta, determine the initial photon direction
// Uses alpha, beta from Cunningham & Bardeen (1973)
void GetMCDirection(Photon *pphot, Real alpha, Real beta) {
  
  MonteCarloBlock *pmcb = pphot->pmy_mcb;
  Real kcon[NCOORD], kcov[NCOORD]; // kcon = k^alpha; kcov = k_alpha
  // kcov is set by alpha, beta, but kcon is what is integrated 
  Real gcon[NCOORD][NCOORD], gcov[NCOORD][NCOORD];
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  Real ucon[NCOORD], bcon[NCOORD];

  // calculate g_{alpha,beta}, g^{a,b} 
  if (pmcb->kerrschild_flag) {
    Metric_KerrSchild(pphot->x, gcov);
    Metric_KerrSchild_Up(pphot->x, gcon);
  } else if (pmcb->boyerlindquist_flag) {
    Metric_BoyerLindquist(pphot->x, gcov);
    Metric_BoyerLindquist_Up(pphot->x, gcon);
  }

  // assumes initial radius r_0 >> alpha, beta 
  kcov[IMC0] = 1.0; // The goal is to have k^t [coordinate frame] to be 1.0 at emission
  kcov[IMC3] = alpha * kcov[IMC0] * sin(pphot->x[IMC2]);
  kcov[IMC2] = beta * kcov[IMC0];

  // use null geodesic equation to determine k_r: k_a * g^(a,b) * k_b = 0
  // solve the quadratic equation k_r^2 * gamma + k_r * zeta + xi = 0
  Real gamma = gcon[IMC1][IMC1];
  Real zeta = kcov[IMC0] * 2. * (gcon[IMC1][IMC0] + gcon[IMC1][IMC2] * beta + 
			    gcon[IMC1][IMC3] * alpha * sin(pphot->x[IMC2]));
  Real xi = SQR(kcov[IMC0]) * (gcon[IMC0][IMC0] + 2. * gcon[IMC0][IMC2] * beta + 
			       2. * gcon[IMC0][IMC3] * alpha * sin(pphot->x[IMC2]) 
			       + gcon[IMC2][IMC2] * SQR(beta)
			       + gcon[IMC3][IMC3] * SQR(alpha) * SQR(sin(pphot->x[IMC2])) + 
			       2. * gcon[IMC2][IMC3] * alpha * beta * sin(pphot->x[IMC2]));
  Real sqrtdis = sqrt(SQR(zeta) - 4. * gamma * xi);
  Real plus = (-zeta + sqrtdis) / (2. * gamma);
  Real minus = (-zeta - sqrtdis) / (2. * gamma);
  //printf("plus: %g  minus %g\n", plus, minus);

  // choose the positive root since dividing by k^t (which is < 0) switches sign
  if ((fabs(plus) < 1.0e-20) && (fabs(minus) < 1.0e-20)) {
    printf("WARNING: both roots of k_r = 0: %g  %g\n", plus, minus);
    pphot->status = DESTROYED;
  } else if ((plus < 0) && (minus > 0)) {
    kcov[IMC1] = minus;
  } else if ((plus > 0) && (minus < 0)) {
    kcov[IMC1] = plus;
  } else if ((plus < 0) && (minus < 0)) {
    printf("WARNING: both roots of k_r < 0: %g  %g\n", plus, minus);
    pphot->status = DESTROYED;
  } else if ((plus > 0) && (minus > 0)) {
    printf("WARNING: both roots of k_r > 0: %g  %g\n", plus, minus);
    pphot->status = DESTROYED;
  }
  // now k_alpha is set, so raise to k^alpha

  //printf("k_r: %g  k_th: %g  k_ph: %g\n", kcov[IMC1], kcov[IMC2], kcov[IMC3]);

  CovToCon(kcov, kcon, gcon); // converts k_alpha to k^alpha
  // normalize k^alpha such that k^t = 1
  for (int i = 0; i < NCOORD; i++) 
    pphot->k[i] = kcon[i] / kcon[IMC0];

  /*
  // set ucon, bcon for constructing the tetrad
  // ucon must be the same as in MonteCarloBlock::TetradTransform, which should be
  // based on the spinning fluid frame (though far from the black hole so close to flat).
  // This should get removed once MonteCarloBlock::RayTrace() is functional
  Real r = pphot->x[IMC1];
  Real omega = 1.0/(pow(r, 3./2.) + a); // circular velocity 

  ucon[IMC0] = sqrt(-1.0/(gcov[IMC0][IMC0] + 2.*gcov[IMC0][IMC3]*omega +
			  SQR(omega)*gcov[IMC3][IMC3]));
  ucon[IMC1] = 0.;
  ucon[IMC2] = 0.;
  ucon[IMC3] = (ucon[IMC0])*omega; 

  ConstructTetrad(ucon, bcon, gcov, econ, ecov);
  
  
  CoordinateToTetrad(kcon, pphot->k, ecov); // converts k^alpha to k^(a)*/

  // pphot->k now properly set and ready for integration
  
  /*
  printf("kcov: %g %g %g %g\n", kcov[IMC0], kcov[IMC1], kcov[IMC2], kcov[IMC3]);
  printf("kcon: %g %g %g %g\n", kcon[IMC0], kcon[IMC1], kcon[IMC2], kcon[IMC3]);
  printf("pphot->k: %g %g %g %g\n", pphot->k[IMC0], pphot->k[IMC1], pphot->k[IMC2],
  pphot->k[IMC3]); */

  return;

}
  
  
