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
  
  // SWD: Commented out to avoid compilation error
  //if (pmcb->kerrschild_flag) {
  //  Metric_BoyerLindquist(pphot->x, gcov0);
  //}

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

  // SWD: This shoudl be removed
  for (i=0;i<NCOORD;i++) {
    for (j=0;j<NCOORD;j++) {
      for (k=0;k<NCOORD;k++) {
	gamma[i][j][k] = 0.;
      }
    }
  }
  
  pmy_mcb->pcoord->Connect(pphot->x, gamma);
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


