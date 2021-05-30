//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mc_geoimg.cpp
//  \brief Problem generator for creating an image with geodesics in kerr spacteime
//
//========================================================================================

#include <iostream> // temporary for testing

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../monte_carlo/montecarlo.hpp"
#include "../monte_carlo/photon.hpp"
#include "../monte_carlo/photonmover.hpp"
#include "../globals.hpp"

static int iphot;
static int nrays;
static int nalpha, nbeta;
static Real *alpha,*beta;
static Real rcam,thcam,phcam; 
static Real r_outer;
static int plane_cross;
static bool backward_integration;
static FILE *input;

static Real spsi,cpsi,szet,czet;

#if !MONTE_CARLO_STATIC
#error "This problem requires monte carlo"
#endif

// User function definitions
void MidplaneCrossing(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover);
void GetMCDirection(Photon *pphot, Real alpha, Real beta);
void GetDirectionTetrad(Photon *pphot, Real alpha, Real beta);
void TransformPhotonAtDisk(MonteCarloBlock *pmcb, Photon *pphot);
void TransformPhotonAtGridEdge(MonteCarloBlock *pmcb, Photon *pphot);

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief monte carlo test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real rideal = 8.314e7;
  Real rho = 1.;
  Real temp = 1.;
  Real gamma = peos->GetGamma(); 


  // Set nominal values for grid, unused
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = rideal*rho*temp/(gamma-1.0);
      }
    }
  }

}

void MonteCarloBlock::InitUserMonteCarloBlockData(ParameterInput *pin){

#ifdef MPI_PARALLEL 
  // Set iphot based on assumption that rays are distributed evenly
  // accross active processes
  int rank = Globals::my_rank;
  int ntot = pin->GetInteger("montecarlo", "nphot");
  if (rank > 0) {
    int nranks = Globals::nranks;
    int myn = ntot/(nranks-1);
    iphot = (rank-1)*myn;
    //if (rank == nranks-1)
    //  myn += ntot % (nranks-1);
    printf("iphot: %d %d %d\n",rank,iphot,myn);
    //nphremain = cadence = myn;
  }
#else
  iphot = 0;
#endif

}

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){

  nuser_var = 4;

  nrays = pin->GetInteger("montecarlo", "nphot");
  nalpha = nbeta = static_cast<int>(sqrt(static_cast<Real>(nrays)));
  backward_integration = pin->GetOrAddBoolean("problem","backward",false);

  // Set r_outer
  Real abh = pcoord->GetSpin();
  Real mbh = pcoord->GetMass();
  r_outer = 1.0 + sqrt(1.0 - SQR(abh));

  Real alpha_min = pin->GetOrAddReal("problem", "alpha_min", -10.);
  Real alpha_max = pin->GetOrAddReal("problem", "alpha_max", 10.);
  Real alpha_range = (alpha_max - alpha_min);
  Real beta_min = pin->GetOrAddReal("problem", "beta_min", -10.);
  Real beta_max = pin->GetOrAddReal("problem", "beta_max", 10.);
  Real beta_range = (beta_max - beta_min);
  alpha = new Real[nalpha];
  beta = new Real[nbeta];
  if (nalpha == 1) 
    alpha[0] = alpha_min;
  else {
    for (int i = 0; i<nalpha; i++) {
      alpha[i] = alpha_min+(static_cast<Real>(i)/static_cast<Real>(nalpha-1))*alpha_range;
      if (fabs(alpha[i]) < 1.e-5) alpha[i] = 1.e-1;
    }
  }
  if (nbeta == 1)
    beta[0] = beta_min;
  else {
    for (int i = 0; i<nbeta; i++) {
      beta[i] = beta_min+(static_cast<Real>(i)/static_cast<Real>(nbeta-1))*beta_range;
      if (fabs(beta[i]) < 1.0e-5) beta[i] = 1.0e-1;
    }
  }
  rcam = pin->GetOrAddReal("problem", "rcam", pin->GetReal("mesh","x1max"));
  thcam = pin->GetOrAddReal("problem", "thcam", 45.) * M_PI / 180.;
  phcam = pin->GetOrAddReal("problem", "phcam", 90.) * M_PI / 180.;
  
  EnrollUserWorkInMove(MidplaneCrossing);
}

void MonteCarloBlock::InitializePhoton(Photon *pphot) {

  MCCoord *pco = pphot->pmy_mcb->pcoord;

  // Set status flag
  pphot->status = EVOLVING;

  // Emit photons from a large radius r >> 1 in units of [GM/c^2]. Ideally this would be 
  // at infinity (or 1/r = 0), but r ~ 1e3 should be fine as long as the region of 
  // interest is small compared to this initial distance. For -10 < alpha, beta < 10, the
  // small angle approximation should still be fine. 
  // To avoid issues with the initial position being on a boundary, I push the position
  // by a small epislon from the ideal starting position.

  Real r, theta, phi, kt, kr, kth, kphi, alpha0, beta0;
  int ia,ib;

  r = rcam * 0.999999;
  theta = thcam;
  phi = phcam + 1.0e-3;
  // Set alpha, beta by either iterating over photons
  ia = iphot / nalpha;
  ib = iphot % nalpha;
  alpha0 = alpha[ia];
  beta0 = beta[ib];
  
  pphot->x[IMC0] = 0.0;
  pphot->x[IMC1] = r;
  pphot->x[IMC2] = theta;
  pphot->x[IMC3] = phi;

  // update the photon's zone indices
  pphot->i1 = -1;
  for(int i=pphot->pmy_mcb->is; i<=pphot->pmy_mcb->ie; i++) {
    if ((pphot->x[IMC1] > pco->x1f(i)) && (pphot->x[IMC1] <= pco->x1f(i+1)))
      pphot->i1 = i;
  }
  if (pphot->i1 < 0) pphot->weight = -1.0;
 
  pphot->i2 = -1;
  for(int i=pphot->pmy_mcb->js; i<=pphot->pmy_mcb->je; i++) {
    if ((pphot->x[IMC2] > pco->x2f(i)) && (pphot->x[IMC2] <= pco->x2f(i+1)))
      pphot->i2 = i;
  }
  if (pphot->i2 < 0) pphot->weight = -1.0;
 
  pphot->i3 = -1;
  for(int i=pphot->pmy_mcb->ks; i<=pphot->pmy_mcb->ke; i++) {
    if ((pphot->x[IMC3] > pco->x3f(i)) && (pphot->x[IMC3] <= pco->x3f(i+1)))
      pphot->i3 = i;
  }
  if (pphot->i3 < 0) pphot->weight = -1.0;

  if (pphot->weight < 0)
    pphot->status = DESTROYED;

  pphot->weight = 1.;

  pphot->user_var[0] = alpha0;
  pphot->user_var[1] = beta0;

  printf("nrays: %d  iphot %d  nalpha: %d  ialpha: %d ",nrays, iphot, nalpha, ia);
  printf("ibeta: %d  alpha: %g beta: %g\n",ib, alpha0, beta0);

  iphot++;

  // Set the initial photon direction using alpha, beta and the position
 
  // Initialize Stokes vector as unpolarized
  pphot->stokes[0] = 1.0;
  pphot->stokes[1] = 0.0;
  pphot->stokes[2] = 0.0;
  pphot->stokes[3] = 0.0;
 
  GetDirectionTetrad(pphot, alpha0, beta0);
  //printf("ktet: %e %e %e %e\n",pphot->k[IMC0],pphot->k[IMC1],pphot->k[IMC2],pphot->k[IMC3]);   
  
  pphot->energy = pphot->k[IMC0];
  pphot->weight = 1.0;
  for (int i=0; i<4; i++)
    pphot->dk[i] = 0.;
  //pphot->PrintPhoton();

  // Set plane crossing flag to zero
  plane_cross = 0;

  if (pphot->weight < 0.0) pphot->status = DESTROYED;

  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = 0.;
  pphot->sct_coef = 0.;
  
}

void MonteCarloBlock::FinalizePhoton(Photon *pphot) {
  
  
  if (pphot->status == DESTROYED) {
    pphot->status = ESCAPED;
    return;
  }
 
  if (backward_integration) {
    TransformPhotonAtGridEdge(this,pphot);
  } else {
    TransformPhotonAtDisk(this,pphot);
  } 
}

void TransformPhotonAtDisk(MonteCarloBlock *pmcb, Photon *pphot) {

  // If r is inside ISCO, do not attempt to transform and instead mark
  // as escaped
  if (pphot->x[IMC1] < r_outer + 1.0e-5) {
    pphot->status = ESCAPED;
    pphot->user_var[2] = 0.;
    pphot->user_var[3] = r_outer;
    return;
  }

  // SWD may not need full tetrad here unless presribed angular dependence
  // to intensity

  // Construct the orthonormal tetrad in comoving frame of circular orbit

  Real gcov[NCOORD][NCOORD];
  pmcb->pcoord->Metric(pphot->x, gcov);

  Real abh = pmcb->pcoord->GetSpin();
  Real mbh = pmcb->pcoord->GetMass();
  Real r = pphot->x[IMC1];
  Real omega = pow(mbh,0.5)/(pow(r, 3./2.) + abh*pow(mbh,3./2.)); // circular velocity 
  Real ucon[NCOORD];
  ucon[IMC0] = sqrt(-1.0/(gcov[IMC0][IMC0] + 2.*gcov[IMC0][IMC3]*omega +
                            SQR(omega)*gcov[IMC3][IMC3]));
  ucon[IMC1] = 0.;
  ucon[IMC2] = 0.;
  ucon[IMC3] = (ucon[IMC0])*omega;

  // create tetrad basis
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  ConstructTetrad(ucon, gcov, econ, ecov);

  // Reverse photon direction to get properties of photon that was emitted
  for (int i = 0; i < NCOORD; i++)
    pphot->k[i] *= -1.;

  //  Transform to comoving tetrad
  Real kcopy[NCOORD];
  for (int i = 0; i < NCOORD; i++)
    kcopy[i] = pphot->k[i];
  CoordinateToTetrad(kcopy, pphot->k, ecov);
  pphot->user_var[2] = pphot->k[IMC0];
    // Get radius at disk crossing
  pphot->user_var[3] = pphot->x[IMC1];

}

void TransformPhotonAtGridEdge(MonteCarloBlock *pmcb, Photon *pphot) {

  // Construct the orthonormal tetrad at edge of simulation grid
  Real ucon[NCOORD];
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  Real gcov[NCOORD][NCOORD], gcon[NCOORD][NCOORD];
  ucon[IMC0] = 1.;
  ucon[IMC1] = 0.;
  ucon[IMC2] = 0.;
  ucon[IMC3] = 0.;
    
  // create tetrad basis
  pmcb->pcoord->Metric(pphot->x, gcov);
  pmcb->pcoord->InverseMetric(pphot->x,gcon);
  Real wcon[NCOORD] = {0,1.,0.,0.}; // Q=1 points along projected BH symmetry axis 
  Real vcov[NCOORD] = {1.,0.,0.,1.};// Make image center point away from origin
  Real vcon[NCOORD]; 
      
  CovToCon(vcov,vcon,gcon);
  //printf("vcon: %e %e %e %e\n",vcon[IMC0],vcon[IMC1],vcon[IMC2],vcon[IMC3]);
  //printf("wcon: %e %e %e %e\n",wcon[IMC0],wcon[IMC1],wcon[IMC2],wcon[IMC3]);
   
  //ConstructTetrad(ucon, pphot->k, gcov, econ, ecov);
  //ConstructTetrad(ucon, pphot->k, vcon, gcov, econ, ecov);
  ConstructTetrad(ucon, vcon, wcon, gcov, econ, ecov);
  std::complex<Real> tcopy[NCOORD][NCOORD];
  
  ComplexCoordinateToTetrad(pphot->polten,tcopy,ecov);
  TensorToStokes(tcopy,pphot->stokes);
 
#ifdef DEBUG
  printf("final: %d\n",iphot-1);
  Real cosi = cos(pphot->x[IMC2]);
  Real pI = SQR(cpsi)+SQR(spsi*cosi);
  Real pQ = SQR(cpsi)-SQR(spsi*cosi);
  Real pU = 2*spsi*cpsi*cosi;
  Real cphi, sphi;
  cphi = cos(pphot->x[IMC3]);
  sphi = sin(pphot->x[IMC3]);
  Real sth, cth, r;
  r = pphot->x[IMC1];
  cth = cos(pphot->x[IMC2]);
  sth = sin(pphot->x[IMC2]);
  printf("prod: %e %e %e %e %e\n",spsi*cphi-cpsi*sphi,cphi,sphi,cpsi,spsi);
  Real Nrr = 2.*SQR(sth)*SQR(spsi*cphi-cpsi*sphi);
  Real Nrth = 2*sth*cth/r*SQR(spsi*cphi-cpsi*sphi);
  Real Nrph = 2./r*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
  Real Nthth = 2.*SQR(cth/r)*SQR(spsi*cphi-cpsi*sphi);
  Real Nthph = -2.*cth/SQR(r)/sth*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
  Real Nphph = 2./SQR(r*sth)*SQR(sphi*spsi+cphi*cpsi);
  printf("Nr: %e %e %e\n",Nrr,Nrth,Nrph);
  printf("Nth: %e %e\n",Nthth,Nthph);
  printf("Nph: %e\n",Nphph);
  Real kr = szet*sth*(cpsi*cphi+spsi*sphi)+czet*cth;
  Real kth = szet*cth/r*(cpsi*cphi+spsi*sphi)-czet*sth/r;
  Real kph = szet/(r*sth)*(spsi*cphi-cpsi*sphi);
  printf("k: %e %e %e\n",kr,kth,kph);
  printf("K: %e %e %e\n",pphot->k[IMC1],pphot->k[IMC2],pphot->k[IMC3]);
  printf("psif: %e %e %e %e %e\n",pI,pQ,pU,pQ/pI,pU/pI);
  printf("x: %g %g %g %g %g %g\n",pphot->x[IMC0],pphot->x[IMC1],pphot->x[IMC2],pphot->x[IMC3],pphot->x[IMC1]*cos(pphot->x[IMC3]),pphot->x[IMC1]*sin(pphot->x[IMC3]));
  //printf("pphot: %e %e %e %e %e\n",pphot->polten[IMC1][IMC1].real(),pphot->polten[IMC1][IMC2].real(),pphot->polten[IMC2][IMC1].real(),
  //         pphot->polten[IMC2][IMC2].real(),pphot->polten[IMC3][IMC3].real());
  printf("tcord[IMC0]: %e %e %e %e\n", pphot->polten[IMC0][IMC0].real(), pphot->polten[IMC0][IMC1].real(), pphot->polten[IMC0][IMC2].real(), pphot->polten[IMC0][IMC3].real());
  printf("tcord[IMC1]: %e %e %e %e\n", pphot->polten[IMC1][IMC0].real(), pphot->polten[IMC1][IMC1].real(), pphot->polten[IMC1][IMC2].real(), pphot->polten[IMC1][IMC3].real());
  printf("tcord[IMC2]: %e %e %e %e\n", pphot->polten[IMC2][IMC0].real(), pphot->polten[IMC2][IMC1].real(), pphot->polten[IMC2][IMC2].real(), pphot->polten[IMC2][IMC3].real());
  printf("tcord[IMC3]: %e %e %e %e\n", pphot->polten[IMC3][IMC0].real(), pphot->polten[IMC3][IMC1].real(), pphot->polten[IMC3][IMC2].real(), pphot->polten[IMC3][IMC3].real());
  printf("ttet[IMC0]: %e %e %e %e\n", tcopy[IMC0][IMC0].real(), tcopy[IMC0][IMC1].real(), tcopy[IMC0][IMC2].real(), tcopy[IMC0][IMC3].real());
  printf("ttet[IMC1]: %e %e %e %e\n", tcopy[IMC1][IMC0].real(), tcopy[IMC1][IMC1].real(), tcopy[IMC1][IMC2].real(), tcopy[IMC1][IMC3].real());
  printf("ttet[IMC2]: %e %e %e %e\n", tcopy[IMC2][IMC0].real(), tcopy[IMC2][IMC1].real(), tcopy[IMC2][IMC2].real(), tcopy[IMC2][IMC3].real());
  printf("ttet[IMC3]: %e %e %e %e\n", tcopy[IMC3][IMC0].real(), tcopy[IMC3][IMC1].real(), tcopy[IMC3][IMC2].real(), tcopy[IMC3][IMC3].real());
  //printf("ttet[11]: %e \n",tcopy[IMC1][IMC1].real());
  Real cosd = spsi*sphi+cpsi*cphi;
  Real sind = spsi*cphi-cpsi*sphi;
  Real I = SQR(cosd)+SQR(sind*cth);
  Real Q = SQR(cosd)-SQR(sind*cth);
  Real U = 2*cth*cosd*sind;
  printf("stopre: %e %e\n",Q/I,U/I);
  printf("stokes: %e %e\n",pphot->stokes[1],pphot->stokes[2]);
  printf("ecov[IMC0]: %e %e %e %e\n", ecov[IMC0][IMC0], ecov[IMC0][IMC1], ecov[IMC0][IMC2], ecov[IMC0][IMC3]);
  printf("ecov[IMC1]: %e %e %e %e\n", ecov[IMC1][IMC0], ecov[IMC1][IMC1], ecov[IMC1][IMC2], ecov[IMC1][IMC3]);
  printf("ecov[IMC2]: %e %e %e %e\n", ecov[IMC2][IMC0], ecov[IMC2][IMC1], ecov[IMC2][IMC2], ecov[IMC2][IMC3]);
  printf("ecov[IMC3]: %e %e %e %e\n", ecov[IMC3][IMC0], ecov[IMC3][IMC1], ecov[IMC3][IMC2], ecov[IMC3][IMC3]);
  printf("gcon: %e %e %e %e\n",gcon[IMC0][IMC0],gcon[IMC1][IMC1],gcon[IMC2][IMC2],gcon[IMC3][IMC3]);
#endif
    
  //printf("x: %g %g %g %g\n",pphot->x[IMC0],pphot->x[IMC1],pphot->x[IMC2],pphot->x[IMC3]);

}

//void LaunchPhotonAtDisk(MonteCarloBlock *pmcb, Photon *pphot) {
//
// 
//}

// Given initial position x^alpha, alpha, beta, determine the initial photon direction
// Uses alpha, beta definitions from Cunningham & Bardeen (1973)
void GetMCDirection(Photon *pphot, Real alpha, Real beta) {
  
  MonteCarloBlock *pmcb = pphot->pmy_mcb;
  //MCCoord *pco = pmcb->pcoord;
  Real kcon[NCOORD], kcov[NCOORD]; // kcon = k^alpha; kcov = k_alpha
  // kcov is set by alpha, beta, but kcon is what is integrated 
  Real gcon[NCOORD][NCOORD], gcov[NCOORD][NCOORD];
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  Real ucon[NCOORD], bcon[NCOORD];

  // calculate g_{alpha,beta}, g^{a,b} 
  pmcb->pcoord->Metric(pphot->x, gcov);
  pmcb->pcoord->InverseMetric(pphot->x, gcon);


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

}
  
  
// Given initial position x^alpha, alpha, beta, determine the initial photon direction
// Equivalent to Cunningham & Bardeen implementation up to minus signs in alpha, beta
void GetDirectionTetrad(Photon *pphot, Real alpha, Real beta) {
  
  // Set metric components
  MCCoord *pcoord = pphot->pmy_mcb->pcoord;
  Real gcov[NCOORD][NCOORD], gcon[NCOORD][NCOORD];
  pcoord->Metric(pphot->x, gcov);
  pcoord->InverseMetric(pphot->x,gcon);

  // Set tetrad vector for camera
  Real ucon[NCOORD] = {0.,0.,0.,1.}; // Static observer
  Real wcon[NCOORD] = {0,-1.,0.,0.}; // Q=1 points along projected BH symmetry axis 
  Real vcov[NCOORD] = {1.,0.,0.,1.};// Make image center point away from origin
  Real vcon[NCOORD];       
  CovToCon(vcov,vcon,gcon);

  // Construct tetrad
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  ConstructTetrad(ucon, vcon, wcon, gcov, econ, ecov);

  // Construct k in the tetrad
  Real ktet[NCOORD];
  Real kx = alpha / pphot->x[IMC1];
  Real ky = beta / pphot->x[IMC1];
  Real knorm = sqrt(1.+SQR(kx)+SQR(ky));
  ktet[IMC0] = -1.; // Photon is moving backward in time
  ktet[IMC1] = kx / knorm;
  ktet[IMC2] = ky / knorm;
  ktet[IMC3] = -1. / knorm; // points along radial direction
    
  TetradToCoordinate(ktet,pphot->k,econ);

}

void MidplaneCrossing(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover) {

 // check if photon has crossed midplane and whether to terminate or keep integrating

  //if (iphot == 2)
  //  pphot->PrintPhoton();

  if (pphot->x[IMC1] < r_outer + 1.0e-5) {
    pphot->status = DESTROYED;
    pphot->user_var[2] = 0.;
    pphot->user_var[3] = r_outer;
    pphot->stokes[1] = 0.;
    pphot->stokes[2] = 0.;
    return;
  }

  bool reverse = false;

  Real rdisk = 1.0e10; // outer radius of accretion disk for plotting
  if (plane_cross == 0) { // photon has not yet crossed the plane
      if (pphot->x[IMC2] >= (M_PI / 2.0)) { // photon has crossed plane for the first time
	if (pphot->x[IMC1] <= rdisk) { // photon is "close" to BH
	  Real step = -(pphot->x[IMC2] - M_PI/2.0) / pphot->k[IMC2];
	  for (int i = 0; i < NCOORD; i++) 
	    pphot->x[i] += pphot->k[i] * step;
	  /*photon_step(pphot, file_output);
	  for (int i=0; i < NCOORD; i++) 
	  pphot->k[i] *= -1;*/
          if (backward_integration)
            reverse = true;
          else
            pphot->status = ESCAPED;
	  plane_cross++;
	} else { // photon far away -> keep integrating
	  plane_cross++;
	}
      } 
    } else if ((plane_cross % 2) == 1) { // photon has crossed plane an odd number
      if (pphot->x[IMC2] <= (M_PI / 2.0)) { 
	if (pphot->x[IMC1] <= rdisk) {
	  pphot->status = ESCAPED;
	} else {
	  plane_cross++;
	}
      }
    } else if ((plane_cross % 2) == 0) { // photon has crossed plane an even number
      if (pphot->x[IMC2] >= (M_PI / 2.0)) {
	if (pphot->x[IMC1] <= rdisk) {
	  pphot->status = ESCAPED;
	} else {
	  plane_cross++;
	}
      }
    }

  if (reverse) {

    // Transform from tetrad frame to comoving frame, initialize stokes, and record
    // emission energy of photon in comoving frame

    Real gcov[NCOORD][NCOORD];
    pmcb->pcoord->Metric(pphot->x, gcov);

    Real r = pphot->x[IMC1];
    Real a = pmcb->pcoord->GetSpin();
    Real mbh = pmcb->pcoord->GetMass();
    Real omega = pow(mbh,0.5)/(pow(r, 3./2.) + a*pow(mbh,3./2.)); // circular velocity 
    // Initialize ucon and vcon (= z unit vector in symmetry plane)
    Real ucon[NCOORD];
    ucon[IMC0] = sqrt(-1.0/(gcov[IMC0][IMC0] + 2.*gcov[IMC0][IMC3]*omega +
                            SQR(omega)*gcov[IMC3][IMC3]));
    ucon[IMC1] = 0.;
    ucon[IMC2] = 0.;
    ucon[IMC3] = (ucon[IMC0])*omega;
    Real vcon[NCOORD];
    vcon[IMC0] = 0.;
    vcon[IMC1] = 0.;
    vcon[IMC2] = -1.;
    vcon[IMC3] = 0.;

    // Reverse photon direction
    for (int i = 0; i < NCOORD; i++) {
      pphot->k[i] *= -1;
      pphot->dk[i] *= -1;
    }

    // create tetrad basis
    Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
    ConstructTetrad(ucon, pphot->k, vcon, gcov, econ, ecov);

    // Initialize and transform Stokes vector
    pphot->stokes[0] = 1.0;
    pphot->stokes[1] = 1.0;
    pphot->stokes[2] = 0.0;
    pphot->stokes[3] = 0.0;
    std::complex<Real> tcopy[NCOORD][NCOORD];
    StokesToTensor(pphot->stokes,tcopy);
    ComplexTetradToCoordinate(tcopy,pphot->polten,econ);

    // Get photon energy in rest frame
    Real kcopy[NCOORD];
    CoordinateToTetrad(pphot->k, kcopy, ecov);
    pphot->user_var[2] = kcopy[IMC0];
    // Get radius at disk crossing
    pphot->user_var[3] = pphot->x[IMC1];
    // set plane crossing to zero
    plane_cross = 0;

#ifdef DEBUG
    printf("init: %d\n",iphot-1);
    Real cphi, sphi;
    Real sth, cth;
    cth = cos(pphot->x[IMC2]);
    sth = sin(pphot->x[IMC2]);
    cphi = cos(pphot->x[IMC3]);
    sphi = sin(pphot->x[IMC3]);
    czet = pphot->k[IMC1]*cth-sth*r*pphot->k[IMC2];
    szet = sqrt(1.-SQR(czet));
    spsi = (sphi*(pphot->k[IMC1]*sth+pphot->k[IMC2]*r*cth)+r*pphot->k[IMC3]*sth*cphi)/szet;
    cpsi = (cphi*(pphot->k[IMC1]*sth+pphot->k[IMC2]*r*cth)-r*pphot->k[IMC3]*sth*sphi)/szet;
    //cpsi = (pphot->k[IMC1]*cphi-r*pphot->k[IMC3]*sin(pphot->x[IMC3]))/szet;
    printf("psi: %e %e\n",spsi,cpsi);
    printf("zeta: %e %e\n",szet,czet);
    printf("th: %e %e\n",cth,sth);
    //printf("k: %e %e %e\n",szet*(cpsi*cphi+spsi*sphi),-czet/r,szet/r*(spsi*cphi-cpsi*sphi));
    //printf("e1: %e %e %e\n",-spsi,cpsi,0.);
    //printf("e2: %e %e %e\n",-czet*cpsi,czet*spsi,szet);
    //printf("e3: %e %e %e\n",szet*cpsi,szet*spsi,czet);
    Real kr = szet*sth*(cpsi*cphi+spsi*sphi)+czet*cth;
    Real kth = szet*cth/r*(cpsi*cphi+spsi*sphi)-czet*sth/r;
    Real kph = szet/(r*sth)*(spsi*cphi-cpsi*sphi);
    printf("k: %e %e %e\n",kr,kth,kph);
    Real Nrr = 2.*SQR(sth)*SQR(spsi*cphi-cpsi*sphi);
    Real Nrth = 2*sth*cth/r*SQR(spsi*cphi-cpsi*sphi);
    Real Nrph = 2./r*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
    Real Nthth = 2.*SQR(cth/r)*SQR(spsi*cphi-cpsi*sphi);
    Real Nthph = -2.*cth/SQR(r)/sth*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
    Real Nphph = 2./SQR(r*sth)*SQR(sphi*spsi+cphi*cpsi);
    printf("Nr: %e %e %e\n",Nrr,Nrth,Nrph);
    printf("Nth: %e %e\n",Nthth,Nthph);
    printf("Nph: %e\n",Nphph);
    printf("prod: %e %e %e %e %e\n",spsi*cphi-cpsi*sphi,cphi,sphi,cpsi,spsi);
    printf("x: %g %g %g %g %g %g\n",pphot->x[IMC0],pphot->x[IMC1],pphot->x[IMC2],pphot->x[IMC3],pphot->x[IMC1]*cos(pphot->x[IMC3]),pphot->x[IMC1]*sin(pphot->x[IMC3]));
    printf("k: %g %g %g %g\n",pphot->k[IMC0],pphot->k[IMC1],pphot->k[IMC2],pphot->k[IMC3]);
    printf("ttet[IMC0]: %e %e %e %e\n", tcopy[IMC0][IMC0].real(), tcopy[IMC0][IMC1].real(), tcopy[IMC0][IMC2].real(), tcopy[IMC0][IMC3].real());
    printf("ttet[IMC1]: %e %e %e %e\n", tcopy[IMC1][IMC0].real(), tcopy[IMC1][IMC1].real(), tcopy[IMC1][IMC2].real(), tcopy[IMC1][IMC3].real());
    printf("ttet[IMC2]: %e %e %e %e\n", tcopy[IMC2][IMC0].real(), tcopy[IMC2][IMC1].real(), tcopy[IMC2][IMC2].real(), tcopy[IMC2][IMC3].real());
    printf("ttet[IMC3]: %e %e %e %e\n", tcopy[IMC3][IMC0].real(), tcopy[IMC3][IMC1].real(), tcopy[IMC3][IMC2].real(), tcopy[IMC3][IMC3].real());
    printf("ttet[11]: %e\n",tcopy[IMC1][IMC1].real());
    printf("tcord[IMC0]: %e %e %e %e\n", pphot->polten[IMC0][IMC0].real(), pphot->polten[IMC0][IMC1].real(), pphot->polten[IMC0][IMC2].real(), pphot->polten[IMC0][IMC3].real());
    printf("tcord[IMC1]: %e %e %e %e\n", pphot->polten[IMC1][IMC0].real(), pphot->polten[IMC1][IMC1].real(), pphot->polten[IMC1][IMC2].real(), pphot->polten[IMC1][IMC3].real());
    printf("tcord[IMC2]: %e %e %e %e\n", pphot->polten[IMC2][IMC0].real(), pphot->polten[IMC2][IMC1].real(), pphot->polten[IMC2][IMC2].real(), pphot->polten[IMC2][IMC3].real());
    printf("tcord[IMC3]: %e %e %e %e\n", pphot->polten[IMC3][IMC0].real(), pphot->polten[IMC3][IMC1].real(), pphot->polten[IMC3][IMC2].real(), pphot->polten[IMC3][IMC3].real());
    //printf("pphot: %e %e %e %e %e\n",pphot->polten[IMC1][IMC1].real(),pphot->polten[IMC1][IMC2].real(),pphot->polten[IMC2][IMC1].real(),
    //printf("gcov[IMC0]: %e %e %e %e\n", gcov[IMC0][IMC0], gcov[IMC0][IMC1], gcov[IMC0][IMC2], gcov[IMC0][IMC3]);
    //printf("gcov[IMC1]: %e %e %e %e\n", gcov[IMC1][IMC0], gcov[IMC1][IMC1], gcov[IMC1][IMC2], gcov[IMC1][IMC3]);
    //printf("gcov[IMC2]: %e %e %e %e\n", gcov[IMC2][IMC0], gcov[IMC2][IMC1], gcov[IMC2][IMC2], gcov[IMC2][IMC3]);
    //printf("gcov[IMC3]: %e %e %e %e\n", gcov[IMC3][IMC0], gcov[IMC3][IMC1], gcov[IMC3][IMC2], gcov[IMC3][IMC3]);
    //Real gcon[4][4];
    //pcoord->InverseMetric(pphot->x,gcon);
    //printf("gcon[IMC0]: %e %e %e %e\n", gcon[IMC0][IMC0], gcon[IMC0][IMC1], gcon[IMC0][IMC2], gcon[IMC0][IMC3]);
    //printf("gcon[IMC1]: %e %e %e %e\n", gcon[IMC1][IMC0], gcon[IMC1][IMC1], gcon[IMC1][IMC2], gcon[IMC1][IMC3]);
    //printf("gcon[IMC2]: %e %e %e %e\n", gcon[IMC2][IMC0], gcon[IMC2][IMC1], gcon[IMC2][IMC2], gcon[IMC2][IMC3]);
    //printf("gcon[IMC3]: %e %e %e %e\n", gcon[IMC3][IMC0], gcon[IMC3][IMC1], gcon[IMC3][IMC2], gcon[IMC3][IMC3]);
    printf("econ[IMC0]: %e %e %e %e\n", econ[IMC0][IMC0], econ[IMC0][IMC1], econ[IMC0][IMC2], econ[IMC0][IMC3]);
    printf("econ[IMC1]: %e %e %e %e\n", econ[IMC1][IMC0], econ[IMC1][IMC1], econ[IMC1][IMC2], econ[IMC1][IMC3]);
    printf("econ[IMC2]: %e %e %e %e\n", econ[IMC2][IMC0], econ[IMC2][IMC1], econ[IMC2][IMC2], econ[IMC2][IMC3]);
    printf("econ[IMC3]: %e %e %e %e\n", econ[IMC3][IMC0], econ[IMC3][IMC1], econ[IMC3][IMC2], econ[IMC3][IMC3]);
#endif
  }
}
