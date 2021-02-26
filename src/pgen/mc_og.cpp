//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mc_og.cpp
//  \brief Problem generator for creating an observer grid near infinity 
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
static Real mu; 
static int plane_cross;
static bool backward_integration;
static FILE *input;

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// User function definitions
void MidplaneCrossing(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover);
void GetMCDirection(Photon *pphot, Real alpha, Real beta);

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
    if (backward_integration) {
      // Read in only parts covered by current rank
      Real dum;
      for(int i=0; i<iphot; i++) {
        fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &dum, &dum, &dum, 
               &dum, &dum, &dum, &dum, &dum, &dum);
      //printf("dum: %f\n",dum);
      }
    }
  }
#else
  iphot = 0;
#endif

  if (!backward_integration)
    EnrollUserWorkInMove(MidplaneCrossing);

}

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){


  nuser_var = 2;

  nrays = pin->GetInteger("montecarlo", "nphot");
  nalpha = nbeta = static_cast<int>(sqrt(static_cast<Real>(nrays)));
  backward_integration = pin->GetOrAddBoolean("problem","backward",false);

  if (backward_integration) {
    // Read input.dat
    input = fopen("input.dat", "r"); //should include starting polarization vectors

  } else {
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
    mu = pin->GetOrAddReal("problem", "mu", 0.5); // 0 < mu < 1
  }
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
  if (backward_integration) {
    fscanf(input, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &r, &theta, &phi, &kt, &kr, &kth, &kphi, &alpha0, &beta0);
  } else {
    r = 1.0e4 - 1.0e-3;
    theta = acos(mu);
    phi = M_PI / 2.0 + 1.0e-3;
    // Set alpha, beta by either iterating over photons
    ia = iphot / nalpha;
    ib = iphot % nalpha;
    alpha0 = alpha[ia];
    beta0 = beta[ib];
  }
  pphot->x[IMC0] = 1.0;
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

  // cweight is a constant weighting factor which accounts for the
  // emissivity of the grid zone in which the photon was emitted
  if (zone_weight_flag) {
    pphot->eweight = 1.;
    pphot->weight = 1.;
  }

  pphot->user_var[0] = alpha0;
  pphot->user_var[1] = beta0;

  printf("nrays: %d  iphot %d  nalpha: %d  ialpha: %d ",nrays, iphot, nalpha, ia);
  printf("ibeta: %d  alpha: %g beta: %g\n",ib, alpha0, beta0);

  iphot++;

  // Set the initial photon direction using alpha, beta and the position
  if (backward_integration) {
    pphot->k[IMC0] = kt;
    pphot->k[IMC1] = kr;
    pphot->k[IMC2] = kth;
    pphot->k[IMC3] = kphi;
    pphot->stokes[0] = 1.0;
    pphot->stokes[1] = 1.0;
    pphot->stokes[2] = 0.0;
    pphot->stokes[3] = 0.0;
    if (pphot->IsNanPhoton()) {
      //printf("\n\n\n\n\n\n\n\n");
      pphot->status = ESCAPED;
    }
    // Transform from tetrad frame to comoving frame
    // Construct the orthonormal tetrad
    Real ucon[NCOORD], vcon[NCOORD];
    Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
    Real kcopy[NCOORD];
    Real gcov[NCOORD][NCOORD];
    Real energy_shift;
    Real kdotu = 0.;
  
    pcoord->Metric(pphot->x, gcov);

    Real r = pphot->x[IMC1];
    Real a = pcoord->GetSpin();
    Real omega = 1.0/(pow(r, 3./2.) + a); // circular velocity 
    // Initialize ucon and vcon (= z unit vector in symmetry plane)
    ucon[IMC0] = sqrt(-1.0/(gcov[IMC0][IMC0] + 2.*gcov[IMC0][IMC3]*omega +
                            SQR(omega)*gcov[IMC3][IMC3]));
    ucon[IMC1] = 0.;
    ucon[IMC2] = 0.;
    ucon[IMC3] = (ucon[IMC0])*omega;

    vcon[IMC0] = 0.;
    vcon[IMC1] = 0.;
    vcon[IMC2] = 1./pphot->x[IMC1];
    vcon[IMC3] = 0.;
    // create tetrad basis
    ConstructTetrad(ucon, vcon, gcov, econ, ecov);

    // Transform to tetrad frame
    for (int i = 0; i < NCOORD; i++) 
      kdotu += pphot->k[i] * ucon[i]; // pphot->k in coordinate frame
    energy_shift = - pphot->k[IMC0] / kdotu; 

    for (int i = 0; i < NCOORD; i++)
      kcopy[i] = pphot->k[i];

    // Transform k
    TetradToCoordinate(kcopy, pphot->k, econ);
    // Invert k to reverse integration
    pphot->k[IMC1] *= -1.;
    pphot->k[IMC2] *= -1.;
    pphot->k[IMC3] *= -1.;
    // construct new tetrad to define stokes parameters
    ConstructTetrad(ucon, pphot->k, gcov, econ, ecov);

    // Initialize and transform Stokes vector
    pphot->stokes[0] = 1.0;
    pphot->stokes[1] = 1.0;
    pphot->stokes[2] = 0.0;
    pphot->stokes[3] = 0.0;
    std::complex<Real> tcopy[NCOORD][NCOORD];
    StokesToTensor(pphot->stokes,tcopy);

    //printf("%g %g %g %g\n",tcopy[IMC1][IMC1],tcopy[IMC1][IMC2],tcopy[IMC2][IMC1],
    //       tcopy[IMC2][IMC2]);
    ComplexTetradToCoordinate(tcopy,pphot->polten,econ);
    //printf("%g %g %g %g\n",pphot->polten[IMC1][IMC1],pphot->polten[IMC1][IMC2],pphot->polten[IMC2][IMC1],
    //       pphot->polten[IMC2][IMC2]);
  } else {
    // Initialize Stokes vector as unpolarized
    pphot->stokes[0] = 1.0;
    pphot->stokes[1] = 0.0;
    pphot->stokes[2] = 0.0;
    pphot->stokes[3] = 0.0;
    GetMCDirection(pphot, alpha0, beta0);
  }
  pphot->energy = pphot->k[IMC0];
  pphot->weight = 1.0;
  pphot->eweight = 1.0;
  for (int i=0; i<4; i++)
    pphot->dk[i] = 0.;
  //pphot->PrintPhoton();

  // Set plane crossing flag to zero
  plane_cross = 0;

  if (pphot->weight < 0.0) pphot->status = DESTROYED;
  if (backward_integration) {
 
    Real a = pcoord->GetSpin();
    Real r_outer = 1.0 + sqrt(1.0 - SQR(a)) + 1.0e-3;
    if ((pphot->x[IMC1] < r_outer) || (pphot->IsNanPhoton())) {
      pphot->status = ESCAPED;
      pphot->stokes[0] = 0.;
      pphot->stokes[1] = 0.;
      pphot->stokes[2] = 0.;
      pphot->stokes[3] = 0.;
    }
  }
  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = 0.;
  pphot->sct_coef = 0.;
  
}

void MidplaneCrossing(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover) {

 // check if photon has crossed midplane and whether to terminate or keep integrating

  //if (iphot == 2)
  //  pphot->PrintPhoton();

  if (pphot->status == DESTROYED) {
     pphot->status = ESCAPED;
     return;
  }

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
}

void MonteCarloBlock::FinalizePhoton(Photon *pphot) {
  
  Real a = pcoord->GetSpin();
  Real r_outer = 1.0 + sqrt(1.0 - SQR(a)) + 1.0e-3;

  if (pphot->status == DESTROYED) {
    pphot->status = ESCAPED;
    return;
  }

  // r is outside ISCO, transform into comoving frame tetrad, assuming
  // circular flow velocity
  if (pphot->x[IMC1] < r_outer + 1.0e-5) {
    pphot->status = ESCAPED;
    return;
  }
  
  if (backward_integration) {
    // Construct the orthonormal tetrad
    Real ucon[NCOORD];
    Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
    Real gcov[NCOORD][NCOORD];
    ucon[IMC0] = 1.;
    ucon[IMC1] = 0.;
    ucon[IMC2] = 0.;
    ucon[IMC3] = 0.;
    
    // create tetrad basis
    pcoord->Metric(pphot->x, gcov);

    Real vcon[NCOORD];
    vcon[IMC0] = 0.;
    vcon[IMC1] = 0.;
    vcon[IMC2] = 1./pphot->x[IMC1];
    vcon[IMC3] = 0.;

    ConstructTetrad(ucon, vcon, gcov, econ, ecov);
    std::complex<Real> tcopy[NCOORD][NCOORD];
    //printf("%g %g %g %g\n",pphot->polten[IMC1][IMC1],pphot->polten[IMC1][IMC2],pphot->polten[IMC2][IMC1],
    //       pphot->polten[IMC2][IMC2]);
    ComplexCoordinateToTetrad(pphot->polten,tcopy,ecov);
    TensorToStokes(tcopy,pphot->stokes);
    //printf("%g %g %g %g\n",tcopy[IMC1][IMC1],tcopy[IMC1][IMC2],tcopy[IMC2][IMC1],
    //       tcopy[IMC2][IMC2]);
  } else {
    // Construct the orthonormal tetrad
    Real ucon[NCOORD], vcon[NCOORD];
    Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
    Real kcopy[NCOORD];
    Real gcov[NCOORD][NCOORD];
    Real energy_shift;
    Real kdotu = 0.;

    pcoord->Metric(pphot->x, gcov);

    Real r = pphot->x[IMC1];
    Real omega = 1.0/(pow(r, 3./2.) + a); // circular velocity 
    
    ucon[IMC0] = sqrt(-1.0/(gcov[IMC0][IMC0] + 2.*gcov[IMC0][IMC3]*omega +
                            SQR(omega)*gcov[IMC3][IMC3]));
    ucon[IMC1] = 0.;
    ucon[IMC2] = 0.;
    ucon[IMC3] = (ucon[IMC0])*omega;
    vcon[IMC0] = 0.;
    vcon[IMC1] = 0.;
    vcon[IMC2] = 1./pphot->x[IMC1];
    vcon[IMC3] = 0.;
    // create tetrad basis
    ConstructTetrad(ucon, vcon, gcov, econ, ecov);

    // Transform to tetrad frame
    for (int i = 0; i < NCOORD; i++) 
      kdotu += pphot->k[i] * ucon[i]; // pphot->k in coordinate frame
    energy_shift = - pphot->k[IMC0] / kdotu; 

    for (int i = 0; i < NCOORD; i++)
      kcopy[i] = pphot->k[i];

    CoordinateToTetrad(kcopy, pphot->k, ecov);
    //printf("final: %g \n",pphot->x[0]);
  } 
}

// Given initial position x^alpha, alpha, beta, determine the initial photon direction
// Uses alpha, beta from Cunningham & Bardeen (1973)
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
  
  
