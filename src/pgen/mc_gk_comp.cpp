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
static Real r0,th0,phi0;
static Real gcov0[4][4];
static Real muk,phik;
static Real rprev;
static bool first;

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// User function definitions
void TurningPointCheck(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover);


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
    printf("iphot: %d %d %d\n",rank,iphot,myn);
  }
#else
  iphot = 0;
#endif

  MCCoord *pcobl = new MCBoyerLindquist(1,1,1,false);
  Real x[4];
  x[IMC0] = 1.; x[IMC1] = r0; x[IMC2] = th0; x[IMC3] = phi0;
  pcobl->SetSpin(pin->GetReal("coord","a"));
  pcobl->SetMass(pin->GetReal("coord","m"));
  pcobl->Metric(x,gcov0);

  EnrollUserWorkInMove(TurningPointCheck);
}

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){

  nuser_var = 6;
  Real signa;
  if(pin->GetOrAddBoolean("problem", "corotating",true))
    signa = 1.;
  else
    signa = -1.;

  Real a = pin->GetReal("coord", "a");
  Real a2 = a*a;
  Real z1 = 1.0 + pow(1.0 - a2, 1./3.) * (pow(1. + a, 1./3.) + pow(1.0-a,1./3.));
  Real z2 = sqrt(3.*a2 + z1*z1);
  Real risco = 3.0 + z2 - signa * sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2)) + 1.0e-3;
  r0 = pin->GetReal("problem", "radius")*risco;
  th0 = pin->GetOrAddReal("problem", "theta",0.5)*M_PI;
  phi0 = pin->GetOrAddReal("problem", "phi", 0.)*M_PI;
  muk = pin->GetOrAddReal("problem", "muk",-0.8);
  phik = pin->GetOrAddReal("problem", "phik", 0.4);

}

void MonteCarloBlock::InitializePhoton(Photon *pphot) {

  MCCoord *pco = pphot->pmy_mcb->pcoord;

  // Set status flag
  pphot->status = EVOLVING;
  pphot->eweight = 1.;
  pphot->weight = 1.;
  
  // Emit photons from a large radius r >> 1 in units of [GM/c^2]. Ideally this would be 
  // at infinity (or 1/r = 0), but r ~ 1e3 should be fine as long as the region of 
  // interest is small compared to this initial distance. For -10 < alpha, beta < 10, the
  // small angle approximation should still be fine. 
  // To avoid issues with the initial position being on a boundary, I push the position
  // by a small epislon from the ideal starting position.

  pphot->x[IMC0] = 1.0;
  pphot->x[IMC1] = r0;
  pphot->x[IMC2] = th0;
  pphot->x[IMC3] = phi0;

  rprev = r0;
  // update the photon's zone indices
  pphot->i1 = -1;
  for(int i=pphot->pmy_mcb->is; i<=pphot->pmy_mcb->ie; i++) {
    if ((pphot->x[IMC1] >= pco->x1f(i)) && (pphot->x[IMC1] <= pco->x1f(i+1)))
      pphot->i1 = i;
  }
  if (pphot->i1 < 0) pphot->weight = -1.0;
 
  pphot->i2 = -1;
  for(int i=pphot->pmy_mcb->js; i<=pphot->pmy_mcb->je; i++) {
    if ((pphot->x[IMC2] >= pco->x2f(i)) && (pphot->x[IMC2] <= pco->x2f(i+1)))
      pphot->i2 = i;
  }
  if (pphot->i2 < 0) pphot->weight = -1.0;
 
  pphot->i3 = -1;
  for(int i=pphot->pmy_mcb->ks; i<=pphot->pmy_mcb->ke; i++) {
    if ((pphot->x[IMC3] >= pco->x3f(i)) && (pphot->x[IMC3] <= pco->x3f(i+1)))
      pphot->i3 = i;
  }
  if (pphot->i3 < 0) pphot->weight = -1.0;

  if (pphot->weight < 0) {
    printf("Warning: photon initial position not found on grid.\n");
    pphot->status = DESTROYED;
    pphot->PrintPhoton();
  }
  // Set the initial photon direction assuming "isotropic" emission
  pphot->energy = 1.0;

  int ith = iphot / 4;
  int iph = iphot % 4;
  Real cth = muk + 0.2 * static_cast<Real>(ith);
  Real phi = (phik + 0.4 * static_cast<Real>(iph)) * M_PI;
  Real sth = sqrt(1.-cth*cth);
 
  iphot++;

  pphot->k[IMC0] = pphot->energy;
  pphot->k[IMC1] = pphot->energy*sth*sin(phi);
  pphot->k[IMC2] = pphot->energy*cth;
  pphot->k[IMC3] = pphot->energy*sth*cos(phi);

  // Initialize Stokes vector as unpolarized
  pphot->stokes[0] = 1.0;
  pphot->stokes[1] = 0.0;
  pphot->stokes[2] = 0.0;
  pphot->stokes[3] = 0.0;

  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = 0.;
  pphot->sct_coef = 0.;

  // Transform to coordinate frame

  Real ucon[NCOORD];
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  Real kcopy[NCOORD];
  Real gcov[NCOORD][NCOORD];
  pco->Metric(pphot->x, gcov);

  Real r = pphot->x[IMC1];
  Real a = pco->GetSpin();
  Real omega = 1.0/(pow(r, 3./2.) + a); // circular velocity 
  // SWD: Eric used for both coordinates -- probably only correct for BL?
  ucon[IMC0] = sqrt(-1.0/(gcov[IMC0][IMC0] + 2.*gcov[IMC0][IMC3]*omega +
                          SQR(omega)*gcov[IMC3][IMC3]));
  ucon[IMC1] = 0.;
  ucon[IMC2] = 0.;
  ucon[IMC3] = (ucon[IMC0])*omega;
    
  // create tetrad basis
  ConstructTetrad(ucon, gcov, econ, ecov);

  // Transform to tetrad frame
  for (int i = 0; i < NCOORD; i++)
    kcopy[i] = pphot->k[i];

  // Transform k
  TetradToCoordinate(kcopy, pphot->k, econ);

  //  Initialize dK
  Real gamma[NCOORD][NCOORD][NCOORD];
  pco->Connect(pphot->x, gamma);

  for (int i = 0; i < 4; i++) {

    pphot->dk[i] = 
      -2.*(pphot->k[0]*(gamma[i][IMC0][IMC1]*pphot->k[IMC1]+gamma[i][IMC0][IMC2]*pphot->k[IMC2]+
                        gamma[i][IMC0][IMC3]*pphot->k[IMC3])+
           pphot->k[IMC1]*(gamma[i][IMC1][IMC2]*pphot->k[IMC2]+gamma[i][IMC1][IMC3]*pphot->k[IMC3])+
           ppot->k[IMC2]*gamma[i][IMC2][IMC3]*pphot->k[IMC3])-
      (gamma[i][IMC0][IMC0]*pphot->k[IMC0]*pphot->k[IMC0]+gamma[i][IMC1][IMC1]*pphot->k[IMC1]*pphot->k[IMC1]+
       gamma[i][IMC2][IMC2]*pphot->k[IMC2]*pphot->k[IMC2]+gamma[i][IMC3][IMC3]*pphot->k[IMC3]*pphot->k[IMC3]);
  }

  // Compute input variables for geokerr and store as user varibles for photon list

  // Geokerr uses BL coordinates so we first transfer from KS to BL and then compute
  // k_\alpha needed to define alpha, beta for geokerr

  Real alpha,beta;
  if (!pphot->pmy_mcb->boyerlindquist_flag) {

    Real delta = SQR(pphot->x[IMC1]) - 2 * pphot->x[IMC1] + SQR(a);
    Real kt0_bl = (pphot->k[IMC0] - 2.*pphot->x[IMC1]/delta*pphot->k[IMC1])*gcov0[IMC0][IMC0] 
      + pphot->k[IMC1]*gcov0[IMC0][IMC1] + pphot->k[IMC2]*gcov0[IMC0][IMC2] 
      + (pphot->k[IMC3] - a/delta*pphot->k[IMC1])*gcov0[IMC0][IMC3];
    Real kth0_bl = (pphot->k[IMC0] - 2.*pphot->x[IMC1]/delta*pphot->k[IMC1])*gcov0[IMC2][IMC0]
      + pphot->k[IMC1]*gcov0[IMC2][IMC1] + pphot->k[IMC2]*gcov0[IMC2][IMC2] 
      + (pphot->k[IMC3] - a/delta*pphot->k[IMC1])*gcov0[IMC2][IMC3];
    Real kphi0_bl = (pphot->k[IMC0] - 2.*pphot->x[IMC1]/delta*pphot->k[IMC1])*gcov0[IMC3][IMC0]
      + pphot->k[IMC1]*gcov0[IMC3][IMC1] + pphot->k[IMC2]*gcov0[IMC3][IMC2] 
      + (pphot->k[IMC3] - a/delta*pphot->k[IMC1])*gcov0[IMC3][IMC3];
    alpha = -kphi0_bl / kt0_bl;
    // Assumes with start at theta = pi/2, then beta^2 = q^2
    beta = sqrt(SQR(kth0_bl)/SQR(kt0_bl));
    if (kth0_bl > 0)
      beta = -beta;
  } else {
    Real kt0 = pphot->k[IMC0]*gcov0[IMC0][IMC0] + pphot->k[IMC1]*gcov0[IMC0][IMC1] + 
      pphot->k[IMC2]*gcov0[IMC0][IMC2] + pphot->k[IMC3]*gcov0[IMC0][IMC3];
    Real kth0 = pphot->k[IMC0]*gcov0[IMC2][IMC0] + pphot->k[IMC1]*gcov0[IMC2][IMC1] + 
      pphot->k[IMC2]*gcov0[IMC2][IMC2] + pphot->k[IMC3]*gcov0[IMC2][IMC3];
    Real kphi0 = pphot->k[IMC0]*gcov0[IMC3][IMC0] + pphot->k[IMC1]*gcov0[IMC3][IMC1] + 
      pphot->k[IMC2]*gcov0[IMC3][IMC2] + pphot->k[IMC3]*gcov0[IMC3][IMC3];
    alpha = -kphi0 / kt0;
    // SWD change?
    // Assumes with start at theta = pi/2, then beta^2 = q^2
    beta = sqrt(SQR(kth0)/SQR(kt0));
  }

  // Geokerr initialization paramters
  pphot->user_var[0] = alpha;
  pphot->user_var[1] = beta;
  pphot->user_var[2] = 1./pphot->x[IMC1]; // ui for 
  first = true; // for setting sign of du/dlambda
  pphot->user_var[4] = 0.; // tpr
  pphot->user_var[5] = cos(th0);


}



void MonteCarloBlock::FinalizePhoton(Photon *pphot) {
  
  Real a = pcoord->GetSpin();
  Real r_outer = 1.0 + sqrt(1.0 - SQR(a)) + 1.0e-3;

  // r is outside ISCO, transform into comoving frame tetrad, assuming
  // circular flow velocity
  if (pphot->x[IMC1] < r_outer + 1.0e-5) {
    pphot->status = ESCAPED;
    return;
  }
  
}

void TurningPointCheck(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover) {

  // Check if r is increasing and set sign of du/dlamda accordingly
  if (first) {
    if (pphot->x[IMC1] > rprev)
      pphot->user_var[3] = -1.;
    else
      pphot->user_var[3] = 1.;
    first = false;
  }
  
  if (pphot->user_var[3] > 0.) 
    if (pphot->x[IMC1] > rprev)
      pphot->user_var[4] = 1.;
  else
    if (pphot->x[IMC1] < rprev)
      pphot->user_var[4] = 1.;
  rprev = pphot->x[IMC1];

}
