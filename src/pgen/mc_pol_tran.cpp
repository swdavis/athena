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


static Real r0,th0,ph0,rfin; 
static Real spsi,cpsi,szet,czet;
static Real Kp[2], polang;
static bool outsphere;

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

void Compute_f(Photon *pphot, Real gcov[4][4],Real K[2], Real a, Real f[4]);
void Compute_K(Photon *pphot, Real f[4], Real a, Real K[2]);
void Compute_delta_gamma(Photon *pphot, Real a, Real delta[3], Real gamma[3]);
void StokesTof(Photon *pphot, Real stokes[4], Real f[4]);
void fToStokes(Photon *pphot, Real f[4], Real stokes[4]);


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


}

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){

  // Determines whether spherical polar predictions are output
  outsphere = pin->GetOrAddBoolean("problem", "outsphere", false);
  //rfin = pin->GetOrAddReal("problem", "rfin", pin->GetReal("mesh","x1max"));
  polang = pin->GetOrAddReal("problem", "polang", 0.) * M_PI / 180.;
  r0 = pin->GetOrAddReal("problem", "r0", 10.);
  th0 = pin->GetOrAddReal("problem", "th0", 45.) * M_PI / 180.;
  ph0 = pin->GetOrAddReal("problem", "ph0", 0.) * M_PI / 180.;
  Real zeta = pin->GetReal("problem","zeta") * M_PI / 180.;
  Real psi = pin->GetReal("problem","psi") * M_PI / 180.;
  cpsi = cos(psi);
  spsi = sin(psi);
  czet = cos(zeta);
  szet = sin(zeta);
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

  pphot->x[IMC0] = 1.0;
  pphot->x[IMC1] = r0;
  pphot->x[IMC2] = th0;
  pphot->x[IMC3] = ph0;

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

  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = 0.;
  pphot->sct_coef = 0.;

  // Arrays for constructing the orthonormal tetrad
  Real ucon[NCOORD], vcon[NCOORD];
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  Real kcopy[NCOORD];
  Real gcov[NCOORD][NCOORD];
  pcoord->Metric(pphot->x, gcov);

 
  Real r = pphot->x[IMC1];
  Real cth = cos(pphot->x[IMC2]);
  Real sth = sin(pphot->x[IMC2]);
  Real cphi = cos(pphot->x[IMC3]);
  Real sphi = sin(pphot->x[IMC3]);
 
  if (outsphere) {
    // For testing, we can assume spherical polar relation.  This can
    // be obtained by setting m = 0 on <coords>.  Alternatively, we could
    // define in tetrad and transform for more general metric.  The
    // following assumes spherical polar:
    // kx = szet*cphi; ky = szet*sphi; kz = czet;
    Real kr = szet*sth*(cpsi*cphi+spsi*sphi)+czet*cth;
    Real kth = szet*cth/r*(cpsi*cphi+spsi*sphi)-czet*sth/r;
    Real kph = szet/(r*sth)*(spsi*cphi-cpsi*sphi);
    printf("k spherical: %e %e %e\n",kr,kth,kph);
  }

  // The following is code uses a tetrad to define k in a general coordinate
  // assuming angles psi and zeta are defined relative to tetrad.  Should reproduce
  // spherical prediction above but will apply even in boyer lindquist
  Real wcon[NCOORD];
  wcon[IMC0] = 0.;
  wcon[IMC1] = sth*sphi;
  wcon[IMC2] = cth*sphi/r;
  wcon[IMC3] = cphi/r/sth;
  ucon[IMC0] = 1.;
  ucon[IMC1] = 0.;
  ucon[IMC2] = 0.;
  ucon[IMC3] = 0.;
  vcon[IMC0] = 0.;
  vcon[IMC1] = cth;
  vcon[IMC2] = -sth;
  vcon[IMC3] = 0.;
  ConstructTetrad(ucon, vcon, wcon, gcov, econ, ecov);
  kcopy[IMC0] = 1.;
  kcopy[IMC1] = szet*cpsi;
  kcopy[IMC2] = szet*spsi;
  kcopy[IMC3] = czet;
  TetradToCoordinate(kcopy, pphot->k, econ);
 
  printf("k tetrad: %e %e %e\n",pphot->k[IMC1],pphot->k[IMC2],pphot->k[IMC3]);
  // Initialize dk to zero
  for (int i=0; i<4; i++)
    pphot->dk[i] = 0.;

  // Initialization for polarization
  // Initialize  Stokes vector
  pphot->stokes[0] = 1.0;
  pphot->stokes[1] = cos(2.*polang);
  pphot->stokes[2] = sin(2.*polang);
  pphot->stokes[3] = 0.0;

  // Construct polarization specific tetrad, which differs from the one 
  // above for specifying initial k direction 
  ucon[IMC0] = 1.;
  ucon[IMC1] = 0.;
  ucon[IMC2] = 0.;
  ucon[IMC3] = 0.;
  vcon[IMC0] = 0.;
  vcon[IMC1] = 0.;
  vcon[IMC2] = -1.;
  vcon[IMC3] = 0.;
 
  // construct tetrad and transform stokes vectors
  ConstructTetrad(ucon, pphot->k, vcon, gcov, econ, ecov);
  std::complex<Real> tcopy[NCOORD][NCOORD];
  StokesToTensor(pphot->stokes,tcopy);
  ComplexTetradToCoordinate(tcopy,pphot->polten,econ);
  
  // Specify polarization basis vectors for comparison
  Real fp[NCOORD];
  printf("pol angle: %e\n",polang);
  for (int i; i< NCOORD; ++i)  {
    // Use tetrad basis vectors
    fp[i] = econ[IMC1][i]*cos(polang)+econ[IMC2][i]*sin(polang);
  }
  Real fp0 = fp[IMC0];
  for (int i; i< NCOORD; ++i)  {
    // Use tetrad basis vectors
    fp[i] -= pphot->k[i]*fp0/pphot->k[IMC0];
  }
  printf("f: %e %e %e %e\n",fp[IMC0],fp[IMC1],fp[IMC2],fp[IMC3]);
  //generate_stokes_to_polarization(takes stokes and k, returns f)
  printf("fnorm: %e\n",DotVec(fp,fp,gcov));
  printf("f.k: %e\n",DotVec(fp,pphot->k,gcov));
  Real abh = pcoord->GetSpin();
  Compute_K(pphot, fp, abh, Kp);
  printf("K1, K2: %e %e\n", Kp[0], Kp[1]);
  Compute_f(pphot, gcov, Kp, abh, fp);
  printf("f2: %e %e %e %e\n",fp[IMC0],fp[IMC1],fp[IMC2],fp[IMC3]);
  printf("f2.k: %e\n",DotVec(fp,pphot->k,gcov));
  if (outsphere) {
    Real czett = pphot->k[IMC1]*cth-sth*r*pphot->k[IMC2];
    Real szett = sqrt(1.-SQR(czet));
    Real spsit = (sphi*(pphot->k[IMC1]*sth+pphot->k[IMC2]*r*cth)+r*pphot->k[IMC3]*sth*cphi)/szet;
    Real cpsit = (cphi*(pphot->k[IMC1]*sth+pphot->k[IMC2]*r*cth)-r*pphot->k[IMC3]*sth*sphi)/szet;

    printf("psi comp: %e %e %e %e\n",spsi,spsit,cpsi,cpsit);
    printf("zeta comp: %e %e %e %e\n",szet,szett,czet,czett);

    Real Nrr = 2.*SQR(sth)*SQR(spsi*cphi-cpsi*sphi);
    Real Nrth = 2*sth*cth/r*SQR(spsi*cphi-cpsi*sphi);
    Real Nrph = 2./r*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
    Real Nthth = 2.*SQR(cth/r)*SQR(spsi*cphi-cpsi*sphi);
    Real Nthph = -2.*cth/SQR(r)/sth*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
    Real Nphph = 2./SQR(r*sth)*SQR(sphi*spsi+cphi*cpsi);
    printf("Nr: %e %e %e\n",Nrr,Nrth,Nrph);
    printf("Nth: %e %e\n",Nthth,Nthph);
    printf("Nph: %e\n",Nphph);
    printf("tcord[IMC1]: %e %e %e\n", pphot->polten[IMC1][IMC1].real(), pphot->polten[IMC1][IMC2].real(), pphot->polten[IMC1][IMC3].real());
    printf("tcord[IMC2]: %e %e\n",pphot->polten[IMC2][IMC2].real(), pphot->polten[IMC2][IMC3].real());
    printf("tcord[IMC3]: %e\n",pphot->polten[IMC3][IMC3].real());
  }
  // output information about tetrad and polarization tensor in coordinate frame
  //printf("tcord[IMC0]: %e %e %e %e\n", pphot->polten[IMC0][IMC0].real(), pphot->polten[IMC0][IMC1].real(), pphot->polten[IMC0][IMC2].real(), pphot->polten[IMC0][IMC3].real());
  //printf("tcord[IMC1]: %e %e %e %e\n", pphot->polten[IMC1][IMC0].real(), pphot->polten[IMC1][IMC1].real(), pphot->polten[IMC1][IMC2].real(), pphot->polten[IMC1][IMC3].real());
  //printf("tcord[IMC2]: %e %e %e %e\n", pphot->polten[IMC2][IMC0].real(), pphot->polten[IMC2][IMC1].real(), pphot->polten[IMC2][IMC2].real(), pphot->polten[IMC2][IMC3].real());
  //printf("tcord[IMC3]: %e %e %e %e\n", pphot->polten[IMC3][IMC0].real(), pphot->polten[IMC3][IMC1].real(), pphot->polten[IMC3][IMC2].real(), pphot->polten[IMC3][IMC3].real());
  //printf("econ[IMC0]: %e %e %e %e\n", econ[IMC0][IMC0], econ[IMC0][IMC1], econ[IMC0][IMC2], econ[IMC0][IMC3]);
  //printf("econ[IMC1]: %e %e %e %e\n", econ[IMC1][IMC0], econ[IMC1][IMC1], econ[IMC1][IMC2], econ[IMC1][IMC3]);
  //printf("econ[IMC2]: %e %e %e %e\n", econ[IMC2][IMC0], econ[IMC2][IMC1], econ[IMC2][IMC2], econ[IMC2][IMC3]);
  //printf("econ[IMC3]: %e %e %e %e\n", econ[IMC3][IMC0], econ[IMC3][IMC1], econ[IMC3][IMC2], econ[IMC3][IMC3]);
    

  // Initialize the absorption and scattering extinction coefficients
  // to the values appropriate in the emitted zone
  pphot->abs_coef = 0.;
  pphot->sct_coef = 0.;

  if (pphot->IsNanPhoton()) {
    pphot->PrintPhoton();
    pphot->status = ESCAPED;
  }

}


void MonteCarloBlock::FinalizePhoton(Photon *pphot) {

  if (pphot->status == DESTROYED) {
    pphot->status = ESCAPED;
    pphot->PrintPhoton();
    return;
  }

  // Construct the orthonormal tetrad
  Real ucon[NCOORD];
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  Real gcov[NCOORD][NCOORD], gcon[NCOORD][NCOORD];
  ucon[IMC0] = 1.;
  ucon[IMC1] = 0.;
  ucon[IMC2] = 0.;
  ucon[IMC3] = 0.;
    
  // create tetrad basis
  pcoord->Metric(pphot->x, gcov);
  pcoord->InverseMetric(pphot->x,gcon);
  Real wcon[NCOORD] = {0,-1.,0.,0.}; // Q=1 points along projected BH symmetry axis 
  Real vcov[NCOORD] = {1.,0.,0.,1.};// Make image center point away from origin
  Real vcon[NCOORD]; 
  
  CovToCon(vcov,vcon,gcon);  
  ConstructTetrad(ucon, vcon, wcon, gcov, econ, ecov);
  std::complex<Real> tcopy[NCOORD][NCOORD];
  ComplexCoordinateToTetrad(pphot->polten,tcopy,ecov);
  TensorToStokes(tcopy,pphot->stokes);
    
  printf("econ[IMC0]: %e %e %e %e\n", econ[IMC0][IMC0], econ[IMC0][IMC1], econ[IMC0][IMC2], econ[IMC0][IMC3]);
  printf("econ[IMC1]: %e %e %e %e\n", econ[IMC1][IMC0], econ[IMC1][IMC1], econ[IMC1][IMC2], econ[IMC1][IMC3]);
  printf("econ[IMC2]: %e %e %e %e\n", econ[IMC2][IMC0], econ[IMC2][IMC1], econ[IMC2][IMC2], econ[IMC2][IMC3]);
  printf("econ[IMC3]: %e %e %e %e\n", econ[IMC3][IMC0], econ[IMC3][IMC1], econ[IMC3][IMC2], econ[IMC3][IMC3]);

  Real fp[NCOORD];
  Real abh = pcoord->GetSpin();

  // Compuate fp from Kp
  Compute_f(pphot, gcov, Kp, abh, fp);

  // Check that fp, k are othonormal
  Real fdotk = DotVec(fp,pphot->k,gcov);
  Real cosb1 = DotVec(fp,econ[IMC1],gcov);
  Real cosb2 = DotVec(fp,econ[IMC2],gcov);
  Real sinb1 = sqrt(1.-SQR(cosb1));
  if (cosb2 < 0.)
    sinb1 *= -1.;
  Real sin2x = 2.*sinb1*cosb1;
  Real cos2x = 2.*SQR(cosb1)-1.;
  Real norm = DotVec(fp,fp,gcov);
  printf("f: %e %e %e %e\n",fp[IMC0],fp[IMC1],fp[IMC2],fp[IMC3]);
  printf("norm f: %e \n",norm);
  printf("k.f: %e\n",fdotk);
  printf("k.k: %e\n",DotVec(fp,pphot->k,gcov));
  printf("fp.e1, fp.e2: %e %e\n",cosb1,cosb2);
  //printf("cos2b,sin2b: %e %e\n",cos2xb,sin2xb);
  if (outsphere) {
    Real cphi = cos(pphot->x[IMC3]);
    Real sphi = sin(pphot->x[IMC3]);
    Real r = pphot->x[IMC1];
    Real cth = cos(pphot->x[IMC2]);
    Real sth = sin(pphot->x[IMC2]);

    Real Nrr = 2.*SQR(sth)*SQR(spsi*cphi-cpsi*sphi);
    Real Nrth = 2*sth*cth/r*SQR(spsi*cphi-cpsi*sphi);
    Real Nrph = 2./r*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
    Real Nthth = 2.*SQR(cth/r)*SQR(spsi*cphi-cpsi*sphi);
    Real Nthph = -2.*cth/SQR(r)/sth*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
    Real Nphph = 2./SQR(r*sth)*SQR(sphi*spsi+cphi*cpsi);

    Real kr = szet*sth*(cpsi*cphi+spsi*sphi)+czet*cth;
    Real kth = szet*cth/r*(cpsi*cphi+spsi*sphi)-czet*sth/r;
    Real kph = szet/(r*sth)*(spsi*cphi-cpsi*sphi);
    printf("k sph pred: %e %e %e\n",kr,kth,kph);
    printf("k actual: %e %e %e\n",pphot->k[IMC1],pphot->k[IMC2],pphot->k[IMC3]);
    //printf("x actual: %e %e %e %e\n",pphot->x[IMC0],pphot->x[IMC1],pphot->x[IMC2],pphot->x[IMC3]);
    printf("inclination: %g\n",pphot->x[IMC2]/M_PI*180.);
    printf("Nr: %e %e %e\n",Nrr,Nrth,Nrph);
    printf("Nth: %e %e\n",Nthth,Nthph);
    printf("Nph: %e\n",Nphph);
    printf("tcord[IMC1]: %e %e %e\n", pphot->polten[IMC1][IMC1].real(), pphot->polten[IMC1][IMC2].real(), pphot->polten[IMC1][IMC3].real());
    printf("tcord[IMC2]: %e %e\n",pphot->polten[IMC2][IMC2].real(), pphot->polten[IMC2][IMC3].real());
    printf("tcord[IMC3]: %e\n",pphot->polten[IMC3][IMC3].real());

    Real cosd = spsi*sphi+cpsi*cphi;
    Real sind = spsi*cphi-cpsi*sphi;
    Real I = SQR(cosd)+SQR(sind*cth);
    Real Q = SQR(cosd)-SQR(sind*cth);
    Real U = 2*cth*cosd*sind;
    printf("stokes sph pred: %e %e\n",Q/I,U/I);
  }

  //printf("tcord[IMC0]: %e %e %e %e\n", pphot->polten[IMC0][IMC0].real(), pphot->polten[IMC0][IMC1].real(), pphot->polten[IMC0][IMC2].real(), pphot->polten[IMC0][IMC3].real());
  //printf("tcord[IMC1]: %e %e %e %e\n", pphot->polten[IMC1][IMC0].real(), pphot->polten[IMC1][IMC1].real(), pphot->polten[IMC1][IMC2].real(), pphot->polten[IMC1][IMC3].real());
  //printf("tcord[IMC2]: %e %e %e %e\n", pphot->polten[IMC2][IMC0].real(), pphot->polten[IMC2][IMC1].real(), pphot->polten[IMC2][IMC2].real(), pphot->polten[IMC2][IMC3].real());
  //printf("tcord[IMC3]: %e %e %e %e\n", pphot->polten[IMC3][IMC0].real(), pphot->polten[IMC3][IMC1].real(), pphot->polten[IMC3][IMC2].real(), pphot->polten[IMC3][IMC3].real());
  //printf("ttet[IMC0]: %e %e %e %e\n", tcopy[IMC0][IMC0].real(), tcopy[IMC0][IMC1].real(), tcopy[IMC0][IMC2].real(), tcopy[IMC0][IMC3].real());
  //printf("ttet[IMC1]: %e %e %e %e\n", tcopy[IMC1][IMC0].real(), tcopy[IMC1][IMC1].real(), tcopy[IMC1][IMC2].real(), tcopy[IMC1][IMC3].real());
  //printf("ttet[IMC2]: %e %e %e %e\n", tcopy[IMC2][IMC0].real(), tcopy[IMC2][IMC1].real(), tcopy[IMC2][IMC2].real(), tcopy[IMC2][IMC3].real());
  //printf("ttet[IMC3]: %e %e %e %e\n", tcopy[IMC3][IMC0].real(), tcopy[IMC3][IMC1].real(), tcopy[IMC3][IMC2].real(), tcopy[IMC3][IMC3].real());


  printf("stokes actual: %e %e\n",pphot->stokes[1],pphot->stokes[2]);
  printf("stokes from fp: %e %e\n",cos2x,sin2x);
  //printf("ecov[IMC0]: %e %e %e %e\n", ecov[IMC0][IMC0], ecov[IMC0][IMC1], ecov[IMC0][IMC2], ecov[IMC0][IMC3]);
  //printf("ecov[IMC1]: %e %e %e %e\n", ecov[IMC1][IMC0], ecov[IMC1][IMC1], ecov[IMC1][IMC2], ecov[IMC1][IMC3]);
  //printf("ecov[IMC2]: %e %e %e %e\n", ecov[IMC2][IMC0], ecov[IMC2][IMC1], ecov[IMC2][IMC2], ecov[IMC2][IMC3]);
  //printf("ecov[IMC3]: %e %e %e %e\n", ecov[IMC3][IMC0], ecov[IMC3][IMC1], ecov[IMC3][IMC2], ecov[IMC3][IMC3]);
  
  
}


void compute_delta_gamma(Photon *pphot, Real a, Real delta[3], Real gamma[3]){
  //compute the r, theta, and phi parts of the delta and gamma constants for the K and f computations
  Real& r = pphot->x[IMC1]; 
  Real& theta = pphot->x[IMC2];
  Real& k_t = pphot->k[IMC0];
  Real& k_r = pphot->k[IMC1];
  Real& k_theta = pphot->k[IMC2];
  Real& k_phi = pphot->k[IMC3];

  delta[0] = r*k_t - r*a*SQR(sin(theta))*k_phi;
  delta[1] = a*a*sin(theta)*cos(theta)*k_t - a*cos(theta)*sin(theta)*(r*r+a*a)*k_phi;
  delta[2] = r*a*SQR(sin(theta))*k_r + a*cos(theta)*sin(theta)*(r*r+a*a)*k_theta;
  gamma[0] = a*cos(theta)*k_t - a*a*cos(theta)*SQR(sin(theta))*k_phi;
  gamma[1] = r*(r*r+a*a)*sin(theta)*k_phi - a*r*sin(theta)*k_t;
  gamma[2] = a*a*cos(theta)*sin(theta)*sin(theta)*k_r - r*(r*r+a*a)*sin(theta)*k_theta;

}

void Compute_K(Photon *pphot, Real f[4], Real a, Real K[2]){
  
  Real delta[3];
  Real gamma[3];
  compute_delta_gamma(pphot, a, delta, gamma);
  
  K[0] = delta[0]*f[IMC1] + delta[1]*f[IMC2] + delta[2]*f[IMC3];
  K[1] = gamma[0]*f[IMC1] + gamma[1]*f[IMC2] + gamma[2]*f[IMC3];

}

void Compute_f(Photon *pphot, Real gcov[4][4], Real K[2], Real a, Real f[4]){

  Real& r = pphot->x[IMC1]; 
  Real& theta = pphot->x[IMC2];
  Real& k_t = pphot->k[IMC0];
  Real& k_r = pphot->k[IMC1];
  Real& k_theta = pphot->k[IMC2];
  Real& k_phi = pphot->k[IMC3]; 

  Real delta[3];
  Real gamma[3];
  compute_delta_gamma(pphot, a, delta, gamma);

  Real g_phiphi = gcov[IMC3][IMC3];
  Real g_thetatheta = gcov[IMC2][IMC2];
  Real g_phit = gcov[IMC3][IMC0];
  Real g_rr = gcov[IMC1][IMC1];

  Real N = (gamma[1]*delta[0] - gamma[0]*delta[1])*g_phiphi*k_phi - (gamma[2]*delta[0] - gamma[0]*delta[2])*g_thetatheta*k_theta 
    + (gamma[1]*delta[0] - gamma[0]*delta[1])*g_phit*k_t + (gamma[2]*delta[1] - gamma[1]*delta[2])*g_rr*k_r;

  f[IMC0] = 0.;
  f[IMC1] = ((gamma[1]*K[0]-delta[1]*K[1])*(g_phiphi*k_phi+g_phit*k_t)-(gamma[2]*K[0]-delta[2]*K[1])*g_thetatheta*k_theta)/N;
  f[IMC2] = -((gamma[0]*K[0]-delta[0]*K[1])*(g_phiphi*k_phi + g_phit*k_t)-(gamma[2]*K[0]-delta[2]*K[1])*g_rr*k_r)/N;
  f[IMC3] = ((gamma[0]*K[0]-delta[0]*K[1])*g_thetatheta*k_theta-(gamma[1]*K[0] - delta[1]*K[1])*g_rr*k_r)/N;

  //Real normf = DotVec(f,f,gcov);
  //for (int i; i< NCOORD; ++i)  
  //  f[i] /= sqrt(normf);
}

Real dot(Real a[4], Real b[4]){
  //3D vector dot product
  return a[IMC1]*b[IMC1] + a[IMC2]*b[IMC2] + a[IMC3]*b[IMC3];
}

void cross(Real a[4], Real b[4], Real c[4]){
  //3D vector cross product (t set to 0)
  c[IMC0] = 0.;
  c[IMC1] = a[IMC2]*b[IMC3] - a[IMC3]*b[IMC2];
  c[IMC2] = a[IMC3]*b[IMC1] - a[IMC1]*b[IMC3];
  c[IMC3] = a[IMC1]*b[IMC2] - a[IMC2]*b[IMC1];
}

Real mag(Real a[4]){
  //magnitude of a vector
  return sqrt(dot(a, a));
}

void proj_ab(Real a[4], Real b[4], Real c[4]){
  //projection of a onto b (3D, t set to 0)
  Real q = dot(a, b)/dot(b, b);
  c[IMC0] = 0.;
  c[IMC1] = q*b[IMC1];
  c[IMC2] = q*b[IMC2];
  c[IMC3] = q*b[IMC3];
}

void rej_ab(Real a[4], Real b[4], Real c[4]){
  //rejection of a from b (3D, t set to 0)
  Real q = dot(a, b)/dot(b, b);
  c[IMC0] = 0.;
  c[IMC1] = a[IMC1] - q*b[IMC1];
  c[IMC2] = a[IMC2] - q*b[IMC2];
  c[IMC3] = a[IMC3] - q*b[IMC3];
}

void orthonormal_basis(Real k_vec[4], Real z[4], Real a[4], Real b[4]){
  rej_ab(z, k_vec, a);
  cross(k_vec, a, b);
  //printf("a %g %g %g %g, a dot k %g\n", a[IMC0], a[IMC1], a[IMC2], a[IMC3], dot(a, k_vec));
  //printf("b %g %g %g %g, b dot k %g\n", b[IMC0], b[IMC1], b[IMC2], b[IMC3], dot(b, k_vec));
  Real a_unit[4] = {a[IMC1]/sqrt(dot(a, a)), a[IMC2]/sqrt(dot(a, a)), a[IMC3]/sqrt(dot(a, a)), a[IMC0]/sqrt(dot(a, a))};
  //printf("a_unit %g %g %g %g, a dot k %g\n", a_unit[IMC0], a_unit[IMC1], a_unit[IMC2], a_unit[IMC3], dot(a_unit, k_vec));
  Real b_unit[4] = {b[IMC1]/sqrt(dot(b, b)), b[IMC2]/sqrt(dot(b, b)), b[IMC3]/sqrt(dot(b, b)), b[IMC0]/sqrt(dot(b, b))};
  //printf("b_unit %g %g %g %g, b dot k %g\n", b_unit[IMC0], b_unit[IMC1], b_unit[IMC2], b_unit[IMC3], dot(b_unit, k_vec));
  a = a_unit;
  b = b_unit;
}

void StokesTof(Photon *pphot, Real stokes[4], Real f[4]){
  //convert a stokes four vector to the polarization 4 vector
  printf("StokesTof:\n");
  Real k_vec[4] ={pphot->k[IMC1], pphot->k[IMC2], pphot->k[IMC3], pphot->k[IMC0]};
  printf("k: %g %g %g %g\n", k_vec[IMC0], k_vec[IMC1], k_vec[IMC2], k_vec[IMC3]);
  Real& Q = stokes[1];
  Real& U = stokes[2];
  Real z[4] = {0., 0., -1., 0.};
  Real a[4];
  Real b[4];
  orthonormal_basis(k_vec, z, a, b);
  printf("U Q U/Q: %g %g %g\n", U, Q, U/Q);
  Real psi = 0.5*atan(U/Q);
  printf("psi: %g\n", psi);
  printf("theta, phi: %g %g\n", pphot->x[IMC2], pphot->x[IMC3]);
  f[IMC0] = cos(psi)*a[IMC0] + sin(psi)*b[IMC0];
  f[IMC1] = cos(psi)*a[IMC1] + sin(psi)*b[IMC1];
  f[IMC2] = cos(psi)*a[IMC2] + sin(psi)*b[IMC2];
  f[IMC3] = cos(psi)*a[IMC3] + sin(psi)*b[IMC3];

}

void fToStokes(Photon *pphot, Real f[4], Real stokes[4]){
  //make so checks for divide by zero and the correct quadrant
  //replace k_vc with rhat
  printf("fToStokes:\n");
  Real k_vec[4] = {pphot->k[IMC1], pphot->k[IMC2], pphot->k[IMC3], pphot->k[IMC0]};
  //Real r[4] = {0., 1., 0., 0.};
  Real z[4] = {0., 0., -1., 0.};
  Real a[4];
  Real b[4];
  Real f_a[4];
  Real f_b[4];
  orthonormal_basis(k_vec, z, a, b);
  proj_ab(f, a, f_a);
  proj_ab(f, b, f_b);
  printf("mag fa, mag fb, magfb/magfa, magfa/magfb: %g %g %g %g\n", mag(f_a), mag(f_b), mag(f_b)/mag(f_a), mag(f_a)/mag(f_b));
  Real psi = 0.5*atan(mag(f_a)/mag(f_b));
  printf("psi (fa/fb): %g\n", psi);
  psi = 0.5*atan(mag(f_b)/mag(f_a));
  //Real psi = 0.;//initial photon angle//acos((f[IMC2]-(b[IMC2]*f[IMC1])/b[IMC1])/(a[IMC2]-(b[IMC2]*a[IMC1])/b[IMC1]));
  Real& theta = pphot->x[IMC2];
  Real& phi = pphot->x[IMC3];
  printf("psi (fb/fa): %g\n", psi);
  printf("theta, phi: %g %g\n", pphot->x[IMC2], pphot->x[IMC3]);
  //Real I = cos(psi-phi)*cos(psi-phi) + sin(psi-phi)*sin(psi-phi)*cos(theta)*cos(theta);
  //Real Q = cos(psi-phi)*cos(psi-phi) - sin(psi-phi)*sin(psi-phi)*cos(theta)*cos(theta);
  //Real U = cos(theta)*sin(2*(psi-phi));
  Real Q = cos(2*psi);
  Real U = sin(2*psi);
  Real I = sqrt(Q*Q + U*U);
  Real V = 0.;
  
  stokes[0] = I;
  stokes[1] = Q;
  stokes[2] = U;
  stokes[3] = V;
}
  
