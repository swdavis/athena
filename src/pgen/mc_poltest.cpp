//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mc_poltest.cpp
//  \brief Problem generator for testing polarized transport
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
#include "../monte_carlo/photonpusher.hpp"
#include "../globals.hpp"

#if !MONTE_CARLO_ENABLED
#error "This problem requires monte carlo"
#endif

namespace {

  // global variables
  Real r0,th0,ph0,rfin;
  Real spsi,cpsi,szet,czet;
  Real Kp[2], polang;
  bool outsphere;
  int i1start,i2start,i3start;

  // user function definitions
  void ComputeDeltaGamma(Real x[4], Real k[4], Real a, Real delta[3], Real gamma[3]);
  void ComputeK(Real x[4], Real k[4], Real f[4], Real a, Real K[2]);
  void Computef(Real x[4], Real k[4], Real gcov[4][4],Real K[2], Real a, Real f[4]);
}

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

//========================================================================================
//! \fn void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin)
//! \brief Analogous to problem generator but used in support of InitializePhoton
//========================================================================================

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {

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

  // set the photon samples 's initial zone indices
  MCCoord *pco = pcoord;
  i1start = -1;
  for(int i=is; i<=ie; i++) {
    if ((r0 >= pcoord->x1f(i)) && (r0 < pcoord->x1f(i+1)))
      i1start = i;
  }
  i2start = -1;
  for(int i=js; i<=je; i++) {
    if ((th0 >= pcoord->x2f(i)) && (th0 < pcoord->x2f(i+1)))
      i2start = i;
  }
  i3start = -1;
  for(int i=ks; i<=ke; i++) {
    if ((ph0 >= pcoord->x3f(i)) && (ph0 < pcoord->x3f(i+1)))
      i3start = i;
  }
  if ((i1start < 0) || (i2start < 0) || (i3start < 0)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in InitUserMonteCarloBlockData" << std::endl
        << "Initial position not found within domain." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}

//========================================================================================
//! \fn void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe, int etype)
//! \brief Initializes Photon packets before integration
//========================================================================================

void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe, int etype) {

  MCCoord *pco = pphot->pmy_mcb->pcoord;

  for (int ip=ips; ip<=ipe; ip++) {

    // Set status flag
    pphot->statp[ip] = EVOLVING;

    // Initialize photon position
    Real x[NCOORD];
    x[IMC0] = pphot->x0p[ip] = 0.;
    x[IMC1] = pphot->x1p[ip] = r0;
    x[IMC2] = pphot->x2p[ip] = th0;
    x[IMC3] = pphot->x3p[ip] = ph0;

    // initialize cell coordinates
    pphot->i1p[ip] = i1start;
    pphot->i2p[ip] = i2start;
    pphot->i3p[ip] = i3start;

    // set weight to 1
    pphot->wp[ip] = 1.;

    // Initialize the absorption and scattering extinction coefficients
    // to the values appropriate in the emitted zone
    pphot->acp[ip] = 0.;
    pphot->scp[ip] = 0.;

    Real r = x[IMC1];
    Real cth = cos(x[IMC2]);
    Real sth = sin(x[IMC2]);
    Real cphi = cos(x[IMC3]);
    Real sphi = sin(x[IMC3]);

    if (outsphere) {
      // For testing, we can assume spherical polar relation.  This can
      // be obtained by setting m = 0 on <coords>.  Alternatively, we could
      // define in tetrad and transform for more general metric.  The
      // following assumes spherical polar:
      // kx = szet*cphi; ky = szet*sphi; kz = czet;
      Real kr = szet*sth*(cpsi*cphi+spsi*sphi)+czet*cth;
      Real kth = szet*cth/r*(cpsi*cphi+spsi*sphi)-czet*sth/r;
      Real kph = szet/(r*sth)*(spsi*cphi-cpsi*sphi);
      printf("Assuming th0=90 deg.\n");
      printf("k spherical: %e %e %e\n",kr,kth,kph);
    }

    // The following code uses a tetrad to define k in a general coordinate
    // assuming angles psi and zeta are defined relative to tetrad.  Should reproduce
    // spherical prediction above but will apply even in boyer lindquist
    // Construct orthonormal tetrad
    Real gcov[4][4];
    pcoord->Metric(x, gcov);

    Real wcon[4];
    wcon[IMC0] = 0.;
    wcon[IMC1] = sth*sphi;
    wcon[IMC2] = cth*sphi/r;
    wcon[IMC3] = cphi/r/sth;
    Real ucon[4];
    ucon[IMC0] = 1.;
    ucon[IMC1] = 0.;
    ucon[IMC2] = 0.;
    ucon[IMC3] = 0.;
    Real vcon[4];
    vcon[IMC0] = 0.;
    vcon[IMC1] = cth;
    vcon[IMC2] = -sth;
    vcon[IMC3] = 0.;
    Real econ[4][4], ecov[4][4];
    ConstructTetrad(ucon, vcon, wcon, gcov, econ, ecov);

    Real kcopy[4];
    kcopy[IMC0] = 1.;
    kcopy[IMC1] = szet*cpsi;
    kcopy[IMC2] = szet*spsi;
    kcopy[IMC3] = czet;
    Real k[4];
    TetradToCoordinate(kcopy, k, econ);
    pphot->k0p[ip] = k[IMC0];
    pphot->k1p[ip] = k[IMC1];
    pphot->k2p[ip] = k[IMC2];
    pphot->k3p[ip] = k[IMC3];

    printf("k tetrad: %d %e %e %e\n",ip,k[IMC1],k[IMC2],k[IMC3]);
    // Initialize dk to zero
    pphot->dk0p[ip] = 0;
    pphot->dk1p[ip] = 0;
    pphot->dk2p[ip] = 0;
    pphot->dk3p[ip] = 0;

    // Initialization for polarization
    // Initialize  Stokes vector
    Real stokes[4];
    stokes[0] = 1.0;
    stokes[1] = cos(2.*polang);
    stokes[2] = sin(2.*polang);
    stokes[3] = 0.0;

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
    ConstructTetrad(ucon, k, vcon, gcov, econ, ecov);

    std::complex<Real> tcopy[4][4];
    StokesToTensor(stokes,tcopy);
    pphot->PolarizationToCoord(tcopy,econ,ip);

    // Specify polarization basis vectors for comparison
    Real fp[4];
    printf("pol angle: %d %e\n",ip,polang);
    for (int i; i< 4; ++i)  {
      // Use tetrad basis vectors
      fp[i] = econ[IMC1][i]*cos(polang)+econ[IMC2][i]*sin(polang);
    }
    Real fp0 = fp[IMC0];
    for (int i; i< 4; ++i)  {
      // Use tetrad basis vectors
      fp[i] -= k[i]*fp0 / k[IMC0];
    }
    printf("f: %e %e %e %e\n",fp[IMC0],fp[IMC1],fp[IMC2],fp[IMC3]);
    //generate_stokes_to_polarization (takes stokes and k, returns f)
    printf("fnorm: %e\n",DotVec(fp,fp,gcov));
    printf("f.k: %e\n",DotVec(fp,k,gcov));
    Real abh = pcoord->GetSpin();
    ComputeK(x, k, fp, abh, Kp);
    printf("K1, K2: %e %e\n", Kp[0], Kp[1]);
    Computef(x, k,  gcov, Kp, abh, fp);
    printf("f2: %e %e %e %e\n",fp[IMC0],fp[IMC1],fp[IMC2],fp[IMC3]);
    printf("f2.k: %e\n",DotVec(fp,k,gcov));
    if (outsphere) {
      Real czett = k[IMC1]*cth-sth*r*k[IMC2];
      Real szett = sqrt(1.-SQR(czet));
      Real spsit = (sphi*(k[IMC1]*sth+k[IMC2]*r*cth)+r*k[IMC3]*sth*cphi)/szet;
      Real cpsit = (cphi*(k[IMC1]*sth+k[IMC2]*r*cth)-r*k[IMC3]*sth*sphi)/szet;

      printf("psi comp: %e %e %e %e\n",spsi,spsit,cpsi,cpsit);
      printf("zeta comp: %e %e %e %e\n",szet,szett,czet,czett);

      Real Nrr = 2.*SQR(sth)*SQR(spsi*cphi-cpsi*sphi);
      Real Nrth = 2*sth*cth/r*SQR(spsi*cphi-cpsi*sphi);
      Real Nrph = 2./r*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
      Real Nthth = 2.*SQR(cth/r)*SQR(spsi*cphi-cpsi*sphi);
      Real Nthph = -2.*cth/SQR(r)/sth*((SQR(cpsi)-SQR(spsi))*sphi*cphi
                                       -(SQR(cphi)-SQR(sphi))*spsi*cpsi);
      Real Nphph = 2./SQR(r*sth)*SQR(sphi*spsi+cphi*cpsi);
      printf("Nr: %e %e %e\n",Nrr,Nrth,Nrph);
      printf("Nth: %e %e\n",Nthth,Nthph);
      printf("Nph: %e\n",Nphph);
      printf("tcord[IMC1]: %e %e %e\n", pphot->polten[IMC1*4+IMC1][ip].real(),
             pphot->polten[IMC1*4+IMC2][ip].real(),pphot->polten[IMC1*4+IMC3][ip].real());
      printf("tcord[IMC2]: %e %e\n",pphot->polten[IMC2*4+IMC2][ip].real(),
             pphot->polten[IMC2*4+IMC3][ip].real());
      printf("tcord[IMC3]: %e\n",pphot->polten[IMC3*4+IMC3][ip].real());
    }


    if (pphot->IsNanPhoton(ip)) {
      pphot->PrintPhoton(ip);
      pphot->statp[ip] = ESCAPED;
    }
  } // loop over ip
}

//========================================================================================
//! \fn void MonteCarloBlock::FinalizePhoton(Photon *pphot, int ip)
//! \brief Complete work at end of photon packets before integration
//========================================================================================

void MonteCarloBlock::FinalizePhoton(Photon *pphot, int ip) {

  if (pphot->statp[ip] == DESTROYED) {
    pphot->statp[ip] = ESCAPED;
    pphot->PrintPhoton(ip);
    return;
  }

  // Construct the orthonormal tetrad
  Real ucon[4];
  ucon[IMC0] = 1.;
  ucon[IMC1] = 0.;
  ucon[IMC2] = 0.;
  ucon[IMC3] = 0.;

  // create tetrad basis
  Real gcov[4][4], gcon[4][4];
  Real x[4];
  x[IMC0] = pphot->x0p[ip];
  x[IMC1] = pphot->x1p[ip];
  x[IMC2] = pphot->x2p[ip];
  x[IMC3] = pphot->x3p[ip];
  pcoord->Metric(x, gcov);
  pcoord->InverseMetric(x,gcon);
  Real wcon[4] = {0,-1.,0.,0.}; // Q=1 points along projected BH symmetry axis
  Real vcov[4] = {1.,0.,0.,1.};// Make image center point away from origin
  Real vcon[4];
  CovToCon(vcov,vcon,gcon);
  Real econ[4][4], ecov[4][4];
  ConstructTetrad(ucon, vcon, wcon, gcov, econ, ecov);
  std::complex<Real> tcopy[4][4];
  pphot->PolarizationToTetrad(tcopy,ecov,ip);
  Real stokes[4];
  TensorToStokes(tcopy,stokes);

  printf("econ[IMC0]: %e %e %e %e\n",econ[IMC0][IMC0],econ[IMC0][IMC1],econ[IMC0][IMC2],
         econ[IMC0][IMC3]);
  printf("econ[IMC1]: %e %e %e %e\n",econ[IMC1][IMC0],econ[IMC1][IMC1],econ[IMC1][IMC2],
         econ[IMC1][IMC3]);
  printf("econ[IMC2]: %e %e %e %e\n",econ[IMC2][IMC0],econ[IMC2][IMC1],econ[IMC2][IMC2],
         econ[IMC2][IMC3]);
  printf("econ[IMC3]: %e %e %e %e\n",econ[IMC3][IMC0],econ[IMC3][IMC1],econ[IMC3][IMC2],
         econ[IMC3][IMC3]);

  Real fp[4];
  Real abh = pcoord->GetSpin();

  Real k[4];
  k[IMC0] = pphot->k0p[ip];
  k[IMC1] = pphot->k1p[ip];
  k[IMC2] = pphot->k2p[ip];
  k[IMC3] = pphot->k3p[ip];

  // Compuate fp from Kp
  Computef(x, k, gcov, Kp, abh, fp);

  // Check that fp, k are othonormal
  Real fdotk = DotVec(fp, k, gcov);
  Real cosb1 = DotVec(fp, econ[IMC1], gcov);
  Real cosb2 = DotVec(fp, econ[IMC2], gcov);
  Real sinb1 = sqrt(1.-SQR(cosb1));
  if (cosb2 < 0.)
    sinb1 *= -1.;
  Real sin2x = 2.*sinb1*cosb1;
  Real cos2x = 2.*SQR(cosb1)-1.;
  Real norm = DotVec(fp,fp,gcov);
  printf("f: %e %e %e %e\n",fp[IMC0],fp[IMC1],fp[IMC2],fp[IMC3]);
  printf("norm f: %e \n",norm);
  printf("k.f: %e\n",fdotk);
  printf("k.k: %e\n",DotVec(fp,k,gcov));
  printf("fp.e1, fp.e2: %e %e\n",cosb1,cosb2);
  //printf("cos2b,sin2b: %e %e\n",cos2xb,sin2xb);
  if (outsphere) {
    Real cphi = cos(x[IMC3]);
    Real sphi = sin(x[IMC3]);
    Real r = x[IMC1];
    Real cth = cos(x[IMC2]);
    Real sth = sin(x[IMC2]);

    Real Nrr = 2.*SQR(sth)*SQR(spsi*cphi-cpsi*sphi);
    Real Nrth = 2*sth*cth/r*SQR(spsi*cphi-cpsi*sphi);
    Real Nrph = 2./r*((SQR(cpsi)-SQR(spsi))*sphi*cphi-(SQR(cphi)-SQR(sphi))*spsi*cpsi);
    Real Nthth = 2.*SQR(cth/r)*SQR(spsi*cphi-cpsi*sphi);
    Real Nthph = -2.*cth/SQR(r)/sth*((SQR(cpsi)-SQR(spsi))*sphi*cphi
                                     -(SQR(cphi)-SQR(sphi))*spsi*cpsi);
    Real Nphph = 2./SQR(r*sth)*SQR(sphi*spsi+cphi*cpsi);

    Real kr = szet*sth*(cpsi*cphi+spsi*sphi)+czet*cth;
    Real kth = szet*cth/r*(cpsi*cphi+spsi*sphi)-czet*sth/r;
    Real kph = szet/(r*sth)*(spsi*cphi-cpsi*sphi);
    printf("k sph pred: %e %e %e\n",kr,kth,kph);
    printf("k actual: %e %e %e\n",k[IMC1],k[IMC2],k[IMC3]);
    printf("inclination: %g\n",x[IMC2]/M_PI*180.);
    printf("Nr: %e %e %e\n",Nrr,Nrth,Nrph);
    printf("Nth: %e %e\n",Nthth,Nthph);
    printf("Nph: %e\n",Nphph);
    printf("tcord[IMC1]: %e %e %e\n",pphot->polten[IMC1*4+IMC1][ip].real(),
           pphot->polten[IMC1*4+IMC2][ip].real(),pphot->polten[IMC1*4+IMC3][ip].real());
    printf("tcord[IMC2]: %e %e\n",pphot->polten[IMC2*4+IMC2][ip].real(),
           pphot->polten[IMC2*4+IMC3][ip].real());
    printf("tcord[IMC3]: %e\n",pphot->polten[IMC3*4+IMC3][ip].real());

    Real cosd = spsi*sphi+cpsi*cphi;
    Real sind = spsi*cphi-cpsi*sphi;
    Real I = SQR(cosd)+SQR(sind*cth);
    Real Q = SQR(cosd)-SQR(sind*cth);
    Real U = 2*cth*cosd*sind;
    printf("stokes sph pred: %e %e\n",Q/I,U/I);
  }

  printf("stokes actual: %e %e\n",stokes[1],stokes[2]);
  printf("stokes from fp: %e %e\n",cos2x,sin2x);

}

namespace {

void ComputeDeltaGamma(Real x[4], Real k[4], Real a, Real delta[3], Real gamma[3]) {
  //compute the r, theta, and phi parts of the delta and gamma constants for the K and
  // f computations
  Real r = x[IMC1];
  Real theta = x[IMC2];
  Real k_t = k[IMC0];
  Real k_r = k[IMC1];
  Real k_theta = k[IMC2];
  Real k_phi = k[IMC3];

  delta[0] = r*k_t - r*a*SQR(sin(theta))*k_phi;
  delta[1] = a*a*sin(theta)*cos(theta)*k_t - a*cos(theta)*sin(theta)*(r*r+a*a)*k_phi;
  delta[2] = r*a*SQR(sin(theta))*k_r + a*cos(theta)*sin(theta)*(r*r+a*a)*k_theta;
  gamma[0] = a*cos(theta)*k_t - a*a*cos(theta)*SQR(sin(theta))*k_phi;
  gamma[1] = r*(r*r+a*a)*sin(theta)*k_phi - a*r*sin(theta)*k_t;
  gamma[2] = a*a*cos(theta)*sin(theta)*sin(theta)*k_r - r*(r*r+a*a)*sin(theta)*k_theta;

}

void ComputeK(Real x[4], Real k[4], Real f[4], Real a, Real K[2]) {

  Real delta[3];
  Real gamma[3];
  ComputeDeltaGamma(x, k, a, delta, gamma);

  K[0] = delta[0]*f[IMC1] + delta[1]*f[IMC2] + delta[2]*f[IMC3];
  K[1] = gamma[0]*f[IMC1] + gamma[1]*f[IMC2] + gamma[2]*f[IMC3];

}

void Computef(Real x[4], Real k[4], Real gcov[4][4], Real K[2], Real a, Real f[4]){

  Real r = x[IMC1];
  Real theta = x[IMC2];
  Real k_t = k[IMC0];
  Real k_r = k[IMC1];
  Real k_theta = k[IMC2];
  Real k_phi = k[IMC3];

  Real delta[3];
  Real gamma[3];
  ComputeDeltaGamma(x, k, a, delta, gamma);

  Real g_phiphi = gcov[IMC3][IMC3];
  Real g_thetatheta = gcov[IMC2][IMC2];
  Real g_phit = gcov[IMC3][IMC0];
  Real g_rr = gcov[IMC1][IMC1];

  Real N = (gamma[1]*delta[0] - gamma[0]*delta[1])*g_phiphi*k_phi - (gamma[2]*delta[0]
           - gamma[0]*delta[2]) * g_thetatheta * k_theta + (gamma[1]*delta[0]
           - gamma[0]*delta[1])*g_phit*k_t + (gamma[2]*delta[1] - gamma[1]*delta[2])
           * g_rr*k_r;

  f[IMC0] = 0.;
  f[IMC1] = ((gamma[1]*K[0]-delta[1]*K[1])*(g_phiphi*k_phi+g_phit*k_t)
             - (gamma[2]*K[0]-delta[2]*K[1])*g_thetatheta*k_theta)/N;
  f[IMC2] = -((gamma[0]*K[0]-delta[0]*K[1])*(g_phiphi*k_phi + g_phit*k_t)
              - (gamma[2]*K[0]-delta[2]*K[1])*g_rr*k_r)/N;
  f[IMC3] = ((gamma[0]*K[0]-delta[0]*K[1])*g_thetatheta*k_theta
             - (gamma[1]*K[0] - delta[1]*K[1])*g_rr*k_r)/N;

}

} //namespace
