//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mc_sphere_lya.cpp
//  \brief Problem generator for monte carlo through a uniform isothermal sphere on a
//         Cartesian grid with lyman alpha scattering
//
//========================================================================================

// C/C++ headers
#include <iostream>
#include <stdexcept>
#include <string>

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
  // Global variables
  Real rad0;
  Real energy0;

  enum UserEstimatorIndex {
    PATH_ENERGY = 0,
    FLUX1,
    FLUX2,
    FLUX3,
    FORCE_MOMENT1,
    FORCE_MOMENT2,
    FORCE_MOMENT3,
    PATH_EXTINCTION,
    NUM_ESTIMATORS
  };

  // function headers
  bool LocateOriginCell(MCCoord *pcoord, int is, int ie, int js, int je, int ks, int ke,
                        int &i1start, int &i2start, int &i3start);
  Real DistanceToSphere(MonteCarloBlock *pmcb, Photon *phot, int ip);
  void AccumulateUserEstimators(MonteCarloBlock *pmcb, Photon *phot,
                                PhotonPusher *ppusher, int ip);

  // Find the cell containing the source.  The lower face is excluded so that
  // an origin on a MeshBlock face belongs to exactly one of the neighboring blocks.
  bool LocateOriginCell(MCCoord *pcoord, int is, int ie, int js, int je, int ks, int ke,
                        int &i1start, int &i2start, int &i3start) {
    i1start = -1;
    for (int i=is; i<=ie; ++i) {
      if ((0. > pcoord->x1f(i)) && (0. <= pcoord->x1f(i+1))) i1start = i;
    }

    i2start = -1;
    for (int i=js; i<=je; ++i) {
      if ((0. > pcoord->x2f(i)) && (0. <= pcoord->x2f(i+1))) i2start = i;
    }

    i3start = -1;
    for (int i=ks; i<=ke; ++i) {
      if ((0. > pcoord->x3f(i)) && (0. <= pcoord->x3f(i+1))) i3start = i;
    }

    return (i1start >= 0) && (i2start >= 0) && (i3start >= 0);
  }

  // Return the forward ray distance from an interior point to the spherical surface.
  Real DistanceToSphere(MonteCarloBlock *, Photon *pphot, int ip) {
    const Real b = pphot->x1p[ip] * pphot->k1p[ip]
                 + pphot->x2p[ip] * pphot->k2p[ip]
                 + pphot->x3p[ip] * pphot->k3p[ip];
    Real discriminant = SQR(b) + SQR(rad0)
                      - SQR(pphot->x1p[ip])
                      - SQR(pphot->x2p[ip])
                      - SQR(pphot->x3p[ip]);
    discriminant = (discriminant > 0.) ? discriminant : 0.;
    const Real distance = -b + sqrt(discriminant);
    return (distance > 0.) ? distance : 0.;
  }
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief monte carlo test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  if (std::string(COORDINATE_SYSTEM) != "cartesian") {
    throw std::runtime_error("mc_sphere_lya requires Cartesian coordinates");
  }

  Real rideal = 8.314e7;
  Real c = 2.99792458e10;
  Real temp = pin->GetReal("problem","temp");
  Real tau = pin->GetReal("problem","tau");
  Real rad0 = pin->GetReal("problem","radius");
  Real vel = pin->GetOrAddReal("problem","velocity",0.);
  Real gamma = peos->GetGamma();
  vel *= c;

  const Real sigma0 = XsecVoigt(MCConstants::nu_lya, temp);
  const Real nH = tau / (rad0 * sigma0);
  const Real rho = nH * MCConstants::mH_cgs;
  printf("sigma0: %g nH: %g rho: %g radius: %g tau: %g\n",
         sigma0, nH, rho, rad0, tau);
  // density is non-zero only in sphere
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = rho*vel;
        phydro->u(IEN,k,j,i) = rideal*rho*temp/(gamma-1.0);
      }
    }
  }

  // add kinetic energy
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
      }}}
}

//========================================================================================
//! \fn void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin)
//! \brief Initializes user data specific to MonteCarlo class
//========================================================================================

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){

  if (pin->DoesParameterExist("problem", "time")) {
    throw std::runtime_error(
        "mc_sphere_lya does not support problem/time; photons evolve until spatial escape");
  }

  // MRW hook
  if (acceleration) {
    throw std::runtime_error(
        "mc_sphere_lya does not support montecarlo/acceleration; MRW resonance "
        "bookkeeping and estimators are not yet validated");
  }

  nuser_var = NUM_ESTIMATORS;
  EnrollUserEscapeDistance(DistanceToSphere);
  EnrollUserWorkInMove(AccumulateUserEstimators);
}

//========================================================================================
//! \fn void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin)
//! \brief Analogous to problem generator but used in support of InitializePhoton
//========================================================================================

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {

  Real x0 = pin->GetReal("problem","x0");
  Real temp = pin->GetReal("problem","temp");
  const Real vth = sqrt(2. * MCConstants::kb_cgs * temp / MCConstants::mH_cgs);
  const Real dopw = MCConstants::nu_lya * vth / MCConstants::c_cgs;
  energy0 = MCConstants::h_cgs * (MCConstants::nu_lya + dopw * x0);

  rad0 = pin->GetReal("problem","radius");

  int i1start, i2start, i3start;
  if (!LocateOriginCell(pcoord, is, ie, js, je, ks, ke,
                        i1start, i2start, i3start)) {
    nphremain = 0;
    nphrun = 0;
  } else {
    // Set number of samples per block because emmision_flag is set to EMISNONE
    nphremain = pin->GetInteger("montecarlo", "nphot");
    nphrun = 0;
  }
}


//========================================================================================
//! \fn void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe, int etype)
//! \brief Initializes Photon packets before integration
//========================================================================================

void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe, int etype) {

  int i1start, i2start, i3start;
  if (!LocateOriginCell(pcoord, is, ie, js, je, ks, ke,
                        i1start, i2start, i3start)) {
    throw std::runtime_error("InitializePhoton called on a MeshBlock without the origin");
  }

  for (int ip=ips; ip<=ipe; ip++) {

    for (int n=0; n<pmy_mc->nuser_var; ++n)
      pphot->user[n][ip] = 0.;
    pphot->nscp[ip] = 0;

    // Set status flag
    pphot->statp[ip] = EVOLVING;
    // Initialize photon at the origin
    pphot->i1p[ip] = i1start;
    pphot->i2p[ip] = i2start;
    pphot->i3p[ip] = i3start;
    pphot->x1p[ip] = 0.;
    pphot->x2p[ip] = 0.;
    pphot->x3p[ip] = 0.;
    pphot->x0p[ip] = 0.; // path length

    // Generate initial angle parameters
    Real phi = 2. * PI * pran->uniform();
    Real cphi = cos(phi);
    Real sphi = sin(phi);
    Real cth = 2. * pran->uniform() - 1.;
    Real sth = sqrt(1. - SQR(cth));
    // Initialize wave vector with isotropic distribution
    pphot->k1p[ip] = sth*cphi;
    pphot->k2p[ip] = sth*sphi;
    pphot->k3p[ip] = cth;
    // k0p is the photon energy (set via ep); the light-travel time bookkeeping
    // below now divides by c explicitly instead of stashing 1/c in k0p.

    // Initialize Photon weights, energy, direction, polarization
    // TODO: Make this an input rather than hardcoding
    Real target_lum = 1.e20;
    pphot->wp[ip] = target_lum / energy0;
    pphot->ep[ip] = energy0;

    // Evolve until the packet crosses the spherical escape surface.
    pphot->dtp[ip] = HUGE_NUMBER;

    // Initialize Stokes vector
    if (pphot->polarized) {
      pphot->sip[ip] = 1.0;
      pphot->sqp[ip] = 0.0;
      pphot->sup[ip] = 0.0;
      pphot->svp[ip] = 0.0;
    }

    // Set status flag
    if (pphot->wp[ip] < 0.0)
      pphot->statp[ip] = DESTROYED;
    else
      pphot->statp[ip] = EVOLVING;

    // Initialize the absorption and scattering extinction coefficients
    // to the values appropriate in the emitted zone
    pphot->acp[ip] = AbsorptionOpacity(this,pphot,ip);
    pphot->scp[ip] = ScatteringOpacity(this,pphot,ip);
    // pphot->PrintPhoton(ip);
  } // loop over ip
}



namespace {

// Accumulate user estimators for the completed movement.  The Cartesian pusher has already
// truncated a final segment at the sphere and marked the packet escaped when applicable.
void AccumulateUserEstimators(MonteCarloBlock *pmcb, Photon *pphot,
                              PhotonPusher *ppusher, int ip) {

  const Real k1 = pphot->k1p[ip];
  const Real k2 = pphot->k2p[ip];
  const Real k3 = pphot->k3p[ip];

  const Real dl = ppusher->dl;

  // Match the lab-frame path-length estimator used by AccumulateMoments.  k0p now aliases
  // ep, so it must not appear as an additional factor in these non-relativistic moments.
  const Real c_cgs = MCConstants::c_cgs;
  const Real dl_cgs = dl * pmcb->l_cgs;
  const Real weight = pphot->wp[ip] * pphot->ep[ip] * dl_cgs / c_cgs;
  const Real extinction = pphot->acp[ip] + pphot->scp[ip];

  pphot->user[PATH_ENERGY][ip] += weight;
  pphot->user[FLUX1][ip] += weight * k1 * c_cgs;
  pphot->user[FLUX2][ip] += weight * k2 * c_cgs;
  pphot->user[FLUX3][ip] += weight * k3 * c_cgs;
  pphot->user[FORCE_MOMENT1][ip] += extinction * weight * k1;
  pphot->user[FORCE_MOMENT2][ip] += extinction * weight * k2;
  pphot->user[FORCE_MOMENT3][ip] += extinction * weight * k3;
  // Together with PATH_ENERGY, this forms the scalar path-weighted mean extinction.
  pphot->user[PATH_EXTINCTION][ip] += extinction * weight;
}

} //namespace
