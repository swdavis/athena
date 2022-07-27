//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file atmos.cpp
//! \brief A hydrostatic atmosphere that tracks ionization using a passive scalar.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../monte_carlo/montecarlo.hpp"
#include "../monte_carlo/photon.hpp"
#include "../monte_carlo/photonmover.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

Real rin,rout;
Real energy0;
int i1,i2,i3;


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Spherical atmosphere in hydrostatic balance
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  //EnrollUserExplicitSourceFunction(TrackIonization);
  // Add user mesh data block
}

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){
  nuser_var = 3;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  Real kb    = 1.380649e-16;
  Real mp    = 1.6726e-24;
  Real clight = 2.997924589e10;
  Real sigma_pi = 6.3e-18;// cm2

  Real GM   = pin->GetOrAddReal("problem", "GM", 0.);
  Real mmw   = pin->GetOrAddReal("problem", "mmw", 1.); // Mean molecular weight of gas
  Real rin   = pin->GetOrAddReal("problem", "rin", 1.e10); // Radius of base of atmosphere
  Real pbase = pin->GetOrAddReal("problem", "pbase", 10.); // Atmospheric pressure at rin
  Real temp  = pin->GetOrAddReal("problem", "temp", 1.e4); // Isothermal temperature
  Real rhobase = pbase * mmw * mp / kb / temp;
  Real gamma = peos->GetGamma();
  Real gm1   = gamma - 1.0;

  // set up ambient medium at equilibrium for an isothermal atmosphere
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      Real column = 0.;
      for (int i=ie+NGHOST; i>=is-NGHOST; --i) {
        Real r = pcoord->x1v(i);
        phydro->u(IDN,k,j,i) = rhobase * exp(GM * mmw * mp / kb / temp * (1./pcoord->x1v(i) - 1./rin));
        Real rho = phydro->u(IDN,k,j,i);

        // Calculate density where n_H = n_p in optically thin limit
        Real Gamma0 = 1. / (6. * 60. * 60.); // s, photoionization rate coefficient
        Real alpha = 4.18e-13;
        Real neq0 = Gamma0 / alpha;

        // Calculate neutral hydrogen column
        Real n_H = std::pow(0.5 * (std::sqrt(neq0 + 4.*rho/mp) - std::sqrt(neq0)), 2.);
        Real dr = pcoord->dx1f(i);
        column += n_H * dr;

        // Calculate ionization rate, initialize ionization fraction scalar
        Real Gamma = Gamma0 / (1. + std::pow(sigma_pi * column, 1.5)); // Trammell et al 2011, Fig 9 powerlaw
        pscalars->s(0,k,j,i) = 0.5 * (std::sqrt(Gamma) * std::sqrt(4.*alpha*rho/mp + Gamma) / alpha
                               - Gamma / alpha) * mp;
        pscalars->r(0,k,j,i) = pscalars->s(0,k,j,i) / rho;

        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = phydro->u(IDN,k,j,i) * kb * temp / mmw / mp / gm1; // Fluid internal energy density
        phydro->u(IEN,k,j,i) += 0.5 * SQR(phydro->u(IM1,k,j,i)) / phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5 * SQR(phydro->u(IM2,k,j,i)) / phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5 * SQR(phydro->u(IM3,k,j,i)) / phydro->u(IDN,k,j,i);
      }
    }
  }
}

//========================================================================================
//! \fn void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin)
//! \brief Analogous to problem generator but used in support of InitializePhoton
//========================================================================================

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {

  Real mp    = 1.6726e-24;
  Real alpha = 4.18e-13;
  Real kb = 1.380649e-16;
  Real c = 2.99792458e10;
  Real nu0 = 2.468e15;
  Real lorw = 6.265e8/(4.*PI);
  Real h = 6.62607015e-27;

  Real x0 = pin->GetReal("problem","x0");
  Real temp = pin->GetReal("problem","temp");
  Real vth = sqrt( 2. * kb * temp / mp);
  Real dopw = nu0 * vth / c;
  Real energy0 = h * (nu0 + dopw * x0);

  // Set variables
  Real rin = pin->GetReal("problem","rin");
  Real rout = pin->GetReal("problem","rout");

  Real Ndot_tot = 0.;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // Integral over recombination in all zones
        Real np = scalars(0,k,j,i)/mp;
        Real vol = pcoord->vol(k,j,i);
        Ndot_tot += alpha * SQR(np) * vol;
  }}}
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real np = scalars(0,k,j,i)/mp;
        Real vol = pcoord->vol(k,j,i);
        emission(k,j,i) = alpha * SQR(np) * vol / Ndot_tot;
  }}}
}

//========================================================================================
//! \fn void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe)
//! \brief Initializes Photon packets before integration
//========================================================================================

void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe) {

  for (int ip=ips; ip<=ipe; ip++) {

    // Set status flag
    pphot->statp[ip] = EVOLVING;

    // Emit the photon
    // random sample position
    Real rin3 = pow(rin, 3.);
    Real rout3 = pow(rout, 3.);
    Real r0 = pow((rout3 - rin3) * pran->uniform() + rin3, 1./3.);

    // random sample angle
    Real phi = 2. * PI * pran->uniform();
    Real cth = 2. * pran->uniform() - 1.;

    // Set coordinates in sphpol
    pphot->x1p[ip] = r0;
    pphot->x2p[ip] = acos(cth);
    pphot->x3p[ip] = phi;
    pphot->x0p[ip] = 0.; //time

    // Get indices of the grid zone the photon is currently in
    i1 = -1;
    for(int i=is; i<=ie; i++) {
      if ((pphot->x1p[ip] > pcoord->x1f(i)) && (pphot->x1p[ip] <= pcoord->x1f(i+1)))
        i1 = i;
    }
    i2 = -1;
    for(int j=js; j<=je; j++) {
      if ((pphot->x2p[ip] > pcoord->x2f(j)) && (pphot->x2p[ip] <= pcoord->x2f(j+1)))
        i2 = j;
    }
    i3 = -1;
    for(int k=ks; k<=ke; k++) {
      if ((pphot->x3p[ip] > pcoord->x3f(k)) && (pphot->x3p[ip] <= pcoord->x3f(k+1)))
        i3 = k;
    }
    if ((i1 < 0) || (i2 < 0) || (i3 < 0)) {
      std::stringstream msg;
      msg << "### FATAL ERROR in InitializePhoton" << std::endl
          << "Photon position not found within domain." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // Initialize wave vector so that it is parallel with r direction
    pphot->k1p[ip] = 1.0;
    pphot->k2p[ip] = 0.;
    pphot->k3p[ip] = 0.;
    pphot->k0p[ip] = 1. / 2.99792458e10;

    // Set weight according to the emission array, which is the relative number
    // of photons per unit time emitted in each cell
    pphot->wp[ip] = emission(i3,i2,i1);
    pphot->ep[ip] = energy0;

    // Set status flag
    if (pphot->wp[ip] < 0.0)
      pphot->statp[ip] = DESTROYED;
    else
      pphot->statp[ip] = EVOLVING;

    // Initialize the absorption and scattering extinction coefficients
    // to the values appropriate in the emitted zone
    pphot->acp[ip] = AbsorptionOpacity(this,pphot,ip);
    pphot->scp[ip] = ScatteringOpacity(this,pphot,ip);
  } // loop over ip

}
