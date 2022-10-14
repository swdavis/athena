//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mc_isoth.cpp
//! \brief Problem generator for monte carlo isothermal atmosphere
//
//========================================================================================

// C++ headers
#include <iostream> // SWD: temporary for testing

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

#if !MONTE_CARLO_ENABLED
#error "This problem requires monte carlo"
#endif

namespace {
  // Global variables
  int i1, i2, i3;
  Real nemit;
  bool tnorm;
  Real logemin, logemax;
  Real DensityProfile(Real x, Real xl, Real xh, Real taul, Real tauh, Real kap);
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief monte carlo test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Determine density via optical depth or constant density
  bool constdens = pin->GetOrAddBoolean("problem","constdens",false);
  Real rho, tau, taumin, taumax;
  if (constdens) {
    rho = pin->GetOrAddReal("problem","dens",-1.);
    tau = pin->GetOrAddReal("problem","tau",-1.);
  } else {
    taumin = pin->GetReal("problem","taumin");
    taumax = pin->GetReal("problem","taumax");
  }

  Real heabund = 0.09; //hardcode for now
  Real mp = 1.6726e-24;
  Real sigmat = 6.65248e-25;
  Real kappaes = sigmat * (1. + 2.*heabund) / (mp * (1.+4.*heabund) );
  if (tau > 0.) {
    Real length;
    if (COORDINATE_SYSTEM == "cartesian") {
      Real xlow = pin->GetReal("mesh","x3min");
      Real xhigh = pin->GetReal("mesh","x3max");
      length = xhigh-xlow;
    } else {
      length = pin->GetReal("mesh","x1max");
    }
    rho = tau / (kappaes * length);
  }

  Real xlow, xhigh;
  if (!constdens) {
    if (COORDINATE_SYSTEM == "cartesian") {
      xlow = pin->GetReal("mesh","x3min");
      xhigh = pin->GetReal("mesh","x3max");
    } else {
      bool radial = pin->GetOrAddBoolean("problem","radial","true");
      if (radial) {
        xlow = pin->GetReal("mesh","x1min");
        xhigh = pin->GetReal("mesh","x1max");
      }
    }
  }

  // Assume constant velocity provided as fraction of speed of light
  Real vel = pin->GetOrAddReal("problem","velocity",0.);
  Real c = 2.99792458e10;
  vel *= c;

  // Assume constant temperature and ideal gas
  Real gamma = peos->GetGamma();
  Real rideal = 8.314e7;
  Real tgas = pin->GetReal("problem","temp");
  // Set initial conditions
  if (COORDINATE_SYSTEM == "cartesian") {
    // density varies in the z direction

    for (int k=ks; k<=ke; k++) {
      Real x1 = pcoord->x3v(k);
      if (!constdens) rho = DensityProfile(x1,xlow,xhigh,taumin,taumax,kappaes);
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = rho*vel;
          phydro->u(IEN,k,j,i) = rideal*rho*tgas/(gamma-1.0);
        }
      }
    }
  } else if  (COORDINATE_SYSTEM == "spherical_polar") {
    bool radial = pin->GetOrAddBoolean("problem","radial","true");

    if (radial) {
      // density varies in the r direction
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            Real x1 = pcoord->x1v(k);
            if (!constdens) rho = DensityProfile(x1,xlow,xhigh,taumin,taumax,kappaes);
            phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = rho*vel;
            phydro->u(IM2,k,j,i) = 0.0;
            phydro->u(IM3,k,j,i) = 0.0;
            phydro->u(IEN,k,j,i) = rideal*rho*tgas/(gamma-1.0);
          }
        }
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
//! \fn void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe)
//! \brief Initializes Photon packets before integration
//========================================================================================

void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe) {

  // Get meshblock dimensions
  Real nx1 = static_cast<Real>(ie-is+1);
  Real nx2 = static_cast<Real>(je-js+1);
  Real nx3 = static_cast<Real>(ke-ks+1);


  for (int ip=ips; ip<=ipe; ip++) {

    if (pmy_mc->emission_eqwt) {
      bool this_zone = false;
      while (!this_zone) {
        if (nemit > 1.) {
          pphot->i1p[ip] = i1+is;
          pphot->i2p[ip] = i2+js;
          pphot->i3p[ip] = i3+ks;
          this_zone = true;
          nemit -= 1.;
        } else if (nemit > 0.) {
          if (pran->uniform() < nemit) {
            pphot->i1p[ip] = i1+is;
            pphot->i2p[ip] = i2+js;
            pphot->i3p[ip] = i3+ks;
            this_zone = true;
          }
          nemit -= 1.;
        } else {
          this_zone = false;
          i3++;
          if (i3 >= nx3) {
            i3 = 0;
            i2++;
            if (i2 >= nx2) {
              i2 = 0;
              i1++;
              if (i1 >= nx1)
                i1 = 0;
            }
          }
          nemit = emission(i3+ks,i2+js,i1+is) / weight;
          //printf("nemit: %g %d %d %d %g %g\n",nemit,i3,i2,i1,emission(i3,i2,i1),weight);
        }
      } // end while (!this_zone)
      pphot->wp[ip] = weight;
    } else {
      // Randomly assign emission zone
      pphot->i1p[ip] = i1 = static_cast<int>(pran->uniform()*nx1)+is;
      pphot->i2p[ip] = i2 = static_cast<int>(pran->uniform()*nx2)+js;
      pphot->i3p[ip] = i3 = static_cast<int>(pran->uniform()*nx3)+ks;

      // Set weight according to the emission array, which is the relative number of
      // photons emitted in each cell
      pphot->wp[ip] = emission(i3,i2,i1);

    }
    // Obtain initial position within zone
    GetZonePosition(pphot,pran,pcoord,ip);

    // Set maximum integration time
    pphot->dtp[ip] = pphot->pmy_mcb->pmy_mc->tmax;

    // Obtain intitial energy, polarization, direction and weight
    // Utilize free-free emission function in emission.cpp
    if(tnorm) {
      Real logtg = log(tgas(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip]));
      PhotonEmitFreeFree(this,pphot,logemin+logtg,logemax+logtg,ip);
    } else{
      PhotonEmitFreeFree(this,pphot,logemin,logemax,ip);
    }

    // Convert k unit vector to k^\alpha
    if (pmy_mc->general_mover_flag) {
      pphot->k0p[ip] = 1.;
      pphot->k2p[ip] /= pphot->x1p[ip];
      pphot->k3p[ip] /= (pphot->x1p[ip]*sin(pphot->x2p[ip]));
      pphot->dk0p[ip] = 0.;
      pphot->dk1p[ip] = 0.;
      pphot->dk2p[ip] = 0.;
      pphot->dk3p[ip] = 0.;
    }
    //printf("%g %g \n",pphot->ep[ip],pphot->wp[ip]);
    // Set status flag
    if (pphot->wp[ip] < 0.0)
      pphot->statp[ip] = DESTROYED;
    else
      pphot->statp[ip] = EVOLVING;

    // initialize scattering number
    pphot->nscp[ip] = 0;

    // Initialize the absorption and scattering extinction coefficients
    // to the values appropriate in the emitted zone
    pphot->acp[ip] = AbsorptionOpacity(this,pphot,ip);
    pphot->scp[ip] = ScatteringOpacity(this,pphot,ip);

    //printf("start: %d %g %g %d\n",pphot->i3p[ip],pphot->wp[ip],pphot->ep[ip],pphot->statp[ip]);
  }
  //pphot->nphot++;

}

//========================================================================================
//! \fn void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin)
//! \brief Analogous to problem generator but used in support of InitializePhoton
//========================================================================================

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {


  // Set the energy boundaries for free-free emission
  tnorm = pin->GetOrAddBoolean("problem","tnorm",false);
  if (tnorm) {
    // interpret as xmin/xmax with x=E/(kb*T)
    const Real kb = 1.380649e-16;
    logemin = log(kb*pin->GetReal("problem", "emin"));
    logemax = log(kb*pin->GetReal("problem", "emax"));
  } else {
    // interpret as emin/emax in eV
    const Real everg = 1.6021772e-12;
    logemin = log(everg*pin->GetReal("problem", "emin"));
    logemax = log(everg*pin->GetReal("problem", "emax"));
  }
  if (pmy_mc->emission_eqwt) {
    i1 = -1; i2 = -1; i3 = -1;
    nemit = 0.;
  }
}

//========================================================================================
//! \fn void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin)
//! \brief Initializes user data specific to MonteCarlo class
//========================================================================================

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin){

  nuser_var = 1;

}

namespace {
Real DensityProfile(Real x, Real xl, Real xh, Real taul, Real tauh, Real kap) {

  Real l0 = (xh-xl) / log(tauh/taul);
  return taul/l0/kap*exp((xh-x)/l0);
}

}
