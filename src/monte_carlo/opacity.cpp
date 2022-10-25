//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!  \file opacity.cpp
//!  \brief contains functions related to setting absorption and scattering opacity

// C++ headers
#include <stdexcept>

// Athena++ headers
#include "montecarlo.hpp"
#include "mcutils.hpp"
#include "../globals.hpp"

// SWD: remove or reconfigure
// Lookup table parameters -- chosen to match GRMONTY (Dolence et al. 2009) defaults
#define MINE 1.e-12
#define MAXE 1.e6
#define MINT 1.e-4
#define MAXT 1.e4
#define NE 220
#define NT 80

// Global variable definitions
Real xsect[NE + 1][NT + 1];
Real dle, dlt, lmine, lmint;
static Real kappad = 1.;
static Real albedo = 1.;
//----------------------------------------------------------------------------------------
//! \fn Real NoOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//  \brief return zero for extinction coeffictent

Real NoOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  return 0.;
}


//----------------------------------------------------------------------------------------
//! \fn Real FreeFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//  \brief calculation extinction coefficient for free-free absorption

Real FreeFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  Real &energy = pphot->ep[ip];
  int &i1 = pphot->i1p[ip];
  int &i2 = pphot->i2p[ip];
  int &i3 = pphot->i3p[ip];

  Real ffnrm = 3.692146e8;
  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.67262192369e-24;
  Real h = 6.62607015e-27;
  Real kb = 1.380649e-16;

  //ffnrm *= 12.;  // Added to match the Athena++ prescription

  Real nh = pmcb->rho(i3,i2,i1) / (mp*(1.+4.*heabund));
  Real nhe = nh*heabund;
  Real ne = nh + 2.*nhe;

  Real nu = energy / h;
  Real tgas = pmcb->tgas(i3,i2,i1);
  Real ehnu = exp(-energy / (kb * tgas) );

  Real aff = ffnrm/sqrt(tgas)/pow(nu,3);

  return ne * (nh + 4. * nhe) * aff * (1. - ehnu);

}

//----------------------------------------------------------------------------------------
//! \fn Real DustAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//! \brief calculation extinction coefficient for dust absorption

Real DustAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  int &i1 = pphot->i1p[ip];
  int &i2 = pphot->i2p[ip];
  int &i3 = pphot->i3p[ip];
  Real &energy = pphot->ep[ip];
  Real kapdust = kappad*1.5e13*(1.-albedo)/albedo;
  return kapdust*pmcb->rho(i3,i2,i1);

}

//----------------------------------------------------------------------------------------
//! \fn Real ThomsonOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//  \brief calculation extinction coefficient for Thomson scattering

Real ThomsonOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  int &i1 = pphot->i1p[ip];
  int &i2 = pphot->i2p[ip];
  int &i3 = pphot->i3p[ip];

  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.67262192369e-24;
  Real sigmat = 6.65248e-25;

  Real kappaes = sigmat * (1. + 2.*heabund) / (mp * (1.+4.*heabund) );
  return kappaes * pmcb->rho(i3,i2,i1);
}

//----------------------------------------------------------------------------------------
//! \fn Real ComptonOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//  \brief Returns compton cross section via lookup table

Real ComptonOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  int &i1 = pphot->i1p[ip];
  int &i2 = pphot->i2p[ip];
  int &i3 = pphot->i3p[ip];
  Real &energy = pphot->ep[ip];

  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.67262192369e-24;
  Real sigmat = 6.65248e-25;
  Real kmec2 = 1.68638e-10;
  Real mec2 = 8.18711e-7;

  Real kappa0 = 1./mp/(1. + 4.*heabund) * (1. + 2.*heabund);
  Real theta = pmcb->tgas(i3,i2,i1)* kmec2;
  Real dens = pmcb->rho(i3,i2,i1);
  Real edim = energy / mec2;

  Real sigma0;
  if (theta <= MINT) {
    if (edim <= MINE)
      sigma0 = sigmat;
    else
      sigma0 = KleinNishina(2.*edim) * sigmat;
    return sigma0 * kappa0 * dens;
  } else if (edim <= MINE) {
    return sigmat * kappa0 * dens;
  } else if ( (edim < MAXE) && (theta < MAXT) ) {
    Real le = log10(edim);
    Real lt = log10(theta);
    Real xi = (le - lmine) / dle;
    Real xj = (lt - lmint) / dlt;
    int i = static_cast<int>(xi);
    int j = static_cast<int>(xj);
    xi -= i;
    xj -= j;

    Real lxsect = (1.-xi) * ( (1.-xj) * xsect[i  ][j] + xj * xsect[i  ][j+1] )
                    + xi  * ( (1.-xj) * xsect[i+1][j] + xj * xsect[i+1][j+1] );

    sigma0 = pow(10., lxsect);
  } else {
    std::cout << "Warning: out of range parameter in Compton Opacity: "
              << edim << " " << theta << std::endl
              << "Assuming Thomson cross section." << std::endl;
    sigma0 = sigmat;
  }

  return sigma0 * kappa0 * dens;

}
//----------------------------------------------------------------------------------------
//! \fn Real ResonanceLineOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//  \brief opacity due resonance line for thermal distribution of atoms
//
Real ResonanceLineOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  int &i1 = pphot->i1p[ip];
  int &i2 = pphot->i2p[ip];
  int &i3 = pphot->i3p[ip];
  Real denscatterers;
  Real &energy = pphot->ep[ip];

  Real h = 6.62607015e-27;
  Real tgas = pmcb->tgas(i3, i2, i1);
  Real mass = 1.660538782e-24;

  if (NSCALARS > 0) {
    denscatterers = pmcb->rho(i3,i2,i1) - pmcb->scalars(i3,i2,i1);
  } else {
    denscatterers = pmcb->rho(i3,i2,i1);
  }
  Real nh = denscatterers / mass;
  //if (pmcb->pmy_block->pmy_mesh->ncycle > 1)
  //  printf("%g %g %g %g\n",nh, nh * XsecVoigt(energy / h, tgas),pmcb->rho(i3,i2,i1),pmcb->scalars(i3,i2,i1));
  return nh * XsecVoigt(energy / h, tgas);
}

//----------------------------------------------------------------------------------------
//! \fn Real DustScatteringOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//! \brief calculation extinction coefficient for dust scattering

Real DustScatteringOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  int &i1 = pphot->i1p[ip];
  int &i2 = pphot->i2p[ip];
  int &i3 = pphot->i3p[ip];
  Real &energy = pphot->ep[ip];
  Real kapdust = kappad*1.5e13;

  return kapdust*pmcb->rho(i3,i2,i1);

}

//----------------------------------------------------------------------------------------
//! \fn void GenerateComptonTable(int io)
//! \brief Generates lookup table used by ComptonOpacity()
//
// Computes look up table for integrated cross section of a Maxwellian
// distribution of electrons. Current version uses direct integration as
// in GRMONTY rather than approximate method of PSS.  Values are chosen to
// match GRMONTY (Dolence et al. 2009) defaults

void GenerateComptonTable(int io) {

  if (io > 0) {
    // generate table from scratch
    dle = log10(MAXE / MINE) / NE;
    dlt = log10(MAXT / MINT) / NT;
    lmine = log10(MINE);
    lmint = log10(MINT);

    for (int i = 0; i <= NE; ++i) {
      Real energy = pow(10.,lmine + i * dle);
      for (int j = 0; j <= NT; ++j) {
        Real theta = pow(10.,lmint + j * dlt);
        xsect[i][j] = log10( ComptonCrossSection(energy,theta) );
      }}
    if (io == 2) {
      // open file for output
      FILE *pfile;
      std::string fname;
      fname.assign("comptontable.out");
      std::stringstream msg;
      if ((pfile = fopen(fname.c_str(),"w")) == NULL) {
        msg << "### FATAL ERROR in function [GenerateComptonTable]"
          <<std::endl<< "Output file '" <<fname<< "' could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }

      // write table
      Real data[NT+1];
      for (int i=0; i<=NE; ++i) {
        for (int j=0; j<=NT; ++j) {
          data[j] = xsect[i][j];
        }
        fwrite(data,sizeof(Real),static_cast<size_t>(NT+1),pfile);
      }
      fclose(pfile);
    }
  } else {
    // open file for output
    FILE *pfile;
    std::string fname;
    //if (Globals::my_rank == 0)
    //  std::cout << "Reading in table for Compton cross section." << std::endl;
    fname.assign("comptontable.out");
    std::stringstream msg;
    if ((pfile = fopen(fname.c_str(),"r")) == NULL) {
      msg << "### FATAL ERROR in function [GenerateComptonTable]"
          <<std::endl<< "Input file '" <<fname<< "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // set global variables
    dle = log10(MAXE / MINE) / NE;
    dlt = log10(MAXT / MINT) / NT;
    lmine = log10(MINE);
    lmint = log10(MINT);
    // read table from file
    for (int i = 0; i <= NE; ++i) {
      for (int j = 0; j <= NT; ++j) {
        Real data;
        fread(&data, sizeof(Real), 1, pfile);
        //if (i == 0) printf("%g ",data);
        xsect[i][j] = data;
      }}

    fclose(pfile);
  }

}
//----------------------------------------------------------------------------------------
//! \fn Real ComptonCrossSection(Real energy, Real theta)
//! \brief Computes compton cross section for GenerateComptonTable

Real ComptonCrossSection(Real energy, Real theta) {

  Real sigmat = 6.65248e-25;

  if (theta < MINT) {
    if (energy < MINE)
      return sigmat;
    else
      return KleinNishina(2.*energy) * sigmat;
  }

  Real dmu = 0.05;
  Real dgam = theta * 0.05;
  Real xsect = 0.;
  for (int i=0; i<40; ++i) {
    Real mu = (0.5 + static_cast<Real>(i)) * dmu - 1.;
    for (int j=0; j<240; ++j) {
      Real gamma = 1.0 + (0.5 + static_cast<Real>(j)) * dgam;
      Real velc = sqrt(1. - 1. / SQR(gamma));
      Real x = 2.0 * energy * gamma * (1. - mu * velc);
      xsect += dmu * dgam * Maxwell(theta,gamma) * KleinNishina(x) * (1. - mu * velc);
    }
  }
  return 0.5 * xsect * sigmat;
}

//----------------------------------------------------------------------------------------
//! \fn Real Maxwell(Real theta, Real gamma)
//! \brief Computes maxwell distribution for total_compton_xsect

Real Maxwell(Real theta, Real gamma)
{
  Real K2exp;
  if (theta < 0.01)
    K2exp = sqrt(0.5 * PI * theta);// * (1. + 1.875 * theta);
  else
    K2exp = BessK(2,1. / theta) * exp(1. / theta);

  return gamma * sqrt(SQR(gamma) - 1.) * exp((1. - gamma) / theta) / (theta * K2exp);
}

//----------------------------------------------------------------------------------------
//! \fn Real KleinNishina(Real x)
//! \brief Computes Klein-Nishina correction for total_compton_xsect

Real KleinNishina(Real x)
{

  if (x < 0.001)
    return (1. - x);
  else
    return 0.75 / x * ( (1. - 4. / x * (1. + 2. / x) ) * log(1. + x) +
                        0.5 + 8. / x - 0.5 / SQR(1. + x) );

}

//----------------------------------------------------------------------------------------
//! \fn void InitializeAccelerationOpacity(MonteCarloBlock *pmcb)
//! \brief Computes mean opacity arrays used by acceleration routines

void InitializeAccelerationOpacity(MonteCarloBlock *pmcb) {

  int ip = 0;
  int nx = 500;
  Real x[500], dx[500];

  // Make grid in x=hnu/kT space from 10^-4 to 10^2
  for(int i=0; i<nx; ++i) {
    Real log10x = static_cast<Real>(i)/(static_cast<Real>(nx)-1.0)*6.0;
    x[i] = pow(10.0,log10x)*1.0e-4;
  }

  // compute dx
  dx[0] = x[1]-x[0];
  dx[nx-1] = x[nx-1]-x[nx-2];
  for(int i=1; i<nx-1; i++) {
    dx[i] = 0.5 * (x[i+1]-x[i-1]);
  }

  // Loop over grid zones
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  Photon *pphot; // SWD temp fix -- will not work below
  Real kb = 1.380649e-16;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        Real temp = pmcb->tgas(k,j,i);
        Real dens = pmcb->rho(k,j,i);
        Real Bint = 0.0;
        Real planck = 0.0;
        Real planck_inv = 0.0;
        for(int l=0; l<nx; ++l) {
          Real energy = x[l] * kb * temp;
          Real abs_coef = pmcb->AbsorptionOpacity(pmcb,pphot,ip);
          Real sct_coef = pmcb->ScatteringOpacity(pmcb,pphot,ip);
          Real Bx = pow(x[l],3)/(exp(x[l])-1.0) *dx[l];
          Bint += Bx;
          planck += Bx * abs_coef;
          planck_inv += Bx / (abs_coef+sct_coef);
        }
        pmcb->planck_inv_opacity(k,j,i) = Bint/planck_inv;
        pmcb->planck_opacity(k,j,i) = planck/Bint;
      }}}
}

// SWD: Useful to retain ability to do different lines but should ensure
// computation only occurs once
//----------------------------------------------------------------------------
//! \fn Real ResLinePre()
//! \brief Species dependent prefactor

Real ResLinePre() {

  Real charge = 4.80320427e-10;
  Real melectron = 9.10938215e-28;
  Real clight = 2.99792458e10;
  Real osc_strength = 0.4164;

  return PI*charge*charge / (melectron*clight) * osc_strength;
}

//----------------------------------------------------------------------------
//! \fn Real xsec_lorentzian(Real nu)
//! \brief lorentzian cross section, used for testing

Real XsecLorentzian(Real nu) {

  Real lorwidth = 6.265e8/(4.*PI);
  Real nu0 = 2.468e15;
  Real lineprofile = lorwidth/ PI / ( SQR(nu-nu0) + SQR(lorwidth) );
  Real sigmatot = ResLinePre() * lineprofile;

  return sigmatot;
}

//-------------------------------------------------------------------------
//! \fn Real XsecDoppler(Real nu, Real tgas)
//! \brief Doppler cross section, used for testing

Real XsecDoppler(Real nu, Real tgas) {

  Real kb = 1.380649e-16;
  Real mass = 1.660538782e-24;
  Real vth = sqrt( 2. * kb * tgas / mass);

  Real clight = 2.99792458e10;
  Real lorwidth = 6.265e8/(4.*PI);
  Real nu0 = 2.468e15;

  Real doppwidth = nu0 * vth / clight;
  Real x = (nu-nu0)/doppwidth;
  Real lineprofile = exp(-x*-x) / sqrt(PI) / doppwidth;

  Real sigmatot = ResLinePre() * lineprofile;

  return sigmatot;

}

//------------------------------------------------------------
//! \fn Real XsecVoigt(Real nu, Real tgas)
//! \brief Voigt cross section

Real XsecVoigt(Real nu, Real tgas) {

  Real kb = 1.380649e-16;
  Real mass = 1.660538782e-24;
  Real vth = sqrt( 2. * kb * tgas / mass);

  Real lorwidth = 6.265e8/(4.*PI);
  Real nu0 = 2.468e15;

  Real clight = 2.99792458e10;

  Real doppwidth = nu0 * vth / clight;
  Real a = lorwidth / doppwidth;
  Real x = (nu-nu0)/doppwidth;

  std::complex<double> z(x, a);
  std::complex<double> f = ZetaVoigt(z);

  Real H = f.imag() / sqrt(PI);

  Real lineprofile = H / sqrt(PI) / doppwidth;

  Real sigmatot = ResLinePre() * lineprofile;

  return sigmatot;
}
