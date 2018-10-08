//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!  \file opacity.cpp
//   \brief contains functions related to setting absorption and scattering opacity

// Athena++ headers
#include "montecarlo.hpp"
#include "mcutils.hpp"

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

//----------------------------------------------------------------------------------------
//! \fn Real NoOpacity(MonteCarloBlock *pmcb, Photon *pphot)
//  \brief return zero for extinction coeffictent

Real NoOpacity(MonteCarloBlock *pmcb, Photon *pphot) {

  return 0.;
}

//----------------------------------------------------------------------------------------
//! \fn Real FreeFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot)
//  \brief calculation extinction coefficient for free-free absorption

Real FreeFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot) {

  Real ffnrm = 3.692146e8;
  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.6726e-24; 
  Real h =  6.6262e-27;
  Real kb = 1.3807e-16;

  Real nhii = pmcb->rho(pphot->i3,pphot->i2,pphot->i1) / (mp*(1.+4.*heabund));
  Real ne = (1. + 2.*heabund) * nhii;

  Real nu = pphot->energy / h;
  Real temp = pmcb->tgas(pphot->i3,pphot->i2,pphot->i1);
  Real ehnu = exp(-pphot->energy / (kb * temp) );

  Real aff = ffnrm/sqrt(temp)/pow(nu,3);
  //printf("opac: %g %g %g %g\n",temp,pmcb->rho(pphot->i3,pphot->i2,pphot->i1),nhii,aff);
  return ne * nhii * aff * (1. - ehnu);

}

//----------------------------------------------------------------------------------------
//! \fn Real ThomsonOpacity(MonteCarloBlock *pmcb, Photon *pphot)
//  \brief calculation extinction coefficient for Thomson scattering

Real ThomsonOpacity(MonteCarloBlock *pmcb, Photon *pphot) {

  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.6726e-24; 
  Real sigmat = 6.65248e-25;

  Real kappaes = sigmat * (1. + 2.*heabund) / (mp * (1.+4.*heabund) );
  return kappaes * pmcb->rho(pphot->i3,pphot->i2,pphot->i1);
}

//----------------------------------------------------------------------------------------
//! \fn Real ComptonOpacity(MonteCarloBlock *pmcb, Photon *pphot)
//  \brief Returns compton cross section via lookup table

Real ComptonOpacity(MonteCarloBlock *pmcb, Photon *pphot) {
  
  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.6726e-24; 
  Real sigmat = 6.65248e-25;
  Real kmec2 = 1.68638e-10;
  Real mec2 = 8.18711e-7;
  
  Real kappa0 = 1./mp/(1. + 4.*heabund) * (1. + 2.*heabund);
  Real theta = pmcb->tgas(pphot->i3,pphot->i2,pphot->i1)* kmec2;
  Real dens = pmcb->rho(pphot->i3,pphot->i2,pphot->i1);
  Real edim = pphot->energy / mec2;

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
//! \fn void GenerateComptonTable(void)
//  \brief Generates lookup table used by ComptonOpacity()
//
// Computes look up table for integrated cross section of a Maxwellian 
// distribution of electrons. Current version uses direct integration as 
// in GRMONTY rather than approximate method of PSS.  Values are chosen to 
// match GRMONTY (Dolence et al. 2009) defaults

void GenerateComptonTable(void) {

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

}
//----------------------------------------------------------------------------------------
//! \fn Real ComptonCrossSection(Real energy, Real theta)
//  \brief Computes compton cross section for GenerateComptonTable

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
//  \brief Computes maxwell distribution for total_compton_xsect

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
//  \brief Computes Klein-Nishina correction for total_compton_xsect

Real KleinNishina(Real x)
{

  if (x < 0.001)
    return (1. - x);
  else
    return 0.75 / x * ( (1. - 4. / x * (1. + 2. / x) ) * log(1. + x) +
                        0.5 + 8. / x - 0.5 / SQR(1. + x) );

}

