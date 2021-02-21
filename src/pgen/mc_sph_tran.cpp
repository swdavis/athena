//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mctest.cpp
//  \brief Test movement through spherical grid 
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


#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

namespace {

void FinalPositionSphericalPolar(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot,
				 Real &rf, Real &thf, Real &phf);

// Global variables
static Real rfp,thfp, phfp;
static Real error_sum = 0.;
static int nerror = 0;
} // namespace

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

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin) {
  nuser_var = 1;
}

void MonteCarloBlock::InitializePhoton(Photon *pphot) {

  MCCoord *pco = pcoord;

  // Set status flag

  pphot->status = EVOLVING;

  // Choose random intial position, weights, energy, and direction
  // for photon emission.  In this version an equal number of photons
  // is emitted in  each grid zone.  The relative emission from each grid 
  // zone is then accounted for by a weighting factor cweight. 

  Real nx1 = static_cast<Real>(ie-is+1);
  Real nx2 = static_cast<Real>(je-js+1);
  Real nx3 = static_cast<Real>(ke-ks+1);

  pphot->i1 = static_cast<int>(pran->uniform()*nx1)+is;
  pphot->i2 = static_cast<int>(pran->uniform()*nx2)+js;
  pphot->i3 = static_cast<int>(pran->uniform()*nx3)+ks;

  // cweight is a constant weighting factor which accounts for the
  // emissivity of the grid zone in which the photon was emitted
  if (zone_weight_flag) {
    pphot->eweight = 1.0;
    pphot->weight = 1.0;
  }


  // Obtain initial position within zone
  GetZonePosition(pphot,pran,pcoord);
  pphot->x[IMC0] = 1.0;
  Real phi = 2. * PI * pran->uniform();
  Real cph = cos(phi);
  Real sph = sin(phi);

  Real cth = 2. * pran->uniform() - 1.;
  Real sth = sqrt(1. - SQR(cth));

  // Initialize wave vector with isotropic distribution
  pphot->k[IMC0] = 1.;
  pphot->k[IMC1] = sth*cph;
  pphot->k[IMC2] = sth*sph;
  pphot->k[IMC3] = cth;

  // Compute cartesian (for testing) SWD: Remove
  // pphot->kcart[0] = pphot->k[0]*sth*cph + pphot->k[1]*cth*cph - pphot->k[2]*sph;
  //pphot->kcart[1] = pphot->k[0]*sth*sph + pphot->k[1]*cth*sph + pphot->k[2]*cph;
  //pphot->kcart[2] = pphot->k[0]*cth - pphot->k[1]*sth;

  // Convert k unit vector to k^\alpha
  if (general_mover_flag) {
    pphot->k[IMC2] /= pphot->x[IMC1];
    pphot->k[IMC3] /= (pphot->x[IMC1]*sin(pphot->x[IMC2]));
  }
  for(int i=0; i<4; i++) pphot->dk[i] = 0.;
  // Initialize Stokes vector
  pphot->stokes[0] = 1.0;
  pphot->stokes[1] = 0.0;
  pphot->stokes[2] = 0.0;

  // Obtain intitial energy, polarization, direction and weight
  // Utilize free-free emission function in emission.cpp
  pphot->energy = 1.;
 

  // initialize kcart
  //pmover->CurvalinearToCartesian(pphot);
  
  if (pphot->weight < 0.0) pphot->status = DESTROYED;

  // Initialize the absorption and scattering extinction coefficients
  pphot->abs_coef = 0.0;
  pphot->sct_coef = 0.0;

  // Compute predicted final photon positions for comparison
  FinalPositionSphericalPolar(this,pco,pphot,rfp,thfp,phfp);
  //pphot->PrintPhoton();
}

void MonteCarloBlock::FinalizePhoton(Photon *pphot) {
 
  if (pphot->status == ESCAPED) {
    Real xp = rfp*sin(thfp)*cos(phfp);
    Real yp = rfp*sin(thfp)*sin(phfp);
    Real zp = rfp*cos(thfp);
    Real xf =  pphot->x[IMC1]*sin(pphot->x[IMC2])*cos(pphot->x[IMC3]);
    Real yf =  pphot->x[IMC1]*sin(pphot->x[IMC2])*sin(pphot->x[IMC3]);
    Real zf =  pphot->x[IMC1]*cos(pphot->x[IMC2]);

    // Save as user variable for photon list
    Real rf = pcoord->x1f(ie+1);
    Real error = sqrt(SQR(xf-xp)+SQR(yf-yp)+SQR(zf-zp))/rf;
    nerror++;
    error_sum += error;
    //printf("error: %g %g %g %g %g %g %g\n",xp,yp,zp,xf,yf,zf,error/rf);
    pphot->user_var[0] = error;
    if (nerror == cadence)
      printf("Mean error: %g\n",error_sum/static_cast<Real>(nerror));
  }

}


namespace {
void FinalPositionSphericalPolar(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot,
				 Real &rf, Real &thf, Real &phf) {

  Real r = pphot->x[IMC1];
  Real cth = cos(pphot->x[IMC2]);
  Real sth = sin(pphot->x[IMC2]);
  Real cph = cos(pphot->x[IMC3]);
  Real sph = sin(pphot->x[IMC3]);

  Real kr, kth, kph;
  if (pmcb->general_mover_flag) {
    kr = pphot->k[IMC1];
    kth = r * pphot->k[IMC2];
    kph = r * sth * pphot->k[IMC3];
  } else {
    kr = pphot->k[IMC1];
    kth = pphot->k[IMC2];
    kph = pphot->k[IMC3];
  }
  // Convert to cartesian
  Real kx = kr * sth * cph + kth * cth*cph - kph * sph;
  Real ky = kr * sth * sph + kth * cth*sph + kph * cph;
  Real kz = kr * cth - kth * sth;
  
  // Outer boundary is r = rf -- find dlr to this boundary
  rf = pco->x1f(pmcb->ie+1);
  Real ndr0 = pphot->x[0] * (sth * (kx * cph + ky * sph) + kz * cth);
  Real det = 1.0 + (SQR(rf) - SQR(pphot->x[0])) / SQR(ndr0);
  Real dlr1 = ndr0 * (sqrt(det) - 1.0);
  Real dlr2 = -ndr0 * (sqrt(det) + 1.0);

  Real dl;
  if (dlr1 > 0.0) {
    if (dlr2 > 0.0) {
      std::cout << "Warning: both roots positive in FinalPositionSphericalPolar: "
		<< dlr1 << " " << dlr2 << std::endl;
    } else {
      dl = dlr1;
    }
  } else if (dlr2 > 0.0) {
    dl = dlr2;
  }

  // Compute other boundary positions
  //theta
  Real zf = r * cth + kz * dl;
  thf = acos(zf / rf);

  //phi
  Real xf = r * sth * cph + kx * dl;
  Real yf = r * sth * sph + ky * dl;
  phf = atan2(yf,xf);
  if (phf < 0.0)
    phf += 2.*PI;
}

} // namespace 
