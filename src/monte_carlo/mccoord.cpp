//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//  \brief implementation of functions in class MCCoord

// SWD: General notes:
// * need to add documentation
// * remove zeroing of metric,connection at top?
// * remove extraneous variable definitions
// * rename X to x to better match athena style
// * add mass parameter to black hole metrics
// * conversion to cgs units for black hole metrics? 

// Athena++ headers
#include "../athena.hpp"
#include "mccoord.hpp"

// constructor
MCCoord::MCCoord(Coordinates *pcoord, MonteCarloBlock *pmcb) {

  x1f.InitWithShallowCopy(pcoord->x1f);
  x2f.InitWithShallowCopy(pcoord->x2f);
  x3f.InitWithShallowCopy(pcoord->x3f);

  // Needed for black hole coordinates
  if (GENERAL_RELATIVITY) {
    bh_mass_ = pcoord->GetMass();
    bh_spin_ = pcoord->GetSpin();
  }

  // Allocate volume array
  int ncells1 = pmcb->nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmcb->nx2 > 1) ncells2 = pmcb->nx2 + 2*(NGHOST);
  if (pmcb->nx3 > 1) ncells3 = pmcb->nx3 + 2*(NGHOST);
  vol.NewAthenaArray(ncells3,ncells2,ncells1);
  // Initialize volume array
  for (int k=pmcb->ks; k<=pmcb->ke; ++k) {
    for (int j=pmcb->js; j<=pmcb->je; ++j) {
      for (int i=pmcb->is; i<=pmcb->ie; ++i) {
        vol(k,j,i) = pcoord->GetCellVolume(k,j,i);
      }}}
  if(pmcb->acceleration) {
    dmin.NewAthenaArray(ncells3,ncells2,ncells1);
    AthenaArray<Real> dw1,dw2,dw3;
    dw1.NewAthenaArray(ncells1);
    dw2.NewAthenaArray(ncells2);
    dw3.NewAthenaArray(ncells3);
    // Initialize dmin array
    for (int k=pmcb->ks; k<=pmcb->ke; ++k) {
      for (int j=pmcb->js; j<=pmcb->je; ++j) {
	pcoord->CenterWidth1(k,j,pmcb->is,pmcb->ie,dw1);
	pcoord->CenterWidth2(k,j,pmcb->is,pmcb->ie,dw2);
	pcoord->CenterWidth3(k,j,pmcb->is,pmcb->ie,dw3);
	for (int i=pmcb->is; i<=pmcb->ie; ++i) {
	  Real dmin0 = std::min(dw1(i),dw2(i));
	  dmin(k,j,i) = std::min(dmin0,dw3(i));
	  //printf("%d %d %d %g\n",k,j,i,dmin(k,j,i));
	}
      }}
    dw1.DeleteAthenaArray();
    dw2.DeleteAthenaArray();
    dw3.DeleteAthenaArray();
  }
}

// constructor
MCCoord::MCCoord(int ncells1, int ncells2, int ncells3, bool acc) {


  x1f.NewAthenaArray(ncells1+1);
  x2f.NewAthenaArray(ncells2+1);
  x3f.NewAthenaArray(ncells3+1);

  vol.NewAthenaArray(ncells3,ncells2,ncells1);
  acceleration = acc;
  if (acceleration) {
    dmin.NewAthenaArray(ncells3,ncells2,ncells1);
  }
}


// destructor
MCCoord::~MCCoord() {

  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();
  vol.DeleteAthenaArray();
  if (acceleration)
    dmin.DeleteAthenaArray();
}

void MCCoord::Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]) {

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      if (i == j) {
	if (i == IMC0)
	  gcov[i][i] = -1.;
	else
	  gcov[i][i] = 1.;
      } else
	gcov[i][j] = 0.;
    }
  }
}

void MCCoord::InverseMetric(Real x[NCOORD], Real gcon[NCOORD][NCOORD]) {

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      if (i == j) {
	if (i == IMC0)
	  gcon[i][i] = -1.;
	else
	  gcon[i][i] = 1.;
      } else
	gcon[i][j] = 0.;
    }
  }
}

void MCCoord::Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]) {

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      for (int k = 0; k < NCOORD; k++) {
	gamma[i][j][k] = 0.;
      }
    }
  }
}

// constructor
MCCartesian::MCCartesian(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

// constructor
MCCartesian::MCCartesian(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

// destructor
MCCartesian::~MCCartesian() {

}

// constructor
MCSphericalPolar::MCSphericalPolar(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

// constructor
MCSphericalPolar::MCSphericalPolar(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

// destructor
MCSphericalPolar::~MCSphericalPolar() {

}

void MCSphericalPolar::Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]) {

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      gcov[i][j] = 0;
    }
  }

  Real sth = sin(x[IMC2]);
  gcov[IMC0][IMC0] = -1; 
  gcov[IMC1][IMC1] = 1;
  gcov[IMC2][IMC2] = x[IMC1] * x[IMC1];
  gcov[IMC3][IMC3] = x[IMC1] * x[IMC1] * sth * sth;

}

void MCSphericalPolar::Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]) {


  for(int i = 0; i < NCOORD; i++) {
    for(int j = 0; j < NCOORD; j++) {
      for(int k = 0; k < NCOORD; k++) {
	gamma[i][j][k]=0;
      }
    }
  }

  Real sth = sin(x[IMC2]);
  Real cth = cos(x[IMC2]);

  gamma[IMC1][IMC2][IMC2] = -x[IMC1];
  gamma[IMC1][IMC3][IMC3] = -x[IMC1]*sth*sth;
  gamma[IMC2][IMC1][IMC2] = 1./x[IMC1];
  gamma[IMC2][IMC2][IMC1] = 1./x[IMC1];
  gamma[IMC2][IMC3][IMC3] = -sth*cth;
  gamma[IMC3][IMC1][IMC3] = 1./x[IMC1];
  gamma[IMC3][IMC2][IMC3] = cth/sth;
  gamma[IMC3][IMC3][IMC1] = 1./x[IMC1];
  gamma[IMC3][IMC3][IMC2] = cth/sth;

}

// constructor
MCCylindrical::MCCylindrical(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

// constructor
MCCylindrical::MCCylindrical(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

// destructor
MCCylindrical::~MCCylindrical() {

}

void MCCylindrical::Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]) {

  for(int i = 0; i < NCOORD; i++) {
    for(int j = 0; j < NCOORD; j++) {
      gcov[i][j] = 0;
    }
  }
  gcov[IMC0][IMC0] = -1;
  gcov[IMC1][IMC1] = 1;
  gcov[IMC2][IMC2] = x[IMC1] * x[IMC1];
  gcov[IMC3][IMC3] = 1;

}

void MCCylindrical::Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]) {

  for(int i = 0; i < NCOORD; i++) {
    for(int j = 0; j < NCOORD; j++) {
      for(int k = 0; k < NCOORD; k++) {
	gamma[i][j][k]=0;
      }
    }
  }

  gamma[IMC1][IMC2][IMC2] = -x[IMC1];
  gamma[IMC2][IMC1][IMC2] = 1./x[IMC1];
  gamma[IMC2][IMC2][IMC1] = 1./x[IMC1];

}


// constructor
MCKerrSchild::MCKerrSchild(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

// constructor
MCKerrSchild::MCKerrSchild(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

// destructor
MCKerrSchild::~MCKerrSchild() {

}

void MCKerrSchild::Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]){

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      gcov[i][j] = 0.;
    }
  }

  Real a = bh_spin_;
  Real r = x[IMC1];
  Real r2 = SQR(r);
  Real th = x[IMC2];
  Real sth = sin(th);
  Real cth = cos(th);
  Real cth2 = SQR(cth);
  Real sth2 = SQR(sth);
  Real a2 = SQR(a);

  Real sigma = r2 + a2 * cth2;
  Real delta = r2 - 2 * r + a2;
  Real A = SQR(r2 + a2) - a2 * delta * sth2;
  
  gcov[IMC0][IMC0] = -1. * (1. - 2. * r / sigma);
  gcov[IMC0][IMC1] = 2 * r / sigma;
  gcov[IMC0][IMC3] = -2. * a * r * sth2 / sigma;

  gcov[IMC1][IMC0] = gcov[IMC0][IMC1];
  gcov[IMC1][IMC1] = 1. + 2. * r / sigma;
  gcov[IMC1][IMC3] = -a * sth2 * (1. + 2. * r / sigma);

  gcov[IMC2][IMC2] = sigma;

  gcov[IMC3][IMC0] = gcov[IMC0][IMC3];
  gcov[IMC3][IMC1] = gcov[IMC1][IMC3];
  gcov[IMC3][IMC3] = A * sth2 / sigma;

}

void MCKerrSchild::InverseMetric(Real x[NCOORD], Real gcon[NCOORD][NCOORD]) {
  
  // equations come from Takahasi (2007) Appendix

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      gcon[i][j] = 0.;
    }
  }

  Real a = bh_spin_;
  Real r = x[IMC1];
  Real r2 = SQR(r);
  Real th = x[IMC2];
  Real sth = sin(th);
  Real cth = cos(th);
  Real cth2 = SQR(cth);
  Real sth2 = SQR(sth);
  Real a2 = SQR(a);

  Real sigma = r2 + a2 * cth2;
  Real delta = r2 - 2 * r + a2;

  gcon[IMC0][IMC0] = -(1. + (2. * r / sigma));
  gcon[IMC0][IMC1] = 2. * r / sigma;

  gcon[IMC1][IMC0] = gcon[IMC0][IMC1];
  gcon[IMC1][IMC1] = delta / sigma;
  gcon[IMC1][IMC3] = a / sigma;

  gcon[IMC2][IMC2] = 1. / sigma;
  
  gcon[IMC3][IMC1] = gcon[IMC1][IMC3];
  gcon[IMC3][IMC3] = 1. / (sigma * sth2);

}

void MCKerrSchild::Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]){

  for(int i = 0; i < NCOORD; i++) {
    for(int j = 0; j < NCOORD; j++) {
      for(int k = 0; k < NCOORD; k++) {
	gamma[i][j][k]=0;
      }
    }
  }

  Real a = bh_spin_;
  Real r = x[IMC1];
  Real r2 = SQR(r);
  Real sth = sin(x[IMC2]);
  Real cth = cos(x[IMC2]);
  Real cth2 = SQR(cth);
  Real sth2 = SQR(sth);
  Real s2th = 2.*sth*cth;
  Real c2th = cth2 - sth2;

  Real a2 = SQR(a);
  Real sigma = r2 + a2 * cth2;
  Real sigma2 = SQR(sigma);
  Real delta = r2 - 2 * r + a2;
  Real A = SQR(r2 + a2) - a2 * delta * sth2;

  gamma[IMC0][IMC0][IMC0] = -2. * r / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC0][IMC0][IMC1] = -1. / sigma * (1. + 2. * r / sigma) * (1. - 2. * r2 / sigma);
  gamma[IMC0][IMC0][IMC2] = -a2 * r * s2th / sigma2;
  gamma[IMC0][IMC0][IMC3] = 2. * a * r * sth2 / sigma2 * (1. - 2. * r2 / sigma);

  gamma[IMC0][IMC1][IMC0] = gamma[IMC0][IMC0][IMC1];
  gamma[IMC0][IMC1][IMC1] = -2. / sigma * (1. + r / sigma) * (1. - 2. * r2 / sigma);
  gamma[IMC0][IMC1][IMC2] = -a2 * r * s2th / sigma2;
  gamma[IMC0][IMC1][IMC3] = a * sth2 / sigma * (1. + 2. * r / sigma) * 
    (1. - 2. * r2 / sigma);

  gamma[IMC0][IMC2][IMC0] = gamma[IMC0][IMC0][IMC2];
  gamma[IMC0][IMC2][IMC1] = gamma[IMC0][IMC1][IMC2];
  gamma[IMC0][IMC2][IMC2] = -2. * r2 / sigma;
  gamma[IMC0][IMC2][IMC3] = a2 * a * r / sigma2 * sth2 * s2th;

  gamma[IMC0][IMC3][IMC0] = gamma[IMC0][IMC0][IMC3];
  gamma[IMC0][IMC3][IMC1] = gamma[IMC0][IMC1][IMC3];
  gamma[IMC0][IMC3][IMC2] = gamma[IMC0][IMC2][IMC3];
  gamma[IMC0][IMC3][IMC3] = -2. * r * sth2 / sigma * (r + a2 * sth2 / sigma * 
						      (1. - 2. * r2 / sigma));


  gamma[IMC1][IMC0][IMC0] = -delta / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC1][IMC0][IMC1] = 1. / sigma * (1. - 2. * r2 / sigma) * (1. - delta / sigma);
  gamma[IMC1][IMC0][IMC2] = 0.;
  gamma[IMC1][IMC0][IMC3] = a * delta * sth2 / sigma2 * (1. - 2. * r2 / sigma);

  gamma[IMC1][IMC1][IMC0] = gamma[IMC1][IMC0][IMC1];
  gamma[IMC1][IMC1][IMC1] = 1. / sigma * (1. - 2. * r2 / sigma) * (2. - delta / sigma);
  gamma[IMC1][IMC1][IMC2] = -a2 / (2. * sigma) * s2th;
  gamma[IMC1][IMC1][IMC3] = a / sigma * sth2 * (r - (1. - 2. * r2 / sigma) * 
						(1. - delta / sigma));

  gamma[IMC1][IMC2][IMC0] = gamma[IMC1][IMC0][IMC2];
  gamma[IMC1][IMC2][IMC1] = gamma[IMC1][IMC1][IMC2];
  gamma[IMC1][IMC2][IMC2] = -r * delta / sigma;
  gamma[IMC1][IMC2][IMC3] = 0.;
 
  gamma[IMC1][IMC3][IMC0] = gamma[IMC1][IMC0][IMC3];
  gamma[IMC1][IMC3][IMC1] = gamma[IMC1][IMC1][IMC3];
  gamma[IMC1][IMC3][IMC2] = gamma[IMC1][IMC2][IMC3];
  gamma[IMC1][IMC3][IMC3] = -delta / sigma * sth2 * (r + a2 * sth2 / sigma * 
						     (1. - 2. * r2 / sigma));
  

  gamma[IMC2][IMC0][IMC0] = -a2 * r * s2th / (sigma2 * sigma);
  gamma[IMC2][IMC0][IMC1] = -a2 * r * s2th / (sigma2 * sigma);
  gamma[IMC2][IMC0][IMC2] = 0.;
  gamma[IMC2][IMC0][IMC3] = a * r * (r2 + a2) * s2th / (sigma2 * sigma);

  gamma[IMC2][IMC1][IMC0] = gamma[IMC2][IMC0][IMC1];
  gamma[IMC2][IMC1][IMC1] = -a2 * r * s2th / (sigma2 * sigma);
  gamma[IMC2][IMC1][IMC2] = r / sigma;
  /*gamma[IMC2][IMC1][IMC3] = a / (2. * sigma) * (1. + 2. * r * (r2 + a2) / sigma2) * s2th;*/
  gamma[IMC2][IMC1][IMC3] = ((a * cth * sth) / (sigma2 * sigma)) * 
    (r2 * r * (r + 2.) + 2. * a2 * r * (r + 1.) * cth2 + a2 * a2 * cth2 * cth2 
    + 2. * a2 * r * sth2); // from Shane's notebook -- not equal to Takahashi+07 

  gamma[IMC2][IMC2][IMC0] = gamma[IMC2][IMC0][IMC2];
  gamma[IMC2][IMC2][IMC1] = gamma[IMC2][IMC1][IMC2];
  gamma[IMC2][IMC2][IMC2] = -a2 * s2th / (2. * sigma);
  gamma[IMC2][IMC2][IMC3] = 0.;

  gamma[IMC2][IMC3][IMC0] = gamma[IMC2][IMC0][IMC3];
  gamma[IMC2][IMC3][IMC1] = gamma[IMC2][IMC1][IMC3];
  gamma[IMC2][IMC3][IMC2] = gamma[IMC2][IMC2][IMC3];
  /*gamma[IMC2][IMC3][IMC3] = -s2th / (2. * sigma) * (delta + 2. * r * 
    SQR((r2 + a2) / sigma));*/
  gamma[IMC2][IMC3][IMC3] = -(cth * sth / (sigma2 * sigma)) * 
    (a2 * a2 * a2 * cth2 * cth2 * cth2 +
     cth2 * cth2 * (3. * a2 * a2 * r2 + a2 * a2 * a2 * sth2) + 
     cth2 * (3. * a2 * r2 * r2 + 2. * a2 * a2 * r2 * sth2) +
     r * (r2 * r2 * r + a2 * r2 * (r + 4.) * sth2 + 2. * a2 * a2 * sth2 * sth2 +
     a2 * a2 * s2th * s2th)); // frome Shane's notebook


  gamma[IMC3][IMC0][IMC0] = -a / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC3][IMC0][IMC1] = -a / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC3][IMC0][IMC2] = -2. * a * r / sigma2 * cth / sth;
  gamma[IMC3][IMC0][IMC3] = a2 * sth2 / sigma2 * (1. - 2. * r2 / sigma);

  gamma[IMC3][IMC1][IMC0] = gamma[IMC3][IMC0][IMC1];
  gamma[IMC3][IMC1][IMC1] = -a / sigma2 * (1. - 2. * r2 / sigma);
  gamma[IMC3][IMC1][IMC2] = -a / sigma * (1. + 2. * r / sigma) * cth / sth;
  gamma[IMC3][IMC1][IMC3] = 1. / sigma * (r + a2 * sth2 / sigma * 
					  (1. - 2. * r2 / sigma));

  gamma[IMC3][IMC2][IMC0] = gamma[IMC3][IMC0][IMC2];
  gamma[IMC3][IMC2][IMC1] = gamma[IMC3][IMC1][IMC2];
  gamma[IMC3][IMC2][IMC2] = -a * r / sigma;
  //gamma[IMC3][IMC2][IMC3] = (1. + 2. * r / sigma * ((r2 + a2) / sigma - 1.)) * cth / sth;
  gamma[IMC3][IMC2][IMC3] = ((1. / 4.) * SQR(a2 + 2. * r2 + a2 * c2th) * cth / sth + 
    a2 * r * s2th) / (sigma2); // from Shane's notebook

  gamma[IMC3][IMC3][IMC0] = gamma[IMC3][IMC0][IMC3];
  gamma[IMC3][IMC3][IMC1] = gamma[IMC3][IMC1][IMC3];
  gamma[IMC3][IMC3][IMC2] = gamma[IMC3][IMC2][IMC3];
  gamma[IMC3][IMC3][IMC3] = -a / sigma * sth2 * (r + a2 * sth2 / sigma * 
						 (1. - 2. * r2 / sigma));

}


// constructor
MCBoyerLindquist::MCBoyerLindquist(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

// constructor
MCBoyerLindquist::MCBoyerLindquist(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

// destructor
MCBoyerLindquist::~MCBoyerLindquist() {

}

void MCBoyerLindquist::Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]) { 

  // equation for the Metric comes from the inside cover of Hartle
  
  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      gcov[i][j] = 0.;
    }
  }

  Real a = bh_spin_;
  Real m = bh_mass_;
  // a in input file is dimensionless spin
  a *= m;
  Real r = x[IMC1];
  Real sth = sin(x[IMC2]);
  Real cth = cos(x[IMC2]);
  Real cth2 = SQR(cth);
  Real sth2 = SQR(sth);

  Real r2 = SQR(r);
  Real a2 = SQR(a);
  Real rho2 = r2 + a2*cth2;
  Real delta = r2 - 2. * m* r + a2;
 
  gcov[IMC0][IMC0] = (-1. + 2.*m*r/rho2);
  gcov[IMC1][IMC1] = rho2/delta;
  gcov[IMC2][IMC2] = rho2;
  gcov[IMC3][IMC3] = (r2 + a2 + 2.*m*r*a2*sth2/rho2)*sth2;

  //gcov[IMC0][IMC3] = -2.*m*a*r*sth2/rho2; SWD <- old, think m is superfluous with a=m*a_*
  gcov[IMC0][IMC3] = -2.*a*r*sth2/rho2;
  gcov[IMC3][IMC0] = gcov[IMC0][IMC3];


}

void MCBoyerLindquist::InverseMetric(Real x[NCOORD], Real gcon[NCOORD][NCOORD]) {

  // Equation comes from ColinsCosmos.com/wiki/boyer-lindquist-coordinates, which
  // sites Frolov & Novikov Section D.1 (but I don't have access to this book)


  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      gcon[i][j] = 0.;
    }
  }


  Real a = bh_spin_;
  Real m = bh_mass_;
  // a in input file is dimensionless spin
  a *= m;
  Real r = x[IMC1];
  Real sth = sin(x[IMC2]);
  Real cth = cos(x[IMC2]);
  Real cth2 = SQR(cth);
  Real sth2 = SQR(sth);

  Real r2 = SQR(r);
  Real a2 = SQR(a);
  Real rho2 = r2 + a2 * cth2;
  Real delta = r2 - 2. * m * r + a2;

  gcon[IMC0][IMC0] = -1. / delta * (r2 + a2 + 2. * r * m * a2 * sth2 / rho2);
  gcon[IMC1][IMC1] = delta / rho2;
  gcon[IMC2][IMC2] = 1. / rho2;
  gcon[IMC3][IMC3] = (delta - a2 * sth2) / (rho2 * delta * sth2);

  gcon[IMC0][IMC3] = -2. * m * r * a / (rho2 * delta);
  gcon[IMC3][IMC0] = gcon[IMC0][IMC3];

}


void MCBoyerLindquist::Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]) {

  // equations for the connection coefficients come from Frutos-Alfaro et al. (2012)

  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      for (int k = 0; k < NCOORD; k++) {
	gamma[i][j][k] = 0.;
      }
    }
  }

  // SWD: Clean this one up
  Real a = bh_spin_;
  Real m = bh_mass_;
  // a has units of mass in Frutos-Alfaro
  a *= m; 
  Real r = x[IMC1];
  Real j = a*m;

  Real sth = sin(x[IMC2]);
  Real cth = cos(x[IMC2]);
  Real cth2 = SQR(cth);
  Real sth2 = 1. - cth2;

  Real r2 = SQR(r);
  Real a2 = SQR(a);
  Real rho2 = r2 + a2*cth2;
  Real rs = 2.*m; 
  Real delta = r2 - rs*r + a2;

  Real rho4 = SQR(rho2);
  Real rho6 = rho4*rho2;

  // Real A = SQR(r2+a2) - a2*delta*sth2;

  gamma[IMC0][IMC0][IMC1] = rs/(2.*rho4*delta)*(r2+a2)*(2.*r2-rho2);
  gamma[IMC0][IMC0][IMC2] = -2.*a*j*r/rho4*sth*cth;
  
  gamma[IMC0][IMC1][IMC0] = gamma[IMC0][IMC0][IMC1];
  gamma[IMC0][IMC1][IMC3] = -j*sth2/(rho4*delta)*(rho2*(r2-a2)+2.*r2*(r2+a2));

  gamma[IMC0][IMC2][IMC3] = 2.*a2*j*r/rho4*cth*sth2*sth;
  gamma[IMC0][IMC3][IMC2] = gamma[IMC0][IMC2][IMC3];

  gamma[IMC1][IMC0][IMC0] = rs*delta/(2.*rho6)*(2.*r2-rho2);
  gamma[IMC1][IMC0][IMC3] = -j*delta/rho6*(2.*r2-rho2)*sth2;
  
  gamma[IMC1][IMC1][IMC1] = 1./(rho2*delta)*(rho2*(rs/2.-r)+r*delta);
  gamma[IMC1][IMC1][IMC2] = -a2/rho2*sth*cth;
  
  gamma[IMC1][IMC2][IMC1] = gamma[IMC1][IMC1][IMC2];
  gamma[IMC1][IMC2][IMC2] = -r*delta/rho2;

  gamma[IMC1][IMC3][IMC0] = gamma[IMC1][IMC0][IMC3];
  gamma[IMC1][IMC3][IMC3] = -delta*sth2/rho6*(r*rho4-a*j*(2.*r2-rho2)*sth2);
  
  gamma[IMC2][IMC0][IMC0] = -2.*a*j*r/rho6*sth*cth;
  gamma[IMC2][IMC0][IMC3] = 2.*j*r/rho6*(r2+a2)*sth*cth;

  gamma[IMC2][IMC1][IMC1] = a2/(rho2*delta)*sth*cth;
  gamma[IMC2][IMC1][IMC2] = r/rho2;
 
  gamma[IMC2][IMC2][IMC1] = gamma[IMC2][IMC1][IMC2];
  gamma[IMC2][IMC2][IMC2] = gamma[IMC1][IMC1][IMC2];

  gamma[IMC2][IMC3][IMC0] = gamma[IMC2][IMC0][IMC3];
  gamma[IMC2][IMC3][IMC3] = -sth*cth/rho6*(rho4*delta+rs*r*SQR(r2+a2));
  //gamma[IMC2][IMC2][IMC3] = -sth*cth/rho6*(A*rho2+(r2+a2)*a2*r*rs*sth2);

  gamma[IMC3][IMC0][IMC1] = j/(rho4*delta)*(2.*r2-rho2);
  gamma[IMC3][IMC0][IMC2] = -2.*j*r*cth/(rho4*sth);

  gamma[IMC3][IMC1][IMC3] = gamma[IMC3][IMC0][IMC1];
  gamma[IMC3][IMC1][IMC3] = 1./(rho4*delta)*(r*rho2*(rho2-rs*r)-a*j*sth2*(2.*r2-rho2));

  gamma[IMC3][IMC2][IMC0] = gamma[IMC3][IMC0][IMC2];
  gamma[IMC3][IMC2][IMC3] = cth/(rho4*sth)*(rho4+2.*a*j*r*sth2);

  gamma[IMC3][IMC3][IMC1] = gamma[IMC3][IMC1][IMC3];
  gamma[IMC3][IMC3][IMC2] = gamma[IMC3][IMC2][IMC3];
}

// constructor
MCMinkowski::MCMinkowski(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

// constructor
MCMinkowski::MCMinkowski(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

// destructor
MCMinkowski::~MCMinkowski() {

}





