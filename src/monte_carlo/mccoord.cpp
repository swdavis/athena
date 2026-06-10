//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//! \brief implementation of functions in class MCCoord

// SWD: General notes:
// * This whole class maybe should be reworked
// * inverse metric for spherical polar and cylindrical
// * remove zeroing of metric,connection at top?
// * remove extraneous variable definitions
// * conversion to cgs units for black hole metrics?

// Athena++ headers
#include "../athena.hpp"
#include "mccoord.hpp"

//----------------------------------------------------------------------------------------
//! MCCoord base class constructor, builds MCCoord from Coord and MonteCarloBlock

MCCoord::MCCoord(Coordinates *pcoord, MonteCarloBlock *pmcb) {

  int nx1 = pcoord->x1f.GetDim1();
  int nx2 = pcoord->x2f.GetDim1();
  int nx3 = pcoord->x3f.GetDim1();

  x1f.NewAthenaArray(nx1);
  x2f.NewAthenaArray(nx2);
  x3f.NewAthenaArray(nx3);

  x1f = pcoord->x1f;
  //for (int i = 0; i < nx1; i++)
  //  x1f(i) *= pmcb->l_cgs;
  x2f = pcoord->x2f;
  x3f = pcoord->x3f;
  //if ( (COORDINATE_SYSTEM == "cartesian") || (COORDINATE_SYSTEM == "minkowski")
  //     || (COORDINATE_SYSTEM == "gr_user") ) {
  //  for (int i = 0; i < nx2; i++)
  //    x2f(i) *= pmcb->l_cgs;
  //  for (int i = 0; i < nx3; i++)
  //    x3f(i) *= pmcb->l_cgs;
  //}
  // Needed for black hole coordinates
  if (GENERAL_RELATIVITY) {
    bh_mass_ = pcoord->GetMass();
    bh_spin_ = pcoord->GetSpin();
  } else {
    // initialize to 0 for flat spacetimes
    bh_mass_ = 0.;
    bh_spin_ = 0.;
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
        // Volume in cgs units
        vol(k,j,i) = pcoord->GetCellVolume(k,j,i) * pow(pmcb->l_cgs,3);
        if (std::isnan(vol(k,j,i)) && (COORDINATE_SYSTEM == "gr_user")) {
          // at the singularity set volume to zero
         vol(k,j,i) = 0.;
        } 
      }}}
  computedmin = pmcb->computedmin;
  if (computedmin) {
    dmin.NewAthenaArray(ncells3,ncells2,ncells1);
    Real dw1,dw2,dw3;
    for (int k=pmcb->ks; k<=pmcb->ke; ++k) {
      for (int j=pmcb->js; j<=pmcb->je; ++j) {
        for (int i=pmcb->is; i<=pmcb->ie; ++i) {
          dw3 = pcoord->dx3f(k);
          dw2 = pcoord->dx2f(j);
          dw1 = pcoord->dx1f(i);
          //dw1 *= pmcb->l_cgs;
          //if ( (COORDINATE_SYSTEM == "cartesian") ||
          //     (COORDINATE_SYSTEM == "minkowski")
          //     || (COORDINATE_SYSTEM == "gr_user") ) {
          //  dw2 *= pmcb->l_cgs;
          //  dw3 *= pmcb->l_cgs;
          //}
          Real dmin0 = std::min(dw1,dw2);
          dmin(k,j,i) = std::min(dmin0,dw3);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! MonteCarlo constructor for processes without own MeshBlock

MCCoord::MCCoord(int ncells1, int ncells2, int ncells3, bool cdmin) {

  x1f.NewAthenaArray(ncells1+1);
  x2f.NewAthenaArray(ncells2+1);
  x3f.NewAthenaArray(ncells3+1);

  vol.NewAthenaArray(ncells3,ncells2,ncells1);
  computedmin = cdmin;
  if (cdmin)
    dmin.NewAthenaArray(ncells3,ncells2,ncells1);
}

//----------------------------------------------------------------------------------------
//! destructor

MCCoord::~MCCoord() {

  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();
  vol.DeleteAthenaArray();
  if (computedmin)
    dmin.DeleteAthenaArray();
}

//----------------------------------------------------------------------------------------
//! \fn void MCCoord::Metric(Real x[4], Real gcov[4][4])
//! \brief compute metric in flat spacetime cartesian

void MCCoord::Metric(Real x[4], Real gcov[4][4]) {

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
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

//----------------------------------------------------------------------------------------
//! \fn void MCCoord::MetricDerivative(Real x[4], Real dgcov[4][4][4])
//! \brief compute metric derivative in flat spacetime cartesian

void MCCoord::MetricDerivative(Real x[4], Real dgcov[4][4][4]) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 4; k++) {
        dgcov[i][j][k]=0;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MCCoord::InverseMetric(Real x[4], Real gcon[4][4])
//! \brief compute inverse metric in flat spacetime cartesian

void MCCoord::InverseMetric(Real x[4], Real gcon[4][4]) {

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
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

//----------------------------------------------------------------------------------------
//! \fn void MCCoord::InverseMetricDerivative(Real x[4], Real dgcon[4][4][4])
//! \brief compute derivative of inverse metric in flat spacetime cartesian

void MCCoord::InverseMetricDerivative(Real x[4], Real dgcon[4][4][4]) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 4; k++) {
        dgcon[i][j][k]=0;
      }
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MCCoord::Connect(Real x[4], Real gamma[4][4][4])
//! \brief compute connection in flat spacetime cartesian

void MCCoord::Connect(Real x[4], Real gamma[4][4][4]) {

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        gamma[i][j][k] = 0.;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MCCoord::Tetrad(Real x[4], Real tetrad[4][4])
//! \brief compute a diagonal tetrad

void MCCoord::Tetrad(Real x[4], Real tetrad[4][4]) {

  for (int l=0; l<4; l++) {
    for (int m=0; m<4; m++) {
      tetrad[l][m] = 0.;
    }
    tetrad[l][l] = 1.;
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MCCoord::InverseTetrad(Real x[4], Real invtet[4][4])
//! \brief compute a diagonal tetrad

void MCCoord::InverseTetrad(Real x[4], Real invtet[4][4]) {

  for (int l=0; l<4; l++) {
    for (int m=0; m<4; m++) {
      invtet[l][m] = 0.;
    }
    invtet[l][l] = 1.;
  }

}

//----------------------------------------------------------------------------------------
//! MCCartesian constructor, built from Coord and MonteCarloBlock

MCCartesian::MCCartesian(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

//----------------------------------------------------------------------------------------
//! MCCartesian constructor for processes without own MeshBlock

MCCartesian::MCCartesian(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

//----------------------------------------------------------------------------------------
//! destructor

MCCartesian::~MCCartesian() {

}

//----------------------------------------------------------------------------------------
//! MCSphericalPolar constructor, built from Coord and MonteCarloBlock

MCSphericalPolar::MCSphericalPolar(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

//----------------------------------------------------------------------------------------
//! MCSphericalPolar constructor for processes without own MeshBlock

MCSphericalPolar::MCSphericalPolar(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

//----------------------------------------------------------------------------------------
//! destructor
MCSphericalPolar::~MCSphericalPolar() {

}

//----------------------------------------------------------------------------------------
//! \fn void MCSphericalPolar::Metric(Real x[4], Real gcov[4][4])
//! \brief compute metric in flat spacetime spherical polar

void MCSphericalPolar::Metric(Real x[4], Real gcov[4][4]) {

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      gcov[i][j] = 0;
    }
  }

  Real r = x[IMC1];
  Real sth = sin(x[IMC2]);
  gcov[IMC0][IMC0] = -1.;
  gcov[IMC1][IMC1] = 1.;
  gcov[IMC2][IMC2] = r * r;
  gcov[IMC3][IMC3] = r * r * sth * sth;

}

//----------------------------------------------------------------------------------------
//! \fn void MCSphericalPolar::InverseMetric(Real x[4], Real gcon[4][4])
//! \brief compute metric in flat spacetime spherical polar

void MCSphericalPolar::InverseMetric(Real x[4], Real gcon[4][4]) {

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      gcon[i][j] = 0;
    }
  }
  Real r = x[IMC1];
  Real sth = sin(x[IMC2]);
  gcon[IMC0][IMC0] = -1;
  gcon[IMC1][IMC1] = 1;
  gcon[IMC2][IMC2] = 1. / (r * r);
  gcon[IMC3][IMC3] = 1. / (r * r * sth *sth);

}

//----------------------------------------------------------------------------------------
//! \fn void MCSphericalPolar::Connect(Real x[4], Real gamma[4][4][4])
//! \brief compute connection in flat spacetime spherical polar

void MCSphericalPolar::Connect(Real x[4], Real gamma[4][4][4]) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 4; k++) {
        gamma[i][j][k]=0;
      }
    }
  }
  Real r = x[IMC1];
  Real sth = sin(x[IMC2]);
  Real cth = cos(x[IMC2]);

  gamma[IMC1][IMC2][IMC2] = -r;
  gamma[IMC1][IMC3][IMC3] = -r * sth * sth;
  gamma[IMC2][IMC1][IMC2] = 1. / r;
  gamma[IMC2][IMC2][IMC1] = gamma[IMC2][IMC1][IMC2];
  gamma[IMC2][IMC3][IMC3] = -sth*cth;
  gamma[IMC3][IMC1][IMC3] = 1. / r;
  gamma[IMC3][IMC2][IMC3] = cth / sth;
  gamma[IMC3][IMC3][IMC1] = gamma[IMC3][IMC1][IMC3];
  gamma[IMC3][IMC3][IMC2] = gamma[IMC3][IMC2][IMC3];

}

//----------------------------------------------------------------------------------------
//! \fn void MCSphericalPolar::InverseMetricDerivative(Real x[4], Real dgcon[4][4][4])
//! \brief compute a diagonal tetrad

void MCSphericalPolar::InverseMetricDerivative(Real x[4], Real dgcon[4][4][4]) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 4; k++) {
        dgcon[i][j][k]=0;
      }
    }
  }

  Real r = x[IMC1];
  Real sth = sin(x[IMC2]);
  Real cth = cos(x[IMC2]);
  dgcon[IMC1][IMC2][IMC2] = -2./(r*r*r);
  dgcon[IMC1][IMC3][IMC3] = -2./(r*r*r)/(sth*sth);
  dgcon[IMC2][IMC3][IMC3] = -2.*cth/(r*r)/(sth*sth*sth);

}

//----------------------------------------------------------------------------------------
//! \fn void MCSphericalPolar::Tetrad(Real x[4], Real tetrad[4][4])
//! \brief compute a diagonal tetrad

void MCSphericalPolar::Tetrad(Real x[4], Real tetrad[4][4]) {

  for (int l=0; l<4; l++) {
    for (int m=0; m<4; m++) {
      tetrad[l][m] = 0.;
    }
    tetrad[l][l] = 1.;
  }
  Real r = x[IMC1];
  Real sth = sin(x[IMC2]);
  tetrad[IMC2][IMC2] = 1. / r;
  tetrad[IMC3][IMC3] = 1. / (r*sth);
}

//----------------------------------------------------------------------------------------
//! \fn void MCSphericalPolar::InverseTetrad(Real x[4], Real invtet[4][4])
//! \brief compute a diagonal tetrad

void MCSphericalPolar::InverseTetrad(Real x[4], Real invtet[4][4]) {

  for (int l=0; l<4; l++) {
    for (int m=0; m<4; m++) {
      invtet[l][m] = 0.;
    }
    invtet[l][l] = 1.;
  }
  Real r = x[IMC1];
  Real sth = sin(x[IMC2]);
  invtet[IMC2][IMC2] = r;
  invtet[IMC3][IMC3] = r * sth;
}
//----------------------------------------------------------------------------------------
//! MCCylindrical constructor, built from Coord and MonteCarloBlock

MCCylindrical::MCCylindrical(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {
}

//----------------------------------------------------------------------------------------
//! MCCylindrical constructor for processes without own MeshBlock

MCCylindrical::MCCylindrical(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {
}

//----------------------------------------------------------------------------------------
//! destructor
MCCylindrical::~MCCylindrical() {

}

//----------------------------------------------------------------------------------------
//! \fn void MCCylindrical::Metric(Real x[4], Real gcov[4][4])
//! \brief compute metric in flat spacetime cylindrical

void MCCylindrical::Metric(Real x[4], Real gcov[4][4]) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      gcov[i][j] = 0;
    }
  }
  gcov[IMC0][IMC0] = -1;
  gcov[IMC1][IMC1] = 1;
  gcov[IMC2][IMC2] = x[IMC1] * x[IMC1];
  gcov[IMC3][IMC3] = 1;

}

//----------------------------------------------------------------------------------------
//! \fn void MCCylindrical::Connect(Real x[4], Real gamma[4][4][4])
//! \brief compute connection in flat spacetime cylindrical

void MCCylindrical::Connect(Real x[4], Real gamma[4][4][4]) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 4; k++) {
        gamma[i][j][k]=0;
      }
    }
  }

  gamma[IMC1][IMC2][IMC2] = -x[IMC1];
  gamma[IMC2][IMC1][IMC2] = 1./x[IMC1];
  gamma[IMC2][IMC2][IMC1] = 1./x[IMC1];
}

//----------------------------------------------------------------------------------------
//! MCKerrSchild constructor, built from Coord and MonteCarloBlock

MCKerrSchild::MCKerrSchild(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

//----------------------------------------------------------------------------------------
//! MCKerrSchild constructor for processes without own MeshBlock

MCKerrSchild::MCKerrSchild(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

//----------------------------------------------------------------------------------------
//! destructor
MCKerrSchild::~MCKerrSchild() {

}

//----------------------------------------------------------------------------------------
//! \fn void MCKerrSchild::Metric(Real x[4], Real gcov[4][4])
//! \brief compute metric in spherical polar Kerr-Schild

void MCKerrSchild::Metric(Real x[4], Real gcov[4][4]) {

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
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

//----------------------------------------------------------------------------------------
//! \fn void MCKerrSchild::InverseMetric(Real x[4], Real gcon[4][4])
//! \brief compute inverse metric in spherical polar Kerr-Schild

void MCKerrSchild::InverseMetric(Real x[4], Real gcon[4][4]) {

  // equations come from Takahasi (2007) Appendix
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
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

//----------------------------------------------------------------------------------------
//! \fn void MCKerrSchild::Connect(Real x[4], Real gamma[4][4][4])
//! \brief compute connection in spherical polar Kerr-Schild

void MCKerrSchild::Connect(Real x[4], Real gamma[4][4][4]) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 4; k++) {
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

  // from Shane's notebook -- not equal to Takahashi+07 (SWD: ?)
  gamma[IMC2][IMC1][IMC3] = ((a * cth * sth) / (sigma2 * sigma)) *
    (r2 * r * (r + 2.) + 2. * a2 * r * (r + 1.) * cth2 + a2 * a2 * cth2 * cth2
    + 2. * a2 * r * sth2);

  gamma[IMC2][IMC2][IMC0] = gamma[IMC2][IMC0][IMC2];
  gamma[IMC2][IMC2][IMC1] = gamma[IMC2][IMC1][IMC2];
  gamma[IMC2][IMC2][IMC2] = -a2 * s2th / (2. * sigma);
  gamma[IMC2][IMC2][IMC3] = 0.;

  gamma[IMC2][IMC3][IMC0] = gamma[IMC2][IMC0][IMC3];
  gamma[IMC2][IMC3][IMC1] = gamma[IMC2][IMC1][IMC3];
  gamma[IMC2][IMC3][IMC2] = gamma[IMC2][IMC2][IMC3];
  /*gamma[IMC2][IMC3][IMC3] = -s2th / (2. * sigma) * (delta + 2. * r *
    SQR((r2 + a2) / sigma));*/
  // from Shane's notebook -- not equal to Takahashi+07 (SWD: ?)
  gamma[IMC2][IMC3][IMC3] = -(cth * sth / (sigma2 * sigma)) *
    (a2 * a2 * a2 * cth2 * cth2 * cth2 +
     cth2 * cth2 * (3. * a2 * a2 * r2 + a2 * a2 * a2 * sth2) +
     cth2 * (3. * a2 * r2 * r2 + 2. * a2 * a2 * r2 * sth2) +
     r * (r2 * r2 * r + a2 * r2 * (r + 4.) * sth2 + 2. * a2 * a2 * sth2 * sth2 +
     a2 * a2 * s2th * s2th));

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
  //gamma[IMC3][IMC2][IMC3] = (1. + 2. * r / sigma * ((r2 + a2) / sigma - 1.)) *
  // cth / sth;
  // from Shane's notebook -- not equal to Takahashi+07 (SWD: ?)
  gamma[IMC3][IMC2][IMC3] = ((1. / 4.) * SQR(a2 + 2. * r2 + a2 * c2th) * cth / sth +
    a2 * r * s2th) / (sigma2);

  gamma[IMC3][IMC3][IMC0] = gamma[IMC3][IMC0][IMC3];
  gamma[IMC3][IMC3][IMC1] = gamma[IMC3][IMC1][IMC3];
  gamma[IMC3][IMC3][IMC2] = gamma[IMC3][IMC2][IMC3];
  gamma[IMC3][IMC3][IMC3] = -a / sigma * sth2 * (r + a2 * sth2 / sigma *
                            (1. - 2. * r2 / sigma));

}

//----------------------------------------------------------------------------------------
//! \fn void MCKerrSchild::InverseMetricDerivative(Real x[4], Real dgcon[4][4][4])
//! \brief compute derivative of the inverse metric in spherical polar Kerr-Schild

void MCKerrSchild::InverseMetricDerivative(Real x[4], Real dgcon[4][4][4]) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 4; k++) {
        dgcon[i][j][k]=0;
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
  Real alts = r2 - a2 * cth2;

  dgcon[IMC1][IMC0][IMC0] = 2. * alts / sigma2;
  dgcon[IMC1][IMC0][IMC1] = -dgcon[IMC1][IMC0][IMC0];
  dgcon[IMC1][IMC1][IMC0] = dgcon[IMC1][IMC0][IMC1];
  dgcon[IMC1][IMC1][IMC1] = 2.*(alts-a2*r*sth2) / sigma2;
  dgcon[IMC1][IMC1][IMC3] = -2*a*r / sigma2;
  dgcon[IMC1][IMC2][IMC2] = -2*r / sigma2;
  dgcon[IMC1][IMC3][IMC1] = dgcon[IMC1][IMC1][IMC3];
  dgcon[IMC1][IMC3][IMC3] = -2*r / sth2 / sigma2;

  dgcon[IMC2][IMC0][IMC0] = -2.*a2*r*s2th / sigma2;
  dgcon[IMC2][IMC0][IMC1] = -dgcon[IMC2][IMC0][IMC0];
  dgcon[IMC2][IMC1][IMC0] = dgcon[IMC2][IMC0][IMC1];
  dgcon[IMC2][IMC1][IMC1] = a2*(a2+r*(r-2.))*s2th / sigma2;
  dgcon[IMC2][IMC1][IMC3] = a*a2*s2th / sigma2;
  dgcon[IMC2][IMC2][IMC2] = a2*s2th / sigma2;
  dgcon[IMC2][IMC3][IMC1] = dgcon[IMC2][IMC1][IMC3];
  dgcon[IMC2][IMC3][IMC3] = -2*(r2+a2*c2th)*cth/sth/sth2 / sigma2;
}

//----------------------------------------------------------------------------------------
//! MCKerrSchild constructor, built from Coord and MonteCarloBlock

MCKerrSchildCartesian::MCKerrSchildCartesian(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {

}

//----------------------------------------------------------------------------------------
//! MCKerrSchild constructor for processes without own MeshBlock

MCKerrSchildCartesian::MCKerrSchildCartesian(int ncells1, int ncells2, int ncells3,
                                             bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {

}

//----------------------------------------------------------------------------------------
//! destructor
MCKerrSchildCartesian::~MCKerrSchildCartesian() {

}

//----------------------------------------------------------------------------------------
//! \fn void MCKerrSchildCartesian::Metric(Real x[4], Real gcov[4][4])
//! \brief compute metric in cartesian Kerr-Schild

void MCKerrSchildCartesian::Metric(Real x[4], Real gcov[4][4]) {

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      gcov[i][j] = 0.;
    }
  }

  // Calculate useful quantities
  Real a = bh_spin_;
  Real m = bh_mass_;
  Real a2 = a * a;
  Real rr2 = x[IMC1] * x[IMC1] + x[IMC2] * x[IMC2] + x[IMC3] * x[IMC3];
  Real r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * a * x[IMC3]));
  Real r = std::sqrt(r2);
  Real f = 2.0 * m * r2 * r / (r2 * r2 + a2 * x[IMC3] * x[IMC3]);

  // Calculate null vector
  Real l_0 = 1.0;
  Real l_1 = (r * x[IMC1] + a * x[IMC2]) / (r2 + a2);
  Real l_2 = (r * x[IMC2] - a * x[IMC1]) / (r2 + a2);
  Real l_3 = x[IMC3] / r;

  // Calculate metric components
  gcov[IMC0][IMC0] = f * l_0 * l_0 - 1.0;
  gcov[IMC0][IMC1] = f * l_0 * l_1;
  gcov[IMC0][IMC2] = f * l_0 * l_2;
  gcov[IMC0][IMC3] = f * l_0 * l_3;

  gcov[IMC1][IMC0] = f * l_1 * l_0;
  gcov[IMC1][IMC1] = f * l_1 * l_1 + 1.0;
  gcov[IMC1][IMC2] = f * l_1 * l_2;
  gcov[IMC1][IMC3] = f * l_1 * l_3;

  gcov[IMC2][IMC0] = f * l_2 * l_0;
  gcov[IMC2][IMC1] = f * l_2 * l_1;
  gcov[IMC2][IMC2] = f * l_2 * l_2 + 1.0;
  gcov[IMC2][IMC3] = f * l_2 * l_3;

  gcov[IMC3][IMC0] = f * l_3 * l_0;
  gcov[IMC3][IMC1] = f * l_3 * l_1;
  gcov[IMC3][IMC2] = f * l_3 * l_2;
  gcov[IMC3][IMC3] = f * l_3 * l_3 + 1.0;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MCKerrSchildCartesian::InverseMetric(Real x[4], Real gcov[4][4])
//! \brief compute metric in cartesian Kerr-Schild

void MCKerrSchildCartesian::InverseMetric(Real x[4], Real gcon[4][4]) {

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      gcon[i][j] = 0.;
    }
  }

  // Calculate useful quantities
  Real a = bh_spin_;
  Real m = bh_mass_;
  Real a2 = a * a;
  Real rr2 = x[IMC1] * x[IMC1] + x[IMC2] * x[IMC2] + x[IMC3] * x[IMC3];
  Real r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * a * x[IMC3]));
  Real r = std::sqrt(r2);
  Real f = 2.0 * m * r2 * r / (r2 * r2 + a2 * x[IMC3] * x[IMC3]);

  // Calculate null vector
  Real l0 = -1.0;
  Real l1 = (r * x[IMC1] + a * x[IMC2]) / (r2 + a2);
  Real l2 = (r * x[IMC2] - a * x[IMC1]) / (r2 + a2);
  Real l3 = x[IMC3] / r;

  // Calculate metric components
  gcon[IMC0][IMC0] = -f * l0 * l0 - 1.0;
  gcon[IMC0][IMC1] = -f * l0 * l1;
  gcon[IMC0][IMC2] = -f * l0 * l2;
  gcon[IMC0][IMC3] = -f * l0 * l3;

  gcon[IMC1][IMC0] = -f * l1 * l0;
  gcon[IMC1][IMC1] = -f * l1 * l1 + 1.0;
  gcon[IMC1][IMC2] = -f * l1 * l2;
  gcon[IMC1][IMC3] = -f * l1 * l3;

  gcon[IMC2][IMC0] = -f * l2 * l0;
  gcon[IMC2][IMC1] = -f * l2 * l1;
  gcon[IMC2][IMC2] = -f * l2 * l2 + 1.0;
  gcon[IMC2][IMC3] = -f * l2 * l3;

  gcon[IMC3][IMC0] = -f * l3 * l0;
  gcon[IMC3][IMC1] = -f * l3 * l1;
  gcon[IMC3][IMC2] = -f * l3 * l2;
  gcon[IMC3][IMC3] = -f * l3 * l3 + 1.0;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MCKerrSchildCartesian::InverseMetricDerivative(Real x[4],
//           Real dgcon[4][4][4])
//! \brief compute derivative of inverse metric in cartesian Kerr-Schild

void MCKerrSchildCartesian::InverseMetricDerivative(Real x[4], Real dgcon[4][4][4]) {

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 4; k++) {
        dgcon[i][j][k]=0;
      }
    }
  }

  // Calculate useful quantities
  Real a = bh_spin_;
  Real m = bh_mass_;
  Real a2 = a * a;
  Real rr2 = x[IMC1] * x[IMC1] + x[IMC2] * x[IMC2] + x[IMC3] * x[IMC3];
  Real r2 = 0.5 * (rr2 - a2 + std::hypot(rr2 - a2, 2.0 * a * x[IMC3]));
  Real r = std::sqrt(r2);
  Real f = 2.0 * m * r2 * r / (r2 * r2 + a2 * x[IMC3] * x[IMC3]);

  // Calculate null vector
  Real l0 = -1.0;
  Real l1 = (r * x[IMC1] + a * x[IMC2]) / (r2 + a2);
  Real l2 = (r * x[IMC2] - a * x[IMC1]) / (r2 + a2);
  Real l3 = x[IMC3] / r;

  // Calculate scalar derivatives
  Real dr_dx = r * x[IMC1] / (2.0 * r2 - rr2 + a2);
  Real dr_dy = r * x[IMC2] / (2.0 * r2 - rr2 + a2);
  Real dr_dz = (r * x[IMC3] + a2 * x[IMC3] / r) / (2.0 * r2 - rr2 + a2);
  Real df_dx = -(r2 * r2 - 3.0 * a2 * x[IMC3] * x[IMC3]) * dr_dx / (r * (r2 * r2 + a2 * x[IMC3] * x[IMC3])) * f;
  Real df_dy = -(r2 * r2 - 3.0 * a2 * x[IMC3] * x[IMC3]) * dr_dy / (r * (r2 * r2 + a2 * x[IMC3] * x[IMC3])) * f;
  Real df_dz =
      -((r2 * r2 - 3.0 * a2 * x[IMC3] * x[IMC3]) * dr_dz + 2.0 * a2 * r * x[IMC3]) / (r * (r2 * r2 + a2 * x[IMC3] * x[IMC3])) * f;

  // Calculate vector derivatives
  Real dl0_dx = 0.0;
  Real dl0_dy = 0.0;
  Real dl0_dz = 0.0;
  Real dl1_dx = ((x[IMC1] - 2.0 * r * l1) * dr_dx + r) / (r2 + a2);
  Real dl1_dy = ((x[IMC1] - 2.0 * r * l1) * dr_dy + a) / (r2 + a2);
  Real dl1_dz = (x[IMC1] - 2.0 * r * l1) * dr_dz / (r2 + a2);
  Real dl2_dx = ((x[IMC2] - 2.0 * r * l2) * dr_dx - a) / (r2 + a2);
  Real dl2_dy = ((x[IMC2] - 2.0 * r * l2) * dr_dy + r) / (r2 + a2);
  Real dl2_dz = (x[IMC2] - 2.0 * r * l2) * dr_dz / (r2 + a2);
  Real dl3_dx = -x[IMC3] / r2 * dr_dx;
  Real dl3_dy = -x[IMC3] / r2 * dr_dy;
  Real dl3_dz = -x[IMC3] / r2 * dr_dz + 1.0 / r;

  // Calculate metric component x-derivatives
  dgcon[IMC1][IMC0][IMC0] = -(df_dx * l0 * l0 + f * dl0_dx * l0 + f * l0 * dl0_dx);
  dgcon[IMC1][IMC0][IMC1] = -(df_dx * l0 * l1 + f * dl0_dx * l1 + f * l0 * dl1_dx);
  dgcon[IMC1][IMC0][IMC2] = -(df_dx * l0 * l2 + f * dl0_dx * l2 + f * l0 * dl2_dx);
  dgcon[IMC1][IMC0][IMC3] = -(df_dx * l0 * l3 + f * dl0_dx * l3 + f * l0 * dl3_dx);
  dgcon[IMC1][IMC1][IMC0] = -(df_dx * l1 * l0 + f * dl1_dx * l0 + f * l1 * dl0_dx);
  dgcon[IMC1][IMC1][IMC1] = -(df_dx * l1 * l1 + f * dl1_dx * l1 + f * l1 * dl1_dx);
  dgcon[IMC1][IMC1][IMC2] = -(df_dx * l1 * l2 + f * dl1_dx * l2 + f * l1 * dl2_dx);
  dgcon[IMC1][IMC1][IMC3] = -(df_dx * l1 * l3 + f * dl1_dx * l3 + f * l1 * dl3_dx);
  dgcon[IMC1][IMC2][IMC0] = -(df_dx * l2 * l0 + f * dl2_dx * l0 + f * l2 * dl0_dx);
  dgcon[IMC1][IMC2][IMC1] = -(df_dx * l2 * l1 + f * dl2_dx * l1 + f * l2 * dl1_dx);
  dgcon[IMC1][IMC2][IMC2] = -(df_dx * l2 * l2 + f * dl2_dx * l2 + f * l2 * dl2_dx);
  dgcon[IMC1][IMC2][IMC3] = -(df_dx * l2 * l3 + f * dl2_dx * l3 + f * l2 * dl3_dx);
  dgcon[IMC1][IMC3][IMC0] = -(df_dx * l3 * l0 + f * dl3_dx * l0 + f * l3 * dl0_dx);
  dgcon[IMC1][IMC3][IMC1] = -(df_dx * l3 * l1 + f * dl3_dx * l1 + f * l3 * dl1_dx);
  dgcon[IMC1][IMC3][IMC2] = -(df_dx * l3 * l2 + f * dl3_dx * l2 + f * l3 * dl2_dx);
  dgcon[IMC1][IMC3][IMC3] = -(df_dx * l3 * l3 + f * dl3_dx * l3 + f * l3 * dl3_dx);

  // Calculate metric component y-derivatives
  dgcon[IMC2][IMC0][IMC0] = -(df_dy * l0 * l0 + f * dl0_dy * l0 + f * l0 * dl0_dy);
  dgcon[IMC2][IMC0][IMC1] = -(df_dy * l0 * l1 + f * dl0_dy * l1 + f * l0 * dl1_dy);
  dgcon[IMC2][IMC0][IMC2] = -(df_dy * l0 * l2 + f * dl0_dy * l2 + f * l0 * dl2_dy);
  dgcon[IMC2][IMC0][IMC3] = -(df_dy * l0 * l3 + f * dl0_dy * l3 + f * l0 * dl3_dy);
  dgcon[IMC2][IMC1][IMC0] = -(df_dy * l1 * l0 + f * dl1_dy * l0 + f * l1 * dl0_dy);
  dgcon[IMC2][IMC1][IMC1] = -(df_dy * l1 * l1 + f * dl1_dy * l1 + f * l1 * dl1_dy);
  dgcon[IMC2][IMC1][IMC2] = -(df_dy * l1 * l2 + f * dl1_dy * l2 + f * l1 * dl2_dy);
  dgcon[IMC2][IMC1][IMC3] = -(df_dy * l1 * l3 + f * dl1_dy * l3 + f * l1 * dl3_dy);
  dgcon[IMC2][IMC2][IMC0] = -(df_dy * l2 * l0 + f * dl2_dy * l0 + f * l2 * dl0_dy);
  dgcon[IMC2][IMC2][IMC1] = -(df_dy * l2 * l1 + f * dl2_dy * l1 + f * l2 * dl1_dy);
  dgcon[IMC2][IMC2][IMC2] = -(df_dy * l2 * l2 + f * dl2_dy * l2 + f * l2 * dl2_dy);
  dgcon[IMC2][IMC2][IMC3] = -(df_dy * l2 * l3 + f * dl2_dy * l3 + f * l2 * dl3_dy);
  dgcon[IMC2][IMC3][IMC0] = -(df_dy * l3 * l0 + f * dl3_dy * l0 + f * l3 * dl0_dy);
  dgcon[IMC2][IMC3][IMC1] = -(df_dy * l3 * l1 + f * dl3_dy * l1 + f * l3 * dl1_dy);
  dgcon[IMC2][IMC3][IMC2] = -(df_dy * l3 * l2 + f * dl3_dy * l2 + f * l3 * dl2_dy);
  dgcon[IMC2][IMC3][IMC3] = -(df_dy * l3 * l3 + f * dl3_dy * l3 + f * l3 * dl3_dy);

  // Calculate metric component z-derivatives
  dgcon[IMC3][IMC0][IMC0] = -(df_dz * l0 * l0 + f * dl0_dz * l0 + f * l0 * dl0_dz);
  dgcon[IMC3][IMC0][IMC1] = -(df_dz * l0 * l1 + f * dl0_dz * l1 + f * l0 * dl1_dz);
  dgcon[IMC3][IMC0][IMC2] = -(df_dz * l0 * l2 + f * dl0_dz * l2 + f * l0 * dl2_dz);
  dgcon[IMC3][IMC0][IMC3] = -(df_dz * l0 * l3 + f * dl0_dz * l3 + f * l0 * dl3_dz);
  dgcon[IMC3][IMC1][IMC0] = -(df_dz * l1 * l0 + f * dl1_dz * l0 + f * l1 * dl0_dz);
  dgcon[IMC3][IMC1][IMC1] = -(df_dz * l1 * l1 + f * dl1_dz * l1 + f * l1 * dl1_dz);
  dgcon[IMC3][IMC1][IMC2] = -(df_dz * l1 * l2 + f * dl1_dz * l2 + f * l1 * dl2_dz);
  dgcon[IMC3][IMC1][IMC3] = -(df_dz * l1 * l3 + f * dl1_dz * l3 + f * l1 * dl3_dz);
  dgcon[IMC3][IMC2][IMC0] = -(df_dz * l2 * l0 + f * dl2_dz * l0 + f * l2 * dl0_dz);
  dgcon[IMC3][IMC2][IMC1] = -(df_dz * l2 * l1 + f * dl2_dz * l1 + f * l2 * dl1_dz);
  dgcon[IMC3][IMC2][IMC2] = -(df_dz * l2 * l2 + f * dl2_dz * l2 + f * l2 * dl2_dz);
  dgcon[IMC3][IMC2][IMC3] = -(df_dz * l2 * l3 + f * dl2_dz * l3 + f * l2 * dl3_dz);
  dgcon[IMC3][IMC3][IMC0] = -(df_dz * l3 * l0 + f * dl3_dz * l0 + f * l3 * dl0_dz);
  dgcon[IMC3][IMC3][IMC1] = -(df_dz * l3 * l1 + f * dl3_dz * l1 + f * l3 * dl1_dz);
  dgcon[IMC3][IMC3][IMC2] = -(df_dz * l3 * l2 + f * dl3_dz * l2 + f * l3 * dl2_dz);
  dgcon[IMC3][IMC3][IMC3] = -(df_dz * l3 * l3 + f * dl3_dz * l3 + f * l3 * dl3_dz);

}

//----------------------------------------------------------------------------------------
//! \fn void MCKerrSchildCartesian::Connect(Real x[4], Real gcov[4][4])
//! \brief compute metric in cartesian Kerr-Schild

void MCKerrSchildCartesian::Connect(Real x[4], Real gamma[4][4][4]) {


}

//----------------------------------------------------------------------------------------
//! MCBoyerLindquist, built from Coord and MonteCarloBlock

MCBoyerLindquist::MCBoyerLindquist(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {
}

//----------------------------------------------------------------------------------------
//! MCBoyerLindquist constructor for processes without own MeshBlock

MCBoyerLindquist::MCBoyerLindquist(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {
}

//----------------------------------------------------------------------------------------
//! destructor
MCBoyerLindquist::~MCBoyerLindquist() {

}

//----------------------------------------------------------------------------------------
//! \fn void MCBoyerLindquist::Metric(Real x[4], Real gcov[4][4])
//! \brief compute metric for Boyer-Lindquist coordinates

void MCBoyerLindquist::Metric(Real x[4], Real gcov[4][4]) {

  // equation for the Metric comes from the inside cover of Hartle
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
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
  gcov[IMC0][IMC3] = -2.*a*r*sth2/rho2;
  gcov[IMC3][IMC0] = gcov[IMC0][IMC3];
}

//----------------------------------------------------------------------------------------
//! \fn void MCBoyerLindquist::InverseMetric(Real x[4], Real gcon[4][4])
//! \brief compute inverse metric for Boyer-Lindquist coordinates

void MCBoyerLindquist::InverseMetric(Real x[4], Real gcon[4][4]) {

  // Equation comes from ColinsCosmos.com/wiki/boyer-lindquist-coordinates, which
  // sites Frolov & Novikov Section D.1 (but I don't have access to this book)

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
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

//----------------------------------------------------------------------------------------
//! \fn void MCBoyerLindquist::Connect(Real x[4], Real gamma[4][4][4])
//! \brief compute connection for Boyer-Lindquist coordinates

void MCBoyerLindquist::Connect(Real x[4], Real gamma[4][4][4]) {

  // equations for the connection coefficients come from Frutos-Alfaro et al. (2012)

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
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

  gamma[IMC3][IMC1][IMC0] = gamma[IMC3][IMC0][IMC1];
  gamma[IMC3][IMC1][IMC3] = 1./(rho4*delta)*(r*rho2*(rho2-rs*r)-a*j*sth2*(2.*r2-rho2));

  gamma[IMC3][IMC2][IMC0] = gamma[IMC3][IMC0][IMC2];
  gamma[IMC3][IMC2][IMC3] = cth/(rho4*sth)*(rho4+2.*a*j*r*sth2);

  gamma[IMC3][IMC3][IMC1] = gamma[IMC3][IMC1][IMC3];
  gamma[IMC3][IMC3][IMC2] = gamma[IMC3][IMC2][IMC3];
}

//----------------------------------------------------------------------------------------
//! MCMinkowski constructor, built from Coord and MonteCarloBlock

MCMinkowski::MCMinkowski(Coordinates *pcoord, MonteCarloBlock *pmcb)
  : MCCoord(pcoord,pmcb) {
}

//----------------------------------------------------------------------------------------
//! MCMinkowski constructor for processes without own MeshBlock
MCMinkowski::MCMinkowski(int ncells1, int ncells2, int ncells3, bool acc)
  : MCCoord(ncells1,ncells2,ncells3,acc) {
}

//----------------------------------------------------------------------------------------
//! destructor
MCMinkowski::~MCMinkowski() {

}
