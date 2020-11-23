//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file frame_transformations.hpp
//  \brief implementation of functions for constructing and transforming tetrads
//         frame

#include <complex>

// Athena++ classes headers
#include "montecarlo.hpp"

#define SMALL_NUMBER 1.e-30
//#define DEBUG

//----------------------------------------------------------------------------------------
//! \fn void ConstructTetrad(Real ucon[NCOORD], Real gcov[NCOORD][NCOORD],
//                           Real ecov[NCOORD][NCOORD], Real econ[NCOORD][NCOORD]
//  \brief construct an orthonormal tetrad using the fluid frame vector ucon with
//         the Gram-Schmidt algorithm

void ConstructTetrad(Real ucon[NCOORD], Real gcov[NCOORD][NCOORD], 
                     Real econ[NCOORD][NCOORD], Real ecov[NCOORD][NCOORD]) {

  // make a trial tetrad where time component is parallel to fluid, the first spatial
  // coordinate is from the magnetic field and the last two are diagonal
  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      if (i == IMC0)
        econ[IMC0][j] = ucon[j]; // time component parallel to fluid frame
      else 
        econ[i][j] = static_cast<Real>(KroneckerDelta(i,j)); // diagonal trial 
    }
  }

#ifdef DEBUG
  printf("\nTrial:\n");
  printf("econ[IMC0]: %g %g %g %g\n", econ[IMC0][IMC0], econ[IMC0][IMC1], econ[IMC0][IMC2], econ[IMC0][IMC3]);
  printf("econ[IMC1]: %g %g %g %g\n", econ[IMC1][IMC0], econ[IMC1][IMC1], econ[IMC1][IMC2], econ[IMC1][IMC3]);
  printf("econ[IMC2]: %g %g %g %g\n", econ[IMC2][IMC0], econ[IMC2][IMC1], econ[IMC2][IMC2], econ[IMC2][IMC3]);
  printf("econ[IMC3]: %g %g %g %g\n", econ[IMC3][IMC0], econ[IMC3][IMC1], econ[IMC3][IMC2], econ[IMC3][IMC3]);
  printf("ucon: %g %g %g %g\n", ucon[IMC0], ucon[IMC1], ucon[IMC2], ucon[IMC3]);
  printf("\n");
#endif

  // begin constructing contravariant tetrad
  NormalizeVec(econ[IMC0], gcov);

  ProjectVecSub(econ[IMC1], econ[IMC0], gcov);
  NormalizeVec(econ[IMC1], gcov);

  ProjectVecSub(econ[IMC2], econ[IMC0], gcov);
  ProjectVecSub(econ[IMC2], econ[IMC1], gcov);
  NormalizeVec(econ[IMC2], gcov);

  ProjectVecSub(econ[IMC3], econ[IMC0], gcov);
  ProjectVecSub(econ[IMC3], econ[IMC1], gcov);
  ProjectVecSub(econ[IMC3], econ[IMC2], gcov);
  NormalizeVec(econ[IMC3], gcov);
  // contravariant tetrad construction complete

  // begin construction covariant tetrad
  for (int i = 0; i < NCOORD; i++) {
    ConToCov(econ[i], ecov[i], gcov);
    if (i == IMC0) {
      for (int j = 0; j < NCOORD; j++) ecov[IMC0][j] *= -1;
    }
  }
  // covariant tetrad construction complete

  // check for orthonormality
#ifdef DEBUG
  printf("--------------------------------------------------\n");
  printf("Othonormality check:\n");
  printf("econ[IMC0]: %g %g %g %g\n", econ[IMC0][IMC0], econ[IMC0][IMC1], econ[IMC0][IMC2], econ[IMC0][IMC3]);
  printf("econ[IMC1]: %g %g %g %g\n", econ[IMC1][IMC0], econ[IMC1][IMC1], econ[IMC1][IMC2], econ[IMC1][IMC3]);
  printf("econ[IMC2]: %g %g %g %g\n", econ[IMC2][IMC0], econ[IMC2][IMC1], econ[IMC2][IMC2], econ[IMC2][IMC3]);
  printf("econ[IMC3]: %g %g %g %g\n", econ[IMC3][IMC0], econ[IMC3][IMC1], econ[IMC3][IMC2], econ[IMC3][IMC3]);
  printf("ecov[IMC0]: %g %g %g %g\n", ecov[IMC0][IMC0], ecov[IMC0][IMC1], ecov[IMC0][IMC2], ecov[IMC0][IMC3]);
  printf("ecov[IMC1]: %g %g %g %g\n", ecov[IMC1][IMC0], ecov[IMC1][IMC1], ecov[IMC1][IMC2], ecov[IMC1][IMC3]);
  printf("ecov[IMC2]: %g %g %g %g\n", ecov[IMC2][IMC0], ecov[IMC2][IMC1], ecov[IMC2][IMC2], ecov[IMC2][IMC3]);
  printf("ecov[IMC3]: %g %g %g %g\n", ecov[IMC3][IMC0], ecov[IMC3][IMC1], ecov[IMC3][IMC2], ecov[IMC3][IMC3]);
  printf("gcov[IMC0]: %g %g %g %g\n", gcov[IMC0][IMC0], gcov[IMC0][IMC1], gcov[IMC0][IMC2], gcov[IMC0][IMC3]);
  printf("gcov[IMC1]: %g %g %g %g\n", gcov[IMC1][IMC0], gcov[IMC1][IMC1], gcov[IMC1][IMC2], gcov[IMC1][IMC3]);
  printf("gcov[IMC2]: %g %g %g %g\n", gcov[IMC2][IMC0], gcov[IMC2][IMC1], gcov[IMC2][IMC2], gcov[IMC2][IMC3]);
  printf("gcov[IMC3]: %g %g %g %g\n", gcov[IMC3][IMC0], gcov[IMC3][IMC1], gcov[IMC3][IMC2], gcov[IMC3][IMC3]);

  Real sum;
  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      sum = 0.;
      for (int k = 0; k < NCOORD; k++) {
	sum += econ[i][k] * ecov[j][k];
      } 
    }
    printf("sum: %g  index: %d\n",sum, i);


  }
  printf("--------------------------------------------------\n");
#endif

}

//----------------------------------------------------------------------------------------
//! \fn void ConstructTetrad(Real ucon[NCOORD], Real kcon[NCOORD], 
//                              Real gcov[NCOORD][NCOORD],
//                              Real ecov[NCOORD][NCOORD], Real econ[NCOORD][NCOORD]
//  \brief construct an orthonormal tetrad using the fluid frame vector ucon and defining
//         e^mu_(3) = k^\mu 


void ConstructTetrad(Real ucon[NCOORD], Real kcon[NCOORD], Real gcov[NCOORD][NCOORD], 
                     Real econ[NCOORD][NCOORD], Real ecov[NCOORD][NCOORD]) {

 // make a trial tetrad where time component is parallel to fluid, the first spatial
  // coordinate is from the magnetic field and the last two are diagonal
  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      if (i == IMC0)
        econ[IMC0][j] = ucon[j]; // time component parallel to fluid frame
      else if (i == IMC3)
        econ[IMC3][j] = kcon[j]; // set by k vector
      else
        econ[i][j] = static_cast<Real>(KroneckerDelta(i,j)); // diagonal trial 
    }
  }

  // begin constructing contravariant tetrad
  NormalizeVec(econ[IMC0], gcov);

  ProjectVecSub(econ[IMC3], econ[IMC0], gcov);
  NormalizeVec(econ[IMC3], gcov);

  ProjectVecSub(econ[IMC1], econ[IMC0], gcov);
  ProjectVecSub(econ[IMC1], econ[IMC3], gcov);
  NormalizeVec(econ[IMC1], gcov);

  ProjectVecSub(econ[IMC2], econ[IMC0], gcov);
  ProjectVecSub(econ[IMC2], econ[IMC3], gcov);
  ProjectVecSub(econ[IMC2], econ[IMC1], gcov);
  NormalizeVec(econ[IMC2], gcov);
  // contravariant tetrad construction complete

  // begin construction covariant tetrad
  for (int i = 0; i < NCOORD; i++) {
    ConToCov(econ[i], ecov[i], gcov);
    if (i == IMC0) {
      for (int j = 0; j < NCOORD; j++) ecov[IMC0][j] *= -1;
    }
  }
#ifdef DEBUG
  Real sum;
  for (int i = 0; i < NCOORD; i++) {
    for (int j = 0; j < NCOORD; j++) {
      sum = 0.;
      for (int k = 0; k < NCOORD; k++) {
	sum += econ[i][k] * ecov[j][k];
      } 
    }
    printf("sum: %g  index: %d\n",sum, i);
  }
  printf("Othonormality check:\n");
  printf("econ[IMC0]: %g %g %g %g\n", econ[IMC0][IMC0], econ[IMC0][IMC1], econ[IMC0][IMC2], econ[IMC0][IMC3]);
  printf("econ[IMC1]: %g %g %g %g\n", econ[IMC1][IMC0], econ[IMC1][IMC1], econ[IMC1][IMC2], econ[IMC1][IMC3]);
  printf("econ[IMC2]: %g %g %g %g\n", econ[IMC2][IMC0], econ[IMC2][IMC1], econ[IMC2][IMC2], econ[IMC2][IMC3]);
  printf("econ[IMC3]: %g %g %g %g\n", econ[IMC3][IMC0], econ[IMC3][IMC1], econ[IMC3][IMC2], econ[IMC3][IMC3]);
  printf("ecov[IMC0]: %g %g %g %g\n", ecov[IMC0][IMC0], ecov[IMC0][IMC1], ecov[IMC0][IMC2], ecov[IMC0][IMC3]);
  printf("ecov[IMC1]: %g %g %g %g\n", ecov[IMC1][IMC0], ecov[IMC1][IMC1], ecov[IMC1][IMC2], ecov[IMC1][IMC3]);
  printf("ecov[IMC2]: %g %g %g %g\n", ecov[IMC2][IMC0], ecov[IMC2][IMC1], ecov[IMC2][IMC2], ecov[IMC2][IMC3]);
  printf("ecov[IMC3]: %g %g %g %g\n", ecov[IMC3][IMC0], ecov[IMC3][IMC1], ecov[IMC3][IMC2], ecov[IMC3][IMC3]);
  printf("gcov[IMC0]: %g %g %g %g\n", gcov[IMC0][IMC0], gcov[IMC0][IMC1], gcov[IMC0][IMC2], gcov[IMC0][IMC3]);
  printf("gcov[IMC1]: %g %g %g %g\n", gcov[IMC1][IMC0], gcov[IMC1][IMC1], gcov[IMC1][IMC2], gcov[IMC1][IMC3]);
  printf("gcov[IMC2]: %g %g %g %g\n", gcov[IMC2][IMC0], gcov[IMC2][IMC1], gcov[IMC2][IMC2], gcov[IMC2][IMC3]);
  printf("gcov[IMC3]: %g %g %g %g\n", gcov[IMC3][IMC0], gcov[IMC3][IMC1], gcov[IMC3][IMC2], gcov[IMC3][IMC3]);
#endif
}


//----------------------------------------------------------------------------------------
//! \fn int KroneckerDelta(int i, int j)
//  \brief Kronecker Delta function

int KroneckerDelta(int i, int j) {

  if (i == j) return 1;
  else return 0;

}

//----------------------------------------------------------------------------------------
//! \fn void ProjectVecSub(Real ucon[NCOORD], Real vcon[NCOORD], 
//                          Real gcov[NCOORD][NCOORD])
//  \brief subtract four-vector projection of ucon onto vcon from ucon in metric gcov

void ProjectVecSub(Real ucon[NCOORD], Real vcon[NCOORD], Real gcov[NCOORD][NCOORD]) {

  Real vcon_dot_vcon = DotVec(vcon, vcon, gcov);
  Real ucon_dot_vcon = DotVec(ucon, vcon, gcov);

  if (fabs(vcon_dot_vcon) < SMALL_NUMBER) {
    printf("Warning: attempted to project out using a zero vector. Vector left as is.\n");
#ifdef DEBUG
    printf("vdotv: %g\n", vcon_dot_vcon);
    printf("vcon: %g %g %g %g\n", vcon[IMC0], vcon[IMC1], vcon[IMC2], vcon[IMC3]);
    printf("gcov[IMC0]: %g %g %g %g\n", gcov[IMC0][IMC0], gcov[IMC0][IMC1], gcov[IMC0][IMC2], gcov[IMC0][IMC3]);
    printf("gcov[IMC1]: %g %g %g %g\n", gcov[IMC1][IMC0], gcov[IMC1][IMC1], gcov[IMC1][IMC2], gcov[IMC1][IMC3]);
    printf("gcov[IMC2]: %g %g %g %g\n", gcov[IMC2][IMC0], gcov[IMC2][IMC1], gcov[IMC2][IMC2], gcov[IMC2][IMC3]);
    printf("gcov[IMC3]: %g %g %g %g\n", gcov[IMC3][IMC0], gcov[IMC3][IMC1], gcov[IMC3][IMC2], gcov[IMC3][IMC3]);
#endif
    return;
  }
  
  for (int i = 0; i < NCOORD; i++) 
    ucon[i] -= vcon[i] * ucon_dot_vcon / vcon_dot_vcon;

}

//----------------------------------------------------------------------------------------
//! \fn Real DotVec(Real ucon[NCOORD], Real vcon[NCOORD], Real gcov[NCOORD][NCOORD])
//  \brief return dot (scalar) product of <ucon, vcon> in metric gcov

Real DotVec(Real ucon[NCOORD], Real vcon[NCOORD], Real gcov[NCOORD][NCOORD]) {

  Real dot = 0;

  for (int i = 0; i < NCOORD; i++) 
    for (int j = 0; j < NCOORD; j++) 
      dot += ucon[i] * gcov[i][j] * vcon[j];

  return dot;

}

//----------------------------------------------------------------------------------------
//! \fn void NormalizeVec(Real ucon[NCOORD], Real gcov[NCOORD][NCOORD])
//  \brief normalize vector ucon in metric gcov

void NormalizeVec(Real ucon[NCOORD], Real gcov[NCOORD][NCOORD]) {

  double mag = sqrt(fabs(DotVec(ucon, ucon, gcov)));

  if (mag < SMALL_NUMBER) {
    printf("Warning: attempted to normalize a zero vector. Vector left as is.\n");
    return;
  }

  for (int i = 0; i < NCOORD; i ++)
    ucon[i] /= mag;

}

//----------------------------------------------------------------------------------------
//! \fn void ConToCov(Real ucon[NCOORD], Real ucov[NCOORD], Real gcov[NCOORD][NCOORD])
//  \brief converts contravariant vector ucon to covariant vector ucov in metric gcov

void ConToCov(Real ucon[NCOORD], Real ucov[NCOORD], Real gcov[NCOORD][NCOORD]) {

  for (int i = 0; i < NCOORD; i++) {
    ucov[i] = 0.;
    for (int j = 0; j < NCOORD; j++) 
      ucov[i] += gcov[i][j] * ucon[j];
  }

}

//----------------------------------------------------------------------------------------
//! \fn void CovToCon(Real ucov[NCOORD], Real ucon[NCOORD], Real gcon[NCOORD][NCOORD])
//  \brief converts covariant vector ucov to contravariant vector ucon in metric gcon

void CovToCon(Real ucov[NCOORD], Real ucon[NCOORD], Real gcon[NCOORD][NCOORD]) {

  for (int i = 0; i < NCOORD; i++) {
    ucon[i] = 0.;
    for (int j = 0; j < NCOORD; j++) 
      ucon[i] += gcon[i][j] * ucov[j];
  }

}

//----------------------------------------------------------------------------------------
//! \fn void CoordinateToTetrad(Real ucoord[NCOORD], Real utet[NCOORD],
//                              Real ecov[NCOORD][NCOORD])
//  \brief transform vector ucoord from coordinate frame to utet in the tetrad frame

void CoordinateToTetrad(Real ucoord[NCOORD], Real utet[NCOORD], 
			Real ecov[NCOORD][NCOORD]) {
  
  for (int i = 0; i < NCOORD; i++) {
    utet[i] = 0.;
    for (int j = 0; j< NCOORD; j++) 
      utet[i] += ecov[i][j] * ucoord[j];
  }

}


//----------------------------------------------------------------------------------------
//! \fn void TetradToCoordinate(Real utet[NCOORD], Real ucoord[NCOORD], 
//                              Real econ[NCOORD][NCOORD])
//  \brief transform vector utet from tetrad frame to ucoord in the coordinate frame

void TetradToCoordinate(Real utet[NCOORD], Real ucoord[NCOORD], 
			Real econ[NCOORD][NCOORD]) {

  for (int i = 0; i < NCOORD; i++) {
    ucoord[i] = 0.;
    for (int j = 0; j < NCOORD; j++) 
      ucoord[i] += econ[j][i] * utet[j];
  }

}

// SWD: make purely real version for linear polarization
//----------------------------------------------------------------------------------------
//! \fn void StokesToTensor(Real stokes[NCOORD],std::complex<Real> tensor[NCOORD][NCOORD])
//  \brief transform stokes vector invariant polarization tensor

void StokesToTensor(Real stokes[NCOORD], std::complex<Real> tensor[NCOORD][NCOORD]) {

  std::complex<Real> I = std::complex<Real>(0.,1.);

  for (int i = 0; i < NCOORD; i++)
    for (int j = 0; j < NCOORD; j++)
      tensor[i][j] = std::complex<Real>(0.,0.);

  //tensor[IMC2][IMC2] = (stokes[0] + stokes[1]);
  //tensor[IMC2][IMC3] = (stokes[2] - I * stokes[3]);
  //tensor[IMC3][IMC2] = (stokes[2] + I * stokes[3]);
  //tensor[IMC3][IMC3] = (stokes[0] - stokes[1]);
  // SWD general for tetrad in MC18
  tensor[IMC1][IMC1] = (stokes[0] + stokes[1]);
  tensor[IMC1][IMC2] = (stokes[2] - I * stokes[3]);
  tensor[IMC2][IMC1] = (stokes[2] + I * stokes[3]);
  tensor[IMC2][IMC2] = (stokes[0] - stokes[1]);
}

//----------------------------------------------------------------------------------------
//! \fn void TensorToStokes(Real stokes[NCOORD],std::complex<Real> tensor[NCOORD][NCOORD])
//  \brief transform invariant polarization tensor to stokes vector

void TensorToStokes(std::complex<Real> tensor[NCOORD][NCOORD], Real stokes[NCOORD]) {

  //stokes[0] = 0.5 * (tensor[IMC2][IMC2] + tensor[IMC3][IMC3]).real();
  //stokes[1] = 0.5 * (tensor[IMC2][IMC2] - tensor[IMC3][IMC3]).real();
  //stokes[2] = 0.5 * (tensor[IMC2][IMC3] + tensor[IMC3][IMC2]).real();
  //stokes[3] = 0.5 * (tensor[IMC3][IMC2] - tensor[IMC2][IMC3]).imag();
  //SWD general for tetrad in MC18
  stokes[0] = 0.5 * (tensor[IMC1][IMC1] + tensor[IMC2][IMC2]).real();
  stokes[1] = 0.5 * (tensor[IMC1][IMC1] - tensor[IMC2][IMC2]).real();
  stokes[2] = 0.5 * (tensor[IMC1][IMC2] + tensor[IMC2][IMC1]).real();
  stokes[3] = 0.5 * (tensor[IMC2][IMC1] - tensor[IMC1][IMC2]).imag();

  Real norm = stokes[0];
  for (int i = 0; i < NCOORD; i++)
    stokes[i] /= norm;

}

//----------------------------------------------------------------------------------------
//! \fn void ComplexCoordinateToTetrad(std::complex<Real> tensor[NCOORD][NCOORD],
//                             std::complex<Real> polten[NCOORD][NCOORD],
//                             Real Ecov[NCOORD][NCOORD])
//  \brief transform complex tensor from coordinate frame to tetrad frame


void ComplexCoordinateToTetrad(std::complex<Real> tcoord[NCOORD][NCOORD], 
                               std::complex<Real> ttet[NCOORD][NCOORD],
                               Real ecov[NCOORD][NCOORD])
{
 
  for (int i = 0; i < NCOORD; i++)
    for (int j = 0; j < NCOORD; j++)
      ttet[i][j] = std::complex<Real>(0.,0.);

  for (int i = 0; i < NCOORD; i++)
    for (int j = 0; j < NCOORD; j++)
      for (int k = 0; k < NCOORD; k++)
	for (int l = 0; l < NCOORD; l++)
	  ttet[i][j] += tcoord[k][l] * ecov[i][k] * ecov[j][l];

}


//----------------------------------------------------------------------------------------
//! \fn void ComplexTetradToCoordinate(std::complex<Real> ttet[NCOORD][NCOORD],
//                                     std::complex<Real> tcoord[NCOORD][NCOORD],
//                                     Real econ[NCOORD][NCOORD])
//  \brief transform complex tensor from tetrad frame to coordinate frame

void ComplexTetradToCoordinate(std::complex<Real> ttet[NCOORD][NCOORD],
                               std::complex<Real> tcoord[NCOORD][NCOORD],
                               Real econ[NCOORD][NCOORD]) {

  for(int i = 0; i < NCOORD; i++)
    for(int j = 0; j < NCOORD; j++)
      tcoord[i][j] = std::complex<Real>(0.,0.);

  for(int i = 0; i < NCOORD; i++)
    for(int j = 0; j < NCOORD; j++)
      for(int k = 0; k < NCOORD; k++)
	for(int l = 0; l < NCOORD; l++)
	  tcoord[i][j] += ttet[k][l] * econ[k][i] * econ[l][j];

}



			     
			    
							   
							    
