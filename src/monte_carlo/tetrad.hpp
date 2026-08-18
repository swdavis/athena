#ifndef MONTE_CARLO_TETRAD_HPP_
#define MONTE_CARLO_TETRAD_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file tetrad.hpp
//! \brief orthonormal frames and the vector algebra that builds them.
//!
//! These are geometric primitives: they act on bare four-vectors and a metric, and know
//! nothing about photons, mesh blocks or radiation moments.  Applying them to Monte Carlo
//! quantities -- projecting a photon into a named frame, transforming accumulated moments
//! between frames -- lives in photon_frames.cpp, which is the layer above this one.
//!
//! Implemented in tetrad.cpp.

// C++ headers
#include <complex>

// Athena++ headers
#include "../athena.hpp"

//! build an orthonormal tetrad by Gram-Schmidt against one, two or three given
//! directions; the timelike leg is taken parallel to ucon and the remaining legs are
//! filled out from coordinate directions.
void ConstructTetrad(Real ucon[4], Real gcov[4][4],
                     Real econ[4][4], Real ecov[4][4]);
void ConstructTetrad(Real ucon[4], Real vcon[4], Real gcov[4][4],
                     Real econ[4][4], Real ecov[4][4]);
void ConstructTetrad(Real ucon[4], Real vcon[4], Real wcon[4],
                     Real gcov[4][4], Real econ[4][4],
                     Real ecov[4][4]);
void InitializeLeviCivita(Real levi[4][4][4][4]);
void ImposeRightHanded(Real econ[4][4], Real gcov[4][4]);
Real KroneckerDelta(int i, int j);

//! vector algebra in a metric
void ProjectVecSub(Real ucon[4], Real vcon[4], Real gcov[4][4]);
Real DotVec(Real ucon[4], Real vcon[4], Real gcov[4][4]);
void NormalizeVec(Real ucon[4], Real gcov[4][4]);
void ConToCov(Real ucon[4], Real ucov[4], Real gcov[4][4]);
void CovToCon(Real ucov[4], Real ucon[4], Real gcon[4][4]);

//! move a four-vector between the coordinate basis and a tetrad
void CoordinateToTetrad(Real ucoord[4], Real utet[4], Real ecov[4][4]);
void TetradToCoordinate(Real utet[4], Real ucoord[4], Real econ[4][4]);

//! polarization, carried as an invariant tensor between scatterings
void StokesToTensor(Real stokes[4], std::complex<Real> tensor[4][4]);
void TensorToStokes(std::complex<Real> tensor[4][4], Real stokes[4]);

void LorentzBoostVector(Real vel[4], Real kold[4]);

#endif // MONTE_CARLO_TETRAD_HPP_
