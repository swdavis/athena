//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mc_snake.cpp
//! \brief Uniform periodic box in sinusoidal ("snake") coordinates, built with -g.
//
// Snake coordinates are flat spacetime written in the sheared frame of White, Stone &
// Gammie (2016), ApJS 225, 22:
//
//     t = t_M,  x = x_M,  y = y_M + a sin(k x_M),  z = z_M
//
// so with beta = a k cos(k x),
//
//     ds^2 = -dt^2 + (1 + beta^2) dx^2 - 2 beta dx dy + dy^2 + dz^2.
//
// This pgen uses the same physical problem as mc_isoth_mink -- a uniform box,
// free-free emission, optionally drifting along x3. Because the spacetime is genuinely flat,
// every frame-invariant scalar must agree between the two runs to round-off.
//
// Units follow mc_isoth_mink exactly: hydro primitives are dimensionless code values,
// converted for the Monte Carlo by rho_cgs and tgas_cgs, because the GR equation of
// state needs p/rho << 1.  See that file for the arithmetic.
//
// Inputs (in addition to the usual mc_isoth ones):
//   <coord>/snake_a   shear amplitude          (default 0 -> degenerates to Minkowski)
//   <coord>/snake_k   shear wavenumber         (default 0)
//   <coord>/m, a      required by GRUser and ignored here; set both to 0

// C++ headers
#include <cmath>
#include <cstring>  // strcmp
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#if !MONTE_CARLO_ENABLED
#error "This problem requires monte carlo"
#endif
#if !GENERAL_RELATIVITY
#error "mc_snake requires general relativity (-g)"
#endif

#include "../monte_carlo/montecarlo.hpp"
#include "../monte_carlo/tetrad.hpp"
#include "../monte_carlo/photon.hpp"
#include "../monte_carlo/photonpusher.hpp"

namespace {
bool tnorm;
Real logemin, logemax;
// Cached so the metric callback, which runs per cell, does not re-parse the input.
// Set in Mesh::InitUserMeshData before the metric is ever evaluated.
Real snake_a = 0.0;
Real snake_k = 0.0;
void SnakeMetric(Real x1, Real x2, Real x3, ParameterInput *pin,
                 AthenaArray<Real> &g, AthenaArray<Real> &g_inv,
                 AthenaArray<Real> &dg_dx1, AthenaArray<Real> &dg_dx2,
                 AthenaArray<Real> &dg_dx3);
void JMeanOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip, int imom,
                  const PhotonFrameState &s);
void AverageEnergy(MonteCarloBlock *pmcb, Photon *pphot, int ip, int imom,
                  const PhotonFrameState &s);
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief enroll the snake metric and check the box is compatible with it
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {

  if (std::strcmp(COORDINATE_SYSTEM, "gr_user") != 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in mc_snake InitUserMeshData" << std::endl
        << "mc_snake supplies its own metric and requires --coord=gr_user, got "
        << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  snake_a = pin->GetOrAddReal("coord", "snake_a", 0.0);
  snake_k = pin->GetOrAddReal("coord", "snake_k", 0.0);

  // The metric depends on x1 alone, through beta = a k cos(k x).  A periodic x1 boundary
  // therefore only makes sense when the shear closes over the box: otherwise the metric
  // is discontinuous across the boundary and the photons see a seam.  x2 and x3 are
  // always safe, since no metric component depends on them.
  std::string ix1 = pin->GetOrAddString("mesh", "ix1_mc_bc", "");
  if (ix1 == "periodic" && snake_a != 0.0 && snake_k != 0.0) {
    Real lx1 = pin->GetReal("mesh", "x1max") - pin->GetReal("mesh", "x1min");
    Real turns = snake_k * lx1 / (2.0 * PI);
    if (std::fabs(turns - std::round(turns)) > 1.0e-10) {
      std::stringstream msg;
      msg << "### FATAL ERROR in mc_snake InitUserMeshData" << std::endl
          << "Periodic x1 requires the shear to close over the box: snake_k * Lx1 must "
          << "be a multiple of 2 pi." << std::endl
          << "Got snake_k = " << snake_k << ", Lx1 = " << lx1
          << ", snake_k*Lx1/(2 pi) = " << turns << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  EnrollUserMetric(SnakeMetric);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief uniform box, optionally moving at a constant velocity along x3
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // dimensionless code values; the Monte Carlo scales them via rho_cgs/tgas_cgs
  Real rho0 = pin->GetOrAddReal("problem", "dens_code", 1.0);
  Real pgas0 = pin->GetOrAddReal("problem", "pgas_code", 1.0e-6);

  // constant velocity as a fraction of the speed of light
  Real beta = pin->GetOrAddReal("problem", "velocity", 0.);
  if (std::fabs(beta) >= 1.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in mc_snake ProblemGenerator" << std::endl
        << "velocity must be given as a fraction of c, got " << beta << std::endl;
    ATHENA_ERROR(msg);
  }
  Real lorentz = 1.0 / std::sqrt(1.0 - SQR(beta));

  // Drift is along x3 deliberately.  g_33 = 1 and no metric component mixes x3 with
  // anything, so utilde^3 = gamma*beta reproduces u^mu = (gamma,0,0,gamma*beta) exactly
  // as it does in Minkowski -- which is what lets a snake run and a mink run be
  // compared directly.  A drift along x1 or x2 would have to account for g_xy.
  Real uu3 = lorentz * beta;

  AthenaArray<Real> bb;
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->w(IDN, k, j, i) = phydro->w1(IDN, k, j, i) = rho0;
        phydro->w(IPR, k, j, i) = phydro->w1(IPR, k, j, i) = pgas0;
        phydro->w(IVX, k, j, i) = phydro->w1(IVX, k, j, i) = 0.0;
        phydro->w(IVY, k, j, i) = phydro->w1(IVY, k, j, i) = 0.0;
        phydro->w(IVZ, k, j, i) = phydro->w1(IVZ, k, j, i) = uu3;
      }
    }
  }
  peos->PrimitiveToConserved(phydro->w, bb, phydro->u, pcoord, is, ie, js, je, ks, ke);

  return;
}

//========================================================================================
//! \fn void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe, int etype)
//! \brief free-free emission, identical to mc_isoth_mink
//========================================================================================

void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe, int etype) {

  SetEmissionCellWeight(pphot, ips, ipe);

  for (int ip = ips; ip <= ipe; ip++) {

    pphot->dtp[ip] = pphot->pmy_mcb->pmy_mc->tmax;

    // position within the selected zone
    GetZonePosition(pphot, pran, pcoord, ip);
    pphot->x0p[ip] = 0.0;

    // energy, direction and weight from the free-free emissivity
    if (tnorm) {
      Real logtg = log(tgas(pphot->i3p[ip], pphot->i2p[ip], pphot->i1p[ip]));
      PhotonEmitFreeFree(this, pphot, logemin + logtg, logemax + logtg, ip);
    } else {
      PhotonEmitFreeFree(this, pphot, logemin, logemax, ip);
    }

    if (pphot->wp[ip] < 0.0)
      pphot->statp[ip] = DESTROYED;
    else
      pphot->statp[ip] = EVOLVING;

    pphot->nscp[ip] = 0;

    pphot->acp[ip] = AbsorptionOpacity(this, pphot, ip);
    pphot->scp[ip] = ScatteringOpacity(this, pphot, ip);
  }

  return;
}

//========================================================================================
//! \fn void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin)
//========================================================================================

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {

  // Set the energy boundaries for free-free emission
  tnorm = pin->GetOrAddBoolean("problem", "tnorm", false);
  if (tnorm) {
    // interpret as xmin/xmax with x = E/(kb*T)
    const Real kb = 1.380649e-16;
    logemin = log(kb * pin->GetReal("problem", "emin"));
    logemax = log(kb * pin->GetReal("problem", "emax"));
  } else {
    // interpret as emin/emax in eV
    const Real everg = 1.6021772e-12;
    logemin = log(everg * pin->GetReal("problem", "emin"));
    logemax = log(everg * pin->GetReal("problem", "emax"));
  }

  return;
}

//========================================================================================
//! \fn void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin)
//========================================================================================

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin) {

  if (!pin->GetOrAddBoolean("montecarlo", "general_pusher", false)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in mc_snake InitUserMonteCarloData" << std::endl
        << "mc_snake requires general_pusher = true" << std::endl;
    ATHENA_ERROR(msg);
  }

  // gr_user cannot say which metric the Monte Carlo should use, so mc_coord must name it
  // and must name this one.  Without the check a snake hydro metric would be paired with
  // a Kerr-Schild photon metric
  if (pin->GetOrAddString("montecarlo", "mc_coord", "") != "snake") {
    std::stringstream msg;
    msg << "### FATAL ERROR in mc_snake InitUserMonteCarloData" << std::endl
        << "mc_snake requires <montecarlo>/mc_coord = snake" << std::endl;
    ATHENA_ERROR(msg);
  }

  nuser_var = 1;
  AllocateUserMoments(2);
  EnrollUserMoment(0, JMeanOpacity, "kapJ");
  EnrollUserMoment(1, AverageEnergy, "eave");

  return;
}

namespace {

//----------------------------------------------------------------------------------------
//! \fn void SnakeMetric(...)
//! \brief the snake metric, its inverse and its coordinate derivatives
//
// Mirrors MCSnake in mccoord.cpp. The spatial block has unit determinant,
// (1+beta^2) - beta^2 = 1, so the inverse is the diagonal swapped with the off-diagonal
// sign flipped, and sqrt(-g) = 1. beta depends on x1 alone, so dg_dx2 and dg_dx3 vanish
// identically.

void SnakeMetric(Real x1, Real x2, Real x3, ParameterInput *pin,
                 AthenaArray<Real> &g, AthenaArray<Real> &g_inv,
                 AthenaArray<Real> &dg_dx1, AthenaArray<Real> &dg_dx2,
                 AthenaArray<Real> &dg_dx3) {

  Real beta = snake_a * snake_k * std::cos(snake_k * x1);
  Real dbeta = -snake_a * SQR(snake_k) * std::sin(snake_k * x1);

  for (int n = 0; n < NMETRIC; ++n) {
    g(n) = 0.0;
    g_inv(n) = 0.0;
    dg_dx1(n) = 0.0;
    dg_dx2(n) = 0.0;
    dg_dx3(n) = 0.0;
  }

  // covariant
  g(I00) = -1.0;
  g(I11) = 1.0 + SQR(beta);
  g(I12) = -beta;
  g(I22) = 1.0;
  g(I33) = 1.0;

  // contravariant
  g_inv(I00) = -1.0;
  g_inv(I11) = 1.0;
  g_inv(I12) = beta;
  g_inv(I22) = 1.0 + SQR(beta);
  g_inv(I33) = 1.0;

  // only x1 derivatives survive
  dg_dx1(I11) = 2.0 * beta * dbeta;
  dg_dx1(I12) = -dbeta;

  return;
}

void JMeanOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip, int imom,
                  const PhotonFrameState &s) {

  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  // energy and path length come from the frame this moment was enrolled with
  Real weight = pphot->wp[ip]*s.e*s.dl/MCConstants::c_cgs;
  pmcb->moments_user(imom,i3,i2,i1) += weight*pphot->acp[ip];
}

void AverageEnergy(MonteCarloBlock *pmcb, Photon *pphot, int ip, int imom,
                  const PhotonFrameState &s) {

  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  // energy and path length come from the frame this moment was enrolled with
  Real weight = pphot->wp[ip]*s.e*s.dl/MCConstants::c_cgs;
  pmcb->moments_user(imom,i3,i2,i1) += weight*s.e;
}

} // namespace
