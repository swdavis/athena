//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarization.cpp
//! \brief supports moving a photon into and out of the polarized scattering basis
//
// See polarization.hpp for the convention these two functions implement.

// C++ headers
#include <cmath>
#include <complex>

// Athena++ headers
#include "../athena.hpp"
#include "polarization.hpp"
#include "montecarlo.hpp"
#include "photon.hpp"
#include "tetrad.hpp"
#include "mccoord.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ToScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//! \brief rotate the photon into the basis the polarized scattering routines assume
//
// On a spherical grid the wavevector is stored in the local orthonormal spherical basis,
// which rotates along the ray; the scattering routines want cartesian components.  Q and
// U are deliberately untouched: they are already referenced to the global cartesian
// meridian, which does not move. Used by non GR coordinates.


static bool MeridianPair(const Real n[3], Real lhat[3], Real rhat[3]);
static void WriteMeridianStokes(Photon *pphot, int ip, Real ecov[4][4],
                                const Real lhat[3], const Real rhat[3]);

void ToScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  if (pmcb->topology != MCTOPO_SPHERICAL) return;

  Real cth = std::cos(pphot->x2p[ip]);
  Real sth = std::sin(pphot->x2p[ip]);
  Real cph = std::cos(pphot->x3p[ip]);
  Real sph = std::sin(pphot->x3p[ip]);

  Real k3[3];
  k3[0] = pphot->k1p[ip];
  k3[1] = pphot->k2p[ip];
  k3[2] = pphot->k3p[ip];

  pphot->k1p[ip] = k3[0]*sth*cph + k3[1]*cth*cph - k3[2]*sph;
  pphot->k2p[ip] = k3[0]*sth*sph + k3[1]*cth*sph + k3[2]*cph;
  pphot->k3p[ip] = k3[0]*cth     - k3[1]*sth;

  Real knorm = std::sqrt(SQR(pphot->k1p[ip]) + SQR(pphot->k2p[ip]) + SQR(pphot->k3p[ip]));
  pphot->k1p[ip] /= knorm;
  pphot->k2p[ip] /= knorm;
  pphot->k3p[ip] /= knorm;
}

//----------------------------------------------------------------------------------------
//! \fn void FromScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//! \brief rotate the photon back out of the scattering basis
//
// Exact inverse of ToScatteringBasis, evaluated at the same position, so the round trip
// is the identity to round-off. Used by non GR coordinates.

void FromScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  if (pmcb->topology != MCTOPO_SPHERICAL) return;

  Real cth = std::cos(pphot->x2p[ip]);
  Real sth = std::sin(pphot->x2p[ip]);
  Real cph = std::cos(pphot->x3p[ip]);
  Real sph = std::sin(pphot->x3p[ip]);

  Real k3[3];
  k3[0] = pphot->k1p[ip];
  k3[1] = pphot->k2p[ip];
  k3[2] = pphot->k3p[ip];

  pphot->k1p[ip] =  k3[0]*sth*cph + k3[1]*sth*sph + k3[2]*cth;
  pphot->k2p[ip] =  k3[0]*cth*cph + k3[1]*cth*sph - k3[2]*sth;
  pphot->k3p[ip] = -k3[0]*sph     + k3[1]*cph;

  Real knorm = std::sqrt(SQR(pphot->k1p[ip]) + SQR(pphot->k2p[ip]) + SQR(pphot->k3p[ip]));
  pphot->k1p[ip] /= knorm;
  pphot->k2p[ip] /= knorm;
  pphot->k3p[ip] /= knorm;
}

//----------------------------------------------------------------------------------------
//! \fn static bool MeridianBasis(...)
//! \brief the scattering basis in the comoving orthonormal frame
//
// Builds the comoving tetrad E -- the same construction TransformToComoving uses, so the
// photon's stored k is expressed in it -- and returns the pair of unit three-vectors the
// scattering routines reference Q and U to:
//
//   n = k-hat in E,   l = normalize(z - (z.n) n),   r = n x l,
//
// with z the third spatial leg of E.  That is the meridian plane containing the frame's
// z axis and the photon direction, which is the convention ScatterThomsonPolarized reads
// off the stored components (see polarization.hpp).  l goes into e_(1) because
// StokesToTensor puts I+Q there, i.e. Q > 0 means polarization along the meridian.
//
// Returns false when n is parallel to z, where the meridian is undefined; the caller
// then leaves the Stokes parameters alone

static bool MeridianBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip,
                          Real econ[4][4], Real ecov[4][4], Real lhat[3], Real rhat[3]) {

  int i1 = pphot->i1p[ip], i2 = pphot->i2p[ip], i3 = pphot->i3p[ip];

  Real x[4];
  x[IMC0] = pphot->x0p[ip];
  x[IMC1] = pphot->x1p[ip];
  x[IMC2] = pphot->x2p[ip];
  x[IMC3] = pphot->x3p[ip];
  Real gcov[4][4];
  pmcb->pcoord->Metric(x, gcov);

  // The frame k is expressed in at a scattering.  The MonteCarloBlock constructor
  // guarantees vel is allocated whenever polarized is set, holding the fluid velocity,
  // the normal observer or the coordinate observer as appropriate, so no test is needed
  // here.
  Real ucon[4];
  for (int m = 0; m < 4; ++m) ucon[m] = pmcb->vel(i3, i2, i1, m);
  ConstructTetrad(ucon, gcov, econ, ecov);

  // k as stored is already in this frame and is a unit direction there
  Real n[3];
  n[0] = pphot->k1p[ip];
  n[1] = pphot->k2p[ip];
  n[2] = pphot->k3p[ip];
  Real nmag = std::sqrt(SQR(n[0]) + SQR(n[1]) + SQR(n[2]));
  if (nmag <= TINY_NUMBER) return false;
  for (int i = 0; i < 3; ++i) n[i] /= nmag;

  return MeridianPair(n, lhat, rhat);
}

//----------------------------------------------------------------------------------------
//! \fn static bool MeridianPair(const Real n[3], Real lhat[3], Real rhat[3])
//! \brief the meridian pair for a unit direction in some orthonormal frame
//
//   l = normalize(z - (z.n) n),   r = n x l,
//
// with z the frame's third spatial leg.  Split out so that the frame can be chosen by the
// caller: the scattering routines want the comoving frame, the outputs want the normal
// observer.  Returns false when n is parallel to z, where the meridian is undefined.

static bool MeridianPair(const Real n[3], Real lhat[3], Real rhat[3]) {

  Real l[3] = {-n[2]*n[0], -n[2]*n[1], 1.0 - n[2]*n[2]};
  Real lmag = std::sqrt(SQR(l[0]) + SQR(l[1]) + SQR(l[2]));
  if (lmag <= TINY_NUMBER) return false;          // photon along the frame z axis
  for (int i = 0; i < 3; ++i) lhat[i] = l[i]/lmag;

  // r = n x l, completing a right-handed triad (l, r, n)
  rhat[0] = n[1]*lhat[2] - n[2]*lhat[1];
  rhat[1] = n[2]*lhat[0] - n[0]*lhat[2];
  rhat[2] = n[0]*lhat[1] - n[1]*lhat[0];

  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void CoherencyToScatteringStokes(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//! \brief coherency tensor -> Stokes, in the scattering basis

void CoherencyToScatteringStokes(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  if (!pmcb->pmy_mc->general_pusher_flag) return;

  Real econ[4][4], ecov[4][4], lhat[3], rhat[3];
  if (!MeridianBasis(pmcb, pphot, ip, econ, ecov, lhat, rhat)) return;

  WriteMeridianStokes(pphot, ip, ecov, lhat, rhat);
}

//----------------------------------------------------------------------------------------
//! \fn static void WriteMeridianStokes(...)
//! \brief project the coherency tensor onto a meridian pair and store the Stokes
//!        parameters the scattering routines and the outputs both read

static void WriteMeridianStokes(Photon *pphot, int ip, Real ecov[4][4],
                                const Real lhat[3], const Real rhat[3]) {

  std::complex<Real> ntet[4][4];
  pphot->PolarizationToTetrad(ntet, ecov, ip);

  // Contract the spatial block onto the meridian pair
  std::complex<Real> nll(0., 0.), nrr(0., 0.), nlr(0., 0.), nrl(0., 0.);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const std::complex<Real> &nij = ntet[i+1][j+1];
      nll += lhat[i]*lhat[j]*nij;
      nrr += rhat[i]*rhat[j]*nij;
      nlr += lhat[i]*rhat[j]*nij;
      nrl += rhat[i]*lhat[j]*nij;
    }
  }

  // Pass the meridian pair to TensorToStokes the(l,r) contraction is the (1,2) block it
  // reads, with l and r standing in for e_(1) and e_(2).
  Real inten = 0.5*(nll + nrr).real();
  if (std::fabs(inten) <= TINY_NUMBER) return;   // TensorToStokes divides by I

  std::complex<Real> nmer[4][4];
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) nmer[a][b] = std::complex<Real>(0., 0.);
  nmer[IMC1][IMC1] = nll;
  nmer[IMC1][IMC2] = nlr;
  nmer[IMC2][IMC1] = nrl;
  nmer[IMC2][IMC2] = nrr;

  Real stokes[4];
  TensorToStokes(nmer, stokes);  // returns I normalised to one

  pphot->sip[ip] = stokes[0];
  pphot->sqp[ip] = stokes[1];
  pphot->sup[ip] = stokes[2];
  pphot->svp[ip] = stokes[3];
}

//----------------------------------------------------------------------------------------
//! \fn void CoherencyToObserverStokes(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//! \brief refresh the stored Stokes parameters from the coherency tensor, in the frame
//!        the outputs are referenced to
//
// Sets sqp/sup/svp to the Stokes parameters referenced to the meridian containing the
// normal observer's z axis and the photon's direction, which is what the outputs write.
// The frame is the normal (Eulerian) observer at the photon's position, n^mu = -alpha
// g^{mu t}, rather than the fluid. Q is positive along the meridian containing the
// frame's third spatial leg

void CoherencyToObserverStokes(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  // The legacy pushers reference the Stokes parameters to a globally constant frame and
  // do not evolve them between scatterings, so they are already current.
  if (!pmcb->pmy_mc->general_pusher_flag) return;

  Real x[4];
  x[IMC0] = pphot->x0p[ip];
  x[IMC1] = pphot->x1p[ip];
  x[IMC2] = pphot->x2p[ip];
  x[IMC3] = pphot->x3p[ip];
  Real gcov[4][4], gcon[4][4];
  pmcb->pcoord->Metric(x, gcov);
  pmcb->pcoord->InverseMetric(x, gcon);

  // normal observer: alpha = 1/sqrt(-g^tt), n^mu = -alpha g^{mu t}
  if (gcon[IMC0][IMC0] >= 0.) return;
  Real alpha = 1.0/std::sqrt(-gcon[IMC0][IMC0]);
  Real ncon[4];
  for (int m = 0; m < 4; ++m) ncon[m] = -alpha*gcon[m][IMC0];

  Real econ[4][4], ecov[4][4];
  ConstructTetrad(ncon, gcov, econ, ecov);

  // The photon is in the coordinate frame here, so its direction has to be projected into
  // the observer's frame before the meridian can be built on it.
  Real kcoord[4], ktet[4];
  pphot->GetFourVector(ip, false, kcoord);
  CoordinateToTetrad(kcoord, ktet, ecov);
  Real nmag = std::sqrt(SQR(ktet[IMC1]) + SQR(ktet[IMC2]) + SQR(ktet[IMC3]));
  if (nmag <= TINY_NUMBER) return;
  Real n[3] = {ktet[IMC1]/nmag, ktet[IMC2]/nmag, ktet[IMC3]/nmag};

  Real lhat[3], rhat[3];
  if (!MeridianPair(n, lhat, rhat)) return;   // photon along the frame z axis

  WriteMeridianStokes(pphot, ip, ecov, lhat, rhat);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//! \fn void ScatteringStokesToCoherency(MonteCarloBlock *pmcb, Photon *pphot, int ip)
//! \brief Stokes -> coherency tensor, undoing CoherencyToScatteringStokes
//
// Rebuilds the basis at the photon's current direction, which after a scattering is the
// outgoing one, since the outgoing Stokes parameters are referenced to the outgoing
// meridian.

void ScatteringStokesToCoherency(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  if (!pmcb->pmy_mc->general_pusher_flag) return;

  Real econ[4][4], ecov[4][4], lhat[3], rhat[3];
  if (!MeridianBasis(pmcb, pphot, ip, econ, ecov, lhat, rhat)) return;

  Real inten = pphot->sip[ip];
  if (std::fabs(inten) <= TINY_NUMBER) return;
  Real q = pphot->sqp[ip]/inten;
  Real u = pphot->sup[ip]/inten;
  Real v = pphot->svp[ip]/inten;

  // StokesToTensor retruns the four block entries, then spread them over the meridian
  // pair: N = N11 l l + N12 l r + N21 r l + N22 r r.
  Real stokes[4] = {1.0, q, u, v};
  std::complex<Real> blk[4][4];
  StokesToTensor(stokes, blk);

  std::complex<Real> ntet[4][4];
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) ntet[a][b] = std::complex<Real>(0., 0.);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ntet[i+1][j+1] = blk[IMC1][IMC1]*lhat[i]*lhat[j] + blk[IMC1][IMC2]*lhat[i]*rhat[j]
                     + blk[IMC2][IMC1]*rhat[i]*lhat[j] + blk[IMC2][IMC2]*rhat[i]*rhat[j];
    }
  }

  pphot->PolarizationToCoord(ntet, econ, ip);
}
