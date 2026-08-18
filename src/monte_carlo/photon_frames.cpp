//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photon_frames.cpp
//! \brief projecting photons and radiation moments between reference frames.
//!
//! Three frames are in play.  MCFRAME_LAB and MCFRAME_COMOVING are orthonormal: in general
//! relativity the tetrads of the normal (Eulerian) observer and of the frame vel is built
//! on, and in flat spacetime the Eulerian frame and its Lorentz boost.  MCFRAME_COORD is
//! the raw coordinate basis.
//!
//! The geometric primitives these are built on -- ConstructTetrad and the vector algebra
//! that supports it -- are the layer below, in tetrad.cpp.  Nothing there knows about
//! photons or mesh blocks; everything here does.
//!
//! Two operations live here and are easy to conflate.  PhotonFrames projects a single
//! photon, per crossing; ComovingFrameMatrix and DeriveComovingMoments transform an
//! accumulated tensor, once per zone at output.  The second is exact -- the transform is
//! the same matrix for every photon in the zone -- which is what lets the comoving moments
//! be derived rather than accumulated.

// C headers

// C++ headers
#include <cmath>
#include <cstring>   // strcmp

// Athena++ headers
#include "photon_frames.hpp"
#include "montecarlo.hpp"
#include "tetrad.hpp"
#include "photon.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

//----------------------------------------------------------------------------------------
//! \fn PhotonFrames::PhotonFrames(...)
//! \brief bind to one photon; frames are projected lazily as they are asked for

PhotonFrames::PhotonFrames(MonteCarloBlock *pmcb, Photon *pphot, int ip, Real dl)
    : pmcb_(pmcb), pphot_(pphot), ip_(ip), dl_(dl) {
  general_ = pmcb->pmy_mc->general_pusher_flag;
  gr_tetrad_ = general_ && GENERAL_RELATIVITY;
  if (general_) pphot->GetFourVector(ip, false, kco_);
  for (int f=0; f<MCFRAME_N; ++f) done_[f] = false;
}

//----------------------------------------------------------------------------------------
//! \fn void PhotonFrames::Fill(MCFrame f)
//! \brief project the photon into one frame

void PhotonFrames::Fill(MCFrame f) {
  PhotonFrameState &s = st_[f];
  const Real ep = pphot_->ep[ip_];
  const int i1 = pphot_->i1p[ip_], i2 = pphot_->i2p[ip_], i3 = pphot_->i3p[ip_];

  if (f == MCFRAME_COORD) {
    // Not orthonormal: n holds k^i/k^0 and is deliberately not normalised.
    s.e = kco_[IMC0];
    for (int i=0; i<3; ++i) s.n[i] = kco_[IMC1+i]/ep;
    s.dl = dl_;
    return;
  }

  if (gr_tetrad_) {
    // One matrix multiply carries the coordinate four-vector to orthonormal components.
    const AthenaArray<Real> &m = (f == MCFRAME_LAB) ? pmcb_->boost_lab : pmcb_->boost_cmv;
    Real p[4];
    for (int a=0; a<4; ++a) {
      p[a] = 0.;
      for (int b=0; b<4; ++b) p[a] += m(i3,i2,i1,a,b) * kco_[b];
    }
    // k is null at the photon to machine precision, but the tetrad is orthonormal with
    // respect to the metric at the zone centre, so the projection is not exactly null.
    // Normalising by the spatial magnitude keeps n a genuine direction and the moment
    // tensor exactly traceless; the energy comes from the time component.
    Real mag = std::sqrt(SQR(p[IMC1]) + SQR(p[IMC2]) + SQR(p[IMC3]));
    s.e = p[IMC0];
    for (int i=0; i<3; ++i) s.n[i] = p[IMC1+i]/mag;
    s.dl = dl_ * s.e / ep;
    return;
  }

  // ---- flat spacetime ----
  if (f == MCFRAME_LAB) {
    if (general_) {
      Real x[4] = {pphot_->x0p[ip_], pphot_->x1p[ip_], pphot_->x2p[ip_], pphot_->x3p[ip_]};
      Real invtet[4][4], kf[4];
      pmcb_->pcoord->InverseTetrad(x, invtet);
      for (int a=0; a<4; ++a) {
        kf[a] = 0.;
        for (int b=0; b<4; ++b) kf[a] += invtet[a][b] * kco_[b];
      }
      // The tetrad time leg is unity for every flat coordinate system, so kf[0] == ep and
      // dividing the spatial parts by ep is the same unit direction the legacy pushers
      // store.  Kept in this form rather than normalised by the spatial magnitude so the
      // arithmetic is unchanged from before this was factored out.
      s.e = kf[IMC0];
      for (int i=0; i<3; ++i) s.n[i] = kf[IMC1+i]/ep;
    } else {
      s.e = ep;
      s.n[0] = pphot_->k1p[ip_];
      s.n[1] = pphot_->k2p[ip_];
      s.n[2] = pphot_->k3p[ip_];
      if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
        // Re-express the direction in the orthonormal basis at the zone centre, which is
        // the basis the moments are accumulated in.
        Real sth = std::sin(pphot_->x2p[ip_]), cth = std::cos(pphot_->x2p[ip_]);
        Real sph = std::sin(pphot_->x3p[ip_]), cph = std::cos(pphot_->x3p[ip_]);
        Real nx = sth*cph*s.n[0] + cth*cph*s.n[1] - sph*s.n[2];
        Real ny = sth*sph*s.n[0] + cth*sph*s.n[1] + cph*s.n[2];
        Real nz = cth*s.n[0] - sth*s.n[1];
        sth = std::sin(pmcb_->pmy_block->pcoord->x2v(i2));
        cth = std::cos(pmcb_->pmy_block->pcoord->x2v(i2));
        sph = std::sin(pmcb_->pmy_block->pcoord->x3v(i3));
        cph = std::cos(pmcb_->pmy_block->pcoord->x3v(i3));
        s.n[0] = sth*cph*nx + sth*sph*ny + cth*nz;
        s.n[1] = cth*cph*nx + cth*sph*ny - sth*nz;
        s.n[2] = -sph*nx + cph*ny;
      }
    }
    s.dl = dl_;
    return;
  }

  // comoving in flat spacetime: boost the lab direction, which is where boost_cmv acts
  const PhotonFrameState &lab = Get(MCFRAME_LAB);
  Real ki[4] = {1., lab.n[0], lab.n[1], lab.n[2]};
  Real kc[4];
  for (int a=0; a<4; ++a) {
    kc[a] = 0.;
    for (int b=0; b<4; ++b) kc[a] += pmcb_->boost_cmv(i3,i2,i1,a,b) * ki[b];
  }
  Real shift = kc[IMC0];
  s.e = ep * shift;
  for (int i=0; i<3; ++i) s.n[i] = kc[IMC1+i]/kc[IMC0];
  s.dl = dl_ * shift;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::ComovingFrameMatrix(...)
//! \brief the lab -> comoving transformation for one zone
//!
//! Both frames are orthonormal at the same event, so this is a Lorentz transformation.
//! In general relativity it must be composed from the two tetrads rather than rebuilt as
//! a pure boost from beta: the Gram-Schmidt legs are grown from the same coordinate trial
//! vectors against different timelike directions, so the map is a boost composed with a
//! rotation whenever the fluid velocity has all three spatial components.  A pure boost
//! would leave Er correct, being a scalar, while silently rotating flux and pressure.

void MonteCarloBlock::ComovingFrameMatrix(int k, int j, int i, const AthenaArray<Real> &g,
                                          const AthenaArray<Real> &gi, Real lam[4][4]) {
  if (!GENERAL_RELATIVITY) {
    // boost_cmv already maps lab-frame components to the comoving frame.
    for (int a=0; a<4; ++a)
      for (int b=0; b<4; ++b) lam[a][b] = boost_cmv(k,j,i,a,b);
    return;
  }

  Real x[4];
  x[IMC0] = 0.;
  x[IMC1] = pmy_block->pcoord->x1v(i);
  x[IMC2] = pmy_block->pcoord->x2v(j);
  x[IMC3] = pmy_block->pcoord->x3v(k);
  Real gcov[4][4];
  pcoord->Metric(x, gcov);

  Real alpha = 1.0/std::sqrt(-gi(I00,i));
  Real ncon[4];
  ncon[IMC0] = -alpha*gi(I00,i);
  ncon[IMC1] = -alpha*gi(I01,i);
  ncon[IMC2] = -alpha*gi(I02,i);
  ncon[IMC3] = -alpha*gi(I03,i);
  Real econL[4][4], ecovL[4][4];
  ConstructTetrad(ncon, gcov, econL, ecovL);

  Real ucon[4];
  for (int m=0; m<4; ++m) ucon[m] = vel(k,j,i,m);
  Real econF[4][4], ecovF[4][4];
  ConstructTetrad(ucon, gcov, econF, ecovF);

  for (int a=0; a<4; ++a) {
    for (int b=0; b<4; ++b) {
      lam[a][b] = 0.;
      for (int m=0; m<4; ++m) lam[a][b] += ecovF[a][m] * econL[b][m];
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::DeriveComovingMoments()
//! \brief fill moments_com by transforming the accumulated lab moments
//!
//! The moments are a rank two tensor and the transformation is the same matrix for every
//! photon in the zone, so sum Lp Lp == L (sum p p) L exactly.  Deriving is therefore not
//! an approximation of accumulating, and it removes the second per-photon projection --
//! measured at 5% of runtime on the general pusher and 12% on the legacy pusher.

void MonteCarloBlock::DeriveComovingMoments() {
  const Real c_cgs = MCConstants::c_cgs;
  AthenaArray<Real> g, gi;
  if (GENERAL_RELATIVITY) {
    g.NewAthenaArray(NMETRIC,ie+1);
    gi.NewAthenaArray(NMETRIC,ie+1);
  }
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      if (GENERAL_RELATIVITY) pmy_block->pcoord->CellMetric(k,j,is,ie,g,gi);
      for (int i=is; i<=ie; ++i) {
        Real lam[4][4];
        ComovingFrameMatrix(k,j,i,g,gi,lam);
        for (int m=0; m<pmy_mc->ntype; ++m) {
          // MomentSlot carries the stored convention, shared with AccumulateMoments
          Real T[4][4];
          for (int a=0; a<4; ++a)
            for (int b=0; b<4; ++b) {
              T[a][b] = moments(m,MomentSlot[a][b],k,j,i);
              if (MomentNeedsC(a,b)) T[a][b] /= c_cgs;
            }

          Real T1[4][4], U[4][4];
          for (int a=0; a<4; ++a)
            for (int b=0; b<4; ++b) {
              T1[a][b] = 0.;
              for (int c=0; c<4; ++c) T1[a][b] += lam[a][c]*T[c][b];
            }
          for (int a=0; a<4; ++a)
            for (int b=0; b<4; ++b) {
              U[a][b] = 0.;
              for (int c=0; c<4; ++c) U[a][b] += T1[a][c]*lam[b][c];
            }

          for (int a=0; a<4; ++a)
            for (int b=a; b<4; ++b) {
              Real v = U[a][b];
              if (MomentNeedsC(a,b)) v *= c_cgs;
              moments_com(m,MomentSlot[a][b],k,j,i) = v;
            }
        }
      }
    }
  }
  if (GENERAL_RELATIVITY) { g.DeleteAthenaArray(); gi.DeleteAthenaArray(); }
}
