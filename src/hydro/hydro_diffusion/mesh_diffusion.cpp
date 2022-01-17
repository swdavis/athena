//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2022 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh_diffusion.cpp
//! \brief functions to calculate mesh diffusion fluxes

// C/C++ headers

// Athena++ headers
#include "../../athena.hpp"        // X[123]DIR
#include "../../athena_arrays.hpp" // AthenaArray
#include "mesh_diffusion.hpp"      // MeshDiffusion

//----------------------------------------------------------------------------------------
//! \fn MeshDiffusion::MeshDiffusion(ParameterInput *pin)
//! \brief constructs an instance of MeshDiffusion

MeshDiffusion::MeshDiffusion(MeshBlock* pmb, ParameterInput* pin) :
    pmb(pmb), nu2mesh{pin->GetOrAddReal("hydro", "nu2mesh", 0.0)} {
}

//----------------------------------------------------------------------------------------
//! \fn void MeshDiffusion::AddFluxes
//! \brief is a wrapper to call individual functions for adding mesh-diffusion fluxes.

void MeshDiffusion::AddFluxes(
    const AthenaArray<Real>& cons, AthenaArray<Real>* flux) const {
  if (nu2mesh > 0.0) AddFluxHyper2(cons, flux);
}

//----------------------------------------------------------------------------------------
//! \fn void MeshDiffusion::AddFluxHyper2
//! \brief calculates and adds mesh hyper-diffusion flux of fourth order.

void MeshDiffusion::AddFluxHyper2(
    const AthenaArray<Real>& cons, AthenaArray<Real>* flux) const {
  const int nvar(cons.GetDim4());
  const int is(pmb->is), js(pmb->js), ks(pmb->ks);
  const int ie(pmb->ie), je(pmb->je), ke(pmb->ke);
  AthenaArray<Real> &x1flux(flux[X1DIR]);
  AthenaArray<Real> &x2flux(flux[X2DIR]);
  AthenaArray<Real> &x3flux(flux[X3DIR]);

  if (ie > is) {
    // Compute fluxes in X1 direction.
    for (int n = 0; n < nvar; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie + 1; ++i)
            x1flux(n,k,j,i) += nu2mesh * (
                cons(n,k,j,i+1) - cons(n,k,j,i-2) -
                3.0 * (cons(n,k,j,i) - cons(n,k,j,i-1)));
  }

  if (je > js) {
    // Compute fluxes in X2 direction.
    for (int n = 0; n < nvar; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je + 1; ++j)
          for (int i = is; i <= ie; ++i)
            x2flux(n,k,j,i) += nu2mesh * (
                cons(n,k,j+1,i) - cons(n,k,j-2,i) -
                3.0 * (cons(n,k,j,i) - cons(n,k,j-1,i)));
  }

  if (ke > ks) {
    // Compute fluxes in X3 direction.
    for (int n = 0; n < nvar; ++n)
      for (int k = ks; k <= ke + 1; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i)
            x3flux(n,k,j,i) += nu2mesh * (
                cons(n,k+1,j,i) - cons(n,k-2,j,i) -
                3.0 * (cons(n,k,j,i) - cons(n,k-1,j,i)));
  }
}
