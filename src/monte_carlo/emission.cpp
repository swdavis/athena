//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!  \file emission.cpp

// Athena++ headers
#include "montecarlo.hpp"
#include "../defs.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn void InitializeEmissionFreefree(MonteCarloBlock *pmbc, MeshBlock *pmb)
//  \brief Initialize emission array for static monte carlo calculation

void InitializeEmissionFreeFree(MonteCarloBlock *pmcb) {
  //void InitializeEmissionFreeFree(MonteCarloBlock *pmcb, MeshBlock *pmb) {

  MeshBlock *pmb = pmcb->pmy_block;

  Real heabund = 0.0; // Should have more general EOS functions
  Real mp = 1.6726e-24;
  Real eta0 = 1.032521e-11;
  Real g = 1.0;

  int il = pmb->is; int iu = pmb->ie;
  int jl = pmb->js; int ju = pmb->je;
  int kl = pmb->ks; int ku = pmb->ke;
  
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        Real temp = 1.e6;  // CHANGE
        Real nhii = pmb->phydro->u(IDN,k,j,i)/mp/(1.0+4.0*heabund);
        Real ne = (1.0+2.0*heabund) * nhii;
        Real vol = pmb->pcoord->GetCellVolume(k,j,i);
        pmcb->emission(k,j,i) = eta0/sqrt(temp)*ne*nhii*g*vol;
      }}}
}
