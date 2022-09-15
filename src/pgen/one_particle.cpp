//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file one_particle.cpp
//  \brief tests one particle.

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../particles/particles.hpp"

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Sets the initial conditions.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Get the (uniform) velocity of the gas.
  Real ux0, uy0, uz0;
  ux0 = pin->GetOrAddReal("problem", "ux0", 0.0);
  uy0 = pin->GetOrAddReal("problem", "uy0", 0.0);
  uz0 = pin->GetOrAddReal("problem", "uz0", 0.0);

  // Set a uniform, steady gas.
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = ux0;
        phydro->u(IM2,k,j,i) = uy0;
        phydro->u(IM3,k,j,i) = uz0;
        phydro->u(IEN,k,j,i) = 1.0 / (peos->GetGamma() - 1.0);
      }
    }
  }

  // Get the position and velocity of the particle.
  Real xp0, yp0, zp0, vpx0, vpy0, vpz0;
  xp0 = pin->GetOrAddReal("problem", "xp0", 0.0);
  yp0 = pin->GetOrAddReal("problem", "yp0", 0.0);
  zp0 = pin->GetOrAddReal("problem", "zp0", 0.0);
  vpx0 = pin->GetOrAddReal("problem", "vpx0", 0.0);
  vpy0 = pin->GetOrAddReal("problem", "vpy0", 0.0);
  vpz0 = pin->GetOrAddReal("problem", "vpz0", 0.0);

  // Check if the particle is in the meshblock.
  bool flag = true;
  if (block_size.nx1 > 1)
    flag = flag && block_size.x1min <= xp0 && xp0 < block_size.x1max;
  if (block_size.nx2 > 1)
    flag = flag && block_size.x2min <= yp0 && yp0 < block_size.x2max;
  if (block_size.nx3 > 1)
    flag = flag && block_size.x3min <= zp0 && zp0 < block_size.x3max;

  // Assign the particle, if any.
  if (flag) {
    ppar->Resize(1);
    ppar->xp[0] = xp0;
    ppar->yp[0] = yp0;
    ppar->zp[0] = zp0;
    ppar->vpx[0] = vpx0;
    ppar->vpy[0] = vpy0;
    ppar->vpz[0] = vpz0;
  }
}
