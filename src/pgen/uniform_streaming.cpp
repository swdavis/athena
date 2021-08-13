//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file uniform_streaming.cpp
//! \brief tests uniform streaming between gas and dust particles.

// C++ standard libraries
#include <cmath>    // round()
#include <sstream>  // stringstream

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../particles/dust_particles.hpp"

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Sets the initial conditions.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Sanity check.
  if (!PARTICLES) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [MeshBlock::ProblemGenerator]" << std::endl
        << "Dust particles were not enabled at configuration. " << std::endl;
    ATHENA_ERROR(msg);
  }

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
      }
    }
  }
  if (NON_BAROTROPIC_EOS) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          phydro->u(IEN,k,j,i) = 1.0 / (peos->GetGamma() - 1.0);
  }

  // Get the dust-to-gas ratio and the velocity of the particles.
  Real dtog, vpx0, vpy0, vpz0;
  dtog = pin->GetOrAddReal("problem", "dtog", 1.0);
  vpx0 = pin->GetOrAddReal("problem", "vpx0", 0.0);
  vpy0 = pin->GetOrAddReal("problem", "vpy0", 0.0);
  vpz0 = pin->GetOrAddReal("problem", "vpz0", 0.0);

  // Find the total number of particles in each direction.
  RegionSize& mesh_size = pmy_mesh->mesh_size;
  int npx1 = (block_size.nx1 > 1) ?
                  pin->GetOrAddInteger("problem", "npx1", mesh_size.nx1) : 1,
      npx2 = (block_size.nx2 > 1) ?
                  pin->GetOrAddInteger("problem", "npx2", mesh_size.nx2) : 1,
      npx3 = (block_size.nx3 > 1) ?
                  pin->GetOrAddInteger("problem", "npx3", mesh_size.nx3) : 1;

  // Find the mass of each particle and the distance between adjacent particles.
  Real vol = mesh_size.x1len * mesh_size.x2len * mesh_size.x3len;
  Real dx1 = mesh_size.x1len / npx1,
       dx2 = mesh_size.x2len / npx2,
       dx3 = mesh_size.x3len / npx3;
  DustParticles::SetOneParticleMass(dtog * vol / (npx1 * npx2 * npx3));

  // Determine number of particles in the block.
  int npx1_loc = static_cast<int>(std::round(block_size.x1len / dx1)),
      npx2_loc = static_cast<int>(std::round(block_size.x2len / dx2)),
      npx3_loc = static_cast<int>(std::round(block_size.x3len / dx3));
  int npar(npx1_loc * npx2_loc * npx3_loc);
  ppar->Resize(npar);

  // Assign the particles.
  int ipar = 0;
  for (int k = 0; k < npx3_loc; ++k) {
    Real zp1 = block_size.x3min + (k + 0.5) * dx3;
    for (int j = 0; j < npx2_loc; ++j) {
      Real yp1 = block_size.x2min + (j + 0.5) * dx2;
      for (int i = 0; i < npx1_loc; ++i) {
        Real xp1 = block_size.x1min + (i + 0.5) * dx1;
        ppar->xp[ipar] = xp1;
        ppar->yp[ipar] = yp1;
        ppar->zp[ipar] = zp1;
        ppar->vpx[ipar] = vpx0;
        ppar->vpy[ipar] = vpy0;
        ppar->vpz[ipar] = vpz0;
        ++ipar;
      }
    }
  }
}
