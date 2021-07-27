//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file dust_particles.cpp
//! \brief implements functions in the DustParticles class

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // nan()

// Athena++ headers
#include "../athena.hpp"
#include "../coordinates/coordinates.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
#include "particle_gravity.hpp"
#include "particles.hpp"

// Class variable initialization
bool DustParticles::initialized(false);
bool DustParticles::backreaction(false);
bool DustParticles::dragforce(true);
bool DustParticles::variable_taus(false);

int DustParticles::iwx = -1, DustParticles::iwy = -1, DustParticles::iwz = -1;
int DustParticles::idpx1 = -1, DustParticles::idpx2 = -1, DustParticles::idpx3 = -1;
int DustParticles::itaus = -1;

Real DustParticles::mass = 1.0, DustParticles::taus0 = 0.0;

//--------------------------------------------------------------------------------------
//! \fn void Particles::FindDensityOnMesh(Mesh *pm, bool include_momentum)
//! \brief finds the mass density of particles on the mesh.  If include_momentum is
//!   true, the momentum density field is also included.
//!
//! \note
//!   Postcondition:
//!   ppm->weight becomes the density in each cell, and if include_momentum
//!   is true, ppm->meshaux(imom1:imom3,:,:,:) becomes the momentum density.

void DustParticles::FindDensityOnMesh(Mesh *pm, bool include_momentum) {
  // Assign the particles onto the mesh.
  Particles::FindDensityOnMesh(pm, include_momentum);

  for (int b = 0; b < pm->nblocal; ++b) {
    ParticleMesh *ppm(pm->my_blocks(b)->ppar->ppm);

    // Find the mass density.
    for (int k = ppm->ks; k <= ppm->ke; ++k)
      for (int j = ppm->js; j <= ppm->je; ++j)
        for (int i = ppm->is; i <= ppm->ie; ++i)
          ppm->weight(k,j,i) *= mass;

    // Find the momentum density.
    if (include_momentum) {
      for (int k = ppm->ks; k <= ppm->ke; ++k)
        for (int j = ppm->js; j <= ppm->je; ++j)
          for (int i = ppm->is; i <= ppm->ie; ++i) {
            ppm->meshaux(imom1,k,j,i) *= mass;
            ppm->meshaux(imom2,k,j,i) *= mass;
            ppm->meshaux(imom3,k,j,i) *= mass;
          }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::Initialize(Mesh *pm, ParameterInput *pin)
//! \brief initializes the class.

void DustParticles::Initialize(Mesh *pm, ParameterInput *pin) {
  // Initialize first the parent class.
  Particles::Initialize(pm, pin);

  if (!initialized) {
    // Add working array at particles for gas velocity/particle momentum change.
    iwx = AddWorkingArray();
    iwy = AddWorkingArray();
    iwz = AddWorkingArray();

    // Define mass.
    mass = pin->GetOrAddReal("particles", "mass", 1.0);

    // Define stopping time.
    variable_taus = pin->GetOrAddBoolean("particles", "variable_taus", variable_taus);
    taus0 = pin->GetOrAddReal("particles", "taus0", taus0);
    if (variable_taus) itaus = AddAuxProperty();

    // Turn on/off back reaction.
    dragforce = taus0 >= 0.0;
    backreaction = pin->GetOrAddBoolean("particles", "backreaction", false);
    if (taus0 == 0.0) backreaction = false;

    if (backreaction) {
      idpx1 = imom1;
      idpx2 = imom2;
      idpx3 = imom3;
    }

    if (SELF_GRAVITY_ENABLED == 2)
      ParticleGravity::Initialize();

    initialized = true;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::SetOneParticleMass(Real new_mass)
//! \brief sets the mass of each particle.

void DustParticles::SetOneParticleMass(Real new_mass) {
  pinput->SetReal("particles", "mass", mass = new_mass);
}

//--------------------------------------------------------------------------------------
//! \fn DustParticles::DustParticles(MeshBlock *pmb, ParameterInput *pin)
//! \brief constructs a DustParticles instance.

static std::vector<Real> dummy_vector(1, std::nan(NULL));

DustParticles::DustParticles(MeshBlock *pmb, ParameterInput *pin)
  : Particles(pmb, pin),
    wx(work[iwx]), wy(work[iwy]), wz(work[iwz]),
    taus(variable_taus ? aux[itaus] : dummy_vector) {
  if (backreaction) {
    // Assign shallow copies for momentum feedback.
    dpx1.InitWithShallowSlice(ppm->meshaux, 4, idpx1, 1);
    dpx2.InitWithShallowSlice(ppm->meshaux, 4, idpx2, 1);
    dpx3.InitWithShallowSlice(ppm->meshaux, 4, idpx3, 1);
  }

  if (SELF_GRAVITY_ENABLED == 2)
    // Activate particle gravity.
    ppgrav = new ParticleGravity(this);
}

//--------------------------------------------------------------------------------------
//! \fn DustParticles::~DustParticles()
//! \brief destroys a DustParticles instance.

DustParticles::~DustParticles() {
  if (backreaction) {
    dpx1.DeleteAthenaArray();
    dpx2.DeleteAthenaArray();
    dpx3.DeleteAthenaArray();
  }

  if (SELF_GRAVITY_ENABLED == 2)
    delete ppgrav;
}

//--------------------------------------------------------------------------------------
//! \fn AthenaArray<Real> DustParticles::GetVelocityField()
//! \brief returns the particle velocity on the mesh.
//!
//! \note
//!   Precondition:
//!   The particle properties on mesh must be assigned using the class method
//!   DustParticles::FindDensityOnMesh().

AthenaArray<Real> DustParticles::GetVelocityField() const {
  AthenaArray<Real> vel(3, ppm->nx3_, ppm->nx2_, ppm->nx1_);
  for (int k = ppm->ks; k <= ppm->ke; ++k)
    for (int j = ppm->js; j <= ppm->je; ++j)
      for (int i = ppm->is; i <= ppm->ie; ++i) {
        Real rho(ppm->weight(k,j,i));
        rho = (rho > 0.0) ? rho : 1.0;
        vel(0,k,j,i) = ppm->meshaux(imom1,k,j,i) / rho;
        vel(1,k,j,i) = ppm->meshaux(imom2,k,j,i) / rho;
        vel(2,k,j,i) = ppm->meshaux(imom3,k,j,i) / rho;
      }
  return vel;
}

//--------------------------------------------------------------------------------------
//! \fn Real DustParticles::NewBlockTimeStep();
//! \brief returns the time step required by particles in the block.

Real DustParticles::NewBlockTimeStep() {
  // Run first the parent class.
  Real dt = Particles::NewBlockTimeStep();

  // Nothing to do for tracer particles.
  if (taus0 <= 0.0) return dt;

  Real epsmax = 0;
  if (backreaction) {
    // Find the maximum local solid-to-gas density ratio.
    Coordinates *pc = pmy_block->pcoord;
    Hydro *phydro = pmy_block->phydro;
    const int is = ppm->is, js = ppm->js, ks = ppm->ks;
    const int ie = ppm->ie, je = ppm->je, ke = ppm->ke;
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          Real epsilon = ppm->weight(k,j,i) / (
                         pc->GetCellVolume(k,j,i) * phydro->u(IDN,k,j,i));
          epsmax = std::max(epsmax, epsilon);
        }
    epsmax *= mass;
  }

  // Return the drag timescale.
  return std::min(dt, static_cast<Real>(cfl_par * taus0 / (1.0 + epsmax)));
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::SourceTerms()
//! \brief adds acceleration to particles.

void DustParticles::SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  if (dragforce) {
    // Interpolate gas velocity onto particles.
    ppm->InterpolateMeshToParticles(meshsrc, IVX, work, iwx, 3);

    // Transform the gas velocity into Cartesian.
    const Coordinates *pc = pmy_block->pcoord;
    for (int k = 0; k < npar; ++k) {
      Real x1, x2, x3;
      //! \todo (ccyang):
      //! - using (xp0, yp0, zp0) is a temporary hack.
      pc->CartesianToMeshCoords(xp0[k], yp0[k], zp0[k], x1, x2, x3);
      pc->MeshCoordsToCartesianVector(x1, x2, x3, wx[k], wy[k], wz[k],
                                                  wx[k], wy[k], wz[k]);
    }

    // Add drag force to particles.
    if (variable_taus) {
      // Variable stopping time
      UserStoppingTime(t, dt, meshsrc);
      for (int k = 0; k < npar; ++k) {
        //! \todo (ccyang):
        //! - This is a temporary hack; to be fixed.
        Real tmpx(vpx[k]), tmpy(vpy[k]), tmpz(vpz[k]);
        //
        Real c = dt / taus[k];
        wx[k] = c * (vpx[k] - wx[k]);
        wy[k] = c * (vpy[k] - wy[k]);
        wz[k] = c * (vpz[k] - wz[k]);
        vpx[k] = vpx0[k] - wx[k];
        vpy[k] = vpy0[k] - wy[k];
        vpz[k] = vpz0[k] - wz[k];
        //
        vpx0[k] = tmpx; vpy0[k] = tmpy; vpz0[k] = tmpz;
        //
      }
    } else if (taus0 > 0.0) {
      // Constant stopping time
      Real c = dt / taus0;
      for (int k = 0; k < npar; ++k) {
        //! \todo (ccyang):
        //! - This is a temporary hack; to be fixed.
        Real tmpx(vpx[k]), tmpy(vpy[k]), tmpz(vpz[k]);
        //
        wx[k] = c * (vpx[k] - wx[k]);
        wy[k] = c * (vpy[k] - wy[k]);
        wz[k] = c * (vpz[k] - wz[k]);
        vpx[k] = vpx0[k] - wx[k];
        vpy[k] = vpy0[k] - wy[k];
        vpz[k] = vpz0[k] - wz[k];
        //
        vpx0[k] = tmpx; vpy0[k] = tmpy; vpz0[k] = tmpz;
        //
      }
    } else if (taus0 == 0.0) {
      // Tracer particles
      for (int k = 0; k < npar; ++k) {
        vpx[k] = wx[k];
        vpy[k] = wy[k];
        vpz[k] = wz[k];
      }
    }
  } else {
    for (int k = 0; k < npar; ++k) {
      //! \todo (ccyang):
      //! - This is a temporary hack; to be fixed.
      Real tmpx(vpx[k]), tmpy(vpy[k]), tmpz(vpz[k]);
      vpx[k] = vpx0[k]; vpy[k] = vpy0[k]; vpz[k] = vpz0[k];
      vpx0[k] = tmpx; vpy0[k] = tmpy; vpz0[k] = tmpz;
    }
  }

  if (SELF_GRAVITY_ENABLED == 2) {
    // Add gravitational force from the Poisson solution.
    ppgrav->FindGravitationalForce(pmy_block->pgrav->phi);
    ppgrav->ExertGravitationalForce(dt);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::UserSourceTerms(Real t, Real dt,
//!                                         const AthenaArray<Real>& meshsrc)
//! \brief adds additional source terms to particles, overloaded by the user.

void __attribute__((weak)) DustParticles::UserSourceTerms(
    Real t, Real dt, const AthenaArray<Real>& meshsrc) {
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::UserStoppingTime(Real t, Real dt,
//!                                          const AthenaArray<Real>& meshsrc)
//! \brief assigns time-dependent stopping time to each particle, overloaded by the user.

void __attribute__((weak)) DustParticles::UserStoppingTime(
    Real t, Real dt, const AthenaArray<Real>& meshsrc) {
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::ReactToMeshAux(
//!              Real t, Real dt, const AthenaArray<Real>& meshsrc)
//! \brief Reacts to meshaux before boundary communications.

void DustParticles::ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  // Nothing to do if no back reaction.
  if (!dragforce || !backreaction) return;

  // Transform the momentum change in mesh coordinates.
  const Coordinates *pc = pmy_block->pcoord;
  for (int k = 0; k < npar; ++k)
    //! \todo (ccyang):
    //! - using (xp0, yp0, zp0) is a temporary hack.
    pc->CartesianToMeshCoordsVector(xp0[k], yp0[k], zp0[k],
        mass * wx[k], mass * wy[k], mass * wz[k], wx[k], wy[k], wz[k]);

  // Assign the momentum change onto mesh.
  ppm->AssignParticlesToMeshAux(work, iwx, idpx1, 3);
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::DepositToMesh(Real t, Real dt,
//!              const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst);
//! \brief Deposits meshaux to Mesh.

void DustParticles::DepositToMesh(
         Real t, Real dt, const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst) {
  if (dragforce && backreaction)
    // Deposit particle momentum changes to the gas.
    ppm->DepositMeshAux(meshdst, idpx1, IM1, 3);
}
