#ifndef PARTICLES_DUST_PARTICLES_HPP_
#define PARTICLES_DUST_PARTICLES_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file dust_particles.hpp
//! \brief defines class DustParticles.
//======================================================================================

// C/C++ Standard Libraries
#include <vector>

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "particles.hpp"

// Forward definitions
class ParticleGravity;

//--------------------------------------------------------------------------------------
//! \class DustParticles
//! \brief defines the class for dust particles that interact with the gas via drag
//!        force.

class DustParticles : public Particles {
friend class MeshBlock;
friend class ParticleGravity;

 public:
  // Class method
  static void FindDensityOnMesh(Mesh *pm, bool include_momentum);
  static void Initialize(Mesh *pm, ParameterInput *pin);
  static void SetOneParticleMass(Real new_mass);
  static bool GetBackReaction();
  static bool GetVariableTaus();
  static Real GetOneParticleMass();
  static Real GetStoppingTime();

  //!Constructor
  DustParticles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  ~DustParticles();

  // Accessors
  AthenaArray<Real> GetMassDensity() const;
  AthenaArray<Real> GetVelocityField() const;

  // Instance method
  Real NewBlockTimeStep();

 private:
  // Class variables
  static bool initialized;    //!> whether or not the class is initialized
  static bool backreaction;   //!> turn on/off back reaction
  static bool dragforce;      //!> turn on/off drag force
  static bool variable_taus;  //!> whether or not the stopping time is variable

  static int iwx, iwy, iwz;         // indices for working arrays
  static int idpx1, idpx2, idpx3;   // indices for momentum change
  static int itaus;                 //!> index for stopping time

  static Real mass;   //!> mass of each particle
  static Real taus0;  //!> constant/default stopping time (in code units)

  // Instance methods.
  void SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void UserSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void UserStoppingTime(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                     AthenaArray<Real>& meshdst);

  // Instance variables
  std::vector<Real> &wx, &wy, &wz;        // shorthand for working arrays
  AthenaArray<Real> dpx1, dpx2, dpx3;     // shorthand for momentum change
  std::vector<Real> &taus;                // shorthand for stopping time
  ParticleGravity *ppgrav;
};

//--------------------------------------------------------------------------------------
//! \fn bool DustParticles::GetBackReaction()
//! \brief returns if the back reaction of the drag is on or off.

inline bool DustParticles::GetBackReaction() {
  return backreaction;
}

//--------------------------------------------------------------------------------------
//! \fn bool DustParticles::GetVariableTaus()
//! \brief returns if the stopping time can be variable or not.

inline bool DustParticles::GetVariableTaus() {
  return variable_taus;
}

//--------------------------------------------------------------------------------------
//! \fn Real DustParticles::GetOneParticleMass()
//! \brief returns the mass of each particle.

inline Real DustParticles::GetOneParticleMass() {
  return mass;
}

//--------------------------------------------------------------------------------------
//! \fn Real DustParticles::GetStoppingTime()
//! \brief returns the stopping time of the drag.

inline Real DustParticles::GetStoppingTime() {
  return taus0;
}

//--------------------------------------------------------------------------------------
//! \fn AthenaArray<Real> DustParticles::GetMassDensity()
//! \brief returns the mass density of particles on the mesh.
//!
//! \note
//!  Precondition:
//!  The particle properties on mesh must be assigned using the class method
//!  DustParticles::FindDensityOnMesh().

inline AthenaArray<Real> DustParticles::GetMassDensity() const {
  return ppm->weight;
}

#endif  // PARTICLES_DUST_PARTICLES_HPP_
