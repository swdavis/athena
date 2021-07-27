#ifndef PARTICLES_PARTICLE_MESH_HPP_
#define PARTICLES_PARTICLE_MESH_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle-mesh.hpp
//! \brief defines ParticleMesh class used for communication between meshblocks needed
//!        by particle-mesh methods.

// C++ standard library
#include <cmath>
#include <vector>

// Athena++ classes headers
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../mesh/mesh.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Particle-mesh constants.
const Real RINF = 1;  //!> radius of influence

//! Define the size of a particle cloud = 2 * RINF + 1
#define NPC 3

// Forward declaration
class Particles;
class ParameterInput;

//--------------------------------------------------------------------------------------
//! \class ParticleMesh
//! \brief defines the class for particle-mesh methods

class ParticleMesh {
friend class Particles;
friend class DustParticles;
friend class ParticleGravity;
friend class OutputType;

 public:
  // Class methods
  static void Initialize(ParameterInput *pin);
  static int AddMeshAux();

  // Constructor and destructor
  explicit ParticleMesh(Particles *ppar);
  ~ParticleMesh();

  // Accessor
  Real FindMaximumWeight() const;

 protected:
  // Class variables
  static int nmeshaux;  //!> number of auxiliaries to the meshblock
  static int iweight;   //!> index to weight in meshaux

  // Instance variables
  AthenaArray<Real> meshaux;   //!> auxiliaries to the meshblock
  int is, ie, js, je, ks, ke;  // beginning and ending indices
  AthenaArray<Real> weight;    //!> shorthand to weight in meshaux

  // Instance methods
  void InterpolateMeshToParticles(
           const AthenaArray<Real>& meshsrc, int ms1,
           std::vector<Real> par[], int p1, int nprop);
  void AssignParticlesToMeshAux(
           const std::vector<Real> par[], int p1, int ma1, int nprop);
  void InterpolateMeshAndAssignParticles(
           const AthenaArray<Real>& meshsrc, int ms1,
           std::vector<Real> pardst[], int pd1, int ni,
           const std::vector<Real> parsrc[], int ps1, int ma1, int na);
  void DepositMeshAux(AthenaArray<Real>& u, int ma1, int mb1, int nprop);

  void ClearBoundary();
  void SendBoundary();
  void StartReceiving();
  bool ReceiveBoundary();

 private:
  struct BoundaryAttributes {
    Real xi1min, xi1max, xi2min, xi2max, xi3min, xi3max;
                               // domain that influences the ghost block
    Real xi1_0, xi2_0, xi3_0;  // origin of the ghost block wrt to the local meshblock
    int ngx1, ngx2, ngx3;      // dimensions of the ghost block
    int ngtot;                 //!> total number of cells in the ghost block
    int irs, ire, jrs, jre, krs, kre;  // beginning/ending indices in meshaux to receive
    int iss, ise, jss, jse, kss, kse;  // beginning/ending indices in meshaux to send
  };

  // Class variables
  static bool initialized_;

  // Instance Variables
  bool active1_, active2_, active3_;  // active dimensions
  Real dxi1_, dxi2_, dxi3_;           // range of influence from a particle cloud
  int nx1_, nx2_, nx3_;               // number of cells in meshaux in each dimension
  int ncells_;                        // total number of cells in meshaux
  int npc1_, npc2_, npc3_;            // size of a particle cloud

  Particles *ppar_;            //!> ptr to my Particles instance
  MeshBlock *pmb_;             //!> ptr to my MeshBlock
  Mesh *pmesh_;                //!> ptr to my Mesh
  BoundaryValues *pbval_;      //!> ptr to my BoundaryValues
  BoundaryData<> bd_;          //!> boundary data
  BoundaryAttributes ba_[56];  //!> ghost block attributes

  // Instance methods
  void InitiateBoundaryData();
  void SetBoundaryAttributes();
  void AssignParticlesToDifferentLevels(
           const std::vector<Real> par[], int p1, int ma1, int nprop);
  int LoadBoundaryBufferSameLevel(Real *buf, const BoundaryAttributes& ba);
  void AddBoundaryBuffer(Real *buf, const BoundaryAttributes& ba);

#ifdef MPI_PARALLEL
  // Variables for MPI
  static MPI_Comm my_comm;  //!> my MPI communicator
#endif
};
#endif  // PARTICLES_PARTICLE_MESH_HPP_
