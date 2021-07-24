#ifndef PARTICLES_PARTICLES_HPP_
#define PARTICLES_PARTICLES_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particles.hpp
//! \brief defines classes for particle dynamics.
//======================================================================================

// C/C++ Standard Libraries
#include <string>
#include <vector>

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../outputs/outputs.hpp"
#include "../parameter_input.hpp"
#include "particle_buffer.hpp"
#include "particle-mesh.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Forward definitions
class ParticleGravity;

//--------------------------------------------------------------------------------------
//! \struct Neighbor
//  \brief defines a structure for links to neighbors

struct Neighbor {
  NeighborBlock *pnb;
  MeshBlock *pmb;
  Neighbor *next, *prev;

  Neighbor() : pnb(NULL), pmb(NULL), next(NULL), prev(NULL) {}
};

//--------------------------------------------------------------------------------------
//! \class Particles
//! \brief defines the base class for all implementations of particles.

class Particles {
friend class MeshBlock;  // Make writing initial conditions possible.
friend class OutputType;
friend class ParticleGravity;
friend class ParticleMesh;

 public:
  // Class methods
  static void AMRCoarseToFine(MeshBlock* pmbc, MeshBlock* pmbf);
  static void AMRFineToCoarse(MeshBlock* pmbf, MeshBlock* pmbc);
  static void Initialize(Mesh *pm, ParameterInput *pin);
  static void PostInitialize(Mesh *pm, ParameterInput *pin);
  static void FindDensityOnMesh(Mesh *pm, bool include_momentum);
  static void FindHistoryOutput(Mesh *pm, Real data_sum[], int pos);
  static void FormattedTableOutput(Mesh *pm, OutputParameters op);
  static void GetHistoryOutputNames(std::string output_names[]);
  static int GetTotalNumber(Mesh *pm);

  // Class constant
  static const int NHISTORY = 7;  //!> number of variables in history output

  // Constructor
  Particles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  virtual ~Particles();

  // Accessor
  Real GetMaximumWeight() const;

  // Instance methods
  void ClearBoundary();
  void ClearNeighbors();
  void Integrate(int step);
  void LinkNeighbors(MeshBlockTree &tree, int64_t nrbx1, int64_t nrbx2, int64_t nrbx3,
                     int root_level);
  void RemoveOneParticle(int k);
  void SendParticleMesh();
  void SendToNeighbors();
  void SetPositionIndices();
  void StartReceiving();
  bool ReceiveFromNeighbors();
  bool ReceiveParticleMesh(int step);
  Real NewBlockTimeStep();

  std::size_t GetSizeInBytes();
  void UnpackParticlesForRestart(char *mbdata, std::size_t &os);
  void PackParticlesForRestart(char *&pdata);

 protected:
  // Class methods
  static int AddIntProperty();
  static int AddRealProperty();
  static int AddAuxProperty();
  static int AddWorkingArray();

  // Class variables
  static bool initialized;  //!> whether or not the class is initialized
  static int nint;          //!> numbers of integer particle properties
  static int nreal;         //!> numbers of real particle properties
  static int naux;          //!> number of auxiliary particle properties
  static int nwork;         //!> number of working arrays for particles

  static int ipid;                 //!> index for the particle ID
  static int ixp, iyp, izp;        // indices for the position components
  static int ivpx, ivpy, ivpz;     // indices for the velocity components

  static int ixp0, iyp0, izp0;     // indices for beginning position components
  static int ivpx0, ivpy0, ivpz0;  // indices for beginning velocity components

  static int ixi1, ixi2, ixi3;     // indices for position indices

  static int imom1, imom2, imom3;  // indices for momentum components on mesh

  static Real cfl_par;  //!> CFL number for particles

  static ParameterInput *pinput;

  // Instance methods
  virtual void AssignShorthands();  //!> Needs to be called everytime
                                    //!> realprop, & auxprop are resized
                                    //!> Be sure to call back when derived.

  void UpdateCapacity(int new_nparmax);  //!> Change the capacity of particle arrays
  void Resize(int new_npar);  //!> Change number of particles

  // Instance variables
  int npar;     //!> number of particles
  int nparmax;  //!> maximum number of particles per meshblock

                               // Data attached to the particles:
  std::vector<int> *intprop;   //!>   integer properties
  AthenaArray<Real> realprop;  //!>   real properties
  AthenaArray<Real> auxprop;   //!>   auxiliary properties (communicated when
                               //!>     particles moving to another meshblock)
  AthenaArray<Real> work;      //!>   working arrays (not communicated)

  ParticleMesh *ppm;  //!> ptr to particle-mesh

                                       // Shorthands:
  std::vector<int>& pid;               //!>   particle ID
  AthenaArray<Real> xp, yp, zp;        //   position
  AthenaArray<Real> vpx, vpy, vpz;     //   velocity
  AthenaArray<Real> xi1, xi2, xi3;     //   position indices in local meshblock
  AthenaArray<Real> xp0, yp0, zp0;     //   beginning position
  AthenaArray<Real> vpx0, vpy0, vpz0;  //   beginning velocity

  MeshBlock* pmy_block;  //!> MeshBlock pointer
  Mesh* pmy_mesh;        //!> Mesh pointer

 private:
  // Class method
  static void ProcessNewParticles(Mesh *pmesh);

  // Instance methods
  virtual void SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) {}
  virtual void UserSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) {}
  virtual void ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc) {}
  virtual void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                             AthenaArray<Real>& meshdst) {}

  int CountNewParticles() const;
  void ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3);
  void EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void FlushReceiveBuffer(ParticleBuffer& recv);
  void GetPositionIndices(int npar,
                          const AthenaArray<Real>& xp,
                          const AthenaArray<Real>& yp,
                          const AthenaArray<Real>& zp,
                          AthenaArray<Real>& xi1,
                          AthenaArray<Real>& xi2,
                          AthenaArray<Real>& xi3);
  void SetNewParticleID(int id);
  void SaveStatus();
  struct Neighbor* FindTargetNeighbor(
      int ox1, int ox2, int ox3, int xi1, int xi2, int xi3);

  // Class variable
  static int idmax;

  // Instance variables
  bool active1_, active2_, active3_;  // active dimensions

  // MeshBlock-to-MeshBlock communication:
  BoundaryValues *pbval_;            //!> ptr to my BoundaryValues
  Neighbor neighbor_[3][3][3];       //!> links to neighbors
  ParticleBuffer recv_[56];          //!> particle receive buffers
  enum BoundaryStatus bstatus_[56];  //!> boundary status
#ifdef MPI_PARALLEL
  static MPI_Comm my_comm;   //!> my MPI communicator
  ParticleBuffer send_[56];  //!> particle send buffers
#endif
};

//--------------------------------------------------------------------------------------
//! \fn Real Particles::GetMaximumWeight()
//! \brief returns the maximum weight on the mesh.

inline Real Particles::GetMaximumWeight() const {
  return ppm->FindMaximumWeight();
}

//--------------------------------------------------------------------------------------
//! \class DustParticles
//! \brief defines the class for dust particles that interact with the gas via drag
//!        force.

class DustParticles : public Particles {
friend class MeshBlock;

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
  void AssignShorthands();
  void SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void UserSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void UserStoppingTime(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                     AthenaArray<Real>& meshdst);

  // Instance variables
  AthenaArray<Real> wx, wy, wz;        // shorthand for working arrays
  AthenaArray<Real> dpx1, dpx2, dpx3;  // shorthand for momentum change
  AthenaArray<Real> taus;              // shorthand for stopping time
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

#endif  // PARTICLES_PARTICLES_HPP_
