#ifndef PARTICLES_PARTICLES_HPP_
#define PARTICLES_PARTICLES_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particles.hpp
//! \brief defines abstract base class Particles.
//======================================================================================

// C/C++ Standard Libraries
#include <string>
#include <vector>

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "particle_buffer.hpp"
#include "particle-mesh.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Forward definitions
class OutputType;
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
friend class DustParticles;

 public:
  // Class methods
  static void AMRCoarseToFine(MeshBlock* pmbc, MeshBlock* pmbf);
  static void AMRFineToCoarse(MeshBlock* pmbf, MeshBlock* pmbc);
  static void Initialize(Mesh *pm, ParameterInput *pin);
  static void PostInitialize(Mesh *pm, ParameterInput *pin);
  static void FindDensityOnMesh(Mesh *pm, bool include_momentum);
  static void FindHistoryOutput(Mesh *pm, Real data_sum[], int pos);
  static void GetHistoryOutputNames(std::string output_names[]);
  static const std::vector<std::string>& GetIntNames() { return ipname; }
  static const std::vector<std::string>& GetRealNames() { return rpname; }
  static int GetNInt() { return nint; }
  static int GetNReal() { return nreal; }
  static int GetTotalNumber(Mesh *pm);

  // Class constant
  static const int NHISTORY = 7;  //!> number of variables in history output

  // Constructor
  Particles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  virtual ~Particles();

  // Accessor
  virtual AthenaArray<Real> GetMassDensity() const;
  virtual AthenaArray<Real> GetVelocityField() const;
  const std::vector<int>* GetIntProps() const { return intprop; }
  const std::vector<Real>* GetRealProps() const { return rp; }
  Real GetMaximumWeight() const;
  int GetNPar() const { return npar; }

  // Instance methods
  void ClearBoundary();
  void ClearNeighbors();
  void Integrate(int stage, Real t, Real dt, Real gamma[]);
  void LinkNeighbors(MeshBlockTree &tree, int64_t nrbx1, int64_t nrbx2, int64_t nrbx3,
                     int root_level);
  void RemoveOneParticle(int k);
  void SendParticleMesh();
  void SendToNeighbors();
  void SetPositionIndices();
  void StartReceiving();
  bool ReceiveFromNeighbors();
  bool ReceiveParticleMesh(Real t, Real dt);
  virtual Real NewBlockTimeStep();

  std::size_t GetSizeInBytes();
  void UnpackParticlesForRestart(char *mbdata, std::size_t &os);
  void PackParticlesForRestart(char *&pdata);

 protected:
  // Class methods
  static int AddIntProperty(const std::string& name);
  static int AddRealProperty(const std::string& name);
  static int AddAuxProperty();
  static int AddWorkingArray();

  // Instance methods
  virtual void SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc);

  // Class variables
  static bool initialized;  //!> whether or not the class is initialized
  static int nint;          //!> numbers of integer particle properties
  static int nreal;         //!> numbers of real particle properties
  static int naux;          //!> number of auxiliary particle properties
  static int nwork;         //!> number of working arrays for particles

  static int ipid;                 //!> index for the particle ID
  static int ixp, iyp, izp;        // indices for the position components
  static int ivpx, ivpy, ivpz;     // indices for the velocity components

  static int ixi1, ixi2, ixi3;     // indices for position indices

  static int imom1, imom2, imom3;  // indices for momentum components on mesh

  static Real cfl_par;  //!> CFL number for particles

  static ParameterInput *pinput;

  // Instance methods
  void Resize(int new_npar);  //!> Change number of particles

  // Instance variables
  int npar;     //!> number of particles

                                     // Data attached to the particles:
  std::vector<int> *intprop;         //!>   integer properties
  std::vector<Real> *rp, *rp1, *drp; //!>   real properties
  std::vector<Real> *aux;            //!>   auxiliary properties (communicated when
                                     //!>     particles moving to another meshblock)
  std::vector<Real> *work;           //!>   working arrays (not communicated)

  ParticleMesh *ppm;  //!> ptr to particle-mesh

                                          // Shorthands:
  std::vector<int> &pid;                  //!>   particle ID
  std::vector<Real> &xp, &yp, &zp;        //   position
  std::vector<Real> &vpx, &vpy, &vpz;     //   velocity
  std::vector<Real> &dxp, &dyp, &dzp;     //   rate of position change
  std::vector<Real> &dvpx, &dvpy, &dvpz;  //   rate of velocity change
  std::vector<Real> &xi1, &xi2, &xi3;     //   position indices in local meshblock

  MeshBlock* pmy_block;  //!> MeshBlock pointer
  Mesh* pmy_mesh;        //!> Mesh pointer

 private:
  // Class method
  static void ProcessNewParticles(Mesh *pmesh);

  // Instance methods
  virtual void UserSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) {}
  virtual void ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc) {}
  virtual void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                             AthenaArray<Real>& meshdst) {}

  int CountNewParticles() const;
  void ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3);
  void FlushReceiveBuffer(ParticleBuffer& recv);
  void GetPositionIndices(
      int ibegin, int iend,
      const std::vector<Real>& xp,
      const std::vector<Real>& yp,
      const std::vector<Real>& zp,
      std::vector<Real>& xi1, std::vector<Real>& xi2, std::vector<Real>& xi3);
  void RealPropCopy(std::vector<Real> *rp1, const std::vector<Real> *rp2);
  void RealPropSwap(std::vector<Real> *rp1, std::vector<Real> *rp2);
  void SetNewParticleID(int id);
  void WeightedAverage(std::vector<Real> *rp_out,
      const std::vector<Real> *rp_in1, const Real weights[]);
  struct Neighbor* FindTargetNeighbor(
      int ox1, int ox2, int ox3, int xi1, int xi2, int xi3);

  // Class variable
  static std::vector<std::string> ipname;  //!> names of integer properties
  static std::vector<std::string> rpname;  //!> names of real properties
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
//! \fn AthenaArray<Real> Particles::GetMassDensity()
//! \brief returns the mass density of particles on the mesh.

inline AthenaArray<Real> Particles::GetMassDensity() const {
  RegionSize& block_size(pmy_block->block_size);
  return AthenaArray<Real>(block_size.nx3, block_size.nx2, block_size.nx1);
}

//--------------------------------------------------------------------------------------
//! \fn AthenaArray<Real> Particles::GetVelocityField()
//! \brief returns the velocity field of particles on the mesh.

inline AthenaArray<Real> Particles::GetVelocityField() const {
  RegionSize& block_size(pmy_block->block_size);
  return AthenaArray<Real>(3, block_size.nx3, block_size.nx2, block_size.nx1);
}

//--------------------------------------------------------------------------------------
//! \fn Real Particles::GetMaximumWeight()
//! \brief returns the maximum weight on the mesh.

inline Real Particles::GetMaximumWeight() const {
  return ppm->FindMaximumWeight();
}

#endif  // PARTICLES_PARTICLES_HPP_
