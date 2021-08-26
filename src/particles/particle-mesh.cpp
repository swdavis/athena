//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle-mesh.cpp
//! \brief implements ParticleMesh class used for operations involved in particle-mesh
//!        methods.

// Standard library
#include <algorithm>
#include <cstring>
#include <sstream>

// Athena++ classes headers
#include "../athena.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../utils/buffer_utils.hpp"
#include "particle-mesh.hpp"
#include "particles.hpp"

// Class variable initialization
bool ParticleMesh::initialized_ = false;
int ParticleMesh::nmeshaux = 0;
int ParticleMesh::iweight = -1;
#ifdef MPI_PARALLEL
MPI_Comm ParticleMesh::my_comm = MPI_COMM_NULL;
#endif

// Local function prototypes.
static Real _WeightFunction(Real dxi);

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::Initialize(ParameterInput *pin)
//! \brief initiates the ParticleMesh class.

void ParticleMesh::Initialize(ParameterInput *pin) {
  if (initialized_) return;

  // Add weight in meshaux.
  iweight = AddMeshAux();

#ifdef MPI_PARALLEL
  // Get my MPI communicator.
  MPI_Comm_dup(MPI_COMM_WORLD, &my_comm);
#endif

  initialized_ = true;
}

//--------------------------------------------------------------------------------------
//! \fn int ParticleMesh::AddMeshAux()
//! \brief adds one auxiliary to the mesh and returns the index.

int ParticleMesh::AddMeshAux() {
  return nmeshaux++;
}

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::ParticleMesh(Particles *ppar, int nmeshaux)
//! \brief constructs a new ParticleMesh instance.

ParticleMesh::ParticleMesh(Particles *ppar) {
  // Save some inputs.
  ppar_ = ppar;
  pmb_ = ppar->pmy_block;
  pmesh_ = pmb_->pmy_mesh;
  pbval_ = pmb_->pbval;

  // Determine active dimensions.
  RegionSize& block_size = pmb_->block_size;
  active1_ = block_size.nx1 > 1;
  active2_ = block_size.nx2 > 1;
  active3_ = block_size.nx3 > 1;

  // Determine the range of a particle cloud.
  dxi1_ = active1_ ? RINF : 0;
  dxi2_ = active2_ ? RINF : 0;
  dxi3_ = active3_ ? RINF : 0;

  // Establish the particle-mesh blocks.
  nx1_ = active1_ ? block_size.nx1 + 2 * NGHOST : 1;
  nx2_ = active2_ ? block_size.nx2 + 2 * NGHOST : 1;
  nx3_ = active3_ ? block_size.nx3 + 2 * NGHOST : 1;
  meshaux.NewAthenaArray(nmeshaux, nx3_, nx2_, nx1_);
  ncells_ = nx1_ * nx2_ * nx3_;

  is = pmb_->is;
  ie = pmb_->ie;
  js = pmb_->js;
  je = pmb_->je;
  ks = pmb_->ks;
  ke = pmb_->ke;

  // Get a shorthand to weights.
  weight.InitWithShallowSlice(meshaux, 4, iweight, 1);

  // Determine the dimensions of each particle cloud.
  npc1_ = active1_ ? NPC : 1;
  npc2_ = active2_ ? NPC : 1;
  npc3_ = active3_ ? NPC : 1;

  // Initialize boundary data.
  bd_.nbmax = 56;
  for (int n = 0; n < bd_.nbmax; n++) {
    bd_.flag[n] = BoundaryStatus::waiting;
    bd_.send[n] = NULL;
    bd_.recv[n] = NULL;
#ifdef MPI_PARALLEL
    bd_.req_recv[n] = MPI_REQUEST_NULL;
    bd_.req_send[n] = MPI_REQUEST_NULL;
#endif
  }
}

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::~ParticleMesh()
//! \brief destructs a ParticleMesh instance.

ParticleMesh::~ParticleMesh() {
  // Destroy the particle meshblock.
  weight.DeleteAthenaArray();
  meshaux.DeleteAthenaArray();

  // Destroy boundary data.
  for (int n = 0; n < bd_.nbmax; n++) {
    if (bd_.recv[n] != NULL) delete [] bd_.recv[n];
#ifdef MPI_PARALLEL
    if (bd_.send[n] != NULL) delete [] bd_.send[n];
    if (bd_.req_recv[n] != MPI_REQUEST_NULL) MPI_Request_free(&bd_.req_recv[n]);
    if (bd_.req_send[n] != MPI_REQUEST_NULL) MPI_Request_free(&bd_.req_send[n]);
#endif
  }
}

//--------------------------------------------------------------------------------------
//! \fn Real ParticleMesh::FindMaximumWeight()
//! \brief returns the maximum weight in the meshblock.

Real ParticleMesh::FindMaximumWeight() const {
  Real wmax = 0.0;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        wmax = std::max(wmax, weight(k,j,i));
  return wmax;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InterpolateMeshToParticles(
//!              const AthenaArray<Real>& meshsrc, int ms1,
//!              std::vector<Real> par[], int p1, int nprop)
//! \brief interpolates meshsrc from property index ms1 to ms1+nprop-1 onto particle
//!     array par (realprop, auxprop, or work in Particles class) from property index p1
//!     to p1+nprop-1.

void ParticleMesh::InterpolateMeshToParticles(
         const AthenaArray<Real>& meshsrc, int ms1,
         std::vector<Real> par[], int p1, int nprop) {
  // Transpose meshsrc.
  int nx1 = meshsrc.GetDim1(), nx2 = meshsrc.GetDim2(), nx3 = meshsrc.GetDim3();
  AthenaArray<Real> u;
  u.NewAthenaArray(nx3,nx2,nx1,nprop);
  for (int n = 0; n < nprop; ++n)
    for (int k = 0; k < nx3; ++k)
      for (int j = 0; j < nx2; ++j)
        for (int i = 0; i < nx1; ++i)
          u(k,j,i,n) = meshsrc(ms1+n,k,j,i);

  // Allocate space for SIMD.
  Real **w1 __attribute__((aligned(64))) = new Real*[npc1_];
  Real **w2 __attribute__((aligned(64))) = new Real*[npc2_];
  Real **w3 __attribute__((aligned(64))) = new Real*[npc3_];
  for (int i = 0; i < npc1_; ++i)
    w1[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc2_; ++i)
    w2[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc3_; ++i)
    w3[i] = new Real[SIMD_WIDTH];
  Real imb1v[SIMD_WIDTH] __attribute__((aligned(64)));
  Real imb2v[SIMD_WIDTH] __attribute__((aligned(64)));
  Real imb3v[SIMD_WIDTH] __attribute__((aligned(64)));

  // Loop over each particle.
  int npar = ppar_->npar;
  for (int k = 0; k < npar; k += SIMD_WIDTH) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Find the domain the particle influences.
      Real xi1(ppar_->xi1[kkk]), xi2(ppar_->xi2[kkk]), xi3(ppar_->xi3[kkk]);
      int imb1 = static_cast<int>(xi1 - dxi1_),
          imb2 = static_cast<int>(xi2 - dxi2_),
          imb3 = static_cast<int>(xi3 - dxi3_);
      xi1 = imb1 + 0.5 - xi1;
      xi2 = imb2 + 0.5 - xi2;
      xi3 = imb3 + 0.5 - xi3;

      imb1v[kk] = imb1;
      imb2v[kk] = imb2;
      imb3v[kk] = imb3;

      // Weigh each cell.
      for (int i = 0; i < npc1_; ++i)
        w1[i][kk] = active1_ ? _WeightFunction(xi1 + i) : 1.0;
      for (int i = 0; i < npc2_; ++i)
        w2[i][kk] = active2_ ? _WeightFunction(xi2 + i) : 1.0;
      for (int i = 0; i < npc3_; ++i)
        w3[i][kk] = active3_ ? _WeightFunction(xi3 + i) : 1.0;
    }

    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Initiate interpolation.
      Real *pd = new Real[nprop];
      for (int i = 0; i < nprop; ++i)
        pd[i] = 0.0;

      int imb1 = imb1v[kk], imb2 = imb2v[kk], imb3 = imb3v[kk];

      for (int ipc3 = 0; ipc3 < npc3_; ++ipc3) {
        for (int ipc2 = 0; ipc2 < npc2_; ++ipc2) {
          for (int ipc1 = 0; ipc1 < npc1_; ++ipc1) {
            Real w = w1[ipc1][kk] * w2[ipc2][kk] * w3[ipc3][kk];

            // Interpolate meshsrc to particles.
            for (int n = 0; n < nprop; ++n)
              pd[n] += w * u(imb3+ipc3,imb2+ipc2,imb1+ipc1,n);
          }
        }
      }

      // Record the final interpolated properties.
      for (int n = 0; n < nprop; ++n)
        par[p1+n][kkk] = pd[n];

      delete [] pd;
    }
  }

  // Release working arrays.
  u.DeleteAthenaArray();
  for (int i = 0; i < npc1_; ++i)
    delete [] w1[i];
  for (int i = 0; i < npc2_; ++i)
    delete [] w2[i];
  for (int i = 0; i < npc3_; ++i)
    delete [] w3[i];
  delete [] w1;
  delete [] w2;
  delete [] w3;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AssignParticlesToMeshAux(
//       const std::vector<Real> par[], int p1, int ma1, int nprop)
//! \brief assigns par (realprop, auxprop, or work in Particles class) from property
//!        index p1 to p1+nprop-1 onto meshaux from property index ma1 and up.

void ParticleMesh::AssignParticlesToMeshAux(
         const std::vector<Real> par[], int p1, int ma1, int nprop) {
  // Zero out meshaux.
  std::fill(&weight(0,0,0), &weight(0,0,0) + ncells_, 0.0);
  std::fill(&meshaux(ma1,0,0,0), &meshaux(ma1+nprop,0,0,0), 0.0);

  // Allocate space for SIMD.
  Real **w1 __attribute__((aligned(64))) = new Real*[npc1_];
  Real **w2 __attribute__((aligned(64))) = new Real*[npc2_];
  Real **w3 __attribute__((aligned(64))) = new Real*[npc3_];
  for (int i = 0; i < npc1_; ++i)
    w1[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc2_; ++i)
    w2[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc3_; ++i)
    w3[i] = new Real[SIMD_WIDTH];
  Real imb1v[SIMD_WIDTH] __attribute__((aligned(64)));
  Real imb2v[SIMD_WIDTH] __attribute__((aligned(64)));
  Real imb3v[SIMD_WIDTH] __attribute__((aligned(64)));

  // Loop over each particle.
  int npar = ppar_->npar;
  for (int k = 0; k < npar; k += SIMD_WIDTH) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Find the domain the particle influences.
      Real xi1(ppar_->xi1[kkk]), xi2(ppar_->xi2[kkk]), xi3(ppar_->xi3[kkk]);
      int imb1 = static_cast<int>(xi1 - dxi1_),
          imb2 = static_cast<int>(xi2 - dxi2_),
          imb3 = static_cast<int>(xi3 - dxi3_);
      xi1 = imb1 + 0.5 - xi1;
      xi2 = imb2 + 0.5 - xi2;
      xi3 = imb3 + 0.5 - xi3;

      imb1v[kk] = imb1;
      imb2v[kk] = imb2;
      imb3v[kk] = imb3;

      // Weigh each cell.
      for (int i = 0; i < npc1_; ++i)
        w1[i][kk] = active1_ ? _WeightFunction(xi1 + i) : 1.0;
      for (int i = 0; i < npc2_; ++i)
        w2[i][kk] = active2_ ? _WeightFunction(xi2 + i) : 1.0;
      for (int i = 0; i < npc3_; ++i)
        w3[i][kk] = active3_ ? _WeightFunction(xi3 + i) : 1.0;
    }

    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Fetch particle properties.
      Real *ps = new Real[nprop];
      for (int i = 0; i < nprop; ++i)
        ps[i] = par[p1+i][kkk];

      int imb1 = imb1v[kk], imb2 = imb2v[kk], imb3 = imb3v[kk];

      for (int ipc3 = 0; ipc3 < npc3_; ++ipc3) {
        for (int ipc2 = 0; ipc2 < npc2_; ++ipc2) {
          for (int ipc1 = 0; ipc1 < npc1_; ++ipc1) {
            Real w = w1[ipc1][kk] * w2[ipc2][kk] * w3[ipc3][kk];

            // Record the weights.
            weight(imb3+ipc3,imb2+ipc2,imb1+ipc1) += w;

            // Assign particles to meshaux.
            for (int n = 0; n < nprop; ++n)
              meshaux(ma1+n,imb3+ipc3,imb2+ipc2,imb1+ipc1) += w * ps[n];
          }
        }
      }

      delete [] ps;
    }
  }

  // Release working array.
  for (int i = 0; i < npc1_; ++i)
    delete [] w1[i];
  for (int i = 0; i < npc2_; ++i)
    delete [] w2[i];
  for (int i = 0; i < npc3_; ++i)
    delete [] w3[i];
  delete [] w1;
  delete [] w2;
  delete [] w3;

  // Treat neighbors of different levels.
  if (pmesh_->multilevel)
    AssignParticlesToDifferentLevels(par, p1, ma1, nprop);
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InterpolateMeshAndAssignParticles(
//!              const AthenaArray<Real>& meshsrc, int ms1,
//!              std::vector<Real> pardst[], int pd1, int ni,
//!              const std::vector<Real> parsrc[], int ps1, int ma1, int na)
//! \brief interpolates meshsrc from property index ms1 to ms1 + ni - 1 onto particle
//!     array pardst from index pd1 to pd1 + ni - 1, and assigns parsrc from property
//!     index ps1 to ps1 + na - 1 onto meshaux from ma1 to ma1 + na - 1.  The arrays
//!     parsrc and pardst can be realprop, auxprop, or work in Particles class.

void ParticleMesh::InterpolateMeshAndAssignParticles(
         const AthenaArray<Real>& meshsrc, int ms1,
         std::vector<Real> pardst[], int pd1, int ni,
         const std::vector<Real> parsrc[], int ps1, int ma1, int na) {
  // Zero out meshaux.
  std::fill(&weight(0,0,0), &weight(0,0,0) + ncells_, 0.0);
  std::fill(&meshaux(ma1,0,0,0), &meshaux(ma1+na,0,0,0), 0.0);

  // Transpose meshsrc.
  int nx1 = meshsrc.GetDim1(), nx2 = meshsrc.GetDim2(), nx3 = meshsrc.GetDim3();
  AthenaArray<Real> u;
  u.NewAthenaArray(nx3,nx2,nx1,ni);
  for (int n = 0; n < ni; ++n)
    for (int k = 0; k < nx3; ++k)
      for (int j = 0; j < nx2; ++j)
        for (int i = 0; i < nx1; ++i)
          u(k,j,i,n) = meshsrc(ms1+n,k,j,i);

  // Allocate space for SIMD.
  Real **w1 __attribute__((aligned(64))) = new Real*[npc1_];
  Real **w2 __attribute__((aligned(64))) = new Real*[npc2_];
  Real **w3 __attribute__((aligned(64))) = new Real*[npc3_];
  for (int i = 0; i < npc1_; ++i)
    w1[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc2_; ++i)
    w2[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc3_; ++i)
    w3[i] = new Real[SIMD_WIDTH];
  Real imb1v[SIMD_WIDTH] __attribute__((aligned(64)));
  Real imb2v[SIMD_WIDTH] __attribute__((aligned(64)));
  Real imb3v[SIMD_WIDTH] __attribute__((aligned(64)));

  // Loop over each particle.
  int npar = ppar_->npar;
  for (int k = 0; k < npar; k += SIMD_WIDTH) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Find the domain the particle influences.
      Real xi1(ppar_->xi1[kkk]), xi2(ppar_->xi2[kkk]), xi3(ppar_->xi3[kkk]);
      int imb1 = static_cast<int>(xi1 - dxi1_),
          imb2 = static_cast<int>(xi2 - dxi2_),
          imb3 = static_cast<int>(xi3 - dxi3_);
      xi1 = imb1 + 0.5 - xi1;
      xi2 = imb2 + 0.5 - xi2;
      xi3 = imb3 + 0.5 - xi3;

      imb1v[kk] = imb1;
      imb2v[kk] = imb2;
      imb3v[kk] = imb3;

      // Weigh each cell.
      for (int i = 0; i < npc1_; ++i)
        w1[i][kk] = active1_ ? _WeightFunction(xi1 + i) : 1.0;
      for (int i = 0; i < npc2_; ++i)
        w2[i][kk] = active2_ ? _WeightFunction(xi2 + i) : 1.0;
      for (int i = 0; i < npc3_; ++i)
        w3[i][kk] = active3_ ? _WeightFunction(xi3 + i) : 1.0;
    }

    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Initiate interpolation and fetch particle properties.
      Real *pd = new Real[ni];
      Real *ps = new Real[na];
      for (int i = 0; i < ni; ++i)
        pd[i] = 0.0;
      for (int i = 0; i < na; ++i)
        ps[i] = parsrc[ps1+i][kkk];

      int imb1 = imb1v[kk], imb2 = imb2v[kk], imb3 = imb3v[kk];

      for (int ipc3 = 0; ipc3 < npc3_; ++ipc3) {
        for (int ipc2 = 0; ipc2 < npc2_; ++ipc2) {
          for (int ipc1 = 0; ipc1 < npc1_; ++ipc1) {
            Real w = w1[ipc1][kk] * w2[ipc2][kk] * w3[ipc3][kk];

            // Record the weights.
            weight(imb3+ipc3,imb2+ipc2,imb1+ipc1) += w;

            // Interpolate meshsrc to particles.
            for (int n = 0; n < ni; ++n)
              pd[n] += w * u(imb3+ipc3,imb2+ipc2,imb1+ipc1,n);

            // Assign particles to meshaux.
            for (int n = 0; n < na; ++n)
              meshaux(ma1+n,imb3+ipc3,imb2+ipc2,imb1+ipc1) += w * ps[n];
          }
        }
      }

      // Record the final interpolated properties.
      for (int n = 0; n < ni; ++n)
        pardst[pd1+n][kkk] = pd[n];

      delete [] pd;
      delete [] ps;
    }
  }

  // Release working array.
  u.DeleteAthenaArray();
  for (int i = 0; i < npc1_; ++i)
    delete [] w1[i];
  for (int i = 0; i < npc2_; ++i)
    delete [] w2[i];
  for (int i = 0; i < npc3_; ++i)
    delete [] w3[i];
  delete [] w1;
  delete [] w2;
  delete [] w3;

  // Treat neighbors of different levels.
  if (pmesh_->multilevel)
    AssignParticlesToDifferentLevels(parsrc, ps1, ma1, na);
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::DepositMeshAux(AthenaArray<Real>& u,
//!                                       int ma1, int mb1, int nprop)
//! \brief deposits data in meshaux from property index ma1 to ma1+nprop-1 to meshblock
//!        data u from property index mb1 and mb1+nprop-1, divided by cell volume.

void ParticleMesh::DepositMeshAux(AthenaArray<Real>& u, int ma1, int mb1, int nprop) {
  Coordinates *pc = pmb_->pcoord;

  for (int n = 0; n < nprop; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          u(mb1+n,k,j,i) += meshaux(ma1+n,k,j,i) / pc->GetCellVolume(k,j,i);
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InitiateBoundaryData()
//! \brief allocates space for boundary data.

void ParticleMesh::InitiateBoundaryData() {
  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    BoundaryAttributes& ba = ba_[n];

    if (bd_.recv[nb.bufid] != NULL) delete [] bd_.recv[nb.bufid];
    int nrecv = (ba.ire - ba.irs + 1) *
                (ba.jre - ba.jrs + 1) *
                (ba.kre - ba.krs + 1) * nmeshaux;
    bd_.recv[nb.bufid] = new Real[nrecv];

#ifdef MPI_PARALLEL
    int nb_rank = nb.snb.rank;
    if (nb_rank != Globals::my_rank) {
      int nsend = ba.ngtot * nmeshaux;
      bd_.send[nb.bufid] = new Real[nsend];
      MPI_Recv_init(bd_.recv[nb.bufid], nrecv, MPI_ATHENA_REAL, nb_rank,
                    (pmb_->lid<<6) | nb.bufid, my_comm, &bd_.req_recv[nb.bufid]);
      MPI_Send_init(bd_.send[nb.bufid], nsend, MPI_ATHENA_REAL, nb_rank,
                    (nb.snb.lid<<6) | nb.targetid, my_comm, &bd_.req_send[nb.bufid]);
    }
#endif
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::SetBoundaryAttributes()
//! \brief initializes or reinitializes attributes for each boundary.

void ParticleMesh::SetBoundaryAttributes() {
  const RegionSize& block_size = pmb_->block_size;
  const Real xi1mid = (is + ie + 1) / 2,
             xi2mid = (js + je + 1) / 2,
             xi3mid = (ks + ke + 1) / 2;
  const int mylevel = pmb_->loc.level;
  const int myfx1 = static_cast<int>(pbval_->loc.lx1 & 1L),
            myfx2 = static_cast<int>(pbval_->loc.lx2 & 1L),
            myfx3 = static_cast<int>(pbval_->loc.lx3 & 1L);
  const int NGH = (NGHOST + 1) / 2, NG2 = 2 * NGHOST;
  const int nx1h = active1_ ? block_size.nx1 / 2 + NGH : 1,
            nx2h = active2_ ? block_size.nx2 / 2 + NGH : 1,
            nx3h = active3_ ? block_size.nx3 / 2 + NGH : 1;

  // Loop over each neighbor block.
  for (int n = 0; n < pbval_->nneighbor; ++n) {
    NeighborBlock& nb = pbval_->neighbor[n];

    // Find the index domain.
    Real xi1min = is, xi1max = ie + 1,
         xi2min = js, xi2max = je + 1,
         xi3min = ks, xi3max = ke + 1;
    Real xi1_0 = xi1min, xi2_0 = xi2min, xi3_0 = xi3min;
    int irs = is, ire = ie, jrs = js, jre = je, krs = ks, kre = ke;
    int iss = is, ise = ie, jss = js, jse = je, kss = ks, kse = ke;

    // Find several depths from the neighbor block.
    Real dxip, dxig;
    int dxir;
    int nblevel = nb.snb.level;
    if (nblevel > mylevel) {
      dxip = 0.5 * RINF;
      dxig = NGHOST;
      dxir = NGH;
    } else if (nblevel < mylevel) {
      dxip = 2 * RINF;
      dxig = 2 * NGH;
      dxir = NG2;
    } else {
      dxip = RINF;
      dxig = dxir = NGHOST;
    }

    // Consider the normal directions.
    NeighborIndexes& ni = nb.ni;
    if (ni.ox1 > 0) {
      xi1min = xi1max - dxip;
      xi1_0 = xi1max;
      irs = ie - dxir + 1;
      iss = ie + 1;
      ise += NGHOST;
    } else if (ni.ox1 < 0) {
      xi1max = xi1min + dxip;
      xi1_0 = xi1min - dxig;
      ire = is + dxir - 1;
      iss -= NGHOST;
      ise = is - 1;
    }

    if (ni.ox2 > 0) {
      xi2min = xi2max - dxip;
      xi2_0 = xi2max;
      jrs = je - dxir + 1;
      jss = je + 1;
      jse += NGHOST;
    } else if (ni.ox2 < 0) {
      xi2max = xi2min + dxip;
      xi2_0 = xi2min - dxig;
      jre = js + dxir - 1;
      jss -= NGHOST;
      jse = js - 1;
    }

    if (ni.ox3 > 0) {
      xi3min = xi3max - dxip;
      xi3_0 = xi3max;
      krs = ke - dxir + 1;
      kss = ke + 1;
      kse += NGHOST;
    } else if (ni.ox3 < 0) {
      xi3max = xi3min + dxip;
      xi3_0 = xi3min - dxig;
      kre = ks + dxir - 1;
      kss -= NGHOST;
      kse = ks - 1;
    }

    // Consider the transverse directions.
    if (nblevel > mylevel) {  // Neighbor block is at a finer level.
      if (ni.type == NeighborConnect::face) {
        if (ni.ox1 != 0) {
          if (active2_) {
            if (ni.fi1) {
              xi2min = xi2mid - dxip;
              xi2_0 = xi2mid;
              jrs = je - nx2h + 1;
            } else {
              xi2max = xi2mid + dxip;
              jre = js + nx2h - 1;
            }
          }
          if (active3_) {
            if (ni.fi2) {
              xi3min = xi3mid - dxip;
              xi3_0 = xi3mid;
              krs = ke - nx3h + 1;
            } else {
              xi3max = xi3mid + dxip;
              kre = ks + nx3h - 1;
            }
          }
        } else if (ni.ox2 != 0) {
          if (active1_) {
            if (ni.fi1) {
              xi1min = xi1mid - dxip;
              xi1_0 = xi1mid;
              irs = ie - nx1h + 1;
            } else {
              xi1max = xi1mid + dxip;
              ire = is + nx1h - 1;
            }
          }
          if (active3_) {
            if (ni.fi2) {
              xi3min = xi3mid - dxip;
              xi3_0 = xi3mid;
              krs = ke - nx3h + 1;
            } else {
              xi3max = xi3mid + dxip;
              kre = ks + nx3h - 1;
            }
          }
        } else {
          if (active1_) {
            if (ni.fi1) {
              xi1min = xi1mid - dxip;
              xi1_0 = xi1mid;
              irs = ie - nx1h + 1;
            } else {
              xi1max = xi1mid + dxip;
              ire = is + nx1h - 1;
            }
          }
          if (active2_) {
            if (ni.fi2) {
              xi2min = xi2mid - dxip;
              xi2_0 = xi2mid;
              jrs = je - nx2h + 1;
            } else {
              xi2max = xi2mid + dxip;
              jre = js + nx2h - 1;
            }
          }
        }
      } else if (ni.type == NeighborConnect::edge) {
        if (ni.ox1 == 0) {
          if (active1_) {
            if (ni.fi1) {
              xi1min = xi1mid - dxip;
              xi1_0 = xi1mid;
              irs = ie - nx1h + 1;
            } else {
              xi1max = xi1mid + dxip;
              ire = is + nx1h - 1;
            }
          }
        } else if (ni.ox2 == 0) {
          if (active2_) {
            if (ni.fi1) {
              xi2min = xi2mid - dxip;
              xi2_0 = xi2mid;
              jrs = je - nx2h + 1;
            } else {
              xi2max = xi2mid + dxip;
              jre = js + nx2h - 1;
            }
          }
        } else {
          if (active3_) {
            if (ni.fi1) {
              xi3min = xi3mid - dxip;
              xi3_0 = xi3mid;
              krs = ke - nx3h + 1;
            } else {
              xi3max = xi3mid + dxip;
              kre = ks + nx3h - 1;
            }
          }
        }
      }
    } else if (nblevel < mylevel) {  // Neighbor block is at a coarser level.
      if (ni.type == NeighborConnect::face) {
        if (ni.ox1 != 0) {
          if (active2_ && myfx2) xi2_0 = xi2min - dxig;
          if (active3_ && myfx3) xi3_0 = xi3min - dxig;
        } else if (ni.ox2 != 0) {
          if (active1_ && myfx1) xi1_0 = xi1min - dxig;
          if (active3_ && myfx3) xi3_0 = xi3min - dxig;
        } else {
          if (active1_ && myfx1) xi1_0 = xi1min - dxig;
          if (active2_ && myfx2) xi2_0 = xi2min - dxig;
        }
      } else if (ni.type == NeighborConnect::edge) {
        if (ni.ox1 == 0) {
          if (active1_ && myfx1) xi1_0 = xi1min - dxig;
        } else if (ni.ox2 == 0) {
          if (active2_ && myfx2) xi2_0 = xi2min - dxig;
        } else {
          if (active3_ && myfx3) xi3_0 = xi3min - dxig;
        }
      }
    }

    // Set the domain that influences the ghost block.
    BoundaryAttributes& ba = ba_[n];
    ba.xi1min = xi1min;  ba.xi1max = xi1max;
    ba.xi2min = xi2min;  ba.xi2max = xi2max;
    ba.xi3min = xi3min;  ba.xi3max = xi3max;

    // Set the origin of the ghost block.
    ba.xi1_0 = xi1_0;
    ba.xi2_0 = xi2_0;
    ba.xi3_0 = xi3_0;

    // Set the dimensions of the ghost block.
    if (nblevel == mylevel) {
      ba.ngx1 = (ni.ox1 == 0) ? block_size.nx1 : NGHOST;
      ba.ngx2 = (ni.ox2 == 0) ? block_size.nx2 : NGHOST;
      ba.ngx3 = (ni.ox3 == 0) ? block_size.nx3 : NGHOST;
    } else if (nblevel < mylevel) {
      ba.ngx1 = (ni.ox1 == 0) ? nx1h : NGH;
      ba.ngx2 = (ni.ox2 == 0) ? nx2h : NGH;
      ba.ngx3 = (ni.ox3 == 0) ? nx3h : NGH;
    } else {
      ba.ngx1 = (ni.ox1 == 0) ? block_size.nx1 : NG2;
      ba.ngx2 = (ni.ox2 == 0) ? block_size.nx2 : NG2;
      ba.ngx3 = (ni.ox3 == 0) ? block_size.nx3 : NG2;
    }
    ba.ngtot = ba.ngx1 * ba.ngx2 * ba.ngx3;

    // Set the indices in meshaux to send and receive.
    ba.irs = irs;  ba.ire = ire;
    ba.jrs = jrs;  ba.jre = jre;
    ba.krs = krs;  ba.kre = kre;

    ba.iss = iss;  ba.ise = ise;
    ba.jss = jss;  ba.jse = jse;
    ba.kss = kss;  ba.kse = kse;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AssignParticlesToDifferentLevels(
//!        const std::vector<Real> par[], int p1, int ma1, int nprop)
//! \brief assigns particle array par from property index p1 to p1+nprop-1 to meshaux
//!        from property index ma1 to ma1+nprop-1 in neighbors of different levels.

void ParticleMesh::AssignParticlesToDifferentLevels(
         const std::vector<Real> par[], int p1, int ma1, int nprop) {
  const int mylevel = pmb_->loc.level;

  // Find neighbor blocks that are on a different level.
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    SimpleNeighborBlock& snb = nb.snb;
    if (snb.level == mylevel) continue;

    // Identify the buffer for assignment.
    Real *pbuf0 = NULL;
    BoundaryData<> *pnbd = NULL;
    if (nb.snb.rank == Globals::my_rank) {
      pnbd = &(pmesh_->FindMeshBlock(snb.gid)->ppar->ppm->bd_);
      pbuf0 = pnbd->recv[nb.targetid];
    }
#ifdef MPI_PARALLEL
    else
      pbuf0 = bd_.send[nb.bufid];
#endif

    // Zero out the buffer.
    BoundaryAttributes& ba = ba_[i];
    Real *pbufw = pbuf0 + iweight * ba.ngtot, *bufw = pbufw;
    for (int j = 0; j < ba.ngtot; ++j)
      *bufw++ = 0.0;

    Real **pbuf = new Real*[nprop];
    Real **buf = new Real*[nprop];
    for (int n = 0; n < nprop; ++n) {
      buf[n] = pbuf[n] = pbuf0 + (ma1 + n) * ba.ngtot;
      for (int j = 0; j < ba.ngtot; ++j)
        *buf[n]++ = 0.0;
    }

    // Find particles that influences the neighbor block.
    for (int k = 0; k < ppar_->npar; ++k) {
      Real xi1(ppar_->xi1[k]), xi2(ppar_->xi2[k]), xi3(ppar_->xi3[k]);
      if ((active1_ && (xi1 <= ba.xi1min || xi1 >= ba.xi1max)) ||
          (active2_ && (xi2 <= ba.xi2min || xi2 >= ba.xi2max)) ||
          (active3_ && (xi3 <= ba.xi3min || xi3 >= ba.xi3max))) continue;

      // Shift and scale the position index of the particle.
      xi1 -= ba.xi1_0;
      xi2 -= ba.xi2_0;
      xi3 -= ba.xi3_0;
      if (snb.level > mylevel) {
        xi1 *= 2;
        xi2 *= 2;
        xi3 *= 2;
      } else {
        xi1 /= 2;
        xi2 /= 2;
        xi3 /= 2;
      }

      // Find the region of the ghost block to assign the particle to.
      int ix1s = std::max(static_cast<int>(xi1 - dxi1_), 0),
          ix1e = std::min(static_cast<int>(xi1 + dxi1_), ba.ngx1-1),
          ix2s = std::max(static_cast<int>(xi2 - dxi2_), 0),
          ix2e = std::min(static_cast<int>(xi2 + dxi2_), ba.ngx2-1),
          ix3s = std::max(static_cast<int>(xi3 - dxi3_), 0),
          ix3e = std::min(static_cast<int>(xi3 + dxi3_), ba.ngx3-1);

      // Stack the properties of the particle and set the pointers to the buffer.
      Real *prop = new Real[nprop];
      int dbuf1 = ix1s + ba.ngx1 * (ix2s + ba.ngx2 * ix3s),
          dbuf2 = ba.ngx1 - ix1e + ix1s - 1,
          dbuf3 = ba.ngx1 * (ba.ngx2 - ix2e + ix2s - 1);

      bufw = pbufw + dbuf1;
      for (int n = 0; n < nprop; ++n) {
        prop[n] = par[p1+n][k];
        buf[n] = pbuf[n] + dbuf1;
      }

      // Assign the particle.
      for (int ix3 = ix3s; ix3 <= ix3e; ++ix3) {
        for (int ix2 = ix2s; ix2 <= ix2e; ++ix2) {
          for (int ix1 = ix1s; ix1 <= ix1e; ++ix1) {
            Real w = (active1_ ? _WeightFunction(ix1 + 0.5 - xi1) : 1.0) *
                     (active2_ ? _WeightFunction(ix2 + 0.5 - xi2) : 1.0) *
                     (active3_ ? _WeightFunction(ix3 + 0.5 - xi3) : 1.0);

            *bufw++ += w;
            for (int n = 0; n < nprop; ++n)
              *buf[n]++ += w * prop[n];
          }

          bufw += dbuf2;
          for (int n = 0; n < nprop; ++n)
            buf[n] += dbuf2;
        }

        bufw += dbuf3;
        for (int n = 0; n < nprop; ++n)
          buf[n] += dbuf3;
      }
      delete [] prop;
    }
    delete [] pbuf;
    delete [] buf;
  }
}

//--------------------------------------------------------------------------------------
//! \fn int ParticleMesh::LoadBoundaryBufferSameLevel(
//!                           Real *buf, const BoundaryAttributes& ba)
//! \brief Fill boundary buffers for sending to a block on the same level

int ParticleMesh::LoadBoundaryBufferSameLevel(Real *buf, const BoundaryAttributes& ba) {
  int p = 0;
  BufferUtility::PackData(meshaux, buf, 0, nmeshaux - 1,
                          ba.iss, ba.ise, ba.jss, ba.jse, ba.kss, ba.kse, p);
  return p;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AddBoundaryBuffer(Real *buf, const BoundaryAttributes& ba)
//! \brief Add boundary buffer from a neighbor block to meshaux.

void ParticleMesh::AddBoundaryBuffer(Real *buf, const BoundaryAttributes& ba) {
  // Add the data to the meshaux.
  for (int n = 0; n < nmeshaux; ++n)
    for (int k = ba.krs; k <= ba.kre; ++k)
      for (int j = ba.jrs; j <= ba.jre; ++j)
        for (int i = ba.irs; i <= ba.ire; ++i)
          meshaux(n,k,j,i) += *buf++;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::ClearBoundary()
//! \brief clears boundary data to neighboring blocks.

void ParticleMesh::ClearBoundary() {
  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank)
      MPI_Wait(&bd_.req_send[nb.bufid], MPI_STATUS_IGNORE);
#endif
    bd_.flag[nb.bufid] = BoundaryStatus::waiting;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::SendBoundary()
//! \brief Send boundary values to neighboring blocks.

void ParticleMesh::SendBoundary() {
  const int mylevel = pmb_->loc.level;

  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    SimpleNeighborBlock& snb = nb.snb;

    // Load boundary values.
    if (snb.rank == Globals::my_rank) {
      BoundaryData<> *pnbd = &(pmesh_->FindMeshBlock(snb.gid)->ppar->ppm->bd_);
      if (snb.level == mylevel)
        LoadBoundaryBufferSameLevel(pnbd->recv[nb.targetid], ba_[n]);
      pnbd->flag[nb.targetid] = BoundaryStatus::arrived;
    } else {
#ifdef MPI_PARALLEL
      if (snb.level == mylevel)
        LoadBoundaryBufferSameLevel(bd_.send[nb.bufid], ba_[n]);
      MPI_Start(&bd_.req_send[nb.bufid]);
#endif
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::StartReceiving()
//! \brief starts receiving meshaux near boundary from neighbor processes.

void ParticleMesh::StartReceiving() {
#ifdef MPI_PARALLEL
  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.snb.rank == Globals::my_rank) continue;
    MPI_Start(&bd_.req_recv[nb.bufid]);
  }
#endif
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::ReceiveBoundary()
//! \brief receives boundary values from neighboring blocks and add to my block and
//!        returns a flag indicating if all receives are completed.

#include <iomanip>

bool ParticleMesh::ReceiveBoundary() {
  bool completed = true;

  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    enum BoundaryStatus& bstatus = bd_.flag[nb.bufid];

    if (bstatus == BoundaryStatus::waiting) {
      if (nb.snb.rank == Globals::my_rank) {
        completed = false;
        continue;
      } else {
#ifdef MPI_PARALLEL
        int flag;
        MPI_Test(&bd_.req_recv[nb.bufid], &flag, MPI_STATUS_IGNORE);
        if (!flag) {
          completed = false;
          continue;
        }
        bstatus = BoundaryStatus::arrived;
#endif
      }
    }

    if (bstatus == BoundaryStatus::arrived) {
      AddBoundaryBuffer(bd_.recv[nb.bufid], ba_[n]);
      bstatus = BoundaryStatus::completed;
    }
  }

  return completed;
}

//--------------------------------------------------------------------------------------
//! \fn Real _WeightFunction(Real dxi)
//! \brief evaluates the weight function given index distance.

Real _WeightFunction(Real dxi) {
  dxi = std::min(std::abs(dxi), static_cast<Real>(1.5));
  return dxi < 0.5 ? 0.75 - dxi * dxi : 0.5 * ((1.5 - dxi) * (1.5 - dxi));
}
