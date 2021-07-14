//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles.cpp
//! \brief implements functions in particle classes

// C++ Standard Libraries
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "particles.hpp"

// Class variable initialization
bool Particles::initialized = false;
int Particles::idmax = 0;
int Particles::nint = 0;
int Particles::nreal = 0;
int Particles::naux = 0;
int Particles::nwork = 0;
int Particles::ipid = -1;
int Particles::ixp = -1, Particles::iyp = -1, Particles::izp = -1;
int Particles::ivpx = -1, Particles::ivpy = -1, Particles::ivpz = -1;
int Particles::ixp0 = -1, Particles::iyp0 = -1, Particles::izp0 = -1;
int Particles::ivpx0 = -1, Particles::ivpy0 = -1, Particles::ivpz0 = -1;
int Particles::ixi1 = -1, Particles::ixi2 = -1, Particles::ixi3 = -1;
int Particles::imom1 = -1, Particles::imom2 = -1, Particles::imom3 = -1;
Real Particles::cfl_par = 1;
ParameterInput* Particles::pinput = NULL;
#ifdef MPI_PARALLEL
MPI_Comm Particles::my_comm = MPI_COMM_NULL;
#endif

// Local function prototypes
static int CheckSide(int xi, int xi1, int xi2);

//--------------------------------------------------------------------------------------
//! \fn void Particles::AMRCoarseToFine(MeshBlock* pmbc, MeshBlock* pmbf)
//! \brief load particles from a coarse meshblock to a fine meshblock.

void Particles::AMRCoarseToFine(MeshBlock* pmbc, MeshBlock* pmbf) {
  // Initialization
  Particles *pparc = pmbc->ppar, *pparf = pmbf->ppar;
  const Real x1min = pmbf->block_size.x1min, x1max = pmbf->block_size.x1max;
  const Real x2min = pmbf->block_size.x2min, x2max = pmbf->block_size.x2max;
  const Real x3min = pmbf->block_size.x3min, x3max = pmbf->block_size.x3max;
  const bool active1 = pparc->active1_,
             active2 = pparc->active2_,
             active3 = pparc->active3_;
  const AthenaArray<Real> &xp = pparc->xp, &yp = pparc->yp, &zp = pparc->zp;
  const Coordinates *pcoord = pmbf->pcoord;

  // Loop over particles in the coarse meshblock.
  for (int k = 0; k < pparc->npar; ++k) {
    Real x1, x2, x3;
    pcoord->CartesianToMeshCoords(xp(k), yp(k), zp(k), x1, x2, x3);
    if ((!active1 || (active1 && x1min <= x1 && x1 < x1max)) &&
        (!active2 || (active2 && x2min <= x2 && x2 < x2max)) &&
        (!active3 || (active3 && x3min <= x3 && x3 < x3max))) {
      // Load a particle to the fine meshblock.
      int npar = pparf->npar;
      if (npar >= pparf->nparmax) pparf->UpdateCapacity(2 * pparf->nparmax);
      for (int j = 0; j < nint; ++j)
        pparf->intprop(j,npar) = pparc->intprop(j,k);
      for (int j = 0; j < nreal; ++j)
        pparf->realprop(j,npar) = pparc->realprop(j,k);
      for (int j = 0; j < naux; ++j)
        pparf->auxprop(j,npar) = pparc->auxprop(j,k);
      ++pparf->npar;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::AMRFineToCoarse(MeshBlock* pmbf, MeshBlock* pmbc)
//! \brief load particles from a fine meshblock to a coarse meshblock.

void Particles::AMRFineToCoarse(MeshBlock* pmbf, MeshBlock* pmbc) {
  // Check the capacity.
  Particles *pparf = pmbf->ppar, *pparc = pmbc->ppar;
  int nparf = pparf->npar, nparc = pparc->npar;
  int npar_new = nparf + nparc;
  if (npar_new > pparc->nparmax) pparc->UpdateCapacity(npar_new);

  // Load the particles.
  for (int j = 0; j < nint; ++j)
    for (int k = 0; k < nparf; ++k)
      pparc->intprop(j,nparc+k) = pparf->intprop(j,k);
  for (int j = 0; j < nreal; ++j)
    for (int k = 0; k < nparf; ++k)
      pparc->realprop(j,nparc+k) = pparf->realprop(j,k);
  for (int j = 0; j < naux; ++j)
    for (int k = 0; k < nparf; ++k)
      pparc->auxprop(j,nparc+k) = pparf->auxprop(j,k);
  pparc->npar = npar_new;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::Initialize(Mesh *pm, ParameterInput *pin)
//! \brief initializes the class.

void Particles::Initialize(Mesh *pm, ParameterInput *pin) {
  if (initialized) return;

  // Add particle ID.
  ipid = AddIntProperty();

  // Add particle position.
  ixp = AddRealProperty();
  iyp = AddRealProperty();
  izp = AddRealProperty();

  // Add particle velocity.
  ivpx = AddRealProperty();
  ivpy = AddRealProperty();
  ivpz = AddRealProperty();

  // Add old particle position.
  ixp0 = AddAuxProperty();
  iyp0 = AddAuxProperty();
  izp0 = AddAuxProperty();

  // Add old particle velocity.
  ivpx0 = AddAuxProperty();
  ivpy0 = AddAuxProperty();
  ivpz0 = AddAuxProperty();

  // Add particle position indices.
  ixi1 = AddWorkingArray();
  ixi2 = AddWorkingArray();
  ixi3 = AddWorkingArray();

  // Initiate ParticleMesh class.
  ParticleMesh::Initialize(pin);
  imom1 = ParticleMesh::AddMeshAux();
  imom2 = ParticleMesh::AddMeshAux();
  imom3 = ParticleMesh::AddMeshAux();

  // Get the CFL number for particles.
  cfl_par = pin->GetOrAddReal("particles", "cfl_par", pm->cfl_number);

  // Remember the pointer to input parameters.
  pinput = pin;

#ifdef MPI_PARALLEL
  // Get my MPI communicator.
  MPI_Comm_dup(MPI_COMM_WORLD, &my_comm);
#endif

  initialized = true;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::PostInitialize(Mesh *pm, ParameterInput *pin)
//! \brief preprocesses the class after problem generator and before the main loop.

void Particles::PostInitialize(Mesh *pm, ParameterInput *pin) {
  // Set particle IDs.
  ProcessNewParticles(pm);

  // Set position indices.
  for (int b = 0; b < pm->nblocal; ++b)
    pm->my_blocks(b)->ppar->SetPositionIndices();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::FindDensityOnMesh(Mesh *pm, bool include_momentum)
//! \brief finds the number density of particles on the mesh.
//!
//!   If include_momentum is true, the momentum density field is also computed,
//!   assuming mass of each particle is unity.
//! \note
//!   Postcondition: ppm->weight becomes the density in each cell, and
//!   if include_momentum is true, ppm->meshaux(imom1:imom3,:,:,:)
//!   becomes the momentum density.

void Particles::FindDensityOnMesh(Mesh *pm, bool include_momentum) {
  // Assign particle properties to mesh and send boundary.
  int nblocks(pm->nblocal);
  for (int b = 0; b < nblocks; ++b) {
    const MeshBlock *pmb(pm->my_blocks(b));
    const Particles *ppar(pmb->ppar);
    ParticleMesh *ppm(ppar->ppm);
    ppm->StartReceiving();
    if (include_momentum) {
      AthenaArray<Real> vp, vp1, vp2, vp3;
      vp.NewAthenaArray(3, ppar->npar);
      vp1.InitWithShallowSlice(vp, 2, 0, 1);
      vp2.InitWithShallowSlice(vp, 2, 1, 1);
      vp3.InitWithShallowSlice(vp, 2, 2, 1);
      const Coordinates *pcoord = pmb->pcoord;
      for (int k = 0; k < ppar->npar; ++k)
        pcoord->CartesianToMeshCoordsVector(ppar->xp(k), ppar->yp(k), ppar->zp(k),
            ppar->vpx(k), ppar->vpy(k), ppar->vpz(k), vp1(k), vp2(k), vp3(k));
      ppm->AssignParticlesToMeshAux(vp, 0, imom1, 3);
    } else {
      ppm->AssignParticlesToMeshAux(ppar->realprop, 0, ppm->iweight, 0);
    }
    ppm->SendBoundary();
  }

  std::vector<bool> completed(nblocks, false);
  bool pending = true;
  while (pending) {
    pending = false;
    for (int i = 0; i < nblocks; ++i) {
      const MeshBlock *pmb(pm->my_blocks(i));
      Coordinates *pc(pmb->pcoord);
      ParticleMesh *ppm(pmb->ppar->ppm);
      if (!completed[i]) {
        // Finalize boundary communications.
        if ((completed[i] = ppm->ReceiveBoundary())) {
          // Convert to densities.
          const int is = ppm->is, ie = ppm->ie;
          const int js = ppm->js, je = ppm->je;
          const int ks = ppm->ks, ke = ppm->ke;
          if (include_momentum) {
            for (int k = ks; k <= ke; ++k)
              for (int j = js; j <= je; ++j)
                for (int i = is; i <= ie; ++i) {
                  Real vol(pc->GetCellVolume(k,j,i));
                  ppm->weight(k,j,i) /= vol;
                  ppm->meshaux(imom1,k,j,i) /= vol;
                  ppm->meshaux(imom2,k,j,i) /= vol;
                  ppm->meshaux(imom3,k,j,i) /= vol;
                }
          } else {
            for (int k = ks; k <= ke; ++k)
              for (int j = js; j <= je; ++j)
                for (int i = is; i <= ie; ++i)
                  ppm->weight(k,j,i) /= pc->GetCellVolume(k,j,i);
          }
          ppm->ClearBoundary();
        } else {
          pending = true;
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::FindHistoryOutput(Mesh *pm, Real data_sum[], int pos)
//! \brief finds the data sums of history output from particles in my process and assign
//!   them to data_sum beginning at index pos.

void Particles::FindHistoryOutput(Mesh *pm, Real data_sum[], int pos) {
  const int NSUM = NHISTORY - 1;

  // Initiate the summations.
  std::int64_t np = 0;
  std::vector<Real> sum(NSUM, 0.0);

  // Sum over each meshblock.
  Real vp1, vp2, vp3;
  for (int b = 0; b < pm->nblocal; ++b) {
    const MeshBlock *pmb(pm->my_blocks(b));
    const Particles *ppar(pmb->ppar);
    np += ppar->npar;
    const Coordinates *pcoord(pmb->pcoord);
    for (int k = 0; k < ppar->npar; ++k) {
      pcoord->CartesianToMeshCoordsVector(ppar->xp(k), ppar->yp(k), ppar->zp(k),
          ppar->vpx(k), ppar->vpy(k), ppar->vpz(k), vp1, vp2, vp3);
      sum[0] += vp1;
      sum[1] += vp2;
      sum[2] += vp3;
      sum[3] += vp1 * vp1;
      sum[4] += vp2 * vp2;
      sum[5] += vp3 * vp3;
    }
  }

  // Assign the values to output variables.
  data_sum[pos++] = static_cast<Real>(np);
  for (int i = 0; i < NSUM; ++i)
    data_sum[pos++] = sum[i];
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::GetHistoryOutputNames(std::string output_names[])
//! \brief gets the names of the history output variables in history_output_names[].

void Particles::GetHistoryOutputNames(std::string output_names[]) {
  output_names[0] = "np";
  output_names[1] = "vp1";
  output_names[2] = "vp2";
  output_names[3] = "vp3";
  output_names[4] = "vp1^2";
  output_names[5] = "vp2^2";
  output_names[6] = "vp3^2";
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::GetTotalNumber(Mesh *pm)
//! \brief returns total number of particles (from all processes).

int Particles::GetTotalNumber(Mesh *pm) {
  int npartot(0);
  for (int b = 0; b < pm->nblocal; ++b)
    npartot += pm->my_blocks(b)->ppar->npar;
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &npartot, 1, MPI_INT, MPI_SUM, my_comm);
#endif
  return npartot;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
//! \brief constructs a Particles instance.

Particles::Particles(MeshBlock *pmb, ParameterInput *pin) {
  // Point to the calling MeshBlock.
  pmy_block = pmb;
  pmy_mesh = pmb->pmy_mesh;
  pbval_ = pmb->pbval;
  nparmax = pin->GetOrAddInteger("particles", "nparmax", 1);
  npar = 0;

  // Check active dimensions.
  active1_ = pmy_mesh->mesh_size.nx1 > 1;
  active2_ = pmy_mesh->mesh_size.nx2 > 1;
  active3_ = pmy_mesh->mesh_size.nx3 > 1;

  // Allocate integer properties.
  intprop.NewAthenaArray(nint,nparmax);

  // Allocate integer properties.
  realprop.NewAthenaArray(nreal,nparmax);

  // Allocate auxiliary properties.
  if (naux > 0) auxprop.NewAthenaArray(naux,nparmax);

  // Allocate working arrays.
  if (nwork > 0) work.NewAthenaArray(nwork,nparmax);

  // Allocate mesh auxiliaries.
  ppm = new ParticleMesh(this);

  // Shallow copy to shorthands.
  AssignShorthands();

  // Initiate ParticleBuffer class.
  ParticleBuffer::SetNumberOfProperties(nint, nreal + naux);
}

//--------------------------------------------------------------------------------------
//! \fn Particles::~Particles()
//! \brief destroys a Particles instance.

Particles::~Particles() {
  // Delete integer properties.
  intprop.DeleteAthenaArray();

  // Delete real properties.
  realprop.DeleteAthenaArray();

  // Delete auxiliary properties.
  if (naux > 0) auxprop.DeleteAthenaArray();

  // Delete working arrays.
  if (nwork > 0) work.DeleteAthenaArray();

  // Clear links to neighbors.
  ClearNeighbors();

  // Delete mesh auxiliaries.
  delete ppm;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ClearBoundary()
//! \brief resets boundary for particle transportation.

void Particles::ClearBoundary() {
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    bstatus_[nb.bufid] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank) {
      ParticleBuffer& recv = recv_[nb.bufid];
      recv.flagn = recv.flagi = recv.flagr = 0;
      send_[nb.bufid].npar = 0;
    }
#endif
  }

  ppm->ClearBoundary();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ClearNeighbors()
//! \brief clears links to neighbors.

void Particles::ClearNeighbors() {
  delete neighbor_[1][1][1].pnb;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        Neighbor *pn = &neighbor_[i][j][k];
        if (pn == NULL) continue;
        while (pn->next != NULL)
          pn = pn->next;
        while (pn->prev != NULL) {
          pn = pn->prev;
          delete pn->next;
          pn->next = NULL;
        }
        pn->pnb = NULL;
        pn->pmb = NULL;
      }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Integrate(int step)
//! \brief updates all particle positions and velocities from t to t + dt.

void Particles::Integrate(int stage) {
  Real t = 0, dt = 0;

  // Determine the integration cofficients.
  switch (stage) {
  case 1:
    t = pmy_mesh->time;
    dt = 0.5 * pmy_mesh->dt;
    SaveStatus();
    break;

  case 2:
    t = pmy_mesh->time + 0.5 * pmy_mesh->dt;
    dt = pmy_mesh->dt;
    break;
  }

  // Conduct one stage of the integration.
  EulerStep(t, dt, pmy_block->phydro->w);
  ReactToMeshAux(t, dt, pmy_block->phydro->w);

  // Update the position index.
  SetPositionIndices();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::LinkNeighbors(MeshBlockTree &tree,
//!         int64_t nrbx1, int64_t nrbx2, int64_t nrbx3, int root_level)
//! \brief fetches neighbor information for later communication.

void Particles::LinkNeighbors(MeshBlockTree &tree,
    int64_t nrbx1, int64_t nrbx2, int64_t nrbx3, int root_level) {
  // Set myself as one of the neighbors.
  Neighbor *pn = &neighbor_[1][1][1];
  pn->pmb = pmy_block;
  pn->pnb = new NeighborBlock;
  pn->pnb->SetNeighbor(Globals::my_rank, pmy_block->loc.level,
      pmy_block->gid, pmy_block->lid, 0, 0, 0, NeighborConnect::none,
      -1, -1, false, false, 0, 0);

  // Save pointer to each neighbor.
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    SimpleNeighborBlock& snb = nb.snb;
    NeighborIndexes& ni = nb.ni;
    Neighbor *pn = &neighbor_[ni.ox1+1][ni.ox2+1][ni.ox3+1];
    while (pn->next != NULL)
      pn = pn->next;
    if (pn->pnb != NULL) {
      pn->next = new Neighbor;
      pn->next->prev = pn;
      pn = pn->next;
    }
    pn->pnb = &nb;
    if (snb.rank == Globals::my_rank) {
      pn->pmb = pmy_mesh->FindMeshBlock(snb.gid);
    } else {
#ifdef MPI_PARALLEL
      send_[nb.bufid].tag = (snb.gid<<8) | (nb.targetid<<2),
      recv_[nb.bufid].tag = (pmy_block->gid<<8) | (nb.bufid<<2);
#endif
    }
  }

  // Collect missing directions from fine to coarse level.
  if (pmy_mesh->multilevel) {
    int my_level = pbval_->loc.level;
    for (int l = 0; l < 3; l++) {
      if (!active1_ && l != 1) continue;
      for (int m = 0; m < 3; m++) {
        if (!active2_ && m != 1) continue;
        for (int n = 0; n < 3; n++) {
          if (!active3_ && n != 1) continue;
          Neighbor *pn = &neighbor_[l][m][n];
          if (pn->pnb == NULL) {
            int nblevel = pbval_->nblevel[n][m][l];
            if (0 <= nblevel && nblevel < my_level) {
              int ngid = tree.FindNeighbor(pbval_->loc, l-1, m-1, n-1)->GetGid();
              for (int i = 0; i < pbval_->nneighbor; ++i) {
                NeighborBlock& nb = pbval_->neighbor[i];
                if (nb.snb.gid == ngid) {
                  pn->pnb = &nb;
                  if (nb.snb.rank == Globals::my_rank)
                    pn->pmb = pmy_mesh->FindMeshBlock(ngid);
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  // Initiate ParticleMesh boundary data.
  ppm->SetBoundaryAttributes();
  ppm->InitiateBoundaryData();

  // Initiate boundary values.
  ClearBoundary();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::RemoveOneParticle(int k)
//! \brief removes particle k in the block.

void Particles::RemoveOneParticle(int k) {
  if (0 <= k && k < npar && --npar != k) {
    xi1(k) = xi1(npar);
    xi2(k) = xi2(npar);
    xi3(k) = xi3(npar);
    for (int j = 0; j < nint; ++j)
      intprop(j,k) = intprop(j,npar);
    for (int j = 0; j < nreal; ++j)
      realprop(j,k) = realprop(j,npar);
    for (int j = 0; j < naux; ++j)
      auxprop(j,k) = auxprop(j,npar);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendParticleMesh()
//! \brief send ParticleMesh meshaux near boundaries to neighbors.

void Particles::SendParticleMesh() {
  if (ppm->nmeshaux > 0)
    ppm->SendBoundary();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendToNeighbors()
//! \brief sends particles outside boundary to the buffers of neighboring meshblocks.

void Particles::SendToNeighbors() {
  const int IS = pmy_block->is;
  const int IE = pmy_block->ie;
  const int JS = pmy_block->js;
  const int JE = pmy_block->je;
  const int KS = pmy_block->ks;
  const int KE = pmy_block->ke;

  for (int k = 0; k < npar; ) {
    // Check if a particle is outside the boundary.
    int xi1i = static_cast<int>(xi1(k)),
        xi2i = static_cast<int>(xi2(k)),
        xi3i = static_cast<int>(xi3(k));
    int ox1 = CheckSide(xi1i, IS, IE),
        ox2 = CheckSide(xi2i, JS, JE),
        ox3 = CheckSide(xi3i, KS, KE);
    if (ox1 == 0 && ox2 == 0 && ox3 == 0) {
      ++k;
      continue;
    }

    // Apply boundary conditions and find the mesh coordinates.
    Real x1, x2, x3;
    ApplyBoundaryConditions(k, x1, x2, x3);

    // Find the neighbor block to send it to.
    if (!active1_) ox1 = 0;
    if (!active2_) ox2 = 0;
    if (!active3_) ox3 = 0;
    Neighbor *pn = FindTargetNeighbor(ox1, ox2, ox3, xi1i, xi2i, xi3i);
    NeighborBlock *pnb = pn->pnb;
    if (pnb == NULL) {
      RemoveOneParticle(k);
      continue;
    }

    // Determine which particle buffer to use.
    ParticleBuffer *ppb = NULL;
    if (pnb->snb.rank == Globals::my_rank) {
      // No need to send if back to the same block.
      if (pnb->snb.gid == pmy_block->gid) {
        pmy_block->pcoord->MeshCoordsToIndices(x1, x2, x3, xi1(k), xi2(k), xi3(k));
        ++k;
        continue;
      }
      // Use the target receive buffer.
      ppb = &pn->pmb->ppar->recv_[pnb->targetid];

    } else {
#ifdef MPI_PARALLEL
      // Use the send buffer.
      ppb = &send_[pnb->bufid];
#endif
    }

    // Check the buffer size.
    if (ppb->npar >= ppb->nparmax)
      ppb->Reallocate((ppb->nparmax > 0) ? 2 * ppb->nparmax : 1);

    // Copy the properties of the particle to the buffer.
    int *pi = ppb->ibuf + ParticleBuffer::nint * ppb->npar;
    for (int j = 0; j < nint; ++j)
      *pi++ = intprop(j,k);
    Real *pr = ppb->rbuf + ParticleBuffer::nreal * ppb->npar;
    for (int j = 0; j < nreal; ++j)
      *pr++ = realprop(j,k);
    for (int j = 0; j < naux; ++j)
      *pr++ = auxprop(j,k);
    ++ppb->npar;

    // Pop the particle from the current MeshBlock.
    RemoveOneParticle(k);
  }

  // Send to neighbor processes and update boundary status.
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    int dst = nb.snb.rank;
    if (dst == Globals::my_rank) {
      Particles *ppar = pmy_mesh->FindMeshBlock(nb.snb.gid)->ppar;
      ppar->bstatus_[nb.targetid] =
          (ppar->recv_[nb.targetid].npar > 0) ? BoundaryStatus::arrived
                                              : BoundaryStatus::completed;
    } else {
#ifdef MPI_PARALLEL
      ParticleBuffer& send = send_[nb.bufid];
      int npsend = send.npar;
      MPI_Send(&npsend, 1, MPI_INT, nb.snb.rank, send.tag, my_comm);
      if (npsend > 0) {
        MPI_Request req = MPI_REQUEST_NULL;
        MPI_Isend(send.ibuf, npsend * ParticleBuffer::nint, MPI_INT,
                  dst, send.tag + 1, my_comm, &req);
        MPI_Request_free(&req);
        MPI_Isend(send.rbuf, npsend * ParticleBuffer::nreal, MPI_ATHENA_REAL,
                  dst, send.tag + 2, my_comm, &req);
        MPI_Request_free(&req);
      }
#endif
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SetPositionIndices()
//! \brief updates position indices of particles.

void Particles::SetPositionIndices() {
  GetPositionIndices(npar, xp, yp, zp, xi1, xi2, xi3);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::StartReceiving()
//! \brief starts receiving ParticleMesh meshaux near boundary from neighbor processes.

void Particles::StartReceiving() {
  ppm->StartReceiving();
}

//--------------------------------------------------------------------------------------
//! \fn bool Particles::ReceiveFromNeighbors()
//! \brief receives particles from neighboring meshblocks and returns a flag indicating
//!        if all receives are completed.

bool Particles::ReceiveFromNeighbors() {
  bool flag = true;

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    enum BoundaryStatus& bstatus = bstatus_[nb.bufid];

#ifdef MPI_PARALLEL
    // Communicate with neighbor processes.
    int nb_rank = nb.snb.rank;
    if (nb_rank != Globals::my_rank && bstatus == BoundaryStatus::waiting) {
      ParticleBuffer& recv = recv_[nb.bufid];
      if (!recv.flagn) {
        // Get the number of incoming particles.
        if (recv.reqi == MPI_REQUEST_NULL)
          MPI_Irecv(&recv.npar, 1, MPI_INT, nb_rank, recv.tag, my_comm, &recv.reqi);
        else
          MPI_Test(&recv.reqi, &recv.flagn, MPI_STATUS_IGNORE);
        if (recv.flagn) {
          if (recv.npar > 0) {
            // Check the buffer size.
            int nprecv = recv.npar;
            if (nprecv > recv.nparmax) {
              recv.npar = 0;
              recv.Reallocate(2 * nprecv - recv.nparmax);
              recv.npar = nprecv;
            }
          } else {
            // No incoming particles.
            bstatus = BoundaryStatus::completed;
          }
        }
      }
      if (recv.flagn && recv.npar > 0) {
        // Receive data from the neighbor.
        if (!recv.flagi) {
          if (recv.reqi == MPI_REQUEST_NULL)
            MPI_Irecv(recv.ibuf, recv.npar * ParticleBuffer::nint, MPI_INT,
                      nb_rank, recv.tag + 1, my_comm, &recv.reqi);
          else
            MPI_Test(&recv.reqi, &recv.flagi, MPI_STATUS_IGNORE);
        }
        if (!recv.flagr) {
          if (recv.reqr == MPI_REQUEST_NULL)
            MPI_Irecv(recv.rbuf, recv.npar * ParticleBuffer::nreal, MPI_ATHENA_REAL,
                      nb_rank, recv.tag + 2, my_comm, &recv.reqr);
          else
            MPI_Test(&recv.reqr, &recv.flagr, MPI_STATUS_IGNORE);
        }
        if (recv.flagi && recv.flagr)
          bstatus = BoundaryStatus::arrived;
      }
    }
#endif

    switch (bstatus) {
      case BoundaryStatus::completed:
        break;

      case BoundaryStatus::waiting:
        flag = false;
        break;

      case BoundaryStatus::arrived:
        ParticleBuffer& recv = recv_[nb.bufid];
        FlushReceiveBuffer(recv);
        bstatus = BoundaryStatus::completed;
        break;
    }
  }

  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn bool Particles::ReceiveParticleMesh(int step)
//! \brief receives ParticleMesh meshaux near boundaries from neighbors and returns a
//!        flag indicating if all receives are completed.

bool Particles::ReceiveParticleMesh(int stage) {
  if (ppm->nmeshaux <= 0) return true;

  // Flush ParticleMesh receive buffers.
  bool flag = ppm->ReceiveBoundary();

  if (flag) {
    // Deposit ParticleMesh meshaux to MeshBlock.
    Hydro *phydro = pmy_block->phydro;
    Real t = 0, dt = 0;

    switch (stage) {
    case 1:
      t = pmy_mesh->time;
      dt = 0.5 * pmy_mesh->dt;
      break;

    case 2:
      t = pmy_mesh->time + 0.5 * pmy_mesh->dt;
      dt = pmy_mesh->dt;
      break;
    }

    DepositToMesh(t, dt, phydro->w, phydro->u);
  }

  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ProcessNewParticles()
//! \brief searches for and books new particles.

void Particles::ProcessNewParticles(Mesh *pmesh) {
  // Count new particles.
  const int nbtotal(pmesh->nbtotal), nblocks(pmesh->nblocal);
  std::vector<int> nnewpar(nbtotal, 0);
  for (int b = 0; b < nblocks; ++b) {
    const MeshBlock *pmb(pmesh->my_blocks(b));
    nnewpar[pmb->gid] = pmb->ppar->CountNewParticles();
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &nnewpar[0], nbtotal, MPI_INT, MPI_MAX, my_comm);
#endif

  // Make the counts cumulative.
  for (int i = 1; i < nbtotal; ++i)
    nnewpar[i] += nnewpar[i-1];

  // Set particle IDs.
  for (int b = 0; b < nblocks; ++b) {
    const MeshBlock *pmb(pmesh->my_blocks(b));
    pmb->ppar->SetNewParticleID(idmax + (pmb->gid > 0 ? nnewpar[pmb->gid - 1] : 0));
  }
  idmax += nnewpar[nbtotal - 1];
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::CountNewParticles()
//! \brief counts new particles in the block.

int Particles::CountNewParticles() const {
  int n = 0;
  for (int i = 0; i < npar; ++i)
    if (pid(i) <= 0) ++n;
  return n;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3)
//! \brief applies boundary conditions to particle k and returns its updated mesh
//!        coordinates (x1,x2,x3).
//! \todo (ccyang):
//! - implement nonperiodic boundary conditions.

void Particles::ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3) {
  bool flag = false;
  RegionSize& mesh_size = pmy_mesh->mesh_size;
  Coordinates *pcoord = pmy_block->pcoord;

  // Find the mesh coordinates.
  Real x10, x20, x30;
  pcoord->IndicesToMeshCoords(xi1(k), xi2(k), xi3(k), x1, x2, x3);
  pcoord->CartesianToMeshCoords(xp0(k), yp0(k), zp0(k), x10, x20, x30);

  // Convert velocity vectors in mesh coordinates.
  Real vp1, vp2, vp3, vp10, vp20, vp30;
  pcoord->CartesianToMeshCoordsVector(xp(k), yp(k), zp(k),
                                      vpx(k), vpy(k), vpz(k), vp1, vp2, vp3);
  pcoord->CartesianToMeshCoordsVector(xp0(k), yp0(k), zp0(k),
                                      vpx0(k), vpy0(k), vpz0(k), vp10, vp20, vp30);

  // Apply periodic boundary conditions in X1.
  if (x1 < mesh_size.x1min) {
    // Inner x1
    x1 += mesh_size.x1len;
    x10 += mesh_size.x1len;
    flag = true;
  } else if (x1 >= mesh_size.x1max) {
    // Outer x1
    x1 -= mesh_size.x1len;
    x10 -= mesh_size.x1len;
    flag = true;
  }

  // Apply periodic boundary conditions in X2.
  if (x2 < mesh_size.x2min) {
    // Inner x2
    x2 += mesh_size.x2len;
    x20 += mesh_size.x2len;
    flag = true;
  } else if (x2 >= mesh_size.x2max) {
    // Outer x2
    x2 -= mesh_size.x2len;
    x20 -= mesh_size.x2len;
    flag = true;
  }

  // Apply periodic boundary conditions in X3.
  if (x3 < mesh_size.x3min) {
    // Inner x3
    x3 += mesh_size.x3len;
    x30 += mesh_size.x3len;
    flag = true;
  } else if (x3 >= mesh_size.x3max) {
    // Outer x3
    x3 -= mesh_size.x3len;
    x30 -= mesh_size.x3len;
    flag = true;
  }

  if (flag) {
    // Convert positions and velocities back in Cartesian coordinates.
    pcoord->MeshCoordsToCartesian(x1, x2, x3, xp(k), yp(k), zp(k));
    pcoord->MeshCoordsToCartesian(x10, x20, x30, xp0(k), yp0(k), zp0(k));
    pcoord->MeshCoordsToCartesianVector(x1, x2, x3,
                                        vp1, vp2, vp3, vpx(k), vpy(k), vpz(k));
    pcoord->MeshCoordsToCartesianVector(x10, x20, x30,
                                        vp10, vp20, vp30, vpx0(k), vpy0(k), vpz0(k));
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc)
//! \brief evolves the particle positions and velocities by one Euler step.

void Particles::EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  // Update positions.
  for (int k = 0; k < npar; ++k) {
    //! \todo (ccyang):
    //! - This is a temporary hack.
    Real tmpx = xp(k), tmpy = yp(k), tmpz = zp(k);
    xp(k) = xp0(k) + dt * vpx(k);
    yp(k) = yp0(k) + dt * vpy(k);
    zp(k) = zp0(k) + dt * vpz(k);
    xp0(k) = tmpx;
    yp0(k) = tmpy;
    zp0(k) = tmpz;
  }

  // Integrate the source terms (e.g., acceleration).
  SourceTerms(t, dt, meshsrc);
  UserSourceTerms(t, dt, meshsrc);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::GetPositionIndices(int npar,
//!                                        const AthenaArray<Real>& xp,
//!                                        const AthenaArray<Real>& yp,
//!                                        const AthenaArray<Real>& zp,
//!                                        AthenaArray<Real>& xi1,
//!                                        AthenaArray<Real>& xi2,
//!                                        AthenaArray<Real>& xi3)
//! \brief finds the position indices of each particle with respect to the local grid.

void Particles::GetPositionIndices(int npar,
                                   const AthenaArray<Real>& xp,
                                   const AthenaArray<Real>& yp,
                                   const AthenaArray<Real>& zp,
                                   AthenaArray<Real>& xi1,
                                   AthenaArray<Real>& xi2,
                                   AthenaArray<Real>& xi3) {
  for (int k = 0; k < npar; ++k) {
    // Convert to the Mesh coordinates.
    Real x1, x2, x3;
    pmy_block->pcoord->CartesianToMeshCoords(xp(k), yp(k), zp(k), x1, x2, x3);

    // Convert to the index space.
    pmy_block->pcoord->MeshCoordsToIndices(x1, x2, x3, xi1(k), xi2(k), xi3(k));
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SetNewParticleID(int id0)
//! \brief searches for new particles and assigns ID, beginning at id + 1.

void Particles::SetNewParticleID(int id) {
  for (int i = 0; i < npar; ++i)
    if (pid(i) <= 0) pid(i) = ++id;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SaveStatus()
//! \brief saves the current positions and velocities for later use.

void Particles::SaveStatus() {
  for (int k = 0; k < npar; ++k) {
    // Save current positions.
    xp0(k) = xp(k);
    yp0(k) = yp(k);
    zp0(k) = zp(k);

    // Save current velocities.
    vpx0(k) = vpx(k);
    vpy0(k) = vpy(k);
    vpz0(k) = vpz(k);
  }
}

//--------------------------------------------------------------------------------------
//! \fn MeshBlock* Particles::FindTargetNeighbor(
//!         int ox1, int ox2, int ox3, int xi1, int xi2, int xi3)
//! \brief finds the neighbor to send a particle to.

struct Neighbor* Particles::FindTargetNeighbor(
    int ox1, int ox2, int ox3, int xi1, int xi2, int xi3) {
  // Find the head of the linked list.
  Neighbor *pn = &neighbor_[ox1+1][ox2+1][ox3+1];

  // Search down the list if the neighbor is at a finer level.
  if (pmy_mesh->multilevel && pn->pnb != NULL &&
      pn->pnb->snb.level > pmy_block->loc.level) {
    RegionSize& bs = pmy_block->block_size;
    int fi[2] = {0, 0}, i = 0;
    if (active1_ && ox1 == 0) fi[i++] = 2 * (xi1 - pmy_block->is) / bs.nx1;
    if (active2_ && ox2 == 0) fi[i++] = 2 * (xi2 - pmy_block->js) / bs.nx2;
    if (active3_ && ox3 == 0) fi[i++] = 2 * (xi3 - pmy_block->ks) / bs.nx3;
    while (pn != NULL) {
      NeighborIndexes& ni = pn->pnb->ni;
      if (ni.fi1 == fi[0] && ni.fi2 == fi[1]) break;
      pn = pn->next;
    }
  }

  // Return the target neighbor.
  return pn;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::FlushReceiveBuffer(ParticleBuffer& recv)
//! \brief adds particles from the receive buffer.

void Particles::FlushReceiveBuffer(ParticleBuffer& recv) {
  // Check the memory size.
  int nprecv = recv.npar;
  if (npar + nprecv > nparmax)
    UpdateCapacity(nparmax + 2 * (npar + nprecv - nparmax));

  // Flush the receive buffers.
  int *pi = recv.ibuf;
  Real *pr = recv.rbuf;
  for (int k = npar; k < npar + nprecv; ++k) {
    for (int j = 0; j < nint; ++j)
      intprop(j,k) = *pi++;
    for (int j = 0; j < nreal; ++j)
      realprop(j,k) = *pr++;
    for (int j = 0; j < naux; ++j)
      auxprop(j,k) = *pr++;
  }

  // Find their position indices.
  AthenaArray<Real> xps, yps, zps, xi1s, xi2s, xi3s;
  xps.InitWithShallowSlice(xp, 1, npar, nprecv);
  yps.InitWithShallowSlice(yp, 1, npar, nprecv);
  zps.InitWithShallowSlice(zp, 1, npar, nprecv);
  xi1s.InitWithShallowSlice(xi1, 1, npar, nprecv);
  xi2s.InitWithShallowSlice(xi2, 1, npar, nprecv);
  xi3s.InitWithShallowSlice(xi3, 1, npar, nprecv);
  GetPositionIndices(nprecv, xps, yps, zps, xi1s, xi2s, xi3s);

  // Clear the receive buffers.
  npar += nprecv;
  recv.npar = 0;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddIntProperty()
//! \brief adds one integer property to the particles and returns the index.

int Particles::AddIntProperty() {
  return nint++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddRealProperty()
//! \brief adds one real property to the particles and returns the index.

int Particles::AddRealProperty() {
  return nreal++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddAuxProperty()
//! \brief adds one auxiliary property to the particles and returns the index.

int Particles::AddAuxProperty() {
  return naux++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddWorkingArray()
//! \brief adds one working array to the particles and returns the index.

int Particles::AddWorkingArray() {
  return nwork++;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::AssignShorthands()
//! \brief assigns shorthands by shallow copying slices of the data.

void Particles::AssignShorthands() {
  pid.InitWithShallowSlice(intprop, 2, ipid, 1);

  xp.InitWithShallowSlice(realprop, 2, ixp, 1);
  yp.InitWithShallowSlice(realprop, 2, iyp, 1);
  zp.InitWithShallowSlice(realprop, 2, izp, 1);
  vpx.InitWithShallowSlice(realprop, 2, ivpx, 1);
  vpy.InitWithShallowSlice(realprop, 2, ivpy, 1);
  vpz.InitWithShallowSlice(realprop, 2, ivpz, 1);

  xp0.InitWithShallowSlice(auxprop, 2, ixp0, 1);
  yp0.InitWithShallowSlice(auxprop, 2, iyp0, 1);
  zp0.InitWithShallowSlice(auxprop, 2, izp0, 1);
  vpx0.InitWithShallowSlice(auxprop, 2, ivpx0, 1);
  vpy0.InitWithShallowSlice(auxprop, 2, ivpy0, 1);
  vpz0.InitWithShallowSlice(auxprop, 2, ivpz0, 1);

  xi1.InitWithShallowSlice(work, 2, ixi1, 1);
  xi2.InitWithShallowSlice(work, 2, ixi2, 1);
  xi3.InitWithShallowSlice(work, 2, ixi3, 1);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::UpdateCapacity(int new_nparmax)
//! \brief changes the capacity of particle arrays while preserving existing data.

void Particles::UpdateCapacity(int new_nparmax) {
  // Increase size of property arrays
  nparmax = new_nparmax;
  intprop.ResizeLastDimension(nparmax);
  realprop.ResizeLastDimension(nparmax);
  if (naux > 0) auxprop.ResizeLastDimension(nparmax);
  if (nwork > 0) work.ResizeLastDimension(nparmax);

  // Reassign the shorthands.
  AssignShorthands();
}

//--------------------------------------------------------------------------------------
//! \fn Real Particles::NewBlockTimeStep();
//! \brief returns the time step required by particles in the block.

Real Particles::NewBlockTimeStep() {
  Coordinates *pc = pmy_block->pcoord;

  // Find the maximum coordinate speed.
  Real dt_inv2_max = 0.0;
  for (int k = 0; k < npar; ++k) {
    Real dt_inv2 = 0.0, vpx1, vpx2, vpx3;
    pc->CartesianToMeshCoordsVector(xp(k), yp(k), zp(k), vpx(k), vpy(k), vpz(k),
                                    vpx1, vpx2, vpx3);
    dt_inv2 += active1_ ? std::pow(vpx1 / pc->dx1f(static_cast<int>(xi1(k))), 2) : 0;
    dt_inv2 += active2_ ? std::pow(vpx2 / pc->dx2f(static_cast<int>(xi2(k))), 2) : 0;
    dt_inv2 += active3_ ? std::pow(vpx3 / pc->dx3f(static_cast<int>(xi3(k))), 2) : 0;
    dt_inv2_max = std::max(dt_inv2_max, dt_inv2);
  }

  // Return the time step constrained by the coordinate speed.
  return dt_inv2_max > 0.0 ? cfl_par / std::sqrt(dt_inv2_max)
                           : std::numeric_limits<Real>::max();
}

//--------------------------------------------------------------------------------------
//! \fn std::size_t Particles::GetSizeInBytes()
//! \brief returns the data size in bytes in the meshblock.

std::size_t Particles::GetSizeInBytes() {
  std::size_t size = sizeof(npar);
  if (npar > 0) size += npar * (nint * sizeof(int) + nreal * sizeof(Real));
  return size;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::UnpackParticlesForRestart()
//! \brief reads the particle data from the restart file.

void Particles::UnpackParticlesForRestart(char *mbdata, std::size_t &os) {
  // Read number of particles.
  std::memcpy(&npar, &(mbdata[os]), sizeof(npar));
  os += sizeof(npar);
  if (nparmax < npar)
    UpdateCapacity(npar);

  if (npar > 0) {
    // Read integer properties.
    std::size_t size = npar * sizeof(int);
    for (int k = 0; k < nint; ++k) {
      std::memcpy(&(intprop(k,0)), &(mbdata[os]), size);
      os += size;
    }

    // Read real properties.
    size = npar * sizeof(Real);
    for (int k = 0; k < nreal; ++k) {
      std::memcpy(&(realprop(k,0)), &(mbdata[os]), size);
      os += size;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::PackParticlesForRestart()
//! \brief pack the particle data for restart dump.

void Particles::PackParticlesForRestart(char *&pdata) {
  // Write number of particles.
  std::memcpy(pdata, &npar, sizeof(npar));
  pdata += sizeof(npar);

  if (npar > 0) {
    // Write integer properties.
    std::size_t size = npar * sizeof(int);
    for (int k = 0; k < nint; ++k) {
      std::memcpy(pdata, &(intprop(k,0)), size);
      pdata += size;
    }
    // Write real properties.
    size = npar * sizeof(Real);
    for (int k = 0; k < nreal; ++k) {
      std::memcpy(pdata, &(realprop(k,0)), size);
      pdata += size;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::FormattedTableOutput()
//! \brief outputs the particle data in tabulated format.

void Particles::FormattedTableOutput(Mesh *pm, OutputParameters op) {
  std::stringstream fname, msg;
  std::ofstream os;

  // Loop over MeshBlocks
  for (int b = 0; b < pm->nblocal; ++b) {
    const MeshBlock *pmb(pm->my_blocks(b));
    const Particles *ppar(pmb->ppar);

    // Create the filename.
    fname << op.file_basename
          << ".block" << pmb->gid << '.' << op.file_id
          << '.' << std::setw(5) << std::right << std::setfill('0') << op.file_number
          << '.' << "par.tab";

    // Open the file for write.
    os.open(fname.str().data());
    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [Particles::FormattedTableOutput]"
          << std::endl << "Output file '" << fname.str() << "' could not be opened"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    // Write the time.
    os << std::scientific << std::showpoint << std::setprecision(18);
    os << "# Athena++ particle data at time = " << pm->time << std::endl;

    // Write the particle data in the meshblock.
    for (int k = 0; k < ppar->npar; ++k)
      os << ppar->pid(k) << "  "
         << ppar->xp(k) << "  " << ppar->yp(k) << "  " << ppar->zp(k) << "  "
         << ppar->vpx(k) << "  " << ppar->vpy(k) << "  " << ppar->vpz(k) << std::endl;

    // Close the file and get the next meshblock.
    os.close();
    fname.str("");
  }
}

//--------------------------------------------------------------------------------------
//! \fn int CheckSide(int xi, nx, int xi1, int xi2)
//! \brief returns -1 if xi < xi1, +1 if xi > xi2, or 0 otherwise.

inline int CheckSide(int xi, int xi1, int xi2) {
  if (xi < xi1) return -1;
  if (xi > xi2) return +1;
  return 0;
}
