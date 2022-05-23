//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//! \brief implementation of functions in class Photon

// C++ Standard Libraries
#include <vector>
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "photon.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// class variable initialization
bool Photon::initialized = false;
bool Photon::polarized = false;
bool Photon::general_mover_flag = false;

int Photon::inscp = -1, Photon::istatp = -1, Photon::itrp = -1;
int Photon::ii1p = -1, Photon::ii2p = -1, Photon::ii3p = -1;
int Photon::ix0p = -1, Photon::ix1p = -1, Photon::ix2p = -1, Photon::ix3p = -1;
int Photon::ik0p = -1, Photon::ik1p = -1, Photon::ik2p = -1, Photon::ik3p = -1;
int Photon::idk0p = -1, Photon::idk1p = -1, Photon::idk2p = -1, Photon::idk3p = -1;
int Photon::iep = -1, Photon::iwp = -1, Photon::iscp = -1, Photon::iacp = -1;
int Photon::isip = -1, Photon::isqp = -1, Photon::isup = -1, Photon::isvp = -1;
int Photon::iuserp = -1, Photon::ipolp = -1;

// Local function prototypes
static int CheckSide(int xi, int xi1, int xi2);
static int nloc = 0;
static int nper = 0;
static int nbuf = 0;
static int nnper = 0;
static int nadj = 0;
static int nmpi = 0;

//----------------------------------------------------------------------------------------
//! Photon constructor

Photon::Photon(MonteCarloBlock *pmcb, ParameterInput *pin)
  : Particles(pmcb->pmy_block, pin),
  // Allocate space for photon data via initialization list
    //user(new std::vector<Real> [pmcb->pmy_mc->nuser_var]),
    //polten(new std::vector<std::complex<Real>> [ncplx]),
    nphot(npar),nscp(intprop[inscp]), statp(intprop[istatp]),
    trp(intprop[itrp]), i1p(intprop[ii1p]), i2p(intprop[ii2p]), i3p(intprop[ii3p]),
    x0p(rp[ix0p]), x1p(rp[ix1p]), x2p(rp[ix2p]), x3p(rp[ix3p]),
    k0p(rp[ik0p]), k1p(rp[ik1p]), k2p(rp[ik2p]), k3p(rp[ik3p]),
    dk0p(rp[idk0p]), dk1p(rp[idk1p]), dk2p(rp[idk2p]),
    dk3p(rp[idk3p]),
    ep(rp[iep]), wp(rp[iwp]), scp(rp[iscp]), acp(rp[iacp]),
  sip(rp[isip]), sqp(rp[isqp]), sup(rp[isup]), svp(rp[isvp]) {

  pmy_mcb = pmcb;
  nphot_limit = pmcb->pmy_mc->max_phots_init;
  nuser_var = pmcb->pmy_mc->nuser_var;
  user = &(rp[iuserp]);
  polten = &(cplxprop[ipolp]);
  npar = 0;


}

//----------------------------------------------------------------------------------------
//! destructor

Photon::~Photon() {

}

//----------------------------------------------------------------------------------------
//! \fn void Photon::PrintPhoton(int ip)
//! \brief print key photon properites

void Photon::PrintPhoton(int ip) {
  // Used primarily for debugging
  std::cout << "----------------------------" << std::endl
            << "Energy, weight: " << ep[ip] << " " << wp[ip] << std::endl
            << "i: " << i1p[ip] << " " << i2p[ip] << " " << i3p[ip] <<std::endl
            << "x: " << x1p[ip] << " " << x2p[ip] << " " << x3p[ip] << " " << x0p[ip]
            << std::endl
            << "k: " << k1p[ip] << " " << k2p[ip] << " " << k3p[ip] << " " << k0p[ip]
            << std::endl;
  if (general_mover_flag) {
    std::cout << "dk: " << dk1p[ip] << " " << dk2p[ip] << " " << dk3p[ip] << " "
              << dk0p[ip] << std::endl;
  }
  if (polarized) {
    std:: cout << "stokes: " << sip[ip] << " " << sqp[ip] << " " << sup[ip] << std::endl
               << "opacity: " << scp[ip] << " " << acp[ip] << std::endl;
  }
  if (nuser_var > 0) {
    std::cout << "User vars:";
      for (int i=0; i<nuser_var; i++) {
        std::cout << " " << user[i][ip];
      }
      std::cout << std::endl;
  }
  if (statp[ip] == EVOLVING)
    std::cout << "EVOLVING" << std::endl;
  else if (statp[ip] == ESCAPED)
    std::cout << "ESCAPED" << std::endl;
  else if (statp[ip] == DESTROYED)
    std::cout << "DESTROYED" << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void Photon::IsNanPhoton(int ip)
//! \brief check for Nan in photon properties

bool Photon::IsNanPhoton(int ip) {

  if (isnan(wp[ip])) return true;
  if (isnan(ep[ip])) return true;
  if (isnan(x0p[ip])) return true;
  if (isnan(x1p[ip])) return true;
  if (isnan(x2p[ip])) return true;
  if (isnan(x3p[ip])) return true;
  if (isnan(k0p[ip])) return true;
  if (isnan(k1p[ip])) return true;
  if (isnan(k2p[ip])) return true;
  if (isnan(k3p[ip])) return true;
  if (polarized) {
    if (isnan(sip[ip])) return true;
    if (isnan(sqp[ip])) return true;
    if (isnan(sup[ip])) return true;
    if (isnan(svp[ip])) return true;
  }
  if (isnan(scp[ip])) return true;
  if (isnan(acp[ip])) return true;

  return false;
}

//--------------------------------------------------------------------------------------
//! \fn void Photon::AllocatePhotons(int nphot)
//! \brief Allocates photons

// SWD: Temporary --> converts protected function to public :(
void Photon::AllocatePhotons(int nphot) {
  // Call Resize function
  Resize(nphot);
}

//----------------------------------------------------------------------------------------
//! \fn void Photon::PolarizationToTetrad(std::complex<Real> ttet[4][4], Real ecov[4][4],
//!                                       const int ip)
//!
//! \brief transform complex tensor from coordinate frame to tetrad frame

void Photon::PolarizationToTetrad(std::complex<Real> ttet[4][4], Real ecov[4][4],
                                  const int ip) {

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      ttet[i][j] = std::complex<Real>(0.,0.);

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      for (int k = 0; k < 4; k++)
        for (int l = 0; l < 4; l++) {
          ttet[i][j] += polten[k*4+l][ip] * ecov[i][k] * ecov[j][l];
        }

}

//----------------------------------------------------------------------------------------
//! \fn void PolarizationToCoord(std::complex<Real> ttet[4][4], Real econ[4][4],
//!                              const int ip)
//!
//! \brief transform complex tensor from tetrad frame to coordinate frame

void Photon::PolarizationToCoord(std::complex<Real> ttet[4][4], Real econ[4][4],
                                 const int ip) {

  for(int i = 0; i < NCOORD; i++)
    for(int j = 0; j < NCOORD; j++)
      polten[i*4+j][ip] = std::complex<Real>(0.,0.);

  for(int i = 0; i < NCOORD; i++)
    for(int j = 0; j < NCOORD; j++)
      for(int k = 0; k < NCOORD; k++)
        for(int l = 0; l < NCOORD; l++) {
          polten[i*4+j][ip] += ttet[k][l] * econ[k][i] * econ[l][j];
        }

}

//--------------------------------------------------------------------------------------
//! \fn Photon::Initialize(MonteCarloBlock *pmcb, ParameterInput *pin)
//! \brief initializes the Photon class.
// SWD: Change name to distinguish with InitializePhoton?
void Photon::Initialize(MonteCarlo *pmc, ParameterInput *pin) {

  // Initialize first the parent class.
  Particles::Initialize(pmc->pmy_mesh, pin);

  if (initialized) return;

  // Add particle ID and status flags, other int parameters.
  inscp = AddIntProperty("nscp");
  istatp = AddIntProperty("statp");
  itrp = AddIntProperty("trp");

  // Add photon position.
  ix0p = AddRealProperty("x0");
  ix1p = AddRealProperty("x1");
  ix2p = AddRealProperty("x2");
  ix3p = AddRealProperty("x3");

  // Add photon momentum.
  ik0p = AddRealProperty("k0");
  ik1p = AddRealProperty("k1");
  ik2p = AddRealProperty("k2");
  ik3p = AddRealProperty("k3");

  if (pmc->general_mover_flag) {
    general_mover_flag = true;
    // Add change in photon momentum.
    idk0p = AddRealProperty("dk0");
    idk1p = AddRealProperty("dk1");
    idk2p = AddRealProperty("dk2");
    idk3p = AddRealProperty("dk3");
  }

  // Add energy, weight, and opacities
  iep = AddRealProperty("ep");
  iwp = AddRealProperty("wp");
  iscp = AddRealProperty("scp");
  iacp = AddRealProperty("acp");

  if (pmc->polarized) {
    polarized = true;
    // Add stokes vectors
    isip = AddRealProperty("sip");
    isqp = AddRealProperty("sqp");
    isup = AddRealProperty("sup");
    isvp = AddRealProperty("svp");
    if (general_mover_flag) {
      // Add complex polarization tensor
      for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
          int idummy = AddComplexProperty("pol"+std::to_string(i)+std::to_string(j));
          if ( (i==0) && (j==0))
            ipolp = idummy;
          }
      }
    }
  }

  // Add particle position indices.
  ii1p = AddIntProperty("i1p");
  ii2p = AddIntProperty("i2p");
  ii3p = AddIntProperty("i3p");

  // Add nuser variables
  for (int i=0; i<pmc->nuser_var; i++) {
    int idummy = AddRealProperty("user"+std::to_string(i));
    if (i==0)
      iuserp = idummy;
  }

#ifdef MPI_PARALLEL
  // Get my MPI communicator.
  MPI_Comm_dup(MPI_COMM_WORLD, &my_comm);
#endif

  initialized = true;
}


//--------------------------------------------------------------------------------------
//! \fn void Photon::SendToNeighbors()
//! \brief sends photons outside boundary to the buffers of neighboring meshblocks.

void Photon::SendToNeighbors() {
  const int IS = pmy_block->is;
  const int IE = pmy_block->ie;
  const int JS = pmy_block->js;
  const int JE = pmy_block->je;
  const int KS = pmy_block->ks;
  const int KE = pmy_block->ke;

  nbuf = nloc = nadj = nper = nnper = 0, nmpi = 0;
  for (int k = npar-1; k >=0; ) {
    if (statp[k] != BUFFERED) {
      --k;
      continue;
    }
    nbuf++;
    // Find which boundary photon has passed beyond
    int ox1 = CheckSide(i1p[k], IS, IE),
        ox2 = CheckSide(i2p[k], JS, JE),
        ox3 = CheckSide(i3p[k], KS, KE);
    if (ox1 == 0 && ox2 == 0 && ox3 == 0) {
      std::cout << "Warning: photon status is BUFFERED but not outside of boundary,"
                << " photon marked destroyed" << std::endl;
      statp[k] = DESTROYED;
      --k;
      continue;
    }

    // Apply periodic boundary conditions and find the mesh coordinates.
    ApplyPeriodicBoundary(x1p[k], x2p[k], x3p[k]);

    // Find the neighbor block to send it to.
    if (!active1_) ox1 = 0;
    if (!active2_) ox2 = 0;
    if (!active3_) ox3 = 0;
    Neighbor *pn = FindTargetNeighbor(ox1, ox2, ox3, i1p[k], i2p[k], i3p[k]);
    NeighborBlock *pnb = pn->pnb;
    if (pnb == nullptr) {
      PrintPhoton(k);
      RemoveOneParticle(k);
      --k;
      std::cout << "[SendToNeighbors] Warning: pnb==nullptr." << std::endl;
      continue;
    }

    // Determine which particle buffer to use.
    ParticleBuffer *ppb = NULL;
    if (pnb->snb.rank == Globals::my_rank) {
      // No need to send if back to the same block.
      if (pnb->snb.gid == pmy_block->gid) {
        Real xi1,xi2,xi3;
        pmy_block->pcoord->MeshCoordsToIndices(x1p[k], x2p[k], x3p[k], xi1, xi2, xi3);
        // Should be positive so no need for floor
        GetPositionIndices(k,k);
        //i1p[k] = static_cast<int>(xi1);
        //i2p[k] = static_cast<int>(xi2);
        //i3p[k] = static_cast<int>(xi3);
        //statp[k] = EVOLVING;
        --k;
        nloc++;
        continue;
      }
      // Use the target receive buffer.
      ppb = &pn->pmb->pmy_mcb->pphot->recv_[pnb->targetid];
      nadj++;
    } else {
#ifdef MPI_PARALLEL
      nmpi++;
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
      *pi++ = intprop[j][k];
    Real *pr(ppb->rbuf + ParticleBuffer::nreal * ppb->npar);
    for (int j = 0; j < nreal; ++j) {
      *pr++ = rp[j][k];
      *pr++ = rp1[j][k];
    }
    for (int j = 0; j < naux; ++j)
      *pr++ = aux[j][k];
    ++ppb->npar;
    // SWDNEW: ADD complex

    // Pop the particle from the current MeshBlock.
    RemoveOneParticle(k);
    --k;
  }

  // Send to neighbor processes and update boundary status.
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    int dst = nb.snb.rank;
    if (dst == Globals::my_rank) {
      Particles *ppar = pmy_mesh->FindMeshBlock(nb.snb.gid)->pmy_mcb->pphot;
      ppar->bstatus_[nb.targetid] =
          (ppar->recv_[nb.targetid].npar > 0) ? BoundaryStatus::arrived
                                              : BoundaryStatus::completed;
    } else {
#ifdef MPI_PARALLEL
      ParticleBuffer& send = send_[nb.bufid];
      int npsend = send.npar;
      MPI_Send(&npsend, 1, MPI_INT, nb.snb.rank, send.tag, my_comm);
      if (npsend > 0) {
        //printf("send: %d %d %d\n",Globals::my_rank,nb.snb.rank,npsend);
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
  printf("send %d %d %d %d %d %d %d\n",Globals::my_rank,nbuf,nloc,nadj,nmpi,nper,nnper);
}

//--------------------------------------------------------------------------------------
//! \fn void Photons::ApplyPeriodicBoundary(Real &x1, Real &x2, Real &x3)
//! \brief applies periodic boundary conditions to photon k and returns its updated mesh
//!        coordinates (x1,x2,x3).

void Photon::ApplyPeriodicBoundary(Real &x1, Real &x2, Real &x3) {
  bool flag = false;
  RegionSize& mesh_size = pmy_mesh->mesh_size;
  Coordinates *pcoord = pmy_block->pcoord;

  // Apply periodic boundary conditions in X1.
  if (x1 <= mesh_size.x1min) {
    // Inner x1
    x1 += mesh_size.x1len;
    flag = true;
  } else if (x1 >= mesh_size.x1max) {
    // Outer x1
    x1 -= mesh_size.x1len;
    flag = true;
  }

  // Apply periodic boundary conditions in X2.
  if (x2 <= mesh_size.x2min) {
    // Inner x2
    x2 += mesh_size.x2len;
    flag = true;
  } else if (x2 >= mesh_size.x2max) {
    // Outer x2
    x2 -= mesh_size.x2len;
    flag = true;
  }

  // Apply periodic boundary conditions in X3.
  if (x3 <= mesh_size.x3min) {
    // Inner x3
    x3 += mesh_size.x3len;
    flag = true;
  } else if (x3 >= mesh_size.x3max) {
    // Outer x3
    x3 -= mesh_size.x3len;
    flag = true;
  }
  if (flag) {
    nper++;
  } else {
    nnper++;
  }
}

//--------------------------------------------------------------------------------------
//! \fn bool Photon::ReceiveFromNeighbors()
//! \brief receives particles from neighboring meshblocks and returns a flag indicating
//!        if all receives are completed.

bool Photon::ReceiveFromNeighbors() {
  bool flag = true;

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    enum BoundaryStatus& bstatus = bstatus_[nb.bufid];

#ifdef MPI_PARALLEL
    // Communicate with neighbor processes.
    int nb_rank = nb.snb.rank;
    //printf("%d %d %d %d\n",Globals::my_rank,i,nb_rank,bstatus);
    if (nb_rank != Globals::my_rank && bstatus == BoundaryStatus::waiting) {
      ParticleBuffer& recv = recv_[nb.bufid];
      if (!recv.mpi_active) {
        // Get the number of incoming particles.
        MPI_Irecv(&recv.npar, 1, MPI_INT, nb_rank, recv.tag, my_comm, &recv.reqi);
        recv.mpi_active = true;
      }
      if (!recv.flagn) {
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
            // Receive data from the neighbor.
            MPI_Irecv(recv.ibuf, recv.npar * ParticleBuffer::nint, MPI_INT,
                      nb_rank, recv.tag + 1, my_comm, &recv.reqi);
            MPI_Irecv(recv.rbuf, recv.npar * ParticleBuffer::nreal, MPI_ATHENA_REAL,
                      nb_rank, recv.tag + 2, my_comm, &recv.reqr);
          } else {
            // No incoming particles.
            bstatus = BoundaryStatus::completed;
          }
        }
      }
      if (recv.flagn && recv.npar > 0) {
        if (!recv.flagi)
          MPI_Test(&recv.reqi, &recv.flagi, MPI_STATUS_IGNORE);
        if (!recv.flagr)
          MPI_Test(&recv.reqr, &recv.flagr, MPI_STATUS_IGNORE);
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
        int nparold = npar;
        FlushReceiveBuffer(recv);
        // Update Photon position indices
        GetPositionIndices(nparold,npar-1);
        //printf("recv %d %d %d\n",Globals::my_rank,nparold,npar-1);
        bstatus = BoundaryStatus::completed;
        break;
    }
  }

  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn void Photon::GetPositionIndices(int ibegin, int iend)
//! \brief finds the position indices of each particle with respect to the local grid.

void Photon::GetPositionIndices(int ibegin, int iend) {

  Real xi1, xi2, xi3;
  for (int k = ibegin; k <= iend; ++k) {
    // Convert to the index space.
    pmy_block->pcoord->MeshCoordsToIndices(x1p[k], x2p[k], x3p[k], xi1, xi2, xi3);
    i1p[k] = static_cast<int>(xi1);
    i2p[k] = static_cast<int>(xi2);
    i3p[k] = static_cast<int>(xi3);
    statp[k] = EVOLVING;
    MonteCarloBlock *pmcb = pmy_mcb;
    if (pmcb->boosts) {
      // Shift photon energy to comoving frame
      Real shift = pmcb->LorentzTransformFrequencyShift(this,k);
      ep[k] *= shift;
      // compute opacities in comoving frame
      acp[k] = pmcb->AbsorptionOpacity(pmcb,this,k);
      scp[k] = pmcb->ScatteringOpacity(pmcb,this,k);
      // Shift energy back to Eulerian frame
      ep[k] /= shift;
      // Shift opacities to Eulerian frame
      acp[k] *= shift;
      scp[k] *= shift;
    } else {
      // No distinction between comovinng frame and eulerian frame
      acp[k] = pmcb->AbsorptionOpacity(pmcb,this,k);
      scp[k] = pmcb->ScatteringOpacity(pmcb,this,k);
    }
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
