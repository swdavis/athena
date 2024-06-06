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
bool Photon::general_pusher_flag = false;

int Photon::inscp = -1, Photon::istatp = -1, Photon::itrp = -1;
int Photon::ii1p = -1, Photon::ii2p = -1, Photon::ii3p = -1;
int Photon::ix0p = -1, Photon::ix1p = -1, Photon::ix2p = -1, Photon::ix3p = -1;
int Photon::ik0p = -1, Photon::ik1p = -1, Photon::ik2p = -1, Photon::ik3p = -1;
int Photon::idk0p = -1, Photon::idk1p = -1, Photon::idk2p = -1, Photon::idk3p = -1;
int Photon::iep = -1, Photon::iwp = -1, Photon::iscp = -1, Photon::iacp = -1;
int Photon::isip = -1, Photon::isqp = -1, Photon::isup = -1, Photon::isvp = -1;
int Photon::iuserp = -1, Photon::ipolp = -1, Photon::idtp = -1;

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
    sip(rp[isip]), sqp(rp[isqp]), sup(rp[isup]), svp(rp[isvp]),
    dtp(rp[idtp]) {

  pmy_mcb = pmcb;
  nphot_limit = pmcb->pmy_mc->max_phots_init;
  nuser_var = pmcb->pmy_mc->nuser_var;
  // SWD: should these be set or controlled by flags?
  user = &(rp[iuserp]);
  polten = &(cplxprop[ipolp]);
  npar = 0;


}

//----------------------------------------------------------------------------------------
//! destructor

Photon::~Photon() {

}

//----------------------------------------------------------------------------------------
//! \fn void Photon::PrintPhoton(std::stringstream msg, int ip)
//! \brief print key photon properites with message

void Photon::PrintPhoton(const std::string &msg, int ip) {
  std::cout << "----------------------------" << std::endl;
  std::cout << "** " << msg << " **" << std::endl;
  PrintPhoton(ip);
}

//----------------------------------------------------------------------------------------
//! \fn void Photon::PrintPhoton(int ip)
//! \brief print key photon properites, primarily for debugging

void Photon::PrintPhoton(int ip) {

  std::cout << "----------------------------" << std::endl
            << "Energy, weight: " << ep[ip] << " " << wp[ip] << std::endl
            << "i: " << i1p[ip] << " " << i2p[ip] << " " << i3p[ip] <<std::endl
            << "x: " << x1p[ip] << " " << x2p[ip] << " " << x3p[ip] << " " << x0p[ip]
            << std::endl
            << "k: " << k1p[ip] << " " << k2p[ip] << " " << k3p[ip] << " " << k0p[ip]
            << std::endl;
  if (general_pusher_flag) {
    std::cout << "dk: " << dk1p[ip] << " " << dk2p[ip] << " " << dk3p[ip] << " "
              << dk0p[ip] << std::endl;
  }
  if (polarized) {
    std:: cout << "stokes: " << sip[ip] << " " << sqp[ip] << " " << sup[ip] << std::endl;
    if (general_pusher_flag) {
      std:: cout << "pol tensor: ";
        for (int k = 0; k < 4; k++) {
          for (int l = 0; l < 4; l++) {
            std:: cout << polten[k*4+l][ip] << " ";
          }
          std::cout << std::endl;
        }
    }
  }
  std::cout << "opacity: " << scp[ip] << " " << acp[ip] << std::endl;
  if (nuser_var > 0) {
    std::cout << "User vars:";
      for (int i=0; i<nuser_var; i++) {
        std::cout << " " << user[i][ip];
      }
      std::cout << std::endl;
  }
  std::cout << "dt: " << dtp[ip] << " ";
  if (statp[ip] == EVOLVING)
    std::cout << "EVOLVING" << std::endl;
  else if (statp[ip] == ESCAPED)
    std::cout << "ESCAPED" << std::endl;
  else if (statp[ip] == DESTROYED)
    std::cout << "DESTROYED" << std::endl;
  else if (statp[ip] == BUFFERED)
    std::cout << "BUFFERED" << std::endl;
  else
    std::cout << std::endl;
  if (nuser_var > 0)
    std::cout << nuser_var << " user vars:";
  for (int i=0; i<nuser_var; i++)
    std::cout << " " << user[i][ip];
  std::cout << std::endl;
  std::cout << "----------------------------" << std::endl;

}

//----------------------------------------------------------------------------------------
//! \fn void Photon::IsNanPhoton(int ip)
//! \brief check for Nan in photon properties

bool Photon::IsNanPhoton(int ip) {

  if (std::isnan(wp[ip])) return true;
  if (std::isnan(ep[ip])) return true;
  if (std::isnan(x0p[ip])) return true;
  if (std::isnan(x1p[ip])) return true;
  if (std::isnan(x2p[ip])) return true;
  if (std::isnan(x3p[ip])) return true;
  if (std::isnan(k0p[ip])) return true;
  if (std::isnan(k1p[ip])) return true;
  if (std::isnan(k2p[ip])) return true;
  if (std::isnan(k3p[ip])) return true;
  if (polarized) {
    if (std::isnan(sip[ip])) return true;
    if (std::isnan(sqp[ip])) return true;
    if (std::isnan(sup[ip])) return true;
  }
  if (std::isnan(scp[ip])) return true;
  if (std::isnan(acp[ip])) return true;

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
// SWD: make trp a user variable

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

  if (pmc->general_pusher_flag) {
    general_pusher_flag = true;
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

  // Add time remaining parameter
  idtp = AddRealProperty("dtp");

  if (pmc->polarized) {
    polarized = true;
    // Add stokes vectors
    isip = AddRealProperty("sip");
    isqp = AddRealProperty("sqp");
    isup = AddRealProperty("sup");
    isvp = AddRealProperty("svp");
    if (general_pusher_flag) {
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
    ApplyPeriodicBoundary(x1p[k], x2p[k], x3p[k], k);
    //printf("%d %d %g %g %g\n",Globals::my_rank,k,x1p[k],x2p[k],x3p[k]);
    // Find the neighbor block to send it to.
    if (!active1_) ox1 = 0;
    if (!active2_) ox2 = 0;
    if (!active3_) ox3 = 0;
    Neighbor *pn = FindTargetNeighbor(ox1, ox2, ox3, i1p[k], i2p[k], i3p[k]);
    NeighborBlock *pnb = pn->pnb;
    if (pnb == nullptr) {
      PrintPhoton("pnb == nullptr",k);
      std::cout << ox1 << " " << ox2 << " " << ox3 << " " << i1p[k] << " "
                << i2p[k] << " " << i3p[k] << std::endl;
      MCCoord *pco = pmy_mcb->pcoord;
      //printf("%d %d %g %g %g %g %g %g\n",Globals::my_rank,k,
      //       pco->x1f(i1p[k]),pco->x1f(i1p[k]+1),
      //       pco->x2f(i2p[k]),pco->x2f(i2p[k]+1),
      //       pco->x3f(i3p[k]),pco->x3f(i3p[k]+1));

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
        GetPositionIndices(k,k);
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
    // copy complex properties
    if (general_pusher_flag && polarized) {
      std::complex<Real> *pc(ppb->cbuf + ParticleBuffer::ncplx * ppb->npar);
      for (int j = 0; j < ncplx; ++j) {
        *pc++ = cplxprop[j][k];
      }
    }
    ++ppb->npar;
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
      //if (npsend > 0)
      //printf("send: %d %d %d %d %d %d %g %g\n",Globals::my_rank,nb.snb.rank,pmy_block->lid,nb.targetid,send.tag+1,npsend,send.rbuf[0],send.rbuf[(npsend-1)*ParticleBuffer::nreal]);
      if (npsend > 0) {
        MPI_Request req = MPI_REQUEST_NULL;
	/*
        MPI_Isend(send.ibuf, npsend * ParticleBuffer::nint, MPI_INT,
                  dst, send.tag + 1, my_comm, &req);
        MPI_Request_free(&req);
        MPI_Isend(send.rbuf, npsend * ParticleBuffer::nreal, MPI_ATHENA_REAL,
                  dst, send.tag + 2, my_comm, &req);
        MPI_Request_free(&req);
	*/
        MPI_Isend(send.rbuf, npsend * ParticleBuffer::nreal, MPI_ATHENA_REAL,
                  dst, send.tag + 1, my_comm, &req);
        MPI_Request_free(&req);
        MPI_Isend(send.ibuf, npsend * ParticleBuffer::nint, MPI_INT,
                  dst, send.tag + 2, my_comm, &req);
        MPI_Request_free(&req);
        // Send complex properties
        if (general_pusher_flag && polarized) {
          MPI_Isend(send.cbuf, npsend * ParticleBuffer::ncplx, MPI_ATHENA_COMPLEX,
                    dst, send.tag + 3, my_comm, &req);
          MPI_Request_free(&req);
        }
      }
#endif
    }
  }
  //if (nbuf > 0)
  //   printf("send %d %d %d %d %d %d %d %d\n",Globals::my_rank,pmy_block->gid,nbuf,nloc,nadj,
  //         nmpi,nper,nnper);
}

//--------------------------------------------------------------------------------------
//! \fn void Photons::ApplyPeriodicBoundary(Real &x1, Real &x2, Real &x3, int k)
//! \brief applies periodic boundary conditions to photon k and returns its updated mesh
//!        coordinates (x1,x2,x3).

void Photon::ApplyPeriodicBoundary(Real &x1, Real &x2, Real &x3, int k) {
  bool flag = false;
  RegionSize& mesh_size = pmy_mesh->mesh_size;
  //MCCoord *pcoord = pmy_mcb->pcoord;
  Real l1cgs, l2cgs = 1., l3cgs = 1.;
  l1cgs = pmy_mcb->l_cgs;
  if ( (COORDINATE_SYSTEM == "cartesian") || (COORDINATE_SYSTEM == "minkowski") ) {
    l2cgs *= pmy_mcb->l_cgs;
    l3cgs *= pmy_mcb->l_cgs;
  }
  Real frac = 1.0e-8;

  // Apply periodic boundary conditions in X1.
  if (x1 <= mesh_size.x1min * l1cgs) {
    // Inner x1
    Real dx = l1cgs * mesh_size.x1min - x1;
    x1 = mesh_size.x1max * l1cgs - dx;
    //x1 = mesh_size.x1max * l1cgs * (1.-frac);
    flag = true;
  } else if (x1 >= mesh_size.x1max * l1cgs) {
    // Outer x1
    Real dx = x1 - l1cgs * mesh_size.x1max;
    x1 = mesh_size.x1min * l1cgs + dx;
    //x1 = mesh_size.x1min * l1cgs * (1.+frac);
    flag = true;
  }

  // Apply periodic boundary conditions in X2.
  if (x2 <= mesh_size.x2min * l2cgs) {
    // Inner x2
    Real dx = l2cgs * mesh_size.x2min - x2;
    x2 = mesh_size.x2max * l2cgs - dx;
    //x2 = mesh_size.x2max * l2cgs * (1.-frac);
    flag = true;
  } else if (x2 >= mesh_size.x2max * l2cgs) {
    // Outer x2
    Real dx = x2 - l2cgs * mesh_size.x2max;
    x2 = mesh_size.x2min * l2cgs + dx;
    //x2 = mesh_size.x2min * l2cgs * (1.+frac);
    flag = true;
  }

  // Apply periodic boundary conditions in X3.
  if (x3 <= mesh_size.x3min * l3cgs) {
    // Inner x3
    Real dx = l3cgs * mesh_size.x3min - x3;
    x3 = mesh_size.x3max * l3cgs - dx;
    //x3 = mesh_size.x3max * l3cgs * (1.-frac);
    flag = true;
  } else if (x3 >= mesh_size.x3max * l3cgs) {
    // Outer x3
    Real dx = x3 - l3cgs * mesh_size.x3max;
    x3 = mesh_size.x3min * l3cgs + dx;
    //x3 = mesh_size.x3min * l3cgs * (1.+frac);
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

    if (nb_rank != Globals::my_rank && bstatus == BoundaryStatus::waiting) {

      ParticleBuffer& recv = recv_[nb.bufid];
      if (!recv.mpi_active) {
        // Get the number of incoming particles.
        MPI_Irecv(&recv.npar, 1, MPI_INT, nb_rank, recv.tag, my_comm, &recv.reqn);
	//printf("r: %d %d %d %d %d %d %d\n",Globals::my_rank,nb_rank,nb.snb.lid,nb.bufid,recv_[nb.bufid].tag,recv.npar,recv.flagn);
        recv.mpi_active = true;
      }
      if (!recv.flagn) {
        MPI_Test(&recv.reqn, &recv.flagn, MPI_STATUS_IGNORE);
        if (recv.flagn) {
          if (recv.npar > 0) {
            // Check the buffer size.
            int nprecv = recv.npar;
            if (nprecv > recv.nparmax) {
	      //printf("buf res: %d %d %d %d %d %d %d\n",Globals::my_rank,nb_rank,nb.snb.lid,nb.bufid,recv_[nb.bufid].tag+1,recv.npar,recv.nparmax);
              recv.npar = 0;
	      //recv.Reallocate(2 * nprecv - recv.nparmax +50);
              recv.Reallocate(2 * nprecv - recv.nparmax);
              recv.npar = nprecv;
            }
            // Receive data from the neighbor.
	    /*
            MPI_Irecv(recv.ibuf, recv.npar * ParticleBuffer::nint, MPI_INT,
                      nb_rank, recv.tag + 1, my_comm, &recv.reqi);
            MPI_Irecv(recv.rbuf, recv.npar * ParticleBuffer::nreal, MPI_ATHENA_REAL,
                      nb_rank, recv.tag + 2, my_comm, &recv.reqr);
	    */
	    int ierr;
	    ierr = MPI_Irecv(recv.rbuf, recv.npar * ParticleBuffer::nreal, MPI_ATHENA_REAL,
                      nb_rank, recv.tag + 1, my_comm, &recv.reqr);

	    /*char err_buffer[MPI_MAX_ERROR_STRING];
	    int resultlen;
	    MPI_Error_string(ierr,err_buffer,&resultlen);
	    printf(err_buffer);
	    printf("\n");*/
            MPI_Irecv(recv.ibuf, recv.npar * ParticleBuffer::nint, MPI_INT,
                      nb_rank, recv.tag + 2, my_comm, &recv.reqi);
	    //int test;
	    //MPI_Status stat;
	    //MPI_Request_get_status(recv.reqr,&test,&stat);
	    //printf("t2: %d %d %d %d %d %d %d %d %d %d\n",Globals::my_rank,nb_rank,nb.snb.lid,nb.bufid,recv.tag+1,pmy_block->lid,test,stat.MPI_SOURCE,stat.MPI_TAG,stat.MPI_ERROR);
            if (general_pusher_flag && polarized) {
              MPI_Irecv(recv.cbuf, recv.npar * ParticleBuffer::ncplx, MPI_ATHENA_COMPLEX,
                        nb_rank, recv.tag + 3, my_comm, &recv.reqc);
            }
          } else {
            // No incoming particles.
            bstatus = BoundaryStatus::completed;
          }
        }
      }

      if (recv.flagn && recv.npar > 0) {
        if (!recv.flagi)
	  MPI_Test(&recv.reqi, &recv.flagi, MPI_STATUS_IGNORE);
        if (!recv.flagr) {
          MPI_Test(&recv.reqr, &recv.flagr, MPI_STATUS_IGNORE);
	}
        if (general_pusher_flag && polarized) {
          if (!recv.flagc)
            MPI_Test(&recv.reqc, &recv.flagc, MPI_STATUS_IGNORE);
          if (recv.flagi && recv.flagr && recv.flagc) {
            bstatus = BoundaryStatus::arrived;
	  }
        } else {
          if (recv.flagi && recv.flagr) {
            bstatus = BoundaryStatus::arrived;
	    //printf("g: %d %d %d %d %d\n",Globals::my_rank,nb_rank,nb.snb.lid,nb.bufid,recv_[nb.bufid].tag+1);
	  } else {
	    // SWD debug
	    //printf("bad: %d %d %d\n",recv.flagn,recv.flagi,recv.flagr);
	    //printf("%g %g %d %d\n",recv_[nb.bufid].rbuf[0],recv_[nb.bufid].rbuf[(recv_[nb.bufid].npar-1)*ParticleBuffer::nreal],recv_[nb.bufid].npar,recv_[nb.bufid].nparmax);
	    //printf("b: %d %d %d %d %d\n",Globals::my_rank,nb_rank,nb.snb.lid,nb.bufid,recv_[nb.bufid].tag+1);
	    //printf("n: %d\n",recv.flagn);
	    MPI_Wait(&recv.reqr, MPI_STATUS_IGNORE);
	    //printf("wait %d\n",Globals::my_rank);
	  }
        }
      }
    }
#endif
    switch (bstatus) {
      case BoundaryStatus::completed:
        break;

      case BoundaryStatus::waiting:
	//printf("w: %d %d %d %d %d %d\n",Globals::my_rank,nb_rank,nb.snb.lid,nb.bufid,recv_[nb.bufid].tag+1,recv_[nb.bufid].npar);
        flag = false;
        break;

      case BoundaryStatus::arrived:
        ParticleBuffer& recv = recv_[nb.bufid];
        int nparold = npar;
        FlushReceiveBuffer(recv);
        // Update Photon position indices
        GetPositionIndices(nparold,npar-1);
        //        printf("recv %d %d %d\n",Globals::my_rank,nparold,npar-1);
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
  int is = pmy_mcb->is, ie = pmy_mcb->ie;
  int js = pmy_mcb->js, je = pmy_mcb->je;
  int ks = pmy_mcb->ks, ke = pmy_mcb->ke;

  Real l1cgs, l2cgs = 1., l3cgs = 1.;
  l1cgs = pmy_mcb->l_cgs;
  if ( (COORDINATE_SYSTEM == "cartesian") || (COORDINATE_SYSTEM == "minkowski") ) {
    l2cgs *= pmy_mcb->l_cgs;
    l3cgs *= pmy_mcb->l_cgs;
  }

  for (int k = ibegin; k <= iend; ++k) {
    // Convert to the index space.
    pmy_block->pcoord->MeshCoordsToIndices(x1p[k]/l1cgs, x2p[k]/l2cgs, x3p[k]/l3cgs,
                                           xi1, xi2, xi3);

    int i1 = i1p[k] = static_cast<int>(xi1);
    //if (i1p[k] < is) {i1p[k] = is;}
    //if (i1p[k] > ie) {i1p[k] = ie;}
    if (i1p[k] = ie+1) {i1p[k] = ie;}
    int i2 = i2p[k] = static_cast<int>(xi2);
    //if (i2p[k] < js) i2p[k] = js;
    //if (i2p[k] > je) i2p[k] = je;
    if (i2p[k] = je+1) i2p[k] = je;
    int i3 =i3p[k] = static_cast<int>(xi3);
    //if (i3p[k] < ks) i3p[k] = ks;
    //if (i3p[k] > ke) i3p[k] = ke;
    if (i3p[k] = ke+1) i3p[k] = ke;

    // MeshCoordsToIndicies can fail for refined grids so check is needed
    MCCoord *pco = pmy_mcb->pcoord;
    bool on_block = true;
    while (x1p[k] > pco->x1f(i1p[k]+1)) {
      i1p[k]++;
      if (i1p[k] > ie) {
        on_block = false;
        break;
      }
    }
    while (x1p[k] < pco->x1f(i1p[k])) {
      i1p[k]--;
      if (i1p[k] < is) {
        on_block = false;
        break;
      }
    }
    while (x2p[k] > pco->x2f(i2p[k]+1)) {
      i2p[k]++;
      if (i2p[k] > je) {
        on_block = false;
        break;
      }
    }
    while (x2p[k] < pco->x2f(i2p[k])) {
      i2p[k]--;
      if (i2p[k] < js) {
        on_block = false;
        break;
      }
    }
    while (x3p[k] > pco->x3f(i3p[k]+1)) {
      i3p[k]++;
      if (i3p[k] > ke) {
        on_block = false;
        break;
      }
    }
    while (x3p[k] < pco->x3f(i3p[k])) {
      i3p[k]--;
      if (i3p[k] < ks) {
        on_block = false;
        break;
      }
    }
    if (on_block)
      statp[k] = EVOLVING;
    else {
      //printf("%d %d %d\n",i1,i2,i3);
      //printf("min/max: %g %g %g %g %g %g\n",pco->x1f(is),pco->x1f(ie+1),pco->x2f(js),
      //         pco->x2f(je+1),pco->x3f(ks),pco->x3f(ke+1));
    //printf("%d %d %d %d %d %d\n",is,ie,js,je,ks,ke);
      PrintPhoton("Warning: [GetPostionIndicies], Photon not on block, destroyed",k);
      statp[k] = DESTROYED;
      continue;
    }

    MonteCarloBlock *pmcb = pmy_mcb;
    if (pmcb->boosts) {
      // Shift photon energy to comoving frame
      //Real shift = pmcb->LorentzTransformFrequencyShift(this,k);
      Real shift = pmcb->FrequencyShiftComoving(this,k);
      if (( std::isinf(shift)) || (std::isnan(shift)) ) {
        printf("shift: %d %d %d %g %g %g %g\n",i1p[k],i2p[k],i3p[k],shift,xi1,xi2,xi3);
      }
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
    if (IsNanPhoton(k)) {
      PrintPhoton("Warning: Nan photon in GetPositionIndicies, destroying",k);
      statp[k] = DESTROYED;
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
