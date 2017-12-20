//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//=======================================================================================
//! \file particle.hpp
//  \brief defines particle classes
//  These classes contain data and functions related to particles
//=======================================================================================

// C++ headers
#include <iostream>   // endl
#include <stdexcept>  // runtime_error
#include <sstream>    // stringstream
#include <cstring>    // memcpy
#include <cstdlib>
#include <cmath>

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "particle.hpp"
#include "../bvals/bvals.hpp"
#include "../globals.hpp"

class BoundaryBase;

//! \class Particle
// \brief data/functions associated with list of particles

void Particle::WaitCounts()
{
  MeshBlock *pbl,*pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  int count_send, count_recv;
  int tag;

// exchange counts
#ifdef MPI_PARALLEL
  for (int n=0; n<pbval->nneighbor; n++){
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank){
      MPI_Wait(&req_cnt_recv[n],MPI_STATUS_IGNORE);
    }
  }
#endif
  return;
}


void Particle::WaitParticle()
{
  MeshBlock *pbl,*pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  int count_send, count_recv;
  int tag;

// exchange particles
#ifdef MPI_PARALLEL
  for (int n=0; n<pbval->nneighbor; n++) {
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank) { // not on the same process
      if (recv_cnt[n]!=0) MPI_Wait(&req_par_recv[n],MPI_STATUS_IGNORE);
    }
  }
#endif

  return;
}

void Particle::SendCounts()
{
  MeshBlock *pbl,*pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  int count_send, count_recv;
  int tag;

  for (int n=0; n<pbval->nneighbor; n++){
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank)
      req_cnt_send[n] = MPI_REQUEST_NULL;
  }

// exchange counts
  for (int n=0; n<pbval->nneighbor; n++) {
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank == Globals::my_rank) { // on the same process
      pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(&(pbl->particle->recv_cnt[nb.targetid]), &(send_cnt[nb.bufid]),sizeof(int));
    } else {
#ifdef MPI_PARALLEL
      tag = BoundaryBase::CreateBvalsMPITag(nb.lid, TAG_PRTL_CNT,nb.targetid);
      MPI_Isend(&send_cnt[n],1,MPI_INT,nb.rank,tag,MPI_COMM_WORLD,&req_cnt_send[n]);
#endif
    }
  }

  return;
}

void Particle::ReceiveCounts()
{
  MeshBlock *pbl,*pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  int count_send, count_recv;
  int tag;

  for (int n=0; n<pbval->nneighbor; n++){
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank)
      req_cnt_recv[n] = MPI_REQUEST_NULL;
  }

// exchange counts
  for (int n=0; n<pbval->nneighbor; n++) {
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank) { // not on the same process
#ifdef MPI_PARALLEL
      tag = BoundaryBase::CreateBvalsMPITag(pmb->lid, TAG_PRTL_CNT, nb.bufid);
      MPI_Irecv(&recv_cnt[n],1,MPI_INT,nb.rank,tag,MPI_COMM_WORLD,&req_cnt_recv[n]);
#endif
    }
  }

  return;
}


void Particle::CheckCnt()
{
  MeshBlock *pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;

#ifdef MPI_PARALLEL
  for (int n=0; n<pbval->nneighbor; n++){
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank){
      MPI_Wait(&req_cnt_send[n],MPI_STATUS_IGNORE);
      MPI_Wait(&req_cnt_recv[n],MPI_STATUS_IGNORE);
    }
  }
#endif


  return;
}

void Particle::SendParticle()
{
  MeshBlock *pbl,*pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  int count_recv, count_send;
  int tag;
  std::stringstream msg;
  for (int n=0; n<pbval->nneighbor; n++){
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank)
      req_par_send[n] = MPI_REQUEST_NULL;
  }

  for (int n=0; n<pbval->nneighbor; n++) {
    NeighborBlock& nb = pbval->neighbor[n];
    count_send=0;
    for (int i=0; i<n; i++)
      count_send+=6*send_cnt[i];
    if (nb.rank == Globals::my_rank) { // on the same process
      pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      count_recv=0;
      for (int i=0; i<nb.targetid; i++)
        count_recv+=6*(pbl->particle->recv_cnt[i]);
      // after calculating offsets, send/receive corresponding pointers
      std::memcpy(pbl->particle->recv_buf+count_recv, send_buf+count_send,send_cnt[nb.bufid]*6*sizeof(Real));
    }
#ifdef MPI_PARALLEL
    else {
     count_send=0;
      for (int i=0; i<n; i++)
        count_send+=6*send_cnt[i];
      tag = BoundaryBase::CreateBvalsMPITag(nb.lid, TAG_PARTICLE, nb.targetid);
      if (send_cnt[n]!=0) MPI_Isend(send_buf+count_send,6*send_cnt[n],MPI_ATHENA_REAL,nb.rank,tag,MPI_COMM_WORLD,&req_par_send[n]);
    }
#endif
  }

  return;
}

void Particle::ReceiveParticle()
{
  MeshBlock *pbl,*pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  int count_recv, count_send;
  int tag;
  std::stringstream msg;
  for (int n=0; n<pbval->nneighbor; n++){
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank)
      req_par_recv[n] = MPI_REQUEST_NULL;
  }

  for (int n=0; n<pbval->nneighbor; n++) {
    NeighborBlock& nb = pbval->neighbor[n];
    count_send=0;
    for (int i=0; i<n; i++)
      count_send+=6*send_cnt[i];
    if (nb.rank != Globals::my_rank) { // not on the same process
#ifdef MPI_PARALLEL
      count_recv=0;
      for (int i=0; i<n; i++)
        count_recv+=6*(recv_cnt[i]);
      tag = BoundaryBase::CreateBvalsMPITag(pmb->lid, TAG_PARTICLE, nb.bufid);
      if (recv_cnt[n]!=0) MPI_Irecv(recv_buf+count_recv,6*recv_cnt[n],MPI_ATHENA_REAL,nb.rank,tag,MPI_COMM_WORLD,&req_par_recv[n]);
    }
#endif
  }

  return;
}


void Particle::ExchangeParticle()
{
  MeshBlock *pbl,*pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  int count_recv, count_send;
  int tag;
  std::stringstream msg;
  for (int n=0; n<pbval->nneighbor; n++){
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank)
      req_par_recv[n] = MPI_REQUEST_NULL;
      req_par_send[n] = MPI_REQUEST_NULL;
  }

  for (int n=0; n<pbval->nneighbor; n++) {
    NeighborBlock& nb = pbval->neighbor[n];
    count_send=0;
    for (int i=0; i<n; i++)
      count_send+=6*send_cnt[i];
    if (nb.rank == Globals::my_rank) { // on the same process
      pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      count_recv=0;
      for (int i=0; i<nb.targetid; i++)
        count_recv+=6*(pbl->particle->recv_cnt[i]);
      // after calculating offsets, send/receive corresponding pointers
      std::memcpy(pbl->particle->recv_buf+count_recv, send_buf+count_send,send_cnt[nb.bufid]*6*sizeof(Real));
    }
#ifdef MPI_PARALLEL
    else {
      count_recv=0;
      for (int i=0; i<n; i++)
        count_recv+=6*(recv_cnt[i]);
      tag = BoundaryBase::CreateBvalsMPITag(pmb->lid, TAG_PARTICLE, nb.bufid);
      if (recv_cnt[n]!=0) MPI_Irecv(recv_buf+count_recv,6*recv_cnt[n],MPI_ATHENA_REAL,nb.rank,tag,MPI_COMM_WORLD,&req_par_recv[n]);   
    }
#endif
  }


  for (int n=0; n<pbval->nneighbor; n++) {
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank) { // not on the same process
      count_send=0;
      for (int i=0; i<n; i++)
        count_send+=6*send_cnt[i];
      tag = BoundaryBase::CreateBvalsMPITag(nb.lid, TAG_PARTICLE, nb.targetid);
      if (send_cnt[n]!=0) MPI_Isend(send_buf+count_send,6*send_cnt[n],MPI_ATHENA_REAL,nb.rank,tag,MPI_COMM_WORLD,&req_par_send[n]);
    }
  }

  for (int n=0; n<pbval->nneighbor; n++) {
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank) { // not on the same process
      if (recv_cnt[n]!=0) MPI_Wait(&req_par_recv[n],MPI_STATUS_IGNORE);
    }
  }

  return;
}

void Particle::ExchangeCounts()
{
  MeshBlock *pbl,*pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  int count_send, count_recv;
  int tag;

  for (int n=0; n<pbval->nneighbor; n++){
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank)
      req_cnt_send[n] = MPI_REQUEST_NULL;
      req_cnt_recv[n] = MPI_REQUEST_NULL;
  }

// exchange counts
  for (int n=0; n<pbval->nneighbor; n++) {
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank == Globals::my_rank) { // on the same process
      pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(&(pbl->particle->recv_cnt[nb.targetid]), &(send_cnt[nb.bufid]),sizeof(int));
    }
#ifdef MPI_PARALLEL
    else {
      tag = BoundaryBase::CreateBvalsMPITag(pmb->lid, TAG_PRTL_CNT, nb.bufid);
      MPI_Irecv(&recv_cnt[n],1,MPI_INT,nb.rank,tag,MPI_COMM_WORLD,&req_cnt_recv[n]);
    }
#endif
  }

#ifdef MPI_PARALLEL
  for (int n=0; n<pbval->nneighbor; n++){
    NeighborBlock& nb = pbval->neighbor[n];
    if (nb.rank != Globals::my_rank){
      tag = BoundaryBase::CreateBvalsMPITag(nb.lid, TAG_PRTL_CNT,nb.targetid);
      MPI_Isend(&send_cnt[n],1,MPI_INT,nb.rank,tag,MPI_COMM_WORLD,&req_cnt_send[n]);
    }
  }

#endif
  return;
}

void Particle::GetPos()
{
  MeshBlock *pmb=pmy_block;
  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;
  Real x2min = pmb->block_size.x2min;
  Real x2max = pmb->block_size.x2max;
  Real x3min = pmb->block_size.x3min;
  Real x3max = pmb->block_size.x3max;
  Real X1MIN = pmb->pmy_mesh->mesh_size.x1min;
  Real X1MAX = pmb->pmy_mesh->mesh_size.x1max;
  Real X2MIN = pmb->pmy_mesh->mesh_size.x2min;
  Real X2MAX = pmb->pmy_mesh->mesh_size.x2max;
  Real X3MIN = pmb->pmy_mesh->mesh_size.x3min;
  Real X3MAX = pmb->pmy_mesh->mesh_size.x3max;
  int ox1, ox2, ox3, if1, if2, if3;
  int fi1=0, fi2=0;
  bool flag = pmb->pmy_mesh->multilevel;
  for (int p=0; p<nparticle; p++) {
    ox1=-1; ox2=-1; ox3=-1;
    if(x1[p] > x1min) ox1++;
    if(x1[p] > x1max) ox1++;
    if(x2[p] > x2min) ox2++;
    if(x2[p] > x2max) ox2++;
    if(x3[p] > x3min) ox3++;
    if(x3[p] > x3max) ox3++;
    if(flag) {
      fi1=0; fi2=0;
      if1 = (x1[p] > x1min + 0.5*(x1max-x1min)) ? 1 : 0;
      if2 = (x2[p] > x2min + 0.5*(x2max-x2min)) ? 1 : 0;
      if3 = (x3[p] > x3min + 0.5*(x3max-x3min)) ? 1 : 0;
      if (ox1==0) {
        if (ox2==0) {
          if (ox3!=0) fi1=if1; fi2=if2;
        } else {
          if (ox3==0) {
            fi1=if1; fi2=if3;
          } else {
            fi1 = if1;
          }
        }
      } else {
        if (ox2==0) {
          if (ox3==0) {
            fi1=if2; fi2=if3;
          } else {
            fi1=if2;
          }
        } else {
          if (ox3==0) fi1=if3;
        }
      }
    }
   
    
    // periodic global boundaries
    if (x1[p] > X1MAX) {
      x1[p] = x1[p] - X1MAX + X1MIN;
    } else if (x1[p] < X1MIN) {
      x1[p] = x1[p] + X1MAX - X1MIN;
    } 
    if (x2[p] > X2MAX) {
      x2[p] = x2[p] - X2MAX + X2MIN;
    } else if (x2[p] < X2MIN) { 
      x2[p] = x2[p] + X2MAX - X2MIN;
    }
    if (x3[p] > X3MAX) {
      x3[p] = x3[p] - X3MAX + X3MIN;
    } else if (x3[p] < X3MIN) {
      x3[p] = x3[p] + X3MAX - X3MIN;
    }

    pos[p] = BoundaryBase::CreateBufferID(ox1,ox2,ox3,fi1,fi2);
  }

  return;
}

void Particle::PackParticle()
{
  MeshBlock *pmb = pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  int offset=0, flag;
  int skip = BoundaryBase::CreateBufferID(0,0,0,0,0);

  for (int n=0; n<pbval->nneighbor; n++) {
    send_cnt[n]=0;
    recv_cnt[n]=0;
  }      
  
  for (int n=0; n<pbval->nneighbor; n++) {
    for (int p=0; p<nparticle; p++) {
      flag=pos[p];
      if (flag==skip) {
        continue;
      } else {
        if (flag == pbval->bufid[n]) {
          send_cnt[n]++;
          *(send_buf+(offset++))=x1[p];
          *(send_buf+(offset++))=x2[p];
          *(send_buf+(offset++))=x3[p];
          *(send_buf+(offset++))=v1[p];
          *(send_buf+(offset++))=v2[p];
          *(send_buf+(offset++))=v3[p];
          nparticle--;
          x1[p]=x1[nparticle];
          x2[p]=x2[nparticle];
          x3[p]=x3[nparticle];
          v1[p]=v1[nparticle];
          v2[p]=v2[nparticle];
          v3[p]=v3[nparticle];
          pos[p]=pos[nparticle];
          p--;
        }
      }
    }
  }

  return;
}

void Particle::UnPackParticle()
{
  MeshBlock *pmb=pmy_block;
  BoundaryValues *pbval=pmb->pbval;
  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;
  Real x2min = pmb->block_size.x2min;
  Real x2max = pmb->block_size.x2max;
  Real x3min = pmb->block_size.x3min;
  Real x3max = pmb->block_size.x3max;

  int total=0;
  for (int i=0; i<pbval->nneighbor; i++) {
    total+=recv_cnt[i];
  }
  int offset=0;
  for (int p=0; p<total; p++) {
    x1[nparticle] = *(recv_buf+(offset++));
    x2[nparticle] = *(recv_buf+(offset++));
    x3[nparticle] = *(recv_buf+(offset++));
    v1[nparticle] = *(recv_buf+(offset++));
    v2[nparticle] = *(recv_buf+(offset++));
    v3[nparticle] = *(recv_buf+(offset++));
    pos[nparticle] = BoundaryBase::CreateBufferID(0,0,0,0,0);
    
    nparticle++;
  }
  return;
}

void Particle::Check()
{
  MeshBlock *pmb=pmy_block;
  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;
  Real x2min = pmb->block_size.x2min;
  Real x2max = pmb->block_size.x2max;
  Real x3min = pmb->block_size.x3min;
  Real x3max = pmb->block_size.x3max;
}


