//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file evals.cpp
//  \brief constructor/destructor and utility functions for ExchangeValues class

// C++ headers
#include <iostream>   // endl
#include <iomanip>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memcpy
#include <cstdlib>
#include <cmath>

// Athena++ classes headers
#include "bvals.hpp"
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh_refinement.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../multigrid/multigrid.hpp"
#include "../gravity/mggravity.hpp"
#include "../coordinates/coordinates.hpp"
#include "../utils/buffer_utils.hpp"
#include "../hybrid/hybrid.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// ExchangeValues constructor - sets functions for the appropriate
// boundary conditions at each of the 6 dirs of a MeshBlock

ExchangeValues::ExchangeValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs)
 : BoundaryBase(pmb->pmy_mesh, pmb->loc, pmb->block_size, input_bcs)
{
  pmy_block_=pmb;
  for(int i=0; i<6; i++)
    ExchangeFunction_[i]=NULL;

// Set exchange functions for each of the 6 boundaries in turn ---------------------------------
  // Inner x1
  nface_=2; nedge_=0;
  switch(block_bcs[INNER_X1]){
    case REFLECTING_BNDRY:
      ExchangeFunction_[INNER_X1] = ReflectInnerExchangeX1;
      break;
    case OUTFLOW_BNDRY:
      ExchangeFunction_[INNER_X1] = OutflowInnerExchangeX1;
      break;
    case BLOCK_BNDRY: // block boundary
    case PERIODIC_BNDRY: // periodic boundary
      ExchangeFunction_[INNER_X1] = NULL;
      break;
    case USER_BNDRY: // user-enrolled BCs
      ExchangeFunction_[INNER_X1] = pmy_mesh_->ExchangeFunction_[INNER_X1];
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in ExchangeValues constructor" << std::endl
          << "Flag ix1_bc=" << block_bcs[INNER_X1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
      break;
   }

  // Outer x1
  switch(block_bcs[OUTER_X1]){
    case REFLECTING_BNDRY:
      ExchangeFunction_[OUTER_X1] = ReflectOuterExchangeX1;
      break;
    case OUTFLOW_BNDRY:
      ExchangeFunction_[OUTER_X1] = OutflowOuterExchangeX1;
      break;
    case BLOCK_BNDRY: // block boundary
    case PERIODIC_BNDRY: // periodic boundary
      ExchangeFunction_[OUTER_X1] = NULL;
      break;
    case USER_BNDRY: // user-enrolled BCs
      ExchangeFunction_[OUTER_X1] = pmy_mesh_->ExchangeFunction_[OUTER_X1];
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in ExchangeValues constructor" << std::endl
          << "Flag ox1_bc=" << block_bcs[OUTER_X1] << " not valid" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

  if (pmb->block_size.nx2 > 1) {
    nface_=4; nedge_=4;
    // Inner x2
    switch(block_bcs[INNER_X2]){
      case REFLECTING_BNDRY:
        ExchangeFunction_[INNER_X2] = ReflectInnerExchangeX2;
        break;
      case OUTFLOW_BNDRY:
        ExchangeFunction_[INNER_X2] = OutflowInnerExchangeX2;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
      case POLAR_BNDRY: // polar boundary
        ExchangeFunction_[INNER_X2] = NULL;
        break;
      case POLAR_BNDRY_WEDGE: //polar boundary with a wedge
        ExchangeFunction_[INNER_X2] = PolarWedgeInnerExchangeX2;
        break;
      case USER_BNDRY: // user-enrolled BCs
        ExchangeFunction_[INNER_X2] = pmy_mesh_->ExchangeFunction_[INNER_X2];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in ExchangeValues constructor" << std::endl
            << "Flag ix2_bc=" << block_bcs[INNER_X2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

    // Outer x2
    switch(block_bcs[OUTER_X2]){
      case REFLECTING_BNDRY:
        ExchangeFunction_[OUTER_X2] = ReflectOuterExchangeX2;
        break;
      case OUTFLOW_BNDRY:
        ExchangeFunction_[OUTER_X2] = OutflowOuterExchangeX2;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
      case POLAR_BNDRY: // polar boundary
        ExchangeFunction_[OUTER_X2] = NULL;
        break;
      case POLAR_BNDRY_WEDGE: //polar boundary with a wedge
        ExchangeFunction_[OUTER_X2] = PolarWedgeOuterExchangeX2;
        break;
      case USER_BNDRY: // user-enrolled BCs
        ExchangeFunction_[OUTER_X2] = pmy_mesh_->ExchangeFunction_[OUTER_X2];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in ExchangeValues constructor" << std::endl
            << "Flag ox2_bc=" << block_bcs[OUTER_X2] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  if (pmb->block_size.nx3 > 1) {
    nface_=6; nedge_=12;
    // Inner x3
    switch(block_bcs[INNER_X3]){
      case REFLECTING_BNDRY:
        ExchangeFunction_[INNER_X3] = ReflectInnerExchangeX3;
        break;
      case OUTFLOW_BNDRY:
        ExchangeFunction_[INNER_X3] = OutflowInnerExchangeX3;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
        ExchangeFunction_[INNER_X3] = NULL;
        break;
      case USER_BNDRY: // user-enrolled BCs
        ExchangeFunction_[INNER_X3] = pmy_mesh_->ExchangeFunction_[INNER_X3];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in ExchangeValues constructor" << std::endl
            << "Flag ix3_bc=" << block_bcs[INNER_X3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
     }

    // Outer x3
    switch(block_bcs[OUTER_X3]){
      case REFLECTING_BNDRY:
        ExchangeFunction_[OUTER_X3] = ReflectOuterExchangeX3;
        break;
      case OUTFLOW_BNDRY:
        ExchangeFunction_[OUTER_X3] = OutflowOuterExchangeX3;
        break;
      case BLOCK_BNDRY: // block boundary
      case PERIODIC_BNDRY: // periodic boundary
        ExchangeFunction_[OUTER_X3] = NULL;
        break;
      case USER_BNDRY: // user-enrolled BCs
        ExchangeFunction_[OUTER_X3] = pmy_mesh_->ExchangeFunction_[OUTER_X3];
        break;
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in ExchangeValues constructor" << std::endl
            << "Flag ox3_bc=" << block_bcs[OUTER_X3] << " not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
  }

  // Count number of blocks wrapping around pole
  if (block_bcs[INNER_X2] == POLAR_BNDRY || block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) {
    if(pmy_mesh_->nrbx3>1 && pmy_mesh_->nrbx3%2!=0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ExchangeValues constructor" << std::endl
          << "Number of MeshBlocks around the pole must be 1 or even." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    int level = pmb->loc.level - pmy_mesh_->root_level;
    num_north_polar_blocks_ = pmy_mesh_->nrbx3 * (1 << level);
  }
  else
    num_north_polar_blocks_ = 0;
  if (block_bcs[OUTER_X2] == POLAR_BNDRY || block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) {
    if(pmy_mesh_->nrbx3>1 && pmy_mesh_->nrbx3%2!=0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ExchangeValues constructor" << std::endl
          << "Number of MeshBlocks around the pole must be 1 or even." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    int level = pmb->loc.level - pmy_mesh_->root_level;
    num_south_polar_blocks_ = pmy_mesh_->nrbx3 * (1 << level);
  } else {
    num_south_polar_blocks_ = 0;
  }

  if (PARTICLE)
    InitExchangeData(bd_mcoup_, BNDRY_MCOUP);

 /* single CPU in the azimuthal direction with the polar boundary*/
  if(pmb->loc.level == pmy_mesh_->root_level &&
     pmy_mesh_->nrbx3 == 1 &&
     (block_bcs[INNER_X2]==POLAR_BNDRY||block_bcs[OUTER_X2]==POLAR_BNDRY||
      block_bcs[INNER_X2]==POLAR_BNDRY_WEDGE||block_bcs[OUTER_X2]==POLAR_BNDRY_WEDGE))
       exc_.NewAthenaArray(pmb->ke+NGHOST+2);

  SearchAndSetNeighbors(pmy_mesh_->tree, pmy_mesh_->ranklist, pmy_mesh_->nslist);
}

// destructor

ExchangeValues::~ExchangeValues()
{
  MeshBlock *pmb=pmy_block_;
  if (PARTICLE) {
    DestroyExchangeData(bd_mcoup_);
  }
  if(pmb->loc.level == pmy_mesh_->root_level &&
     pmy_mesh_->nrbx3 == 1 &&
     (block_bcs[INNER_X2]==POLAR_BNDRY||block_bcs[OUTER_X2]==POLAR_BNDRY||
      block_bcs[INNER_X2]==POLAR_BNDRY_WEDGE||block_bcs[OUTER_X2]==POLAR_BNDRY_WEDGE))
       exc_.DeleteAthenaArray();
}


//----------------------------------------------------------------------------------------
//! \fn void ExchangeValues::InitExchangeData(BoundaryData &bd, enum BoundaryType type)
//  \brief Initialize BoundaryData structure
void ExchangeValues::InitExchangeData(BoundaryData &bd, enum BoundaryType type)
{
  MeshBlock *pmb=pmy_block_;
  bool multilevel=pmy_mesh_->multilevel;
  int cng=pmb->cnghost, cng1=0, cng2=0, cng3=0;
  if(pmb->block_size.nx2>1) cng1=cng, cng2=cng;
  if(pmb->block_size.nx3>1) cng3=cng;
  int f2d=0, f3d=0;
  if(pmb->block_size.nx2 > 1) f2d=1;
  if(pmb->block_size.nx3 > 1) f3d=1;
  int size;
  bd.nbmax=maxneighbor_;

  for(int n=0;n<bd.nbmax;n++) {
    // Clear flags and requests
    bd.flag[n]=BNDRY_WAITING;
    bd.send[n]=NULL;
    bd.recv[n]=NULL;
#ifdef MPI_PARALLEL
    bd.req_send[n]=MPI_REQUEST_NULL;
    bd.req_recv[n]=MPI_REQUEST_NULL;
#endif

    // Allocate buffers
    // calculate the buffer size
    switch(type) {
      case BNDRY_MCOUP: {
        size=((ExchangeValues::ni[n].ox1==0)?pmb->block_size.nx1:NGHOST)
            *((ExchangeValues::ni[n].ox2==0)?pmb->block_size.nx2:NGHOST)
            *((ExchangeValues::ni[n].ox3==0)?pmb->block_size.nx3:NGHOST);
        if(multilevel) {
          int f2c=((ExchangeValues::ni[n].ox1==0)?((pmb->block_size.nx1+1)/2):NGHOST)
                 *((ExchangeValues::ni[n].ox2==0)?((pmb->block_size.nx2+1)/2):NGHOST)
                 *((ExchangeValues::ni[n].ox3==0)?((pmb->block_size.nx3+1)/2):NGHOST);
          int c2f=((ExchangeValues::ni[n].ox1==0)?((pmb->block_size.nx1+1)/2+cng1):cng)
                 *((ExchangeValues::ni[n].ox2==0)?((pmb->block_size.nx2+1)/2+cng2):cng)
                 *((ExchangeValues::ni[n].ox3==0)?((pmb->block_size.nx3+1)/2+cng3):cng);
          size=std::max(size,c2f);
          size=std::max(size,f2c);
        }
        size*=NMCOUP;
      }
      break;
      default: {
        std::stringstream msg;
        msg << "### FATAL ERROR in InitExchangeData" << std::endl
            << "Invalid exchange type is specified." << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      break;
    }
    bd.send[n]=new Real [size];
    bd.recv[n]=new Real [size];
  }
}


//----------------------------------------------------------------------------------------
//! \fn void ExchangeValues::DestroyExchangeData(BoundaryData &bd)
//  \brief Destroy BoundaryData structure
void ExchangeValues::DestroyExchangeData(BoundaryData &bd)
{
  for(int n=0;n<bd.nbmax;n++) {
    delete [] bd.send[n];
    delete [] bd.recv[n];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ExchangeValues::Initialize(void)
//  \brief Initialize MPI requests

void ExchangeValues::Initialize(void)
{
  MeshBlock* pmb=pmy_block_;
  int myox1, myox2, myox3;
  int tag;
  int cng, cng1, cng2, cng3;
  int ssize, rsize;
  cng=cng1=pmb->cnghost;
  cng2=(pmb->block_size.nx2>1)?cng:0;
  cng3=(pmb->block_size.nx3>1)?cng:0;
  long int &lx1=pmb->loc.lx1;
  long int &lx2=pmb->loc.lx2;
  long int &lx3=pmb->loc.lx3;
  int &mylevel=pmb->loc.level;
  myox1=((int)(lx1&1L));
  myox2=((int)(lx2&1L));
  myox3=((int)(lx3&1L));
  int f2d=0, f3d=0;
  if(pmb->block_size.nx2 > 1) f2d=1;
  if(pmb->block_size.nx3 > 1) f3d=1;


  // count the number of the fine meshblocks contacting on each edge
  int eid=0;
  if(pmb->block_size.nx2 > 1) {
    for(int ox2=-1;ox2<=1;ox2+=2) {
      for(int ox1=-1;ox1<=1;ox1+=2) {
        int nis, nie, njs, nje;
        nis=std::max(ox1-1,-1), nie=std::min(ox1+1,1);
        njs=std::max(ox2-1,-1), nje=std::min(ox2+1,1);
        int nf=0, fl=mylevel;
        for(int nj=njs; nj<=nje; nj++) {
          for(int ni=nis; ni<=nie; ni++) {
            if(nblevel[1][nj+1][ni+1] > fl)
              fl++, nf=0;
            if(nblevel[1][nj+1][ni+1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
  }
  if(pmb->block_size.nx3 > 1) {
    for(int ox3=-1;ox3<=1;ox3+=2) {
      for(int ox1=-1;ox1<=1;ox1+=2) {
        int nis, nie, nks, nke;
        nis=std::max(ox1-1,-1), nie=std::min(ox1+1,1);
        nks=std::max(ox3-1,-1), nke=std::min(ox3+1,1);
        int nf=0, fl=mylevel;
        for(int nk=nks; nk<=nke; nk++) {
          for(int ni=nis; ni<=nie; ni++) {
            if(nblevel[nk+1][1][ni+1] > fl)
              fl++, nf=0;
            if(nblevel[nk+1][1][ni+1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
    for(int ox3=-1;ox3<=1;ox3+=2) {
      for(int ox2=-1;ox2<=1;ox2+=2) {
        int njs, nje, nks, nke;
        njs=std::max(ox2-1,-1), nje=std::min(ox2+1,1);
        nks=std::max(ox3-1,-1), nke=std::min(ox3+1,1);
        int nf=0, fl=mylevel;
        for(int nk=nks; nk<=nke; nk++) {
          for(int nj=njs; nj<=nje; nj++) {
            if(nblevel[nk+1][nj+1][1] > fl)
              fl++, nf=0;
            if(nblevel[nk+1][nj+1][1]==fl)
              nf++;
          }
        }
        edge_flag_[eid]=(fl==mylevel);
        nedge_fine_[eid++]=nf;
      }
    }
  }

#ifdef MPI_PARALLEL
  // Initialize non-polar neighbor communications to other ranks
  for(int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if(nb.rank!=Globals::my_rank) {
      if(nb.level==mylevel) { // same
        ssize=rsize=((nb.ox1==0)?pmb->block_size.nx1:NGHOST)
                   *((nb.ox2==0)?pmb->block_size.nx2:NGHOST)
                   *((nb.ox3==0)?pmb->block_size.nx3:NGHOST);
      }
      else if(nb.level<mylevel) { // coarser
        ssize=((nb.ox1==0)?((pmb->block_size.nx1+1)/2):NGHOST)
             *((nb.ox2==0)?((pmb->block_size.nx2+1)/2):NGHOST)
             *((nb.ox3==0)?((pmb->block_size.nx3+1)/2):NGHOST);
        rsize=((nb.ox1==0)?((pmb->block_size.nx1+1)/2+cng1):cng1)
             *((nb.ox2==0)?((pmb->block_size.nx2+1)/2+cng2):cng2)
             *((nb.ox3==0)?((pmb->block_size.nx3+1)/2+cng3):cng3);
      }
      else { // finer
        ssize=((nb.ox1==0)?((pmb->block_size.nx1+1)/2+cng1):cng1)
             *((nb.ox2==0)?((pmb->block_size.nx2+1)/2+cng2):cng2)
             *((nb.ox3==0)?((pmb->block_size.nx3+1)/2+cng3):cng3);
        rsize=((nb.ox1==0)?((pmb->block_size.nx1+1)/2):NGHOST)
             *((nb.ox2==0)?((pmb->block_size.nx2+1)/2):NGHOST)
             *((nb.ox3==0)?((pmb->block_size.nx3+1)/2):NGHOST);
      }
      if (HYBRID){
        ssize*=NMCOUP; rsize*=NMCOUP;
        // specify the offsets in the view point of the target block: flip ox? signs
        tag=CreateBvalsMPITag(nb.lid, TAG_MCOUP, nb.targetid);
        MPI_Send_init(bd_mcoup_.send[nb.bufid],ssize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&(bd_mcoup_.req_send[nb.bufid]));
        tag=CreateBvalsMPITag(pmb->lid, TAG_MCOUP, nb.bufid);
        MPI_Recv_init(bd_mcoup_.recv[nb.bufid],rsize,MPI_ATHENA_REAL,
                      nb.rank,tag,MPI_COMM_WORLD,&(bd_mcoup_.req_recv[nb.bufid]));
      }
    }
  }

#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ExchangeValues::CheckBoundary(void)
//  \brief checks if the boundary conditions are correctly enrolled

void ExchangeValues::CheckBoundary(void)
{
  MeshBlock *pmb=pmy_block_;
  for(int i=0;i<nface_;i++) {
    if(block_bcs[i]==USER_BNDRY) {
      if(ExchangeFunction_[i]==NULL) {
        std::stringstream msg;
        msg << "### FATAL ERROR in ExchangeValues::CheckBoundary" << std::endl
            << "A user-defined boundary is specified but the exchange boundary function "
            << "is not enrolled in direction " << i  << "." << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ExchagneValues::StartReceivingForInit(bool flag)
//  \brief initiate MPI_Irecv for initialization

void ExchangeValues::StartReceivingForInit(bool flag)
{
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_block_;
  for(int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if(nb.rank!=Globals::my_rank) { 
      if (flag) {  // normal case
        if (HYBRID)
          MPI_Start(&(bd_mcoup_.req_recv[nb.bufid]));
      }        
    }
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ExchangeValues::StartReceivingAll(void)
//  \brief initiate MPI_Irecv for all the sweeps

void ExchangeValues::StartReceivingAll(void)
{
  firsttime_=true;
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;
  for(int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if(nb.rank!=Globals::my_rank) { 
      if (HYBRID)
        MPI_Start(&(bd_mcoup_.req_recv[nb.bufid]));
    }
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ExchangeValues::ClearExchangeForInit(void)
//  \brief clean up the exchange flags for initialization

void ExchangeValues::ClearExchangeForInit(bool flag)
{
  MeshBlock *pmb=pmy_block_;

  // Note step==0 corresponds to initial exchange of conserved variables, while step==1
  // corresponds to primitives sent only in the case of GR with refinement
  for(int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if (HYBRID)
      bd_mcoup_.flag[nb.bufid] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank) {
      if (flag) {  // normal case
        if (HYBRID)
          MPI_Wait(&(bd_mcoup_.req_send[nb.bufid]),MPI_STATUS_IGNORE); // Wait for Isend
       }
     }
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void ExchangeValues::ClearExchangeAll(void)
//  \brief clean up the exchange flags after each loop

void ExchangeValues::ClearExchangeAll(void)
{
  MeshBlock *pmb=pmy_block_;

  // Clear non-polar boundary communications
  for(int n=0;n<nneighbor;n++) {
    NeighborBlock& nb = neighbor[n];
    if (HYBRID) {
      bd_mcoup_.flag[nb.bufid] = BNDRY_WAITING;
    }
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank) {
      if (HYBRID)
        MPI_Wait(&(bd_mcoup_.req_send[nb.bufid]),MPI_STATUS_IGNORE); // Wait for Isend
    }
#endif
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void ExchangeValues::ApplyExchangePhysicalBoundaries(AthenaArray<Real> &pdst,
//           AthenaArray<Real> &cdst, const Real time, const Real dt)
//  \brief Apply all the exchange boundary conditions for hybrid coupling array

void ExchangeValues::ApplyExchangePhysicalBoundaries(AthenaArray<Real> &pdst,
                     const Real time, const Real dt)
{
  MeshBlock *pmb=pmy_block_;
  Coordinates *pco=pmb->pcoord;
  int bis=pmb->is-NGHOST, bie=pmb->ie+NGHOST, bjs=pmb->js, bje=pmb->je,
      bks=pmb->ks, bke=pmb->ke;
  if(ExchangeFunction_[INNER_X2]==NULL && pmb->block_size.nx2>1) bjs=pmb->js-NGHOST;
  if(ExchangeFunction_[OUTER_X2]==NULL && pmb->block_size.nx2>1) bje=pmb->je+NGHOST;
  if(ExchangeFunction_[INNER_X3]==NULL && pmb->block_size.nx3>1) bks=pmb->ks-NGHOST;
  if(ExchangeFunction_[OUTER_X3]==NULL && pmb->block_size.nx3>1) bke=pmb->ke+NGHOST;
 
  if(pmb->block_size.nx3>1) { // 3D
    bjs=pmb->js-NGHOST;
    bje=pmb->je+NGHOST;

    // Apply exchange function on inner-x3
    if (ExchangeFunction_[INNER_X3] != NULL) {
      ExchangeFunction_[INNER_X3](pmb, pco, pdst, time, dt,
                                  bis, bie, bjs, bje, pmb->ks, pmb->ke);
    }

    // Apply exchange function on outer-x3
    if (ExchangeFunction_[OUTER_X3] != NULL) {
      ExchangeFunction_[OUTER_X3](pmb, pco, pdst, time, dt,
                                  bis, bie, bjs, bje, pmb->ks, pmb->ke);
    }
  }
  // Apply exchange function on inner-x1
  if (ExchangeFunction_[INNER_X1] != NULL) {
    ExchangeFunction_[INNER_X1](pmb, pco, pdst, time, dt,
                                pmb->is, pmb->ie, bjs, bje, bks, bke);
  }

  // Apply exchange function on outer-x1
  if (ExchangeFunction_[OUTER_X1] != NULL) {
    ExchangeFunction_[OUTER_X1](pmb, pco, pdst, time, dt,
                                pmb->is, pmb->ie, bjs, bje, bks, bke);
  }
  if(pmb->block_size.nx2>1) { // 2D or 3D

    // Apply exchange function on inner-x2
    if (ExchangeFunction_[INNER_X2] != NULL) {
      ExchangeFunction_[INNER_X2](pmb, pco, pdst, time, dt,
                                  bis, bie, pmb->js, pmb->je, bks, bke);
    }

    // Apply exchange function on outer-x2
    if (ExchangeFunction_[OUTER_X2] != NULL) {
      ExchangeFunction_[OUTER_X2](pmb, pco, pdst, time, dt,
                                  bis, bie, pmb->js, pmb->je, bks, bke);
    }
  }

  return;
}
