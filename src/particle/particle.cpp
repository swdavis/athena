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
#include <sstream>    //stringstream
#include <stdlib.h>

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "particle.hpp"
#include "../bvals/bvals.hpp"
#include "../globals.hpp"
#include "../utils/utils.hpp"

class BoundaryBase;

//! \class Particle
// \brief data/functions associated with list of particles

Particle::Particle(MeshBlock *pmb, ParameterInput *pin)
{
  std::stringstream msg;

  pmy_block = pmb;
  max_prtl = pin->GetInteger("particle","max_prtl");
  buffer_size = pin->GetInteger("particle","buffer_size");
//  x1.NewAthenaArray(max_prtl);
//  x2.NewAthenaArray(max_prtl);
//  x3.NewAthenaArray(max_prtl);
//  v1.NewAthenaArray(max_prtl);
//  v2.NewAthenaArray(max_prtl);
//  v3.NewAthenaArray(max_prtl);
//  pos.NewAthenaArray(max_prtl);

  posix_memalign((void **)&x1, CACHELINE_BYTES, max_prtl*sizeof(Real));
  posix_memalign((void **)&x2, CACHELINE_BYTES, max_prtl*sizeof(Real));
  posix_memalign((void **)&x3, CACHELINE_BYTES, max_prtl*sizeof(Real));
  posix_memalign((void **)&v1, CACHELINE_BYTES, max_prtl*sizeof(Real));
  posix_memalign((void **)&v2, CACHELINE_BYTES, max_prtl*sizeof(Real));
  posix_memalign((void **)&v3, CACHELINE_BYTES, max_prtl*sizeof(Real));
  posix_memalign((void**)&pos, CACHELINE_BYTES, max_prtl*sizeof(short int));

  for (int i=0;i<56;i++){
    send_cnt[i]=0; recv_cnt[i]=0;
  }
  send_buf = new Real[6*buffer_size];
  recv_buf = new Real[6*buffer_size];
  int ncells1 = pmb->block_size.nx1;
  int ncells2 = pmb->block_size.nx2;
  int ncells3 = pmb->block_size.nx3;
  nparticle = ncells1*ncells2*ncells3 * pin->GetInteger("particle","ppc0");

  srand(1231);
  long int idum=-1;
  for (int p=0; p<nparticle; p++) {
    Real r1 = 1.0*rand()/RAND_MAX;
    Real r2 = 1.0*rand()/RAND_MAX;
    Real r3 = 1.0*rand()/RAND_MAX;
//    x1(p)= (pmb->block_size.x1max-pmb->block_size.x1min)*r1 + pmb->block_size.x1min;
//    x2(p)= (pmb->block_size.x2max-pmb->block_size.x2min)*r2 + pmb->block_size.x2min;
//    x3(p)= (pmb->block_size.x3max-pmb->block_size.x3min)*r3 + pmb->block_size.x3min;
//    v1(p)=0.01;
//    v2(p)=0.1;
//    v3(p)=-0.05;
//    pos(p) = BoundaryBase::CreateBufferID(0,0,0,0,0);

    x1[p]= (pmb->block_size.x1max-pmb->block_size.x1min)*r1 + pmb->block_size.x1min;
    x2[p]= (pmb->block_size.x2max-pmb->block_size.x2min)*r2 + pmb->block_size.x2min;
    x3[p]= (pmb->block_size.x3max-pmb->block_size.x3min)*r3 + pmb->block_size.x3min;
    v1[p]= ran_gaussian(&idum);
    v2[p]= ran_gaussian(&idum);
    v3[p]= ran_gaussian(&idum);
    pos[p] = BoundaryBase::CreateBufferID(0,0,0,0,0);
  }



  if (nparticle > max_prtl) {
    msg << "### Overflow while initializing particle array" << std::endl
        << "nparticle = " << nparticle << std::endl
        << "max_prtl  = " << max_prtl << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  for (int n=0; n<pmb->pbval->nneighbor; n++){
    flag_par[n] = BNDRY_WAITING;
    flag_cnt[n] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    req_par_send[n] = MPI_REQUEST_NULL;
    req_par_recv[n] = MPI_REQUEST_NULL;
    req_cnt_send[n] = MPI_REQUEST_NULL;
    req_cnt_recv[n] = MPI_REQUEST_NULL;
#endif
  }
}

Particle::~Particle()
{
//  x1.DeleteAthenaArray();
//  x2.DeleteAthenaArray();
//  x3.DeleteAthenaArray();
//  v1.DeleteAthenaArray();
//  v2.DeleteAthenaArray();
//  v3.DeleteAthenaArray();
//  pos.DeleteAthenaArray();

  free(x1);
  free(x2);
  free(x3);
  free(v1);
  free(v2);
  free(v3);
  free(pos);

  delete send_buf;
  delete recv_buf;
}


double Particle::ran_gaussian (long int *idum)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  if (*idum < 0) iset = 0;
  if (iset == 0) {
    do {
      v1 = 2.0 * ran2(idum) - 1.0;
      v2 = 2.0 * ran2(idum) - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >=1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1 * fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  } 
}
