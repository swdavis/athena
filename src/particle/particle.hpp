#ifndef PARTICLE_HPP
#define PARTICLE_HPP
//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//=======================================================================================
//! \file particle.hpp
//  \brief defines particle classes
//  These classes contain data and functions related to particles
//=======================================================================================

// TODO:
// 1. Initialization of particles
// 2. Add boundary structure to particle class
// 3. Physical boundaries
// 4. DepositBuffer
// 5. Add pointers to move and deposit functions (and set initialization of those functions)
// 6. Add vtk output of mcoup


// Athena headers
#include "../bvals/bvals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#define CACHELINE_BYTES 64
#define SIMD_WIDTH 16

class MeshBlock;
class ParameterInput;

//! \struct ParticleBoundaryData
//  \brief structure storing information for particle exchange
typedef struct ParticleBoundaryData {
  int nbmax;
  enum BoundaryStatus flag[56];

#ifdef MPI_PARALLEL
  MPI_Request req_recv[56];
#endif
} ParticleBoundaryData;

//! \class Particle
// \brief data/functions associated with list of particles

class Particle {
public:
  Particle(MeshBlock *pmb, ParameterInput *pin);   // create particle array
  ~Particle();                                     // destroy particle array

  MeshBlock* pmy_block;                            // MeshBlock pointer

// Properties of particle array
//---------------------------------------------------------------------------------------
  //AthenaArray<Real> x1,x2,x3,v1,v2,v3;             // particle positions and velocities
  long int nparticle, max_prtl;                    // number of particles, maximum number
                                                   // of particles per meshblock
  //AthenaArray<Real> prop;                          // array of particle properties

  Real *x1, *x2, *x3, *v1, *v2, *v3;
  Real *prop;
  static double ran_gaussian(long int *idum);
// Buffers and functions for particle exchange
//---------------------------------------------------------------------------------------
  long int buffer_size;                            // buffer size for particle exchange
  Real *send_buf, *recv_buf;                       // exchange buffers
  int send_cnt[56], recv_cnt[56];                  // numbers of particles to send to
                                                   // neighboring meshblocks, numbers
                                                   // of particles to receive
  enum BoundaryStatus flag_par[56];
  enum BoundaryStatus flag_cnt[56];
#ifdef MPI_PARALLEL
  MPI_Request req_par_send[56];
  MPI_Request req_par_recv[56];
  MPI_Request req_cnt_send[56];
  MPI_Request req_cnt_recv[56];
#endif

// timers
  Real timer_move, timer_deposit, timer_movedeposit,
       timer_exchange, timer_wait;

  void SendParticle();
  void ReceiveParticle();
  void SendCounts();
  void ReceiveCounts();
  void WaitCounts();
  void WaitParticle();
  void CheckCnt();
  void Check();
  void Move();
  void MoveDeposit();
  void Deposit();
  void DepositBuffer();

//  AthenaArray<short int> pos;                      // particle positions on the grid
  short int *pos;

  void GetPos();                                   // function to set pos array
  void PackParticle();                             // fill send_buf
  void ExchangeCounts();                           // send/receive particle counts
  void UnPackParticle();
  void Test();
  void ExchangeParticle();                         // send/receive particles

  void Initialize();                               // initialize particle array
};

#endif
