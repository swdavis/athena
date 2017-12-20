#ifndef TASK_LIST_HYBRID_HPP
#define TASK_LIST_HYBRID_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//!   \file task_list_hybrid.hpp
//    \brief provides functionality to control dynamic execution using tasks in hybrid-PIC

#include <stdint.h>

// Athena++ headers
#include "../athena.hpp"
#include "task_list.hpp"

// forward declarations
class Mesh;
class MeshBlock;
class TaskList;

//----------------------------------------------------------------------------------------
//! \class TimeIntegratorTaskListHybrid
//  \brief data and function definitions for TimeIntegratorTaskListHybrid derived class

class TimeIntegratorTaskListHybrid : public TaskList {
public:
  TimeIntegratorTaskListHybrid(ParameterInput *pin, Mesh *pm);
  ~TimeIntegratorTaskListHybrid() {};

  void AddTimeIntegratorTask(uint64_t id, uint64_t dep);

  // functions
  enum TaskStatus StartAllReceive(MeshBlock *pmb, int step);
  enum TaskStatus ClearAllBoundary(MeshBlock *pmb, int step);

  enum TaskStatus CalculateEMF(MeshBlock *pmb, int step);

  enum TaskStatus FieldIntegrate(MeshBlock *pmb, int step);

  enum TaskStatus ExchangeSend(MeshBlock *pmb, int step);
  enum TaskStatus FieldSend(MeshBlock *pmb, int step);

  enum TaskStatus ExchangeReceive(MeshBlock *pmb, int step);
  enum TaskStatus FieldReceive(MeshBlock *pmb, int step);

  enum TaskStatus Push(MeshBlock *pmb, int step);
  enum TaskStatus Deposit(MeshBlock *pmb, int step);
  enum TaskStatus DepositBuffer(MeshBlock *pmb, int step);

  enum TaskStatus CountsSend(MeshBlock *pmb, int step);
  enum TaskStatus CountsReceive(MeshBlock *pmb, int step);
  enum TaskStatus CountsWait(MeshBlock *pmb, int step);

  enum TaskStatus PackParticle(MeshBlock *pmb, int step);
  enum TaskStatus UnPackParticle(MeshBlock *pmb, int step);

  enum TaskStatus ParticleSend(MeshBlock *pmb, int step);
  enum TaskStatus ParticleReceive(MeshBlock *pmb, int step);
  enum TaskStatus ParticleWait(MeshBlock *pmb, int step);

  enum TaskStatus UserWork(MeshBlock *pmb, int step);
  enum TaskStatus NewBlockTimeStep(MeshBlock *pmb, int step);

};


//----------------------------------------------------------------------------------------
// 64-bit integers with "1" in different bit positions used to ID each hybrid task.

namespace HybridIntegratorTaskNames {
  const uint64_t NONE=0;
  const uint64_t START_ALLRECV=1LL<<0;
  const uint64_t CLEAR_ALLBND=1LL<<1;

  // exchange deposited quantities
  const uint64_t SEND_EXCH=1LL<<2;
  const uint64_t RECV_EXCH=1LL<<3;

  // calculate electric field for CT
  const uint64_t CALC_EFLD=1LL<<4;
  // integrate magnetic field half timestep
  const uint64_t INT_FLD=1LL<<5;

  // MHD boundary values
  const uint64_t SEND_FLD=1LL<<6;
  const uint64_t RECV_FLD=1LL<<7;

  // Move() and MoveDeposit()
  const uint64_t PUSH=1LL<<8;

  // particle exchange
  const uint64_t PACK_PRTL=1LL<<9;
  const uint64_t SEND_CNTS=1LL<<10;
  const uint64_t RECV_CNTS=1LL<<11;
  const uint64_t SEND_PRTL=1LL<<12;
  const uint64_t RECV_PRTL=1LL<<13;

  // deposit
  const uint64_t DEPOSIT=1LL<<14;
  const uint64_t DEPOSIT_BUFF=1LL<<15;
  const uint64_t UNPACK_PRTL=1LL<<16;

  const uint64_t NEW_DT  =1LL<<17;
  const uint64_t USERWORK=1LL<<18;

  const uint64_t WAIT_CNTS=1LL<<19;
  const uint64_t WAIT_PRTL=1LL<<20;
};

#endif // TASK_LIST_HYBRID_HPP
