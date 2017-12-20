//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file time_integrator.cpp
//  \brief derived class for time integrator task list.  Can create task lists for one
//  of many different time integrators (e.g. van Leer, RK2, RK3, etc.)

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include "time.h"
#include <sys/time.h>
// Athena++ classes headers
#include "task_list.hpp"
#include "task_list_hybrid.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../hybrid/hybrid.hpp"
#include "../particle/particle.hpp"
#include "../globals.hpp"
//----------------------------------------------------------------------------------------
//  TimeIntegratorTaskListHybrid constructor

TimeIntegratorTaskListHybrid::TimeIntegratorTaskListHybrid(ParameterInput *pin, Mesh *pm)
  : TaskList(pm)
{

  // Now assemble list of tasks for each step of time integrator
  {using namespace HybridIntegratorTaskNames;

    nsub_steps=2;
 
    AddTimeIntegratorTask(START_ALLRECV,NONE);

    // compute electric field for constrained trasport
    AddTimeIntegratorTask(CALC_EFLD,NONE);
    // integrate field half timestep
    AddTimeIntegratorTask(INT_FLD,CALC_EFLD);
    // MHD boundary values
    AddTimeIntegratorTask(SEND_FLD,INT_FLD);
    AddTimeIntegratorTask(RECV_FLD,START_ALLRECV);
    // push particles
    AddTimeIntegratorTask(PUSH,RECV_FLD);
    // physical boundaries
    // AddTimeInregratorTask(PHYS_BNDR,PUSH);
    // predictor step: outflow - add particles to ghost zones
    // corrector step: outflow - remove particles from ghost zones
    //                 reflecting - reflect particles from ghost zones
    //                 periodic - shift positions
    //                 shearing - shift positions and velocities
    // add GetPos to this task
    // make GetPos function a pointer which depends on bvals

    // exchange particles
    AddTimeIntegratorTask(PACK_PRTL,PUSH);
    AddTimeIntegratorTask(SEND_CNTS,PACK_PRTL);
    AddTimeIntegratorTask(RECV_CNTS,NONE);
    AddTimeIntegratorTask(WAIT_CNTS,RECV_CNTS);
    AddTimeIntegratorTask(RECV_PRTL,WAIT_CNTS);
    AddTimeIntegratorTask(SEND_PRTL,WAIT_CNTS); 
    AddTimeIntegratorTask(DEPOSIT,PACK_PRTL);
    AddTimeIntegratorTask(WAIT_PRTL,RECV_PRTL);
    // deposit
    AddTimeIntegratorTask(DEPOSIT_BUFF,(DEPOSIT|WAIT_PRTL)); 
    AddTimeIntegratorTask(UNPACK_PRTL,DEPOSIT_BUFF);

    // exchange moments of distribution function
    AddTimeIntegratorTask(SEND_EXCH,DEPOSIT_BUFF);
    AddTimeIntegratorTask(RECV_EXCH,START_ALLRECV);

    // everything else
    AddTimeIntegratorTask(USERWORK,RECV_EXCH);
    AddTimeIntegratorTask(NEW_DT,USERWORK);
    AddTimeIntegratorTask(CLEAR_ALLBND,NEW_DT);

  } // end of using namespace block
}

//----------------------------------------------------------------------------------------//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  

void TimeIntegratorTaskListHybrid::AddTimeIntegratorTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace HybridIntegratorTaskNames;
  switch((id)) {
    case (START_ALLRECV):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::StartAllReceive);
      break;
    case (CLEAR_ALLBND):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::ClearAllBoundary);
      break;

    case (CALC_EFLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::CalculateEMF);
      break;

    case (INT_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::FieldIntegrate);
      break;

    case (SEND_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::FieldSend);
      break;
    case (RECV_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::FieldReceive);
      break;

    case (SEND_EXCH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::ExchangeSend);
      break;
    case (RECV_EXCH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::ExchangeReceive);
      break;

    case (PACK_PRTL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::PackParticle);
      break;

    case (UNPACK_PRTL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::UnPackParticle);
      break;

    case (SEND_CNTS):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::CountsSend);
      break;
    case (RECV_CNTS):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::CountsReceive);
      break;

    case (SEND_PRTL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::ParticleSend);
      break;
    case (RECV_PRTL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::ParticleReceive);
      break;

    case (WAIT_CNTS):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::CountsWait);
      break;
    case (WAIT_PRTL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::ParticleWait);
      break;


    case (PUSH):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::Push);
      break;

    case (DEPOSIT):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::Deposit);
      break;

    case (DEPOSIT_BUFF):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::DepositBuffer);
      break;

     case (USERWORK):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::UserWork);
      break;
    case (NEW_DT):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskListHybrid::NewBlockTimeStep);
      break;

    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTimeIntegratorTaskHybrid" << std::endl
          << "Invalid Task "<< id << " is specified" << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief

//----------------------------------------------------------------------------------------
// Functions to start/end MPI communication

enum TaskStatus TimeIntegratorTaskListHybrid::StartAllReceive(MeshBlock *pmb, int step)
{
  //pmb->pbval->StartReceivingAll();
  //pmb->peval->StartReceivingAll();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::ClearAllBoundary(MeshBlock *pmb, int step)
{
  //pmb->pbval->ClearBoundaryAll();
  //pmb->peval->ClearExchangeAll();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to calculates fluxes

enum TaskStatus TimeIntegratorTaskListHybrid::CalculateEMF(MeshBlock *pmb, int step)
{

  return TASK_NEXT;
}


//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

enum TaskStatus TimeIntegratorTaskListHybrid::FieldIntegrate(MeshBlock *pmb, int step)
{

  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to communicate conserved variables between MeshBlocks


enum TaskStatus TimeIntegratorTaskListHybrid::FieldSend(MeshBlock *pmb, int step)
{
 /* if(step == 1) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b1);
  } else if(step == 2) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }*/
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskListHybrid::FieldReceive(MeshBlock *pmb, int step)
{
  /*bool ret;
  if(step == 1) {
    ret=pmb->pbval->ReceiveFieldBoundaryBuffers(pmb->pfield->b1);
  } else if(step == 2) {
    ret=pmb->pbval->ReceiveFieldBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }
  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }*/
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions for everything else

enum TaskStatus TimeIntegratorTaskListHybrid::UserWork(MeshBlock *pmb, int step)
{
  //if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  //pmb->UserWorkInLoop();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::NewBlockTimeStep(MeshBlock *pmb, int step)
{
  //if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::Push(MeshBlock *pmb, int step)
{
  if (step==1) {
    struct timeval tv_init,tv_curr;
    gettimeofday(&tv_init,NULL);
    
    pmb->particle->MoveDeposit();

    pmb->phybrid->TransposeMCOUP();

    gettimeofday(&tv_curr,NULL);
    Real step_time = (double)(tv_curr.tv_sec - tv_init.tv_sec)
                   + 1.0e-6*(double)(tv_curr.tv_usec - tv_init.tv_usec);
    pmb->particle->timer_movedeposit = step_time;
  } else if (step==2) {
    struct timeval tv_init,tv_curr;
    gettimeofday(&tv_init,NULL);

    pmb->particle->Move();

    gettimeofday(&tv_curr,NULL);
    Real step_time = (double)(tv_curr.tv_sec - tv_init.tv_sec)
                   + 1.0e-6*(double)(tv_curr.tv_usec - tv_init.tv_usec);
    pmb->particle->timer_move = step_time; 
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::Deposit(MeshBlock *pmb, int step)
{
  if (step==2) { 
    struct timeval tv_init,tv_curr;
    gettimeofday(&tv_init,NULL);
    
    pmb->particle->Deposit();
    
    gettimeofday(&tv_curr,NULL);
    Real step_time = (double)(tv_curr.tv_sec - tv_init.tv_sec)
                   + 1.0e-6*(double)(tv_curr.tv_usec - tv_init.tv_usec);
    pmb->particle->timer_deposit = step_time;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::DepositBuffer(MeshBlock *pmb, int step)
{
  if (step==2){
    struct timeval tv_init,tv_curr;
    gettimeofday(&tv_init,NULL);

    pmb->particle->DepositBuffer();

    pmb->phybrid->TransposeMCOUP();

    gettimeofday(&tv_curr,NULL);
    Real step_time = (double)(tv_curr.tv_sec - tv_init.tv_sec)
                   + 1.0e-6*(double)(tv_curr.tv_usec - tv_init.tv_usec);
    pmb->particle->timer_deposit += step_time;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::CountsSend(MeshBlock *pmb, int step)
{
  if (step==2) {
    pmb->particle->SendCounts();
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::CountsReceive(MeshBlock *pmb, int step)
{
  if (step==2) {
    pmb->particle->ReceiveCounts();
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::CountsWait(MeshBlock *pmb, int step)
{
  if (step==2) {
    struct timeval tv_init,tv_curr;
    gettimeofday(&tv_init,NULL);
 
    pmb->particle->WaitCounts();
 
    gettimeofday(&tv_curr,NULL);
    Real step_time = (double)(tv_curr.tv_sec - tv_init.tv_sec)
                   + 1.0e-6*(double)(tv_curr.tv_usec - tv_init.tv_usec);
    pmb->particle->timer_wait = step_time;
  }

  return TASK_SUCCESS;
}
enum TaskStatus TimeIntegratorTaskListHybrid::ParticleSend(MeshBlock *pmb, int step)
{
  if (step==2) {
    pmb->particle->SendParticle();
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::ParticleReceive(MeshBlock *pmb, int step)
{
  if (step==2) {
    pmb->particle->ReceiveParticle();
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::ParticleWait(MeshBlock *pmb, int step)
{
  if (step==2) {
    struct timeval tv_init,tv_curr;
    gettimeofday(&tv_init,NULL);

    pmb->particle->WaitParticle();

    gettimeofday(&tv_curr,NULL);
    Real step_time = (double)(tv_curr.tv_sec - tv_init.tv_sec)
                   + 1.0e-6*(double)(tv_curr.tv_usec - tv_init.tv_usec);
    pmb->particle->timer_wait += step_time;
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::ExchangeSend(MeshBlock *pmb, int step)
{
  pmb->peval->SendExchangeBuffers(pmb->phybrid->mcoup_,HYBRID_MCOUP);

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::ExchangeReceive(MeshBlock *pmb, int step)
{
  pmb->peval->ReceiveExchangeBuffersWithWait(pmb->phybrid->mcoup_,HYBRID_MCOUP);

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::PackParticle(MeshBlock *pmb, int step)
{
  if (step==2) {

    struct timeval tv_init,tv_curr;
    gettimeofday(&tv_init,NULL);

    pmb->particle->GetPos();
    pmb->particle->PackParticle();

    gettimeofday(&tv_curr,NULL);
    Real step_time = (double)(tv_curr.tv_sec - tv_init.tv_sec) 
              + 1.0e-6*(double)(tv_curr.tv_usec - tv_init.tv_usec);
    pmb->particle->timer_exchange = step_time;
  }

  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskListHybrid::UnPackParticle(MeshBlock *pmb, int step)
{
  if (step==2) {
    pmb->particle->UnPackParticle();
  }

  return TASK_SUCCESS;
}
