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

// Athena++ classes headers
#include "task_list.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../gravity/gravity.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"

//----------------------------------------------------------------------------------------
//  TimeIntegratorTaskList constructor

TimeIntegratorTaskList::TimeIntegratorTaskListHybrid(ParameterInput *pin, Mesh *pm)
  : TaskList(pm)
{

  {using namespace HybridIntegratorTaskNames;

    // start receiving for exchange bvals and mhd bvals
    AddTimeIntegratorTask(START_ALLRECV,NONE);

    // deposit particles for predictor step
    AddTimeIntegratorTask(DEPOSIT,START_ALLRECV);
    
    // exchange deposited quantities
    AddTimeIntegratorTask(SEND_EXCH,DEPOSIT);
    AddTimeIntegratorTask(RECV_EXCH,DEPOSIT);

    // calculate electric field
    AddTimeIntegratorTask(CALC_EFLD,RECV_EXCH);

    // advance magnetic field by half timestep
    AddTimeIntegratorTask(INT_FLD,CALC_EFLD);

    // move particles for predictor step and deposit
    AddTimeIntegratorTask(MOVEDEPOSIT,INT_FLD);

    // exchange deposited quantities
    AddTimeIntegratorTask(SEND_EXCH,MOVEDEPOSIT);
    AddTimeIntegratorTask(RECV_EXCH,MOVEDEPOSIT);

    // calculate electric field
    AddTimeIntegratorTask(CALC_EFLD,RECV_EXCH);

    // advance magnetic field by half timestep
    AddTimeIntegratorTask(INT_FLD,CALC_EFLD);

    // move particles for corrector step
    AddTimeIntegratorTask(MOVE,INT_FLD);

    // exchange particles
    AddTimeIntegratorTask(PACK_PRTL,MOVE)
    AddTimeIntegratorTask(SEND_CNTS,PACK_PRTL);
    AddTimeIntegratorTask(RECV_CNTS,PACK_PRTL);
    AddTimeIntegratorTask(SEND_PRTL,RECV_CNTS);
    AddTimeIntegratorTask(RECV_PRTL,RECV_CNTS);
    AddTimeIntegratorTask(UNPACK_PRTL,RECV_PRTL);

  } // end of using namespace block
}

//----------------------------------------------------------------------------------------//! \fn
//  \brief Sets id and dependency for "ntask" member of task_list_ array, then iterates
//  value of ntask.  

void TimeIntegratorTaskList::AddTimeIntegratorTask(uint64_t id, uint64_t dep)
{
  task_list_[ntasks].task_id=id;
  task_list_[ntasks].dependency=dep;

  using namespace HydroIntegratorTaskNames;
  switch((id)) {
    case (START_ALLRECV):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::StartAllReceive);
      break;
    case (CLEAR_ALLBND):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ClearAllBoundary);
      break;

    case (CALC_HYDFLX):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateFluxes);
      break;
    case (CALC_FLDFLX):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateEMF);
      break;

    case (SEND_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FluxCorrectSend);
      break;
    case (SEND_FLDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFCorrectSend);
      break;

    case (RECV_HYDFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FluxCorrectReceive);
      break;
    case (RECV_FLDFLX):
      task_list_[ntasks].TaskFunc= 
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::EMFCorrectReceive);
      break;

    case (INT_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroIntegrate);
      break;
    case (INT_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldIntegrate);
      break;

    case (SRCTERM_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroSourceTerms);
      break;

    case (SEND_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroSend);
      break;
    case (SEND_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldSend);
      break;

    case (RECV_HYD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::HydroReceive);
      break;
    case (RECV_FLD):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::FieldReceive);
      break;

    case (PROLONG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Prolongation);
      break;
    case (CON2PRIM):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Primitives);
      break;
    case (PHY_BVAL):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::PhysicalBoundary);
      break;
    case (USERWORK):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::UserWork);
      break;
    case (NEW_DT):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::NewBlockTimeStep);
      break;
    case (AMR_FLAG):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CheckRefinement);
      break;

    case (SOLV_GRAV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::GravSolve);
      break;
    case (SEND_GRAV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::GravSend);
      break;
    case (RECV_GRAV):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::GravReceive);
      break;
    case (CORR_GFLX):
      task_list_[ntasks].TaskFunc=
        static_cast<enum TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::GravFluxCorrection);
      break;


    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in AddTimeIntegratorTask" << std::endl
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

enum TaskStatus TimeIntegratorTaskList::StartAllReceive(MeshBlock *pmb, int step)
{
  pmb->pbval->StartReceivingAll();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::ClearAllBoundary(MeshBlock *pmb, int step)
{
  pmb->pbval->ClearBoundaryAll();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to calculates fluxes

enum TaskStatus TimeIntegratorTaskList::CalculateFluxes(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;

  if((step == 1) && (integrator == "vl2")) {
    phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, 1);
    return TASK_NEXT;
  }

  if((step == 1) && (integrator == "rk2")) {
    phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, 2);
    return TASK_NEXT;
  } 

  if(step == 2) {
    phydro->CalculateFluxes(phydro->w1, pfield->b1, pfield->bcc1, 2);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::CalculateEMF(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pfield->ComputeCornerE(pmb->phydro->w,  pmb->pfield->bcc);
    return TASK_NEXT;
  }

  if(step == 2) {
    pmb->pfield->ComputeCornerE(pmb->phydro->w1, pmb->pfield->bcc1);
    return TASK_NEXT;
  } 

  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to communicate fluxes between MeshBlocks for flux correction step with AMR

enum TaskStatus TimeIntegratorTaskList::FluxCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendFluxCorrection(FLUX_HYDRO);
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectSend(MeshBlock *pmb, int step)
{
  pmb->pbval->SendEMFCorrection();
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive fluxes between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::FluxCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveFluxCorrection(FLUX_HYDRO) == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::EMFCorrectReceive(MeshBlock *pmb, int step)
{
  if(pmb->pbval->ReceiveEMFCorrection() == true) {
    return TASK_NEXT;
  } else {
    return TASK_FAIL;
  }
}

//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

enum TaskStatus TimeIntegratorTaskList::HydroIntegrate(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;

  if(step == 1) {
    ph->AddFluxDivergenceToAverage(ph->u,ph->u,ph->w,pf->bcc,step_wghts[0],ph->u1);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "vl2")) {
    ph->AddFluxDivergenceToAverage(ph->u,ph->u,ph->w1,pf->bcc1,step_wghts[1],ph->u);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "rk2")) {
   ph->AddFluxDivergenceToAverage(ph->u,ph->u1,ph->w1,pf->bcc1,step_wghts[1],ph->u);
   return TASK_NEXT;
  }

  return TASK_FAIL;
}

enum TaskStatus TimeIntegratorTaskList::FieldIntegrate(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pfield->CT(pmb->pfield->b, pmb->pfield->b, step_wghts[0], pmb->pfield->b1);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "vl2")) {
    pmb->pfield->CT(pmb->pfield->b, pmb->pfield->b, step_wghts[1], pmb->pfield->b);
    return TASK_NEXT;
  }

  if((step == 2) && (integrator == "rk2")) {
    pmb->pfield->CT(pmb->pfield->b, pmb->pfield->b1, step_wghts[1], pmb->pfield->b);
    return TASK_NEXT;
  }

  return TASK_FAIL;
}

//----------------------------------------------------------------------------------------
// Functions to add source terms

enum TaskStatus TimeIntegratorTaskList::HydroSourceTerms(MeshBlock *pmb, int step)
{
  Hydro *ph=pmb->phydro;
  Field *pf=pmb->pfield;

  // return if there are no source terms to be added
  if (ph->psrc->hydro_sourceterms_defined == false) return TASK_NEXT;

  Real dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
  Real time;
  // *** this must be changed for the RK3 integrator
  if(step == 1) {
    time=pmb->pmy_mesh->time;
    ph->psrc->AddHydroSourceTerms(time,dt,ph->flux,ph->w,pf->bcc,ph->u1);
  } else if(step == 2) {
    if      (integrator == "vl2") time=pmb->pmy_mesh->time + 0.5*pmb->pmy_mesh->dt;
    else if (integrator == "rk2") time=pmb->pmy_mesh->time +     pmb->pmy_mesh->dt;
    ph->psrc->AddHydroSourceTerms(time,dt,ph->flux,ph->w1,pf->bcc1,ph->u);
  } else {
    return TASK_FAIL;
  }

  return TASK_NEXT;
}

//----------------------------------------------------------------------------------------
// Functions to communicate conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroSend(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pbval->SendCellCenteredBoundaryBuffers(pmb->phydro->u1, HYDRO_CONS);
  } else if(step == 2) {
    pmb->pbval->SendCellCenteredBoundaryBuffers(pmb->phydro->u, HYDRO_CONS);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::FieldSend(MeshBlock *pmb, int step)
{
  if(step == 1) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b1);
  } else if(step == 2) {
    pmb->pbval->SendFieldBoundaryBuffers(pmb->pfield->b);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

enum TaskStatus TimeIntegratorTaskList::HydroReceive(MeshBlock *pmb, int step)
{
  bool ret;
  if(step == 1) {
    ret=pmb->pbval->ReceiveCellCenteredBoundaryBuffers(pmb->phydro->u1, HYDRO_CONS);
  } else if(step == 2) {
    ret=pmb->pbval->ReceiveCellCenteredBoundaryBuffers(pmb->phydro->u, HYDRO_CONS);
  } else {
    return TASK_FAIL;
  }
  if(ret==true) {
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::FieldReceive(MeshBlock *pmb, int step)
{
  bool ret;
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
  }
}

//----------------------------------------------------------------------------------------
// Functions for everything else

enum TaskStatus TimeIntegratorTaskList::Prolongation(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;

  if(step == 1) {
    dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
    pbval->ProlongateBoundaries(phydro->w1, phydro->u1, pfield->b1, pfield->bcc1,
                                pmb->pmy_mesh->time+dt, dt);
  } else if(step == 2) {
    dt=pmb->pmy_mesh->dt;
    pbval->ProlongateBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                pmb->pmy_mesh->time+dt, dt);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::Primitives(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  if(pbval->nblevel[1][1][0]!=-1) is-=NGHOST;
  if(pbval->nblevel[1][1][2]!=-1) ie+=NGHOST;
  if(pbval->nblevel[1][0][1]!=-1) js-=NGHOST;
  if(pbval->nblevel[1][2][1]!=-1) je+=NGHOST;
  if(pbval->nblevel[0][1][1]!=-1) ks-=NGHOST;
  if(pbval->nblevel[2][1][1]!=-1) ke+=NGHOST;

  if(step == 1) {
    pmb->peos->ConservedToPrimitive(phydro->u1, phydro->w, pfield->b1,
                                    phydro->w1, pfield->bcc1, pmb->pcoord,
                                    is, ie, js, je, ks, ke);
  } else if(step == 2) {
    pmb->peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b,
                                    phydro->w, pfield->bcc, pmb->pcoord,
                                    is, ie, js, je, ks, ke);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::PhysicalBoundary(MeshBlock *pmb, int step)
{
  Hydro *phydro=pmb->phydro;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  Real dt;
  if(step == 1) {
    dt = (step_wghts[(step-1)].c)*(pmb->pmy_mesh->dt);
    pbval->ApplyPhysicalBoundaries(phydro->w1, phydro->u1, pfield->b1, pfield->bcc1,
                                   pmb->pmy_mesh->time+dt, dt);
  } else if(step == 2) {
    dt=pmb->pmy_mesh->dt;
    pbval->ApplyPhysicalBoundaries(phydro->w,  phydro->u,  pfield->b,  pfield->bcc,
                                   pmb->pmy_mesh->time+dt, dt);
  } else {
    return TASK_FAIL;
  }
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::UserWork(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->UserWorkInLoop();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::NewBlockTimeStep(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->phydro->NewBlockTimeStep();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::CheckRefinement(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->pmr->CheckRefinementCondition();
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::GravSolve(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_NEXT; // only do on last sub-step

  pmb->pgrav->Solver(pmb->phydro->u);
  return TASK_NEXT;
}

enum TaskStatus TimeIntegratorTaskList::GravSend(MeshBlock *pmb, int step)
{
//  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  pmb->pbval->SendGravityBoundaryBuffers(pmb->pgrav->phi);
  return TASK_SUCCESS;
}

enum TaskStatus TimeIntegratorTaskList::GravReceive(MeshBlock *pmb, int step)
{
//  if (step != nsub_steps) return TASK_SUCCESS; // only do on last sub-step

  if(pmb->pbval->ReceiveGravityBoundaryBuffers(pmb->pgrav->phi) == true){
    return TASK_SUCCESS;
  } else {
    return TASK_FAIL;
  }
}

enum TaskStatus TimeIntegratorTaskList::GravFluxCorrection(MeshBlock *pmb, int step)
{
  if (step != nsub_steps) return TASK_NEXT; // only do on last sub-step

  pmb->phydro->CorrectGravityFlux();
  return TASK_NEXT;
}
