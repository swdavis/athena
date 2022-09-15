//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//! \brief implementation of functions in class MonteCarlo, MCRandom

// C++ headers
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "montecarlo.hpp"
#include "../globals.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../utils/buffer_utils.hpp"
#include "../scalars/scalars.hpp"

// GSL library
#if RAN3 == 0
#include <gsl/gsl_randist.h>
#endif

//----------------------------------------------------------------------------------------
//! MonteCarlo constructor, builds monte carlo using parameters in input file

MonteCarlo::MonteCarlo(ParameterInput *pin, Mesh *pmesh) {

  MonteCarloBlock *pfirst;
  std::stringstream msg;

  pmy_mesh = pmesh;

  UserWorkInMove=nullptr;
  InitEmission=nullptr;
  GetTemperature=nullptr;

  // Set flags that control emission, absorption and scattering
  emission_meth = GetEmissionFlag(pin->GetOrAddString("montecarlo","emission","none"));
  if (emission_meth == EMISNONE) {
    InitEmission = nullptr; // left unset
    emission_array_flag = false; // do not allocate memory for array
  } else if (emission_meth ==  EMISUSER) {
    InitEmission = nullptr; // must be set in InitUserMonteCarloData
    emission_array_flag = true; // allocate memory for array
  } else if (emission_meth ==  EMISFF) {
    InitEmission = InitializeEmissionFreeFree;
    emission_array_flag = true; // allocate memory for array
  }

  // read bc flags for each of the 6 physical boundaries.
  mc_bcs[BoundaryFace::inner_x1] = GetMCBoundaryFlag(pin->GetString("mesh","ix1_mc_bc"));
  mc_bcs[BoundaryFace::outer_x1] = GetMCBoundaryFlag(pin->GetString("mesh","ox1_mc_bc"));
  mc_bcs[BoundaryFace::inner_x2] = GetMCBoundaryFlag(pin->GetString("mesh","ix2_mc_bc"));
  mc_bcs[BoundaryFace::outer_x2] = GetMCBoundaryFlag(pin->GetString("mesh","ox2_mc_bc"));
  mc_bcs[BoundaryFace::inner_x3] = GetMCBoundaryFlag(pin->GetString("mesh","ix3_mc_bc"));
  mc_bcs[BoundaryFace::outer_x3] = GetMCBoundaryFlag(pin->GetString("mesh","ox3_mc_bc"));

  // intitialize boundary functions
  for (int dir=0; dir<6; dir++)
    BoundaryFunction_[dir]=nullptr;

  boosts = pin->GetOrAddBoolean("montecarlo","boosts",false);
  coupled = pin->GetOrAddBoolean("montecarlo","coupled",false);
  polarized = pin->GetOrAddBoolean("montecarlo","polarized",false);
  acceleration = pin->GetOrAddBoolean("montecarlo","acceleration",false);
  time_acc = pin->GetOrAddBoolean("montecarlo","time_acc",false);
  raytrace_flag = pin->GetOrAddBoolean("montecarlo", "raytrace", false);
  if (raytrace_flag)
    general_mover_flag = true;
  else
    general_mover_flag = pin->GetOrAddBoolean("montecarlo","general_mover",false);
  scattering_meth = GetScatteringFlag(pin->GetOrAddString("montecarlo","scattering",
                                                          "none"));
                                                          
  nuser_var = 0; // Initialize photon user variables to zero

  // Set photon numbers
  nphtot = pin->GetInteger("montecarlo","nphot");
  cadence = pin->GetOrAddInteger("montecarlo","cadence",nphtot);
  nout = nphtot/cadence;
  max_phots_init = 10000;
#ifdef MPI_PARALLEL
  max_list_size = cadence/(Globals::nranks)+1;
#else
  max_list_size = cadence;
#endif

  checkscat = pin->GetOrAddInteger("montecarlo","checkscat",10000);
  checkmove = pin->GetOrAddInteger("montecarlo","checkmove",10000);

  // Initialize user MonteCarlo data before initializing MonteCarloBlocks
  // Should be caleld before Output constuctor
  InitUserMonteCarloData(pin);

  // Initialize output
  pmcout = new MCOutput(this,pin);

  // Create and intitialize randon number generator
  iseed = pin->GetInteger("montecarlo","iseed");

  // Initialize arrays sizes for photon instances
  Photon::Initialize(this,pin);

  // Initialize ncells and broadcast
  if (Globals::my_rank == 0) {
    ncells = pmesh->GetTotalCells();
  }
#ifdef MPI_PARALLEL
  // then broadcasts it
  MPI_Bcast(&ncells, sizeof(int64_t), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  // Create MonteCarloBlock for each MeshBlock for this process
  nblocal = pmy_mesh->nblocal;
  nbtotal = pmy_mesh->nbtotal;
  int root_level = pmy_mesh->root_level;
  my_blocks.NewAthenaArray(nblocal);
  for (int i=0; i<nblocal; i++) {
    MeshBlock *pmb = pmy_mesh->my_blocks(i);
    my_blocks(i) = new MonteCarloBlock(pmb, NULL, this, pin);
    pmb->pmy_mcb = my_blocks(i);
    // Set neighbors for photon class
    // SWD: ideally moved to mesh constructor, but awkward
    int nrbx1 = pmy_mesh->mesh_size.nx1/pmb->block_size.nx1;
    int nrbx2 = pmy_mesh->mesh_size.nx2/pmb->block_size.nx2;
    int nrbx3 = pmy_mesh->mesh_size.nx3/pmb->block_size.nx3;
    my_blocks(i)->pphot->LinkNeighbors(pmy_mesh->tree, nrbx1, nrbx2, nrbx3, root_level);
  }
}

//----------------------------------------------------------------------------------------
//! destructor

MonteCarlo::~MonteCarlo() {

  delete pmcout;

}

//----------------------------------------------------------------------------------------
//! \fn enum AbsorptionOpacityFlag GetAbsorptionOpacityFlag(std::string input_string)
//! \brief set absorption opacity flag

enum AbsorptionOpacityFlag GetAbsorptionOpacityFlag(std::string input_string) {
  if (input_string == "user") {
    return ABSUSER;
  } else if (input_string == "none") {
    return ABSNONE;
  } else if (input_string == "freefree") {
    return ABSFF;
  } else if (input_string == "dust") {
    return ABSDUST;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetAbsorptionOpacityFlag" << std::endl
        << "Input string=" << input_string << " not valid absorption opacity"
        << std::endl;
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
//! \fn enum AbsorptionMethodFlag GetAbsorptionMethodFlag(std::string input_string)
//! \brief set absorption method flag

enum AbsorptionMethodFlag GetAbsorptionMethodFlag(std::string input_string) {
  if (input_string == "weight") {
    return ABSWEIGHT;
  } else if (input_string == "prob") {
    return ABSPROB;
  } else if (input_string == "tau") {
    return ABSTAU;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetAbsorptionMethodFlag" << std::endl
        << "Input string=" << input_string << " not valid absorption method" << std::endl;
    ATHENA_ERROR(msg);
  }

}

//----------------------------------------------------------------------------------------
//! \fn enum ScatteringFlag GetScatteringFlag(std::string input_string)
//! \brief set scatering flag

enum ScatteringFlag GetScatteringFlag(std::string input_string) {
  if (input_string == "user") {
    return SCATUSER;
  } else if (input_string == "none") {
    return SCATNONE;
  } else if (input_string == "isotropic") {
    return SCATISO;
  } else if (input_string == "thomson") {
    return SCATTHOM;
  } else if (input_string == "compton") {
    return SCATCOMP;
  } else if (input_string == "resonance") {
    return SCATRES;
  } else if (input_string == "dust") {
    return SCATDUST;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetScatteringFlag" << std::endl
        << "Input string=" << input_string << " not valid scattering type" << std::endl;
    ATHENA_ERROR(msg);
  }

}

//----------------------------------------------------------------------------------------
//! \fn enum EmissionFlag GetEmissionFlag(std::string input_string)
//! \brief set emission flag

enum EmissionFlag GetEmissionFlag(std::string input_string) {

  if (input_string == "none") {
    return EMISNONE;
  } else if (input_string == "user") {
    return EMISUSER;
  } else if (input_string == "freefree") {
    return EMISFF;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetEmissionFlag" << std::endl
        << "Input string=" << input_string << " not valid emission type" << std::endl;
    ATHENA_ERROR(msg);
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserMCBoundaryFunction(enum BoundaryFace dir,
//!       BValHydro_t my_bc)
//!  \brief Enroll a user-defined monte carlo boundary function

void MonteCarlo::EnrollUserMCBoundaryFunction(enum BoundaryFace dir, MCBValFunc_t my_bc) {
  std::stringstream msg;
  if (dir<0 || dir>5) {
    msg << "### FATAL ERROR in EnrollMCBoundaryCondition function" << std::endl
        << "dirName = " << dir << " not valid" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (mc_bcs[dir]!=MC_USER_BNDRY) {
    msg << "### FATAL ERROR in EnrollUserMCBoundaryFunction" << std::endl
        << "The boundary condition flag must be set to the string 'user' in the "
        << " <mesh> block in the input file to use user-enrolled BCs" << std::endl;
    ATHENA_ERROR(msg);
  }
  BoundaryFunction_[dir]=my_bc;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserEmissionInitialization(EmisFunc_t emissfunc)
//! \brief Enroll a user-defined function for initializing emission methods

void MonteCarlo::EnrollUserEmissionInitialization(EmisFunc_t emissfunc) {

  InitEmission = emissfunc;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserGetTemperature(TempFunc_t tempfunc)
//! \brief Enroll a user-defined function for computing temperature

void MonteCarlo::EnrollUserGetTemperature(TempFunc_t tempfunc) {

  GetTemperature = tempfunc;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserWorkInMove(UserMoveFunc_t userfunc)
//! \brief Enroll a user-defined condition to be called during photon moves

void MonteCarlo::EnrollUserWorkInMove(UserMoveFunc_t userfunc) {

  UserWorkInMove = userfunc;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserScatteringFunction(ScatFunc_t scatfunc)
//! \brief Enroll a user-defined scattering function

void MonteCarlo::EnrollUserScatteringFunction(ScatFunc_t scatfunc) {

  UserScattering = scatfunc;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserOpacityFunction(OpacFunc_t opacfunc, bool abs)
//! \brief Enroll a user-defined opacity function

void MonteCarlo::EnrollUserOpacityFunction(OpacFunc_t opacfunc, bool abs) {

  if (abs)
    UserAbsorptionOpacity = opacfunc;
  else
    UserScatteringOpacity = opacfunc;
}

//----------------------------------------------------------------------------------------
//! \fn enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string)
//! \brief set boundary flag

enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string) {

  if (input_string == "periodic") {
    return MC_PERIODIC_BNDRY;
  } else if (input_string == "escape") {
    return MC_ESCAPE_BNDRY;
  } else if (input_string == "absorb") {
    return MC_ABSORB_BNDRY;
  } else if (input_string == "destroy") {
    return MC_DESTROY_BNDRY;
  } else if (input_string == "polar") {
    return MC_POLAR_BNDRY;
  } else if (input_string == "reflecting") {
    return MC_REFLECT_BNDRY;
  } else if (input_string == "user") {
    return MC_USER_BNDRY;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetMCBoundaryFlag" << std::endl
        << "Input string=" << input_string << " not valid boundary type" << std::endl;
    ATHENA_ERROR(msg);
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetDensity(MonteCarloBlock *pmcb)
//! \brief Make hard copy of density from MeshBlock to MonteCarloBlock.
//  Uses hard copy so that rho is always in cgs units

void MonteCarlo::GetDensity(MonteCarloBlock *pmcb) {

  // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        pmcb->rho(k,j,i) = pmcb->codetocgs_rho * pmcb->pmy_block->phydro->u(IDN,k,j,i);
      }}}
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetScalars(MonteCarloBlock *pmcb)
//! \brief Make hard copy of scalars from MeshBlock to MonteCarloBlock.

void MonteCarlo::GetScalars(MonteCarloBlock *pmcb) {

  // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        pmcb->scalars(k,j,i) = pmcb->pmy_block->pscalars->s(0,k,j,i);
      }}}
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetVelocities(MonteCarloBlock *pmcb)
//! \brief Make hard copy of velocites from MeshBlock to MonteCarloBlock.
//  Uses hard copy so that velocities is always fraction of speed of light

void MonteCarlo::GetVelocity(MonteCarloBlock *pmcb) {

  // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        Real rho = pmcb->pmy_block->phydro->u(IDN,k,j,i);
        pmcb->vel(0,k,j,i) = pmcb->codetocgs_vel *
          pmcb->pmy_block->phydro->u(IM1,k,j,i) / rho;
        pmcb->vel(1,k,j,i) = pmcb->codetocgs_vel *
          pmcb->pmy_block->phydro->u(IM2,k,j,i) / rho;
        pmcb->vel(2,k,j,i) = pmcb->codetocgs_vel *
          pmcb->pmy_block->phydro->u(IM3,k,j,i) / rho;
      }}}
}

//----------------------------------------------------------------------------------------
//! \fn void DefaultGetTemperature(MonteCarloBlock *pmcb)
//! \brief default function for computing temperature if no user function provided.
//  Assumes EOS of from P=RTd.

void DefaultGetTemperature(MonteCarloBlock *pmcb) {

  Real rideal = 8.314e7;
  Hydro* phydro = pmcb->pmy_block->phydro;

   // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  // get pressure
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        pmcb->tgas(k,j,i) = pmcb->codetocgs_tgas * phydro->w(IEN,k,j,i) /
                            phydro->w(IDN,k,j,i)/rideal;
      }}}
}

/*
//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::SendMonteCarloBlocks(int dest)
//! \brief send all monte carlo blocks to another process

void MonteCarlo::SendMonteCarloBlocks(int dest) {

#ifdef MPI_PARALLEL
  // Count number of blocks
  int blcnt=0;
  int head_buf[10];
  blcnt = 0;
  MonteCarloBlock *pmcb = pblock;
  while (pmcb != nullptr) {
    blcnt++;
    pmcb = pmcb->next;
  }
  // Send number of monte carlo block and dimensions
  head_buf[0] = blcnt;
  head_buf[1] = pblock->nx1;
  head_buf[2] = pblock->nx2;
  head_buf[3] = pblock->nx3;
  head_buf[4] = pblock->is;
  head_buf[5] = pblock->ie;
  head_buf[6] = pblock->js;
  head_buf[7] = pblock->je;
  head_buf[8] = pblock->ks;
  head_buf[9] = pblock->ke;

  MPI_Request send_rq;
  unsigned int tag = 0; //temporary
  MPI_Isend(head_buf,10,MPI_INT,dest,tag,MPI_COMM_WORLD,&send_rq);
  MPI_Wait(&send_rq, MPI_STATUS_IGNORE);
#endif

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::ReveiveMonteCarloBlocks(ParameterInput *pin, int source)
//! \brief initialize monte carlo blocks from another process

void MonteCarlo::ReceiveMonteCarloBlocks(ParameterInput *pin, int source) {

#ifdef MPI_PARALLEL
  // Receive number and dimensions of monte carlo blocks
  int head_buf[10];
  MPI_Request recv_rq;
  unsigned int tag = 0; //temporary
  MPI_Irecv(head_buf,10,MPI_INT,source,tag,MPI_COMM_WORLD,&recv_rq);
  MPI_Wait(&recv_rq, MPI_STATUS_IGNORE);

  int blcnt = head_buf[0];
  MCBlockSize blocksize;
  blocksize.nx1 = head_buf[1];
  blocksize.nx2 = head_buf[2];
  blocksize.nx3 = head_buf[3];
  blocksize.is = head_buf[4];
  blocksize.ie = head_buf[5];
  blocksize.js = head_buf[6];
  blocksize.je = head_buf[7];
  blocksize.ks = head_buf[8];
  blocksize.ke = head_buf[9];

  // create monte carlo blocks to receive data
  MonteCarloBlock *pfirst, *plast=nullptr;

  for(int i=0; i<blcnt; ++i) {
    pblock = new MonteCarloBlock(nullptr, &blocksize, this, pin);
    if (plast == nullptr)
      pfirst = pblock;
    else
      plast->next = pblock;
    plast = pblock;
  }
  pblock = pfirst;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::SendMonteCarloData(int dest)
//! \brief send all monte carlo blocks to another process

void MonteCarlo::SendMonteCarloData(int dest) {
#ifdef MPI_PARALLEL
  // Send data for each block
  MonteCarloBlock *pmcb = pblock;
  Real *send_buf;
  int size = 3; //tgas,rho,vol
  if (boosts) size+=3;
  if (computedmin) size+=1; //dmin array
  if (NSCALARS > 0) size+=1; //scalars array
  size *= (pmcb->nx1*pmcb->nx2*pmcb->nx3); // all blocks have same size
  size += pmcb->nx1+1; size += pmcb->nx2+1; size+= pmcb->nx3+1;
  send_buf = new Real[size];
  MPI_Request send_rq;
  unsigned int tag = 1;
  while (pmcb != nullptr) {
    int p=0;
    BufferUtility::PackData(pmcb->rho,send_buf,pmcb->is,pmcb->ie,pmcb->js,pmcb->je,
                            pmcb->ks,pmcb->ke,p);
    BufferUtility::PackData(pmcb->tgas,send_buf,pmcb->is,pmcb->ie,pmcb->js,pmcb->je,
                            pmcb->ks,pmcb->ke,p);
    if (boosts)
      BufferUtility::PackData(pmcb->vel,send_buf,0,2,pmcb->is,pmcb->ie,pmcb->js,pmcb->je,
                              pmcb->ks,pmcb->ke,p);
    BufferUtility::PackData(pmcb->pcoord->vol,send_buf,pmcb->is,pmcb->ie,pmcb->js,
                            pmcb->je,pmcb->ks,pmcb->ke,p);
    if (computedmin)
      BufferUtility::PackData(pmcb->pcoord->dmin,send_buf,pmcb->is,pmcb->ie,pmcb->js,
                              pmcb->je,pmcb->ks,pmcb->ke,p);
    if (NSCALARS > 0)
      BufferUtility::PackData(pmcb->scalars,send_buf,pmcb->is,pmcb->ie,pmcb->js,pmcb->je,
                              pmcb->ks,pmcb->ke,p);
    for (int i=pmcb->is; i<=pmcb->ie+1; ++i)
      send_buf[p++] = pmcb->pcoord->x1f(i);
     for (int i=pmcb->js; i<=pmcb->je+1; ++i)
      send_buf[p++] = pmcb->pcoord->x2f(i);
    for (int i=pmcb->ks; i<=pmcb->ke+1; ++i)
      send_buf[p++] = pmcb->pcoord->x3f(i);
    MPI_Isend(send_buf,size,MPI_ATHENA_REAL,dest,tag++,MPI_COMM_WORLD,&send_rq);
    MPI_Wait(&send_rq, MPI_STATUS_IGNORE);
    pmcb=pmcb->next;
  }
  delete send_buf;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::ReceiveMonteCarloData(int source)
//! \brief initialize monte carlo data from another process

void MonteCarlo::ReceiveMonteCarloData(int source) {

#ifdef MPI_PARALLEL
  MonteCarloBlock *pmcb=pblock;
  Real *recv_buf;
  int size = 3; //tgas,rho,vol
  if (boosts) size+=3;
  if (computedmin) size+=1; //dmin array
  if (NSCALARS > 0) size+=1; //scalars array
  size *= (pmcb->nx1*pmcb->nx2*pmcb->nx3); // all blocks have same size
  size += pmcb->nx1+1; size += pmcb->nx2+1; size+= pmcb->nx3+1;
  recv_buf = new Real[size];
  MPI_Request recv_rq;
  unsigned int tag = 1;
  while (pmcb != nullptr) {
    MPI_Irecv(recv_buf,size,MPI_ATHENA_REAL,source,tag++,MPI_COMM_WORLD,&recv_rq);
    MPI_Wait(&recv_rq, MPI_STATUS_IGNORE);
    int p=0;
    BufferUtility::UnpackData(recv_buf, pmcb->rho, pmcb->is, pmcb->ie, pmcb->js,
                                pmcb->je, pmcb->ks, pmcb->ke, p);
    BufferUtility::UnpackData(recv_buf, pmcb->tgas, pmcb->is, pmcb->ie, pmcb->js,
                                pmcb->je, pmcb->ks, pmcb->ke, p);
    if (boosts)
      BufferUtility::UnpackData(recv_buf, pmcb->vel, 0, 2, pmcb->is, pmcb->ie, pmcb->js,
                                  pmcb->je, pmcb->ks, pmcb->ke, p);
    BufferUtility::UnpackData(recv_buf, pmcb->pcoord->vol, pmcb->is, pmcb->ie, pmcb->js,
                                pmcb->je, pmcb->ks, pmcb->ke, p);
    if (computedmin)
      BufferUtility::UnpackData(recv_buf, pmcb->pcoord->dmin, pmcb->is, pmcb->ie,
                                pmcb->js, pmcb->je, pmcb->ks, pmcb->ke, p);
    if (NSCALARS > 0)
      BufferUtility::UnpackData(recv_buf, pmcb->scalars, pmcb->is, pmcb->ie,
                                pmcb->js, pmcb->je, pmcb->ks, pmcb->ke, p);
    for (int i=pmcb->is; i<=pmcb->ie+1; ++i)
      pmcb->pcoord->x1f(i) = recv_buf[p++];
    for (int i=pmcb->js; i<=pmcb->je+1; ++i)
      pmcb->pcoord->x2f(i) = recv_buf[p++];
    for (int i=pmcb->ks; i<=pmcb->ke+1; ++i)
      pmcb->pcoord->x3f(i) = recv_buf[p++];
    // initialize emission array
    if (InitEmission != nullptr)
      pmcb->minweight *= InitEmission(pmcb);
    if (acceleration && !(pmcb->coherent_scattering) && !(scattering_meth == SCATRES))
      InitializeAccelerationOpacity(pmcb);
    pmcb=pmcb->next;
  }
  delete recv_buf;
#endif
}
*/
//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::SendMonteCarloSpectra(int dest)
//! \brief send all monte carlo spectra to another process

void MonteCarlo::SendMonteCarloSpectra(int dest) {
#ifdef MPI_PARALLEL
  Spectrum *pspec = pmcout->pspec;

  int maxsize = 0;
  while (pspec != nullptr) {
    int ne = pspec->range.ne;
    int ncth = pspec->range.ncth;
    int nphi = pspec->range.nphi;
    int size = 2;
    if (pspec->polarized) {
      size += 4;
    }
    size *= (ne*ncth*nphi);
    maxsize = (size > maxsize) ? size : maxsize;
    pspec = pspec->next;
  }
  Real *send_buf;
  send_buf = new Real[maxsize];
  MPI_Request send_rq;
  unsigned int tag = 100;

  pspec = pmcout->pspec;
  while (pspec != nullptr) {
    int p=0;
    int ne = pspec->range.ne;
    int ncth = pspec->range.ncth;
    int nphi = pspec->range.nphi;
    int size = 2;
    if (pspec->polarized)
      size += 4;
    size *= (ne*ncth*nphi);
    ne--; ncth--; nphi--;
    BufferUtility::PackData(pspec->intensity,send_buf,0,ne,0,ncth,0,nphi,p);
    BufferUtility::PackData(pspec->intensity_sq,send_buf,0,ne,0,ncth,0,nphi,p);
    if (pspec->polarized) {
      BufferUtility::PackData(pspec->stokesq,send_buf,0,ne,0,ncth,0,nphi,p);
      BufferUtility::PackData(pspec->stokesq_sq,send_buf,0,ne,0,ncth,0,nphi,p);
      BufferUtility::PackData(pspec->stokesu,send_buf,0,ne,0,ncth,0,nphi,p);
      BufferUtility::PackData(pspec->stokesu_sq,send_buf,0,ne,0,ncth,0,nphi,p);
    }
    MPI_Isend(send_buf,size,MPI_ATHENA_REAL,dest,tag++,MPI_COMM_WORLD,&send_rq);
    MPI_Wait(&send_rq, MPI_STATUS_IGNORE);
    pspec = pspec->next;
  }
  delete send_buf;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::MakeOutputs()
//! \brief write MonteCarlo outputs

void MonteCarlo::MakeOutputs() {

  nphrun = nphtot;
  pmcout->OutputSpectrum(this);
  pmcout->OutputPhotonList(nblock);
  pmcout->OutputTrajectoryList();

}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::ReceiveMonteCarloSpectra(int dest)
//! \brief receive monte carlo spectra from another process

void MonteCarlo::ReceiveMonteCarloSpectra(int source) {
#ifdef MPI_PARALLEL
  Spectrum *pspec = pmcout->pspec;
  int maxsize = 0;
  while (pspec != nullptr) {
    int ne = pspec->range.ne;
    int ncth = pspec->range.ncth;
    int nphi = pspec->range.nphi;
    int size = 2;
    if (pspec->polarized)
      size += 4;
    size *= (ne*ncth*nphi);
    maxsize = (size > maxsize) ? size : maxsize;
    pspec = pspec->next;
  }
  Real *recv_buf;
  recv_buf = new Real[maxsize];
  MPI_Request recv_rq;
  unsigned int tag = 100; // temporary

  pspec = pmcout->pspec;
  while (pspec != nullptr) {
    int ne = pspec->range.ne;
    int ncth = pspec->range.ncth;
    int nphi = pspec->range.nphi;
    int size = 2;
    if (pspec->polarized)
      size += 4;
    size *= (ne*ncth*nphi);
    ne--; ncth--; nphi--;
    MPI_Irecv(recv_buf,size,MPI_ATHENA_REAL,MPI_ANY_SOURCE,tag++,MPI_COMM_WORLD,&recv_rq);
    MPI_Wait(&recv_rq, MPI_STATUS_IGNORE);
    Spectrum *ptemp = new Spectrum(pspec);
    int p=0;
    BufferUtility::UnpackData(recv_buf,ptemp->intensity,0,ne,0,ncth,0,nphi,p);
    BufferUtility::UnpackData(recv_buf,ptemp->intensity_sq,0,ne,0,ncth,0,nphi,p);
    if (pspec->polarized) {
      BufferUtility::UnpackData(recv_buf,ptemp->stokesq,0,ne,0,ncth,0,nphi,p);
      BufferUtility::UnpackData(recv_buf,ptemp->stokesq_sq,0,ne,0,ncth,0,nphi,p);
      BufferUtility::UnpackData(recv_buf,ptemp->stokesu,0,ne,0,ncth,0,nphi,p);
      BufferUtility::UnpackData(recv_buf,ptemp->stokesu_sq,0,ne,0,ncth,0,nphi,p);
    }
    pspec->AddSpectrum(ptemp);
    pspec = pspec->next;
    delete ptemp;
  }
  delete recv_buf;
#endif
}

/*
//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::CollectMoments(void)
//! \brief collect moments from other processes for output

void MonteCarlo::CollectMoments(void) {

  if (origin < 0) {
    // Set Moments to zero on origin blocks
    MonteCarloBlock *pmcb=pblock;
    while (pmcb != nullptr) {
      pmcb->ResetMoments();
      pmcb=pmcb->next;
    }
    // Retrieve moments from destination processes
    for(int i=0; i<nderv; ++i) {
      ReceiveMoments(derv[i],true);
    }
    pmcb=pblock;
    while (pmcb != nullptr) {
      pmcb->nphdone = nphlist[pmcb->myblockid];
      pmcb->NormalizeMoments(true);
      pmcb=pmcb->next;
    }
  } else {
    // Return moments to origin
    SendMoments(origin);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::SendMoments(int dest)
//! \brief receive momdents

void MonteCarlo::SendMoments(int dest) {
#ifdef MPI_PARALLEL
  // Send data for each block
  MonteCarloBlock *pmcb = pblock;
  Real *send_buf;
  int size = (NMOM-3) * (pmcb->nx1*pmcb->nx2*pmcb->nx3);
  send_buf = new Real[size];
  MPI_Request send_rq;
  unsigned int tag = 1000; // temporary
  while (pmcb != nullptr) {
    int p=0;
    BufferUtility::PackData(pmcb->moments,send_buf,0,(NMOM-4),pmcb->is,pmcb->ie,pmcb->js,
                            pmcb->je,pmcb->ks,pmcb->ke,p);
    MPI_Isend(send_buf,size,MPI_ATHENA_REAL,dest,tag,MPI_COMM_WORLD,&send_rq);
    MPI_Wait(&send_rq, MPI_STATUS_IGNORE);
    pmcb=pmcb->next;
  }
  delete send_buf;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::ReceiveMoments(int source, bool sum_moments)
//! \brief receive moments from source

void MonteCarlo::ReceiveMoments(int source, bool sum_moments) {
#ifdef MPI_PARALLEL
  // Receive data from each block
  MonteCarloBlock *pmcb=pblock;
  Real *recv_buf;
  int size = (NMOM-3) * (pmcb->nx1*pmcb->nx2*pmcb->nx3);
  recv_buf = new Real[size];
  MPI_Request recv_rq;
  unsigned int tag = 1000; //temporary
  while (pmcb != nullptr) {
    MPI_Irecv(recv_buf,size,MPI_ATHENA_REAL,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&recv_rq);
    MPI_Wait(&recv_rq, MPI_STATUS_IGNORE);
    int p=0;
    if (sum_moments)
      BufferUtility::UnpackDataSum(recv_buf, pmcb->moments,0,(NMOM-4),pmcb->is,pmcb->ie,
                                   pmcb->js,pmcb->je,pmcb->ks,pmcb->ke,p);
    else
      BufferUtility::UnpackData(recv_buf, pmcb->moments,0,(NMOM-4),pmcb->is,pmcb->ie,
                                pmcb->js,pmcb->je,pmcb->ks,pmcb->ke,p);
    pmcb=pmcb->next;
  }
  delete recv_buf;
#endif
}
*/

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::Initialize(ParameterInput *pinput)
//! \brief initialize grid data in each monte carlo block

void MonteCarlo::Initialize(ParameterInput *pin) {

  if (MONTE_CARLO_STATIC) {
    // initialize timing parameters if static calculation
    dt = pin->GetOrAddReal("montecarlo","dt",1.);
    tmax = pin->GetOrAddReal("montecarlo","tmax",HUGE_NUMBER);
  } else if (MONTE_CARLO_DYNAMIC) {
    dt = tmax = pmy_mesh->dt;
  }

  // Set number of photons per montecarloblock
  nblock = nphtot / nbtotal;
  nphtot = nblock * nbtotal; // reset nblock to be multiple of nbtotal
  int nchunk = nblock / 2;
  if (GetTemperature == nullptr)
    GetTemperature = DefaultGetTemperature;

  // Initialize monte carlo blocks
  for (int i=0; i<nblocal; i++) {
    MonteCarloBlock *pmcb = my_blocks(i);
    // Initialize variables over all blocks
    GetDensity(pmcb);
    GetTemperature(pmcb);
    if (boosts) GetVelocity(pmcb);
    if (InitEmission != nullptr) {
      pmcb->minweight *= InitEmission(pmcb);
    }
    // set photons to be computed and counters
    pmcb->nphremain = nblock;
    pmcb->nphdone = 0;
    pmcb->nchunk = nchunk;
    pmcb->nscat = pmcb->nesc = pmcb->nabs = pmcb->ndes = 0;
    // Call problem generators for Monte Carlo
    pmcb->MonteCarloProblemGenerator(pin);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::RunStaticMonteCarlo(Outputs *pouts, Mesh *pmesh,
//!                                          ParameterInput *pinput)
//! \brief Finish Initialization of MonteCarloBlocks and run steady-state MC calculation

void MonteCarlo::RunStaticMonteCarlo(Outputs *pouts, Mesh *pmesh,
                                     ParameterInput *pinput) {

  // initialize counter
  int ntot = 0;

  bool photons_remain = true; // True if photons on any process
  while(photons_remain) {

    for(int nb=0; nb<nblocal; ++nb){
      my_blocks(nb)->TransferPhotonsOnBlock();
    }
    photons_remain = CheckAndBroadCastPhotonsRemaining();
  }

  // Report diagnostic results from all blocks
  int nesc = 0, nabs = 0, ndes = 0, nscat = 0;
  for(int nb=0; nb<nblocal; ++nb){
    MonteCarloBlock *pmcb = my_blocks(nb);
    nesc += pmcb->nesc;
    nabs += pmcb->nabs;
    ndes += pmcb->ndes;
    nscat += pmcb->nscat;
    ntot += pmcb->nphdone;
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,&nesc,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&nabs,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&ndes,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&nscat,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&ntot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (Globals::my_rank == 0) {
    std::cout  << "ntot, nesc, nabs, ndes, nscat: "
               << ' ' << ntot << ' ' << nesc
               << ' ' << nabs << ' ' << ndes << ' ' << nscat << ' ';
    if (ntot > 0)
      std::cout << static_cast<Real>(nscat)/static_cast<Real>(ntot) << std::endl;
    else
      std::cout << 0. << std::endl;
  }

  if (pmcout->moments) {
    // normalize moments for output
    for(int nb=0; nb<nblocal; ++nb){
      MonteCarloBlock *pmcb = my_blocks(nb);
      pmcb->NormalizeMoments(true,static_cast<Real>(ntot));
    }
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::RunStaticMonteCarlo()
//! \brief Checks if photons remain on any process and sends message to rank 0
//!        process if none remaining.

bool MonteCarlo::CheckAndBroadCastPhotonsRemaining() {

  // Send photons from all blocks
  for(int nb=0; nb<nblocal; ++nb){
    my_blocks(nb)->pphot->SendToNeighbors();
  }

  // Receive photons from all blocks
  bool complete = false;
  while(!complete) {
    complete = true;
    for(int nb=0; nb<nblocal; ++nb) {
      bool success = my_blocks(nb)->pphot->ReceiveFromNeighbors();
      if (!success)
        complete = false;
    }
  }
  // Clear Boundaries
  for(int nb=0; nb<nblocal; ++nb)
    my_blocks(nb)->pphot->ClearBoundary();

  // Check if photons have completed
  int nremain=0,nprop=0;
  for(int nb=0; nb<nblocal; ++nb){
    MonteCarloBlock *pmcb = my_blocks(nb);
    nremain += pmcb->nphremain;
    nprop += pmcb->pphot->nphot;
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,&nprop,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&nremain,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  bool active;
  if ((nremain > 0) || (nprop > 0)) {
    active = true;
  } else {
    active = false;
  }
  return active;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::RunDynamicMonteCarlo(Outputs *pouts, Mesh *pmesh,
//!                                          ParameterInput *pinput)
//! \brief Run one step in dynamic MC calculation

void MonteCarlo::RunDynamicMonteCarlo(Outputs *pouts, Mesh *pmesh,
                                     ParameterInput *pinput) {

  // initialize counter
  int ntot = 0;

  dt = tmax = pmy_mesh->dt;

  // Rest monte carlo blocks
  for (int i=0; i<nblocal; i++) {
    MonteCarloBlock *pmcb = my_blocks(i);
    // Initialize variables over all blocks
    GetDensity(pmcb);
    GetTemperature(pmcb);
    if (boosts) GetVelocity(pmcb);
    if (InitEmission != nullptr) {
      pmcb->minweight *= InitEmission(pmcb);
    }
    // set photons to be computed and counters
    pmcb->nphremain = nblock;
    pmcb->nphdone = 0;
    pmcb->nscat = pmcb->nesc = pmcb->nabs = pmcb->ndes = 0;
  }

  if (pmcout->moments) {
    // reset moments at start of new timestep
    for(int nb=0; nb<nblocal; ++nb){
      MonteCarloBlock *pmcb = my_blocks(nb);
      pmcb->ResetMoments();
    }
  }

  bool photons_remain = true; // True if photons on any process
  while(photons_remain) {

    for(int nb=0; nb<nblocal; ++nb){
      my_blocks(nb)->TransferPhotonsOnBlock();
    }
    photons_remain = CheckAndBroadCastPhotonsRemaining();
  }

  // Report diagnostic results from all blocks
  int nesc = 0, nabs = 0, ndes = 0, nscat = 0;
  for(int nb=0; nb<nblocal; ++nb){
    MonteCarloBlock *pmcb = my_blocks(nb);
    nesc += pmcb->nesc;
    nabs += pmcb->nabs;
    ndes += pmcb->ndes;
    nscat += pmcb->nscat;
    ntot += pmcb->nphdone;
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,&nesc,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&nabs,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&ndes,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&nscat,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&ntot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  if (Globals::my_rank == 0) {
    std::cout  << "ntot, nesc, nabs, ndes, nscat: "
               << ' ' << ntot << ' ' << nesc
               << ' ' << nabs << ' ' << ndes << ' ' << nscat << ' ';
    if (ntot > 0)
      std::cout << static_cast<Real>(nscat)/static_cast<Real>(ntot) << std::endl;
    else
      std::cout << 0. << std::endl;
  }

  if (pmcout->moments) {
    // normalize moments for output
    for(int nb=0; nb<nblocal; ++nb){
      MonteCarloBlock *pmcb = my_blocks(nb);
      pmcb->NormalizeMoments(true,static_cast<Real>(ntot));
    }
  }

}


//----------------------------------------------------------------------------------------
//! MCRandom constructor, builds Athena++ random number generator
//  current implementation is wrapper for gsl function

MCRandom::MCRandom(int iseed) {

#if RAN3
  r3seed = static_cast<long>(-iseed);
  ran3(&r3seed);
#else
  dev = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(dev, iseed);
#endif

}

//----------------------------------------------------------------------------------------
//! destructor

MCRandom::~MCRandom() {

}

Real MCRandom::uniform() {

#if RAN3
  return ran3(&r3seed);
#else
  return static_cast<Real>(gsl_rng_uniform(dev));
#endif
}

Real MCRandom::chisquare(int n) {
#if RAN3
  std::stringstream msg;
  msg << "### FATAL ERROR in MCRandom::chisquare" << std::endl
      << "Not supported with RAN3 random number generator" << std::endl;
  ATHENA_ERROR(msg);
#else
  return static_cast<Real>(gsl_ran_chisq(dev,n));
#endif

}

// Used by ran3
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

Real MCRandom::ran3(long *idum) {

  static int inext,inextp;
  static long ma[56];       // The value 56 (range ma[1..55]) is special and
  static int iff=0;         // should not be modified; see Knuth.

  long mj,mk;
  int i,ii,k;

  if (*idum < 0 || iff == 0) {    //Initialization.
    iff=1;
    mj=labs(MSEED-labs(*idum));   // Initialize ma[55] using the seed idum
    mj %= MBIG;                   // and the large number MSEED.
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {         //  Now initialize the rest of the table,
      ii=(21*i) % 55;             //  in a slightly random order,with
      ma[ii]=mk;                  //  numbers that are not especially random.
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)    // We randomize them by "warming up the generator."
      for (i=1;i<=55;i++) {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;     // Prepare indices for our first generated number.
    inextp=31;   //  The constant 31 is special; see Knuth.
    *idum=1;
  }
  // Here is where we start, except on initialization.
  if (++inext == 56) inext=1;     // Increment inext and inextp, wrapping
  if (++inextp == 56) inextp=1;   // around  56 to 1.

  mj=ma[inext]-ma[inextp];   // Generate a new random number subtractively.
  if (mj < MZ) mj += MBIG;   // Be sure that it is in range.
  ma[inext]=mj;              // Store it,
  return static_cast<Real>(mj*FAC);             // and output the derived uniform deviate.
}
