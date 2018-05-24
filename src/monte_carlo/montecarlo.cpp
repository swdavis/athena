//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//  \brief implementation of functions in class MonteCarlo, MCRandom

#include <gsl/gsl_randist.h>
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "montecarlo.hpp"

#include "../globals.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../utils/buffer_utils.hpp"

// constructor, initializes data structures and parameters

MonteCarlo::MonteCarlo(ParameterInput *pin, Mesh *pmesh) {

  MonteCarloBlock *pfirst;
  std::stringstream msg;

  pmy_mesh = pmesh;

  InitEmission=NULL;
  GetTemperature=NULL;
 
  // Set flags that control emission, absorption and scattering
  emission_meth = GetEmissionFlag(pin->GetOrAddString("montecarlo","emission","error"));
  if (emission_meth ==  EMISFF) {
    InitEmission = InitializeEmissionFreeFree;
  }
  absorption_meth = GetAbsorptionFlag(pin->GetOrAddString("montecarlo","absorption",
                                                          "error"));
  scattering_meth = GetScatteringFlag(pin->GetOrAddString("montecarlo","scattering",
                                                          "error"));
  // read bc flags for each of the 6 boundaries.
  mc_bcs[INNER_X1] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ix1_mc_bc","escape"));
  mc_bcs[OUTER_X1] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ox1_mc_bc","escape"));
  mc_bcs[INNER_X2] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ix2_mc_bc","escape"));
  mc_bcs[OUTER_X2] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ox2_mc_bc","escape"));
  mc_bcs[INNER_X3] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ix3_mc_bc","escape"));
  mc_bcs[OUTER_X3] = GetMCBoundaryFlag(pin->GetOrAddString("mesh","ox3_mc_bc","escape"));

  // Initialize output
  pmcout = new MCOutput(this,pin);
  //pmcout->CheckFace(mc_bcs);

  lorentz_transform = pin->GetOrAddBoolean("montecarlo","lorentz_transform",false);
  polarized = pin->GetOrAddBoolean("montecarlo","polarized",false);

  // Create and intitialize randon number generator
  iseed = pin->GetInteger("montecarlo","iseed");

  // Initialize ncells and broadcast
  if (Globals::my_rank == 0) {
    ncells = pmesh->GetTotalCells();
#ifdef MPI_PARALLEL
    // then broadcasts it
    MPI_Bcast(&ncells, sizeof(int64_t), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
  }
  // Initialize all montecarlo blocks to correspond to mesh blocks
  mcranks = pin->GetOrAddInteger("montecarlo","mcranks",0);
  int nmesh = (Globals::nranks)-mcranks;
  if(Globals::my_rank < nmesh) {
    // use mesh blocks initialized on my process
    int myblockid = pmesh->nslist[Globals::my_rank];
    MeshBlock *pmb = pmesh->pblock;
    pblock = new MonteCarloBlock(pmb, NULL, this, pin);
    //pblock->myblockid = myblockid;
    //pblock->nphremain = nphlist[myblockid++]; 
    pfirst = pblock;
    pmb=pmb->next;
    while (pmb != NULL)  {
      pblock->next = new MonteCarloBlock(pmb, NULL, this, pin);
      //pblock->myblockid = myblockid;
      //pblock->nphremain = nphlist[myblockid++]; 
      pblock = pblock->next;
      pmb=pmb->next;
    }
    pblock = pfirst;
    // set list of destination processes
    source = -1;
    ndest = 0;
    for(int i=nmesh; i<Globals::nranks; ++i) {
      if ((i % nmesh)==Globals::my_rank) 
        ndest++;
    }
    dest = new int[ndest];
    int j = 0;
    for(int i=nmesh; i<Globals::nranks; ++i) {
      if ((i % nmesh)==Globals::my_rank) {
        dest[j++] = i;
        SendMonteCarloBlocks(i);
      }
    }
  } else {
    // no mesh blocks on my process, copy from another process and source
    ndest = 0;
    dest = NULL;
    source = Globals::my_rank % nmesh;
    ReceiveMonteCarloBlocks(pin,source);
  }

  // set number of photons for each block
  nphot = pin->GetInteger("montecarlo","nphot");
  int nbtotal = pmesh->nbtotal;
  nphlist = new int[nbtotal];
  int nperb = nphot / nbtotal;
  for (int i=0; i<nbtotal; ++i)
    nphlist[i] = nperb;
  // ensure equal number of photons per mesh block
  nphot = 0;
  for (int i=0; i<nbtotal; ++i)
    nphot += nphlist[i];
  if(nphot < nbtotal) {
    msg << "### FATAL ERROR Monte Carlo Constructor" << std::endl
        << "Number of photons < number of mesh blocks= " << nphot << " "
        << nbtotal << std::endl;
    throw std::runtime_error(msg.str().c_str());
  } 
  // Divide photons between ranks and set blockid
  int *nranks;
  nranks = new int[nmesh];
  for(int i=0; i<nmesh; ++i)
    nranks[i]=0;
  for(int i=0; i<Globals::nranks; ++i) {
    nranks[i % nmesh]++;
  }
  int my_mesh_rank = Globals::my_rank % nmesh;
  int myblockid = pmesh->nslist[my_mesh_rank];
  MonteCarloBlock *pmcb = pblock;
  while (pmcb != NULL) {
    pmcb->myblockid = myblockid;
    pmcb->nphremain = nphlist[myblockid++]/nranks[my_mesh_rank];
    pmcb = pmcb->next;
  }
}

// destructor

MonteCarlo::~MonteCarlo() {

  delete pmcout;
  if (pblock != NULL) {
    while(pblock->next != NULL)
      delete pblock->next;
    delete pblock;
  }
  delete nphlist;

}

//----------------------------------------------------------------------------------------
//! \fn enum AbsorptionFlag GetAbsorptionFlag(std::string input_string)
//  \brief set absorption flag

enum AbsorptionFlag GetAbsorptionFlag(std::string input_string) {
  if (input_string == "user") {
    return ABSUSER;
  } else if (input_string == "none") {
    return ABSNONE;
  } else if (input_string == "freefree") {
    return ABSFF;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetAbsorptionFlag" << std::endl
        << "Input string=" << input_string << " not valid absorption type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}

//----------------------------------------------------------------------------------------
//! \fn enum ScatteringFlag GetScatteringFlag(std::string input_string)
//  \brief set scatering flag

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
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetAbsorptionFlag" << std::endl
        << "Input string=" << input_string << " not valid scattering type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}

//----------------------------------------------------------------------------------------
//! \fn enum EmissionFlag GetEmissionFlag(std::string input_string)
//  \brief set emission flag

enum EmissionFlag GetEmissionFlag(std::string input_string) {
  if (input_string == "user") {
    return EMISUSER;
  } else if (input_string == "freefree") {
    return EMISFF;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetEmissionFlag" << std::endl
        << "Input string=" << input_string << " not valid emission type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}

//----------------------------------------------------------------------------------------
//! \fn enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string)
//  \brief set boundary flag

enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string) {

  if (input_string == "periodic") {
    return MC_PERIODIC_BNDRY;
  } else if (input_string == "escape") {
    return MC_ESCAPE_BNDRY;
  } else if (input_string == "absorb") {
    return MC_ABSORB_BNDRY;
  } else if (input_string == "polar") {
    return MC_POLAR_BNDRY;
  } else if (input_string == "reflect") {
    return MC_REFLECT_BNDRY;
  } else if (input_string == "user") {
    return MC_USER_BNDRY;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetMCBoundaryFlag" << std::endl
        << "Input string=" << input_string << " not valid boundary type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}
//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetDensity(MonteCarloBlock *pmcb)
//  \brief Make hard copy of density from MeshBlock to MonteCarloBlock. Uses hard copy
//  so that rho is always in cgs units

void MonteCarlo::GetDensity(MonteCarloBlock *pmcb) {

  
  // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        pmcb->rho(k,j,i) = pmcb->codetocgs_rho * pmcb->pmy_block->phydro->u(IDN,k,j,i);
      }}}
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetVelocities(MonteCarloBlock *pmcb)
//  \brief Make hard copy of velocites from MeshBlock to MonteCarloBlock. Uses hard copy
//  so that velocities is always fraction of speed of light

void MonteCarlo::GetVelocity(MonteCarloBlock *pmcb) {

  
  // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        Real rho = pmcb->pmy_block->phydro->u(IDN,k,j,i);
        pmcb->vel(0,k,j,i) = pmcb->codetoc_vel*pmcb->pmy_block->phydro->u(IM1,k,j,i)/rho;
        pmcb->vel(1,k,j,i) = pmcb->codetoc_vel*pmcb->pmy_block->phydro->u(IM2,k,j,i)/rho;
        pmcb->vel(2,k,j,i) = pmcb->codetoc_vel*pmcb->pmy_block->phydro->u(IM3,k,j,i)/rho;     
        // transform to cartesian if not cartesian
      }}}
}

//----------------------------------------------------------------------------------------
//! \fn void DefaultGetTemperature(MonteCarloBlock *pmcb)
//  \brief default function for computing temperature if no user function provided.
//  Assumes that code values correspond to cgs with simple equation of state.

void DefaultGetTemperature(MonteCarloBlock *pmcb) {

  Real rideal = 8.314e7;
  Hydro* phydro = pmcb->pmy_block->phydro;

   // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu+1; ++i) {
        pmcb->tgas(k,j,i) = phydro->w(IEN,k,j,i)/phydro->w(IDN,k,j,i)/rideal;

      }}}

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserEmissionInitialization(EmisFunc_t emissfunc)
//  \brief Enroll a user-defined function for initializing emission methods

void MonteCarlo::EnrollUserEmissionInitialization(EmisFunc_t emissfunc) {
  InitEmission = emissfunc;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::EnrollUserGetTemperature(TempFunc_t tempfunc)
//  \brief Enroll a user-defined function for computing temperature

void MonteCarlo::EnrollUserGetTemperature(TempFunc_t tempfunc) {
  GetTemperature = tempfunc;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::SendMonteCarloBlocks(int dest)
//  \brief send all monte carlo blocks to another process

void MonteCarlo::SendMonteCarloBlocks(int dest) {

#ifdef MPI_PARALLEL 
  // Count number of blocks
  int blcnt=0;
  int head_buf[10];
  blcnt = 0;
  MonteCarloBlock *pmcb = pblock;
  while (pmcb != NULL) {
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
  unsigned int tag = 0;
  MPI_Isend(head_buf,10,MPI_INT,dest,tag,MPI_COMM_WORLD,&send_rq);
  MPI_Wait(&send_rq, MPI_STATUS_IGNORE);
#endif

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::ReveiveMonteCarloBlocks(ParameterInput *pin, int source)
//  \brief initialize monte carlo blocks from another process

void MonteCarlo::ReceiveMonteCarloBlocks(ParameterInput *pin, int source) {

#ifdef MPI_PARALLEL 
  // Receive number and dimensions of monte carlo blocks
  int head_buf[10];
  MPI_Request recv_rq;
  unsigned int tag = 0;
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

  // creat monte carlo blocks to receive data
  MonteCarloBlock *pfirst, *plast=NULL;
 
  for(int i=0; i<blcnt; ++i) {
    pblock = new MonteCarloBlock(NULL, &blocksize, this, pin);
    if (plast == NULL) 
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
//  \brief send all monte carlo blocks to another process

void MonteCarlo::SendMonteCarloData(int dest) {
#ifdef MPI_PARALLEL 
  // Send data for each block
  MonteCarloBlock *pmcb = pblock;
  Real *send_buf;
  int size = 3; //tgas,rho,vol
  if (lorentz_transform) size+=3;
  size *= (pmcb->nx1*pmcb->nx2*pmcb->nx3); // all blocks have same size
  size += pmcb->nx1+1; size += pmcb->nx2+1; size+= pmcb->nx3+1;
  send_buf = new Real[size];
  MPI_Request send_rq;
  unsigned int tag = 1;
  while (pmcb != NULL) {
    int p=0;
    BufferUtility::Pack3DData(pmcb->rho,send_buf,pmcb->is,pmcb->ie,pmcb->js,pmcb->je,pmcb->ks,pmcb->ke,p);
    BufferUtility::Pack3DData(pmcb->tgas,send_buf,pmcb->is,pmcb->ie,pmcb->js,pmcb->je,pmcb->ks,pmcb->ke,p);
    if (lorentz_transform)
      BufferUtility::Pack4DData(pmcb->vel,send_buf,0,2,pmcb->is,pmcb->ie,pmcb->js,pmcb->je,pmcb->ks,pmcb->ke,p);
    BufferUtility::Pack3DData(pmcb->pcoord->vol,send_buf,pmcb->is,pmcb->ie,pmcb->js,pmcb->je,pmcb->ks,pmcb->ke,p);
    for (int i=pmcb->is; i<=pmcb->ie; ++i) 
      send_buf[p++] = pmcb->pcoord->x1f(i);
     for (int i=pmcb->js; i<=pmcb->je; ++i) 
      send_buf[p++] = pmcb->pcoord->x2f(i);
    for (int i=pmcb->ks; i<=pmcb->ke; ++i) 
      send_buf[p++] = pmcb->pcoord->x3f(i);
    MPI_Isend(send_buf,size,MPI_ATHENA_REAL,dest,tag++,MPI_COMM_WORLD,&send_rq);
    MPI_Wait(&send_rq, MPI_STATUS_IGNORE);
    pmcb=pmcb->next;
  }
#endif
}
  
//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::ReceiveMonteCarloData(int source)
//  \brief initialize monte carlo data from another process

void MonteCarlo::ReceiveMonteCarloData(int source) {

#ifdef MPI_PARALLEL
  MonteCarloBlock *pmcb=pblock;
  Real *recv_buf;
  int size = 3; //tgas,rho,vol
  if (lorentz_transform) size+=3;
  size *= (pmcb->nx1*pmcb->nx2*pmcb->nx3); // all blocks have same size
  size += pmcb->nx1+1; size += pmcb->nx2+1; size+= pmcb->nx3+1;
  recv_buf = new Real[size];
  MPI_Request recv_rq;
  unsigned int tag = 1;
  while (pmcb != NULL) {
    MPI_Irecv(recv_buf,size,MPI_ATHENA_REAL,source,tag++,MPI_COMM_WORLD,&recv_rq);
    MPI_Wait(&recv_rq, MPI_STATUS_IGNORE);
    int p=0;
    BufferUtility::Unpack3DData(recv_buf, pmcb->rho, pmcb->is, pmcb->ie, pmcb->js, 
                                pmcb->je, pmcb->ks, pmcb->ke, p);
    BufferUtility::Unpack3DData(recv_buf, pmcb->tgas, pmcb->is, pmcb->ie, pmcb->js, 
                                pmcb->je, pmcb->ks, pmcb->ke, p);
    if (lorentz_transform)
      BufferUtility::Unpack4DData(recv_buf, pmcb->vel, 0, 2, pmcb->is, pmcb->ie, pmcb->js, 
                                  pmcb->je, pmcb->ks, pmcb->ke, p);
    BufferUtility::Unpack3DData(recv_buf, pmcb->pcoord->vol, pmcb->is, pmcb->ie, pmcb->js, 
                                pmcb->je, pmcb->ks, pmcb->ke, p);
    for (int i=pmcb->is; i<=pmcb->ie; ++i) 
      pmcb->pcoord->x1f(i) = recv_buf[p++];
    for (int i=pmcb->js; i<=pmcb->je; ++i) 
      pmcb->pcoord->x2f(i) = recv_buf[p++];
    for (int i=pmcb->ks; i<=pmcb->ke; ++i) 
      pmcb->pcoord->x3f(i) = recv_buf[p++];
    pmcb=pmcb->next;
  }

#endif
}

//----------------------------------------------------------------------------------------
//! \fn unsigned int MonteCarlo::CreateMCMPITag(int bid)
//  \brief calculate an MPI tag for monte carlo communications

unsigned int MonteCarlo::CreateMCMPITag(int bid) {
  return bid;
}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::InitializeMonteCarloBlocks(void)
//  \brief initialize grid data in each monte carlo block

void MonteCarlo::InitializeMonteCarloBlocks(void) {

  MonteCarloBlock *pmcb = pblock;

  // Check/set function pointers
  if (InitEmission == NULL) {
    std::stringstream msg;
    msg << "### FATAL ERROR in RunStaticMonteCarlo()" << std::endl
        << "InitEmission function pointer not set." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (GetTemperature == NULL)
    GetTemperature = DefaultGetTemperature;

  if (source < 0) {
    // Initialize variables over all blocks
    GetDensity(pmcb);
    GetTemperature(pmcb);
    //(pmcb->*(pmcb->GetTemperature2))();
    if (lorentz_transform) GetVelocity(pmcb);
    InitEmission(pmcb);
    pmcb = pmcb->next;
    while (pmcb != NULL) {
      GetDensity(pmcb);
      GetTemperature(pmcb);
      if (lorentz_transform) GetVelocity(pmcb);
      InitEmission(pmcb);
      pmcb = pmcb->next;
    }
    for(int i=0; i<ndest; ++i) {
      SendMonteCarloData(dest[i]);
    }
  } else {
    // Get data from another process
    ReceiveMonteCarloData(source);
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::RunStaticMonteCarlo(void)
//  \brief start evolving photons in each monte carlo block

void MonteCarlo::RunStaticMonteCarlo(void) {
 
  InitializeMonteCarloBlocks();

  MonteCarloBlock *pmcb = pblock;
  // transfer photons over all blocks
  pmcb = pblock;
  pmcb->TransferPhotons(pmcb->nphremain);
  pmcb = pmcb->next;
  while (pmcb != NULL) {
    pmcb->TransferPhotons(pmcb->nphremain);
    pmcb = pmcb->next;
  }

  return;
}

// constructor
MCCoord::MCCoord(Coordinates *pcoord, MonteCarloBlock *pmcb) {

  x1f.InitWithShallowCopy(pcoord->x1f);
  x2f.InitWithShallowCopy(pcoord->x2f);
  x3f.InitWithShallowCopy(pcoord->x3f);

  // Allocate volume array
  int ncells1 = pmcb->nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmcb->nx2 > 1) ncells2 = pmcb->nx2 + 2*(NGHOST);
  if (pmcb->nx3 > 1) ncells3 = pmcb->nx3 + 2*(NGHOST);
  vol.NewAthenaArray(ncells3,ncells2,ncells1);
  // Initialize volume array
  for (int k=pmcb->ks; k<=pmcb->ke; ++k) {
    for (int j=pmcb->js; j<=pmcb->je; ++j) {
      for (int i=pmcb->is; i<=pmcb->ie; ++i) {
        vol(k,j,i) = pcoord->GetCellVolume(k,j,i);
      }}}
}

// constructor
MCCoord::MCCoord(int ncells1, int ncells2, int ncells3) {

  x1f.NewAthenaArray(ncells1+1);
  x2f.NewAthenaArray(ncells2+1);
  x3f.NewAthenaArray(ncells3+1);

  vol.NewAthenaArray(ncells3,ncells2,ncells1);

}


// destructor
MCCoord::~MCCoord() {

  x1f.DeleteAthenaArray();
  x2f.DeleteAthenaArray();
  x3f.DeleteAthenaArray();
  vol.DeleteAthenaArray();
}

// constructor

MCRandom::MCRandom(int iseed) {
  dev = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(dev, iseed);
}

// destructor

MCRandom::~MCRandom() {

}

Real MCRandom::uniform() {

  return static_cast<Real>(gsl_rng_uniform(dev));
}
