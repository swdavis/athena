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
#include "../scalars/scalars.hpp"

// GSL library
#if RAN3 == 0
#include <gsl/gsl_randist.h>
#endif

//----------------------------------------------------------------------------------------
//! MonteCarlo constructor, builds monte carlo using parameters in input file

MonteCarlo::MonteCarlo(ParameterInput *pin, Mesh *pmesh) {

  std::stringstream msg;

  pmy_mesh = pmesh;

  UserWorkInMove=nullptr;
  GetEmission=nullptr;
  GetTemperature=nullptr;

  // Set flags that control emission, absorption and scattering
  InitializeEmissionFlags(pin);

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

  // SWD: replace dynamic/coupled with single method flag?
  dynamic = pin->GetOrAddBoolean("montecarlo","dynamic",false);
  coupled = pin->GetOrAddBoolean("montecarlo","coupled",false);
  boosts = pin->GetOrAddBoolean("montecarlo","boosts",false);
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

  // Set mininmum weight if using weighting for absorption
  weightratio = pin->GetOrAddReal("montecarlo","minweight",1.0e-20);

  // Set photon integration totals
  nsamp = pin->GetInteger("montecarlo","nphot");
  nout = pin->GetOrAddInteger("montecarlo","nout",1);

  // Set default size parameters
  max_phots_init = pin->GetOrAddInteger("montecarlo","max_phots_init",10000);
  list_size_init = pin->GetOrAddInteger("montecarlo","list_size_init",10000);
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
  for (int i=0; i<nblocal; i++)
    delete my_blocks(i);
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
//! \fn void MonteCarlo::EnrollUserEmissionFunction(EmisFunc_t emissfunc)
//! \brief Enroll a user-defined function for computing emission array

void MonteCarlo::EnrollUserEmissionFunction(EmisFunc_t emissfunc) {

  GetEmission = emissfunc;
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
        pmcb->rho(k,j,i) = pmcb->rho_cgs * pmcb->pmy_block->phydro->u(IDN,k,j,i);
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
        pmcb->vel(0,k,j,i) = pmcb->vel_cgs *
          pmcb->pmy_block->phydro->u(IM1,k,j,i) / rho;
        pmcb->vel(1,k,j,i) = pmcb->vel_cgs *
          pmcb->pmy_block->phydro->u(IM2,k,j,i) / rho;
        pmcb->vel(2,k,j,i) = pmcb->vel_cgs *
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

  Real tconv;
  if (pmcb->tgas_cgs <= 0.)
    tconv = 1. / rideal;
  else
    tconv = pmcb->tgas_cgs;

  // compute temperature from pressure and density
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        Real tgas = tconv * phydro->w(IEN,k,j,i) / phydro->w(IDN,k,j,i);
        // apply temperature floor
        pmcb->tgas(k,j,i) = (tgas > pmcb->tfloor_cgs) ? tgas : pmcb->tfloor_cgs;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::Initialize(ParameterInput *pinput)
//! \brief initialize grid data in each monte carlo block

void MonteCarlo::Initialize(ParameterInput *pin) {

  if (dynamic) {
    tmax = pin->GetOrAddReal("montecarlo","tmax",-1.);
    if (tmax < 0.)
      tmax = pmy_mesh->dt;
    tint = pmy_mesh->dt;
  } else {
    // initialize timing parameters if static calculation
    tint = pin->GetOrAddReal("montecarlo","tint",1.);
    tmax = pin->GetOrAddReal("montecarlo","tmax",HUGE_NUMBER);
  }
  // convert to cgs units
  Real vel_cgs = pin->GetOrAddReal("problem","vel_cgs",1.);
  Real l_cgs = pin->GetOrAddReal("problem","l_cgs",1.);
  time_cgs = pin->GetOrAddReal("problem","time_cgs",l_cgs/vel_cgs);
  tint *= time_cgs;
  tmax *= time_cgs;

  if (GetTemperature == nullptr)
    GetTemperature = DefaultGetTemperature;

  // Initialize monte carlo blocks
  for (int i=0; i<nblocal; i++) {
    MonteCarloBlock *pmcb = my_blocks(i);
    // Initialize variables over all blocks
    GetDensity(pmcb);
    GetTemperature(pmcb);
    if (boosts) GetVelocity(pmcb);
    if (NSCALARS > 0) GetScalars(pmcb);

    // initialize counters to zero
    pmcb->nscat = pmcb->nesc = pmcb->nabs = pmcb->ndes = 0;
    pmcb->loop_max_size = pin->GetOrAddInteger("montecarlo","loop_max_size",1000);

    // Call problem generators for Monte Carlo
    pmcb->MonteCarloProblemGenerator(pin);
  }

  // Initialize emission arrays, if needed
  //ComputeEmission();

}

//----------------------------------------------------------------------------------------
//! \fn void void MonteCarlo::InitializeEmissionFlags(ParameterInput *pin)
//! \brief initialize flags that control emission

void MonteCarlo::InitializeEmissionFlags(ParameterInput *pin) {
  // This function sets flags that determine how emission will be handled
  // Options currently include none, user, and free-free
  // "none" means that no arrays are set up to assist photon initialization
  // "user" means that meshblock arrays will be set up but the user will provide the
  //        emissivity function
  // "freefree" mean that the default freefree emissivity function will be used for
  //            the emission array
  emission_flag = GetEmissionFlag(pin->GetOrAddString("montecarlo","emission","none"));

  // Set emmisivity functions and flag for determining emission array
  if (emission_flag == EMISNONE) {
    GetEmission = nullptr; // left unset
    emission_array = false; // do not allocate memory for array
  } else if (emission_flag ==  EMISUSER) {
    GetEmission = nullptr; // must be set in InitUserMonteCarloData
    emission_array = true; // allocate memory for array
  } else if (emission_flag ==  EMISFF) {
    GetEmission = GetEmissionFreeFree;
    emission_array = true; // allocate memory for array
  }

  // In addition the user has a choice of how this emission array will be utilized
  // either it will determine number of photon per zone with constant initial weights
  // or it will be used to set the initial weight
  emission_eqwt = pin->GetOrAddBoolean("montecarlo","equal_weight",false);

}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::ComputeEmission()
//! \brief compute the emission on all blocks

void MonteCarlo::ComputeEmission() {

  if (emission_flag == EMISNONE) {
    // Do nothing.  nphremain needs to be set in the problem generator

    // Don't compute emission array -- just set nphremain on meshblocks
    // to all be the same value, resetting nsamp if needed    
    /*
    int nphblock = nsamp / nbtotal;
    nsamp = nphblock * nbtotal; // adjust nsamp if needed

    for (int nb=0; nb<nblocal; nb++) {
      my_blocks(nb)->nphremain = nphblock;
      my_blocks(nb)->nphrun = 0;
    }
    */

    return;
  } else if (GetEmission == nullptr) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ComputeEmission" << std::endl
        << "emission method is not none, but GetEmission is not set."
        << std::endl;
    ATHENA_ERROR(msg);
  }

  // compute emission on all blocks
  Real emm_min = SQR(HUGE_NUMBER), emm_max = -HUGE_NUMBER, emm_tot = 0.;
  Real *tot_block = new Real[nblocal];

  for (int nb=0; nb<nblocal; nb++) {
    Real min_block, max_block;
    my_blocks(nb)->ComputeEmissionArray(min_block,max_block,tot_block[nb]);
    emm_tot += tot_block[nb];
    emm_min = (emm_min < min_block) ? emm_min : min_block;
    emm_max = (emm_max > max_block) ? emm_max : max_block;
  }

  // Compute emmision properties overall processes
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,&emm_min,1,MPI_ATHENA_REAL,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&emm_max,1,MPI_ATHENA_REAL,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&emm_tot,1,MPI_ATHENA_REAL,MPI_SUM,MPI_COMM_WORLD);
#endif

  if (emission_eqwt) {
    // emmision weights are all equal
    Real ave_weight = emm_tot/static_cast<Real>(nsamp);
    int nsampnew = 0;
    for (int nb=0; nb<nblocal; nb++) {
      Real ntarget = tot_block[nb] / ave_weight;
      int nsblock = static_cast<int>(ntarget);
      Real diff = ntarget - static_cast<Real>(nsblock);
      if (my_blocks(nb)->pran->uniform() < diff)
        nsblock += 1;
      if (nsblock > 0) {
        my_blocks(nb)->emiss_to_weight = tot_block[nb] / static_cast<Real>(nsblock);
      } else {
        my_blocks(nb)->emiss_to_weight = 0.;
      }
      my_blocks(nb)->minweight = weightratio * ave_weight;
      my_blocks(nb)->nphremain = nsblock;
      my_blocks(nb)->nphrun = 0;
      nsampnew += nsblock;
    }
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE,&nsampnew,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
    if (nsamp != nsampnew) {
      if (Globals::my_rank == 0)
        std::cout << "Updating nsample to " << nsampnew << std::endl;
      nsamp = nsampnew;
    }
  } else {
    // emission weights will just be equal to emmisivity
    // Determine number of photons per block per step assuming each block is equal
    int nphblock = nsamp / nbtotal;
    if ((nphblock * nbtotal != nsamp) && (Globals::my_rank == 0))
      std::cout << "Updating nsample to " << nphblock * nbtotal << std::endl;
    nsamp = nphblock * nbtotal; // adjust nsamp if needed

    for (int nb=0; nb<nblocal; nb++) {
      my_blocks(nb)->nphremain = nphblock;
      my_blocks(nb)->nphrun = 0;
      my_blocks(nb)->minweight = weightratio * emm_min;
      my_blocks(nb)->emiss_to_weight = static_cast<Real>(ncells)/static_cast<Real>(nsamp);
    }

  }
  // Report emissivity ranges
  if (Globals::my_rank == 0) {
    std::cout << "Emission array range (min, max), total: " << emm_min << " "
              << emm_max << " " << emm_tot << std::endl;
  }


}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarlo::RunStaticMonteCarlo(Outputs *pouts, Mesh *pmesh,
//!                                          ParameterInput *pinput)
//! \brief Finish Initialization of MonteCarloBlocks and run steady-state MC calculation

void MonteCarlo::RunStaticMonteCarlo(Outputs *pouts, Mesh *pmesh,
                                     ParameterInput *pinput) {

  tint /= static_cast<Real>(nout);
  for (int i=0; i<nout; i++) {

    ComputeEmission();

    // Clear Boundary buffers for photons
    for(int nb=0; nb<nblocal; ++nb)
      my_blocks(nb)->pphot->ClearBoundary();

    // initialize monte carlo counter
    nphrun = 0;
    bool photons_remain = true; // True if photons on any process
    while(photons_remain) {

      for(int nb=0; nb<nblocal; ++nb){
        if (raytrace_flag)
          my_blocks(nb)->RayTracePhotonsOnBlock();
        else
          my_blocks(nb)->TransferPhotonsOnBlock();
      }
      photons_remain = CheckAndBroadCastPhotonsRemaining();
    }

    // Report diagnostic results from all blocks
    int64_t nesc = 0, nabs = 0, ndes = 0, nscat = 0;
    for(int nb=0; nb<nblocal; ++nb){
      MonteCarloBlock *pmcb = my_blocks(nb);
      nesc += pmcb->nesc;
      nabs += pmcb->nabs;
      ndes += pmcb->ndes;
      nscat += pmcb->nscat;
      nphrun += pmcb->nphrun;
    }
    pmcout->UpdateOutputCount(nphrun);

    int ntot = nphrun;
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE,&nesc,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&nabs,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&ndes,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&nscat,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&ntot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
    if (Globals::my_rank == 0) {
      std::cout  << "ntot, nesc, nabs, ndes, nscat/ntot: "
                 << ' ' << ntot << ' ' << nesc
                 << ' ' << nabs << ' ' << ndes << ' ';
      if (ntot > 0)
        std::cout << static_cast<Real>(nscat)/static_cast<Real>(ntot) << std::endl;
      else
        std::cout << 0. << std::endl;
    }

    Real norm_mom = static_cast<Real>(i+1);
    // normalize moments for output
    for(int nb=0; nb<nblocal; ++nb) {
      MonteCarloBlock *pmcb = my_blocks(nb);
      if (pmcb->moments_rad)
        pmcb->NormalizeMoments(true,norm_mom);
      if (pmcb->call_srcterms)
        pmcb->NormalizeSourceTerms(true,norm_mom);
    }

    // Write outputs
    pouts->MakeOutputs(pmesh,this,pinput,true);

    // unnormalize moments after output
    for(int nb=0; nb<nblocal; ++nb) {
      MonteCarloBlock *pmcb = my_blocks(nb);
      if (pmcb->moments_rad)
        pmcb->NormalizeMoments(false,norm_mom);
      if (pmcb->call_srcterms)
        pmcb->NormalizeSourceTerms(false,norm_mom);
    }

  } // end loop over nout

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
  //if (Globals::my_rank == 0)
  //  printf("here\n");
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
                                     ParameterInput *pin) {

  // initialize counter
  nphrun = 0;
  tmax = pin->GetOrAddReal("montecarlo","tmax",-1.);
  if (tmax < 0.)
    tmax = pmy_mesh->dt;
  tint = pmy_mesh->dt;


  // Reset monte carlo blocks
  for (int i=0; i<nblocal; i++) {
    MonteCarloBlock *pmcb = my_blocks(i);
    // Initialize variables over all blocks
    GetDensity(pmcb);
    GetTemperature(pmcb);
    if (boosts) GetVelocity(pmcb);
    if (NSCALARS > 0) GetScalars(pmcb); //scalars
    // reset counters
    pmcb->nscat = pmcb->nesc = pmcb->nabs = pmcb->ndes = 0;
  }

  // update emission arrays, if needed
  ComputeEmission();

  // reset moments/sourcterms for start of new timestep
  for(int nb=0; nb<nblocal; ++nb) {
    MonteCarloBlock *pmcb = my_blocks(nb);
    if (pmcb->call_moments)
      pmcb->ResetMoments();
    if (pmcb->call_srcterms)
      pmcb->ResetSourceTerms();
  }

  bool photons_remain = true; // True if photons on any process
  while(photons_remain) {

    for(int nb=0; nb<nblocal; ++nb){
      my_blocks(nb)->TransferPhotonsOnBlock();
    }
    photons_remain = CheckAndBroadCastPhotonsRemaining();
  }

  // Report diagnostic results from all blocks
  int64_t nesc = 0, nabs = 0, ndes = 0, nscat = 0;
  for(int nb=0; nb<nblocal; ++nb){
    MonteCarloBlock *pmcb = my_blocks(nb);
    nesc += pmcb->nesc;
    nabs += pmcb->nabs;
    ndes += pmcb->ndes;
    nscat += pmcb->nscat;
    nphrun += pmcb->nphrun;
  }
  pmcout->UpdateOutputCount(nphrun);

  int ntot = nphrun;
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

  for(int nb=0; nb<nblocal; ++nb) {
    MonteCarloBlock *pmcb = my_blocks(nb);
    if (pmcb->moments_rad)
      pmcb->NormalizeMoments(true,1.);
    if (pmcb->call_srcterms)
      pmcb->NormalizeSourceTerms(true,1.);

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
