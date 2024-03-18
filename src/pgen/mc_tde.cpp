//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file from_array.cpp
//! \brief Problem generator for initializing with preexisting array from HDF5 input

// C headers

// C++ headers
#include <algorithm>  // max()
#include <string>     // c_str(), string

// Athena++ headers
#include "../athena.hpp"              // Real
#include "../athena_arrays.hpp"       // AthenaArray
#include "../field/field.hpp"         // Field
#include "../globals.hpp"             // Globals
#include "../hydro/hydro.hpp"         // Hydro
#include "../eos/eos.hpp"                  // EquationOfState
#include "../inputs/hdf5_reader.hpp"  // HDF5ReadRealArray()
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"     // ParameterInput
#include "../monte_carlo/montecarlo.hpp"
#include "../monte_carlo/photon.hpp"

namespace {
  // Global variables
  bool tnorm;
  Real logemin, logemax;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief monte carlo test problem generator
//! Inputs:
//! - pin: parameters
//! Outputs: (none)
//! Notes:
//! - reads in array using a slightly modified verions of from_array.cpp
//!   - NHYDRO
//!   - total number of MeshBlocks
//!   - MeshBlock/nx3
//!   - MeshBlock/nx2
//!   - MeshBlock/nx1

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // Determine locations of initial values
  std::string input_filename = pin->GetString("problem", "input_filename");
  std::string dataset = pin->GetString("problem", "dataset");
  bool conserved = pin->GetBoolean("problem", "conserved");
  int index_dens = pin->GetInteger("problem", "index_dens");
  int index_mom1 = pin->GetInteger("problem", "index_mom1");
  int index_mom2 = pin->GetInteger("problem", "index_mom2");
  int index_mom3 = pin->GetInteger("problem", "index_mom3");
  int index_etot = pin->GetInteger("problem", "index_etot");
  std::string dataset_b1, dataset_b2, dataset_b3;
  if (MAGNETIC_FIELDS_ENABLED) {
    dataset_b1 = pin->GetString("problem", "dataset_b1");
    dataset_b2 = pin->GetString("problem", "dataset_b2");
    dataset_b3 = pin->GetString("problem", "dataset_b3");
  }
  // Set conserved array selections
  int start_cons_file[5];
  start_cons_file[1] = gid;
  start_cons_file[2] = 0;
  start_cons_file[3] = 0;
  start_cons_file[4] = 0;
  int start_cons_indices[5];
  start_cons_indices[IDN] = index_dens;
  start_cons_indices[IM1] = index_mom1;
  start_cons_indices[IM2] = index_mom2;
  start_cons_indices[IM3] = index_mom3;
  start_cons_indices[IEN] = index_etot;
  int count_cons_file[5];
  count_cons_file[0] = 1;
  count_cons_file[1] = 1;
  count_cons_file[2] = block_size.nx3;
  count_cons_file[3] = block_size.nx2;
  count_cons_file[4] = block_size.nx1;
  int start_cons_mem[4];
  start_cons_mem[1] = ks;
  start_cons_mem[2] = js;
  start_cons_mem[3] = is;
  int count_cons_mem[4];
  count_cons_mem[0] = 1;
  count_cons_mem[1] = block_size.nx3;
  count_cons_mem[2] = block_size.nx2;
  count_cons_mem[3] = block_size.nx1;

  // Set conserved values from file
  for (int n = 0; n < NHYDRO; ++n) {
    start_cons_file[0] = start_cons_indices[n];
    start_cons_mem[0] = n;
    //std::cout << gid << " " << input_filename << std::endl;
  
    if (conserved)
      HDF5ReadRealArray(input_filename.c_str(), dataset.c_str(), 5, start_cons_file,
                        count_cons_file, 4, start_cons_mem,
                        count_cons_mem, phydro->u, true);
    else
      HDF5ReadRealArray(input_filename.c_str(), dataset.c_str(), 5, start_cons_file,
                        count_cons_file, 4, start_cons_mem,
                        count_cons_mem, phydro->w, true);
  }

  // Set field array selections
  int start_field_file[4];
  start_field_file[0] = gid;
  start_field_file[1] = 0;
  start_field_file[2] = 0;
  start_field_file[3] = 0;
  int count_field_file[4];
  count_field_file[0] = 1;
  int start_field_mem[3];
  start_field_mem[0] = ks;
  start_field_mem[1] = js;
  start_field_mem[2] = is;
  int count_field_mem[3];

  // Set magnetic field values from file
  if (MAGNETIC_FIELDS_ENABLED) {
    // Set B1
    count_field_file[1] = block_size.nx3;
    count_field_file[2] = block_size.nx2;
    count_field_file[3] = block_size.nx1 + 1;
    count_field_mem[0] = block_size.nx3;
    count_field_mem[1] = block_size.nx2;
    count_field_mem[2] = block_size.nx1 + 1;
    HDF5ReadRealArray(input_filename.c_str(), dataset_b1.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, pfield->b.x1f, true);

    // Set B2
    count_field_file[1] = block_size.nx3;
    count_field_file[2] = block_size.nx2 + 1;
    count_field_file[3] = block_size.nx1;
    count_field_mem[0] = block_size.nx3;
    count_field_mem[1] = block_size.nx2 + 1;
    count_field_mem[2] = block_size.nx1;
    HDF5ReadRealArray(input_filename.c_str(), dataset_b2.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, pfield->b.x2f, true);

    // Set B3
    count_field_file[1] = block_size.nx3 + 1;
    count_field_file[2] = block_size.nx2;
    count_field_file[3] = block_size.nx1;
    count_field_mem[0] = block_size.nx3 + 1;
    count_field_mem[1] = block_size.nx2;
    count_field_mem[2] = block_size.nx1;
    HDF5ReadRealArray(input_filename.c_str(), dataset_b3.c_str(), 4, start_field_file,
                      count_field_file, 3, start_field_mem,
                      count_field_mem, pfield->b.x3f, true);
  }

  // Make no-op collective reads if using MPI and ranks have unequal numbers of blocks
#ifdef MPI_PARALLEL
  {
    int num_blocks_this_rank = pmy_mesh->nblist[Globals::my_rank];
    if (lid == num_blocks_this_rank - 1) {
      int block_shortage_this_rank = 0;
      for (int rank = 0; rank < Globals::nranks; ++rank) {
        block_shortage_this_rank =
            std::max(block_shortage_this_rank,
                     pmy_mesh->nblist[rank] - num_blocks_this_rank);
      }
      for (int block = 0; block < block_shortage_this_rank; ++block) {
        for (int n = 0; n < NHYDRO; ++n) {
          start_cons_file[0] = start_cons_indices[n];
          start_cons_mem[0] = n;
          if (conserved)
            HDF5ReadRealArray(input_filename.c_str(), dataset.c_str(), 5,
                              start_cons_file, count_cons_file, 4,
                              start_cons_mem, count_cons_mem,
                              phydro->u, true, true);
          else
            HDF5ReadRealArray(input_filename.c_str(), dataset.c_str(), 5,
                              start_cons_file, count_cons_file, 4,
                              start_cons_mem, count_cons_mem,
                              phydro->w, true, true);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          count_field_file[1] = block_size.nx3;
          count_field_file[2] = block_size.nx2;
          count_field_file[3] = block_size.nx1 + 1;
          count_field_mem[0] = block_size.nx3;
          count_field_mem[1] = block_size.nx2;
          count_field_mem[2] = block_size.nx1 + 1;
          HDF5ReadRealArray(input_filename.c_str(), dataset_b1.c_str(), 4,
                            start_field_file, count_field_file, 3,
                            start_field_mem, count_field_mem,
                            pfield->b.x1f, true, true);
          count_field_file[1] = block_size.nx3;
          count_field_file[2] = block_size.nx2 + 1;
          count_field_file[3] = block_size.nx1;
          count_field_mem[0] = block_size.nx3;
          count_field_mem[1] = block_size.nx2 + 1;
          count_field_mem[2] = block_size.nx1;
          HDF5ReadRealArray(input_filename.c_str(), dataset_b2.c_str(), 4,
                            start_field_file, count_field_file, 3,
                            start_field_mem, count_field_mem,
                            pfield->b.x2f, true, true);
          count_field_file[1] = block_size.nx3 + 1;
          count_field_file[2] = block_size.nx2;
          count_field_file[3] = block_size.nx1;
          count_field_mem[0] = block_size.nx3 + 1;
          count_field_mem[1] = block_size.nx2;
          count_field_mem[2] = block_size.nx1;
          HDF5ReadRealArray(input_filename.c_str(), dataset_b3.c_str(), 4,
                            start_field_file, count_field_file, 3,
                            start_field_mem, count_field_mem,
                            pfield->b.x3f, true, true);
        }
      }
    }
  }
#endif

  if (!conserved) {
    // Set index bounds
    int il = is - NGHOST;
    int iu = ie + NGHOST;
    int jl = js;
    int ju = je;
    if (block_size.nx2 > 1) {
      jl -= NGHOST;
      ju += NGHOST;
    }
    int kl = ks;
    int ku = ke;
    if (block_size.nx3 > 1) {
      kl -= NGHOST;
      ku += NGHOST;
    }
    // Calculate cell-centered magnetic field
    if (MAGNETIC_FIELDS_ENABLED)
      pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, il, iu, jl, ju,
                                         kl, ku);
    // Initialize conserved
    peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, il, iu, jl, ju,
                               kl, ku);
  }

}


//========================================================================================
//! \fn void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin)
//! \brief Analogous to problem generator but used in support of InitializePhoton
//========================================================================================

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {

  tnorm = pin->GetOrAddBoolean("problem","tnorm",false);
  if (tnorm) {
    // interpret as xmin/xmax with x=E/(kb*T)
    Real kb = 1.380649e-16;
    logemin = log(kb*pin->GetReal("problem", "emin"));
    logemax = log(kb*pin->GetReal("problem", "emax"));
  } else {
    // Set the energy boundaries for free-free emission
    Real everg = 1.6021772e-12;
    logemin = log(everg*pin->GetReal("problem", "emin"));
    logemax = log(everg*pin->GetReal("problem", "emax"));
  }

}

//========================================================================================
//! \fn void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe)
//! \brief Initializes Photon packets before integration
//========================================================================================

void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe) {


  // Set initial cells and emission weights for all photon samples
  SetEmissionCellWeight(pphot,ips,ipe);

  for (int ip=ips; ip<=ipe; ip++) {

    // Set status flag
    pphot->statp[ip] = EVOLVING;

    // Obtain initial position within zone
    GetZonePosition(pphot,pran,pcoord,ip);

    // Set maximum integration time
    pphot->dtp[ip] = pphot->pmy_mcb->pmy_mc->tmax;

    // Obtain intitial energy, polarization, direction and weight
    // Utilize free-free emission function in emission.cpp
    if(tnorm) {
      Real logtg = log(tgas(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip]));
      PhotonEmitFreeFree(this,pphot,logemin+logtg,logemax+logtg,ip);
    } else{
      PhotonEmitFreeFree(this,pphot,logemin,logemax,ip);
    }

    if (pphot->wp[ip] < 0.0)
      pphot->statp[ip] = DESTROYED;
    else
      pphot->statp[ip] = EVOLVING;

    // initialize scattering number
    pphot->nscp[ip] = 0;

    // Initialize the absorption and scattering extinction coefficients
    // to the values appropriate in the emitted zone
    pphot->acp[ip] = AbsorptionOpacity(this,pphot,ip);
    pphot->scp[ip] = ScatteringOpacity(this,pphot,ip);

  } // loop over ip

}
