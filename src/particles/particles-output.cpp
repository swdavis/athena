//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2021 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles-output.cpp
//! \brief implements functions for class ParticlesOutput and its derived classes.

// C/C++ Standard Libraries
#include <cstring>   // size_t, memcpy()
#include <fstream>   // ofstream
#include <iomanip>   // setprecision(), setw()
#include <iostream>  // <<, endl, ios, scientific, showpoint
#include <limits>    // numeric_limits<T>
#include <sstream>   // ostringstream
#include <string>    // string
#include <vector>    // vector<T>

// Athena++ headers
#include "../defs.hpp"       // ATHENA_ERROR()
#include "../globals.hpp"    // my_rank, nranks
#include "../mesh/mesh.hpp"  // MeshBlock
#include "particles.hpp"     // Particles
#include "particles-output.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Constants
static const std::size_t SIZE_OF_INT(sizeof(int));
static const std::size_t SIZE_OF_REAL(sizeof(Real));

// Local functions
static char* add_data(char* pbuf, int n);
static char* add_data(char* pbuf, Real n);
static char* add_data(char* pbuf, const std::string& s);

//--------------------------------------------------------------------------------------
//! \fn std::string ParticlesOutput::ComposeFileName() const
//! \brief composes the output file name without the extension.

std::string ParticlesOutput::ComposeFileName() const {
  std::ostringstream fname;
  fname << op.file_basename << ".pout."
        << std::setw(5) << std::right << std::setfill('0') << op.file_number;
  return fname.str();
}

//--------------------------------------------------------------------------------------
//! \fn std::string ParticlesOutput::ComposeFileName(int block_id) const
//! \brief composes the output file name without the extension, given block_id.

std::string ParticlesOutput::ComposeFileName(int block_id) const {
  std::ostringstream fname;
  fname << op.file_basename << ".block" << block_id << ".pout."
        << std::setw(5) << std::right << std::setfill('0') << op.file_number;
  return fname.str();
}

//--------------------------------------------------------------------------------------
//! \fn void ParticlesOutput::SetNextOutput(ParameterInput* pin)
//! \brief advances the timer for next output.

void ParticlesOutput::SetNextOutput(ParameterInput* pin) {
  op.next_time += op.dt;
  pin->SetReal(op.block_name, "next_time", op.next_time);
  pin->SetInteger(op.block_name, "file_number", ++op.file_number);
}

//--------------------------------------------------------------------------------------
//! \fn POutBinaries::POutBinaries(const OutputParameters& op)
//! \brief is a constructor.

POutBinaries::POutBinaries(const OutputParameters& op)
: ParticlesOutput(op) {
  // Compute the data size for the header.
  header_size = SIZE_OF_REAL + 4 * SIZE_OF_INT;
  for (int j = 0; j < nint; ++j)
    header_size += ipname[j].size() + 1;
  for (int j = 0; j < nreal; ++j)
    header_size += rpname[j].size() + 1;

  // Compute the data size for each particle.
  psize = nint * SIZE_OF_INT + nreal * SIZE_OF_REAL;
}

//--------------------------------------------------------------------------------------
//! \fn void POutBinaries::WriteOutputFile(Mesh *pm)
//! \brief outputs the particle data in raw binaries.

void POutBinaries::WriteOutputFile(const Mesh *pm) {
  // Count total number of particles in each process.
  int my_npar(0), npar_tot(0);
  for (int b = 0; b < pm->nblocal; ++b)
    my_npar += pm->my_blocks(b)->ppar->GetNPar();
#ifdef MPI_PARALLEL
  int *npar_in_rank(new int[Globals::nranks]);
  MPI_Allgather(&my_npar, 1, MPI_INT, npar_in_rank, 1, MPI_INT, MPI_COMM_WORLD);
  for (int i = 0; i < Globals::nranks; ++i)
    npar_tot += npar_in_rank[i];
#else // MPI_PARALLEL
  npar_tot = my_npar;
#endif // MPI_PARALLEL

  // Create the output file.
  bool flag;
  const std::string fname(ComposeFileName() + ".dat");
#ifdef MPI_PARALLEL
  MPI_File fhandle;
  flag = MPI_File_open(MPI_COMM_WORLD, fname.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY,
      MPI_INFO_NULL, &fhandle) != MPI_SUCCESS;
#else // MPI_PARALLEL
  std::ofstream os(fname, std::ios::out|std::ios::binary);
  flag = !os.is_open();
#endif // MPI_PARALLEL
  if (flag) {
    std::ostringstream msg;
    msg << "### FATAL ERROR in function [POutBinaries::WriteOutputFile]\n"
        << "Output file '" << fname << "' could not be opened.\n";
    ATHENA_ERROR(msg);
  }

  // Write the header.
  char *buf, *pbuf;
  if (Globals::my_rank == 0) {
    pbuf = buf = new char[header_size];
    pbuf = add_data(pbuf, static_cast<int>(SIZE_OF_REAL));
    pbuf = add_data(pbuf, pm->time);
    pbuf = add_data(pbuf, nint);
    pbuf = add_data(pbuf, nreal);
    for (int j = 0; j < nint; ++j)
      pbuf = add_data(pbuf, ipname[j]);
    for (int j = 0; j < nreal; ++j)
      pbuf = add_data(pbuf, rpname[j]);
    pbuf = add_data(pbuf, npar_tot);
#ifdef MPI_PARALLEL
    if (MPI_File_write(fhandle, buf, pbuf - buf, MPI_CHAR,
                       MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "### FATAL ERROR in function [POutBinaries::WriteOutputFile]\n"
          << "Unable to write the header.\n";
      ATHENA_ERROR(msg);
    }
#else // MPI_PARALLEL
    os.write(buf, pbuf - buf);
#endif // MPI_PARALLEL
    delete [] buf;
  }

  // Write the particle data.
  pbuf = buf = new char[my_npar * psize];
  for (int b = 0; b < pm->nblocal; ++b) {
    const Particles *ppar(pm->my_blocks(b)->ppar);
    const std::vector<int>* intprop(ppar->GetIntProps());
    const std::vector<Real>* realprop(ppar->GetRealProps());
    for (int k = 0; k < ppar->GetNPar(); ++k) {
      for (int j = 0; j < nint; ++j)
        pbuf = add_data(pbuf, intprop[j][k]);
      for (int j = 0; j < nreal; ++j)
        pbuf = add_data(pbuf, realprop[j][k]);
    }
  }
#ifdef MPI_PARALLEL
  int npar_prev(0);
  for (int i = 0; i < Globals::my_rank; ++i)
    npar_prev += npar_in_rank[i];
  MPI_Offset offset(header_size + npar_prev * psize);
  if (MPI_File_write_at_all(fhandle, offset, buf, pbuf - buf, MPI_CHAR,
                     MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    std::ostringstream msg;
    msg << "### FATAL ERROR in function [POutBinaries::WriteOutputFile]\n"
        << "Unable to write the particle data.\n";
    ATHENA_ERROR(msg);
  }
#else // MPI_PARALLEL
  os.write(buf, pbuf - buf);
#endif // MPI_PARALLEL
  delete [] buf;

  // Close the file.
#ifdef MPI_PARALLEL
  MPI_File_close(&fhandle);
#else // MPI_PARALLEL
  os.close();
#endif // MPI_PARALLEL

#ifdef MPI_PARALLEL
  delete [] npar_in_rank;
#endif // MPI_PARALLEL
}

//--------------------------------------------------------------------------------------
//! \fn void POutFormattedTable::WriteOutputFile(Mesh *pm)
//! \brief outputs the particle data in tabulated format.

void POutFormattedTable::WriteOutputFile(const Mesh *pm) {
  const int iprec(std::numeric_limits<int>::digits10);
  const int rprec(std::numeric_limits<Real>::max_digits10);
  const int wi(iprec+2);
  const int wr(rprec+8);
  std::ostringstream msg;
  std::ofstream os;

  // Loop over MeshBlocks
  for (int b = 0; b < pm->nblocal; ++b) {
    const MeshBlock *pmb(pm->my_blocks(b));
    const Particles *ppar(pmb->ppar);

    // Create the filename.
    std::string fname(ComposeFileName(pmb->gid) + ".tab");

    // Open the file for write.
    os.open(fname);
    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [POutFormattedTable::WriteOutputFile]"
          << std::endl << "Output file '" << fname << "' could not be opened"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    // Write the time and numbers of properties.
    os << std::scientific << std::showpoint << std::setprecision(rprec);
    os << "# Athena++ particle data at time = " << pm->time
       << " with nint = " << nint << " and nreal = " << nreal << std::endl;

    // Write the column head.
    os << '#';
    for (int j = 0; j < nint; ++j)
       os << std::setw(wi) << ipname[j];
    for (int j = 0; j < nreal; ++j)
       os << std::setw(wr) << rpname[j];
    os << '\n';

    // Write the particle data in the meshblock.
    const std::vector<int>* intprop(ppar->GetIntProps());
    const std::vector<Real>* realprop(ppar->GetRealProps());
    for (int k = 0; k < ppar->GetNPar(); ++k) {
      os << ' ';
      for (int j = 0; j < nint; ++j)
        os << std::setw(wi) << intprop[j][k];
      for (int j = 0; j < nreal; ++j)
        os << std::setw(wr) << realprop[j][k];
      os << '\n';
    }

    // Close the file.
    os.close();
  }
}

//--------------------------------------------------------------------------------------
//! \fn char* add_data(char* pbuf, int n)
//! \brief writes an integer to the buffer and advance the pointer.

inline char* add_data(char* pbuf, int n) {
  std::memcpy(pbuf, &n, SIZE_OF_INT);
  return pbuf + SIZE_OF_INT;
}

//--------------------------------------------------------------------------------------
//! \fn char* add_data(char* pbuf, Real x)
//! \brief writes a real to the buffer and advance the pointer.

inline char* add_data(char* pbuf, Real x) {
  std::memcpy(pbuf, &x, SIZE_OF_REAL);
  return pbuf + SIZE_OF_REAL;
}

//--------------------------------------------------------------------------------------
//! \fn char* add_data(char* pbuf, const std::string& s)
//! \brief writes a string to the buffer and advance the pointer.

inline char* add_data(char* pbuf, const std::string& s) {
  const std::size_t size(s.size() + 1);
  std::memcpy(pbuf, s.c_str(), size);
  return pbuf + size;
}
