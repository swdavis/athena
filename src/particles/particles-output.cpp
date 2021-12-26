//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2021 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles-output.cpp
//! \brief implements functions for class ParticlesOutput and its derived classes.

// C/C++ Standard Libraries
#include <fstream>   // ofstream
#include <iomanip>   // setprecision(), setw()
#include <iostream>  // <<, endl, scientific, showpoint
#include <limits>    // numeric_limits<T>
#include <sstream>   // ostringstream
#include <vector>    // vector<T>

// Athena++ headers
#include "../defs.hpp"       // ATHENA_ERROR()
#include "../mesh/mesh.hpp"  // MeshBlock
#include "particles.hpp"     // Particles
#include "particles-output.hpp"

//--------------------------------------------------------------------------------------
//! \fn void ParticlesOutput::SetNextOutput(ParameterInput* pin)
//! \brief advances the timer for next output.

void ParticlesOutput::SetNextOutput(ParameterInput* pin) {
  op.next_time += op.dt;
  pin->SetReal(op.block_name, "next_time", op.next_time);
  pin->SetInteger(op.block_name, "file_number", ++op.file_number);
}

//--------------------------------------------------------------------------------------
//! \fn void POutFormattedTable::WriteOutputFile(Mesh *pm)
//! \brief outputs the particle data in tabulated format.

void POutFormattedTable::WriteOutputFile(const Mesh *pm) {
  const int iprec(std::numeric_limits<int>::digits10);
  const int rprec(std::numeric_limits<Real>::max_digits10);
  const int wi(iprec+2);
  const int wr(rprec+8);
  std::ostringstream fname, msg;
  std::ofstream os;

  // Loop over MeshBlocks
  for (int b = 0; b < pm->nblocal; ++b) {
    const MeshBlock *pmb(pm->my_blocks(b));
    const Particles *ppar(pmb->ppar);

    // Create the filename.
    fname << op.file_basename << ".block" << pmb->gid << ".pout."
          << std::setw(5) << std::right << std::setfill('0') << op.file_number << ".tab";

    // Open the file for write.
    os.open(fname.str().data());
    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [POutFormattedTable::WriteOutputFile]"
          << std::endl << "Output file '" << fname.str() << "' could not be opened"
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

    // Close the file and get the next meshblock.
    os.close();
    fname.str("");
  }
}

