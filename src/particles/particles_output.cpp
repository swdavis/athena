//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2021 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles_output.cpp
//! \brief implements functions for class ParticlesOutput and its derived classes.

// C/C++ Standard Libraries
#include <fstream>   // ofstream
#include <iomanip>   // setprecision(), setw()
#include <iostream>  // <<, endl, scientific, showpoint
#include <limits>    // numeric_limits<T>
#include <sstream>   // ostringstream

// Athena++ headers
#include "../defs.hpp"       // ATHENA_ERROR()
#include "../mesh/mesh.hpp"  // MeshBlock
#include "particles.hpp"     // Particles
#include "particles_output.hpp"

//--------------------------------------------------------------------------------------
//! \fn void POutFormattedTable::WriteOutputFile(Mesh *pm)
//! \brief outputs the particle data in tabulated format.

void POutFormattedTable::WriteOutputFile(const Mesh *pm) {
  const int nint(Particles::nint), nreal(Particles::nreal);
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
    fname << op.file_basename
          << ".block" << pmb->gid << '.' << op.file_id
          << '.' << std::setw(5) << std::right << std::setfill('0') << op.file_number
          << '.' << "par.tab";

    // Open the file for write.
    os.open(fname.str().data());
    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [POutFormattedTable::WriteOutputFile]"
          << std::endl << "Output file '" << fname.str() << "' could not be opened"
          << std::endl;
      ATHENA_ERROR(msg);
    }

    // Write the time.
    os << std::scientific << std::showpoint << std::setprecision(rprec);
    os << "# Athena++ particle data at time = " << pm->time << std::endl;

    // Write the column head.
    os << '#';
    for (int j = 0; j < nint; ++j)
       os << std::setw(wi) << Particles::ipname[j];
    for (int j = 0; j < nreal; ++j)
       os << std::setw(wr) << Particles::rpname[j];
    os << '\n';

    // Write the particle data in the meshblock.
    for (int k = 0; k < ppar->npar; ++k) {
      os << ' ';
      for (int j = 0; j < nint; ++j)
        os << std::setw(wi) << ppar->intprop[j][k];
      for (int j = 0; j < nreal; ++j)
        os << std::setw(wr) << ppar->rp[j][k];
      os << '\n';
    }

    // Close the file and get the next meshblock.
    os.close();
    fname.str("");
  }
}

