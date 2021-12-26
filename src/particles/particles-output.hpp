#ifndef PARTICLES_PARTICLES_OUTPUT_HPP_
#define PARTICLES_PARTICLES_OUTPUT_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2021 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particles-output.hpp
//! \brief defines class ParticlesOutput
//======================================================================================

// C/C++ Standard Libraries
#include <string>  // string
#include <vector>  // vector<T>

// Athena headers
#include "../mesh/mesh.hpp"                 // Mesh
#include "../outputs/output_parameters.hpp" // OutputParameters
#include "../parameter_input.hpp"           // ParameterInput
#include "particles.hpp"                    // Particles

//--------------------------------------------------------------------------------------
//! \class ParticlesOutput
//! \brief defines abstract class for output of particles data.

class ParticlesOutput {
 public:
  // Constructors
  explicit ParticlesOutput(const OutputParameters& op_in) :
      op(op_in), nint(Particles::GetNInt()), nreal(Particles::GetNReal()),
      ipname(Particles::GetIntNames()), rpname(Particles::GetRealNames()) {}
  virtual ~ParticlesOutput() {}

  // Accessors
  std::string ComposeFileName() const;
  std::string ComposeFileName(int block_id) const;
  bool CheckTimer(const Mesh *pm) const;

  // Mutators
  void SetNextOutput(ParameterInput* pin);

  // Virtual function to write output file.
  virtual void WriteOutputFile(const Mesh *pm) = 0;

 protected:
  // Instance Variables
  OutputParameters op;   //!> output parameters
  const int nint, nreal; //!> numbers of int and Real properties for each particle
  const std::vector<std::string>
      &ipname, &rpname;  //!> names of the properties
};

//--------------------------------------------------------------------------------------
//! \fn bool ParticlesOutput::CheckTimer(const Mesh *pm) const
//! \brief checks if it is time for next output.

inline bool ParticlesOutput::CheckTimer(const Mesh *pm) const {
  return (pm->time == pm->start_time) ||
         (pm->time >= op.next_time) ||
         (pm->time >= pm->tlim);
}

//--------------------------------------------------------------------------------------
//! \class POutBinaries
//! \brief defines derived class for writing a binary file.

class POutBinaries : public ParticlesOutput {
 public:
  // Constructors
  explicit POutBinaries(const OutputParameters& op) : ParticlesOutput(op) {}

  // Function to write raw data to a binary file.
  void WriteOutputFile(const Mesh *pm);
};

//--------------------------------------------------------------------------------------
//! \class POutFormattedTable
//! \brief defines derived class for writing an ASCII file.

class POutFormattedTable : public ParticlesOutput {
 public:
  // Constructors
  explicit POutFormattedTable(const OutputParameters& op) : ParticlesOutput(op) {}

  // Function to write a formatted table to an ASCII file.
  void WriteOutputFile(const Mesh *pm);
};

#endif  // PARTICLES_PARTICLES_OUTPUT_HPP_
