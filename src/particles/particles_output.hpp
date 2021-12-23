#ifndef PARTICLES_PARTICLES_OUTPUT_HPP_
#define PARTICLES_PARTICLES_OUTPUT_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2021 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particles_output.hpp
//! \brief defines class ParticlesOutput
//======================================================================================

// C/C++ Standard Libraries

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
      op(op_in), nint(Particles::GetNInt()), nreal(Particles::GetNReal()) {}
  ~ParticlesOutput() {}

  // Mutators
  void SetNextOutput(ParameterInput* pin);

  // Virtual function to write output file.
  virtual void WriteOutputFile(const Mesh *pm) = 0;

 protected:
  // Instance Variables
  OutputParameters op;  //!> output parameters
  int nint, nreal;      //!> numbers of int and Real properties for each particle
};

//--------------------------------------------------------------------------------------
//! \class POutFormattedTable
//! \brief defines derived class for writing an ASCII file.

class POutFormattedTable : public ParticlesOutput {
 public:
  // Constructors
  explicit POutFormattedTable(const OutputParameters& op_in) : ParticlesOutput(op_in) {}

  // Function to write a formatted table to an ASCII file.
  void WriteOutputFile(const Mesh *pm);
};

#endif  // PARTICLES_PARTICLES_OUTPUT_HPP_
