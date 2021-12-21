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
#include "../outputs/output_parameters.hpp"
#include "particles.hpp"

//--------------------------------------------------------------------------------------
//! \class ParticlesOutput
//! \brief defines the class for output of particles data.

class ParticlesOutput {
 public:
  // Constructors
  explicit ParticlesOutput(const OutputParameters& op_in) : op(op_in) {}
  ~ParticlesOutput() {}

 protected:
  // Instance Variables
  OutputParameters op;
};

#endif  // PARTICLES_PARTICLES_OUTPUT_HPP_
