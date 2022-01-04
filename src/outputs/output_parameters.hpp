#ifndef OUTPUTS_OUTPUT_PARAMETERS_HPP_
#define OUTPUTS_OUTPUT_PARAMETERS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file output_parameters.hpp
//! \brief provides the structure OutputParameters.

// C/C++ headers
#include <string>

//----------------------------------------------------------------------------------------
//! \struct OutputParameters
//! \brief  container for parameters read from `<output>` block in the input file

struct OutputParameters {
  int block_number;
  std::string block_name;
  std::string file_basename;
  std::string file_id;
  std::string variable;
  std::string file_type;
  std::string data_format;
  Real next_time, dt;
  int dcycle;
  int file_number;
  bool output_slicex1, output_slicex2, output_slicex3;
  bool output_sumx1, output_sumx2, output_sumx3;
  bool include_ghost_zones, cartesian_vector;
  bool orbital_system_output;
  int islice, jslice, kslice;
  Real x1_slice, x2_slice, x3_slice;
  // TODO(felker): some of the parameters in this class are not initialized in constructor
  OutputParameters() : block_number(0), next_time(0.0), dt(0.0), file_number(0),
                       output_slicex1(false),output_slicex2(false),output_slicex3(false),
                       output_sumx1(false), output_sumx2(false), output_sumx3(false),
                       include_ghost_zones(false), cartesian_vector(false),
                       islice(0), jslice(0), kslice(0) {}
};

#endif // OUTPUTS_OUTPUT_PARAMETERS_HPP_
