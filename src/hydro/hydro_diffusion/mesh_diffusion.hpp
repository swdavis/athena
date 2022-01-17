#ifndef HYDRO_HYDRO_DIFFUSION_MESH_DIFFUSION_HPP_
#define HYDRO_HYDRO_DIFFUSION_MESH_DIFFUSION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2022 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh_diffusion.hpp
//! \brief defines class MeshDiffusion

// C/C++ headers

// Athena++ headers
#include "../../athena_arrays.hpp"   // AthenaArray
#include "../../mesh/mesh.hpp"       // MeshBlock
#include "../../parameter_input.hpp" // ParameterInput

//----------------------------------------------------------------------------------------
//! \class MeshDiffusion
//! \brief data and functions for mesh diffusion

class MeshDiffusion {
 public:
  // Constructors
  MeshDiffusion(MeshBlock* pmb, ParameterInput* pin);

  // Diffusion operations
  void AddFluxes(const AthenaArray<Real>& cons, AthenaArray<Real>* flux) const;

 private:
  // Instance Variables
  MeshBlock* pmb; // pointer to my mesh-block
  Real nu2mesh;   // 4th-order diffusion coefficient

  // Diffusion operations
  void AddFluxHyper2(const AthenaArray<Real> &cons, AthenaArray<Real> *flux) const;
};

#endif // HYDRO_HYDRO_DIFFUSION_MESH_DIFFUSION_HPP_
