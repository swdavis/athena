//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file default_pgen.cpp
//  \brief Provides default (empty) versions of all functions in problem generator files
//  This means user does not have to implement these functions if they are not needed.
//
// The attribute "weak" is used to ensure the loader selects the user-defined version of
// functions rather than the default version given here.

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../particle/particle.hpp"

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Should be used to set initial conditions.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  if (MAGNETIC_FIELDS_ENABLED){
    for (int k=ks; k<=ke+1; ++k){
      for (int j=js; j<=je+1; ++j){
        for (int i=is; i<=ie+1; ++i){
          pfield->b.x1f(k,j,i) = 0.0;
          pfield->b.x2f(k,j,i) = 0.0;
          pfield->b.x3f(k,j,i) = 1.0;
        }
      }
    }
  }

   return;
}
