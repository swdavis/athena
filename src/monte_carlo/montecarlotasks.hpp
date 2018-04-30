#ifndef MONTECARLO_TASKS_HPP
#define MONTECARLO_TASKS_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo_tasks.hpp
//  \brief definitions for MonteCarloTasks class

// Athena++ classes headers
#include "../athena.hpp"
//#include "../athena_arrays.hpp"

class Mesh;
class ParameterInput;
class MonteCarlo;

// Current design focusses on implementing simple static post-processing so this class
// implementation will evolve.

//! \class MonteCarloTasks
//  \brief monte carlo task list

class MonteCarloTasks {
public:
  MonteCarloTasks(Mesh *pm, ParameterInput *pin);
  ~MonteCarloTasks();

};


#endif // MONTECARLO_TASKS_HPP
