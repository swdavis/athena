#ifndef MCCOORD_HPP
#define MCCOORD_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mccoord.hpp
//  \brief definitions for MCCoord class

// Athena++ classes headers
#include "../athena.hpp"
#include "montecarlo.hpp"

//----------------------------------------------------------------------------------------
//! \class MCCoord
//  \brief monte carlo specific coordinate values

class MCCoord {
public:
  MCCoord(Coordinates *pcoord, MonteCarloBlock *pmcb);
  MCCoord(int ncells1, int ncells2, int ncells3, bool acc);
  ~MCCoord();

  bool acceleration;

  AthenaArray<Real> x1f, x2f, x3f; // face  positions
  AthenaArray<Real> vol;
  AthenaArray<Real> dmin;

  Real GetMass() const {return bh_mass_;}
  Real GetSpin() const {return bh_spin_;}

private:
 // GR-specific variables
  Real bh_mass_;
  Real bh_spin_;

};



#endif // MCCOORD_HPP
