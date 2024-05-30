#ifndef MCCOORD_HPP
#define MCCOORD_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mccoord.hpp
//! \brief definitions for MCCoord base class and derived classes

// Athena++ classes headers
#include "../athena.hpp"
#include "montecarlo.hpp"

#define NCOORD 4

//----------------------------------------------------------------------------------------
//! \class MCCoord
//! \brief monte carlo specific coordinate value, base class

class MCCoord {
public:
  MCCoord(Coordinates *pcoord, MonteCarloBlock *pmcb);
  MCCoord(int ncells1, int ncells2, int ncells3, bool acc);
  ~MCCoord();

  bool computedmin;

  AthenaArray<Real> x1f, x2f, x3f; // face  positions
  AthenaArray<Real> vol;
  AthenaArray<Real> dmin;

  virtual void Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
  virtual void InverseMetric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
  virtual void Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);

  Real GetMass() const {return bh_mass_;}
  Real GetSpin() const {return bh_spin_;}
  void SetMass(Real m) {bh_mass_ = m;}
  void SetSpin(Real a) {bh_spin_ = a;}

protected:
 // GR-specific variables
  Real bh_mass_;
  Real bh_spin_;

};

//----------------------------------------------------------------------------------------
//! \class MCCartesian
//! \brief derived class for Cartesian coordinates

class MCCartesian : public MCCoord {
public:
  MCCartesian(Coordinates *pcoord, MonteCarloBlock *pmcb);
  MCCartesian(int ncells1, int ncells2, int ncells3, bool acc);
  ~MCCartesian();
};

//----------------------------------------------------------------------------------------
//! \class MCSphericalPolar
//! \brief derived class for Spherical coordinates

class MCSphericalPolar : public MCCoord {
public:
  MCSphericalPolar(Coordinates *pcoord, MonteCarloBlock *pmcb);
  MCSphericalPolar(int ncells1, int ncells2, int ncells3, bool acc);
  ~MCSphericalPolar();

  // functions
  void Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
  void Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
};

//----------------------------------------------------------------------------------------
//! \class MCKerrSchild
//! \brief derived class for Kerr-Schild coordinates

class MCKerrSchild: public MCCoord {
public:
  MCKerrSchild(Coordinates *pcoord, MonteCarloBlock *pmcb);
  MCKerrSchild(int ncells1, int ncells2, int ncells3, bool acc);
  ~MCKerrSchild();

  // functions
  void Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
  void Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
  void InverseMetric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);

};

//----------------------------------------------------------------------------------------
//! \class MCKerrSchildCartesian
//! \brief derived class for cartesian Kerr-Schild coordinates

class MCKerrSchildCartesian: public MCCoord {
public:
  MCKerrSchildCartesian(Coordinates *pcoord, MonteCarloBlock *pmcb);
  MCKerrSchildCartesian(int ncells1, int ncells2, int ncells3, bool acc);
  ~MCKerrSchildCartesian();

  // functions
  void Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
  void Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
  void InverseMetric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);

};

//----------------------------------------------------------------------------------------
//! \class MCCylindrical
//! \brief derived class for Cylindrical coordinates

class MCCylindrical : public MCCoord {
public:
  MCCylindrical(Coordinates *pcoord, MonteCarloBlock *pmcb);
  MCCylindrical(int ncells1, int ncells2, int ncells3, bool acc);
  ~MCCylindrical();

  // functions
  void Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
  void Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);

};

//----------------------------------------------------------------------------------------
//! \class MCBoyerLindquist
//! \brief derived class for Boyer-Lindquist coordinates

class MCBoyerLindquist : public MCCoord {
public:
  MCBoyerLindquist(Coordinates *pcoord, MonteCarloBlock *pmcb);
  MCBoyerLindquist(int ncells1, int ncells2, int ncells3, bool acc);
  ~MCBoyerLindquist();

  // functions
  void Connect(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
  void Metric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
  void InverseMetric(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
};

//----------------------------------------------------------------------------------------
//! \class MCMinkowski
//! \brief derived class for Minkowski coordinates

class MCMinkowski : public MCCoord {
public:
  MCMinkowski(Coordinates *pcoord, MonteCarloBlock *pmcb);
  MCMinkowski(int ncells1, int ncells2, int ncells3, bool acc);
  ~MCMinkowski();
};

#endif // MCCOORD_HPP
