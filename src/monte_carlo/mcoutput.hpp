#ifndef MCOUTPUT_HPP
#define MCOUTPUT_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photon.hpp
//  \brief definitions for MCOutput class

// Athena++ classes headers
#include "../athena.hpp"
#include "montecarlo.hpp"

class Spectrum {
public:
  Spectrum(Real emin, Real emax, int nfreq, int nmu, int phi, bool polarized);
  ~Spectrum();

  int nfreq, nmu, nphi;
  Real emin, emax;
  bool polarized;

  AthenaArray<Real> energies;
  AthenaArray<Real> intensity;
  AthenaArray<Real> intensity_sq;
  AthenaArray<Real> stokesq;
  AthenaArray<Real> stokesq_sq;
  AthenaArray<Real> stokesu;
  AthenaArray<Real> stokesu_sq;

  //functions
  void BuildFrequencyGrid(Real emin, Real emax, int nfreq);
  void UpdateSpectrum(Photon *pphot);
  int GetEbin(Real energy);

};

class MCOutput {
public:
  MCOutput(MonteCarlo *pmc, ParameterInput *pin);
  ~MCOutput();

  MonteCarlo *pmy_mc;
  void OutputSpectra(MonteCarlo *pmc);
  void OutputSpectrum(Spectrum *pspec, Real norm, std::string outfile);

};

#endif //MC_OUTPUT
