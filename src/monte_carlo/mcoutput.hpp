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


//----------------------------------------------------------------------------------------
//! \struct MomentumRange
//  \brief physical ranges of photon momentum grid

typedef struct MomentumRange {
  int ne, ncth, nphi;
  Real emin, emax;
  Real phimin, phimax;
  Real cthmin, cthmax;

} MomentumRange;

//----------------------------------------------------------------------------------------
//! \class Spectrum
//  \brief spectral bins for outputs

class Spectrum {
public:
  Spectrum(MomentumRange input_range, bool polarized);
  Spectrum(Spectrum *pspec);
  ~Spectrum();

  MomentumRange range;
  bool polarized;
  Spectrum *next;  // next spectrum
  enum BoundaryFace face;

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
  void SetSurface(std::string input_face);

};

//----------------------------------------------------------------------------------------
//! \class MCOutput
//  \brief class for handling monte carlo specific spectral outputs

class MCOutput {
public:
  MCOutput(MonteCarlo *pmc, ParameterInput *pin);
  ~MCOutput();

  MonteCarlo *pmy_mc;
  Spectrum *pspec;
  bool moments;

  //functions
  void OutputSpectra(MonteCarlo *pmc);
  void OutputSpectrum(Spectrum *pspec, Real norm, std::string outfile);

};

#endif //MC_OUTPUT
