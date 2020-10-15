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

class Photon;

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
  Spectrum(MomentumRange input_range, bool polarized, bool logarithmic);
  Spectrum(Spectrum *pspec);
  ~Spectrum();

  std::string base_name;
  MomentumRange range;
  bool polarized;
  bool cartesian_axis;
  bool coordinates;
  bool logarithmic;
  bool pathbin; // Replace energy bin with path length bin
  bool radbin; // Replace energy bin with radius bin
  bool legacy; // Output format to be used (temporarily retained for testing)
  Spectrum *next;  // next spectrum
  enum BoundaryFace face;
  int id;
  int output_number;// current output number
  Real x1min,x1max,x2min,x2max,x3min,x3max;

  AthenaArray<Real> energies;
  AthenaArray<Real> intensity;
  AthenaArray<Real> intensity_sq;
  AthenaArray<Real> stokesq;
  AthenaArray<Real> stokesq_sq;
  AthenaArray<Real> stokesu;
  AthenaArray<Real> stokesu_sq;

  //functions
  void BuildEnergyGrid(Real emin, Real emax, int nen, bool xlog);
  void UpdateSpectrum(Photon *pphot);
  int EnergyBin(Real energy);
  bool AngleBinsCartesian(Photon *pphot, int &Phibin, int &mubin);
  bool AngleBinsSphericalPolar(Photon *pphot, int &Phibin, int &mubin);
  void SetSurface(std::string input_face);
  bool ScreenCoordinates(Photon *pphot);
  void ResetSpectrum();
  void AddSpectrum(Spectrum *pspec);
  void WriteSpectrum(std::string filename, int ntot);
  void WriteSpectrumLegacy(std::string outfile, Real norm);
};

//----------------------------------------------------------------------------------------
//! \class PhotonList
//  \brief List of Photon properties (usually escaping photons)

class PhotonList {
public:
  PhotonList(int list_mem_size, bool pol, bool rel, int nuser_out);
  ~PhotonList();

  std::string base_name;

  int length; // number of occupied elements
  int nparams; // number of properties for each photon in list
  int output_number;// current output number
  int nuser_out;
  bool polarized;
  bool relativistic;
  AthenaArray<Real> photons;  // array of photon properies

  //functions
  void AddPhoton(Photon *pphot);
  void WriteList(std::string filename, int ntot);

private:
  int max_len;  // number of photons allowed with current allocated memory
  void ResizeList(int new_size); 

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
  PhotonList *pphlist;

  bool moments;

  //functions
  void CollectSpectrum(MonteCarlo *pmc);
  void OutputSpectrum(MonteCarlo *pmc);
  void OutputPhotonList(int nphtot);
};

#endif //MC_OUTPUT
