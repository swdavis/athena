#ifndef MCOUTPUT_HPP
#define MCOUTPUT_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file photon.hpp
//! \brief definitions for MCOutput and output classes

// Athena++ classes headers
#include "../athena.hpp"
#include "montecarlo.hpp"

class Photon;

//----------------------------------------------------------------------------------------
//! \struct MomentumRange
//! \brief physical ranges of photon momentum grid

typedef struct MomentumRange {
  int ne, ncth, nphi;
  Real emin, emax;
  Real phimin, phimax;
  Real cthmin, cthmax;

} MomentumRange;

//----------------------------------------------------------------------------------------
//! \class Spectrum
//! \brief spectral output

class Spectrum {
public:
  Spectrum(MomentumRange input_range, bool polarized, bool logarithmic);
  Spectrum(Spectrum *pspec);
  ~Spectrum();

  MonteCarlo *pmy_mc;
  int nsrun;  // total number of photons samples run for this spectrum
  std::string base_name;
  MomentumRange range;
  bool polarized;
  bool polar_axis;
  bool coordinates;
  bool logarithmic;

  Spectrum *next;  // next spectrum
  enum BoundaryFace face;
  int id; // spectrum id -- maybe multiple spectra requested
  int output_number;// current output number
  Real x1min,x1max,x2min,x2max,x3min,x3max;
  Real dt; // targe integration time for this spectrum
  Real last_time;

  AthenaArray<Real> energies;
  AthenaArray<Real> count;
  AthenaArray<Real> intensity;
  AthenaArray<Real> intensity_sq;
  AthenaArray<Real> stokesq;
  AthenaArray<Real> stokesq_sq;
  AthenaArray<Real> stokesu;
  AthenaArray<Real> stokesu_sq;

  //functions
  void BuildEnergyGrid(Real emin, Real emax, int nen, bool xlog);
  void UpdateSpectrum(Photon *pphot, int ip);
  int EnergyBin(Real energy);
  int EnergyBinUniform(Real energy, bool loge);
  bool AngleBinsCartesian(Real k[4], int &Phibin, int &mubin);
  bool AngleBinsSphericalPolar(Real k[4], int &Phibin, int &mubin);
  void SetSurface(std::string input_face);
  enum BoundaryFace GetPhotonFace(Photon *phot, int ip);
  bool ScreenCoordinates(Photon *pphot, int ip);
  void ResetSpectrum();
  void AddSpectrum(Spectrum *pspec);
  void WriteSpectrum(std::string filename);

};

//----------------------------------------------------------------------------------------
//! \class PhotonList
//! \brief List of output Photon properties

class PhotonList {
public:
  PhotonList(int list_size_init, bool pol, int nuser);
  ~PhotonList();

  MonteCarlo *pmy_mc;
  std::string base_name;

  int nsrun;  // total number of photons samples run for this list
  int length; // number of occupied elements
  int nparams; // number of properties for each photon in list
  int output_number;// current output number
  int nuser_out;
  bool polarized;
  Real dt; // targe integration time for this spectrum
  Real last_time;
  AthenaArray<Real> photons;  // array of photon properies

  //functions
  void AddPhoton(Photon *pphot, int ip);
  void WriteList(std::string filename);
  void ResetList();

private:
  int len_limit;  // number of photons allowed with current allocated memory
  void ResizeList(int new_size);

};

//----------------------------------------------------------------------------------------
//! \class PhotonTrajectoryList
//  \brief List of photon trajectories

class PhotonTrajectoryList {
public:
  PhotonTrajectoryList(int init_len_limit, int init_step_limit, int nuser);
  ~PhotonTrajectoryList();

  std::string base_name;

  int length; // number of trajectories
  int maxstep;
  int nparams; // number of properties for each photon trajectory
  int output_number;// current output number
  int nuser_out;

  int *nsteps;  // step number for each trajectory
  AthenaArray<Real> trajectories;  // array of photon properies

  //functions
  void InitializeTrajectory(int itraj);
  void CompleteTrajectory(int itraj);
  void AddToTrajectory(Photon *pphot, int ip);
  void WriteList(std::string filename);
  //void ResetList();

private:
  int len_limit;  // number of trajectories allowed with current allocated memory
  int step_limit; // maximum number of steps
  void ResizeList(int new_len_limit, int new_step_limit);

};

//----------------------------------------------------------------------------------------
//! \class Image
//  \brief Image created using ray traced photons (SWD: in progress)

class Image {
public:
  Image(int list_mem_size, bool pol, bool rel, int nuser);
  ~Image();

  std::string base_name;

  int nx; // number of horizontal pixels
  int ny; // number of vertical pixels
  Real xcam[4]; // position of camera (tetrad)
  Real kcam[4]; // Camera direction
  Real xmin, xmax; // horizontal angular extent of image
  Real ymin, ymax; // vertical angular extent of image
  int nparams; // number of properties stored for each pixel
  int output_number;// current output number
  int nuser_out;
  bool polarized;
  bool relativistic;
  AthenaArray<Real> image;  // pixel array

  //functions
  void WriteImage(std::string filename);

  //private:

};

// SWD: Will need to be modifed with new parallelization scheme
//----------------------------------------------------------------------------------------
//! \class MCOutput
//! \brief class for managing monte carlo specific spectral outputs

class MCOutput {
public:
  MCOutput(MonteCarlo *pmc, ParameterInput *pin);
  ~MCOutput();

  MonteCarlo *pmy_mc;
  Spectrum *pspec;
  PhotonList *pphlist;
  PhotonTrajectoryList *ptraj;

  bool mom_flag_lab;
  bool mom_flag_com;
  bool mom_flag_src;
  bool mom_flag_usr;

  //functions
  void OutputSpectrum(bool wtflag);
  void SendMonteCarloSpectrum(Spectrum *spect, int dest);
  void ReceiveMonteCarloSpectrum(Spectrum *spect, bool add);
  void OutputPhotonList(bool wtflag);
  void OutputTrajectoryList();
  void UpdateOutputCount(int nph);
  void MakeOutputs(bool wtflag);
};

#endif //MC_OUTPUT
