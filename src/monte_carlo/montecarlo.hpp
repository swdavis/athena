#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file montecarlo.hpp
//  \brief definitions for MonteCarlo class
//
// Current design focusses on implementing static post-processing so these class
// implementations will evolve.

#include <sstream>
#include <gsl/gsl_randist.h>

// Athena++ classes headers
#include "../athena.hpp"
//#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "photon.hpp"
#include "mcbvals.hpp"
#include "mcoutput.hpp"

class Mesh;
class MeshBlock;
class MonteCarloBlock;
class ParameterInput;
class Photon;
class PhotonMover;
class MCRandom;
class MCBoundaryValues;
class MCOutoupt;

// Flags for controlling monte carlo emission, scattering, absorption, bcs
enum EmissionFlag {EMISUSER = 0, EMISFF = 1};
enum AbsorptionFlag {ABSUSER = 0, ABSNONE = 1, ABSFF = 2};
enum ScatteringFlag {SCATUSER = 0, SCATNONE =1, SCATISO = 2, SCATTHOM = 3, SCATCOMP =4};
enum MCBoundaryFlag {MC_PERIODIC_BNDRY = 0, MC_ESCAPE_BNDRY = 1, MC_ABSORB_BNDRY = 2,
                     MC_POLAR_BNDRY = 3, MC_REFLECT_BNDRY = 4, MC_USER_BNDRY = 5};
// Array indices for monte carlo radiation moments
enum {MCIER=0, MCIFR1=1, MCIFR2=2, MCIFR3=3, MCIPR11=4, MCIPR22=5, MCIPR33=6,
      MCIPR12=7, MCIPR13=8, MCIPR23=9, MCIPR21=10, MCIPR31=11, MCIPR32=12};

enum {TOEUL=0, TOCOM=1};

//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime
typedef void (*EmisFunc_t)(MonteCarloBlock *pmcb);
typedef void (*TempFunc_t)(MonteCarloBlock *pmcb);
//typedef void (MonteCarloBlock::*TempFunc2_t)(void);
typedef Real (*OpacFunc_t)(MonteCarloBlock *pmcb, Photon *phot);
typedef void (*ScatFunc_t)(MonteCarloBlock *pmcb, Photon *phot);
typedef void (*GetZonePos_t)(Photon *phot, MCRandom *pran, Coordinates *pco);

//---------------------- prototypes for provided functions -------------------------------
void InitializeEmissionFreeFree(MonteCarloBlock *pmcb);
void DefaultGetTemperature(MonteCarloBlock *pmcb);
void PhotonEmitFreeFree(MonteCarloBlock *pmcb, Photon *pphot);
//--------------------- prototypes for opacity.cpp functions -----------------------------
Real NoOpacity(MonteCarloBlock *pmcb, Photon *pphot);
Real FreeFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot);
Real ThomsonOpacity(MonteCarloBlock *pmcb, Photon *pphot);
Real ComptonOpacity(MonteCarloBlock *pmcb, Photon *pphot);
void GenerateComptonTable(void);
Real ComptonCrossSection(Real energy, Real theta);
Real Maxwell(Real theta, Real gamma);
Real KleinNishina(Real x);
//--------------------- prototypes for scatter.cpp functions -----------------------------
void ScatterIsotropic(MonteCarloBlock *pmcb, Photon *pphot);
void ScatterThomsonPolarized(MonteCarloBlock *pmcb, Photon *pphot);
void ScatterComptonUnpolarized(MonteCarloBlock *pmcb, Photon *pphot);
Real Bigy(Real x, Real xp);
Real SigmaHat(Real x);
Real ElectronDist(Real tgas, MCRandom *pran);
//--------------------- prototypes for emission.cpp functions ----------------------------
void GetZonePositionCartesian(Photon *pphot, MCRandom *pran, Coordinates *pco);
void GetZonePositionSphericalPolar(Photon *pphot, MCRandom *pran, Coordinates *pco);
//---------------------- prototypes for setting flags ------------------------------------
enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string);
enum EmissionFlag GetEmissionFlag(std::string input_string);
enum AbsorptionFlag GetAbsorptionFlag(std::string input_string);
enum ScatteringFlag GetScatteringFlag(std::string input_string);


//! \class MCRandom
//  \brief monte carlo random number generator

class MCRandom {
public:
  MCRandom(int iseed);
  ~MCRandom();

  gsl_rng *dev;
  
  Real uniform();
};

//! \class MonteCarlo
//  \brief monte carlo functions and data

class MonteCarlo {
public:
  MonteCarlo(ParameterInput *pin, Mesh *pmesh);
  ~MonteCarlo();

  // data
  Mesh *pmy_mesh;
  MCOutput *pmcout;
  MonteCarloBlock *pblock;
  int nphot;         // total number of photons to integrate
  int iseed;  // seed to initialized random number generator(s)
  enum EmissionFlag emission_meth;
  enum AbsorptionFlag absorption_meth;
  enum ScatteringFlag scattering_meth;
  enum MCBoundaryFlag mc_bcs[6];
  bool moments_flag; // Compute/output moments
  bool lorentz_trans_flag;  // Compute lorentz transformations

  EmisFunc_t InitEmission;
  TempFunc_t GetTemperature;

  // functions
  void RunStaticMonteCarlo();
  void InitUserMonteCarloData(ParameterInput *pin);
  void EnrollUserEmissionInitialization(EmisFunc_t emissfunc);
  void EnrollUserGetTemperature(TempFunc_t tempfunc);

private:

  // functions
  void GetDensity(MonteCarloBlock *pmcb);
  void GetVelocity(MonteCarloBlock *pmcb);
  void CountCellsOnBlocks(void);

};

//! \class MonteCarloBlock
//  \brief monte carlo functions and data contained on each mesh block

class MonteCarloBlock {
public:
  MonteCarloBlock(MeshBlock *pmb, MonteCarlo *pmc, ParameterInput *pin);
  ~MonteCarloBlock();

  // data
  MonteCarlo* pmy_mc; // MonteCarlo
  MeshBlock* pmy_block;    // MeshBlock corresponding to this MonteCarlo
  Coordinates *pmy_coord;

  Photon* pphoton; // ptr to photon packet
  PhotonMover* pmover; // ptr to photon mover
  MCRandom *pran; // ptr to random number generator
  MCBoundaryValues *pbval; // ptr to MC boundary values
  Spectrum *pspec; // ptr to spectrum

  MonteCarloBlock *prev, *next;

  enum EmissionFlag emission_meth;
  enum AbsorptionFlag absorption_meth;
  enum ScatteringFlag scattering_meth;
  enum MCBoundaryFlag mcb_bcs[6];

  // function pointers
  //TempFunc2_t GetTemperature2;
  GetZonePos_t GetZonePosition;
  OpacFunc_t AbsorptionOpacity;
  OpacFunc_t ScatteringOpacity;
  ScatFunc_t Scatter;

  int nphot;  // total number of photons for this block;
  int ncells;
  int is,ie,js,je,ks,ke;
  int nfreq, nmu, nphi, nsurf;

  bool zone_weight_flag; // flag for zone weighting
  bool weighted_absorption; // flag controling how absorption is handled
  bool moments_flag; // Compute/output moments
  bool emission_array_flag;  // Compute and save zone emissivities
  bool lorentz_trans_flag;  // Compute lorentz transformations
  bool coherent_scattering; // photon does notchange energy after scattering
  bool polarized; // track photon polarization

  Real codetocgs_rho, codetoc_vel;
  Real emin, emax, elog, eminlog;

  AthenaArray<Real> emission;
  AthenaArray<Real> moments;
  AthenaArray<Real> rho;
  AthenaArray<Real> tgas;
  AthenaArray<Real> vel;


  // functions
  void MonteCarloProblemGenerator(ParameterInput *pin);
  void TransferPhotons();  // Transfer photons on this block
  void InitializePhoton(Photon *pphot);
  void DefaultGetTemperature();
  void UpdateMoments(Photon *pphot, Real dl);
  void NormalizeMoments(bool normalize);
  //void GetPhotonsFromNeighbors();
  //void SendPhotonsToNeighbors();

private:
  
 
};



#endif // MONTECARLO_HPP
