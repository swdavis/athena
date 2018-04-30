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
#include "photon.hpp"

class Mesh;
class MeshBlock;
class MonteCarloBlock;
class ParameterInput;
class Photon;
class PhotonMover;

// Flags for controlling monte carlo emission, scattering, and absorption
enum EmissionFlag {EMISUSER = 0, EMISFF = 1};
enum AbsorptionFlag {ABSUSER = 0, ABSNONE = 1, ABSFF = 2};
enum ScatteringFlag {SCATUSER = 0, SCATNONE =1, SCATISO = 2, SCATTHOM = 3, SCATCOMP =4};
enum {TOEUL=0, TOCOM=1};

//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime
typedef void (*EmisFunc_t)(MonteCarloBlock *pmcb);
//typedef void (*EmisFunc_t)(MonteCarloBlock *pmcb, MeshBlock *pmb);
typedef Real (*OpacFunc_t)(MeshBlock *pmb);
typedef void (*ScatFunc_t)(MeshBlock *pmb);

//---------------------- prototypes for default functions --------------------------------
//void InitializeEmissionFreeFree(MonteCarloBlock *pmcb, MeshBlock *pmb);
void InitializeEmissionFreeFree(MonteCarloBlock *pmcb);

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
  MonteCarloBlock *pblock;
  int ntot;         // total number of photons to integrate
  enum EmissionFlag emission_meth;
  EmisFunc_t InitEmission;

  // functions
  void LaunchPhotons();
  void EnrollUserEmissionInitialization(EmisFunc_t emissfunc);

private:
  enum EmissionFlag GetEmissionFlag(std::string input_string);
};

//! \class MonteCarloBlock
//  \brief monte carlo functions and data contained on each mesh block

class MonteCarloBlock {
public:
  MonteCarloBlock(MeshBlock *pmb, MonteCarlo *pmc, ParameterInput *pin);
  ~MonteCarloBlock();

  // data
  MonteCarlo* pmy_mc; // MonteCarlo
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this MonteCarlo
  Coordinates *pmy_coord;

  Photon* pphoton; // ptr to photon packet
  PhotonMover* pmover; // ptr to photon mover

  MonteCarloBlock *prev, *next;

  enum EmissionFlag emission_meth;
  enum AbsorptionFlag absorption_meth;
  enum ScatteringFlag scattering_meth;

  int ntot;  // total number of photons for this block;
  int is,ie,js,je,ks,ke;

  bool zone_weight_flag; // flag for zone weighting
  bool moments_flag; // Compute/output moments
  bool emission_array_flag;  // Compute and save zone emissivities
  bool lorentz_trans_flag;  // Compute lorentz transformations

  MCRandom *pran;

  AthenaArray<Real> emission;
  AthenaArray<Real> moments;
  AthenaArray<Real> rho;
  AthenaArray<Real> tgas;
  AthenaArray<Real> vel;

  // functions
  void MonteCarloProblemGenerator(ParameterInput *pin);
  void TransferPhotons();  // Transfer photons on this block
  void InitializePhoton(MeshBlock *pmb, Photon *pphot);
  //void EnrollUserEmissionInitialization(EmisFunc_t emissfunc);
  //void InitializeEmissionFreeFree();

private:

  //enum EmissionFlag GetEmissionFlag(std::string input_string);
  enum AbsorptionFlag GetAbsorptionFlag(std::string input_string);
  enum ScatteringFlag GetScatteringFlag(std::string input_string);

};



#endif // MONTECARLO_HPP
