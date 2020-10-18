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
#include "../coordinates/coordinates.hpp"
#include "../outputs/outputs.hpp"
#include "photon.hpp"
#include "mcbvals.hpp"
#include "mcoutput.hpp"
#include "mccoord.hpp"

class Mesh;
class MeshBlock;
class MonteCarloBlock;
class ParameterInput;
class Photon;
class PhotonMover;
class MCRandom;
class MCBoundaryValues;
class MCOutoupt;
class MCCoord;

// SWD: Make into a general MACRO set by configure?
#define NCOORD 4

// Flags for controlling monte carlo emission, scattering, absorption, bcs
enum EmissionFlag {EMISUSER = 0, EMISFF = 1};
enum AbsorptionFlag {ABSUSER = 0, ABSNONE = 1, ABSFF = 2};
enum ScatteringFlag {SCATUSER = 0, SCATNONE =1, SCATISO = 2, SCATTHOM = 3, SCATCOMP =4};
enum MCBoundaryFlag {MC_PERIODIC_BNDRY = 0, MC_ESCAPE_BNDRY = 1, MC_ABSORB_BNDRY = 2,
                     MC_POLAR_BNDRY = 3, MC_REFLECT_BNDRY = 4, MC_USER_BNDRY = 5,
                     MC_BLOCK_BNDRY = 6};
// Array indices for monte carlo radiation moments
enum {MCIER=0, MCIFR1=1, MCIFR2=2, MCIFR3=3, MCIPR11=4, MCIPR22=5, MCIPR33=6,
      MCIPR12=7, MCIPR13=8, MCIPR23=9, MCIEN = 10, MCIPR21=11, MCIPR31=12, 
      MCIPR32=13};

//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime
typedef void (*EmisFunc_t)(MonteCarloBlock *pmcb);
typedef void (*TempFunc_t)(MonteCarloBlock *pmcb);
typedef void (*MCBValFunc_t)(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot);
typedef Real (*OpacFunc_t)(MonteCarloBlock *pmcb, Photon *phot);
typedef void (*ScatFunc_t)(MonteCarloBlock *pmcb, Photon *phot);
typedef void (*UserMoveFunc_t)(MonteCarloBlock *pmcb, Photon *phot, PhotonMover *pmover);
typedef void (*GetZonePos_t)(Photon *phot, MCRandom *pran, MCCoord *pco);
typedef void (*ConnectFunc_t)(Real *x, Real gamma[NCOORD][NCOORD][NCOORD]);
typedef void (*MCMetricFunc_t)(Real *x,  Real gcov[NCOORD][NCOORD]);

//---------------------- prototypes for provided functions -------------------------------
void DefaultGetTemperature(MonteCarloBlock *pmcb);
//--------------------- prototypes for opacity.cpp functions -----------------------------
Real NoOpacity(MonteCarloBlock *pmcb, Photon *pphot);
Real FreeFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot);
Real ThomsonOpacity(MonteCarloBlock *pmcb, Photon *pphot);
Real ComptonOpacity(MonteCarloBlock *pmcb, Photon *pphot);
void GenerateComptonTable(int io);
Real ComptonCrossSection(Real energy, Real theta);
Real Maxwell(Real theta, Real gamma);
Real KleinNishina(Real x);
void InitializeAccelerationOpacity(MonteCarloBlock *pmcb);
//--------------------- prototypes for scatter.cpp functions -----------------------------
void NoScatter(MonteCarloBlock *pmcb, Photon *pphot);
void ScatterIsotropic(MonteCarloBlock *pmcb, Photon *pphot);
void ScatterThomsonPolarized(MonteCarloBlock *pmcb, Photon *pphot);
void ScatterThomsonUnpolarized(MonteCarloBlock *pmcb, Photon *pphot);
void ScatterComptonUnpolarized(MonteCarloBlock *pmcb, Photon *pphot);
void ScatterComptonPolarized(MonteCarloBlock *pmcb, Photon *pphot);
Real Bigy(Real x, Real xp);
Real SigmaHat(Real x);
Real ElectronDistOld(Real tgas, MCRandom *pran);
Real ElectronDist(Real tgas, MCRandom *pran);
//--------------------- prototypes for emission.cpp functions ----------------------------
void InitializeEmissionFreeFree(MonteCarloBlock *pmcb);
void PhotonEmitFreeFree(MonteCarloBlock *pmcb, Photon *pphot);
Real PlanckDist(Real temp,MCRandom *pran);
void GetZonePositionCartesian(Photon *pphot, MCRandom *pran, MCCoord *pco);
void GetZonePositionSphericalPolar(Photon *pphot, MCRandom *pran, MCCoord *pco);
void GetZonePositionSphericalPolarGR(Photon *pphot, MCRandom *pran, MCCoord *pcoord);
void GetZonePositionCylindrical(Photon *pphot, MCRandom *pran, MCCoord *pcoord);
void GetZonePositionCylindricalGR(Photon *pphot, MCRandom *pran, MCCoord *pcoord);
//--------------------- protoypes for grmover.cpp functions ------------------------------
void Metric_KerrSchild(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
void Metric_KerrSchild_Up(Real x[NCOORD], Real gcon[NCOORD][NCOORD]);
void Metric_BoyerLindquist(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
void Metric_BoyerLindquist_Up(Real x[NCOORD], Real gcon[NCOORD][NCOORD]);
void Metric_Cartesian(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
void Metric_SphericalPolar(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
void Metric_Cylindrical(Real x[NCOORD], Real gcov[NCOORD][NCOORD]);
void Connect_KerrSchild(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
void Connect_BoyerLindquist(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
void Connect_Cartesian(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
void Connect_SphericalPolar(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
void Connect_Cylindrical(Real x[NCOORD], Real gamma[NCOORD][NCOORD][NCOORD]);
void GetMCDirection(Photon *pphot, Real alpha, Real beta);
//------------------ prototypes for frame_transformations.cpp functions ------------------
void ConstructTetrad(Real ucon[NCOORD], Real gcov[NCOORD][NCOORD],
		      Real econ[NCOORD][NCOORD], Real ecov[NCOORD][NCOORD]);
int KroneckerDelta(int i, int j);
void ProjectVecSub(Real ucon[NCOORD], Real vcon[NCOORD], Real gcov[NCOORD][NCOORD]);
Real DotVec(Real ucon[NCOORD], Real vcon[NCOORD], Real gcov[NCOORD][NCOORD]);
void NormalizeVec(Real ucon[NCOORD], Real gcov[NCOORD][NCOORD]);
void ConToCov(Real ucon[NCOORD], Real ucov[NCOORD], Real gcov[NCOORD][NCOORD]);
void CovToCon(Real ucov[NCOORD], Real ucon[NCOORD], Real gcon[NCOORD][NCOORD]);
void CoordinateToTetrad(Real ucoord[NCOORD], Real utet[NCOORD], Real ecov[NCOORD][NCOORD]);
void TetradToCoordinate(Real utet[NCOORD], Real ucoord[NCOORD], Real econ[NCOORD][NCOORD]);
//---------------------- prototypes for setting flags ------------------------------------
enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string);
enum EmissionFlag GetEmissionFlag(std::string input_string);
enum AbsorptionFlag GetAbsorptionFlag(std::string input_string);
enum ScatteringFlag GetScatteringFlag(std::string input_string);

//----------------------------------------------------------------------------------------
//! \struct MCBlockSize
//  \brief physical size of monte carlo block

typedef struct MCBlockSize {
  int nx1,nx2,nx3;
  int is,ie,js,je,ks,ke;

} MCBlockSize;

//----------------------------------------------------------------------------------------
//! \class MCRandom
//  \brief monte carlo random number generator

class MCRandom {
public:
  MCRandom(int iseed);
  ~MCRandom();

  gsl_rng *dev;
  
  Real uniform();
  Real chisquare(int n);
};

//----------------------------------------------------------------------------------------
//! \class MonteCarlo
//  \brief monte carlo functions and data

class MonteCarlo {
  friend class MCBoundaryValues;

public:
  MonteCarlo(ParameterInput *pin, Mesh *pmesh);
  ~MonteCarlo();

  // data
  Mesh *pmy_mesh;
  MCOutput *pmcout;
  MonteCarloBlock *pblock;

  int nphtot;  // total number of photons to integrate
  int cadence; // number of photons per output
  int nout;  // number of outputs
  int *nphlist; // number of photons per block
  int64_t ncells; // total number of cells in mesh
  int iseed;  // seed to initialized random number generator(s)
  int mcranks; // number of monte carlo only ranks
  int nderv; // number of derivative processes
  int *derv; // pointer to array of derivative processes
  int origin; // process with origin mesh block
  int blocksize;
  int nphrun; // number of photons run thus far
  int max_list_size; // maximum number of photons run per output on any process
  int nuser_var;

  enum EmissionFlag emission_meth;
  enum AbsorptionFlag absorption_meth;
  enum ScatteringFlag scattering_meth;
  enum MCBoundaryFlag mc_bcs[6];

  bool boosts;  // Compute lorentz transformations
  bool emission_array_flag;  // Compute and save zone emissivities
  bool polarized;// track photon polarization
  bool acceleration;  // use MRW acceleration
  bool time_acc;  // use MRW acceleration with time limit
  bool raytrace_flag; // Will trace photons rather than scatter

  EmisFunc_t InitEmission;
  TempFunc_t GetTemperature;

  // functions
  void RunStaticMonteCarlo(Outputs *pouts, Mesh *pmesh, ParameterInput *pinput);
  void RunStaticMonteCarloNew(void);
  void InitUserMonteCarloData(ParameterInput *pin);
  void EnrollUserEmissionInitialization(EmisFunc_t emissfunc);
  void EnrollUserGetTemperature(TempFunc_t tempfunc);
  void SendMonteCarloSpectra(int dest);
  void ReceiveMonteCarloSpectra(int source);
  void CollectMoments(void);
  void EnrollUserMCBoundaryFunction(enum BoundaryFace dir, MCBValFunc_t my_bc);

private:

  // functions
  MCBValFunc_t BoundaryFunction_[6];

  void GetDensity(MonteCarloBlock *pmcb);
  void GetVelocity(MonteCarloBlock *pmcb);
  void DistributePhotonsToBlocks(void);
  void SendMonteCarloBlocks(int dest);
  void ReceiveMonteCarloBlocks(ParameterInput *pin, int source);
  void SendMonteCarloData(int dest);
  void ReceiveMonteCarloData(int source);
  unsigned int CreateMCMPITag(int bid);
  void InitializeMonteCarloBlocks(void);
  void SendMoments(int dest);
  void ReceiveMoments(int source, bool sum);

};

//----------------------------------------------------------------------------------------
//! \class MonteCarloBlock
//  \brief monte carlo functions and data contained on each mesh block

class MonteCarloBlock {
public:
  MonteCarloBlock(MeshBlock *pmb, MCBlockSize *pblsize, MonteCarlo *pmc, 
                  ParameterInput *pin);
  ~MonteCarloBlock();

  // data
  MonteCarlo* pmy_mc; // MonteCarlo
  MeshBlock* pmy_block;    // MeshBlock corresponding to this MonteCarlo
  MCCoord *pcoord;

  Photon* pphoton; // ptr to photon packet
  PhotonMover* pmover; // ptr to photon mover
  MCRandom *pran; // ptr to random number generator
  MCBoundaryValues *pbval; // ptr to MC boundary values
  Spectrum *pspec; // ptr to spectrum
  PhotonList *pphlist; // ptr to photon list

  MonteCarloBlock *next;

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
  //UserMoveFunc_t UserWorkInMove;
  ConnectFunc_t Connection;
  MCMetricFunc_t Metric;

  int nphdone; // Photons integrated thus far
  int nphremain; // total number of photons to integrate
  int myblockid;
  int nx1,nx2,nx3;
  int is,ie,js,je,ks,ke;
  int nfreq, nmu, nphi, nsurf;
  int cadence;

  bool zone_weight_flag; // flag for zone weighting
  bool weighted_absorption; // flag controling how absorption is handled
  bool moments_flag; // Compute/output moments
  bool emission_array_flag;  // Compute and save zone emissivities
  bool boosts;  // Compute lorentz transformations
  bool coherent_scattering; // photon does notchange energy after scattering
  //bool polarized; // track photon polarization
  bool acceleration;  // use MRW acceleration
  bool time_acc;  // use MRW acceleration with time limit

  // Associated with general mover
  // SWD some of these should be eliminated others moved to MonteCarlo
  bool general_mover_flag; // use general integration for mover
  bool kerrschild_flag; // use KerrSchild coordinates and BH test
  bool boyerlindquist_flag; // use Boyer-Lindquist coordinates 
  bool orthotet_flag; // use orthonormal tetrad for TransferPhotons()
  bool varystep_flag; // use variable (true) or constant (false) step

  Real codetocgs_rho, codetoc_vel;
  Real emin, emax, elog, eminlog;
  // SWD: used by general mover, move/eliminate 
  Real stepsize, a, velocity;

  AthenaArray<Real> emission;
  AthenaArray<Real> moments;
  AthenaArray<Real> rho;
  AthenaArray<Real> tgas;
  AthenaArray<Real> vel;
  AthenaArray<Real> planck_opacity; // for acceleration
  AthenaArray<Real> planck_inv_opacity; // for acceleration

  // functions
  void InitUserMonteCarloBlockData(ParameterInput *pin);
  void MonteCarloProblemGenerator(ParameterInput *pin);
  void RayTracePhotons(int nphtot); // Ray trace photon on this block
  void TransferPhotons(int nphtot); // Transfer photons on this block
  void LorentzTransform(Photon *pphot, const Real sign);
  Real LorentzTransformFrequencyShift(Photon *pphot);
  void InitializePhoton(Photon *pphot);
  void FinalizePhoton(Photon *pphot);
  //void DefaultGetTemperature();
  void UpdateMoments(Photon *pphot, Real dl);
  void NormalizeMoments(bool normalize);
  void ResetMoments();
  //void GetPhotonsFromNeighbors();
  //void SendPhotonsToNeighbors();
  void EnrollUserWorkInMove(UserMoveFunc_t userfunc);
  void TetradTransform(Photon *pphot, const Real sign);
  Real TetradTransformFrequencyShift(Photon *pphot);

private:
   void SetBoundaryValues(enum MCBoundaryFlag *input_bcs);
 
};



#endif // MONTECARLO_HPP
