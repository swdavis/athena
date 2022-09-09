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
#include <complex>

// Athena++ classes headers
#include "../athena.hpp"
#include "../coordinates/coordinates.hpp"
#include "../outputs/outputs.hpp"
#include "photon.hpp"
#include "mcbvals.hpp"
#include "mcoutput.hpp"
#include "mccoord.hpp"

// GSL library
#if RAN3 == 0
#include <gsl/gsl_randist.h>
#endif

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
#define NMOM 16

// Flags for controlling monte carlo emission, scattering, absorption, bcs
enum EmissionFlag {EMISUSER = 0, EMISNONE = 1, EMISFF = 2};
enum AbsorptionOpacityFlag {ABSUSER = 0, ABSNONE = 1, ABSFF = 2, ABSDUST =3};
enum AbsorptionMethodFlag {ABSWEIGHT = 0, ABSPROB = 1, ABSTAU = 2};
enum ScatteringFlag {SCATUSER = 0, SCATNONE =1, SCATISO = 2, SCATTHOM = 3, SCATCOMP =4,
                     SCATRES = 5, SCATDUST = 6};
enum MCBoundaryFlag {MC_PERIODIC_BNDRY = 0, MC_ESCAPE_BNDRY = 1, MC_ABSORB_BNDRY = 2,
                     MC_DESTROY_BNDRY = 3, MC_POLAR_BNDRY = 4, MC_REFLECT_BNDRY = 5,
                     MC_USER_BNDRY = 6, MC_BLOCK_BNDRY = 7};
// Array indices for monte carlo radiation moments
enum {MCIER=0, MCIFR1=1, MCIFR2=2, MCIFR3=3, MCIPR11=4, MCIPR22=5, MCIPR33=6,
      MCIPR12=7, MCIPR13=8, MCIPR23=9, MCINET = 10, MCIEN = 11, MCIKJ = 12, MCIPR21=13,
      MCIPR31=14, MCIPR32=15};
//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime
typedef Real (*EmisFunc_t)(MonteCarloBlock *pmcb);
typedef void (*TempFunc_t)(MonteCarloBlock *pmcb);
typedef void (*MCBValFunc_t)(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot, int ip);
typedef Real (*OpacFunc_t)(MonteCarloBlock *pmcb,  Photon *pphot, int ip);
typedef void (*ScatFunc_t)(MonteCarloBlock *pmcb, Photon *phot, int ips, int ipe);
typedef void (*UserMoveFunc_t)(MonteCarloBlock *pmcb, Photon *phot, PhotonMover *pmover,
                               int ip);
typedef void (*GetZonePos_t)(Photon *phot, MCRandom *pran, MCCoord *pco, int ip);

//---------------------- prototypes for provided functions -------------------------------
void DefaultGetTemperature(MonteCarloBlock *pmcb);
//--------------------- prototypes for opacity.cpp functions -----------------------------
Real NoOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip);
Real FreeFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip);
Real DustAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip);
Real ThomsonOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip);
Real ComptonOpacity(MonteCarloBlock *pmcb,  Photon *pphot, int ip);
Real ResonanceLineOpacity(MonteCarloBlock *pmcb,  Photon *pphot, int ip);
Real DustScatteringOpacity(MonteCarloBlock *pmcb,  Photon *pphot, int ip);
void GenerateComptonTable(int io);
Real ComptonCrossSection(Real energy, Real theta);
Real Maxwell(Real theta, Real gamma);
Real KleinNishina(Real x);
Real ResLinePre();
Real XsecLorentzian(Real nu);
Real XsecDoppler(Real nu, Real tgas);
Real XsecVoigt(Real nu, Real tgas);
void InitializeAccelerationOpacity(MonteCarloBlock *pmcb);
//--------------------- prototypes for scatter.cpp functions -----------------------------
void NoScatter(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);
void ScatterIsotropic(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);
void ScatterThomsonPolarized(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);
void ScatterThomsonUnpolarized(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);
void ScatterComptonUnpolarized(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);
void ScatterComptonPolarized(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);
void ScatterResonanceLine(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);
void ScatterDust(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);
Real Bigy(Real x, Real xp);
Real SigmaHat(Real x);
Real ElectronDistPozdnyakov(Real tgas, MCRandom *pran);
Real ElectronDist(Real tgas, MCRandom *pran);
void SampleDipole(Real theta_in, Real phi_in, Real &theta_out, Real &phi_out,
                  MCRandom *pran);
Real SampleVelocityParallel(Real a, Real x_in, MCRandom *pran);
//--------------------- prototypes for emission.cpp functions ----------------------------
Real InitializeEmissionFreeFree(MonteCarloBlock *pmcb);
void PhotonEmitFreeFree(MonteCarloBlock *pmcb, Photon *pphot, Real lemin, Real lemax,
                        int ip);
Real PlanckDist(Real temp,MCRandom *pran);
void GetZonePositionCartesian(Photon *pphot, MCRandom *pran, MCCoord *pco, int ip);
void GetZonePositionSphericalPolar(Photon *pphot, MCRandom *pran, MCCoord *pco, int ip);
void GetZonePositionCylindrical(Photon *pphot, MCRandom *pran, MCCoord *pco, int ip);

//------------------ prototypes for frame_transformations.cpp functions ------------------
// SWD:  Add these to MCCoord class, utils, keep here?
void ConstructTetrad(Real ucon[4], Real gcov[4][4],
                     Real econ[4][4], Real ecov[4][4]);
void ConstructTetrad(Real ucon[4], Real vcon[4], Real gcov[4][4],
                     Real econ[4][4], Real ecov[4][4]);
void ConstructTetrad(Real ucon[4], Real vcon[4], Real wcon[4],
                     Real gcov[4][4], Real econ[4][4],
                     Real ecov[4][4]);
void InitializeLeviCivita(Real levi[4][4][4][4]);
void ImposeRightHanded(Real econ[4][4], Real gcov[4][4]);
Real KroneckerDelta(int i, int j);
void ProjectVecSub(Real ucon[4], Real vcon[4], Real gcov[4][4]);
Real DotVec(Real ucon[4], Real vcon[4], Real gcov[4][4]);
void NormalizeVec(Real ucon[4], Real gcov[4][4]);
void ConToCov(Real ucon[4], Real ucov[4], Real gcov[4][4]);
void CovToCon(Real ucov[4], Real ucon[4], Real gcon[4][4]);
void CoordinateToTetrad(Real ucoord[4],Real utet[4],Real ecov[4][4]);
void TetradToCoordinate(Real utet[4],Real ucoord[4],Real econ[4][4]);
void StokesToTensor(Real stokes[4], std::complex<Real> tensor[4][4]);
void TensorToStokes(std::complex<Real> tensor[4][4], Real stokes[4]);

//---------------------- prototypes for setting flags ------------------------------------
enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string);
enum EmissionFlag GetEmissionFlag(std::string input_string);
enum AbsorptionOpacityFlag GetAbsorptionOpacityFlag(std::string input_string);
enum AbsorptionMethodFlag GetAbsorptionMethodFlag(std::string input_string);
enum ScatteringFlag GetScatteringFlag(std::string input_string);


//----------------------------------------------------------------------------------------
//! \struct MCBlockSize
//! \brief physical size of monte carlo block

typedef struct MCBlockSize {
  int nx1,nx2,nx3;
  int is,ie,js,je,ks,ke;

} MCBlockSize;

//----------------------------------------------------------------------------------------
//! \class MCRandom
//! \brief monte carlo random number generator

class MCRandom {
public:
  MCRandom(int iseed);
  ~MCRandom();


  Real uniform();
  Real chisquare(int n);

private:
  long r3seed;
#if RAN3 == 0
  gsl_rng *dev;
#endif
  Real ran3(long *idum);
};

//----------------------------------------------------------------------------------------
//! \class MonteCarlo
//! \brief monte carlo functions and data

class MonteCarlo {
  friend class MCBoundaryValues;

public:
  MonteCarlo(ParameterInput *pin, Mesh *pmesh);
  ~MonteCarlo();

  // data
  Mesh *pmy_mesh;
  MCOutput *pmcout;
  AthenaArray<MonteCarloBlock*> my_blocks;

  Real dt;     // Monte Carlo timestep
  Real tmax;   // Maximum evolution time
  int nphtot;  // total number of photons to integrate
  int nphdone; // total photons completed accross all blocks on all processes
  int nblock;  // number of photons per step per block
  int nblocal; // number of montecarloblocks on this process
  int nbtotal; // total number of montecarloblocks
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
  int max_phots_init; // maximum number of photon elements
  int nuser_var;
  int checkmove,checkscat;

  enum EmissionFlag emission_meth;
  enum MCBoundaryFlag mc_bcs[6];

  bool boosts;  // Compute lorentz transformations
  bool coupled; // is monte carlo evolution coupled to hydro
  bool emission_array_flag;  // Compute and save zone emissivities
  bool polarized;// track photon polarization
  bool acceleration;  // use MRW acceleration
  bool time_acc;  // use MRW acceleration with time limit
  bool raytrace_flag; // Will trace photons rather than scatter
  bool general_mover_flag; // Use integration for photon movement

  // function pointers
  UserMoveFunc_t UserWorkInMove;
  EmisFunc_t InitEmission;
  TempFunc_t GetTemperature;
  ScatFunc_t UserScattering;
  OpacFunc_t UserScatteringOpacity;
  OpacFunc_t UserAbsorptionOpacity;

  // functions
  // SWD: some of these functions could/should be private
  void RunStaticMonteCarlo(Outputs *pouts, Mesh *pmesh, ParameterInput *pinput);
  void RunDynamicMonteCarlo(Outputs *pouts, Mesh *pmesh, ParameterInput *pinput);
  bool CheckAndBroadCastPhotonsRemaining();
  void InitUserMonteCarloData(ParameterInput *pin);
  // Enroll User functions
  void EnrollUserMCBoundaryFunction(enum BoundaryFace dir, MCBValFunc_t my_bc);
  void EnrollUserEmissionInitialization(EmisFunc_t emissfunc);
  void EnrollUserGetTemperature(TempFunc_t tempfunc);
  void EnrollUserWorkInMove(UserMoveFunc_t userfunc);
  void EnrollUserOpacityFunction(OpacFunc_t opacfunc, bool abs);
  void EnrollUserScatteringFunction(ScatFunc_t scatfunc);
  void SendMonteCarloSpectra(int dest);
  void ReceiveMonteCarloSpectra(int source);
  //void CollectMoments(void);
  void Initialize(ParameterInput *pinput);
  void MakeOutputs();

private:

  // functions
  MCBValFunc_t BoundaryFunction_[6];

  void GetDensity(MonteCarloBlock *pmcb);
  void GetVelocity(MonteCarloBlock *pmcb);
  //void SendMonteCarloBlocks(int dest);
  //void ReceiveMonteCarloBlocks(ParameterInput *pin, int source);
  //void SendMonteCarloData(int dest);
  //void ReceiveMonteCarloData(int source);
  //void SendMoments(int dest);
  //void ReceiveMoments(int source, bool sum);

};

//----------------------------------------------------------------------------------------
//! \class MonteCarloBlock
//! \brief monte carlo functions and data contained on each mesh block

class MonteCarloBlock {
public:
  MonteCarloBlock(MeshBlock *pmb, MCBlockSize *pblsize, MonteCarlo *pmc,
                  ParameterInput *pin);
  ~MonteCarloBlock();

  // data
  MonteCarlo* pmy_mc; // MonteCarlo
  MeshBlock* pmy_block;    // MeshBlock corresponding to this MonteCarloBlock
  MonteCarloBlock *next;
  MCCoord *pcoord;

  Photon* pphot; // ptr to photon packet
  PhotonMover* pmover; // ptr to photon mover
  MCRandom *pran; // ptr to random number generator
  MCBoundaryValues *pbval; // ptr to MC boundary values

  Spectrum *pspec; // ptr to spectrum
  PhotonList *pphlist; // ptr to photon list
  PhotonTrajectoryList *ptraj;

  enum MCBoundaryFlag mcb_bcs[6];

  // function pointers
  GetZonePos_t GetZonePosition;
  OpacFunc_t AbsorptionOpacity;
  OpacFunc_t ScatteringOpacity;
  ScatFunc_t Scatter;

  int nphdone; // Photons integrated thus far
  int nphremain; // total number of photons to integrate
  int nabs, nesc, ndes, nscat;
  int nchunk;
  int lid;
  int nx1,nx2,nx3;
  int is,ie,js,je,ks,ke;
  int nfreq, nmu, nphi, nsurf;
  int cadence;

  bool weighted_absorption; // flag controling how absorption is handled
  bool moments_flag; // Compute/output moments
  bool moments_comoving; // Compute in comoving frame
  bool emission_array_flag;  // Compute and save zone emissivities
  bool boosts;  // Compute lorentz transformations
  bool coupled; // Whether time dependent code is coupled to hydro
  bool coherent_scattering; // photon does notchange energy after scattering
  bool acceleration;  // use MRW acceleration
  bool time_acc;  // use MRW acceleration with time limit

  // Set flags
  //enum EmissionFlag emission_meth;
  enum AbsorptionMethodFlag absorption_meth;
  enum AbsorptionOpacityFlag absorption_opac;
  enum ScatteringFlag scattering_meth;

  // Associated with general mover
  // SWD some of these should be eliminated others moved to MonteCarlo?
  bool boyerlindquist_flag; // use Boyer-Lindquist coordinates
  bool orthotet_flag; // use orthonormal tetrad for TransferPhotons()
  bool varystep_flag; // use variable (true) or constant (false) step

  Real codetocgs_rho, codetocgs_vel, codetocgs_tgas;
  Real stepsize;
  Real minweight;

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
  void TransferPhotonsStatic(); // Transfer photons on this block
  void TransferPhotonsDynamic(); // time dependent transfer on this block
  void CoupleMonteCarloToFluid(Real dt);  // couple monte carlo to mesh
  void CoordinateToComoving(Photon *pphot, int ips, int ipe);
  void ComovingToCoordinate(Photon *pphot, int ips, int ipe);
  void LorentzTransform(Photon *pphot, const Real sign, int ips, int ipe);
  Real LorentzTransformFrequencyShift(Photon *pphot, int ip);
  void TetradTransform(Photon *pphot, const Real sign, int ips, int ipe);
  void InitializePhoton(Photon *pphot, int ips, int ipe);
  void FinalizePhoton(Photon *pphot, int ip);
  void UpdateMoments(Photon *pphot, Real dl, Real etau, int ip);
  void NormalizeMoments(bool normalize, Real norm);
  void ResetMoments();
  void UpdateCooling(Photon *pphot, Real energy0, Real weight0, int ip);
  //void GetPhotonsFromNeighbors();
  //void SendPhotonsToNeighbors();

private:
   void SetBoundaryValues(enum MCBoundaryFlag *input_bcs);
};

#endif // MONTECARLO_HPP
