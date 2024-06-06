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
class PhotonPusher;
class MCRandom;
class MCBoundaryValues;
class MCOutoupt;
class MCCoord;

// SWD: Make into a general MACRO set by configure?
// SWD: or make a paramter that is set
#define NMOM 15

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
      MCIPR12=7, MCIPR13=8, MCIPR23=9, MCIPR21=10, MCIPR31=11, MCIPR32=12};
enum SourceTermFlag {MCRS0 = 0, MCRS1=1, MCRS2=2, MCRS3=3, MCRSP0=4, MCRSP1=5,
                     MCRSP2=6, MCRSP3=7};
//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime
typedef Real (*EmisFunc_t)(MonteCarloBlock *pmcb, int k, int j, int i);
typedef void (*TempFunc_t)(MonteCarloBlock *pmcb);
typedef void (*MCBValFunc_t)(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot, int ip);
typedef Real (*OpacFunc_t)(MonteCarloBlock *pmcb,  Photon *pphot, int ip);
typedef void (*ScatFunc_t)(MonteCarloBlock *pmcb, Photon *phot, int ips, int ipe);
typedef void (*UserMoveFunc_t)(MonteCarloBlock *pmcb, Photon *phot, PhotonPusher *ppusher,
                               int ip);
typedef void (*GetZonePos_t)(Photon *phot, MCRandom *pran, MCCoord *pco, int ip);
typedef void (*UserMomentFunc_t)(MonteCarloBlock *pmcb, Photon *phot, Real dl, int ip,
                                 int imom);
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
Real GetEmissionFreeFree(MonteCarloBlock *pmcb, int k, int j, int i);
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

  Real tint;     // Monte Carlo timestep
  Real tmax;   // Maximum evolution time
  Real time_cgs; // conversion of time to cgs units
  Real weightratio; // used for setting minimum weight for absorption
  int64_t nsamp;  // total number of photons to integrate per timestep
  int64_t nphrun;  // number of photons completed
  int nblocal; // number of montecarloblocks on this process
  int nbtotal; // total number of montecarloblocks
  int nout;  // number of outputs
  int64_t ncells; // total number of cells in mesh
  int iseed;  // seed to initialized random number generator(s)

  int list_size_init; // maximum number of photons run per output on any process
  int max_phots_init; // maximum number of photon elements
  int nuser_var, nuser_mom;
  int checkmove,checkscat;

  enum EmissionFlag emission_flag;
  enum MCBoundaryFlag mc_bcs[6];
  enum ScatteringFlag scattering_meth;

  bool dynamic; // is monte carlo evolving with time
  bool coupled; // is monte carlo evolution coupled to hydro
  bool boosts;  // Compute lorentz transformations
  bool tetrads; // convert from coordinate frame
  bool emission_array;  // Compute and save zone emissivities
  bool emission_eqwt; // Set initial weights equal
  bool polarized;// track photon polarization
  bool acceleration;  // use MRW acceleration
  bool computedmin;
  bool time_acc;  // use MRW acceleration with time limit
  bool raytrace_flag; // Will trace photons rather than scatter
  bool general_pusher_flag; // Use integration for photon movement

  // function pointers
  UserMoveFunc_t UserWorkInMove;
  EmisFunc_t GetEmission;
  TempFunc_t UserGetTemperature;
  ScatFunc_t UserScattering;
  OpacFunc_t UserScatteringOpacity;
  OpacFunc_t UserAbsorptionOpacity;

  std::string *user_moment_names;
  UserMomentFunc_t *user_moment_func;

  // functions
  // SWD: some of these functions could/should be private
  void RunStaticMonteCarlo(Outputs *pouts, Mesh *pmesh, ParameterInput *pinput);
  void RunDynamicMonteCarlo(Outputs *pouts, Mesh *pmesh, ParameterInput *pinput);
  bool CheckAndBroadCastPhotonsRemaining();
  void InitUserMonteCarloData(ParameterInput *pin);
  // Enroll User functions
  void EnrollUserMCBoundaryFunction(enum BoundaryFace dir, MCBValFunc_t my_bc);
  void EnrollUserEmissionFunction(EmisFunc_t emissfunc);
  void EnrollUserGetTemperature(TempFunc_t tempfunc);
  void EnrollUserWorkInMove(UserMoveFunc_t userfunc);
  void EnrollUserOpacityFunction(OpacFunc_t opacfunc, bool abs);
  void AllocateUserMoments(int n);
  void EnrollUserMoment(int i, UserMomentFunc_t my_func, const char *name);
  void EnrollUserScatteringFunction(ScatFunc_t scatfunc);
  void Initialize(ParameterInput *pinput);
  void InitializeEmissionFlags(ParameterInput *pinput);
  void ComputeEmission();

private:

  // functions
  MCBValFunc_t BoundaryFunction_[6];

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
  PhotonPusher* ppusher; // ptr to photon pusher
  MCRandom *pran; // ptr to random number generator
  MCBoundaryValues *pbval; // ptr to MC boundary values

  // ouput pointers
  Spectrum *pspec; // ptr to spectrum
  PhotonList *pphlist; // ptr to photon list
  PhotonTrajectoryList *ptraj; // ptr to traj list

  enum MCBoundaryFlag mcb_bcs[6];

  // function pointers
  GetZonePos_t GetZonePosition;
  OpacFunc_t AbsorptionOpacity;
  OpacFunc_t ScatteringOpacity;
  ScatFunc_t Scatter;

  int nphrun; // Photons initialized thus far
  int nphremain; // total number of photons to integrate
  int64_t nabs, nesc, ndes, nscat;
  int loop_max_size;
  int nx1,nx2,nx3;
  int is,ie,js,je,ks,ke;
  int nsrc, nmom; // # of elements in sourcterm, moments arrays

  bool weighted_absorption; // flag controling how absorption is handled
  bool mom_flag_lab; // Compute/output moments
  bool mom_flag_com; // Compute moments in comoving frame
  bool mom_flag_src; // Compute source terms for output
  bool mom_flag_usr; // Compute user defined monte carlo moments
  bool call_moments;
  bool call_srcterms;

  bool boosts;  // Compute lorentz transformations
  bool tetrads; // Compute tetrads
  bool coupled; // Whether time dependent code is coupled to hydro
  bool coherent_scattering; // photon does notchange energy after scattering
  bool acceleration;  // use MRW acceleration
  bool computedmin;
  bool time_acc;  // use MRW acceleration with time limit

  // Set flags
  enum AbsorptionMethodFlag absorption_meth;
  enum AbsorptionOpacityFlag absorption_opac;
  enum ScatteringFlag scattering_meth;

  // Associated with general pusher
  // SWD some of these should be eliminated others moved to MonteCarlo?
  bool boyerlindquist_flag; // use Boyer-Lindquist coordinates
  bool orthotet_flag; // use orthonormal tetrad for TransferPhotons()
  bool varystep_flag; // use variable (true) or constant (false) step

  Real rho_cgs, vel_cgs, tgas_cgs, tfloor_cgs, l_cgs;
  Real betamax;
  Real stepsize;
  Real minweight;
  Real emiss_to_weight; // used relate weight to emission array

  AthenaArray<Real> emission;
  AthenaArray<Real> moments;
  AthenaArray<Real> moments_com;
  AthenaArray<Real> moments_user;
  AthenaArray<Real> sourceterms;
  AthenaArray<Real> scalars;
  AthenaArray<Real> rho;
  AthenaArray<Real> tgas;
  AthenaArray<Real> vel;
  AthenaArray<Real> tran_cmv;
  AthenaArray<Real> tran_crd;
  AthenaArray<Real> planck_opacity; // for acceleration
  AthenaArray<Real> planck_inv_opacity; // for acceleration

  // functions
  void InitUserMonteCarloBlockData(ParameterInput *pin);
  void MonteCarloProblemGenerator(ParameterInput *pin);
  void RayTracePhotonsOnBlock(); // Ray trace photon on this block
  void TransferPhotonsOnBlock(); // Transfer photons on this block
  void CoupleMonteCarloToFluid(Real dt);  // couple monte carlo to mesh
  void LorentzTransform(Photon *pphot, const Real sign, int ips, int ipe);
  Real LorentzTransformFrequencyShift(Photon *pphot, int ip);
  void TetradTransform(Photon *pphot, const Real sign, int ips, int ipe);
  void InitializePhoton(Photon *pphot, int ips, int ipe);
  void FinalizePhoton(Photon *pphot, int ip);
  void UpdateMoments(Photon *pphot, Real dl, Real etau, int ip);
  void UpdateMoments(Photon *pphot, Real dl, int ip);
  void UpdateMomentsAcceleration(Photon *pphot, Real dl, Real pl, Real k1, Real k2,
                                 Real k3,Real etau, int ip);
  void UpdateMomentsOld(Photon *pphot, Real dl, Real pl, Real k1, Real k2, Real k3,
                        Real etau, int ip);
  void NormalizeMoments(bool normalize);
  void ResetMoments();
  void UpdateSourceTerms(Photon *pphot, Real energy0, Real weight0,
                         Real k1p0, Real k2p0, Real k3p0, int ip);
  void NormalizeSourceTerms(bool normalize);
  void ResetSourceTerms();
  // Functions for handling distributed emission over cells
  void ComputeEmissionArray(Real &emm_min, Real &emm_max, Real &emm_tot);
  void SetEmissionCellWeight(Photon *pphot, int ips, int ipe);
  void GetDensity();
  void GetScalars();
  void GetVelocity();
  void GetTemperature();
  void ComputeTransformations();
  void TransformToComoving(Photon *pphot, int ips, int ipe);
  void TransformToCoordinate(Photon *pphot, int ips, int ipe);
  Real FrequencyShiftComoving(Photon *pphot, int ips);
  void  FrequencyAngleShiftComoving(Photon *pphot, int ip, Real &shift,
                                    Real &k1, Real &k2, Real &k3);
  Real FrequencyShiftCoordinate(Photon *pphot, int ips);

private:
  int i1_, i2_, i3_; // used for emission
  Real nemit_; // used for emission

  void SetBoundaryValues(enum MCBoundaryFlag *input_bcs);
};

#endif // MONTECARLO_HPP
