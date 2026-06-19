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
#include <random>
// Athena++ classes headers
#include "../athena.hpp"
#include "../coordinates/coordinates.hpp"
#include "../outputs/outputs.hpp"
#include "../field/field.hpp"
#include "photon.hpp"
#include "mcbvals.hpp"
#include "mcoutput.hpp"
#include "mccoord.hpp"

// GSL library
#if GSL
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
enum EmissionFlag {EMISUSER = 0, EMISNONE = 1, EMISFF = 2, EMISBB = 3, MULTI = 4};
enum EmissionGeometry {EMISVOL = 0, EMISAREA = 1, EMISGNONE = 2};
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
enum SourceTermFlag {MCRS0 = 0, MCRS1=1, MCRS2=2, MCRS3=3, MCRF0=4, MCRF1=5,
                     MCRF2=6, MCRF3=7, MCNABS=8};
//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime
typedef Real (*EmisFunc_t)(MonteCarloBlock *pmcb, int k, int j, int i, int etype);
typedef void (*DensFunc_t)(MonteCarloBlock *pmcb);
typedef void (*TempFunc_t)(MonteCarloBlock *pmcb);
typedef void (*NumbFunc_t)(MonteCarloBlock *pmcb);
typedef void (*MCBValFunc_t)(MonteCarloBlock *pmcb, MCCoord *pco, Photon *pphot, int ip);
typedef Real (*OpacFunc_t)(MonteCarloBlock *pmcb,  Photon *pphot, int ip);
typedef void (*ScatFunc_t)(MonteCarloBlock *pmcb, Photon *phot, int ips, int ipe);
typedef void (*UserMoveFunc_t)(MonteCarloBlock *pmcb, Photon *phot, PhotonPusher *ppusher,
                               int ip);
typedef void (*GetZonePos_t)(Photon *phot, MCRandom *pran, MCCoord *pco, int ip);
typedef void (*UserMomentFunc_t)(MonteCarloBlock *pmcb, Photon *phot, Real dl, int ip,
                                 int imom);
typedef void (*UserSourcetermFunc_t)(MonteCarloBlock *pmcb, Photon *pphot,
                                    Real energy0, Real weight0, Real k1p0,
                                    Real k2p0, Real k3p0, int ip);
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
Real GetEmissionFreeFree(MonteCarloBlock *pmcb, int k, int j, int i, int etype);
void PhotonEmitFreeFree(MonteCarloBlock *pmcb, Photon *pphot, Real lemin, Real lemax,
                        int ip);
Real GetEmissionBlackbody(MonteCarloBlock *pmcb, int k, int j, int i, int etype);
void PhotonEmitBlackbody(MonteCarloBlock *pmcb, Photon *pphot, BoundaryFace face, int ip);
Real PlanckDist(Real temp, MCRandom *pran);
void GetZonePositionCartesian(Photon *pphot, MCRandom *pran, MCCoord *pco, int ip);
void GetZonePositionSphericalPolar(Photon *pphot, MCRandom *pran, MCCoord *pco, int ip);
void GetZonePositionCylindrical(Photon *pphot, MCRandom *pran, MCCoord *pco, int ip);
void GetZonePositionCartesianFace(Photon *pphot, MCRandom *pran, MCCoord *pcoord,
                                  BoundaryFace face, int ip);
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
void LorentzBoostVector(Real vel[4], Real kold[4]);

//---------------------- prototypes for setting flags ------------------------------------
enum MCBoundaryFlag GetMCBoundaryFlag(std::string input_string);
enum EmissionFlag GetEmissionFlag(std::string input_string);
enum EmissionGeometry GetEmissionGeometry(std::string input_string);
enum BoundaryFace SetEmissionSurface(std::string input_face);
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
  Real chisquare(Real nu);
  int binomial(unsigned int n, Real p);
  void SampleMultinomial(int n, int m, Real *prob, int *counts);

private:

#if GSL
  gsl_rng *dev;
#endif
  std::mt19937 gen;
  std::uniform_real_distribution<Real> uniform_dist;

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

  Real tint;   // Monte Carlo timestep
  Real tmax;   // Maximum evolution time
  Real time_cgs; // conversion of time to cgs units
  Real weightratio; // used for setting minimum weight for absorption

  int ntype; // number of emission types 
  int64_t nsamp;  // total number of photons to integrate per timestep/output
  int64_t *nsamptype; // number of sample per type
  int nblocal; // number of montecarloblocks on this process
  int nbtotal; // total number of montecarloblocks
  int nout;  // number of outputs
  int64_t ncells; // total number of cells in mesh
  int iseed;  // seed to initialized random number generator(s)
 
  int list_size_init; // maximum number of photons run per output on any process
  int max_phots_init; // maximum number of photon elements
  int nuser_var, nuser_mom;
  int checkmove,checkscat;
  int emission_method;
  int *emission_geometry;
  BoundaryFace *emission_face;

  enum EmissionFlag emission_flag;

  enum MCBoundaryFlag mc_bcs[6];
  enum ScatteringFlag scattering_meth;

  bool dynamic; // is monte carlo evolving with time
  bool coupled; // is monte carlo evolution coupled to hydro
  bool boosts;  // Compute lorentz transformations
  bool using_bfield; // set magnetic fields
  bool tetrads; // convert from coordinate frame
  bool emission_array;  // Compute and save zone emissivities
  bool *emission_eqwt; // Set initial weights equal
  
  bool polarized;// track photon polarization
  bool acceleration;  // use MRW acceleration
  bool computedmin;
  bool time_acc;  // use MRW acceleration with time limit
  bool raytrace_flag; // Will trace photons rather than scatter
  bool general_pusher_flag; // Use integration for photon movement
  bool verbose; // print out more information during run

  // function pointers
  UserMoveFunc_t UserWorkInMove;
  EmisFunc_t *GetEmission; // array of function pointers
  DensFunc_t UserGetDensity;
  TempFunc_t UserGetTemperature;
  NumbFunc_t UserGetNumberDensity;
  ScatFunc_t UserScattering;
  OpacFunc_t UserScatteringOpacity;
  OpacFunc_t UserAbsorptionOpacity;
  std::string *user_moment_names;
  UserMomentFunc_t *user_moment_func;
  UserSourcetermFunc_t UserSourcetermFunc;

  // functions
  // SWD: some of these functions could/should be private
  void RunMonteCarlo(Outputs *pouts, Mesh *pmesh, ParameterInput *pinput);
  bool CheckAndBroadCastPhotonsRemaining();
  void InitUserMonteCarloData(ParameterInput *pin);
  // Enroll User functions
  void EnrollUserMCBoundaryFunction(enum BoundaryFace dir, MCBValFunc_t my_bc);
  void EnrollUserEmissionFunction(EmisFunc_t emissfunc);
  void EnrollUserEmissionFunction(EmisFunc_t emissfunc, int etype);
  void EnrollUserGetDensity(DensFunc_t densfunc);
  void EnrollUserGetTemperature(TempFunc_t tempfunc);
  void EnrollUserGetNumberDensity(NumbFunc_t numbunc);
  void EnrollUserWorkInMove(UserMoveFunc_t userfunc);
  void EnrollUserOpacityFunction(OpacFunc_t opacfunc, bool abs);
  void AllocateUserMoments(int n);
  void EnrollUserMoment(int i, UserMomentFunc_t my_func, const char *name);
  void EnrollUserSourcetermUpdate(UserSourcetermFunc_t my_func);
  void EnrollUserScatteringFunction(ScatFunc_t scatfunc);
  void Initialize(ParameterInput *pinput);
  void InitializeEmission(ParameterInput *pin);
  void DistributeSamples(int etype);
  void NormalizeDomainOutputs(bool normalize);
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

  int64_t nphrun; // Photons initialized thus far
  int64_t nphremain; // total number of photons to integrate
  int64_t nabs, nesc, ndes, nscat, nrem; // counters
  int loop_max_size;
  int nx1,nx2,nx3;
  int is,ie,js,je,ks,ke;
  int nsrc, nmom; // # of elements in sourcterm, moments arrays
  int nf_scat;

  bool weighted_absorption; // flag controling how absorption is handled
  bool mom_flag_lab; // Compute/output moments
  bool mom_flag_com; // Compute moments in comoving frame
  bool mom_flag_src; // Compute source terms for output
  bool mom_flag_usr; // Compute user defined monte carlo moments
  bool mom_flag_scat; // Compute scattering source terms
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

  Real rho_cgs, vel_cgs, tgas_cgs, tfloor_cgs, tceiling_cgs, l_cgs;
  Real betamax;
  Real stepsize;
  Real minweight;
  Real emiss_to_weight; // used relate weight to emission array
  Real emin_scat, emax_scat, dloge_scat; // min/max energy for scattering moments

  AthenaArray<Real> emission;
  AthenaArray<Real> moments;
  AthenaArray<Real> moments_com;
  AthenaArray<Real> moments_user;
  AthenaArray<Real> moments_scat;
  AthenaArray<Real> energy_scat;
  AthenaArray<Real> freq_scat_mid;
  AthenaArray<Real> sourceterms;
  AthenaArray<Real> scalars;
  AthenaArray<Real> rho;
  AthenaArray<Real> nel;
  AthenaArray<Real> nion;
  AthenaArray<Real> tgas;
  AthenaArray<Real> vel;
  AthenaArray<Real> bcc;
  AthenaArray<Real> boost_cmv;
  AthenaArray<Real> boost_lab;
  AthenaArray<Real> planck_opacity; // for acceleration
  AthenaArray<Real> planck_inv_opacity; // for acceleration

  // functions
  void InitUserMonteCarloBlockData(ParameterInput *pin);
  void MonteCarloProblemGenerator(ParameterInput *pin);
  void RayTracePhotonsOnBlock(int etype); // Ray trace photon on this block
  void TransferPhotonsOnBlock(int etype); // Transfer photons on this block
  void CoupleMonteCarloToFluid(Real dt);  // couple monte carlo to mesh
  void LorentzTransform(Photon *pphot, const Real sign, int ips, int ipe);
  Real LorentzTransformFrequencyShift(Photon *pphot, int ip);
  void TetradTransform(Photon *pphot, const Real sign, int ips, int ipe);
  void InitializePhoton(Photon *pphot, int ips, int ipe, int etype);
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
  void ComputeEmissionArray(int etype, Real &emm_min, Real &emm_max, Real &emm_tot);
  void ComputeEmissionSampleArray();
  //void ComputeEmissionSampleArray(BoundaryFace face);
  void SetEmissionCellWeight(Photon *pphot, int ips, int ipe);
  void SetEmissionCellWeightArea(Photon *pphot, BoundaryFace face, int ips, int ipe);
  void GetDensity();
  void GetNumberDensity();
  void GetScalars();
  void GetVelocity();
  void GetBField();
  void GetTemperature();
  void ComputeTransformations();
  void TransformToComoving(Photon *pphot, int ips, int ipe);
  void TransformToCoordinate(Photon *pphot, int ips, int ipe);
  Real FrequencyShiftComoving(Photon *pphot, int ips);
  void  FrequencyAngleShiftComoving(Photon *pphot, int ip, Real &shift,
                                    Real &k1, Real &k2, Real &k3);
  Real FrequencyShiftCoordinate(Photon *pphot, int ips);
  void UserWorkAfterTransfer(int etype);

private:
  int i1_, i2_, i3_; // used for emission
  AthenaArray<int> emit_count_; // used for emission
  void SetBoundaryValues(enum MCBoundaryFlag *input_bcs);
};

#endif // MONTECARLO_HPP
