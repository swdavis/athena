//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code
// contributors Licensed under the 3-clause BSD License, see LICENSE file for
// details
//========================================================================================
//! \file mc_hotjupiter.cpp
//! \brief Plane-parallel source of radiation from a star incident on the atomic H layer
//! of an exoplanet atmosphere, including photoelectric heating, recombination, and
//! resonance scattering of Lya

// C headers

// C++ headers
#include <cmath>
#include <cstdio>  // fopen(), fprintf(), freopen()
#include <cstring> // strcmp()
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../monte_carlo/montecarlo.hpp"
#include "../monte_carlo/photon.hpp"
#include "../monte_carlo/photonpusher.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../inputs/hdf5_reader.hpp"

#if NSCALARS < 1
#error "This problem generator requires scalars to track neutral hydrogen fraction"
#endif

namespace {

enum emission_type {ION_STR= 0, LYA_STR = 1, LYA_REC = 2};

// GLOBAL VARIABLES
// ----------------

Real kb_cgs, mp_cgs, c_cgs, h_cgs;
Real rin, rout; // inner and outer limits of coordinate system
Real energy0, nu0, threshold;
Real gm_planet, gm_star, sep, psi;
Real temp0, nbase;
Real lya_flux;
Real ion_flux;
Real edot_lya;
Real edot_ion;
Real linewidth, linewidth_cutoff, linewidth_cutoff_energy;
Real stddev;
Real numin, numax, powa, mean_nu, chromo_temp, sigmamin;
Real numin_erg;
Real ev_to_erg;
Real user_dt;
Real mdot;
static int iset = 0; // Flag for performance-saving random deviate sampling (internal use only - DO NOT SET)
static Real gset; // Variable that contains extra deviate from sampling using Box-Muller method
Real pfloor;
Real dfloor;

// Ghost zones
static Real rho_gz[2];
static Real P_gz[2];
static Real nfrac_gz[2];

bool flag_incident_from_z;
bool flag_zero_opacity;
bool flag_escape_after_scatter;
bool flag_lya_volume_emis, flag_lya_surface_emis, flag_ion_surface_emis;
bool flag_point_mass = false;
int  flag_tidal_gravity;
bool flag_dynamic = false;
bool flag_wind = false;
bool flag_pow_law = false;

enum PhotonType {NO_TYPE=0, LYA=1, IONIZING=2};
enum EmissionGeometryLya {NO_GEOMETRY=0, VOLUME=1, SURFACE=2};

// Function for sampling photon frequencies
void gasdev(MeshBlock *pmb, Real mean, Real stddev, Real &samp);

void ForceEscape(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);

Real GetIsowindVelocity(Real x);

// Function to calculate the emission array for Lya emitted from recombinations
// and from electron impact excitation
//Real MultipleEmissivities(MonteCarloBlock *pmcb, int k, int j, int i);
Real VolumeEmissivityLya(MonteCarloBlock *pmcb, int k, int j, int i, int etype);
Real SurfaceEmissivityLya(MonteCarloBlock *pmcb, int k, int j, int i, int etype);
Real SurfaceEmissivityIonizing(MonteCarloBlock *pmcb, int k, int j, int i, int etype);

void F1pos(MonteCarloBlock *pmcb, Photon *pphot, Real dl, int ip, int imom);
void F1neg(MonteCarloBlock *pmcb, Photon *pphot, Real dl, int ip, int imom);
void a1pos(MonteCarloBlock *pmcb, Photon *pphot, Real dl, int ip, int imom);
void a1neg(MonteCarloBlock *pmcb, Photon *pphot, Real dl, int ip, int imom);

// Tidal gravity source term
void GetTidalAcceleration(Real r, Real th, Real ph, Real rho, Real &a_r, Real &a_th, Real &a_ph);
void TwoPointMass(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);
void HillTidalGravity(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);
void ThirdOrderTidalGravity(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);

// Explicit heating function, test
void ExplicitEUVHeating(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);


// INITIALIZATION FUNCTIONS
// ========================

void ResonantScattering(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe);
Real ResonantScatteringOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip);
Real BoundFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip);


Real ConstantTimestep(MeshBlock *pmb) {
  if (user_dt > TINY_NUMBER) {
    return user_dt;
  } else {
    return std::numeric_limits<Real>().max();
  }
}
void GetIonizationTemperature(MonteCarloBlock *pmcb);
//void EscapeCoords(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover, int ip);

void UpdateSourceTerms(MonteCarloBlock *pmcb, Photon *pphot,Real energy0, Real weight0,
                  Real k1p0, Real k2p0, Real k3p0, int ip);

} // end namespace

// SWD: Why aren't these inside anonymous namespace -- move

// Boundary Conditions
void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void StaticInflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh);



// SETUP
// =====

void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Initialize file-scope variables with values read in from the input file
  // Enroll the user gravity source term
  h_cgs = 6.62607015e-27;
  kb_cgs = 1.380649e-16;
  mp_cgs = 1.660538782e-24;
  c_cgs = 2.99792458e10;
  psi = pin->GetOrAddReal("problem", "psi", 0.0)*PI; // radians, defaults to star in -x direction

  // Parameters for Lya emission
  nu0 = 2.4660675e+15;
  energy0 = h_cgs * nu0;

  // H-Ionizing threshold energy at 912 angstrom (13.6 eV)
  threshold = h_cgs*c_cgs/9.12e-6;

  // Gravity source term configuration
  gm_planet = pin->GetReal("problem", "GM");
  gm_star = pin->GetReal("problem", "gm_star");
  sep = pin->GetReal("problem", "orbit_sep")*1.49597871e+13; // cm
  flag_tidal_gravity = pin->GetOrAddInteger("problem", "tidal_gravity", 0);
  if (flag_tidal_gravity == 1) {
    EnrollUserExplicitSourceFunction(HillTidalGravity);
  }
  if (flag_tidal_gravity == 2) {
    EnrollUserExplicitSourceFunction(ThirdOrderTidalGravity);
  }

  // Parameters to set initial atmospheric conditions
  temp0 = pin->GetOrAddReal("problem", "isothermal_temp", 1e4);
  nbase = pin->GetOrAddReal("problem", "nbase", 1.0e10); // base hydrogen nH density in cm^-3 at rin
  flag_wind = pin->GetBoolean("problem", "wind");


  // Sources of photons
  flag_incident_from_z = pin->GetBoolean("problem", "incident_from_z");
  flag_lya_volume_emis = pin->GetBoolean("problem", "lya_volume_emis");
  flag_lya_surface_emis = pin->GetBoolean("problem", "lya_surface_emis");
  flag_ion_surface_emis = pin->GetBoolean("problem", "ion_surface_emis");
  lya_flux = pin->GetReal("problem", "lya_flux");
  ion_flux = pin->GetReal("problem", "ion_flux");
  edot_lya = PI * SQR(mesh_size.x1max) * lya_flux;
  edot_ion = PI * SQR(mesh_size.x1max) * ion_flux;
  linewidth = pin->GetReal("problem", "linewidth");
  linewidth_cutoff = pin->GetReal("problem", "linewidth_cutoff");

  // Debug flag to turn off opacity entirely
  flag_zero_opacity = pin->GetBoolean("problem", "zero_opacity");
  flag_escape_after_scatter = pin->GetBoolean("problem", "escape_after_scatter");

  // Normal distribution for sampling initial photon frequencies
  // Using a standard deviation of 67 km/s
  stddev = linewidth*1.0e5 / c_cgs * energy0;
  linewidth_cutoff_energy = linewidth_cutoff*1.0e5 / c_cgs * energy0;
  sigmamin = 6.3e-18; // absorption edge cross section

  // Time evolution
  // Set a constant timestep (best for photoionization equilibrium convergence tests)
  flag_dynamic = pin->GetBoolean("montecarlo", "dynamic");
  user_dt = pin->GetReal("problem", "user_dt");
  EnrollUserTimeStepFunction(ConstantTimestep);

  // Boundary conditions
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, StaticInflowInnerX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OutflowOuterX1);

  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  int nx1 = pmy_mesh->mesh_size.nx1+2*NGHOST;
  int nx2 = pmy_mesh->mesh_size.nx2+2*NGHOST;
  int nx3 = pmy_mesh->mesh_size.nx3+2*NGHOST;

  AllocateRealUserMeshBlockDataField(3);

  // Array for photon emissivities
  ruser_meshblock_data[0].NewAthenaArray(4,nx3,nx2,nx1);

  // Array for each component of the tidal gravity acceleration
  ruser_meshblock_data[1].NewAthenaArray(3,nx3,nx2,nx1);

  // Array for data loaded in from external file
  ruser_meshblock_data[2].NewAthenaArray(6,ncells3,ncells2,ncells1);

  AllocateUserOutputVariables(6);
  SetUserOutputVariableName(0, "vol_emis");
  SetUserOutputVariableName(1, "surf_emis");
  SetUserOutputVariableName(2, "gravsrc_r");
  SetUserOutputVariableName(3, "gravsrc_th");
  SetUserOutputVariableName(4, "gravsrc_ph");
  SetUserOutputVariableName(5, "collisional_cooling");
  return;
}

// set user output variables to corresponding meshblock arrays
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  AthenaArray<Real> &uov = user_out_var;
  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++) {
        // arrays to store tidal gravity components
        uov(2,k,j,i) = ruser_meshblock_data[1](0,k,j,i);
        uov(3,k,j,i) = ruser_meshblock_data[1](1,k,j,i);
        uov(4,k,j,i) = ruser_meshblock_data[1](2,k,j,i);
      }
    }
  }
}

// Read in input params that define the attributes of the Monte Carlo photons
void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin) {

  //AllocateUserMoments(4);
  //EnrollUserMoment(0, F1pos, "F1pos");
  //EnrollUserMoment(1, F1neg, "F1neg");
  //EnrollUserMoment(2, a1pos, "a1pos");
  //EnrollUserMoment(3, a1neg, "a1neg");

  // set up emission methods
  nsamp = nsamptype[0] = pin->GetInteger("problem", "nion");
  emission_eqwt[0] = pin->GetOrAddBoolean("problem", "ion_eqwt",false);
  EnrollUserEmissionFunction(SurfaceEmissivityIonizing,0);
  emission_geometry[0] = EMISAREA;
  emission_face[0] = SetEmissionSurface("outer_x1");
  if (ntype > 1) {
    nsamp += nsamptype[1] = pin->GetInteger("problem", "nlyastr");
    emission_eqwt[1] = pin->GetOrAddBoolean("problem", "lys_eqwt",false);
    emission_geometry[1] = EMISAREA;
    EnrollUserEmissionFunction(SurfaceEmissivityLya,1);
    emission_face[1] = SetEmissionSurface("outer_x1");
  }
  if (ntype == 3) {
    nsamp += nsamptype[2] = pin->GetInteger("problem", "nlyarec");
    emission_eqwt[2] = pin->GetOrAddBoolean("problem", "lyr_eqwt",false);
    emission_geometry[2] = EMISVOL;
    EnrollUserEmissionFunction(VolumeEmissivityLya,2);
    emission_face[2] = SetEmissionSurface("none");
  }

  EnrollUserSourcetermUpdate(UpdateSourceTerms);
  // Photon user variables for incident and outgoing direction vector and position
  nuser_var = 1;

  //EnrollUserWorkInMove(EscapeCoords);
  EnrollUserOpacityFunction(ResonantScatteringOpacity, false);
  EnrollUserOpacityFunction(BoundFreeAbsorptionOpacity, true);
  EnrollUserScatteringFunction(ResonantScattering);
  //EnrollUserEmissionFunction(MultipleEmissivities);
  EnrollUserGetTemperature(GetIonizationTemperature);
}

// INITIALIZATION
// ==============

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Set initial conditions before the main loop starts. Note that any variables
  // or arrays initialized here are *not* saved in restart files and are
  // inaccessible except at initialization

  // If initializing passive scalar from an existing file, load it here
  bool flag_initialize_scalar_from_file = pin->GetBoolean("problem", "initialize_scalar_from_file");
  bool flag_initialize_pressure_from_file = pin->GetBoolean("problem", "initialize_pressure_from_file");
  bool flag_initialize_velocity_from_file = pin->GetBoolean("problem", "initialize_velocity_from_file");
  bool flag_initialize_density_from_file = pin->GetBoolean("problem", "initialize_density_from_file");

  // BEGIN INITIALIZATION FILE I/O
  if ((flag_initialize_scalar_from_file || flag_initialize_pressure_from_file)) {
    std::string input_filename = pin->GetString("problem", "input_filename");
    std::string dataset_prim = pin->GetString("problem", "dataset_prim");
    int start_prim_file[5];
    start_prim_file[1] = gid;
    start_prim_file[2] = 0;
    start_prim_file[3] = 0;
    start_prim_file[4] = 0;
    int start_prim_indices[6];
    start_prim_indices[0] = 0;
    start_prim_indices[1] = 1;
    start_prim_indices[2] = 2;
    start_prim_indices[3] = 3;
    start_prim_indices[4] = 4;
    start_prim_indices[5] = 5;
    int count_prim_file[5];
    count_prim_file[0] = 1;
    count_prim_file[1] = 1;
    count_prim_file[2] = block_size.nx3;
    count_prim_file[3] = block_size.nx2;
    count_prim_file[4] = block_size.nx1;
    int start_prim_mem[4];
    start_prim_mem[1] = ks;
    start_prim_mem[2] = js;
    start_prim_mem[3] = is;
    int count_prim_mem[4];
    count_prim_mem[0] = 1;
    count_prim_mem[1] = block_size.nx3;
    count_prim_mem[2] = block_size.nx2;
    count_prim_mem[3] = block_size.nx1;


    for (int n=0; n<6; ++n) {
      start_prim_file[0] = start_prim_indices[n];
      start_prim_mem[0] = n;
      HDF5ReadRealArray(input_filename.c_str(), dataset_prim.c_str(), 5, start_prim_file,
                        count_prim_file, 4, start_prim_mem,
                        count_prim_mem, ruser_meshblock_data[2], true);
    }
  }
  // END INITIALIZATION FILE I/O


  rin = pmy_mesh->mesh_size.x1min;       // base radius
  rout = pmy_mesh->mesh_size.x1max;

  Real mmw0 = 1.0;                      // mean molecular weight at base
  Real a2 = kb_cgs * temp0 / mmw0 / mp_cgs;     // sound speed squared at base
  Real a = std::sqrt(a2);
  Real lambda = gm_planet / (rin * a2);
  Real rsonic = gm_planet / (2.0 * a2);
  //printf("lambda: %g\n", lambda);
  //printf("gm_planet, rmin, asound: %g %g %g\n", gm_planet, rin, a);

  Real gamma = peos->GetGamma();
  Real invgm1 = 1.0/(gamma - 1.0);

  Real Gamma0 = 4.e-5; // 1/s, photoionization rate coefficient
  Real alpha = 2.6e-13;
  Real neq0 = Gamma0 / alpha;
  Real sigmapi = 6.e-18;

  Real vbase = 0.0;
  Real vel1 = 0.0;
  Real vel2 = 0.0;
  Real vel3 = 0.0;
  if (flag_wind) {
    vbase = a * GetIsowindVelocity(2.0/lambda);
    mdot = 4.*PI * SQR(rin) * nbase * mp_cgs * vbase;
  }

  //for (int n=0; n<6; ++n) {
  //  printf("Input data (%d, 4, 4, 4): %g\n", n, ruser_meshblock_data[2](n,4,4,4));
  //}
  //printf("\n");

  // density and pressure floors
  Real float_min = std::numeric_limits<float>::min();
  dfloor = pin->GetOrAddReal("hydro", "dfloor", (1024*(float_min)));
  pfloor = pin->GetOrAddReal("hydro", "pfloor", (1024*(float_min)));

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is-NGHOST; i<=ie+NGHOST; i++) {

        Real r = pcoord->x1v(i);
        Real th = pcoord->x2v(j);
        Real ph = pcoord->x3v(k);
        Real H = rin / lambda * (SQR(r) / SQR(rin));
        Real nH = nbase * std::exp(lambda * (rin/r - 1.));
        Real ion_rate_atten = Gamma0 / (1.0 + std::pow(nH*sigmapi*H, 1.5));
        Real np = std::sqrt(ion_rate_atten * nH / alpha);
        Real neutral_frac = nH / (nH + 2.*np);
        Real mean_mol_weight = (1. + neutral_frac)/2.;
        // ^ This will be wrong for the isothermal wind, but there is no analytic solution

        Real rho;
        if (flag_wind) {
          vel1 = a * GetIsowindVelocity(r/rsonic);
          rho = mdot / (4.0 * PI * r * r * vel1);
          // Re-assign neutral and ion densities based on new density profile
          // Doing it this way so that nH + np = rho
          // If we kept the nH from HSE, we may have places where nH > new rho
          nH = neutral_frac * rho / mean_mol_weight / mp_cgs;
          np = (1.-mean_mol_weight) * rho / mean_mol_weight / mp_cgs;
        } else {
          vel1 = 0.0;
          rho = mp_cgs * (nH + np);
        }

        Real mmw = (nH + np)/(2.*np + nH);
        Real a2_mmw = kb_cgs * temp0 / mmw / mp_cgs;
        //Real P = rho * a2;
        Real P = rho * a2_mmw;

        // Set global ghost zone constant values for inner boundary
        if ((r <= rin) && (i < NGHOST)) {
          rho_gz[i] = rho;
          P_gz[i] = P;
        }

        if (flag_initialize_pressure_from_file) {
          P = ruser_meshblock_data[2](1,k,j,i);
        }

        if (flag_initialize_density_from_file) {
          rho = ruser_meshblock_data[2](0,k,j,i);
        }

        if (flag_initialize_velocity_from_file) {
          //printf("INITIALIZED FROM FILE: %g     INITIALIZED FROM PROBLEM: %g\n", ruser_meshblock_data[2](2,k,j,i), vel1);
          vel1 = ruser_meshblock_data[2](2,k,j,i);
          vel2 = ruser_meshblock_data[2](3,k,j,i);
          vel3 = ruser_meshblock_data[2](4,k,j,i);
        } else {
          vel2 = 0.0;
          vel3 = 0.0;
        }

        // Set corresponding conserved fluid variables for the above
        phydro->u(IDN,k,j,i) = rho;             // Density
        phydro->u(IM1,k,j,i) = rho * vel1;             // X1 Momentum
        phydro->u(IM2,k,j,i) = rho * vel2;             // X2 Momentum
        phydro->u(IM3,k,j,i) = rho * vel3;             // X3 Momentum
        phydro->u(IEN,k,j,i) = P*invgm1;        // Internal energy
        phydro->u(IEN,k,j,i) += 0.5 * SQR(phydro->u(IM1,k,j,i)) / phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5 * SQR(phydro->u(IM2,k,j,i)) / phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5 * SQR(phydro->u(IM3,k,j,i)) / phydro->u(IDN,k,j,i);
        //printf("radius, density, velocity, pressure: %e, %e, %e, %e\n", r, rho, vel1, P);

        //phydro->w(IDN,k,j,i) = rho;
        //phydro->w(IVX,k,j,i) = vel1;
        //phydro->w(IVY,k,j,i) = vel2;
        //phydro->w(IVZ,k,j,i) = vel3;
        //phydro->w(IPR,k,j,i) = P;

        // Ionization state of the gas
        if (flag_initialize_scalar_from_file) {
          //printf("INITIALIZED FROM FILE: %g     INITIALIZED FROM PROBLEM: %g\n", ruser_meshblock_data[2](0,k,j,i), rho);
          pscalars->s(0,k,j,i) = ruser_meshblock_data[2](5,k,j,i)*rho;
          //pscalars->r(0,k,j,i) = ruser_meshblock_data[2](5,k,j,i); // concentration (0 <= r <= 1)

        } else {
          // Initialize presuming optically-thin photoionization rate equilibrium
					// This passive scalar tracks mass-fraction: s = rho_H, r = rho_H / rho = n_H / (n_H + n_p)
          pscalars->s(0,k,j,i) = nH * mp_cgs;
          //pscalars->r(0,k,j,i) = nH / (nH + np); // concentration (0 <= r <= 1)
        }
				//printf("ionization fraction: %e\n", pscalars->r(0,k,j,i));

        if (!flag_dynamic) {
          Real a_r, a_th, a_ph;
          GetTidalAcceleration(r, th, ph, rho, a_r, a_th, a_ph);

          // Update the user output quantities
          ruser_meshblock_data[1](0,k,j,i) = a_r;
          ruser_meshblock_data[1](1,k,j,i) = a_th;
          ruser_meshblock_data[1](2,k,j,i) = a_ph;
        }
      }

      //for (int i=ie+NGHOST; i>=is-NGHOST; i--) {
      //  if (((k==ks) && (j==js)) && i>=is) {
      //    Real tgas = temp0;//pmy_mcb->tgas(k,j,i);
      //    Real chi = pscalars->s(0,k,j,i)/mp_cgs * XsecVoigt(energy0/h_cgs, tgas);
      //    lc_tau += chi * pcoord->dx1f(i);
      //    printf("i=%d   nh=%g   xsec=%g   dx1f=%g   lc_tau = %g\n", i, pscalars->s(0,k,j,i)/mp_cgs, XsecVoigt(energy0/h_cgs, tgas), pcoord->dx1f(i), lc_tau);
      //  }
      //}
    }
  }
  //printf("LINE CENTER OPTICAL DEPTH THROUGH DOMAIN: %g\n", lc_tau);
}


void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {
  numin = pin->GetReal("problem", "numin");
  numin_erg = numin * h_cgs; // absorption edge energy
  numax = pin->GetReal("problem", "numax");
  powa = pin->GetReal("problem", "powa");
	if (powa <= 2.) {
		mean_nu = std::log(numax/numin) / (std::pow(numax, 1-powa) - std::pow(numin, 1-powa)) * (1-powa);
	} else {
		mean_nu = (std::pow(numax, 2-powa) - std::pow(numin, 2-powa)) * (1-powa)
			/ (std::pow(numax, 1-powa) - std::pow(numin, 1-powa)) / (2-powa);
	}
	//printf("mean freq: %g\n", mean_nu);
	flag_pow_law = pin->GetBoolean("problem", "pow_law");
	chromo_temp = pin->GetReal("problem", "chromo_temp");
}

void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe, int etype) {
  // Function called each time a photon is initialized

  std::stringstream msg;
  AthenaArray<Real> emis = pmy_block->ruser_meshblock_data[0];
  PhotonType phot_type = NO_TYPE;
  EmissionGeometryLya emis_geometry = NO_GEOMETRY;

  // 1. Determine from which cell the photon is emitted
  // --------------------------------------------------
  // Sets initial cells and emission weights for all photon samples.
  // This function is unmodified from the single photon population case since
  // it doesn't need to know which emissivity the photon came from, only the
  // total number

  // Note: on the RunPhotoionization pass, the emissivity array will contain
  // only the emissivity from the ionizing radiation source
  if (etype == ION_STR) {
    SetEmissionCellWeightArea(pphot,pmy_mc->emission_face[0],ips,ipe);
  } else if (etype == LYA_STR) {
    SetEmissionCellWeightArea(pphot,pmy_mc->emission_face[1],ips,ipe);
  } else if (etype == LYA_REC) {
    SetEmissionCellWeight(pphot,ips,ipe);
  }


  // Calculate numerical prefactors for sampling so we don't repeat math operations
  Real numinpow = std::pow(numin, 1.0-powa);
  Real numaxpow = std::pow(numax, 1.0-powa);
  Real nuexp = 1.0 / (1.0 - powa);

  // Loop over photon index
  for (int ip=ips; ip<=ipe; ip++) {

    pphot->dtp[ip] = pphot->pmy_mcb->pmy_mc->tmax;

    int i = pphot->i1p[ip];
    int j = pphot->i2p[ip];
    int k = pphot->i3p[ip];

    // 2. Calculate probability of photon being emitted from each source
    // -----------------------------------------------------------------
    // SetEmissionCellWeight has chosen a zone for this photon, so now we look at the
    // emissivity arrays for that zone and sample a type for the photon.
    // Could technically do this calculation outside of the loop over photons, might
    // save a small amount of time
    Real vol_lya = emis(0,k,j,i) + emis(1,k,j,i);
    Real sur_lya = emis(2,k,j,i);
    Real sur_ion = emis(3,k,j,i);
    Real emis_tot;

    if (!flag_lya_volume_emis) vol_lya = 0.0;
    if (!flag_lya_surface_emis) sur_lya = 0.0;
    if (!flag_ion_surface_emis) sur_ion = 0.0;

    // 3. Sample emission probabilities to determine photon type
    // ---------------------------------------------------------
    // Sample random numbers, compare with probability of photon belonging to each
    // source of emissivity

    if (pmy_mc->emission_method == ION_STR) {
      // photoionization pass
      emis_tot = sur_ion;
      emis_geometry = SURFACE;
      phot_type = IONIZING;
    } else { // main monte carlo simulation
      emis_tot = sur_lya + vol_lya;
      Real frac_vol;
      if (emis_tot == 0.0) {
        printf("TOTAL EMISSIVITY IS ZERO FOR LYMAN ALPHA - SHOULD NOT HAPPEN");
      } else {
        frac_vol = vol_lya / emis_tot;
      }
      if (pran->uniform() <= frac_vol) {
        emis_geometry = VOLUME;
        phot_type = LYA;
      } else {
        emis_geometry = SURFACE;
        phot_type = LYA;
      }
    }

    // 4. Calculate emission geometry
    // ------------------------------
    // Set zone index, position within zone, initial direction vector for the case of
    // either volume emission or surface emission
    switch (emis_geometry) {

      case NO_GEOMETRY: {
        msg << "### FATAL ERROR in function MonteCarloBlock::InitializePhoton"
            << std::endl << "Photon has no emission geometry defined" <<std::endl;
        ATHENA_ERROR(msg);
        break;
      } // END CASE NO_GEOMETRY

      case VOLUME: {
        // Obtain initial position within zone
        GetZonePosition(pphot,pran,pcoord,ip);

        // Sample isotropic angle in orthonormal basis
        Real mu = 2.*pran->uniform()-1.0;
        Real stheta = std::sqrt(1.0-mu*mu);
        Real phi = 2.*PI*pran->uniform();
        pphot->k0p[ip] = 1. / c_cgs;
        pphot->k1p[ip] = stheta * std::cos(phi);
        pphot->k2p[ip] = stheta * std::sin(phi);
        pphot->k3p[ip] = mu;

        // Resize direction vector to sphpol code coords
        //if (pphot->general_pusher_flag) {
        //  pphot->k2p[ip] /= pphot->x1p[ip];
        //  pphot->k3p[ip] /= pphot->x1p[ip] * sin(pphot->x2p[ip]);
        //}

        break;
      } // END CASE VOLUME

      case SURFACE: {

        // Set initial position on the outer radial surface of this zone
        // Use a rejection method for theta and phi within the zone's bounds

        Real ph, th;
        Real phmax = pmy_block->pcoord->x3f(k+1);
        Real phmin = pmy_block->pcoord->x3f(k);
        Real thmax = pmy_block->pcoord->x2f(j+1);
        Real thmin = pmy_block->pcoord->x2f(j);
        Real sinphmax = std::sin(phmax);
        Real sinphmin = std::sin(phmin);
        Real sinthmax = std::sin(thmax);
        Real sinthmin = std::sin(thmin);

        if (flag_incident_from_z) {

          Real sthp = std::sin(thmax);
          Real sthm = std::sin(thmin);
          Real unif_1 = pran->uniform();
          Real sth_samp = std::sqrt(unif_1*(SQR(sthp) - SQR(sthm)) + SQR(sthm));
          Real th = std::asin(sth_samp);
          if (th < 0.0) {
            printf("WARNING: arcsin domain - sampled theta was %g\n", th);
            th += PI;
          }

          Real unif_2 = pran->uniform();
          Real ph_samp = unif_2 * (phmax - phmin) + phmin;

          // r position - R
          pphot->x1p[ip] = pmy_block->pmy_mesh->mesh_size.x1max;

          // theta position - [0, pi]
          pphot->x2p[ip] = th;

          // phi position - [0, 2pi]
          pphot->x3p[ip] = ph_samp;

          // Set direction vector - parallel to star-planet separation vector
          pphot->k0p[ip] = 1.;
          pphot->k1p[ip] = -std::cos(th);
          pphot->k2p[ip] = std::sin(th);
          pphot->k3p[ip] = 0.0;

          //Real sth = std::sin(th);
          //Real cth = std::cos(th);
          //Real sph = std::sin(ph_samp);
          //Real cph = std::cos(ph_samp);

          //Real xhat = sth*cph*pphot->k1p[ip] + cth*cph*pphot->k2p[ip] - sph*pphot->k3p[ip];
          //Real yhat = sth*sph*pphot->k1p[ip] + cth*sph*pphot->k2p[ip] + cph*pphot->k3p[ip];
          //Real zhat = cth*pphot->k1p[ip] - sth*pphot->k2p[ip];
          //printf("xhat = %g      yhat = %g       zhat = %g\n", xhat, yhat, zhat);

        } else { // NOT INCIDENT_FROM_Z

         // // Sample phi from probability distribution within meshblock bounds
         // Real unif_1 = pran->uniform();
         // Real cph_samp = unif_1 * (std::cos(phmax) - std::cos(phmin)) + std::cos(phmin);
         // Real ph_samp = std::atan2(std::sqrt(1.0 - SQR(cph_samp)), cph_samp); // TODO: Does this work?
         //                                                                      // This assumes sin phi is always positive

         // // Return domain to 0 < ph < 2 PI
         // if (ph_samp < 0) {
         //   ph_samp += 2.0 * PI;
         // }
         // ph = ph_samp;
         // Real sph = std::sin(ph);

         // // Use rejection method to sample theta
         // bool reject = true;
         // int nreject = 0;
         // while (reject) {

         //   // Sample uniformly within the bounding box
         //   Real th_sample = pran->uniform() * (thmax - thmin) + thmin;

         //   // Evaluate theta probability at the sampled theta
         //   Real sth_sample = sin(th_sample);
         //   Real th_prob = SQR(sth_sample) * sph;

         //   // Draw a random number between 0 and 1 (uniform)
         //   Real unif_2 = pran->uniform();

         //   nreject++;
         //   // Keep if uniform sample is under the curve of the theta probability
         //   // distribution
         //   if (unif_2 < th_prob) {
         //     reject = false;
         //     th = th_sample;
         //   }
         // }
         // //printf("Rejection efficiency: %g\n", 1.0/static_cast<Real>(nreject));

         // CMF: sampling functions for phi and theta within meshblock boundaries
         // Assumes star is in the -x direction (psi = 0)
          Real unif_ph = pran->uniform();
          ph = PI - std::asin(unif_ph * (sinphmax-sinphmin) + sinphmin);

          Real prob0 = 0.25 * (2*(thmax-thmin) + std::sin(2*thmin) - std::sin(2*thmax));
          Real probmax;
          if (thmax <= PI/2) {
            probmax = SQR(sinthmax) / prob0;
          } else if (thmin >= PI/2) {
            probmax = SQR(sinthmin) / prob0;
          } else {
            probmax = 1/prob0;
          }

          bool reject = true;
          int step = 0;
          while (reject) {
            Real unif_th = pran->uniform();
            Real th_samp = unif_th * (thmax-thmin) + thmin;
            Real sth_samp = std::sin(th_samp);
            Real th_prob = SQR(sth_samp)/prob0;
            Real unif_2 = pran->uniform();
            if (unif_2 <= th_prob) {
              reject = false;
              th = th_samp;
            }
            step++;
          }

          // r position - R
          pphot->x1p[ip] = pmy_block->pmy_mesh->mesh_size.x1max;

          // theta position - [0, pi]
          pphot->x2p[ip] = th;

          // phi position - [0, 2pi]
          pphot->x3p[ip] = ph;

          Real sth = std::sin(th);
          Real cth = std::cos(th);
          Real sph = std::sin(ph);
          Real cph = std::cos(ph);

          // Set direction vector - parallel to star-planet separation vector
          pphot->k0p[ip] = 1.;
          pphot->k1p[ip] = sth*cph;
          pphot->k2p[ip] = cth*cph;
          pphot->k3p[ip] = -sph;

         // Real xhat = sth*cph*pphot->k1p[ip] + cth*cph*pphot->k2p[ip] - sph*pphot->k3p[ip];
         // Real yhat = sth*sph*pphot->k1p[ip] + cth*sph*pphot->k2p[ip] + cph*pphot->k3p[ip];
         // Real zhat = cth*pphot->k1p[ip] - sth*pphot->k2p[ip];
         // printf("photon x, y, z directions: %g, %g, %g\n", xhat, yhat, zhat);

        } // END IF INCIDENT_FROM_Z

        // Resize direction vector to sphpol code coords
        //if (pphot->general_pusher_flag) {
        //  pphot->k2p[ip] /= pphot->x1p[ip];
        //  pphot->k3p[ip] /= pphot->x1p[ip] * sin(pphot->x2p[ip]);
        //}

        // Set photon position indices
        //Real xi1, xi2, xi3;
        //pmy_block->pcoord->MeshCoordsToIndices(pphot->x1p[ip], pphot->x2p[ip],
        //                                       pphot->x3p[ip], xi1, xi2, xi3);
        //pphot->i1p[ip] = static_cast<int>(xi1);
        //pphot->i2p[ip] = static_cast<int>(xi2);
        //pphot->i3p[ip] = static_cast<int>(xi3);

        break;
      } // END CASE SURFACE

    } // END SWITCH EMIS_GEOMETRY

    bool dayside = true;
    if (flag_incident_from_z) {
      Real pth = pphot->x2p[ip];
      if (pth >= PI/2) {
        dayside = false;
      }
    } else {
      Real pph = pphot->x3p[ip];
      if (pph <= PI/2 || pph >= 3*PI/2) {
        dayside = false;
      }
    }

    // 5. Set photon type
    // ------------------
    // Set the photon's energy and absorption / scattering opacities
    switch (phot_type) {
      case NO_TYPE: {
        msg << "### FATAL ERROR in function MonteCarloBlock::InitializePhoton"
            << std::endl << "Photon has no type defined" <<std::endl;
        ATHENA_ERROR(msg);
        break;
      }
      case LYA: {
        if (emis_geometry == SURFACE) {
          Real energy;
          gasdev(pmy_block, energy0, stddev, energy);
          while (std::fabs(energy - energy0) > linewidth_cutoff_energy) {
            gasdev(pmy_block, energy0, stddev, energy);
          }
					pphot->ep[ip] = energy;
          //printf("ip: %d\n initial energy: %g\n", ip, energy);
        } else {
          // Photon is emitted via recombination - should have line center energy
          pphot->ep[ip] = energy0;
        }
        break;
      }
      case IONIZING: {
        Real nu_phot;
        if (!flag_pow_law) {
          nu_phot = numin + kb_cgs * chromo_temp / h_cgs * std::log(1./(1-pran->uniform()));
        } else {
          nu_phot = pow(pran->uniform() * (numaxpow - numinpow) + numinpow, nuexp);
        }
        pphot->ep[ip] = h_cgs * nu_phot;
        break;
      }
    }

    // Set status flag
    if (pphot->wp[ip] < 0.0)
      pphot->statp[ip] = DESTROYED;
    else
      pphot->statp[ip] = EVOLVING;

    // CMF: quick and dirty killer of nightside photons
    if (pphot->wp[ip] == 0.0) {
      printf("Warning: photon with zero weight found.\n");
      pphot->PrintPhoton(ip);
      pphot->statp[ip] = ABSORBED;
    }

    // Initialize the absorption and scattering extinction coefficients
    // to the values in the emitted zone
    pphot->acp[ip] = AbsorptionOpacity(this,pphot,ip);
    pphot->scp[ip] = ScatteringOpacity(this,pphot,ip);

  } // end loop over ip
}


// Function called each time a Lya photon escapes
/*void MonteCarloBlock::FinalizePhoton(Photon *pphot, int ip) {
  // For the pure absorption test, set user variable to be the cross section
  int i = pphot->i1p[ip];
  int j = pphot->i2p[ip];
  int k = pphot->i3p[ip];
  Real nH = pmy_block->pscalars->s(0,k,j,i)/mp_cgs;
  Real sigma = pphot->scp[ip] / nH;
  pphot->user[0][ip] = sigma;
  return;
}*/


/*
 * Bottom (r_inner) boundary condition: static inflow diode
 * If radial velocity is positive (outwards), copy it
 * If radial velocity is negative (inwards), set it to zero
 * Copy any theta and phi velocities
 * Density and pressure are fixed to initial values
*/
void StaticInflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
  FaceField &b, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku,
  int ngh) {

  Real gamma = pmb->peos->GetGamma();
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {

        // fix initial density and pressure
        prim(IDN,k,j,il-i) = rho_gz[il-i];
        //prim(IPR,k,j,il-i) = pres_gz[il-i];
        prim(IPR,k,j,il-i) = prim(IPR,k,j,il) * std::pow(rho_gz[il-i]/prim(IDN,k,j,il), gamma);

        // outflow diode for velocity
        prim(IVY,k,j,il-i) = prim(IVY,k,j,il);
        prim(IVZ,k,j,il-i) = prim(IVZ,k,j,il);
        if (prim(IVX,k,j,il) >= 0.0) {
          prim(IVX,k,j,il-i) = prim(IVX,k,j,il);
        } else {
          prim(IVX,k,j,il-i) = 0.0;
        }

        // set neutral fraction to fully neutral
        pmb->pscalars->r(0,k,j,il-i) = 1.0;
        pmb->pscalars->s(0,k,j,il-i) = prim(IDN,k,j,il-i);
      }
    }
  }
  return;
}


/*
 * Bottom (r_inner) boundary condition: outflow diode
 * If radial velocity is positive (outwards), copy it
 * If radial velocity is negative (inwards), set it to zero
 * Density, pressure, ionization are fixed to initial values
 * Copy all other vars
*/
void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
  FaceField &b, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku,
  int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        if (prim(IVX,k,j,il) >= 0.0) {
          prim(IVX,k,j,il-i) = prim(IVX,k,j,il);
        } else {
          prim(IVX,k,j,il-i) = 0.0;
        }
        prim(IVY,k,j,il-i) = prim(IVY,k,j,il);
        prim(IVZ,k,j,il-i) = prim(IVZ,k,j,il);

        prim(IPR,k,j,il-i) = P_gz[ngh-i];
        prim(IDN,k,j,il-i) = rho_gz[ngh-i];

        pmb->pscalars->r(0,k,j,il-i) = 1.0;
        pmb->pscalars->s(0,k,j,il-i) = prim(IDN,k,j,il-i);
      }
    }
  }
  return;
}


/*
 * Top (r_outer) boundary condition: outflow diode
 * If velocity is positive (outwards), copy velocity, density and pressure
 * Otherwise, if velocity is negative (inwards), set all to zero
 * Copy all other vars
*/
void OutflowOuterX1 (MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
  FaceField &b, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku,
  int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        Real dr = pco->x1v(iu+i) - pco->x1v(iu);
        if (prim(IVX,k,j,iu) >= 0.0) {
          prim(IVX,k,j,iu+i) = prim(IVX,k,j,iu);
          prim(IPR,k,j,iu+i) = prim(IPR,k,j,iu);
          prim(IDN,k,j,iu+i) = prim(IDN,k,j,iu);
        } else {
          prim(IVX,k,j,iu+i) = 0.0;
          prim(IPR,k,j,iu+i) = pfloor;
          prim(IDN,k,j,iu+i) = dfloor;
        }
        prim(IVY,k,j,iu+i) = prim(IVY,k,j,iu);
        prim(IVZ,k,j,iu+i) = prim(IVZ,k,j,iu);
      }
    }
  }
  return;
}

void MonteCarloBlock::UserWorkAfterTransfer(int etype) {

  // only update ionization after transfer of ionizing photons
  if (etype != ION_STR)
    return;

  Real dt = pmy_block->pmy_mesh->dt;
  // for checking ionization
  //dt = 1.e10;
  Real tint = pmy_mc->tint;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        Real norm = 1./ (tint * pcoord->vol(k,j,i));
        // Get temperatures to calculate recombination and impact excitation terms
        Real tempo1e4K = tgas(k,j,i) / 1.e4;
        Real invtemp = 1/tempo1e4K;

        // Do an implicit update of the neutral fraction to solve for the ionization state at the end of this step
        Real rho = pmy_block->phydro->u(IDN,k,j,i); // SWD: Why not use MCBlock rho?
        //Real neutral_frac = pmy_block->pscalars->s(0,k,j,i)/rho;// Neutral fraction
        //Real mean_mol_weight = (1. + neutral_frac)/2.;
        //printf("mean mol weight: %g\n", mean_mol_weight);
        const Real mass = 1.660538782e-24;
        Real nh = pmy_block->pscalars->s(0,k,j,i)/mass;
        Real np = (rho/mass) - nh;
        Real na = nh + np;

        // Calculate the photoionization rate, Gamma, from the number of photons absorbed
        // per cell per time
        // absweight already normalized by (vol*dt) in NormalizeSourceTerms
        Real absweight = sourceterms(MCNABS,k,j,i) *= norm;
        Real vol = pcoord->vol(k,j,i);
        Real Gamma = absweight / nh;


        // Calculate the recombination rate, alpha, from the temperature of this cell
        Real alpha = 2.54e-13 * std::pow(tempo1e4K, -0.8164-0.0208*std::log(tempo1e4K));

        Real nC, nR;
        nR = 1. / alpha / dt;
        nC = Gamma / alpha;

        // Calculate the update to the neutral H number density
        //Real discriminant = SQR(nR) + 4*na*(nR+nC) + 2*nR*nC + SQR(nC) - 4*nh*nR;
        //Real update = na + 0.5*(nR + nC - std::sqrt(discriminant));
        //Real update = 0.5*(nR+2*na+nC) - 0.5*std::sqrt(discriminant);

        //CMF: see numerical recipes 5.6
        Real bb = 2*na + nC + nR;
        Real cc = nh*nR + SQR(na);
        Real dd = SQR(bb) - 4*cc;
        Real qq = 0.5*(bb + std::sqrt(dd));
        Real update = cc/qq;

        // Double-check that the implicit update gives a neutral fraction between 0 and 1 --- NOT guaranteed if dt is large
        Real neutral_frac = update/na;
        if (neutral_frac > 1.0) {
          if ((k==ks)&&(j==js)&&(i==14))
            printf("(Block %d) UpdateIonizationFraction: neutral fraction %g is greater than 1.0\n", pmy_block->lid, neutral_frac);
          neutral_frac = 1.0;
        } else if (neutral_frac < 0.0) {
          if ((k==ks)&&(j==js)&&(i==14))
            printf("(Block %d) UpdateIonizationFraction: neutral fraction %g is less than 0.0\n", pmy_block->lid, neutral_frac);
          neutral_frac = 0.0;
        }

        // check result
       // if (absweight > 0.0) {
       //   printf("nh'=%g, na=%g, nh0=%g, alpha=%g, Gamma=%g, dt=%g\n", update, na, nh, alpha, Gamma, dt);
       // }

        nh = neutral_frac * na;
        np = na - nh;

        pmy_block->pscalars->s(0,k,j,i) = nh*mass;

        // calculate impact excitation cooling
        // see Christie, Arras, Li, 2013, eq.8 and Table 2
        Real c1s2s = 1.21e-8 * std::pow(invtemp, 0.455) * std::exp(-11.84/tempo1e4K);
        Real c1s2p = 1.71e-8 * std::pow(invtemp, 0.077) * std::exp(-11.84/tempo1e4K);
        Real ctot = c1s2s + c1s2p;

        const Real ev_to_erg = 1.602176634e-12;
        Real cool = ctot*nh*np*10.2*ev_to_erg;
        //printf("cool, ctot, nh, np: %g, %g, %g, %g\n", cool, ctot, nh, np);
        //sourceterms(MCRS0,k,j,i) -= cool;
        // CMF: modify sourceterms in unnormalized way, so that final result can be normalized in RunDynamic
        sourceterms(MCRS0,k,j,i) -= cool*vol*dt;
        pmy_block->user_out_var(5,k,j,i) += cool;
        //if ((sourceterms(MCRS0,k,j,i) > 0.0) || (sourceterms(MCRS0,k,j,i) < 0.0))
          //printf("from UpdateIonizationFraction: %g, %g, %d, %d, %d\n", sourceterms(MCRS0,k,j,i), cool, k, j, i);
      }
    }
  }

}

// Definitions for namespace functions
namespace { // begin namespace
            //
void F1pos(MonteCarloBlock *pmcb, Photon *pphot, Real dl, int ip, int imom) {

  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  const Real c_cgs = 2.99792458e10;
  if (pphot->k1p[ip] >= 0.0) {
    Real weight = pphot->ep[ip]*pphot->wp[ip]*dl*pphot->k1p[ip]/c_cgs;
    pmcb->moments_user(imom,i3,i2,i1) += weight;
  }
}

void F1neg(MonteCarloBlock *pmcb, Photon *pphot, Real dl, int ip, int imom) {

  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  const Real c_cgs = 2.99792458e10;
  if (pphot->k1p[ip] < 0.0) {
    Real weight = pphot->ep[ip]*pphot->wp[ip]*dl*pphot->k1p[ip]/c_cgs;
    pmcb->moments_user(imom,i3,i2,i1) += weight;
  }
}

void a1pos(MonteCarloBlock *pmcb, Photon *pphot, Real dl, int ip, int imom) {

  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  const Real c_cgs = 2.99792458e10;
  if (pphot->k1p[ip] >= 0.0) {
    Real weight = pphot->ep[ip]*pphot->wp[ip]*dl*pphot->k1p[ip]/c_cgs * (pphot->scp[ip] + pphot->acp[ip]);
    pmcb->moments_user(imom,i3,i2,i1) += weight;
  }
}

void a1neg(MonteCarloBlock *pmcb, Photon *pphot, Real dl, int ip, int imom) {

  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  const Real c_cgs = 2.99792458e10;
  if (pphot->k1p[ip] <= 0.0) {
    Real weight = pphot->ep[ip]*pphot->wp[ip]*dl*pphot->k1p[ip]/c_cgs * (pphot->scp[ip] + pphot->acp[ip]);
    pmcb->moments_user(imom,i3,i2,i1) += weight;
  }
}

void GetTidalAcceleration(Real r, Real th, Real ph, Real rho, Real &a_r, Real &a_th, Real &a_ph) {

  Real sth = std::sin(th);
  Real cth = std::cos(th);//std::sqrt(1.0 - SQR(sth));

  Real phm = ph - psi;
  Real sphm = std::sin(phm);
  Real cphm = std::cos(phm);//std::sqrt(1.0 - SQR(sphm));

  Real prefac = gm_star*std::pow(SQR(r)+2.0*sep*r*sth*cphm+SQR(sep),-1.5);
  a_r = - prefac*(r+sep*sth*cphm)
        - gm_planet/(SQR(r))
        + r*SQR(sth)*(gm_star+gm_planet)/(sep*sep*sep)
        + sth*cphm*gm_star/(SQR(sep));
  a_th= - prefac*(sep*cth*cphm)
        + r*sth*cth*(gm_star+gm_planet)/(sep*sep*sep)
        + cth*cphm*gm_star/(SQR(sep));
  a_ph= - prefac*(sep*sphm)
        - sphm*gm_star/(SQR(sep));
  return;
}

void TwoPointMass(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar) {

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real r = pmb->pcoord->x1v(i);
        Real th = pmb->pcoord->x2v(j);
        Real ph = pmb->pcoord->x3v(k);
        Real rho = prim(IDN,k,j,i);

				// Tidal gravity + Centrifugal
        Real a_r, a_th, a_ph;
        GetTidalAcceleration(r, th, ph, rho, a_r, a_th, a_ph);

        // Coriolis
        // Omega = (GM/a^3)(costh e_r - sinth e_th), a_cor = -2Omega cross vel
        Real omega = std::sqrt((gm_star+gm_planet)/(sep*sep*sep));
        Real sinth = std::sin(th);
        Real costh = std::cos(th);
        Real vel_r = prim(IVX,k,j,i);
        Real vel_th = prim(IVY,k,j,i);
        Real vel_ph = prim(IVZ,k,j,i);
        a_r += 2*omega * vel_ph*sinth;
        a_th += 2*omega * vel_ph*costh;
        a_ph -= 2*omega * (vel_r*sinth + vel_th*costh);

        Real src1 = a_r*rho*dt;
        Real src2 = a_th*rho*dt;
        Real src3 = a_ph*rho*dt;
        cons(IM1,k,j,i) += src1;
        cons(IM2,k,j,i) += src2;
        cons(IM3,k,j,i) += src3;

        // Update the user meshblock quantities
        pmb->ruser_meshblock_data[1](0,k,j,i) = a_r;
        pmb->ruser_meshblock_data[1](1,k,j,i) = a_th;
        pmb->ruser_meshblock_data[1](2,k,j,i) = a_ph;

        // Update conserved gas energy
        cons(IEN,k,j,i) += src1*prim(IVX,k,j,i)
                         + src2*prim(IVY,k,j,i)
                         + src3*prim(IVZ,k,j,i);
      }
    }
  }
}

// Assumes the star is in the -x direction, i.e. psi = 0
void HillTidalGravity(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
    AthenaArray<Real> &cons_scalar) {
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real rho = prim(IDN,k,j,i);
        Real r = pmb->pcoord->x1v(i);
        Real th = pmb->pcoord->x2v(j);
        Real ph = pmb->pcoord->x3v(k);
        Real vel_r = prim(IVX,k,j,i);
        Real vel_th = prim(IVY,k,j,i);
        Real vel_ph = prim(IVZ,k,j,i);
        Real cosph = std::cos(ph);
        Real sinph = std::sin(ph);
        Real costh = std::cos(th);
        Real sinth = std::sin(th);
        Real x = r * cosph * sinth;
        //Real y = r * sinph * sinth;
        Real z = r * costh;

        // Tidal gravity, second order
        Real prefac = gm_star / (sep*sep*sep);
        Real a_r = prefac * (3*x*sinth*cosph - z*costh);
        Real a_th = prefac * (3*x*costh*cosph + z*sinth);
        Real a_ph = prefac * (-3*x*sinph);

        // Coriolis, second order
        Real omega = std::sqrt((gm_star+gm_planet)/(sep*sep*sep));
        a_r += 2*omega * vel_ph*sinth;
        a_th += 2*omega * vel_ph*costh;
        a_ph -= 2*omega * (vel_r*sinth + vel_th*costh);

        Real src1 = a_r*rho*dt;
        Real src2 = a_th*rho*dt;
        Real src3 = a_ph*rho*dt;
        cons(IM1,k,j,i) += src1;
        cons(IM2,k,j,i) += src2;
        cons(IM3,k,j,i) += src3;

        // Update the user meshblock quantities
        pmb->ruser_meshblock_data[1](0,k,j,i) = a_r;
        pmb->ruser_meshblock_data[1](1,k,j,i) = a_th;
        pmb->ruser_meshblock_data[1](2,k,j,i) = a_ph;

        // Update conserved gas energy
        cons(IEN,k,j,i) += src1*prim(IVX,k,j,i)
                         + src2*prim(IVY,k,j,i)
                         + src3*prim(IVZ,k,j,i);
      }
    }
  }
}

// Assumes the star is in the -x direction, i.e. psi = 0
void ThirdOrderTidalGravity(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
  AthenaArray<Real> &cons_scalar) {
	Real prefac2 = gm_star / (sep*sep*sep);
	Real prefac3 = 3*gm_star / (2*sep*sep*sep*sep);
	Real omega = std::sqrt((gm_star+gm_planet)/(sep*sep*sep));

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
		Real ph = pmb->pcoord->x3v(k);
		Real cosph = std::cos(ph);
		Real sinph = std::sin(ph);

    for (int j=pmb->js; j<=pmb->je; ++j) {
			Real th = pmb->pcoord->x2v(j);
			Real costh = std::cos(th);
			Real sinth = std::sin(th);

			// unit vector conversions
			Real unit_xr = sinth*cosph;
			Real unit_xt = costh*cosph;
			Real unit_xp = -sinph;

			Real unit_yr = sinth*sinph;
			Real unit_yt = costh*sinph;
			Real unit_yp = cosph;

			Real unit_zr = costh;
			Real unit_zt = -sinth;

      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real r = pmb->pcoord->x1v(i);
				Real rho = prim(IDN,k,j,i);
				Real vel_r  = prim(IVX,k,j,i);
				Real vel_th = prim(IVY,k,j,i);
				Real vel_ph = prim(IVZ,k,j,i);
				Real x = r * sinth * cosph;
				Real y = r * sinth * sinph;
				Real z = r * costh;

				// Tidal gravity, second order
        Real a_r  = prefac2 * (3*x*unit_xr - z*unit_zr);
				Real a_th = prefac2 * (3*x*unit_xt + z*unit_zt);
				Real a_ph = prefac2 * (3*x*unit_xp);

				// Tidal gravity, third order
				a_r  += prefac3 * ((3*r*r - x*x)*unit_xr + 4*x*y*unit_yr + 4*x*z*unit_zr);
				a_th += prefac3 * ((3*r*r - x*x)*unit_xt + 4*x*y*unit_yt + 4*x*z*unit_zt);
				a_ph += prefac3 * ((3*r*r - x*x)*unit_xp + 4*x*y*unit_yp);

				// Coriolis
				a_r  += 2*omega * vel_ph*sinth;
				a_th += 2*omega * vel_ph*costh;
				a_ph -= 2*omega * (vel_r*sinth + vel_th*costh);

        Real src1 = a_r*rho*dt;
        Real src2 = a_th*rho*dt;
        Real src3 = a_ph*rho*dt;
        cons(IM1,k,j,i) += src1;
        cons(IM2,k,j,i) += src2;
        cons(IM3,k,j,i) += src3;

        // Update the user meshblock quantities
        pmb->ruser_meshblock_data[1](0,k,j,i) = a_r;
        pmb->ruser_meshblock_data[1](1,k,j,i) = a_th;
        pmb->ruser_meshblock_data[1](2,k,j,i) = a_ph;

        // Update conserved gas energy
        cons(IEN,k,j,i) += src1*prim(IVX,k,j,i)
                         + src2*prim(IVY,k,j,i)
                         + src3*prim(IVZ,k,j,i);
      }
    }
  }
}


// Real MultipleEmissivities(MonteCarloBlock *pmcb, int k, int j, int i) {
// // Wrapper function for combining each source of emissivity in the problem

//   Real emis_vol_recomb = VolumeEmissivityLya(pmcb, k, j, i);
//   Real emis_surf_lya = SurfaceEmissivityLya(pmcb, k, j, i);
//   Real emis_surf_ionizing = SurfaceEmissivityIonizing(pmcb, k, j, i);

//   if (!flag_lya_volume_emis) emis_vol_recomb = 0.0;
//   if (!flag_lya_surface_emis) emis_surf_lya = 0.0;
//   if (!flag_ion_surface_emis) emis_surf_ionizing = 0.0;

//   // If this is the photoionization pass, return only that emissivity
//   if (pmcb->pmy_mc->photoionization) {
//     return emis_surf_ionizing;
//   }
//   // Otherwise, return the combined emissivity of the non-ionizing sources
//   return emis_vol_recomb + emis_surf_lya;
// }


Real VolumeEmissivityLya(MonteCarloBlock *pmcb, int k, int j, int i, int etype) {
// Sets the value of the emission array which determines where Lya photons are
// likely to be emitted. The associated cooling of the gas is handled in
// UpdateSourceTerms, not here - this is just where the emissivity is calculated

   Real rho = pmcb->pmy_block->phydro->u(IDN,k,j,i);
   Real tempo1e4K = pmcb->tgas(k,j,i) / 1.0e4;
   Real invtemp = 1.0/tempo1e4K;
   //const Real mp = 1.660538782e-24;

   // LYA COOLING TERMS
   // c.f. Christie Arras & Li (2013), Table 2
   // ----------------------------------------

   // Recombination rate as a function of temperature
   // Ionized hydrogen and free electrons combine to produce Lya
   // Energy per volume per time = 10.2 eV * alpha * np * ne
   Real alpha = 2.54e-13 * std::pow(tempo1e4K, -0.8164-0.0208*std::log(tempo1e4K));

   // Electron impact excitation rate as a function of temperature
   // We presume thermal electrons excite neutral H(1s) to the n=2 state (either 2s or 2p)
   // which then immediately de-excite, producing Lya
   // Energy per volume per time = 10.2 eV * ctot * nH * ne
   Real c1s2s = 1.21e-8 * std::pow(invtemp, 0.455) * std::exp(-11.84/tempo1e4K);
   Real c1s2p = 1.71e-8 * std::pow(invtemp, 0.077) * std::exp(-11.84/tempo1e4K);
   Real ctot = c1s2s + c1s2p;

   // Number density of neutral H
	 //Real neutral_frac = pmcb->pmy_block->pscalars->s(0,k,j,i) / rho;
	 //Real mean_mol_weight = (1. + neutral_frac) / 2.;
   //Real nH = neutral_frac * rho / mean_mol_weight / mp_cgs;
	 Real nH = pmcb->pmy_block->pscalars->s(0,k,j,i)/mp_cgs;

   // Number density of protons (= number density of electrons)
   //Real np = (1.-mean_mol_weight) * rho / mean_mol_weight /	mp_cgs;
	 Real np = (rho/mp_cgs) - nH;

   Real recombination = alpha*SQR(np);
   Real impact = ctot*nH*np;
   Real emis = recombination + impact;

  //Real emis = 0.0;
  //Real mmw0 = 1.0;                      // mean molecular weight at base
  //Real a2 = kb_cgs * temp0 / mmw0 / mp_cgs;     // sound speed squared at base
  //Real lambda = gm_planet / (rin * a2);

  //Real r = pmcb->pmy_block->pcoord->x1v(i);  // r coordinate
  //Real H = rin / lambda * (SQR(r) / SQR(rin)); // scale height
  //Real nH = pmcb->scalars(k,j,i) / mp_cgs;

  //Real Gamma0 = 4.e-5; // 1/s, photoionization rate coefficient
  //Real sigmapi = 6.e-18;

  //Real ion_rate_atten = Gamma0 / (1.0 + std::pow(nH*sigmapi*H, 1.5));
  //emis = ion_rate_atten * nH;

  pmcb->pmy_block->user_out_var(0,k,j,i) = emis;

  pmcb->pmy_block->ruser_meshblock_data[0](0,k,j,i) = recombination;
  pmcb->pmy_block->ruser_meshblock_data[0](1,k,j,i) = impact;
  return emis;
}


Real SurfaceEmissivityLya(MonteCarloBlock *pmcb, int k, int j, int i, int etype) {

  Coordinates *pco = pmcb->pmy_block->pcoord;

  Real emis = 0.0;
  if (pco->x1f(i+1) >= rout) {
    if (flag_incident_from_z) {
      if (pco->x2v(j) >= PI/2.)
        return 0.0;
      Real cthm = std::cos(pco->x2f(j));
      Real cthp = std::cos(pco->x2f(j+1));
      emis = 0.5*edot_lya/PI*(pco->x3f(k+1)-pco->x3f(k))*(SQR(cthm)-SQR(cthp))/energy0;
    } else {
      Real thm = pco->x2f(j);
      Real thp = pco->x2f(j+1);
      Real phm = pco->x3f(k);
      Real php = pco->x3f(k+1);

     // if (phm < 0.) {
     //   phm += 2.0 * PI;
     // } else if (phm >= 2.0 * PI) {
     //   phm -= 2.0 * PI;
     // }

     // if (php < 0) {
     //   php += 2.0 * PI;
     // } else if (php >= 2.0 * PI) {
     //   php -= 2.0 * PI;
     // }

     // bool dayside = true;
     // if ((phm < PI/2.) || (php > 3.*PI/2.)) {
     //   dayside = false;
     // }

      //printf("k=%d  j=%d  i=%d     phm=%g     php=%g    dayside=%d\n", k,j,i,phm,php,dayside);
      // Zero emissivity if phi coordinate does not fall within dayside bounds
      if ((phm < PI/2.) || (php > 3.*PI/2.)) {
        return 0.0;
      }

      Real sthm = std::sin(thm);
      Real cthm = std::cos(thm);
      Real sthp = std::sin(thp);
      Real cthp = std::cos(thp);
      Real sphm = std::sin(phm);
      Real sphp = std::sin(php);

      emis = -0.5*edot_lya/PI*(thp - thm - (std::sin(thp - thm)*std::cos(thp + thm)))*(sphp-sphm)/energy0;
    }
  }

  Real vol = pco->GetCellVolume(k,j,i);
//	if (emis > 0.)
//	printf("lya emis, i, j, k: %g, %d, %d, %d\n", emis/vol, i, j, k);

  pmcb->pmy_block->user_out_var(1,k,j,i) = emis/vol;
  pmcb->pmy_block->ruser_meshblock_data[0](2,k,j,i) = emis/vol;
  return emis/vol;
}


Real SurfaceEmissivityIonizing(MonteCarloBlock *pmcb, int k, int j, int i, int etype) {
  Coordinates *pco = pmcb->pmy_block->pcoord;


  Real emis = 0.0;
  if (pco->x1f(i+1) >= rout) {
    if (flag_incident_from_z) {
      if (pco->x2v(j) >= PI/2.)
        return 0.0;
      Real cthm = std::cos(pco->x2f(j));
      Real cthp = std::cos(pco->x2f(j+1));
      Real r2 = pco->x1f(i+1)*pco->x1f(i+1);
      Real nflux = r2 * ion_flux/(h_cgs*mean_nu)/pco->GetFace1Area(k,j,i+1);
      emis = 0.5*nflux*(pco->x3f(k+1)-pco->x3f(k))*(SQR(cthm)-SQR(cthp));

    } else {
      Real phm = pco->x3f(k);
      Real php = pco->x3f(k+1);
      if ((phm < PI/2.) || (php > 3.*PI/2.)) {
        return 0.0;
      }
      Real thm = pco->x2f(j);
      Real thp = pco->x2f(j+1);

      Real sthm = std::sin(thm);
      Real cthm = std::cos(thm);
      Real sthp = std::sin(thp);
      Real cthp = std::cos(thp);
      Real sphm = std::sin(phm);
      Real sphp = std::sin(php);
      Real r2 = pco->x1f(i+1)*pco->x1f(i+1);
      Real nflux = r2*ion_flux/(h_cgs*mean_nu) / pco->GetFace1Area(k,j,i+1);
      emis = -0.5*nflux*(thp-thm-(std::sin(thp - thm)*std::cos(thp + thm)))*(sphp-sphm);
    }
  }


	// CMF note: currently no user_out_var space for EUV only?
  pmcb->pmy_block->ruser_meshblock_data[0](3,k,j,i) = emis;

  return emis;
}


void ResonantScattering(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe) {
  if (!flag_escape_after_scatter) {
    ScatterResonanceLine(pmcb, pphot, ips, ipe);
  } else {
    ForceEscape(pmcb, pphot, ips, ipe);
  }
}

Real BoundFreeAbsorptionOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {
  Real opac = 0.0;
	if (!flag_zero_opacity) {
		if (pphot->ep[ip] >= numin_erg) {
			int i1 = pphot->i1p[ip];
			int i2 = pphot->i2p[ip];
			int i3 = pphot->i3p[ip];
			Real energy = pphot->ep[ip];
			Real xsec = sigmamin * std::pow((energy / h_cgs) / (numin), -3.0);
			//Real rho = pmcb->pmy_block->phydro->u(IDN,i3,i2,i1);
			//Real neutral_frac = pmcb->pmy_block->pscalars->s(0,i3,i2,i1) / rho;
			//Real mean_mol_weight = (1. + neutral_frac) / 2.;
			//Real nH = neutral_frac * rho / mean_mol_weight / mp_cgs;
			Real nH = pmcb->pmy_block->pscalars->s(0,i3,i2,i1)/mp_cgs;
			opac = xsec * nH;
			//printf("In BoundFreeAbsorptionOpacity; energy, nH, mfp: %g, %g, %g\n", energy/1.6e-12, nH, 1./opac);
		}
	}
  return opac;
}

Real ResonantScatteringOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {
  Real opac = 0.0;
  Real erg_away_from_lc = std::fabs(pphot->ep[ip] - energy0);
  if (erg_away_from_lc <= linewidth_cutoff_energy) {
    if (!flag_zero_opacity) {
      opac = ResonanceLineOpacity(pmcb, pphot, ip);
    }
  }
  return opac;
}


void gasdev(MeshBlock *pmb, Real mean, Real stddev, Real &samp) {
  Real fac, rsq, v1, v2;

  if (iset == 0) {
    do {
      v1 = 2.0f * pmb->pmy_mcb->pran->uniform() - 1.0f;
      v2 = 2.0f * pmb->pmy_mcb->pran->uniform() - 1.0f;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0f || rsq == 0.0f);

    fac = std::sqrt(-2.0f * std::log(rsq) / rsq);
    gset = v1 * fac * stddev + mean;
    samp = v2 * fac * stddev + mean;
    iset = 1;
  } else {
    samp = gset;
    iset = 0;
  }
  return;
}

void ForceEscape(MonteCarloBlock *pmcb, Photon *pphot, int ips, int ipe) {
  for (int ip=ips; ip<=ipe; ip++) {
    pphot->statp[ip] = ESCAPED;
  }
}

Real GetIsowindVelocity(Real x) {
  Real y = (x <= 1.0) ? std::exp(1.5-2.0/x) : 1.1;
  int l = 0;
  while (l<1000) {
    l++;
    Real f = 0.5*y*y - std::log(y*x*x) - 2.0/x + 1.5;
    Real dfdy = y - 1./y;
    if (std::abs(f) < 1e-10)
      break;
    y -= f / dfdy;
    l++;
  }
  return y;
}


void GetIonizationTemperature(MonteCarloBlock *pmcb) {
	Hydro* phydro = pmcb->pmy_block->phydro;

  int il = pmcb->is; int iu = pmcb->ie;
  int jl = pmcb->js; int ju = pmcb->je;
  int kl = pmcb->ks; int ku = pmcb->ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
				Real rho = phydro->u(IDN,k,j,i);
				//Real neutral_frac = pmcb->pmy_block->pscalars->s(0,k,j,i)/rho;
				Real mean_mol_weight = 1./(2. - pmcb->pmy_block->pscalars->r(0,k,j,i));
        Real tgas = phydro->w(IPR,k,j,i) * mean_mol_weight * mp_cgs / rho / kb_cgs;

        // apply temperature floor
        pmcb->tgas(k,j,i) = (tgas > pmcb->tfloor_cgs) ? tgas : pmcb->tfloor_cgs;
      }
    }
  }
}


void ExplicitEUVHeating(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
	Real rmax = pmb->pmy_mesh->mesh_size.x1max;
	Real total_vol = (4./3.)*PI*std::pow(rmax,3);
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
		//	if (pmb->pcoord->x2v(j) > PI/2) {
		//		continue;
		//	}
      for (int i = pmb->is; i <= pmb->ie; ++i) {
				Real rad = pmb->pcoord->x1v(i);
				Real the = pmb->pcoord->x2v(j);
				//Real cell_vol = pmb->pcoord->GetCellVolume(k,j,i);
				Real energy = dt*edot_ion/total_vol*(rmax/rad)*std::cos(the);
				//printf("rad, the, energy = %g, %g, %g\n", rad, the, energy);
        cons(IEN,k,j,i) += energy;
      }
    }
  }
  return;
}

void UpdateSourceTerms(MonteCarloBlock *pmcb, Photon *pphot, Real energy0, Real weight0,
                  Real k1p0, Real k2p0, Real k3p0, int ip) {

  // if continuous absorptioin, handle source terms in UpdateMoments()
  if (pmcb->absorption_meth == ABSTAU)
    return;

  Real hplanck = 6.62607015e-27;
  Real threshold = 3.28808816e+15 * hplanck;

  // Update soucterms for ionizing radiation
  if (energy0 > threshold) {
    Real heat = weight0 * (energy0 - threshold);
    //printf("energy0, weight0: %g %g %g\n", energy0, weight0);
    if ((std::isinf(heat)) || (std::isnan(heat))) {
      std::cout << "Warning: UpdateSourceTerms heating is : " << heat << std::endl;
      pphot->PrintPhoton(ip);
      pphot->statp[ip] = DESTROYED;
    } else if (heat >= TINY_NUMBER) {
      int &i = pphot->i1p[ip];
      int &j = pphot->i2p[ip];
      int &k = pphot->i3p[ip];
      pmcb->sourceterms(MCRS0,k,j,i) += heat;
      pmcb->sourceterms(MCNABS,k,j,i) += weight0;
    }
  }

  // update momentum source terms as usual
  Real k1 = pphot->k1p[ip];
  Real k2 = pphot->k2p[ip];
  Real k3 = pphot->k3p[ip];

  // Normalize k vector if using general mover in spherical polar coords
  if ((COORDINATE_SYSTEM == "spherical_polar") && (pphot->general_pusher_flag)) {
    k2 *= pphot->x1p[ip];
    k3 *= pphot->x1p[ip] * sin(pphot->x2p[ip]);
    k2p0 *= pphot->x1p[ip];
    k3p0 *= pphot->x1p[ip] * sin(pphot->x2p[ip]);
  }
  Real norm0 = sqrt(SQR(k1p0) + SQR(k2p0) + SQR(k3p0));
  if ((fabs(norm0-1.) > 1.0e-8) && (norm0 > 1.0e-8)) {
    k1p0 /= norm0;
    k2p0 /= norm0;
    k3p0 /= norm0;
  }
  Real norm = sqrt(SQR(k1) + SQR(k2) + SQR(k3));
  if ((fabs(norm-1.) > 1.0e-8) && (norm > 1.0e-8)) {
    k1 /= norm;
    k2 /= norm;
    k3 /= norm;
  }
  //norm0 = sqrt(SQR(k1p0) + SQR(k2p0) + SQR(k3p0));
  //printf("norm0: %f\n", norm0);
  //norm = sqrt(SQR(k1) + SQR(k2) + SQR(k3));
  //printf("norm: %f\n", norm);

  Real c_cgs = 2.99792458e10;
  // Components of momentum change --- assumes orthonormal basis
  Real dp1p = pphot->wp[ip] * k1 * pphot->ep[ip] / c_cgs
              - weight0 * k1p0 * energy0 / c_cgs;
  Real dp2p = pphot->wp[ip] * k2 * pphot->ep[ip] / c_cgs
              - weight0 * k2p0 * energy0 / c_cgs;
  Real dp3p = pphot->wp[ip] * k3 * pphot->ep[ip] / c_cgs
              - weight0 * k3p0 * energy0 / c_cgs;

  Real cool = (pphot->wp[ip] * pphot->ep[ip]) - (weight0 * energy0);
  //if (energy0 == 0.0)
  //  printf("weight, cool: %g %g\n",pphot->weight,cool);

  if ((std::isinf(cool)) || (std::isnan(cool))) {
    std::cout << "Warning: UpdateSourceTerms cooling is : " << cool << std::endl;
    pphot->PrintPhoton(ip);
    pphot->statp[ip] = DESTROYED;
  } else if ((std::isinf(dp1p)) || (std::isnan(dp1p))) {
    std::cout << "Warning: UpdateSourceTerms momentum change (k1p) is : "
              << dp1p << std::endl;
    pphot->PrintPhoton(ip);
  } else if ((std::isinf(dp2p)) || (std::isnan(dp2p))) {
    std::cout << "Warning: UpdateSourceTerms momentum change (k2p) is : "
              << dp2p << std::endl;
    pphot->PrintPhoton(ip);
  } else if ((std::isinf(dp3p)) || (std::isnan(dp3p))) {
    std::cout << "Warning: UpdateSourceTerms momentum change (k3p) is : "
              << dp3p << std::endl;
    pphot->PrintPhoton(ip);
    pphot->statp[ip] = DESTROYED;
    std::cout << "Warning: UpdateCooling cooling is : " << cool << std::endl;
  } else {
    int &i = pphot->i1p[ip];
    int &j = pphot->i2p[ip];
    int &k = pphot->i3p[ip];
    //sourceterms(MCRS0,k,j,i) -= cool; // BCM: We do not want the photons to cool here - handle cooling separately via explicit source terms
    pmcb->sourceterms(MCRS1,k,j,i) -= dp1p;
    pmcb->sourceterms(MCRS2,k,j,i) -= dp2p;
    pmcb->sourceterms(MCRS3,k,j,i) -= dp3p;
  }


}
//void EscapeCoords(MonteCarloBlock *pmcb, Photon *pphot, PhotonMover *pmover, int ip) {
//  if (pphot->statp[ip] == ESCAPED || pphot->statp[ip] == ABSORBED) {
//    pphot->user[6][ip] = pphot->x1p[ip];
//    pphot->user[7][ip] = pphot->x2p[ip];
//    pphot->user[8][ip] = pphot->x3p[ip];
//    pphot->user[9][ip] = pphot->k1p[ip];
//    pphot->user[10][ip] = pphot->k2p[ip];
//    pphot->user[11][ip] = pphot->k3p[ip];
//  }
//}

} // end namespace
