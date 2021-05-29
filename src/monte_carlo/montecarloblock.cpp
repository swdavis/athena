//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file montecarloblock.cpp
//  \brief implementation of functions in class MonteCarloBlock

#include <iostream>
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "montecarlo.hpp"
#include "photon.hpp"
#include "photonmover.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../globals.hpp"

// SWD: Eliminate these
#define MAXSCAT 10000

// constructor, initializes data structures and parameters

MonteCarloBlock::MonteCarloBlock(MeshBlock *pmb,  MCBlockSize *pblsize, MonteCarlo *pmc, 
                                 ParameterInput *pin) {
  
  pmy_mc = pmc;

  // Set related meshblock, coordinate
  pmy_block = pmb;

  // Set pointer to this monte carlo block in pmb if not NULL
  if (pmb != NULL) {
    pmb->pmy_mcb = this;
  }

  // Construct pointer to photon 
  pphoton  = new Photon(this,pmy_mc->nuser_var); // Currently one photon per block

  // Initialize to NULL and set below
  pmover = NULL;
  pcoord = NULL;

  // Initialize input parameters and flags
  weighted_absorption = pin->GetOrAddBoolean("montecarlo","abs_weight",true);
  //polarized = pin->GetOrAddBoolean("montecarlo","polarized",false);

  // get seed and intitialize randon number generator
  int rank = Globals::my_rank;
  int iseed = pmy_mc->iseed + rank *100;  // temporary solution
  pran = new MCRandom(iseed);

  next=NULL;

  // Set flags (initialized in MonteCarlo class)
  emission_meth = pmy_mc->emission_meth;
  absorption_meth = pmy_mc->absorption_meth;
  scattering_meth = pmy_mc->scattering_meth;
  boosts = pmy_mc->boosts;
  emission_array_flag = pmy_mc->emission_array_flag;
  moments_flag = pmy_mc->pmcout->moments; // set in mcoutput
  moments_comoving = pmy_mc->pmcout->moments_comoving;
  acceleration = pmy_mc->acceleration;
  time_acc = pmy_mc->time_acc;

  // *currently* assumes all block boundaries are physical
  SetBoundaryValues(pmy_mc->mc_bcs);

  // Initialize pbval after mcb_bcs is set
  pbval = new MCBoundaryValues(this,pin);

  // Setup output spectra
  Spectrum *pfirst = NULL, *plast;
  Spectrum *psmcout = pmy_mc->pmcout->pspec;
  // Loop over output spectra and make local equivalent for each
  while (psmcout != NULL) {
    pspec = new Spectrum(psmcout);
    //pspec = new Spectrum(pmy_mc->pmcout->pspec);
    if (pfirst == NULL)
      pfirst = pspec;
    else
      plast->next = pspec;
    plast = pspec;
    psmcout = psmcout->next;
  }
  pspec = pfirst;

  // Setup output photon list and trajectory list
  pphlist = pmy_mc->pmcout->pphlist;
  ptraj = pmy_mc->pmcout->ptraj;

  // set local mesh parameters to correspond to mesh block
  if (pmb != NULL) {
    is = pmb->is; ie = pmb->ie;
    js = pmb->js; je = pmb->je;
    ks = pmb->ks; ke = pmb->ke;
    nx1 = pmb->block_size.nx1;
    nx2 = pmb->block_size.nx2;
    nx3 = pmb->block_size.nx3;
  } else {
    if (pblsize == NULL) {
      std::stringstream msg;
      msg << "### FATAL ERROR Monte Carlo Block Constructor" << std::endl
          << "Both input Mesh Block and Block Size are NULL." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    } else {
      is = pblsize->is; ie = pblsize->ie;
      js = pblsize->js; je = pblsize->je;
      ks = pblsize->ks; ke = pblsize->ke;
      nx1 = pblsize->nx1;
      nx2 = pblsize->nx2;
      nx3 = pblsize->nx3;
      //pcoord = new MCCoord(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),acceleration);
    }
  }

  codetocgs_rho = 1.; codetocgs_vel = 1., codetocgs_tgas = 1.;  // default cgs for code units
  // SWD:  stepsize control needs to be modified
  stepsize = pin->GetOrAddReal("montecarlo","stepsize",1.0e-3);
  minweight = pin->GetOrAddReal("montecarlo","minweight",1.0e-20);

  // Flags for handling photon movement
  general_mover_flag = pin->GetOrAddBoolean("montecarlo","general_mover",false);
  boyerlindquist_flag = pin->GetOrAddBoolean("montecarlo","boyerlindquist",false);
  orthotet_flag = pin->GetOrAddBoolean("montecarlo", "orthotet", false);
  varystep_flag = pin->GetOrAddBoolean("montecarlo", "varystep", false);

  // Set up photon movement and initialization methods
  if (COORDINATE_SYSTEM == "cartesian") {
    GetZonePosition = GetZonePositionCartesian;
    if (general_mover_flag) {
      pmover = new GeneralMover(this);
      if (pmb != NULL)
        pcoord = new MCCartesian(pmb->pcoord,this);
      else
        pcoord = new MCCartesian(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                 acceleration);
    } else {
      pmover = new CartesianMover(this);
      if (pmb != NULL)
        pcoord = new MCCoord(pmb->pcoord,this);
      else
        pcoord = new MCCoord(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                             acceleration);
    }
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    GetZonePosition = GetZonePositionSphericalPolar;
    if (general_mover_flag) {
      pmover = new GeneralMover(this);
      if (pmb != NULL)
        pcoord = new MCSphericalPolar(pmb->pcoord,this);
      else
        pcoord = new MCSphericalPolar(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                      acceleration);
    } else {
      pmover = new SphericalPolarMover(this);
      if (pmb != NULL)
        pcoord = new MCCoord(pmb->pcoord,this);
      else
        pcoord = new MCCoord(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                             acceleration);
    } 
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    GetZonePosition = GetZonePositionCylindrical;
    pmover = new GeneralMover(this);
    if (pmb != NULL)
      pcoord = new MCCylindrical(pmb->pcoord,this);
    else
      pcoord = new MCCylindrical(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                 acceleration);
  } else if (COORDINATE_SYSTEM == "kerr-schild") {
    GetZonePosition = GetZonePositionSphericalPolar;//approximate
    pmover = new GeneralMover(this);
    if (boyerlindquist_flag) {
     if (pmb != NULL)
       pcoord = new MCBoyerLindquist(pmb->pcoord,this);
     else {
       pcoord = new MCBoyerLindquist(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                     acceleration);
       pcoord->SetSpin(pin->GetReal("coord", "a"));
       pcoord->SetMass(pin->GetReal("coord", "m"));
     }
    } else {
      if (pmb != NULL)
        pcoord = new MCKerrSchild(pmb->pcoord,this);
      else {
        pcoord = new MCKerrSchild(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                  acceleration);
        pcoord->SetSpin(pin->GetReal("coord", "a"));
        pcoord->SetMass(pin->GetReal("coord", "m"));
      }
    }
  } else if (COORDINATE_SYSTEM == "minkowski") {
    GetZonePosition = GetZonePositionCartesian;
    pmover = new GeneralMover(this);
    if (pmb != NULL)
      pcoord = new MCMinkowski(pmb->pcoord,this);
    else
      pcoord = new MCMinkowski(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                               acceleration);
  } else {
      std::stringstream msg;
      msg << "### ERROR in MonteCarloBlock constructor" << std::endl
          << COORDINATE_SYSTEM
          << "coordinates not currently supported with Monte Carlo" 
	  << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  // Set pcoord in pmover
  pmover->pcoord = pcoord;
  // Set absorption opacity
  if (absorption_meth == ABSUSER) {
    AbsorptionOpacity = pmy_mc->UserAbsorptionOpacity;
  } else if (absorption_meth == ABSNONE) {
    AbsorptionOpacity = NoOpacity;
  } else if (absorption_meth == ABSFF) {
    AbsorptionOpacity = FreeFreeAbsorptionOpacity;
  }
 
  // Set scattering opacity and method
  if (scattering_meth == SCATUSER) {
    ScatteringOpacity = pmy_mc->UserScatteringOpacity;
    Scatter = pmy_mc->UserScattering;
  } else if (scattering_meth == SCATNONE) {
    ScatteringOpacity = NoOpacity;
    Scatter = NoScatter;  // should not be called
    coherent_scattering = true;
  } else if (scattering_meth == SCATISO) {
    if (pmy_mc->polarized) {
      std::stringstream msg;
      msg << "### ERROR in MonteCarloBlock constructor" << std::endl
          << "Istropic scattering not suppored for polarized = " 
	  << pmy_mc->polarized << std::endl;
      throw std::runtime_error(msg.str().c_str());
    } else {
      ScatteringOpacity = ThomsonOpacity;
      Scatter = ScatterIsotropic;
      coherent_scattering = true;
    }
  } else if (scattering_meth == SCATTHOM) {
    ScatteringOpacity = ThomsonOpacity;
    if (pmy_mc->polarized) {
      Scatter = ScatterThomsonPolarized;
    } else
      Scatter = ScatterThomsonUnpolarized;
    coherent_scattering = true;
  } else if (scattering_meth == SCATCOMP) {
    int comptonio = pin->GetOrAddInteger("montecarlo","comptonio",1);
    GenerateComptonTable(comptonio);
    ScatteringOpacity = ComptonOpacity;
    if (pmy_mc->polarized) {
      Scatter = ScatterComptonPolarized;
    } else {
      Scatter = ScatterComptonUnpolarized;
    }
    coherent_scattering = false;
  } else if (scattering_meth == SCATRES) {
    ScatteringOpacity = ResonanceLineOpacity;
    if (pmy_mc->polarized) {
      std::stringstream msg;
      msg << "### ERROR in MonteCarloBlock constructor" << std::endl
          << "Lyman alpha scattering not suppored for polarized = " 
	  << pmy_mc->polarized << std::endl;
      throw std::runtime_error(msg.str().c_str());
    } else {
      Scatter = ScatterResonanceLine;
      coherent_scattering = false;
    }
  }


  // Allocate (/initialize) variable arrays needed for evolution/output
  int ncells1 = nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (nx2 > 1) ncells2 = nx2 + 2*(NGHOST);
  if (nx3 > 1) ncells3 = nx3 + 2*(NGHOST);
  rho.NewAthenaArray(ncells3,ncells2,ncells1);
  tgas.NewAthenaArray(ncells3,ncells2,ncells1);
  if (boosts) vel.NewAthenaArray(3,ncells3,ncells2,ncells1);
  // moments is 1 (Er) + 3 (Fr) + 9 (Pr) + 1 (Eave) + 1 (net cool)
  if (moments_flag) moments.NewAthenaArray(NMOM,ncells3,ncells2,ncells1);
  if (emission_array_flag) emission.NewAthenaArray(ncells3,ncells2,ncells1);
  if (acceleration && !(coherent_scattering)) {
    planck_opacity.NewAthenaArray(ncells3,ncells2,ncells1);
    planck_inv_opacity.NewAthenaArray(ncells3,ncells2,ncells1);
  }

  // Create user monte carlo block data
  InitUserMonteCarloBlockData(pin);

}


// destructor

MonteCarloBlock::~MonteCarloBlock() {

  delete pphoton;
  delete pmover;
  delete pran;
  delete pspec;
  delete pphlist;
  delete ptraj;

  rho.DeleteAthenaArray();
  tgas.DeleteAthenaArray();
  if (boosts) vel.DeleteAthenaArray();
  if (moments_flag) moments.DeleteAthenaArray();
  if (emission_array_flag) emission.DeleteAthenaArray();
  if (acceleration && !(coherent_scattering)) {
    planck_opacity.DeleteAthenaArray();
    planck_inv_opacity.DeleteAthenaArray();
  }
}

/*void MonteCarloBlock::DefaultGetTemperature() {

  Real rideal = 8.314e7;
  Hydro* phydro = pmy_block->phydro;

   // MonteCarloBlock ranges should always match MeshBlock ranges
  int il = is; int iu = ie;
  int jl = js; int ju = je;
  int kl = ks; int ku = ke;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        tgas(k,j,i) = phydro->w(IEN,k,j,i)/phydro->w(IDN,k,j,i)/rideal;

      }}}

      }*/

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::RayTracePhotons()
//  \brief Integrate photons to termination condtion without scattering

void MonteCarloBlock::RayTracePhotons(int nphot) {

    Real const to_comv = 1.0;
    Real const to_eulr = -1.0; 
    int nscat = 0, nesc = 0, nabs = 0;
    int nprop = (nphot > nphremain) ? nphremain : nphot;
 
    for (int i=0; i<nprop; ++i) {
  
      // user definied photon initialization
      InitializePhoton(pphoton);
 
      // Photon initialized in coordinate frame
      // move photon until  stopping condition
      pmover->Move(pphoton);
      if (ptraj != NULL) ptraj->CompleteTrajectory();
      // User defined completion work
      FinalizePhoton(pphoton);
      if (pphoton->status == ESCAPED) {
	// loop over spectra and update
	Spectrum *pspect = pspec;
	while (pspect != NULL) {
          pspect->UpdateSpectrum(pphoton);
	  pspect = pspect->next;
	}
        if (pphlist != NULL) {
          pphlist->AddPhoton(pphoton);
        }
	nesc++;
      } else
	nabs++;
    }

    std::cout  << "nesc, nabs: " << nesc << ' ' << nabs << ' ' << Globals::my_rank 
               << std::endl;
    std::cout << "nscat: " << nscat << ' ' << Globals::my_rank << std::endl;


    return;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::TransferPhotons()
//  \brief perform radiation transfer nphtot photons

void MonteCarloBlock::TransferPhotons(int nphot) {

  Real const to_comv = 1.0;
  Real const to_eulr = -1.0; 
  int nscat = 0, nesc = 0, nabs = 0;
  int nprop = (nphot > nphremain) ? nphremain : nphot;
  for(int i=0; i<nprop; ++i) {

    // user definied photon initialization
    InitializePhoton(pphoton);
    
    // Lorentz transform E, k to Eulerian frame and update opacities.
    if (boosts) {
      // SWD: Change this to a function pointer ?
      if (orthotet_flag) {
        TetradTransform(pphoton, to_eulr);  
      } else {
        LorentzTransform(pphoton,to_eulr);
      }
    }
    if (moments_flag)
      UpdateCooling(pphoton,0.,0.);

    // move photon to next scattering/absorption or to boundary
    pmover->Move(pphoton);
    int iscat = 0;
    Real xmax = 0.;
    while (pphoton->status == EVOLVING) {
      
      // Account for absorption
      Real weight0 = pphoton->weight;
      if (weighted_absorption) {
        pphoton->weight *= (pphoton->sct_coef / (pphoton->sct_coef+pphoton->abs_coef));
        if(pphoton->weight <= minweight) {
          pphoton->status = ABSORBED;
        }
      } else {
        if (pran->uniform() > (pphoton->sct_coef / (pphoton->sct_coef+pphoton->abs_coef)) )
          pphoton->status = ABSORBED;
      }
      if (moments_flag)
        UpdateCooling(pphoton,0.,weight0);
   
      // Scatter the photon packet
      if (pphoton->status == EVOLVING) {
        Real e_pre_scat = pphoton->energy;
	// Lorentz transform to comoving frame for scattering
	if (boosts) {
          if (orthotet_flag) {
            TetradTransform(pphoton, to_comv);
          } else { 
            LorentzTransform(pphoton,to_comv);
          }
	}
        Scatter(this,pphoton);
	iscat++;
      
	if (iscat %  MAXSCAT == 0) {
	  // Check for possible infinite loop due to NaN in photon
	  if (pphoton->IsNanPhoton()) {
	    pphoton->status = DESTROYED;
	    std::cout << "Warning: IsNanPhoton() returned true, photon destroyed" 
		      << std::endl;
	    pphoton->PrintPhoton();
	  }
	}
	// Update the absorption and scattering extinction coefficients
	// with the new energy.
	if (!coherent_scattering) {
	  pphoton->abs_coef = AbsorptionOpacity(this,pphoton);
	  pphoton->sct_coef = ScatteringOpacity(this,pphoton);
	}
	// Lorentz transform to Eulerian frame and shift opacities
	if (boosts) {
          if (orthotet_flag) {
            TetradTransform(pphoton, to_eulr);
          } else {
            LorentzTransform(pphoton,to_eulr);
          }
        }
        if (moments_flag)
          UpdateCooling(pphoton,e_pre_scat,0.);
      }

      // move photon to next scattering/absorption or to boundary
      pmover->Move(pphoton);
    }
    if (ptraj != NULL) ptraj->CompleteTrajectory();
    if (pphoton->status == ESCAPED) {
      // User defined completion work
      FinalizePhoton(pphoton);
      // loop over spectra and update
      Spectrum *pspect = pspec;
      while (pspect != NULL) {
        pspect->UpdateSpectrum(pphoton);
	pspect = pspect->next;
      }
      if (pphlist != NULL) {
        pphlist->AddPhoton(pphoton);
      }
      nesc++;
    } else
      nabs++;
    nscat += iscat;
  }
  
  nphdone += nprop;
  
  std::cout  << "nesc, nabs: " << nesc << ' ' << nabs << ' ' << Globals::my_rank << std::endl;
  std::cout << "nscat: " << nscat << ' ' << Globals::my_rank << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::LorentzTransform(Photon *pphot, const Real sign)
//  \brief Lorentz transform photon packet
//
// Does not transform stokes vectors but this seems
// to be correct -- the plane of polarization is invariant under lorentz
// transformation as discussed in Cocke & Holm (1972) Nature letter.
// weight is not transformed either as weight represents number of photons
// in the packet which is invariant.
// to_comv: sign = 1.0;
// to_eulr: sign = -1.0; 

void MonteCarloBlock::LorentzTransform(Photon *pphot, const Real sign) {

  Real &k1 = pphot->k[0];
  Real &k2 = pphot->k[1];
  Real &k3 = pphot->k[2];
  int i1 = pphot->i1, i2 = pphot->i2, i3 = pphot->i3;
  
  Real beta[3];
  for (int i=0; i<3; ++i) {
    beta[i] = sign * vel(i,i3,i2,i1) / 2.9979e10;
  }
  Real beta2= SQR(beta[0]) + SQR(beta[1]) + SQR(beta[2]);

  if(beta2 > 0.) {
    // SWD: pretabulate gamma for each zone?
    Real gamma = 1. / sqrt(1. - beta2); // assumes v^2 < c^2 checked elsewhere
    Real bdk = k1 * beta[0] + k2 * beta[1] + k3 * beta[2];
    Real gonembdk = gamma * (1. - bdk);
    //Real aber = (gamma-1.) * bdk / beta2 - gamma;
    Real aber = gamma*(1.-gamma*bdk/(gamma+1.));
    
    pphot->energy *= gonembdk;
    //printf("%g %g %g %g\n",gonembdk,gamma,bdk,pphot->energy/1.6021772e-12);
    Real kz = k3;
    k1 = (k1 - aber * beta[0]) / gonembdk;
    k2 = (k2 - aber * beta[1]) / gonembdk;
    k3 = (k3 - aber * beta[2]) / gonembdk;
    //printf("%g;",k3);
    //if ((k3 < 1./32)&&(k3 > 0.)) printf("kz: %g %g %g %g %g\n",kz,k3,gonembdk,aber*beta[2],aber);
    
    // Transform opacities
    // Must be performed even when transforming to comoving frame because inverse process 
    // is performed to go back to Eulerian frame in cases where scattering is coherent
    pphot->abs_coef /= gonembdk;
    pphot->sct_coef /= gonembdk;
    
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real MonteCarloBlock::LorentzTransformFrequencyShift(Photon *pphot)
//  \brief Returns frequency shift from Lorentz Trransformation

Real MonteCarloBlock::LorentzTransformFrequencyShift(Photon *pphot) {

  Real k1 = pphot->k[0];
  Real k2 = pphot->k[1];
  Real k3 = pphot->k[2];
  int i1 = pphot->i1, i2 = pphot->i2, i3 = pphot->i3;
  
  Real beta[3];
  for (int i=0; i<3; ++i) {
    beta[i] = vel(i,i3,i2,i1) / 2.9979e10;
  }
  Real beta2= SQR(beta[0]) + SQR(beta[1]) + SQR(beta[2]);

  Real gonembdk;
  if(beta2 > 0.) {
    Real gamma = 1. / sqrt(1. - beta2); // assumes v^2 < c^2 checked elsewhere
    Real bdk = k1 * beta[0] + k2 * beta[1] + k3 * beta[2];
    gonembdk = gamma * (1. - bdk);
  } else {
    gonembdk = 1.;
  }
  // Always called from lab frame so returna nu'/nu
  return gonembdk;
  
}

// SWD: This is an untested modification of Eric's original method
//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::TetradTransform(Photon *pphot, const Real sign)
//  \brief Tetrad transform photon packet
//
// Does not transform stokes vectors..
// weight is not transformed either as weight represents number of photons
// in the packet which is invariant.
// to_comv: sign = 1.0;
// to_eulr: sign = -1.0; 

void MonteCarloBlock::TetradTransform(Photon *pphot, const Real sign) {

  // Get velocity of cell
  int i1 = pphot->i1, i2 = pphot->i2, i3 = pphot->i3;
  Real beta[3];
  for (int i=0; i<3; ++i) {
    beta[i] = sign * vel(i,i3,i2,i1) / 2.9979e10;
  }
  Real beta2= SQR(beta[0]) + SQR(beta[1]) + SQR(beta[2]);
  Real gamma = 1. / sqrt(1. - beta2);

  // Define velocity four vector for cell
  Real ucon[NCOORD];
  ucon[IMC0] = gamma;
  ucon[IMC1] = gamma * beta[0];
  ucon[IMC2] = gamma * beta[1];
  ucon[IMC3] = gamma * beta[2];
  

  // get metric values for current position
  Real gcov[NCOORD][NCOORD];
  pcoord->Metric(pphot->x, gcov);

  // create tetrad basis
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  ConstructTetrad(ucon, gcov, econ, ecov);

  if (sign > 0) { // tranforming to comoving frame

    // SWD: mirrors Lorentz transformation but redundant -> should use CoordinateToTetrad
    Real kdotu;
    for (int i = 0; i < NCOORD; i++) 
      kdotu += pphot->k[i] * ucon[i]; // pphot->k in coordinate frame
    Real energy_shift = - pphot->k[IMC0] / kdotu; 


    Real kcopy[NCOORD];
    for (int i = 0; i < NCOORD; i++)
      kcopy[i] = pphot->k[i];
    CoordinateToTetrad(kcopy, pphot->k, ecov); // updates pphot->k 

    // transform energy and extinction coefficients
    pphot->energy *= energy_shift; 
    pphot->abs_coef *= energy_shift;
    pphot->sct_coef *= energy_shift;


  } else { // transforming to coordinate frame

    Real kcopy[NCOORD];
    for (int i = 0; i < NCOORD; i++)
      kcopy[i] = pphot->k[i];
    TetradToCoordinate(kcopy, pphot->k, econ); // updates pphot->k

    // Eric's implementation -- needs to be updated
    Real kdotu = DotVec(pphot->k, ucon, gcov);
    if (fabs(kdotu) < 1.0e-30) {
      printf("warning: kdotu = %g\n", kdotu);
      kdotu = 1.0e-30;
    }
    Real energy_shift = - kdotu / pphot->k[IMC0]; // new calculation
    
    // transform energy and opacities
    pphot->energy *= energy_shift;
    pphot->abs_coef *= energy_shift;
    pphot->sct_coef *= energy_shift;
    
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl)
//  \brief add contribution to radiation moments in current zone

void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl) {
  
  // SWD: needs to be modifed for non Cartesian kvectors
  Real k1,k2,k3;
  Real energy, abs_coef, step;
  if (moments_comoving) {
    // boost relevant quanitities to comoving frame
    k1 = pphot->k[0];
    k2 = pphot->k[1];
    k3 = pphot->k[2];
    energy = pphot->energy;
    int i1 = pphot->i1, i2 = pphot->i2, i3 = pphot->i3;
    Real beta[3];
    for (int i=0; i<3; ++i) {
      beta[i] = vel(i,i3,i2,i1) / 2.99792458e10;
    }
    Real beta2= SQR(beta[0]) + SQR(beta[1]) + SQR(beta[2]);
    
    if(beta2 > 0.) {
      Real gamma = 1. / sqrt(1. - beta2); // assumes v^2 < c^2 checked elsewhere
      Real bdk = k1 * beta[0] + k2 * beta[1] + k3 * beta[2];
      Real gonembdk = gamma * (1. - bdk);
      Real aber = gamma*(1.-gamma*bdk/(gamma+1.));
    
      energy *= gonembdk;
      k1 = (k1 - aber * beta[0]) / gonembdk;
      k2 = (k2 - aber * beta[1]) / gonembdk;
      k3 = (k3 - aber * beta[2]) / gonembdk;
      abs_coef = pphot->abs_coef /= gonembdk;
      step = dl * gonembdk;
    }
  } else {
    // Use eulerian values
    k1 = pphot->k[0];
    k2 = pphot->k[1];
    k3 = pphot->k[2];
    energy = pphot->energy;
    abs_coef = pphot->abs_coef;
    step = dl;
  }

  // Weight moments by time spent in domain
  Real weight = pphot->weight * energy * step / 2.99792458e10;
  if ((isinf(weight)) || (isnan(weight))) {
    std::cout << "Warning: UpdateMoments weight is : " << weight << std::endl;
  } else {
    // Higher order moments are weighted by curvalinear coorindates k
    Real weight1 = weight * k1;
    Real weight2 = weight * k2;
    Real weight3 = weight * k3;

    int i = pphot->i1;
    int j = pphot->i2;
    int k = pphot->i3;

    // Add contribution to corresponding moments
    // Energy density
    // SWD: Modify this appropriately
    //if (general_mover_flag) 
    //  weight *= pphot->k[IMC0];
    moments(MCIER,k,j,i) += weight;
    // Flux
    moments(MCIFR1,k,j,i) += weight1 * 2.99792458e10;
    moments(MCIFR2,k,j,i) += weight2 * 2.99792458e10;
    moments(MCIFR3,k,j,i) += weight3 * 2.99792458e10;
    // Radiation Pressure
    //Real weightp = weight1 * pphot->k[0];
    Real weightp = weight1 * 2.99792458e10;
    moments(MCIPR11,k,j,i) += weightp;
    weightp = weight2 * k2;
    moments(MCIPR22,k,j,i) += weightp;
    weightp = weight3 * k3;
    moments(MCIPR33,k,j,i) += weightp;
    weightp = weight1 * k2;
    moments(MCIPR12,k,j,i) += weightp;
    weightp = weight1 * k3;
    moments(MCIPR13,k,j,i)  += weightp;
    weightp = weight2 * k3;
    moments(MCIPR23,k,j,i) += weightp;
    // Photon mean energy
    moments(MCIEN,k,j,i) += weight * energy;
    // Jmean opacity
    moments(MCIKJ,k,j,i) += weight * abs_coef;
  }

}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::NormalizeMoments(bool normalize)
//  \brief (un)normalized moments for output and copy symmetric elements

void MonteCarloBlock::NormalizeMoments(bool normalize) {

  // Normalize all moments by number of photons emitted
  Real normall = static_cast<Real>(nphdone);
 
  if (normalize) {
   // Normalize energy density weighted averages first
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (moments(MCIER,k,j,i) != 0.) {
            moments(MCIKJ,k,j,i) /= moments(MCIER,k,j,i);
            moments(MCIEN,k,j,i) /= moments(MCIER,k,j,i);
          }
        }}}
    // Normalize remaining moments by volume and global norm (counts)
    for (int n=0; n<11; ++n) {
      Real norm = normall;
      for (int k=ks; k<=ke; ++k) {
	for (int j=js; j<=je; ++j) {
	  for (int i=is; i<=ie; ++i) {
	    moments(n,k,j,i) /= (pcoord->vol(k,j,i) * norm);
	  }}}
    }
    // Copy normalized moments to symmetric elements
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
	for (int i=is; i<=ie; ++i) {
	  moments(MCIPR21,k,j,i) = moments(MCIPR12,k,j,i);
	  moments(MCIPR31,k,j,i) = moments(MCIPR13,k,j,i);
	  moments(MCIPR32,k,j,i) = moments(MCIPR23,k,j,i);
	}}}
  } else {
    // Undo normalization for continuing evolution
    for (int n=0; n<11; ++n) {
      Real norm = normall;
      for (int k=ks; k<=ke; ++k) {
	for (int j=js; j<=je; ++j) {
	  for (int i=is; i<=ie; ++i) {
	    moments(n,k,j,i) *= (pcoord->vol(k,j,i) * norm);
	  }}}
    }
    // Unnormalize energy density weighted averages after moments
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (moments(MCIER,k,j,i) != 0.) {
            moments(MCIKJ,k,j,i) *= moments(MCIER,k,j,i);
            moments(MCIEN,k,j,i) *= moments(MCIER,k,j,i);
          }
        }}}
 }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::ResetMoments()
//  \brief set moments to zero on origin blocks

void MonteCarloBlock::ResetMoments() {

    // set moments to zero
  for (int n=0; n<11; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
	for (int i=is; i<=ie; ++i) {
	  moments(n,k,j,i) = 0.;
	}}}
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateCooling(Photon *pphot, Real energy0, Real weight0)
//  \brief compute net photon cooling rate

void MonteCarloBlock::UpdateCooling(Photon *pphot, Real energy0, Real weight0) {
  
  Real cool = (pphot->weight - weight0) * (pphot->energy - energy0);
  if ((isinf(cool)) || (isnan(cool))) {
    std::cout << "Warning: UpdateCooling cooling is : " << cool << std::endl;
    pphot->PrintPhoton();
  } else {
    int i = pphot->i1;
    int j = pphot->i2;
    int k = pphot->i3;
    moments(MCINET,k,j,i) -= cool;
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs)
//  \brief set boundary values on monte carlo block

void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs) {

  // set x1 boundaries
  mcb_bcs[INNER_X1] = input_bcs[INNER_X1];
  mcb_bcs[OUTER_X1] = input_bcs[OUTER_X1];

  // set x2 boundaries
  mcb_bcs[INNER_X2] = input_bcs[INNER_X2];
  mcb_bcs[OUTER_X2] = input_bcs[OUTER_X2];

  // set x3 boundaries
  mcb_bcs[INNER_X3] = input_bcs[INNER_X3];
  mcb_bcs[OUTER_X3] = input_bcs[OUTER_X3];
  
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs)
//  \brief set boundary values on monte carlo block

/*void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs) {

  // set x1 boundaries
  if (pmy_block->pbval->block_bcs[INNER_X1] == BLOCK_BNDRY) {
    mcb_bcs[INNER_X1] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[INNER_X1] = input_bcs[INNER_X1];
  }
  if (pmy_block->pbval->block_bcs[OUTER_X1] == BLOCK_BNDRY) {
    mcb_bcs[OUTER_X1] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[OUTER_X1] = input_bcs[OUTER_X1];
  }

  // set x2 boundaries
  if (pmy_block->pbval->block_bcs[INNER_X2] == BLOCK_BNDRY) {
    mcb_bcs[INNER_X2] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[INNER_X2] = input_bcs[INNER_X2];
  }
  if (pmy_block->pbval->block_bcs[OUTER_X2] == BLOCK_BNDRY) {
    mcb_bcs[OUTER_X2] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[OUTER_X2] = input_bcs[OUTER_X2];
  }

  // set x3 boundaries
  if (pmy_block->pbval->block_bcs[INNER_X3] == BLOCK_BNDRY) {
    mcb_bcs[INNER_X3] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[INNER_X3] = input_bcs[INNER_X3];
  }
  if (pmy_block->pbval->block_bcs[OUTER_X3] == BLOCK_BNDRY) {
    mcb_bcs[OUTER_X3] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[OUTER_X3] = input_bcs[OUTER_X3];
  }
  }*/

