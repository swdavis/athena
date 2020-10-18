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

#define MINWEIGHT 1.0e-30
#define MAXSCAT 10000
//#define MINWEIGHT 1.0e-30

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
  zone_weight_flag = pin->GetOrAddBoolean("montecarlo","zone_weight",true);
  weighted_absorption = pin->GetOrAddBoolean("montecarlo","abs_weight",true);
  //polarized = pin->GetOrAddBoolean("montecarlo","polarized",false);

  // get seed and intitialize randon number generator
  int rank = Globals::my_rank;
  int iseed = pmy_mc->iseed + rank *100;  // temporary solution
  pran = new MCRandom(iseed);

  next=NULL;

  // Set energy range in ergs (input assumed in eV)
  Real everg = 1.6021772e-12;
  emin = everg * pin->GetReal("montecarlo","emin");
  emax = everg * pin->GetReal("montecarlo","emax");
  elog = log10(emax/emin);
  eminlog = log10(emin);

  // Set flags (initialized in MonteCarlo class)
  emission_meth = pmy_mc->emission_meth;
  absorption_meth = pmy_mc->absorption_meth;
  scattering_meth = pmy_mc->scattering_meth;
  boosts = pmy_mc->boosts;
  emission_array_flag = pmy_mc->emission_array_flag;
  moments_flag = pmy_mc->pmcout->moments; // set in mcoutput
  acceleration = pmy_mc->acceleration;
  time_acc = pmy_mc->time_acc;

  // *currently** assumes all block boundaries are physical
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

  // Setup output photon list
  pphlist = pmy_mc->pmcout->pphlist;

  // set local mesh parameters to correspond to mesh block
  if (pmb != NULL) {
    is = pmb->is; ie = pmb->ie;
    js = pmb->js; je = pmb->je;
    ks = pmb->ks; ke = pmb->ke;
    nx1 = pmb->block_size.nx1;
    nx2 = pmb->block_size.nx2;
    nx3 = pmb->block_size.nx3;
    pcoord = new MCCoord(pmb->pcoord,this);
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
      pcoord = new MCCoord(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),acceleration);
    }
  }

  codetocgs_rho = 1.0; codetoc_vel = 1.0;  // default cgs for code units

  stepsize = pin->GetOrAddReal("montecarlo","stepsize",1.0e-3);
  // SWD: a needs to go
  a = pin->GetOrAddReal("montecarlo", "spin", 0.0);
  velocity = pin->GetOrAddReal("problem", "velocity", 0.0);


  // Set up photon movement and initialization methods
  general_mover_flag = pin->GetOrAddBoolean("montecarlo","general_mover",false);
  kerrschild_flag = pin->GetOrAddBoolean("montecarlo","kerrschild",false);
  boyerlindquist_flag = pin->GetOrAddBoolean("montecarlo","boyerlindquist",false);
  orthotet_flag = pin->GetOrAddBoolean("montecarlo", "orthotet", false);
  varystep_flag = pin->GetOrAddBoolean("montecarlo", "varystep", false);

 if (general_mover_flag) {
    pmover = new GeneralMover(this);
    if (COORDINATE_SYSTEM == "cartesian") {
      GetZonePosition = GetZonePositionCartesian;
      if (kerrschild_flag) {
	Connection = Connect_KerrSchild;
	Metric = Metric_KerrSchild;
      }
      else {
	Connection = Connect_Cartesian;
	Metric = Metric_Cartesian;
      }
    } else if ((COORDINATE_SYSTEM == "spherical_polar") ||
               (COORDINATE_SYSTEM == "kerr-schild")){
      GetZonePosition = GetZonePositionSphericalPolar;
      if (kerrschild_flag) {
	Connection = Connect_KerrSchild;
	Metric = Metric_KerrSchild;
      } else if (boyerlindquist_flag) {
	Connection = Connect_BoyerLindquist;
	Metric = Metric_BoyerLindquist;
      } else {
	Connection = Connect_SphericalPolar;
	Metric = Metric_SphericalPolar;
      }
    } else if (COORDINATE_SYSTEM == "cylindrical") {
      GetZonePosition = GetZonePositionCylindricalGR;
      Connection = Connect_Cylindrical;
      Metric = Metric_Cylindrical;
    }
  } else {
    if (COORDINATE_SYSTEM == "cartesian") {
      pmover = new CartesianMover(this);
      GetZonePosition = GetZonePositionCartesian;
    } else if (COORDINATE_SYSTEM == "spherical_polar") {
      pmover = new SphericalPolarMover(this);
      GetZonePosition = GetZonePositionSphericalPolar;
    } else if (COORDINATE_SYSTEM == "cylindrical") {
      std::stringstream msg;
      msg << "### ERROR in MonteCarloBlock constructor" << std::endl
          << "cylindircal coordinates only suppored with general photon mover" 
	  << std::endl;
      throw std::runtime_error(msg.str().c_str());
      //pmover = new CylindricalMover(this);
      //GetZonePosition = GetZonePositionCylindrical;
    }
  }

  // Set absorption opacity
  if (absorption_meth == ABSUSER) {
    AbsorptionOpacity = NULL;
  } else if (absorption_meth == ABSNONE) {
    AbsorptionOpacity = NoOpacity;
  } else if (absorption_meth == ABSFF) {
    AbsorptionOpacity = FreeFreeAbsorptionOpacity;
  }
 
  // Set scattering opacity and method
  if (scattering_meth == SCATUSER) {
    ScatteringOpacity = NULL;
    Scatter = NULL;
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
  }

  // Allocate (/initialize) variable arrays needed for evolution/output
  int ncells1 = nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (nx2 > 1) ncells2 = nx2 + 2*(NGHOST);
  if (nx3 > 1) ncells3 = nx3 + 2*(NGHOST);
  rho.NewAthenaArray(ncells3,ncells2,ncells1);
  tgas.NewAthenaArray(ncells3,ncells2,ncells1);
  if (boosts) vel.NewAthenaArray(3,ncells3,ncells2,ncells1);
  if (moments_flag) moments.NewAthenaArray(14,ncells3,ncells2,ncells1);
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
      // SWD: Change this to a function pointer
      if (orthotet_flag) {
        TetradTransform(pphoton, to_eulr);  
      } else {
        LorentzTransform(pphoton,to_eulr);
      }
    }

    // move photon to next scattering/absorption or to boundary
    pmover->Move(pphoton);
    int iscat = 0;
    while (pphoton->status == EVOLVING) {
      
      // Account for absorption
      if (weighted_absorption) {
        pphoton->weight *= (pphoton->sct_coef / (pphoton->sct_coef+pphoton->abs_coef));
        if(pphoton->weight <= MINWEIGHT)
          pphoton->status = DESTROYED;
      } else {
        if (pran->uniform() > (pphoton->sct_coef / (pphoton->sct_coef+pphoton->abs_coef)) )
          pphoton->status = DESTROYED;
      }
      
      // Scatter the photon packet
      if (pphoton->status == EVOLVING) {
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
            LorentzTransform(pphoton,to_eulr);
          } else {
            TetradTransform(pphoton, to_eulr);
          }
        }
      }

      // move photon to next scattering/absorption or to boundary
      pmover->Move(pphoton);
    }

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
    Real gamma = 1. / sqrt(1. - beta2); // assumes v^2 < c^2 checked elsewhere
    Real bdk = k1 * beta[0] + k2 * beta[1] + k3 * beta[2];
    Real gonembdk = gamma * (1. - bdk);
    //Real aber = (gamma-1.) * bdk / beta2 - gamma;
    Real aber = gamma*(1.-gamma*bdk/(gamma+1.));
    
    pphot->energy *= gonembdk;
    //printf("%g %g %g %g %g\n",pphot->eweight,gonembdk,gamma,bdk,pphot->energy/1.6021772e-12);
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

  // Construct the orthonormal tetrad
  Real ucon[NCOORD];
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  Real kcopy[NCOORD];
  Real gcov[NCOORD][NCOORD];
  Real energy_shift;
  Real kdotu = 0.;
  
  Metric(pphot->x, gcov);

  if ((kerrschild_flag) or (boyerlindquist_flag)) {
    Real r = pphot->x[IMC1];
    Real omega = 1.0/(pow(r, 3./2.) + a); // circular velocity 

    ucon[IMC0] = sqrt(-1.0/(gcov[IMC0][IMC0] + 2.*gcov[IMC0][IMC3]*omega +
			    SQR(omega)*gcov[IMC3][IMC3]));
    Real del = SQR(r)+SQR(a)-2.*r;
    Real aeq = SQR(r)+SQR(a);
    aeq = SQR(aeq)-SQR(a)*del;
    Real enu = sqrt(del/aeq)*r;
    Real vg = (SQR(r)-2.*a*sqrt(r)+SQR(a))/sqrt(del)*omega;
    Real gam = 1./sqrt(1.-SQR(vg));
    Real comp = gam/enu;
    //printf("u0: %g %g\n",comp,ucon[IMC0]);
    ucon[IMC1] = 0.;
    ucon[IMC2] = 0.;
    ucon[IMC3] = (ucon[IMC0])*omega;
  } else {
    Real beta2 = SQR(velocity);
    Real gamma = 1. / sqrt(1. - beta2);

    ucon[IMC0] = gamma;
    ucon[IMC1] = gamma * velocity;
    ucon[IMC2] = 0.;
    ucon[IMC3] = 0.;
  }
  //printf("ucon: %g %g %g %g\n", ucon[IMC0], ucon[IMC1], ucon[IMC2], ucon[IMC3]);

  // create tetrad basis
  ConstructTetrad(ucon, gcov, econ, ecov);

  if (sign > 0) { // tranforming to comoving frame

    for (int i = 0; i < NCOORD; i++) 
      kdotu += pphot->k[i] * ucon[i]; // pphot->k in coordinate frame

    energy_shift = - pphot->k[IMC0] / kdotu; 

    for (int i = 0; i < NCOORD; i++)
      kcopy[i] = pphot->k[i];
    CoordinateToTetrad(kcopy, pphot->k, ecov); // updates pphot->k 

    Real klow[NCOORD];
    ConToCov(kcopy,klow,gcov);
    Real r = pphot->x[IMC1];
    Real omega = 1.0/(pow(r, 3./2.) + a); // circular velocity 
    Real del = SQR(r)+SQR(a)-2.*r;
    Real aeq = SQR(r)+SQR(a);
    aeq = SQR(aeq)-SQR(a)*del;
    Real enu = sqrt(del/aeq)*r;
    Real vg = (SQR(r)-2.*a*sqrt(r)+SQR(a))/sqrt(del)*omega;
    Real gam = 1./sqrt(1.-SQR(vg));
    Real comp1 = (1.+omega*klow[IMC3]/klow[IMC0])*gam/enu;
    Real comp2 = -DotVec(kcopy,ucon,gcov);
    //printf("kt: %g %g %g\n",pphot->k[IMC0],comp1,comp2);
     
    //  update pphot->energy and coefficients
    pphot->energy *= energy_shift; // now in tetrad frame
    pphot->abs_coef *= energy_shift;
    pphot->sct_coef *= energy_shift;
    //printf("energy_shift to comoving frame: %g\n", energy_shift);

  } else { // transforming to coordinate frame

    TetradToCoordinate(kcopy, pphot->k, econ); // updates pphot->k

    /*for (int i = 0; i < NCOORD; i++) 
      kdotu += pphot->k[i] * ucon[i]; // pphot->k in coordinate frame*/
    kdotu = DotVec(pphot->k, ucon, gcov);

    if (fabs(kdotu) < 1.0e-30) {
      printf("warning: kdotu = %g\n", kdotu);
      kdotu = 1.0e-30;
    }

    //energy_shift = pphot->k[IMC0] / kdotu;
    energy_shift = - kdotu / pphot->k[IMC0]; // new calculation
    //printf("energy_shift to coordinate frame: %g\n", energy_shift);
    /*printf("energy_shift: %g  k^t_tet = %g  k^t_coord = %g  tet/coord = %g\n",
      energy_shift, kcopy[IMC0], pphot->k[IMC0], kcopy[IMC0]/pphot->k[IMC0]);*/
    
    // update pphot->energy and opacities
    pphot->energy *= energy_shift;
    pphot->abs_coef *= energy_shift;
    pphot->sct_coef *= energy_shift;
    
  }

}

//----------------------------------------------------------------------------------------
//! \fn Real MonteCarloBlock::TetradTransformFrequencyShift(Photon *pphot)
//  \brief return energy shift

Real MonteCarloBlock::TetradTransformFrequencyShift(Photon *pphot) {

  // Construct the orthonormal tetrad
  Real r = pphot->x[IMC1];
  Real ucon[NCOORD];
  Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
  Real kcopy[NCOORD];
  Real gcov[NCOORD][NCOORD];
  Real energy_shift;
  Real omega;
  Real kdotu = 0.;
  Real gamma, beta2;

  omega = 1.0/(pow(r, 3./2.) + a);
  
  Metric(pphot->x, gcov);

  if (COORDINATE_SYSTEM == "spherical_polar") {
    ucon[IMC0] = sqrt(-1.0/(gcov[IMC0][IMC0] + 2.*gcov[IMC0][IMC3]*omega +
			    SQR(omega)*gcov[IMC3][IMC3]));
    ucon[IMC1] = 0.;
    ucon[IMC2] = 0.;
    ucon[IMC3] = (ucon[IMC0])*omega;
  } else if (COORDINATE_SYSTEM == "cartesian") {
    beta2 = SQR(velocity);
    gamma = 1. / sqrt(1. - velocity);
    ucon[IMC0] = gamma;
    ucon[IMC1] = gamma * velocity;
    ucon[IMC2] = gamma * 0.;
    ucon[IMC3] = gamma * 0.;
  }
  
  for (int i = 0; i < NCOORD; i++) 
    kdotu += pphot->k[i] * ucon[i];

  energy_shift = kdotu / pphot->k[IMC0];

  return energy_shift;

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl)
//  \brief add contribution to radiation moments in current zone

void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl) {
  
  Real weight = pphot->eweight * pphot->weight * pphot->energy * dl;
  if ((isinf(weight)) || (isnan(weight))) {
    std::cout << "Warning: UpdateMoments weight is : " << weight << std::endl;
  } else {
    // Higher order moments are weighted by curvalinear coorindates k
    Real weight1 = weight * pphot->k[0];
    Real weight2 = weight * pphot->k[1];
    Real weight3 = weight * pphot->k[2];

    int i = pphot->i1;
    int j = pphot->i2;
    int k = pphot->i3;

    // Add contribution to corresponding moments
    // Energy density
    if (general_mover_flag) 
      weight *= pphot->k[IMC0];
    moments(MCIER,k,j,i) += weight;
    // Flux
    moments(MCIFR1,k,j,i) += weight1;
    moments(MCIFR2,k,j,i) += weight2;
    moments(MCIFR3,k,j,i) += weight3;
    // Radiation Pressure
    Real weightp = weight1 * pphot->k[0];
    moments(MCIPR11,k,j,i) += weightp;
    weightp = weight2 * pphot->k[1];
    moments(MCIPR22,k,j,i) += weightp;
    weightp = weight3 * pphot->k[2];
    moments(MCIPR33,k,j,i) += weightp;
    weightp = weight1 * pphot->k[1];
    moments(MCIPR12,k,j,i) += weightp;
    weightp = weight1 * pphot->k[2];
    moments(MCIPR13,k,j,i)  += weightp;
    weightp = weight2 * pphot->k[2];
    moments(MCIPR23,k,j,i) += weightp;
    // Photon mean energy
    moments(MCIEN,k,j,i) += weight * pphot->energy;
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::NormalizeMoments(bool normalize)
//  \brief (un)normalized moments for output and copy symmetric elements

void MonteCarloBlock::NormalizeMoments(bool normalize) {

  if (normalize) {
    // Normalize moments
    for (int n=0; n<11; ++n) {
      //Real norm = static_cast<Real>(nphdone)*pmy_mc->normalization;
      Real norm = static_cast<Real>(nphdone);
      if ((n == 0) || (n >= 4))
	norm *= 2.9979e10;
      for (int k=ks; k<=ke; ++k) {
	for (int j=js; j<=je; ++j) {
	  for (int i=is; i<=ie; ++i) {
	    moments(n,k,j,i) /= (pcoord->vol(k,j,i) * norm);
	  }}}
    }
    // Copy noramilzed moments to symmetric elements
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
      //Real norm = static_cast<Real>(nphdone)*pmy_mc->normalization;
      Real norm = static_cast<Real>(nphdone);
      if ((n == 0) || (n >= 4))
	norm *= 2.9979e10;
      for (int k=ks; k<=ke; ++k) {
	for (int j=js; j<=je; ++j) {
	  for (int i=is; i<=ie; ++i) {
	    moments(n,k,j,i) *= (pcoord->vol(k,j,i) * norm);
	  }}}
    }
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

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::EnrollUserWorkInMove(UserMoveFunc_t userfunc)
//  \brief Enroll a user-defined condition to be called during photon moves

void MonteCarloBlock::EnrollUserWorkInMove(UserMoveFunc_t userfunc) {

  pmover->UserWorkInMove = userfunc;

}
