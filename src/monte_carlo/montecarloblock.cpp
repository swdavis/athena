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

#define MINWEIGHT 1.0e-15

// constructor, initializes data structures and parameters

MonteCarloBlock::MonteCarloBlock(MeshBlock *pmb, MonteCarlo *pmc, ParameterInput *pin) {
  
  pmy_mc = pmc;

  // Set related meshblock, coordinate
  pmy_block = pmb;
  // Set pointer to this monte carlo block in pmb (***temporary***)
  pmb->pmy_mcb = this;
  
  pmy_coord = pmb->pcoord;

  // Construct pointer to photon 
  pphoton  = new Photon(this); // Currently one photon per block (will change)

  // Set photon mover based on coordinate system
  if (COORDINATE_SYSTEM == "cartesian") {
    pmover = new CartesianMover(this);
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    pmover = new SphericalPolarMover(this);
  }

  // Initialize input parameters and flags
  nphot = pin->GetInteger("montecarlo","nphot");
  zone_weight_flag = pin->GetOrAddBoolean("montecarlo","zone_weight",true);
  weighted_absorption = pin->GetOrAddBoolean("montecarlo","abs_weight",true);
  polarized = pin->GetOrAddBoolean("montecarlo","polarized",false);

  // get seed and intitialize randon number generator
  int rank = Globals::my_rank;
  int iseed = pmy_mc->iseed + rank *100;  // temporary solution
  pran = new MCRandom(iseed);

  prev=NULL;
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
  moments_flag = pmy_mc->moments_flag;
  
  // *currently** assumes all block boundaries are physical
  for (int i=0; i<6; ++i) {
    mcb_bcs[i] = pmy_mc->mc_bcs[i];
  }
  // Initialize pbval after mcb_bcs is set
  pbval = new MCBoundaryValues(this,pin);

  // Setup output spectrum **currently* assumes single spectrum
  int nfreq = pin->GetInteger("montecarlo","nfreq");
  int nmu = pin->GetInteger("montecarlo","nmu");
  int nphi = pin->GetInteger("montecarlo","nphi");
  pspec = new Spectrum(emin,emax,nfreq,nmu,nphi,polarized);

  // set local mesh parameters to correspond to mesh block
  is = pmb->is; ie = pmb->ie;
  js = pmb->js; je = pmb->je;
  ks = pmb->ks; ke = pmb->ke;
  ncells = (ie-is+1)*(je-js+1)*(ke-ks+1);
  codetocgs_rho = 1.0; codetoc_vel = 1.0;  // default cgs for code units

  // Allocate memory for emissivity (if used) and radiation moments
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  // allocate emission array for user function only if requested in input block
  // all current non-user methods require emission array
  if (emission_meth == EMISUSER) {
    emission_array_flag = pin->GetBoolean("montecarlo","emiss_array");
  } else if (emission_meth ==  EMISFF) {
    emission_array_flag = true;
  }
  if (emission_array_flag) emission.NewAthenaArray(ncells3,ncells2,ncells1);

  // Allocate (/initialize) variable arrays needed for evolution/output
  rho.NewAthenaArray(ncells3,ncells2,ncells1);
  tgas.NewAthenaArray(ncells3,ncells2,ncells1);
  if (pmy_mc->lorentz_trans_flag) vel.NewAthenaArray(3,ncells3,ncells2,ncells1);
  if (moments_flag) moments.NewAthenaArray(13,ncells3,ncells2,ncells1);

  // Set function pointers
  if (COORDINATE_SYSTEM == "cartesian") {
      GetZonePosition = GetZonePositionCartesian;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    GetZonePosition = GetZonePositionSphericalPolar;
  }
  if (absorption_meth == ABSUSER) {
    AbsorptionOpacity = NULL;
  } else if (absorption_meth == ABSNONE) {
    AbsorptionOpacity = NoOpacity;
  } else if (absorption_meth == ABSFF) {
    AbsorptionOpacity = FreeFreeAbsorptionOpacity;
  }
 
  if (scattering_meth == SCATUSER) {
    ScatteringOpacity = NULL;
    Scatter = NULL;
  } else if (scattering_meth == SCATNONE) {
    ScatteringOpacity = NoOpacity;
    Scatter = NULL;  // should not be called
    coherent_scattering = true;
  } else if (scattering_meth == SCATISO) {
    if (polarized) {
      std::stringstream msg;
      msg << "### ERROR in MonteCarloBlock constructor" << std::endl
          << "Istropic scattering but polarized = " << polarized << std::endl;
    throw std::runtime_error(msg.str().c_str());
    } else {
      ScatteringOpacity = ThomsonOpacity;
      Scatter = ScatterIsotropic;
      coherent_scattering = true;
    }
  } else if (scattering_meth == SCATTHOM) {
    ScatteringOpacity = ThomsonOpacity;
    if (polarized) {
      Scatter = ScatterThomsonPolarized;
    }
    coherent_scattering = true;
  } else if (scattering_meth == SCATCOMP) {
    GenerateComptonTable();
    ScatteringOpacity = ComptonOpacity;
    if (!polarized) {
      Scatter = ScatterComptonUnpolarized;
    } 
  }

}

// destructor

MonteCarloBlock::~MonteCarloBlock() {

  delete pphoton;
  delete pmover;
  delete pran;
  delete pspec;

  rho.DeleteAthenaArray();
  tgas.DeleteAthenaArray();
  if (lorentz_trans_flag) vel.DeleteAthenaArray();
  if (emission_array_flag) emission.DeleteAthenaArray();
  if (moments_flag) moments.DeleteAthenaArray();
}

void MonteCarloBlock::DefaultGetTemperature() {

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

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::TransferPhotons()
//  \brief perform radiation transfer of all photons on block

void MonteCarloBlock::TransferPhotons() {

  int nscat = 0, nesc = 0, nabs =0;
  for(int i=0; i<nphot; ++i) {

    // user definied photon initialization
    InitializePhoton(pphoton);
    
    //std::cout << pphoton->weight << std::endl;
    // Lorentz transform E, k to Eulerian frame and update opacities.
    //if (lorenz_transform)
    //  photon->lorentz_transform(pmy_bplo,TOEUL);

    // move photon to next scattering/absorption or to boundary
    pmover->Move(pphoton);
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
        
      // Lorentz transform to comoving frame for scattering
      //lorentz_transform(&Packet,pG,to_comv);

      // Scatter the photon packet
      if (pphoton->status == EVOLVING) {
        Scatter(this,pphoton);
	//pmover->CartesianToCurvalinear(pphoton);
	nscat++;
      }
      // Update the absorption and scattering extinction coefficients
      // with the new energy.
      if (!coherent_scattering) {
        pphoton->abs_coef = AbsorptionOpacity(this,pphoton);
        pphoton->sct_coef = ScatteringOpacity(this,pphoton);
      }
        
        // Lorentz transform to Eulerian frame and shift opacities
        //lorentz_transform(&Packet,pG,to_eulr);

      // move photon to next scattering/absorption or to boundary
      pmover->Move(pphoton);
    }
    if (pphoton->status == ESCAPED) {
      pspec->UpdateSpectrum(pphoton);
      nesc++;
    } else
      nabs++;
    
  }
  // Normalize moments for output
  if (moments_flag)
    NormalizeMoments(true);
  
  std::cout  << "nesc, nabs: " << nesc << ' ' << nabs << std::endl;
  std::cout << "nscat: " << nscat << std::endl;
}



//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateMoments(MonteCarloBlock *pmcb, Photon *pphot, Real dl)
//  \brief add contribution to radiation moments in current zone

void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl) {

  int& i = pphot->i1;
  int& j = pphot->i2;
  int& k = pphot->i3;
  
  Real wght = pphot->eweight * pphot->weight * pphot->energy * dl;
  Real wght1 = wght * pphot->k[0]; // curvalinear coorindates k
  Real wght2 = wght * pphot->k[1];
  Real wght3 = wght * pphot->k[2];

  moments(MCIER,k,j,i) += wght;
  moments(MCIFR1,k,j,i) += wght1;
  moments(MCIFR2,k,j,i) += wght2;
  moments(MCIFR3,k,j,i) += wght3;

  wght = wght1 * pphot->k[0];
  moments(MCIPR11,k,j,i) += wght;
  wght = wght2 * pphot->k[1];
  moments(MCIPR22,k,j,i) += wght;
  wght = wght3 * pphot->k[2];
  moments(MCIPR33,k,j,i) += wght;
  wght = wght1 * pphot->k[1];
  moments(MCIPR12,k,j,i) += wght;
  wght = wght1 * pphot->k[2];
  moments(MCIPR13,k,j,i)  += wght;
  wght = wght2 * pphot->k[2];
  moments(MCIPR23,k,j,i) += wght;
  
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::NormalizeMoments(bool normalize)
//  \brief (un)normalized moments for output and copy symmetric elements

void MonteCarloBlock::NormalizeMoments(bool normalize) {

  Coordinates *pco = pmy_coord;

  if (normalize) {
    // Normalize moments
    for (int n=0; n<10; ++n) {
      Real norm = static_cast<Real>(nphot)/static_cast<Real>(ncells);
      if ((n == 0) || (n >= 4))
	norm *= 2.9979e10;
      for (int k=ks; k<=ke; ++k) {
	for (int j=js; j<=je; ++j) {
	  for (int i=is; i<=ie; ++i) {
	    Real vol = pco->GetCellVolume(k,j,i);
	    moments(n,k,j,i) /= (vol * norm);
	  }}}}
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
    for (int n=0; n<10; ++n) {
      Real norm = static_cast<Real>(nphot)/static_cast<Real>(ncells);
      if ((n == 0) || (n >= 4))
	norm *= 2.9979e10;
      for (int k=ks; k<=ke; ++k) {
	for (int j=js; j<=je; ++j) {
	  for (int i=is; i<=ie; ++i) {
	    Real vol = pco->GetCellVolume(k,j,i);
	    moments(n,k,j,i) *= (vol * norm);
	  }}}}
  }
}
