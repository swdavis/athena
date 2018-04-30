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

// constructor, initializes data structures and parameters

MonteCarloBlock::MonteCarloBlock(MeshBlock *pmb, MonteCarlo *pmc, ParameterInput *pin) {
  
  pmy_mc = pmc;

  // Set related meshblock, coordinate
  pmy_block = pmb;
  pmy_coord = pmb->pcoord;

  // Construct pointer to photon 
  pphoton  = new Photon(this);

  // Set photon mover based on coordinate system
  if (COORDINATE_SYSTEM == "cartesian") {
    pmover = new CartesianMover(this);
  }

  // Initialize input parameters and flags
  ntot = pin->GetInteger("montecarlo","nphot");
  zone_weight_flag = pin->GetOrAddBoolean("montecarlo","zone_weight",true);
  moments_flag = pin->GetOrAddBoolean("montecarlo","moments",true);
  lorentz_trans_flag = pin->GetOrAddBoolean("montecarlo","lorentz_trans",true);

  // Create and intitialize randon number generator
  int iseed = pin->GetInteger("montecarlo","iseed");
  pran = new MCRandom(iseed);

  prev=NULL;
  next=NULL;

  // read emission, absorption and scattering flags
  emission_meth = pmy_mc->emission_meth;
  //emission_meth = GetEmissionFlag(pin->GetOrAddString("montecarlo","emission","error"));
  absorption_meth = GetAbsorptionFlag(pin->GetOrAddString("montecarlo","absorption","error"));
  scattering_meth = GetScatteringFlag(pin->GetOrAddString("montecarlo","scattering","error"));

  // set local mesh parameters to correspond to mesh block
  is = pmb->is; ie = pmb->ie;
  js = pmb->js; je = pmb->je;
  ks = pmb->ks; ke = pmb->ke;

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
    //InitEmission = InitializeEmissionFreeFree;
    emission_array_flag = true;
  }
  if (emission_array_flag) emission.NewAthenaArray(ncells3,ncells2,ncells1);

  // Allocate variable arrays needed for evolution/output
  rho.NewAthenaArray(ncells3,ncells2,ncells1);
  tgas.NewAthenaArray(ncells3,ncells2,ncells1);
  if (lorentz_trans_flag) vel.NewAthenaArray(3,ncells3,ncells2,ncells1);
  if (moments_flag) moments.NewAthenaArray(4,ncells3,ncells2,ncells1);


}

// destructor

MonteCarloBlock::~MonteCarloBlock() {

  delete pphoton;
  delete pmover;
  delete pran;

  rho.DeleteAthenaArray();
  tgas.DeleteAthenaArray();
  if (lorentz_trans_flag) vel.DeleteAthenaArray();
  if (emission_array_flag) emission.DeleteAthenaArray();
  if (moments_flag) moments.DeleteAthenaArray();
}


enum AbsorptionFlag MonteCarloBlock::GetAbsorptionFlag(std::string input_string) {
  if (input_string == "user") {
    return ABSUSER;
  } else if (input_string == "none") {
    return ABSNONE;
  } else if (input_string == "freefree") {
    return ABSFF;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetAbsorptionFlag" << std::endl
        << "Input string=" << input_string << " not valid absorption type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}

enum ScatteringFlag MonteCarloBlock::GetScatteringFlag(std::string input_string) {
  if (input_string == "user") {
    return SCATUSER;
  } else if (input_string == "none") {
    return SCATNONE;
  } else if (input_string == "isotropic") {
    return SCATISO;
  } else if (input_string == "thomson") {
    return SCATTHOM;
  } else if (input_string == "compton") {
    return SCATCOMP;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetAbsorptionFlag" << std::endl
        << "Input string=" << input_string << " not valid scattering type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

}


void MonteCarloBlock::TransferPhotons() {

  for(int i=0; i<ntot; ++i) {

    // user definied photon initialization
    InitializePhoton(pmy_block,pphoton);

    // Lorentz transform E, k to Eulerian frame and update opacities.
    //if (lorenz_transform)
    //  photon->lorentz_transform(pmy_bplo,TOEUL);

    // Move photon to the first scattering/absorption event
    if (pphoton->status != DESTROYED) {
    //MovePhoton(pG,&Packet,pOut);
    }
    
    // Propogate photon packet until it leaves calculation domain or
    // it is absorbed.          
    /*    while(!(Packet.escape) && !(Packet.absorb)) {
#ifdef PHOTON_WEIGHT // Default method                    
      Packet.weight *= (Packet.sigma/(Packet.alpha+Packet.sigma));
#else
      if (random_mcgrid() > (Packet.sigma/(Packet.alpha+Packet.sigma)))
        Packet.weight = 0.0;
#endif
      pOut->nscat += 1;
      if(Packet.weight>WMIN) {
        // Scatter photon and reduce weight
#ifdef VELOCITIES
        // Lorentz transform to comoving frame for scattering
        lorentz_transform(&Packet,pG,to_comv);
#endif
        // Scatter the photon packet
        scatter(pG,&Packet);

        // Update the absorption and scattering extinction coefficients
        // with the new energy.
        ind=indexi(Packet.ix1,Packet.ix2,Packet.ix3,pG);
        Packet.alpha = absopac(pG->temp[ind],pG->dens[ind],Packet.energy);
        Packet.sigma = sctopac(pG->temp[ind],pG->dens[ind],Packet.energy);
        
#ifdef VELOCITIES
        // Lorentz transform to Eulerian frame and shift opacities
        lorentz_transform(&Packet,pG,to_eulr);
#endif

        // Move photon to next scattering/absorption event
        transfer(pG,&Packet,pOut);
        //if (errflag) return;
        
      } else {
        pOut->nabs++;
        Packet.absorb=1;
      }
      }*/

  }
  std::cout  << ntot << std::endl;
}

