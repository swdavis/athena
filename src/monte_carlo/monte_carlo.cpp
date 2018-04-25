//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file monte_carlo.cpp
//  \brief implementation of functions in class MonteCarlo

#include <iostream> // temporary for testing

#include <gsl/gsl_randist.h>

// Athena++ headers
#include "monte_carlo.hpp"
#include "photon.hpp"
#include "photonmover.hpp"
#include "../athena.hpp"
#include "../parameter_input.hpp"

// constructor, initializes data structures and parameters

MonteCarlo::MonteCarlo(MeshBlock *pmb, ParameterInput *pin) {
  
  pmy_block = pmb;

  // Construct pointer to photon 
  pphoton  = new Photon(this);

  // Set photon mover
  if (COORDINATE_SYSTEM == "cartesian") {
    pmover = new CartesianMover(this);
  }
  ntot = pin->GetInteger("montecarlo","nphot");
  zone_weight = pin->GetOrAddBoolean("montecarlo","zone_weight",true);

  int iseed = pin->GetInteger("montecarlo","iseed");
  pran = new MCRandom(iseed);
}

// destructor

MonteCarlo::~MonteCarlo() {

  delete pphoton;
}

void MonteCarlo::TransferPhotons() {

  for(int i=0; i<ntot; ++i) {

    // user definied photon initialization
    pphoton->InitializePhoton(pmy_block);

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

// constructor

MCRandom::MCRandom(int iseed) {
  dev = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(dev, iseed);
}

// destructor

MCRandom::~MCRandom() {

}

Real MCRandom::uniform() {

  return static_cast<Real>(gsl_rng_uniform(dev));
}
