//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file montecarloblock.cpp
//! \brief implementation of functions in class MonteCarloBlock

// C++ headers
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

// SWD: remove these
static Real test;
static bool first = true;

//----------------------------------------------------------------------------------------
//! MonteCarloBlock constructor, builds MonteCarloBlock from parameter input

MonteCarloBlock::MonteCarloBlock(MeshBlock *pmb,  MCBlockSize *pblsize, MonteCarlo *pmc,
                                 ParameterInput *pin) {

  pmy_mc = pmc;

  // Set related meshblock, coordinate
  pmy_block = pmb;

  // Set pointer to this monte carlo block in pmb if not nullptr
  if (pmb != nullptr) {
    pmb->pmy_mcb = this;
  }

  // Construct pointer to photon
  pphot  = new Photon(this,pmy_mc->nuser_var,pmy_mc->max_phots_init);

  // Initialize to nullptr and set below
  pmover = nullptr;
  pcoord = nullptr;

  // get seed and intitialize randon number generator
  int rank = Globals::my_rank;
  int iseed = pmy_mc->iseed + rank *100;  // temporary solution
  pran = new MCRandom(iseed);

  next=nullptr;

  // SWD: eliminate some or all of these?
  // set local flags based on monte_carlo
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
  Spectrum *pfirst = nullptr, *plast;
  Spectrum *psmcout = pmy_mc->pmcout->pspec;
  // Loop over output spectra and make local equivalent for each
  while (psmcout != nullptr) {
    pspec = new Spectrum(psmcout);
    //pspec = new Spectrum(pmy_mc->pmcout->pspec);
    if (pfirst == nullptr)
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
  if (pmb != nullptr) {
    is = pmb->is; ie = pmb->ie;
    js = pmb->js; je = pmb->je;
    ks = pmb->ks; ke = pmb->ke;
    nx1 = pmb->block_size.nx1;
    nx2 = pmb->block_size.nx2;
    nx3 = pmb->block_size.nx3;
  } else {
    if (pblsize == nullptr) {
      std::stringstream msg;
      msg << "### FATAL ERROR Monte Carlo Block Constructor" << std::endl
          << "Both input Mesh Block and Block Size are nullptr." << std::endl;
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
  // default cgs conversion from code units
  codetocgs_rho = 1.;
  codetocgs_vel = 1.;
  codetocgs_tgas = 1.;

  // SWD:  stepsize control needs to be modified
  stepsize = pin->GetOrAddReal("montecarlo","stepsize",1.0e-3);
  minweight = pin->GetOrAddReal("montecarlo","minweight",1.0e-20);

  // Flags for handling photon movement
  general_mover_flag = pin->GetOrAddBoolean("montecarlo","general_mover",false);
  sphpol_alt_flag = pin->GetOrAddBoolean("montecarlo", "sphpol_alt",false);
  boyerlindquist_flag = pin->GetOrAddBoolean("montecarlo","boyerlindquist",false);
  orthotet_flag = pin->GetOrAddBoolean("montecarlo", "orthotet", false);
  varystep_flag = pin->GetOrAddBoolean("montecarlo", "varystep", false);

  // Scattering
  scattering_meth = GetScatteringFlag(pin->GetOrAddString("montecarlo","scattering",
                                                          "none"));
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
  } else if (scattering_meth == SCATDUST) {
    ScatteringOpacity = DustScatteringOpacity;
    if (pmy_mc->polarized) {
      Scatter = ScatterDust;
      coherent_scattering = true;
    } else {
      std::stringstream msg;
      msg << "### ERROR in MonteCarloBlock constructor" << std::endl
          << "Dust scattering not suppored for polarized = "
          << pmy_mc->polarized << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }

  // Set up photon movement and initialization methods
  computedmin = false;
  if ((acceleration)||(sphpol_alt_flag))
    computedmin = true;
  pmy_mc->computedmin = computedmin;
  if (COORDINATE_SYSTEM == "cartesian") {
    GetZonePosition = GetZonePositionCartesian;
    if (general_mover_flag) {
      pmover = new GeneralMover(this);
      if (pmb != nullptr)
        pcoord = new MCCartesian(pmb->pcoord,this);
      else
        pcoord = new MCCartesian(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                 computedmin);
    } else {
      pmover = new CartesianMover(this);
      if (pmb != nullptr)
        pcoord = new MCCoord(pmb->pcoord,this);
      else
        pcoord = new MCCoord(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                             computedmin);
    }
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    GetZonePosition = GetZonePositionSphericalPolar;
    if (general_mover_flag) {
      pmover = new GeneralMover(this);
      if (pmb != nullptr)
        pcoord = new MCSphericalPolar(pmb->pcoord,this);
      else
        pcoord = new MCSphericalPolar(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                      computedmin);
    } else if (sphpol_alt_flag) {
      pmover = new SphericalPolarAltMover(this);
      if (pmb != nullptr)
        pcoord = new MCCoord(pmb->pcoord,this);
      else
        pcoord = new MCCoord(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                             computedmin);
    } else {
      pmover = new SphericalPolarMover(this);
      if (pmb != nullptr)
        pcoord = new MCCoord(pmb->pcoord,this);
      else
        pcoord = new MCCoord(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                             computedmin);
    }
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    GetZonePosition = GetZonePositionCylindrical;
    pmover = new GeneralMover(this);
    if (pmb != nullptr)
      pcoord = new MCCylindrical(pmb->pcoord,this);
    else
      pcoord = new MCCylindrical(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                 computedmin);
  } else if (COORDINATE_SYSTEM == "kerr-schild") {
    GetZonePosition = GetZonePositionSphericalPolar;//approximate
    pmover = new GeneralMover(this);
    if (boyerlindquist_flag) {
     if (pmb != nullptr)
       pcoord = new MCBoyerLindquist(pmb->pcoord,this);
     else {
       pcoord = new MCBoyerLindquist(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                     computedmin);
       pcoord->SetSpin(pin->GetReal("coord", "a"));
       pcoord->SetMass(pin->GetReal("coord", "m"));
     }
    } else {
      if (pmb != nullptr)
        pcoord = new MCKerrSchild(pmb->pcoord,this);
      else {
        pcoord = new MCKerrSchild(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                  computedmin);
        pcoord->SetSpin(pin->GetReal("coord", "a"));
        pcoord->SetMass(pin->GetReal("coord", "m"));
      }
    }
  } else if (COORDINATE_SYSTEM == "minkowski") {
    GetZonePosition = GetZonePositionCartesian;
    pmover = new GeneralMover(this);
    if (pmb != nullptr)
      pcoord = new MCMinkowski(pmb->pcoord,this);
    else
      pcoord = new MCMinkowski(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                               computedmin);
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

  // Set absorption opacity and method
  absorption_meth = GetAbsorptionMethodFlag(pin->GetOrAddString("montecarlo","abs_method",
                                                                "weight"));
  absorption_opac = GetAbsorptionOpacityFlag(pin->GetOrAddString("montecarlo",
                                                                 "absorption","none"));
  if (absorption_opac == ABSUSER) {
    AbsorptionOpacity = pmy_mc->UserAbsorptionOpacity;
  } else if (absorption_opac == ABSNONE) {
    AbsorptionOpacity = NoOpacity;
  } else if (absorption_opac == ABSFF) {
    AbsorptionOpacity = FreeFreeAbsorptionOpacity;
  } else if (absorption_opac == ABSDUST) {
    AbsorptionOpacity = DustAbsorptionOpacity;
  }

  // Allocate (/initialize) variable arrays needed for evolution/output
  int ncells1 = nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (nx2 > 1) ncells2 = nx2 + 2*(NGHOST);
  if (nx3 > 1) ncells3 = nx3 + 2*(NGHOST);
  rho.NewAthenaArray(ncells3,ncells2,ncells1);
  tgas.NewAthenaArray(ncells3,ncells2,ncells1);
  if (boosts) vel.NewAthenaArray(3,ncells3,ncells2,ncells1);
  if (NSCALARS > 0) scalars.NewAthenaArray(NSCALARS,ncells3,ncells2,ncells1);
  // moments is 1 (Er) + 3 (Fr) + 9 (Pr) + 1 (Eave) + 1 (net cool)
  if (moments_flag) moments.NewAthenaArray(NMOM,ncells3,ncells2,ncells1);
  if (emission_array_flag) emission.NewAthenaArray(ncells3,ncells2,ncells1);
  if (acceleration && !(coherent_scattering) && !(scattering_meth == SCATRES)) {
    planck_opacity.NewAthenaArray(ncells3,ncells2,ncells1);
    planck_inv_opacity.NewAthenaArray(ncells3,ncells2,ncells1);
  }

  // Create user monte carlo block data
  InitUserMonteCarloBlockData(pin);

}

//----------------------------------------------------------------------------------------
//! destructor

MonteCarloBlock::~MonteCarloBlock() {

  delete pphot;
  delete pmover;
  delete pran;
  delete pspec;
  delete pphlist;
  delete ptraj;

  rho.DeleteAthenaArray();
  tgas.DeleteAthenaArray();
  if (boosts) vel.DeleteAthenaArray();
  if (NSCALARS > 0) scalars.DeleteAthenaArray();
  if (moments_flag) moments.DeleteAthenaArray();
  if (emission_array_flag) emission.DeleteAthenaArray();
  if (acceleration && !(coherent_scattering) && !(scattering_meth == SCATRES)) {
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
//! \brief Integrate photons to termination condtion without scattering

void MonteCarloBlock::RayTracePhotons(int nphot) {

    Real const to_comv = 1.0;
    Real const to_eulr = -1.0;
    int nscat = 0, nesc = 0, nabs = 0, ndes = 0;
    int ntodo = (nphot > nphremain) ? nphremain : nphot;

    int nloop = 100;
    int nprop = ntodo;

    while (nprop > 0) {

      // Emit photons to replace those that left meshblock or were terminated
      nloop = (nloop > nprop) ? nprop : nloop;
      int nold = pphot->nphot;
      pphot->Resize(nloop);
      // user definied photon initialization
      InitializePhoton(pphot,nold,pphot->nphot-1);
      if (ptraj != nullptr) {
        for (int ip=nold; ip < pphot->nphot; ip++)
          ptraj->InitializeTrajectory(pphot->trp[ip]);
      }
      // Photon initialized in coordinate frame
      // move photon until  stopping condition
      pmover->Move(pphot,0,pphot->nphot-1);
      if (ptraj != nullptr) {
        for (int ip=nold; ip < pphot->nphot; ip++)
          ptraj->CompleteTrajectory(pphot->trp[ip]);
      }

      for (int ip=0; ip<pphot->nphot; ip++) {
        // User defined completion work
        FinalizePhoton(pphot,ip);
        if (pphot->statp[ip] == ESCAPED) {
          // loop over spectra and update
          Spectrum *pspect = pspec;
          while (pspect != nullptr) {
            pspect->UpdateSpectrum(pphot,ip);
            pspect = pspect->next;
          }
          if (pphlist != nullptr) {
            pphlist->AddPhoton(pphot,ip);
          }
          nesc++;
        } else if (pphot->statp[ip] == ABSORBED) {
          nabs++;
        } else if (pphot->statp[ip] == DESTROYED) {
          ndes++;
        }
        pphot->RemoveOneParticle(ip);
        nprop--;
      } // end loop over ip
    } // while nprop > 0

    std::cout  << "rank, nesc, nabs, ndes, nscat: " << Globals::my_rank << ' ' << nesc
               << ' ' << nabs << ' ' << ndes << ' '
               << static_cast<Real>(nscat)/static_cast<Real>(nprop)
               << std::endl;

    return;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::TransferPhotons()
//! \brief perform radiation transfer nphtot photons

void MonteCarloBlock::TransferPhotons(int nphot) {

  Real const to_comv = 1.0;
  Real const to_eulr = -1.0;
  int nscat = 0, nesc = 0, nabs = 0, ndes = 0;
  int ntodo = (nphot > nphremain) ? nphremain : nphot;
  nphdone += ntodo;

  int nloop = 100;
  int nprop = ntodo;
  while(nprop > 0) {

    // Emit photons to replace those that left meshblock or were terminated
    nloop = (nloop > nprop) ? nprop : nloop;
    int nold = pphot->nphot;
    pphot->Resize(nloop);

    //printf("nold: %d %d %d %d\n",nprop,nold,nprop,nloop);
    // user definied photon initialization
    InitializePhoton(pphot,nold,pphot->nphot-1);
    if (ptraj != nullptr) {
      for (int ip=nold; ip < pphot->nphot; ip++)
        ptraj->InitializeTrajectory(pphot->trp[ip]);
    }
    // Lorentz transform E, k to Eulerian frame and update opacities
    // only for newly emitted photons
    if (boosts) {
      LorentzTransform(pphot,to_eulr,nold,pphot->nphot-1);
    }
    if (moments_flag) {
      // Update cooling to relect newly emitted photons
      for (int ip=nold; ip<pphot->nphot; ip++) {
        UpdateSourceTerms(pphot,0.,0.,ip,0.,0.,0.);
      }
    }

    // move all photons to next interaction or boundary
    pmover->Move(pphot,0,pphot->nphot-1);

    for (int ip=0; ip<pphot->nphot; ip++) {

      if (pphot->statp[ip] == EVOLVING) {
        // Account for absorption
        Real weight0 = pphot->wp[ip];
        Real energy0 = pphot->ep[ip];
        Real k1p0 = pphot->k1p[ip];
        Real k2p0 = pphot->k2p[ip];
        Real k3p0 = pphot->k3p[ip];
        if (absorption_meth == ABSWEIGHT) {
          pphot->wp[ip] *= (pphot->scp[ip]/(pphot->scp[ip]+pphot->acp[ip]));
          if(pphot->wp[ip] <= minweight) {
            pphot->statp[ip] = ABSORBED;
          }
        } else if (absorption_meth == ABSPROB) {
          if (pran->uniform() > (pphot->scp[ip]/(pphot->scp[ip]+pphot->acp[ip])) )
            pphot->wp[ip] = 0.;
          pphot->statp[ip] = ABSORBED;
        } else if (absorption_meth == ABSTAU) {
          if(pphot->wp[ip] <= minweight) {
            pphot->statp[ip] = ABSORBED;
          }
        }
        if (moments_flag) {
          UpdateSourceTerms(pphot,energy0,weight0,ip,k1p0,k2p0,k3p0);
        }
      } // status == evolving

      if (pphot->statp[ip] == EVOLVING) {
        // Scatter the photon
        Real e_pre_scat = pphot->ep[ip];
        Real weight_pre_scat = pphot->wp[ip];
        Real k1p_pre_scat = pphot->k1p[ip];
        Real k2p_pre_scat = pphot->k2p[ip];
        Real k3p_pre_scat = pphot->k3p[ip];
        // Lorentz transform to comoving frame for scattering
        if (boosts) {
          LorentzTransform(pphot,to_comv,ip,ip);
        }
        Scatter(this,pphot,ip,ip);
        nscat++;
        pphot->nscp[ip]++;
        if (pphot->nscp[ip] %  pmy_mc->checkscat == 0) {
          // Check for possible infinite loop due to NaN in photon
          if (pphot->IsNanPhoton(ip)) {
            pphot->statp[ip] = DESTROYED;
            std::cout << "Warning: IsNanPhoton() returned true, photon destroyed"
                      << std::endl;
            pphot->PrintPhoton(ip);
          }
        }

        // Update the absorption and scattering extinction coefficients
        // with the new energy.
        if (!coherent_scattering) {
          pphot->acp[ip] = AbsorptionOpacity(this,pphot,ip);
          pphot->scp[ip] = ScatteringOpacity(this,pphot,ip);
        }
        // Lorentz transform to Eulerian frame and shift opacities
        if (boosts) {
          LorentzTransform(pphot,to_eulr,ip,ip);
        }
        if (moments_flag) {
            UpdateSourceTerms(pphot,e_pre_scat,weight_pre_scat,ip,
                              k1p_pre_scat,k2p_pre_scat,k3p_pre_scat);
        }
      } // status == evolving

    } // End loop over ip

    //for (int ip=0; ip<pphot->nphot; ip++) {
    for (int ip=pphot->nphot-1; ip >= 0; ip--) {
      if (pphot->statp[ip] != EVOLVING) {

        if (pphot->statp[ip] == ESCAPED) {
          if (ptraj != nullptr) {
            ptraj->CompleteTrajectory(pphot->trp[ip]);
          }
          // User defined completion work
          FinalizePhoton(pphot,ip);
          // SWD: temporary, needed for output
          //pphot->VectorsToWorkingArrays(ip);

          // loop over spectra and update
          Spectrum *pspect = pspec;
          while (pspect != nullptr) {
            pspect->UpdateSpectrum(pphot,ip);
            pspect = pspect->next;
          }
          if (pphlist != nullptr) {
            pphlist->AddPhoton(pphot,ip);
          }
          nesc++;
        } else if (pphot->statp[ip] == ABSORBED) {
          nabs++;
        } else if (pphot->statp[ip] == DESTROYED) {
          ndes++;
        }
        pphot->RemoveOneParticle(ip);
        nprop--;
      }
    } // End loop over ip
  } // while nprop > 0

  std::cout  << "rank, nesc, nabs, ndes, nscat: " << Globals::my_rank << ' ' << nesc
             << ' ' << nabs << ' ' << ndes << ' '
             << static_cast<Real>(nscat)/static_cast<Real>(ntodo) << std::endl;
}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::TransferPhotons()
//!  \brief perform radiation transfer nphtot photons

void MonteCarloBlock::TransferPhotonsOld(int nphot) {

  Real const to_comv = 1.0;
  Real const to_eulr = -1.0;
  int nscat = 0, nesc = 0, nabs = 0, ndes = 0;
  int nprop = (nphot > nphremain) ? nphremain : nphot;

  int nremain = nprop;
  pphot->Resize(10);
  while(nremain > 0) {
    //while (pphot->nphot < pphot->nphot_limit) {
      // user definied photon initialization
    InitializePhoton(pphot,0,pphot->nphot);
      //}

    // Lorentz transform E, k to Eulerian frame and update opacities.
    if (boosts) {
      LorentzTransform(pphot,to_eulr,0,pphot->nphot);
    }

    if (moments_flag) {
      for (int ip=0; ip<pphot->nphot; ip++) {
        UpdateSourceTerms(pphot,0.,0.,ip,0.,0.,0.);
      }
    }
    // move photon to next scattering/absorption or to boundary
    pmover->Move(pphot,0,pphot->nphot);

    for (int ip=0; ip<pphot->nphot; ip++) {

    int iscat = 0;
    while (pphot->statp[ip] == EVOLVING) {

      // Account for absorption
      Real weight0 = pphot->wp[ip];
      if (absorption_meth == ABSWEIGHT) {
        pphot->wp[ip] *= (pphot->scp[ip]/(pphot->scp[ip]+pphot->acp[ip]));
        if(pphot->wp[ip] <= minweight) {
          pphot->statp[ip] = ABSORBED;
        }
      } else if (absorption_meth == ABSPROB) {
        if (pran->uniform() > (pphot->scp[ip]/(pphot->scp[ip]+pphot->acp[ip])) )
          pphot->wp[ip] = 0.;
          pphot->statp[ip] = ABSORBED;
      } else if (absorption_meth == ABSTAU) {
        if(pphot->wp[ip] <= minweight) {
          pphot->statp[ip] = ABSORBED;
        }
      }
      if (moments_flag) {
          UpdateSourceTerms(pphot,0.,weight0,ip,0.,0.,0.);
      }

      // Scatter the photon packet
      if (pphot->statp[ip] == EVOLVING) {
        Real e_pre_scat = pphot->ep[ip];
        // Lorentz transform to comoving frame for scattering
        if (boosts) {
          LorentzTransform(pphot,to_comv,ip,ip);
        }
        Scatter(this,pphot,ip,ip);
        iscat++;

        if (iscat %  pmy_mc->checkscat == 0) {
          // Check for possible infinite loop due to NaN in photon
          if (pphot->IsNanPhoton()) {
            pphot->statp[ip] = DESTROYED;
            std::cout << "Warning: IsNanPhoton() returned true, photon destroyed"
                      << std::endl;
            pphot->PrintPhoton(ip);
          }
        }
        // Update the absorption and scattering extinction coefficients
        // with the new energy.
        if (!coherent_scattering) {
          pphot->acp[ip] = AbsorptionOpacity(this,pphot,ip);
          pphot->scp[ip] = ScatteringOpacity(this,pphot,ip);
        }
        // Lorentz transform to Eulerian frame and shift opacities
        if (boosts) {
          LorentzTransform(pphot,to_eulr,ip,ip);
        }
        if (moments_flag) {
            UpdateSourceTerms(pphot,e_pre_scat,0.,ip,0.,0.,0.);
        }
      }

      // move photon to next scattering/absorption or to boundary
      pmover->Move(pphot,ip,ip);

    }
    if (ptraj != nullptr) {
      ptraj->CompleteTrajectory(pphot->trp[ip]);
    }
    nscat += iscat;

    } // End loop over ip

    for (int ip=0; ip<pphot->nphot; ip++) {

    if (pphot->statp[ip] == ESCAPED) {
      pphot->VectorsToWorkingArrays(ip);
      // User defined completion work
      FinalizePhoton(pphot,ip);
      // loop over spectra and update
      Spectrum *pspect = pspec;
      while (pspect != nullptr) {
        pspect->UpdateSpectrum(pphot,ip);
        pspect = pspect->next;
      }
      if (pphlist != nullptr) {
        pphlist->AddPhoton(pphot,ip);
      }
      nesc++;
    } else if (pphot->statp[ip] == ABSORBED) {
      nabs++;
    } else if (pphot->statp[ip] == DESTROYED) {
      ndes++;
    }
    if (pphot->statp[ip] != EVOLVING)
      nremain--;

    } // End loop over ip
  }

  nphdone += nprop;
  std::cout  << "rank, nesc, nabs, ndes, nscat: " << Globals::my_rank << ' ' << nesc
             << ' ' << nabs << ' ' << ndes << ' '
             << static_cast<Real>(nscat)/static_cast<Real>(nprop) << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::ComovingToCoordinate(Photon *pphot, int ips, int ipe)
//! \brief Transform photon sample to coordinate/Eulerian frame

void MonteCarloBlock::ComovingToCoordinate(Photon *pphot, int ips, int ipe) {

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::ComovingToCoordinate(Photon *pphot, int ips, int ipe)
//! \brief Transform photon sample to comoving frame

void MonteCarloBlock::CoordinateToComoving(Photon *pphot, int ips, int ipe) {

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::LorentzTransform(Photon *pphot, const Real sign, int ips,
//!                                            int ipe)
//! \brief Lorentz transform photon sample
//
// Does not transform stokes vectors but this seems
// to be correct -- the plane of polarization is invariant under lorentz
// transformation as discussed in Cocke & Holm (1972) Nature letter.
// weight is not transformed either as weight represents number of photons
// in the packet which is invariant.
// to_comv: sign = 1.0;
// to_eulr: sign = -1.0;

void MonteCarloBlock::LorentzTransform(Photon *pphot, const Real sign, int ips,
                                       int ipe) {

  for (int ip=ips; ip<=ipe; ip++) {

    Real &k1 = pphot->k1p[ip];
    Real &k2 = pphot->k2p[ip];
    Real &k3 = pphot->k3p[ip];

    int i1 = pphot->i1p[ip];
    int i2 = pphot->i2p[ip];
    int i3 = pphot->i3p[ip];

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
      Real aber = gamma*(1.-gamma*bdk/(gamma+1.));
      pphot->ep[ip] *= gonembdk;
      k1 = (k1 - aber * beta[0]) / gonembdk;
      k2 = (k2 - aber * beta[1]) / gonembdk;
      k3 = (k3 - aber * beta[2]) / gonembdk;

      // Transform opacities
      // Must be performed even when transforming to comoving frame because inverse
      // process is performed to go back to Eulerian frame in cases where scattering
      // is coherent
      pphot->acp[ip] /= gonembdk;
      pphot->scp[ip] /= gonembdk;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real MonteCarloBlock::LorentzTransformFrequencyShift(Photon *pphot, int ip)
//!  \brief Returns frequency shift from Lorentz Trransformation

Real MonteCarloBlock::LorentzTransformFrequencyShift(Photon *pphot, int ip) {

  Real k1 = pphot->k1p[ip];
  Real k2 = pphot->k2p[ip];
  Real k3 = pphot->k3p[ip];
  int i1 = pphot->i1p[ip], i2 = pphot->i2p[ip], i3 = pphot->i3p[ip];

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

// SWD: This is an untested modification of Eric's original method, unfinished
//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::TetradTransform(Photon *pphot, const Real sign, int ips,
//!                                           int ipe)
//!  \brief Tetrad transform photon packet
//
// Does not transform stokes vectors..
// weight is not transformed either as weight represents number of photons
// in the packet which is invariant.
// to_comv: sign = 1.0;
// to_eulr: sign = -1.0;

void MonteCarloBlock::TetradTransform(Photon *pphot, const Real sign, int ips, int ipe) {

  for (int ip=ips; ip<=ipe; ip++) {
    // Get velocity of cell
    int i1 = pphot->i1p[ip], i2 = pphot->i2p[ip], i3 = pphot->i3p[ip];
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
    Real x[NCOORD];
    x[IMC0] = pphot->x0p[ip];
    x[IMC1] = pphot->x1p[ip];
    x[IMC2] = pphot->x2p[ip];
    x[IMC3] = pphot->x3p[ip];
    pcoord->Metric(x, gcov);

    // create tetrad basis
    Real econ[NCOORD][NCOORD], ecov[NCOORD][NCOORD];
    ConstructTetrad(ucon, gcov, econ, ecov);

    if (sign > 0) { // tranforming to comoving frame

      Real kcopy[NCOORD];
      kcopy[IMC0] = pphot->k0p[ip];
      kcopy[IMC1] = pphot->k1p[ip];
      kcopy[IMC2] = pphot->k2p[ip];
      kcopy[IMC3] = pphot->k3p[ip];
      Real k[NCOORD];
      CoordinateToTetrad(kcopy, k, ecov);
      pphot->k0p[ip] = k[IMC0];
      pphot->k1p[ip] = k[IMC1];
      pphot->k2p[ip] = k[IMC2];
      pphot->k3p[ip] = k[IMC3];

      Real energy_shift = kcopy[IMC0] / k[IMC0];
      // transform energy and extinction coefficients
      pphot->ep[ip] *= energy_shift;
      pphot->acp[ip] *= energy_shift;
      pphot->scp[ip] *= energy_shift;

    } else { // transforming to coordinate frame

      Real kcopy[NCOORD];
      kcopy[IMC0] = pphot->k0p[ip];
      kcopy[IMC1] = pphot->k1p[ip];
      kcopy[IMC2] = pphot->k2p[ip];
      kcopy[IMC3] = pphot->k3p[ip];
      Real k[NCOORD];
      TetradToCoordinate(kcopy, k, econ); // updates pphot->k

      Real energy_shift = kcopy[IMC0] / pphot->k[IMC0]; // new calculation

      // transform energy and opacities
      pphot->ep[ip] *= energy_shift;
      pphot->acp[ip] *= energy_shift;
      pphot->scp[ip] *= energy_shift;
    }
  } // loop over ip
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl, Real etau, int ip)
//! \brief add contribution to radiation moments in current zone

void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl, Real etau, int ip) {
  // SWD: needs to be modifed for non general mover kvectors
  
  Real k1 = pphot->k1p[ip];
  Real k2 = pphot->k2p[ip];
  Real k3 = pphot->k3p[ip];

  // Normalize k vector if using general mover in spherical polar coords
  if ((COORDINATE_SYSTEM == "spherical_polar") && (general_mover_flag)) {
    k2 *= pphot->x1p[ip];
    k3 *= pphot->x1p[ip] * sin(pphot->x2p[ip]);
  }

  Real energy, abs_coef, step;
  if (moments_comoving) {
    // boost relevant quanitities to comoving frame
    energy = pphot->ep[ip];
    int i1 = pphot->i1p[ip], i2 = pphot->i2p[ip], i3 = pphot->i3p[ip];
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
      abs_coef = pphot->acp[ip] / gonembdk;
      step = dl * gonembdk;
    }
  } else {
    // Use eulerian values
    energy = pphot->ep[ip];
    abs_coef = pphot->acp[ip];
    step = dl;
  }
  // Account for attenuation along ray
  Real leff;
  if (absorption_meth == ABSTAU) {
    if (fabs(1.-etau) < TINY_NUMBER) {
      leff = step;
    } else {
      leff = (1.-etau)/abs_coef;
      //printf("%g %g %g %g\n",etau,leff/step,abs_coef,step);
    }
  } else {
    leff = step;
  }
  // Weight moments by time spent in domain
  Real weight = pphot->wp[ip] * energy * leff / 2.99792458e10;
  if ((std::isinf(weight)) || (std::isnan(weight))) {
    std::cout << "Warning: UpdateMoments weight is : " << weight << std::endl;
  } else {
    // Higher order moments are weighted by curvalinear coordinates k
    Real weight1 = weight * k1;
    Real weight2 = weight * k2;
    Real weight3 = weight * k3;

    int i = pphot->i1p[ip];
    int j = pphot->i2p[ip];
    int k = pphot->i3p[ip];

    // SWD: Modify this appropriately
    //if (general_mover_flag)
    //  weight *= pphot->k0p[ip]

    // Add contribution to corresponding moments
    // Energy density
    moments(MCIER,k,j,i) += weight;
    // Flux
    moments(MCIFR1,k,j,i) += weight1 * 2.99792458e10;
    moments(MCIFR2,k,j,i) += weight2 * 2.99792458e10;
    moments(MCIFR3,k,j,i) += weight3 * 2.99792458e10;
    // Radiation Pressure

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
//! \brief (un)normalized moments for output and copy symmetric elements

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
    for (int n=0; n<14; ++n) {
      Real norm = normall;
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
//            if ((n==11) && (moments(n,k,j,i) != 0)) printf("moments: %g    vol: %g     norm: %g    i: %d   j: %d   k: %d    val: %g\n", moments(n,k,j,i), pcoord->vol(k,j,i), norm, i, j, k, moments(n,k,j,i)/pcoord->vol(k,j,i)/norm);
//            if ((n==11) && (moments(n,k,j,i) != 0)) printf("%g\n", moments(n,k,j,i)/pcoord->vol(k,j,i)/norm);
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
    for (int n=0; n<14; ++n) {
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
//! \brief set moments to zero on origin blocks

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
//! \fn void MonteCarloBlock::UpdateSourceTerms(Photon *pphot, Real energy0, Real weight0,
//                                          int ip)
//! \brief compute net photon cooling rate

void MonteCarloBlock::UpdateSourceTerms(Photon *pphot, Real energy0, Real weight0, int ip, Real k1p0, Real k2p0, Real k3p0) {
  Real c = 2.99792458e10;

  Real k1 = pphot->k1p[ip];
  Real k2 = pphot->k2p[ip];
  Real k3 = pphot->k3p[ip];

  // Normalize k vector if using general mover in spherical polar coords
  if ((COORDINATE_SYSTEM == "spherical_polar") && (general_mover_flag)) {
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

  // Components of momentum change --- assumes orthonormal basis
  Real dp1p = pphot->wp[ip] * k1 * pphot->ep[ip] / c - weight0 * k1p0 * energy0 / c;
  Real dp2p = pphot->wp[ip] * k2 * pphot->ep[ip] / c - weight0 * k2p0 * energy0 / c;
  Real dp3p = pphot->wp[ip] * k3 * pphot->ep[ip] / c - weight0 * k3p0 * energy0 / c;

  Real cool = (pphot->wp[ip] * pphot->ep[ip]) - (weight0 * energy0);
  //if (energy0 == 0.0)
  //  printf("weight, cool: %g %g\n",pphot->weight,cool);
  if ((std::isinf(cool)) || (std::isnan(cool))) {
    std::cout << "Warning: UpdateSourceTerms cooling is : " << cool << std::endl;
    pphot->PrintPhoton(ip);
  } else if ((std::isinf(dp1p)) || (std::isnan(dp1p))) {
    std::cout << "Warning: UpdateSourceTerms momentum change (k1p) is : " << dp1p << std::endl;
    pphot->PrintPhoton(ip);
  } else if ((std::isinf(dp2p)) || (std::isnan(dp2p))) {
    std::cout << "Warning: UpdateSourceTerms momentum change (k2p) is : " << dp2p << std::endl;
    pphot->PrintPhoton(ip);
  } else if ((std::isinf(dp3p)) || (std::isnan(dp3p))) {
    std::cout << "Warning: UpdateSourceTerms momentum change (k3p) is : " << dp3p << std::endl;
    pphot->PrintPhoton(ip);
  } else {
    int &i = pphot->i1p[ip];
    int &j = pphot->i2p[ip];
    int &k = pphot->i3p[ip];
    moments(MCINET,k,j,i) -= cool;
    moments(MCIP1,k,j,i) -= dp1p;
    moments(MCIP2,k,j,i) -= dp2p;
    moments(MCIP3,k,j,i) -= dp3p;
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs)
//! \brief set boundary values on monte carlo block

void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs) {

  // set x1 boundaries
  mcb_bcs[BoundaryFace::inner_x1] = input_bcs[BoundaryFace::inner_x1];
  mcb_bcs[BoundaryFace::outer_x1] = input_bcs[BoundaryFace::outer_x1];

  // set x2 boundaries
  mcb_bcs[BoundaryFace::inner_x2] = input_bcs[BoundaryFace::inner_x2];
  mcb_bcs[BoundaryFace::outer_x2] = input_bcs[BoundaryFace::outer_x2];

  // set x3 boundaries
  mcb_bcs[BoundaryFace::inner_x3] = input_bcs[BoundaryFace::inner_x3];
  mcb_bcs[BoundaryFace::outer_x3] = input_bcs[BoundaryFace::outer_x3];

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs)
//! \brief set boundary values on monte carlo block

/*void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs) {

  // set x1 boundaries
  if (pmy_block->pbval->block_bcs[BoundaryFace::inner_x1] == BLOCK_BNDRY) {
    mcb_bcs[BoundaryFace::inner_x1] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[BoundaryFace::inner_x1] = input_bcs[BoundaryFace::inner_x1];
  }
  if (pmy_block->pbval->block_bcs[BoundaryFace::outer_x1] == BLOCK_BNDRY) {
    mcb_bcs[BoundaryFace::outer_x1] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[BoundaryFace::outer_x1] = input_bcs[BoundaryFace::outer_x1];
  }

  // set x2 boundaries
  if (pmy_block->pbval->block_bcs[BoundaryFace::inner_x2] == BLOCK_BNDRY) {
    mcb_bcs[BoundaryFace::inner_x2] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[BoundaryFace::inner_x2] = input_bcs[BoundaryFace::inner_x2];
  }
  if (pmy_block->pbval->block_bcs[BoundaryFace::outer_x2] == BLOCK_BNDRY) {
    mcb_bcs[BoundaryFace::outer_x2] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[BoundaryFace::outer_x2] = input_bcs[BoundaryFace::outer_x2];
  }

  // set x3 boundaries
  if (pmy_block->pbval->block_bcs[BoundaryFace::inner_x3] == BLOCK_BNDRY) {
    mcb_bcs[BoundaryFace::inner_x3] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[BoundaryFace::inner_x3] = input_bcs[BoundaryFace::inner_x3];
  }
  if (pmy_block->pbval->block_bcs[BoundaryFace::outer_x3] == BLOCK_BNDRY) {
    mcb_bcs[BoundaryFace::outer_x3] = MC_BLOCK_BNDRY;
  } else {
    mcb_bcs[BoundaryFace::outer_x3] = input_bcs[BoundaryFace::outer_x3];
  }
  }*/
