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
#include "photonpusher.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../globals.hpp"
#include "../scalars/scalars.hpp"

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
  //pphot  = new Photon(this,pmy_mc->nuser_var,pmy_mc->max_phots_init);
  pphot  = new Photon(this,pin);

  // Initialize to nullptr and set below
  ppusher = nullptr;
  pcoord = nullptr;

  // get seed and intitialize randon number generator
  int rank = Globals::my_rank;
  int iseed = pmy_mc->iseed+pmy_block->gid*10;  // temporary solution
  //printf(" MonteCarloBlock gid %d rank %d iseed %d\n",pmy_block->gid,rank,iseed); 

  pran = new MCRandom(iseed);

  next=nullptr;

  // SWD: eliminate some or all of these?
  // set local flags based on monte_carlo
  // set in monte carlo
  boosts = pmy_mc->boosts;
  coupled = pmy_mc->coupled;
  acceleration = pmy_mc->acceleration;
  time_acc = pmy_mc->time_acc;
  // set in mcoutput if output requested
  mom_flag_lab = pmy_mc->pmcout->mom_flag_lab;
  mom_flag_com = pmy_mc->pmcout->mom_flag_com;
  if (mom_flag_com && !boosts) {
    std::stringstream msg;
    msg << "FATAL ERROR: comoving frame moments requested but booosts set to false."
        << std::endl;
    ATHENA_ERROR(msg);
  }
  mom_flag_src = pmy_mc->pmcout->mom_flag_src;
  mom_flag_usr = pmy_mc->pmcout->mom_flag_usr || (pmy_mc->nuser_mom > 0);
  mom_flag_scat = pmy_mc->pmcout->mom_flag_scat;

  call_srcterms = coupled || mom_flag_src;
  call_moments = mom_flag_lab || mom_flag_com || call_srcterms || mom_flag_usr;
  // Set boundary values for this block
  SetBoundaryValues(pmy_mc->mc_bcs);

  // Initialize pbval after mcb_bcs is set
  pbval = new MCBoundaryValues(this,pin);

  // Setup outputs
  pspec = pmy_mc->pmcout->pspec;
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
      std::cout << "Warning: comoving frame moments requested but booosts set to false.\n"
                << "S" << std::endl;
    }
  }
  mom_flag_src = pmy_mc->pmcout->mom_flag_src;
  mom_flag_usr = pmy_mc->pmcout->mom_flag_usr || (pmy_mc->nuser_mom > 0);

  call_srcterms = coupled || mom_flag_src;
  call_moments = mom_flag_lab || mom_flag_com || call_srcterms || mom_flag_usr;
  // Set boundary values for this block
  SetBoundaryValues(pmy_mc->mc_bcs);

  // Initialize pbval after mcb_bcs is set
  pbval = new MCBoundaryValues(this,pin);

  // Setup outputs
  pspec = pmy_mc->pmcout->pspec;
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
      ATHENA_ERROR(msg);
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
  rho_cgs = pin->GetOrAddReal("problem","rho_cgs",1.);
  vel_cgs = pin->GetOrAddReal("problem","vel_cgs",1.);
  tgas_cgs = pin->GetOrAddReal("problem","tgas_cgs",-1.);
  tfloor_cgs = pin->GetOrAddReal("problem","tfloor_cgs",0.);
  tceiling_cgs = pin->GetOrAddReal("problem","tceiling_cgs",HUGE_NUMBER);
  l_cgs = pin->GetOrAddReal("problem","l_cgs",1.);
  betamax = pin->GetOrAddReal("problem","betamax",0.999);

  // SWD:  stepsize control needs to be modified
  stepsize = pin->GetOrAddReal("montecarlo","stepsize",1.0e-3);

  // Flags for handling photon steps
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
      ATHENA_ERROR(msg);
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
    if (pmy_block->gid == 0) {
      if (comptonio > 0)
        std::cout << "Creating table for Compton cross section." << std::endl;
      else
        std::cout << "Reading in table for Compton cross section." << std::endl;
    }
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
      ATHENA_ERROR(msg);
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
      ATHENA_ERROR(msg);
    }
  }

  // Set up photon movement and initialization methods
  computedmin = false;
  if (acceleration)
    computedmin = true;
  pmy_mc->computedmin = computedmin;
  tetrads = true;
  if (COORDINATE_SYSTEM == "cartesian") {
    tetrads = false;
    GetZonePosition = GetZonePositionCartesian;
    if (pmy_mc->general_pusher_flag) {
      ppusher = new GeneralPusher(this);
      if (pmb != nullptr)
        pcoord = new MCCartesian(pmb->pcoord,this);
      else
        pcoord = new MCCartesian(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                 computedmin);
    } else {
      ppusher = new CartesianPusher(this);
      if (pmb != nullptr)
        pcoord = new MCCoord(pmb->pcoord,this);
      else
        pcoord = new MCCoord(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                             computedmin);
    }
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    GetZonePosition = GetZonePositionSphericalPolar;
    if (pmy_mc->general_pusher_flag) {
      ppusher = new GeneralPusher(this);
      if (pmb != nullptr)
        pcoord = new MCSphericalPolar(pmb->pcoord,this);
      else
        pcoord = new MCSphericalPolar(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                      computedmin);
    } else {
      tetrads = false;
      ppusher = new SphericalPolarPusher(this);
      if (pmb != nullptr)
        pcoord = new MCCoord(pmb->pcoord,this);
      else
        pcoord = new MCCoord(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                             computedmin);
    }
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    GetZonePosition = GetZonePositionCylindrical;
    ppusher = new GeneralPusher(this);
    if (pmb != nullptr)
      pcoord = new MCCylindrical(pmb->pcoord,this);
    else
      pcoord = new MCCylindrical(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                                 computedmin);
  } else if (COORDINATE_SYSTEM == "kerr-schild") {
    GetZonePosition = GetZonePositionSphericalPolar;//approximate
    ppusher = new GeneralPusher(this);
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
    ppusher = new GeneralPusher(this);
    if (pmb != nullptr)
      pcoord = new MCMinkowski(pmb->pcoord,this);
    else
      pcoord = new MCMinkowski(nx1+2*(NGHOST),nx2+2*(NGHOST),nx3+2*(NGHOST),
                               computedmin);
  } else if (COORDINATE_SYSTEM == "gr_user") {
    GetZonePosition = GetZonePositionCartesian;
    ppusher = new GeneralPusher(this);
    if (pmb != nullptr)
        pcoord = new MCKerrSchildCartesian(pmb->pcoord,this);
    pcoord->SetSpin(pin->GetReal("coord", "a"));
    pcoord->SetMass(pin->GetReal("coord", "m"));
  } else {
      std::stringstream msg;
      msg << "### ERROR in MonteCarloBlock constructor" << std::endl
          << COORDINATE_SYSTEM
          << " coordinates not currently supported with Monte Carlo"
          << std::endl;
      ATHENA_ERROR(msg);
  }
  pmy_mc->tetrads = tetrads;

  // Set pcoord in ppusher
  ppusher->pcoord = pcoord;

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
  nel.NewAthenaArray(ncells3,ncells2,ncells1);
  nion.NewAthenaArray(ncells3,ncells2,ncells1);
  tgas.NewAthenaArray(ncells3,ncells2,ncells1);
  if (boosts || tetrads) {
    boost_cmv.NewAthenaArray(ncells3,ncells2,ncells1,4,4);
    boost_lab.NewAthenaArray(ncells3,ncells2,ncells1,4,4);
  }
  if (boosts) vel.NewAthenaArray(ncells3,ncells2,ncells1,4);
  if (NSCALARS > 0) scalars.NewAthenaArray(ncells3,ncells2,ncells1);
  // moments is 1 (Er) + 3 (Fr) + 9 (Pr) + 1 (Eave) + 1 (net cool)
  nmom = 13;
  if (mom_flag_lab) moments.NewAthenaArray(nmom,ncells3,ncells2,ncells1);
  if (mom_flag_com) moments_com.NewAthenaArray(nmom,ncells3,ncells2,ncells1);
  if (pmy_mc->nuser_mom > 0)
    moments_user.NewAthenaArray(pmy_mc->nuser_mom,ncells3,ncells2,ncells1);
  nsrc = 8;
  if (call_srcterms) sourceterms.NewAthenaArray(nsrc,ncells3,ncells2,ncells1);
  if (mom_flag_scat) {
    nf_scat = pin->GetInteger("montecarlo","nf_scat");
    Real everg = 1.602176634e-12;
    emin_scat = pin->GetReal("montecarlo","emin_scat") * everg;
    emax_scat = pin->GetReal("montecarlo","emax_scat") * everg;
    moments_scat.NewAthenaArray(nf_scat,ncells3,ncells2,ncells1);
    dloge_scat = (std::log10(emax_scat/emin_scat))/static_cast<Real>(nf_scat);
    energy_scat.NewAthenaArray(nf_scat+1);
    freq_scat_mid.NewAthenaArray(nf_scat);
    Real h_cgs = 6.62607015e-27;
    for (int i=0; i<=nf_scat; i++) {
      energy_scat(i) = std::log10(emin_scat) + static_cast<Real>(i)*dloge_scat; // keep log
      if (i > 0)
        freq_scat_mid(i-1) = 0.5*( pow(10.,energy_scat(i-1)) + pow(10.,energy_scat(i)) )/h_cgs;
;
    }
    if ((Globals::my_rank ==0) && (pmy_block->lid == 0)) {
      FILE *pfile;
      std::stringstream msg;
      if ((pfile = fopen("sourceterm_frequencies.txt","w")) == NULL) {
        msg << "### FATAL ERROR in MonteCarloBlock Constructor" << std::endl
            << "Output file sourceterm_frequencies.txt could not be opened";
        throw std::runtime_error(msg.str().c_str());
      }
      fprintf(pfile,"%d\n", nf_scat);
      for (int i=0; i< nf_scat; i++) {
        fprintf(pfile,"%d %e %e %e\n",i,freq_scat_mid(i),pow(10.,energy_scat(i))/h_cgs,
                pow(10.,energy_scat(i+1))/h_cgs);
      }
      fclose(pfile);
    }
  }
  if (pmy_mc->emission_array) emission.NewAthenaArray(ncells3,ncells2,ncells1);
  if (pmy_mc->emission_eqwt[0]) emit_count_.NewAthenaArray(ncells3,ncells2,ncells1);
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
  delete ppusher;
  delete pcoord;
  delete pbval;
  delete pran;
  //delete pspec;
  //delete pphlist;
  //delete ptraj;

  rho.DeleteAthenaArray();
  nel.DeleteAthenaArray();
  nion.DeleteAthenaArray();
  tgas.DeleteAthenaArray();
  if (boosts || tetrads) {
    boost_cmv.DeleteAthenaArray();
    boost_lab.DeleteAthenaArray();
  }
  if (boosts) vel.DeleteAthenaArray();
  if (NSCALARS > 0) scalars.DeleteAthenaArray();
  if (mom_flag_lab) moments.DeleteAthenaArray();
  if (mom_flag_com) moments_com.DeleteAthenaArray();
  if (pmy_mc->nuser_mom > 0) moments_user.DeleteAthenaArray();
  if (call_srcterms) sourceterms.DeleteAthenaArray();
  if (pmy_mc->emission_array) emission.DeleteAthenaArray();
  if (pmy_mc->emission_eqwt[0]) emit_count_.DeleteAthenaArray();
  if (acceleration && !(coherent_scattering) && !(scattering_meth == SCATRES)) {
    planck_opacity.DeleteAthenaArray();
    planck_inv_opacity.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::RayTracePhotonsOnBlock()
//! \brief Integrate photons to termination condtion without scattering

void MonteCarloBlock::RayTracePhotonsOnBlock() {

  Real const to_comv = 1.0;
  Real const to_eulr = -1.0;
  int nbuf = 0;

  printf("remain: %d \n",nphremain);
  // Emit photons to replace those that left meshblock or were terminated
  // Limit ntodo to number of remaining photons on block
  int ntodo = (loop_max_size > nphremain) ? nphremain : loop_max_size;

  // if photons remain to transfer, make space for new photons
  if (ntodo > 0) {
    int nold = pphot->nphot;
    pphot->AllocatePhotons(nold+ntodo);
    nphremain -= ntodo;
    nphrun += ntodo;

    // user definied photon initialization
    InitializePhoton(pphot,nold,pphot->nphot-1);
    if (ptraj != nullptr) {
      for (int ip=nold; ip < pphot->nphot; ip++)
        ptraj->InitializeTrajectory(pphot->trp[ip]);
    }
  }
  int ntot = pphot->nphot;

  // Photon initialized in coordinate frame
  // move photon until  stopping condition
  ppusher->Move(pphot,0,pphot->nphot-1);

  for (int ip=pphot->nphot-1; ip >= 0; ip--) {
    if (pphot->statp[ip] != EVOLVING) {

      if (pphot->statp[ip] != BUFFERED) {
        // User defined completion work
        FinalizePhoton(pphot,ip);

        if (ptraj != nullptr) {
          ptraj->CompleteTrajectory(pphot->trp[ip]);
        }
      }
      if (pphot->statp[ip] == ESCAPED) {
        // loop over outputs for escaping photons and update
        Spectrum *pspect = pspec;
        while (pspect != nullptr) {
          pspect->UpdateSpectrum(pphot,ip);
          pspect = pspect->next;
        }
        if (pphlist != nullptr) {
          pphlist->AddPhoton(pphot,ip);
        }
        nesc++;
        pphot->RemoveOneParticle(ip);
      } else if (pphot->statp[ip] == ABSORBED) {
        nabs++;
        pphot->RemoveOneParticle(ip);
      } else if (pphot->statp[ip] == DESTROYED) {
        ndes++;
        pphot->PrintPhoton("destroy in ray tracing",ip);
        pphot->RemoveOneParticle(ip);
      } else if (pphot->statp[ip] == BUFFERED) {
        nbuf++;
      }
    }
  } // end loop over ip
  //if ((pphot->nphot > 0) && (Globals::my_rank == 0))
  //  pphot->PrintPhoton(pphot->nphot-1);
  /*std::cout  << "rank, ntot, nnew, nesc, nabs, ndes, nbuf, nscat: " << Globals::my_rank
             << ' ' << ntot << ' ' << ntodo << ' ' << nesc
             << ' ' << nabs << ' ' << ndes << ' ' << nbuf
             << ' ' << nscat << std::endl;*/
}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::TransferPhotonsOnBlock()
//! \brief perform radiation transfer

void MonteCarloBlock::TransferPhotonsOnBlock() {

  //int nbuf = 0;
  int nold = pphot->nphot;
  int ntot = nold + nphremain;

  // Emit photons to replace those that left meshblock or were terminated
  // limit ntot < loop_max_size unless nold is larger than loop_max_size
  ntot = (loop_max_size > ntot) ? ntot : loop_max_size;
  ntot = (nold > ntot) ? nold : ntot;
  int nnew = ntot - nold;

  if (ntot == 0) // nothing to do
    return;

  // if photons remain to transfer, make space for new photons
  if (nnew > 0) {
    pphot->AllocatePhotons(ntot);
    nphremain -= nnew;
    nphrun += nnew;

    // user definied photon initialization
    InitializePhoton(pphot,nold,pphot->nphot-1);

    // Lorentz transform E, k to Eulerian frame and update opacities
    // only for newly emitted samples
    if (boosts || tetrads) {
      TransformToCoordinate(pphot,nold,pphot->nphot-1);
    }
  
    // Update the absorption and scattering extinction coefficients
    if (call_srcterms) {
      // Update source terms to reflect newly emitted samples
      for (int ip=nold; ip<pphot->nphot; ip++) {
        UpdateSourceTerms(pphot,0.,0.,0.,0.,0.,ip);
      }
    }
  }

  // move all samples to next interaction or boundary
  ppusher->Move(pphot,0,pphot->nphot-1);

  // perform all absorption and scattering related tasks for all samples
  for (int ip=0; ip<pphot->nphot; ip++) {
    // record initial weight and direction
    Real weight0 = pphot->wp[ip];
    Real k1p0 = pphot->k1p[ip];
    Real k2p0 = pphot->k2p[ip];
    Real k3p0 = pphot->k3p[ip];

    // account for absorption
    if (pphot->statp[ip] == EVOLVING) {
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
    } // status == evolving

    // account for scattering
    if (pphot->statp[ip] == EVOLVING) {
      Real e_pre_scat = pphot->ep[ip];
      // Lorentz transform to comoving frame for scattering
      if (boosts || tetrads) {
        TransformToComoving(pphot,ip,ip);
      }
      // call scattering function and update counters
      Scatter(this,pphot,ip,ip);
      nscat++;
      pphot->nscp[ip]++;
      if (pphot->nscp[ip] %  pmy_mc->checkscat == 0) {
        // Check for possible infinite loop due to NaN in photon
        if (pphot->IsNanPhoton(ip)) {
          pphot->statp[ip] = DESTROYED;
          if (pmy_mc->verbose) {
            pphot->PrintPhoton("Warning: Nan encounterd in TransferPhotons(),"
                               " photon destroyed",ip);
            }
        }
      }

      // Update the absorption and scattering extinction coefficients
      // if scattering can change sample energy.
      if (!coherent_scattering) {
        pphot->acp[ip] = AbsorptionOpacity(this,pphot,ip);
        pphot->scp[ip] = ScatteringOpacity(this,pphot,ip);
      }
      // Lorentz transform to Eulerian frame and shift opacities
      if (boosts || tetrads) {
        TransformToCoordinate(pphot,ip,ip);
      }
      // Update moments that compute radiation force and net heating/cooling
      if (call_srcterms) {
        UpdateSourceTerms(pphot,e_pre_scat,weight0,k1p0,k2p0,k3p0,ip);
      }
    } // status == evolving

  } // End loop over ip

  // Perform tasks for samples that have left meshblock or otherwise terminated
  // their evolution.  Loop is reversed because of way particles (photon samples)
  // are popped
  for (int ip=pphot->nphot-1; ip >= 0; ip--) {
    if (pphot->statp[ip] != EVOLVING) {
      if (pphot->statp[ip] != BUFFERED) {
        // User defined completion work
        FinalizePhoton(pphot,ip);
      }

      if (pphot->statp[ip] == ESCAPED) {
       // loop over outputs for escaping photons and update
        Spectrum *pspect = pspec;
        while (pspect != nullptr) {
          pspect->UpdateSpectrum(pphot,ip);
          pspect = pspect->next;
        }
        if (pphlist != nullptr) {
          pphlist->AddPhoton(pphot,ip);
        }
        nesc++;
        pphot->RemoveOneParticle(ip);
      } else if (pphot->statp[ip] == ABSORBED) {
        nabs++;
        pphot->RemoveOneParticle(ip);
      } else if (pphot->statp[ip] == DESTROYED) {
        pphot->RemoveOneParticle(ip);
        ndes++;
      //} else if (pphot->statp[ip] == BUFFERED) {
      //  nbuf++;
      }
    }
  } // End loop over ip

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::CoupleMonteCarloToFluid(Real dt)
//! \brief update hydro momentum and energy based on radiative cooling and forces from MC

void MonteCarloBlock::CoupleMonteCarloToFluid(Real dt) {

  if (!coupled) return;

  Real edot_cgs = pmy_mc->time_cgs / (rho_cgs * SQR(vel_cgs));
  Real pdot_cgs = pmy_mc->time_cgs / (rho_cgs * vel_cgs);

  MeshBlock *pmb = pmy_block;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          pmb->phydro->u(IEN,k,j,i) += dt * edot_cgs * sourceterms(MCRS0,k,j,i);
          pmb->phydro->u(IM1,k,j,i) += dt * pdot_cgs * sourceterms(MCRF1,k,j,i)
                                          * pmb->phydro->u(IDN,k,j,i);
          pmb->phydro->u(IM2,k,j,i) += dt * pdot_cgs * sourceterms(MCRF2,k,j,i)
                                          * pmb->phydro->u(IDN,k,j,i);
          pmb->phydro->u(IM3,k,j,i) += dt * pdot_cgs * sourceterms(MCRF3,k,j,i)
                                          * pmb->phydro->u(IDN,k,j,i);
        }
      }
  }
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
      beta[i] = sign * vel(i3,i2,i1,i+1)/vel(i3,i2,i1,0);
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
      // is computed

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
    beta[i] = vel(i3,i2,i1,i+1)/vel(i3,i2,i1,0);
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
      beta[i] = sign * vel(i3,i2,i1,i+1);
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

      Real energy_shift = kcopy[IMC0] / pphot->k0p[ip]; // new calculation

      // transform energy and opacities
      pphot->ep[ip] *= energy_shift;
      pphot->acp[ip] *= energy_shift;
      pphot->scp[ip] *= energy_shift;
    }
  } // loop over ip
}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl, Real etau, int ip)
//! \brief Overload for UpdateMoments with additional wait time argument

void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl, Real etau, int ip) {

  // Account for attenuation along ray
  Real leff;
  if (absorption_meth == ABSTAU) {
    if (fabs(1.-etau) < TINY_NUMBER) {
      leff = dl;
    } else {
      leff = (1.-etau)/pphot->acp[ip];;
    }
  } else {
    leff = dl;
  }
  UpdateMoments(pphot, leff, ip);
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl, int ip)
//! \brief add contribution to radiation moments in current zone

void MonteCarloBlock::UpdateMoments(Photon *pphot, Real dl, int ip) {

  // SWD Needs to add support for moments in three bases:
  // comoving frame (supported already)
  // tetrad frame  (currently lab frame)
  // coordinate frame

  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  Real k0,k1,k2,k3,weight;
  const Real c_cgs = 2.99792458e10;
  if (pmy_mc->general_pusher_flag) {
    Real ki[4];
    ki[0] = pphot->k0p[ip];
    ki[1] = pphot->k1p[ip];
    ki[2] = pphot->k2p[ip];
    ki[3] = pphot->k3p[ip];
    Real x[4];
    x[0] = pphot->x0p[ip];
    x[1] = pphot->x1p[ip];
    x[2] = pphot->x2p[ip];
    x[3] = pphot->x3p[ip];
    Real invtet[4][4], kf[4];
    pcoord->InverseTetrad(x,invtet);
    for (int j=0; j<4; j++) {
      kf[j] = 0.;
      for (int i=0; i<4; i++) {
        kf[j] += invtet[j][i] * ki[i];
      }
    }
    k0 = kf[0];
    k1 = kf[1];
    k2 = kf[2];
    k3 = kf[3];
    // Weight moments by time spent in domain
    weight = pphot->wp[ip] * pphot->ep[ip] / k0 * dl * l_cgs / c_cgs;
  } else {
    k0 = pphot->k0p[ip];
    k1 = pphot->k1p[ip];
    k2 = pphot->k2p[ip];
    k3 = pphot->k3p[ip];
    // Weight moments by time spent in domain
    weight = pphot->wp[ip] * pphot->ep[ip] * dl * l_cgs / c_cgs;
  }

  // Normalize k vector if using general pusher in spherical polar coords
  //if ((COORDINATE_SYSTEM == "spherical_polar") && (pphot->general_pusher_flag)) {
  //  k2 *= pphot->x1p[ip];
  //  k3 *= pphot->x1p[ip] * sin(pphot->x2p[ip]);
  //}

  if (mom_flag_lab) {
    if (std::isinf(weight) || std::isnan(weight) || std::isnan(k0) || std::isnan(k1) || 
        std::isnan(k2) || std::isnan(k3)) {
      pphot->statp[ip] = DESTROYED;
      if (pmy_mc->verbose) {
        pphot->PrintPhoton ("Warning: Nan/Inf encountered in UpdateMoments(),"
                            " photon destroyed",ip);
      }
      return;
    } else {
      // Add contribution to corresponding moments
      // Energy density
      moments(MCIER,i3,i2,i1) += weight * k0 * k0;
      // Flux
      moments(MCIFR1,i3,i2,i1) += weight * k0 * k1 * c_cgs;
      moments(MCIFR2,i3,i2,i1) += weight * k0 * k2 * c_cgs;
      moments(MCIFR3,i3,i2,i1) += weight * k0 * k3 * c_cgs;
      // Radiation Pressure
      moments(MCIPR11,i3,i2,i1) += weight * k1 * k1;
      moments(MCIPR22,i3,i2,i1) += weight * k2 * k2;
      moments(MCIPR33,i3,i2,i1) += weight * k3 * k3;
      moments(MCIPR12,i3,i2,i1) += weight * k1 * k2;
      moments(MCIPR13,i3,i2,i1) += weight * k1 * k3;
      moments(MCIPR23,i3,i2,i1) += weight * k2 * k3;
    }
  }

  // add contribution to scattering source terms
  // SWD: Ultimately want comoving frame values
  if (mom_flag_scat) {
    Real loge = std::log10(pphot->ep[ip]);
    Real log10 = 2.302585092994046;
    int n = std::floor((loge-energy_scat(0))/dloge_scat);
    if (n >= 0 && n < nf_scat) {
      Real norm = c_cgs/(4.*PI*freq_scat_mid(n)*dloge_scat*log10);
      moments_scat(n,i3,i2,i1) += norm * pphot->scp[ip] * weight * k0 * k0;
    }
  }

  if (mom_flag_com) {
    // boost relevant quanitities to comoving frame
    //FrequencyAngleShiftComoving(pphot,ip,shift,k1c,k2c,k3c);

    Real ki[4],kc[4];
    ki[0] = k0;
    ki[1] = k1;
    ki[2] = k2;
    ki[3] = k3;
    for (int j=0; j<4; j++) {
      kc[j] = 0.;
      for (int i=0; i<4; i++) {
        kc[j] += boost_cmv(i3,i2,i1,j,i) * ki[i];
       }
    }
    Real k0c,k1c,k2c,k3c,weight;
    if (pmy_mc->general_pusher_flag) {
      k0c = kc[0];
      k1c = kc[1];
      k2c = kc[2];
      k3c = kc[3];
      weight = pphot->wp[ip] * pphot->ep[ip] / k0c * dl * l_cgs / c_cgs;
    } else {
      Real shift = kc[0]/k0;
      k1c = kc[1]/kc[0];
      k2c = kc[2]/kc[0];
      k3c = kc[3]/kc[0];
      weight = pphot->wp[ip] * pphot->ep[ip] * dl * SQR(shift) * l_cgs / c_cgs;
    }

    //Real dlcom = dl * shift;
    //Real ecom = pphot->ep[ip] * shift;

    //Real weight = pphot->wp[ip] * ecom * dlcom / c_cgs;
    if (std::isinf(weight) || std::isnan(weight) ||
        std::isinf(k1c) || std::isnan(k1c) ||
        std::isinf(k2c) || std::isnan(k2c) ||
        std::isinf(k3c) || std::isnan(k3c) ) {
      pphot->statp[ip] = DESTROYED;
      if (pmy_mc->verbose) {
        pphot->PrintPhoton("Warning: Nan/Inf encountered in UpdateMoments(),"
                           " comoving frame",ip);
      }
      return;
    } else {
      // Add contribution to corresponding moments
      // Energy density
      moments_com(MCIER,i3,i2,i1) += weight * k0c * k0c;
      // Flux
      moments_com(MCIFR1,i3,i2,i1) += weight * k0c * k1c * c_cgs;
      moments_com(MCIFR2,i3,i2,i1) += weight * k0c * k2c * c_cgs;
      moments_com(MCIFR3,i3,i2,i1) += weight * k0c * k3c * c_cgs;
      // Radiation Pressure
      moments_com(MCIPR11,i3,i2,i1) += weight * k1c * k1c;
      moments_com(MCIPR22,i3,i2,i1) += weight * k2c * k2c;
      moments_com(MCIPR33,i3,i2,i1) += weight * k3c * k3c;
      moments_com(MCIPR12,i3,i2,i1) += weight * k1c * k2c;
      moments_com(MCIPR13,i3,i2,i1) += weight * k1c * k3c;
      moments_com(MCIPR23,i3,i2,i1) += weight * k2c * k3c;
    }
  }

  if (mom_flag_usr) {
    for (int i=0; i<pmy_mc->nuser_mom; i++) {
      pmy_mc->user_moment_func[i](this,pphot,dl,ip,i);
    }
  }

  if (call_srcterms) {
    // Radiative Acceleration from flux (always lab frame)
    Real weight = pphot->wp[ip] * pphot->ep[ip] * dl * l_cgs / c_cgs;
    Real abs_coef = pphot->acp[ip];
    Real sct_coef = pphot->scp[ip];
    sourceterms(MCRF1,i3,i2,i1) += (sct_coef+abs_coef) * weight * k1;
    sourceterms(MCRF2,i3,i2,i1) += (sct_coef+abs_coef) * weight * k2;
    sourceterms(MCRF3,i3,i2,i1) += (sct_coef+abs_coef) * weight * k3;
  }


}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateMomentsAcceleration(Photon *pphot, Real dl, Real pl,
//        Real k1, Real k2, Real k3,Real etau, int ip)
//! \brief add contribution to radiation moments in current zone for acceleration

void MonteCarloBlock::UpdateMomentsAcceleration(Photon *pphot, Real dl, Real pl, Real k1,
                                                Real k2, Real k3, Real etau, int ip) {

  const Real c_cgs = 2.99792458e10;;
  Real k1p = pphot->k1p[ip];
  Real k2p = pphot->k2p[ip];
  Real k3p = pphot->k3p[ip];

  // Normalize k vector if using general pusher in spherical polar coords

  if ((COORDINATE_SYSTEM == "spherical_polar") && (pphot->general_pusher_flag)) {
    k2p *= pphot->x1p[ip];
    k3p *= pphot->x1p[ip] * sin(pphot->x2p[ip]);
  }

  Real energy, abs_coef, sct_coef, step;
  // BCM: Comoving moments currently do not work with code acceleration
  if (mom_flag_com) {
    // boost relevant quanitities to comoving frame
    energy = pphot->ep[ip];
    int i1 = pphot->i1p[ip], i2 = pphot->i2p[ip], i3 = pphot->i3p[ip];
    Real beta[3];
    for (int i=0; i<3; ++i) {
      beta[i] = vel(i3,i2,i1,i+1);
    }
    Real beta2= SQR(beta[0]) + SQR(beta[1]) + SQR(beta[2]);

    if(beta2 > 0.) {
      Real gamma = 1. / sqrt(1. - beta2); // assumes v^2 < c^2 checked elsewhere
      Real bdk = k1p * beta[0] + k2p * beta[1] + k3p * beta[2];
      Real gonembdk = gamma * (1. - bdk);
      Real aber = gamma*(1.-gamma*bdk/(gamma+1.));

      energy *= gonembdk;
      k1p = (k1p - aber * beta[0]) / gonembdk;
      k2p = (k2p - aber * beta[1]) / gonembdk;
      k3p = (k3p - aber * beta[2]) / gonembdk;
      abs_coef = pphot->acp[ip] / gonembdk;
      sct_coef = pphot->scp[ip] / gonembdk;
      step = dl * gonembdk;
    }
  } else {
    // Use eulerian values
    energy = pphot->ep[ip];
    abs_coef = pphot->acp[ip];
    sct_coef = pphot->scp[ip];
    step = dl;
  }
  // Account for attenuation along ray
  Real leff;
  if (absorption_meth == ABSTAU) {
    if (fabs(1.-etau) < TINY_NUMBER) {
      leff = step;
    } else {
      leff = (1.-etau)/abs_coef;
    }
  } else {
    leff = step;
  }
  // Weight moments by time spent in domain

  Real weight = pphot->wp[ip] * energy * leff / c_cgs;
  Real path_weight = weight * (pl / dl);

  if ((std::isinf(weight)) || (std::isnan(weight))) {
    pphot->statp[ip] = DESTROYED;
    if (pmy_mc->verbose) {
      pphot->PrintPhoton("Warning: Nan/Inf encountered in UpdateMoments(),"
                         " photon destroyed",ip);
    }
  } else {
    // Higher order moments are weighted by displacement direction vector k
    Real weight1 = weight * k1;
    Real weight2 = weight * k2;
    Real weight3 = weight * k3;

    int i = pphot->i1p[ip];
    int j = pphot->i2p[ip];
    int k = pphot->i3p[ip];

    if (mom_flag_lab) {
      // Add contribution to corresponding moments
      // Energy density
      moments(MCIER,k,j,i) += path_weight;
      // Flux
      moments(MCIFR1,k,j,i) += weight1 * c_cgs;
      moments(MCIFR2,k,j,i) += weight2 * c_cgs;
      moments(MCIFR3,k,j,i) += weight3 * c_cgs;

      // Radiation Pressure
      Real weightp = weight1 * k1p;
      moments(MCIPR11,k,j,i) += weightp;
      weightp = weight2 * k2p;
      moments(MCIPR22,k,j,i) += weightp;
      weightp = weight3 * k3p;
      moments(MCIPR33,k,j,i) += weightp;
      weightp = weight1 * k2p;
      moments(MCIPR12,k,j,i) += weightp;
      weightp = weight1 * k3p;
      moments(MCIPR13,k,j,i)  += weightp;
      weightp = weight2 * k3p;
      moments(MCIPR23,k,j,i) += weightp;
      // Photon mean energy
      //moments(MCIEN,k,j,i) += weight * energy;
      // Jmean opacity
      //moments(MCIKJ,k,j,i) += weight * abs_coef;
    }

    if (call_srcterms) {
      // Radiative Acceleration from flux
      sourceterms(MCRF1,k,j,i) += (sct_coef+abs_coef) * weight1;
      sourceterms(MCRF2,k,j,i) += (sct_coef+abs_coef) * weight2;
      sourceterms(MCRF3,k,j,i) += (sct_coef+abs_coef) * weight3;
    }
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateMomentsOld(Photon *pphot, Real dl, Real etau, int ip)
//! \brief add contribution to radiation moments in current zone
// SWD: remove this!
void MonteCarloBlock::UpdateMomentsOld(Photon *pphot, Real dl, Real pl, Real k1, Real k2,
                                    Real k3, Real etau, int ip) {

  const Real c_cgs = 2.99792458e10;
  Real k1p = pphot->k1p[ip];
  Real k2p = pphot->k2p[ip];
  Real k3p = pphot->k3p[ip];

  // Normalize k vector if using general pusher in spherical polar coords

  if ((COORDINATE_SYSTEM == "spherical_polar") && (pphot->general_pusher_flag)) {
    k2p *= pphot->x1p[ip];
    k3p *= pphot->x1p[ip] * sin(pphot->x2p[ip]);
  }

  Real energy, abs_coef, sct_coef, step;
  // BCM: Comoving moments currently do not work with code acceleration
  if (mom_flag_com) {
    // boost relevant quanitities to comoving frame
    energy = pphot->ep[ip];
    int i1 = pphot->i1p[ip], i2 = pphot->i2p[ip], i3 = pphot->i3p[ip];
    Real beta[3];
    for (int i=0; i<3; ++i) {
      beta[i] = vel(i3,i2,i1,i+1);
    }
    Real beta2= SQR(beta[0]) + SQR(beta[1]) + SQR(beta[2]);

    if(beta2 > 0.) {
      Real gamma = 1. / sqrt(1. - beta2); // assumes v^2 < c^2 checked elsewhere
      Real bdk = k1p * beta[0] + k2p * beta[1] + k3p * beta[2];
      Real gonembdk = gamma * (1. - bdk);
      Real aber = gamma*(1.-gamma*bdk/(gamma+1.));

      energy *= gonembdk;
      k1p = (k1p - aber * beta[0]) / gonembdk;
      k2p = (k2p - aber * beta[1]) / gonembdk;
      k3p = (k3p - aber * beta[2]) / gonembdk;
      abs_coef = pphot->acp[ip] / gonembdk;
      sct_coef = pphot->scp[ip] / gonembdk;
      step = dl * gonembdk;
    }
  } else {
    // Use eulerian values
    energy = pphot->ep[ip];
    abs_coef = pphot->acp[ip];
    sct_coef = pphot->scp[ip];
    step = dl;
  }
  // Account for attenuation along ray
  Real leff;
  if (absorption_meth == ABSTAU) {
    if (fabs(1.-etau) < TINY_NUMBER) {
      leff = step;
    } else {
      leff = (1.-etau)/abs_coef;
    }
  } else {
    leff = step;
  }
  // Weight moments by time spent in domain

  Real weight = pphot->wp[ip] * energy * leff / c_cgs;
  Real path_weight = weight * (pl / dl);

  if ((std::isinf(weight)) || (std::isnan(weight))) {
    pphot->statp[ip] = DESTROYED;
    if (pmy_mc->verbose) {
      pphot->PrintPhoton("Warning: Nan/Inf encountered in UpdateMoments(),"
                         " photon destroyed",ip);
    }
                        } else {
    // Higher order moments are weighted by displacement direction vector k
    Real weight1 = weight * k1;
    Real weight2 = weight * k2;
    Real weight3 = weight * k3;

    int i = pphot->i1p[ip];
    int j = pphot->i2p[ip];
    int k = pphot->i3p[ip];

    if (mom_flag_lab) {
      // Add contribution to corresponding moments
      // Energy density
      moments(MCIER,k,j,i) += path_weight;
      // Flux
      moments(MCIFR1,k,j,i) += weight1 * c_cgs;
      moments(MCIFR2,k,j,i) += weight2 * c_cgs;
      moments(MCIFR3,k,j,i) += weight3 * c_cgs;

      // Radiation Pressure
      Real weightp = weight1 * k1p;
      moments(MCIPR11,k,j,i) += weightp;
      weightp = weight2 * k2p;
      moments(MCIPR22,k,j,i) += weightp;
      weightp = weight3 * k3p;
      moments(MCIPR33,k,j,i) += weightp;
      weightp = weight1 * k2p;
      moments(MCIPR12,k,j,i) += weightp;
      weightp = weight1 * k3p;
      moments(MCIPR13,k,j,i)  += weightp;
      weightp = weight2 * k3p;
      moments(MCIPR23,k,j,i) += weightp;

      // Photon mean energy
      //moments(MCIEN,k,j,i) += weight * energy;
      // Jmean opacity
      //moments(MCIKJ,k,j,i) += weight * abs_coef;
    }

    if (call_srcterms) {
      // Radiative Acceleration from flux
      sourceterms(MCRF1,k,j,i) += (sct_coef+abs_coef) * weight1;
      sourceterms(MCRF2,k,j,i) += (sct_coef+abs_coef) * weight2;
      sourceterms(MCRF3,k,j,i) += (sct_coef+abs_coef) * weight3;
    }
  }

}

//------------------------------------------------- --------------------------------------
//! \fn void MonteCarloBlock::NormalizeMoments(bool normalize)
//! \brief (un)normalized moments for output and copy symmetric elements

void MonteCarloBlock::NormalizeMoments(bool normalize) {

  // Get integration time
  // Fix for dynamic MC
  //Real dt = pmy_mc->tint;
  Real dt = pmy_mc->pmy_mesh->time;
  Real norm;

  if (mom_flag_lab) {
    for (int n=0; n<nmom-3; ++n) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            if (normalize) {
              norm = 1./ (dt * pcoord->vol(k,j,i));
            } else {
              norm = dt * pcoord->vol(k,j,i);
            }
            moments(n,k,j,i) *= norm;
          }
        }
      }
    }
    // Copy normalized moments to symmetric elements
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          moments(MCIPR21,k,j,i) = moments(MCIPR12,k,j,i);
          moments(MCIPR31,k,j,i) = moments(MCIPR13,k,j,i);
          moments(MCIPR32,k,j,i) = moments(MCIPR23,k,j,i);
        }
      }
    }
  } // end if (mom_flag_lab)
  if (mom_flag_com) {
    for (int n=0; n<nmom-3; ++n) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            if (normalize)
              norm = 1./ (dt * pcoord->vol(k,j,i));
            else
              norm = dt * pcoord->vol(k,j,i);
            moments_com(n,k,j,i) *= norm;
          }
        }
      }
    }
    // Copy normalized moments to symmetric elements
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          moments_com(MCIPR21,k,j,i) = moments_com(MCIPR12,k,j,i);
          moments_com(MCIPR31,k,j,i) = moments_com(MCIPR13,k,j,i);
          moments_com(MCIPR32,k,j,i) = moments_com(MCIPR23,k,j,i);
        }
      }
    }
  }
  if (mom_flag_scat) {
    for (int n=0; n<nf_scat; ++n) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            if (normalize)
              norm = 1./ (dt * pcoord->vol(k,j,i));
            else
              norm = dt * pcoord->vol(k,j,i);
            moments_scat(n,k,j,i) *= norm;
          }
        }
      }
    }
  }
  if (mom_flag_usr) {
    for (int n=0; n<pmy_mc->nuser_mom; ++n) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            if (normalize)
              norm = 1./ (dt * pcoord->vol(k,j,i));
            else
              norm = dt * pcoord->vol(k,j,i);
            moments_user(n,k,j,i) *= norm;
          }
        }
      }
    }
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::ResetMoments()
//! \brief set moments to zero on block

void MonteCarloBlock::ResetMoments() {

    // set moments to zero
  for (int n=0; n<NMOM-5; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          moments(n,k,j,i) = 0.;
        }
      }
    }
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::UpdateSourceTerms(Photon *pphot, Real energy0,
//                                              Real weight0, int ip, Real k1p0,
//                                              Real k2p0, Real k3p0)
//! \brief compute net photon cooling rate and momentum change

void MonteCarloBlock::UpdateSourceTerms(Photon *pphot, Real energy0,
                                        Real weight0, Real k1p0,
                                        Real k2p0, Real k3p0, int ip) {

  // Updates
  Real c_cgs = 2.99792458e10;
  Real k1 = pphot->k1p[ip];
  Real k2 = pphot->k2p[ip];
  Real k3 = pphot->k3p[ip];

  // Normalize k vector if using general pusher in spherical polar coords
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

  // Components of momentum change --- assumes orthonormal basis
  Real dp1p = pphot->wp[ip] * k1 * pphot->ep[ip] / c_cgs
              - weight0 * k1p0 * energy0 / c_cgs;
  Real dp2p = pphot->wp[ip] * k2 * pphot->ep[ip] / c_cgs
              - weight0 * k2p0 * energy0 / c_cgs;
  Real dp3p = pphot->wp[ip] * k3 * pphot->ep[ip] / c_cgs
              - weight0 * k3p0 * energy0 / c_cgs;

  Real cool = (pphot->wp[ip] * pphot->ep[ip]) - (weight0 * energy0);

  if ((std::isinf(cool)) || (std::isnan(cool))) {
    std::cout << "Warning: UpdateSourceTerms cooling is : " << cool << std::endl;
    pphot->PrintPhoton(ip);
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
    sourceterms(MCRS0,k,j,i) -= cool;
    sourceterms(MCRS1,k,j,i) -= dp1p;
    sourceterms(MCRS2,k,j,i) -= dp2p;
    sourceterms(MCRS3,k,j,i) -= dp3p;
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::NormalizeSourceTerms(bool normalize)
//! \brief (un)normalized source terms for output

void MonteCarloBlock::NormalizeSourceTerms(bool normalize) {

  // Get integration time
  Real dt = pmy_mc->tint;

  // Normalize sourcterms
  for (int n=0; n<nsrc; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (normalize)
            sourceterms(n,k,j,i) /= (dt * pcoord->vol(k,j,i));
          else
            sourceterms(n,k,j,i) *= (dt * pcoord->vol(k,j,i));
        }
      }
    }
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::ResetSourceTerms()
//! \brief set sourceterms to zero on block

void MonteCarloBlock::ResetSourceTerms() {

  // set sourceterms to zero
  for (int n=0; n<8; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          sourceterms(n,k,j,i) = 0.;
        }
      }
    }
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs)
//! \brief set boundary values on monte carlo block

void MonteCarloBlock::SetBoundaryValues(enum MCBoundaryFlag *input_bcs) {

  // set x1 boundaries
  if(pmy_block->pbval->block_bcs[BoundaryFace::inner_x1] == BoundaryFlag::block)
    mcb_bcs[BoundaryFace::inner_x1] = MC_BLOCK_BNDRY;
  else
    mcb_bcs[BoundaryFace::inner_x1] = input_bcs[BoundaryFace::inner_x1];

  if(pmy_block->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::block)
    mcb_bcs[BoundaryFace::outer_x1] = MC_BLOCK_BNDRY;
  else
    mcb_bcs[BoundaryFace::outer_x1] = input_bcs[BoundaryFace::outer_x1];

  // set x2 boundaries
  if(pmy_block->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::block)
    mcb_bcs[BoundaryFace::inner_x2] = MC_BLOCK_BNDRY;
  else
    mcb_bcs[BoundaryFace::inner_x2] = input_bcs[BoundaryFace::inner_x2];

  if(pmy_block->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::block)
    mcb_bcs[BoundaryFace::outer_x2] = MC_BLOCK_BNDRY;
  else
    mcb_bcs[BoundaryFace::outer_x2] = input_bcs[BoundaryFace::outer_x2];

  // set x3 boundaries
  if(pmy_block->pbval->block_bcs[BoundaryFace::inner_x3] == BoundaryFlag::block)
    mcb_bcs[BoundaryFace::inner_x3] = MC_BLOCK_BNDRY;
  else
    mcb_bcs[BoundaryFace::inner_x3] = input_bcs[BoundaryFace::inner_x3];

  if(pmy_block->pbval->block_bcs[BoundaryFace::outer_x3] == BoundaryFlag::block)
    mcb_bcs[BoundaryFace::outer_x3] = MC_BLOCK_BNDRY;
  else
    mcb_bcs[BoundaryFace::outer_x3] = input_bcs[BoundaryFace::outer_x3];

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::ComputeEmissionArray(int etype, Real &em_min, Real &em_max, Real
//!                                                &em_tot)
//! \brief compute emission array

void MonteCarloBlock::ComputeEmissionArray(int etype, Real &em_min, Real &em_max, Real &em_tot) {


  Real dt = pmy_mc->tint;
  em_min = SQR(HUGE_NUMBER);
  em_max = -HUGE_NUMBER;
  em_tot = 0.;

  EmisFunc_t GetEmission = pmy_mc->GetEmission[etype];
  if (pmy_mc->emission_geometry[etype] == EMISVOL) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real vol = pcoord->vol(k,j,i);
          emission(k,j,i) = GetEmission(this,k,j,i,etype) * dt * vol;
          em_tot += emission(k,j,i);
          if (emission(k,j,i) > em_max) em_max = emission(k,j,i);
          if (emission(k,j,i) < em_min) em_min = emission(k,j,i);
          if (std::isnan(emission(k,j,i)))
            printf("emission[%d %d %d]: %g %g %g\n",i,j,k,vol,emission(k,j,i),pmy_mc->GetEmission[etype](this,k,j,i,etype));
        }
      }
    }
  } else if (pmy_mc->emission_geometry[etype] == EMISAREA) {
    // emmision array is reused so loop over block and reset to zero
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          emission(k,j,i) = 0.;
        }
      }
    }
    int il=is, iu=ie, jl=js, ju=je, kl=ks, ku=ke;
    int ip=0, jp=0, kp=0;
    int iface;
    BoundaryFace face = pmy_mc->emission_face[etype];
    // Check if we are on a physical boundary corresponding to emission face
    // Could use the mesh information, but this will generalize to arbitrary
    // surface in future
    Mesh *pm = pmy_block->pmy_mesh;
    Real tol = 0.001;
    bool physical_boundary = false;
    switch(face) {
      case BoundaryFace::inner_x1: {
        Real diff = std::fabs(pm->mesh_size.x1min-pcoord->x1f(is));
        Real dx = std::fabs(pcoord->x1f(is+1)-pcoord->x1f(is));
        if (diff < tol*dx) {
          physical_boundary = true;
          il = iu = is;
          iface = 0;
        } 
        break;
      }
      case BoundaryFace::outer_x1: {
        Real diff = std::fabs(pm->mesh_size.x1max-pcoord->x1f(ie+1));
        Real dx = std::fabs(pcoord->x1f(ie+1)-pcoord->x1f(ie));
        if (diff < tol*dx) {
          physical_boundary = true;
          il = iu = ie;
          iface = 3;
        }
        break;
      }
      case BoundaryFace::inner_x2: {
        Real diff = std::fabs(pm->mesh_size.x2min-pcoord->x2f(js));
        Real dx = std::fabs(pcoord->x2f(js+1)-pcoord->x2f(js));
        if (diff < tol*dx) {
          physical_boundary = true;
          jl = ju = js;
          iface = 1;
        }
        break;
      }
      case BoundaryFace::outer_x2: {
        Real diff = std::fabs(pm->mesh_size.x2max-pcoord->x2f(je+1));
        Real dx = std::fabs(pcoord->x2f(je+1)-pcoord->x2f(je));
        if (diff < tol*dx) {
          physical_boundary = true;
          jl = ju = je;
          iface = 4;
        }
        break;
      }
      case BoundaryFace::inner_x3: {
        Real diff = std::fabs(pm->mesh_size.x3min-pcoord->x3f(ks));
        Real dx = std::fabs(pcoord->x3f(ks+1)-pcoord->x3f(ks));
        if (diff < tol*dx) {
          physical_boundary = true;
          kl = ku = ks;
          iface = 2;
        }
        break;
      }
      case BoundaryFace::outer_x3: {
        Real diff = std::fabs(pm->mesh_size.x3max-pcoord->x3f(ke+1));
        Real dx = std::fabs(pcoord->x3f(ke+1)-pcoord->x3f(ke));
        if (diff < tol*dx) {
          physical_boundary = true;
          kl = iu = ke;
          iface = 5;
        }
        break;
      }
      default:
        std::stringstream msg;
        msg << "### FATAL ERROR in function [MonteCarloBlock::SetEmissionCellWeightArea]"
            << std::endl << "Face not valid" << std::endl;
        throw std::runtime_error(msg.str().c_str());
        break;
    }
    if (physical_boundary) {
      Coordinates *pbcoord = pmy_block->pcoord; // SWD: Maybe should improve this
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
          for (int i=il; i<=iu; ++i) {
            Real area;
            if (iface == 0)
              area = pbcoord->GetFace1Area(k,j,i);
            else if (iface == 1)
              area = pbcoord->GetFace2Area(k,j,i);
            else if (iface == 2)
              area = pbcoord->GetFace3Area(k,j,i);
            else if (iface == 3)
              area = pbcoord->GetFace1Area(k,j,i+1);
            else if (iface == 4)
              area = pbcoord->GetFace2Area(k,j+1,i);
            else if (iface == 5)
              area = pbcoord->GetFace3Area(k+1,j,i);
            area *= l_cgs*l_cgs;
            emission(k,j,i) = GetEmission(this,k,j,i,etype) * dt * area;
            em_tot += emission(k,j,i);
            if (emission(k,j,i) > em_max) em_max = emission(k,j,i);
            if (emission(k,j,i) < em_min) em_min = emission(k,j,i);
            if (std::isnan(emission(k,j,i)))
              printf("emission[%d %d %d]: %g %g %g\n",i,j,k,area,emission(k,j,i),
                     GetEmission(this,k,j,i,etype));
          }
        }
      }
    } else {
      // not a physical boundary so there is emission
      em_tot = em_max = em_min = 0.;
    }
  } // end if EMISAREA
  // if using equal weight scheme, intialize variables for SetEmissionCellWeight
  i1_ = -1; i2_= -1; i3_ = -1;

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::ComputeEmissionSampleArray()
//! \brief compute emission array for equal weight scheme

void MonteCarloBlock::ComputeEmissionSampleArray() {

  int ncells = nx1 * nx2 * nx3;
  Real prob[ncells];
  int count[ncells];

  // contruct probability array
  Real total_emission = 0.;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        int n = (k-ks)*nx2*nx1 + (j-js)*nx1 + i-is;
        prob[n] = emission(k,j,i);
        total_emission += emission(k,j,i);
      }
    }
  }

  for (int i=0; i<ncells; ++i) {
    prob[i] /= total_emission;
  }
  // sample multinomial distribution
  pran->SampleMultinomial(nphremain,ncells,prob,count);
  // set counts in emit_count_ array
  for (int k=ks; k<=ke; ++k) {
    int sum =0;
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        int n = (k-ks)*nx2*nx1 + (j-js)*nx1 + i-is;
        emit_count_(k,j,i) = count[n];
        sum += count[n];
      }
    }
  }

}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::SetEmissionCellWeight(Photon *pphot, int ips, int ipe)
//! \brief set emission cell and weight for photon via emission array

void MonteCarloBlock::SetEmissionCellWeight(Photon *pphot, int ips, int ipe) {

  if (pmy_mc->emission_eqwt[0]) {
    // Set intial zone based on probability within zone
    
    for (int ip=ips; ip<=ipe; ip++) {
      bool this_zone = false;
      while (!this_zone) {
        int i = i1_ + is;
        int j = i2_ + js;
        int k = i3_ + ks;
        if (emit_count_(k,j,i) > 0) {
          pphot->i1p[ip] = i;
          pphot->i2p[ip] = j;
          pphot->i3p[ip] = k;
          this_zone = true;
          emit_count_(k,j,i) -= 1;
        } else {
          // Update zone
          this_zone = false;
          i3_++;
          if (i3_ >= nx3) {
            i3_ = 0;
            i2_++;
            if (i2_ >= nx2) {
              i2_ = 0;
              i1_++;
              if (i1_ >= nx1)
                i1_ = 0;
            }
          }
        }
      } // end while (!this_zone)
      // Set weight to constant value for all photons
      pphot->wp[ip] = emiss_to_weight;
      
    } // end loop over ip
  } else {
    for (int ip=ips; ip<=ipe; ip++) {
      // Randomly assign emission zone
      pphot->i1p[ip] = static_cast<int>(pran->uniform()*nx1)+is;
      pphot->i2p[ip] = static_cast<int>(pran->uniform()*nx2)+js;
      pphot->i3p[ip] = static_cast<int>(pran->uniform()*nx3)+ks;

      // Set weight according to the emission array, which is the relative number of
      // photons emitted in each cell
      pphot->wp[ip] = emission(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip])
        * emiss_to_weight;
    } // end loop over ip
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::SetEmissionCellWeightArea(Photon *pphot, BoundaryFace face,
//        int ips, int ipe)
//! \brief set emission cell and weight for photon via emission array

void MonteCarloBlock::SetEmissionCellWeightArea(Photon *pphot, BoundaryFace face, int ips,
                                                int ipe) {

  if (pmy_mc->emission_eqwt[0]) {
    // Set intial zone based on probability within zone
    
    for (int ip=ips; ip<=ipe; ip++) {
      bool i1flag = true;
      bool i2flag = true;
      bool i3flag = true;
      switch(face) {
        case BoundaryFace::inner_x1:
          i1_ = 0;
          i1flag = false;
          break;
        case BoundaryFace::outer_x1:
          i1_ = nx1-1;
          i1flag = false;
          break;
        case BoundaryFace::inner_x2:
          i2_ = 0;
          i2flag = false;
          break;
        case BoundaryFace::outer_x2:
          i2_ = nx2-1;
          i2flag = false;
          break;
        case BoundaryFace::inner_x3:
          i3_ = 0;
          i3flag = false;
          break;
        case BoundaryFace::outer_x3:
          i3_ = nx3-1;
          i3flag = false;
          break;
        default:
          std::stringstream msg;
          msg << "### FATAL ERROR in function [MonteCarloBlock::SetEmissionCellWeightArea]"
              << std::endl << "Face not valid" << std::endl;
          throw std::runtime_error(msg.str().c_str());
          break;
      }
      bool this_zone = false;
      while (!this_zone) {
        int i = i1_ + is;
        int j = i2_ + js;
        int k = i3_ + ks;
        if (emit_count_(k,j,i) > 0) {
          pphot->i1p[ip] = i;
          pphot->i2p[ip] = j;
          pphot->i3p[ip] = k;
          this_zone = true;
          emit_count_(k,j,i) -= 1;
        } else {
          // Update zone
          this_zone = false;
          if (i3flag) i3_++;
          if (i3_ >= nx3) {
            i3_ = 0;
            if (i2flag) i2_++;
            if (i2_ >= nx2) {
              i2_ = 0;
              if (i1flag) i1_++;
              if (i1_ >= nx1)
                i1_ = 0;
            }
          }
        }
      } // end while (!this_zone)
      // Set weight to constant value for all photons
      pphot->wp[ip] = emiss_to_weight;
      
    } // end loop over ip
  } else {

    for (int ip=ips; ip<=ipe; ip++) {
      // Randomly assign emission zone
      Real weight_reduce;
      switch(face) {
        case BoundaryFace::inner_x1:
          pphot->i1p[ip] = is;
          pphot->i2p[ip] = static_cast<int>(pran->uniform()*nx2)+js;
          pphot->i3p[ip] = static_cast<int>(pran->uniform()*nx3)+ks;
          weight_reduce = 1./static_cast<Real>(nx1);
          break;
        case BoundaryFace::outer_x1:
          pphot->i1p[ip] = ie;
          pphot->i2p[ip] = static_cast<int>(pran->uniform()*nx2)+js;
          pphot->i3p[ip] = static_cast<int>(pran->uniform()*nx3)+ks;
          weight_reduce = 1./static_cast<Real>(nx1);
          break;
        case BoundaryFace::inner_x2:
          pphot->i1p[ip] = static_cast<int>(pran->uniform()*nx1)+is;
          pphot->i2p[ip] = js;
          pphot->i3p[ip] = static_cast<int>(pran->uniform()*nx3)+ks;
          weight_reduce = 1./static_cast<Real>(nx2);
          break;
        case BoundaryFace::outer_x2:
          pphot->i1p[ip] = static_cast<int>(pran->uniform()*nx1)+is;
          pphot->i2p[ip] = je;
          pphot->i3p[ip] = static_cast<int>(pran->uniform()*nx3)+ks;
          weight_reduce = 1./static_cast<Real>(nx2);
          break;
        case BoundaryFace::inner_x3:
          pphot->i1p[ip] = static_cast<int>(pran->uniform()*nx1)+is;
          pphot->i2p[ip] = static_cast<int>(pran->uniform()*nx2)+js;
          pphot->i3p[ip] = ks;
          weight_reduce = 1./static_cast<Real>(nx3);
          break;
        case BoundaryFace::outer_x3:
          pphot->i1p[ip] = static_cast<int>(pran->uniform()*nx1)+is;
          pphot->i2p[ip] = static_cast<int>(pran->uniform()*nx2)+js;
          pphot->i3p[ip] = ke;
          weight_reduce = 1./static_cast<Real>(nx3);
          break;
        default:
          std::stringstream msg;
          msg << "### FATAL ERROR in function [MonteCarloBlock::SetEmissionCellWeightArea]"
              << std::endl << "Face not valid" << std::endl;
          throw std::runtime_error(msg.str().c_str());
          break;
      }

      // Set weight according to the emission array, which is the relative number of
      // photons emitted in each cell. Note that emiss_to_weight assumes all cells
      // in block contribute so needs to be receduced by number of cell in direction
      // normal to the face
      pphot->wp[ip] = emission(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip])
        * emiss_to_weight * weight_reduce;
    } // end loop over ip
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetDensity()
//! \brief Make hard copy of density from MeshBlock to MonteCarloBlock.
//  Uses hard copy so that rho is always in cgs units

void MonteCarloBlock::GetDensity() {

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        rho(k,j,i) = rho_cgs * pmy_block->phydro->u(IDN,k,j,i);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetNumberDensity()
//! \brief default function for computing number densities if no user function provided.

void MonteCarloBlock::GetNumberDensity() {

  if (pmy_mc->UserGetNumberDensity != nullptr) {
    pmy_mc->UserGetNumberDensity(this);
    return;
  }

  Real heabund = 0.09; //hardcode for now (should be parameter)
  Real mp = 1.67262192369e-24;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real nh = rho(k,j,i) / (mp*(1.+4.*heabund));
        Real nhe = nh*heabund;
        nion(k,j,i) = nh + 4. * nhe;
        nel(k,j,i) = nh + 2. * nhe;
      }
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetScalars(MonteCarloBlock *pmcb)
//! \brief Make a hard copy of scalars from MeshBlock to MonteCarloBlock.

void MonteCarloBlock::GetScalars() {

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        scalars(k,j,i) = pmy_block->pscalars->s(0,k,j,i);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetVelocities()
//! \brief Make hard copy of velocites from MeshBlock to MonteCarloBlock.
//  Uses hard copy so that velocities is always fraction of speed of light

void MonteCarloBlock::GetVelocity() {

  if (GENERAL_RELATIVITY) {

    AthenaArray<Real> g, gi;
    g.NewAthenaArray(NMETRIC,ie+1);
    gi.NewAthenaArray(NMETRIC,ie+1);
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pmy_block->pcoord->CellMetric(k,j,is,ie,g,gi);
        for (int i=is; i<=ie; ++i) {
          Real alpha = 1.0/std::sqrt(-gi(I00,i));
          Real uu1 = pmy_block->phydro->u(IVX,k,j,i); 
          Real uu2 = pmy_block->phydro->u(IVY,k,j,i);
          Real uu3 = pmy_block->phydro->u(IVZ,k,j,i);
        
          Real gamma2 = 1. + g(I11,i)*uu1*uu1 + g(I22,i)*uu2*uu2 + g(I33,i)*uu3*uu3 +
                        2.0*g(I12,i)*uu1*uu2 + 2.*g(I13,i)*uu1*uu3 + 2.*g(I23,i)*uu2*uu3;
          Real gamma = std::sqrt(gamma2);
          
          vel(k,j,i,0) = -gamma*alpha*gi(I00,i);
          vel(k,j,i,1) = uu1 - gamma*alpha*gi(I01,i);
          vel(k,j,i,2) = uu2 - gamma*alpha*gi(I02,i);
          vel(k,j,i,3) = uu3 - gamma*alpha*gi(I03,i);
          Real beta0 = std::sqrt(SQR(vel(k,j,i,1)/vel(k,j,i,0))+SQR(vel(k,j,i,2)/vel(k,j,i,0))+SQR(vel(k,j,i,3)/vel(k,j,i,0)));
          if (beta0 >= betamax)
            printf("beta > betamax: %g\n", beta0);
        }
      }
    }
  } else {
    Real c_cgs = 2.99792458e10;
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real rho = pmy_block->phydro->u(IDN,k,j,i);
          vel(k,j,i,1) = vel_cgs * pmy_block->phydro->u(IM1,k,j,i) / (rho * c_cgs);
          vel(k,j,i,2) = vel_cgs * pmy_block->phydro->u(IM2,k,j,i) / (rho * c_cgs);
          vel(k,j,i,3) = vel_cgs * pmy_block->phydro->u(IM3,k,j,i) / (rho * c_cgs);
          Real beta0 = std::sqrt(SQR(vel(k,j,i,1))+SQR(vel(k,j,i,2))+SQR(vel(k,j,i,3)));
          Real beta = (beta0 > betamax) ? betamax : beta0;
          Real gamma = 1. / std::sqrt(1. - beta*beta);
          vel(k,j,i,0) = gamma;
          //gamma = 1;
          if (beta0 > 0.) {
            for (int l=1; l<4; ++l) {
              vel(k,j,i,l) *= gamma * beta / beta0 ;
            }
          }
        }
      }
    }
  } // end if (GENERAL_RELATIVITY) else
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::GetTemperature()
//! \brief default function for computing temperature if no user function provided.
//  Assumes EOS of the form P=RTd.

void MonteCarloBlock::GetTemperature() {

  if (pmy_mc->UserGetTemperature != nullptr) {
    pmy_mc->UserGetTemperature(this);
    return;
  }

  Real rideal = 8.314e7;
  Hydro* phydro = pmy_block->phydro;

  Real tconv;
  if (tgas_cgs <= 0.)
    tconv = 1. / rideal;
  else
    tconv = tgas_cgs;

  // compute temperature from pressure and density
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real temp = tconv * phydro->w(IEN,k,j,i) / phydro->w(IDN,k,j,i);
        // apply temperature floor
        temp = (temp > tfloor_cgs) ? temp : tfloor_cgs;
        temp = (temp < tceiling_cgs) ? temp : tceiling_cgs;
        tgas(k,j,i) = temp;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::ComputeTransformations()
//! \brief compute transformation matrices between comoving frame and lab frame

void MonteCarloBlock::ComputeTransformations() {

  // loop over all cells on block
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        boost_cmv(k,j,i,0,0) = vel(k,j,i,0);
        boost_lab(k,j,i,0,0) = vel(k,j,i,0);
        for (int m=1; m<4; m++) {
          boost_cmv(k,j,i,0,m) = -vel(k,j,i,m);
          //if (std::isnan(boost_cmv(k,j,i,0,m)))
          //    printf("boost: %d %d %d %g\n",k,j,i,vel(k,j,i,m));
          boost_lab(k,j,i,0,m) = vel(k,j,i,m);
        }
        for (int l=1; l<4; l++) {
          boost_cmv(k,j,i,l,0) = -vel(k,j,i,l);
          boost_lab(k,j,i,l,0) = vel(k,j,i,l);
          for (int m=1; m<4; m++) {
            boost_cmv(k,j,i,l,m) = vel(k,j,i,l)*vel(k,j,i,m)/(1.+vel(k,j,i,0));
            boost_lab(k,j,i,l,m) = boost_cmv(k,j,i,l,m);
          }
          boost_cmv(k,j,i,l,l) += 1.;
          boost_lab(k,j,i,l,l) += 1.;
        }
        /*for (int l=0; l<4; l++) {
          for (int m=0; m<4; m++) {
            Real sum = 0.;
            for (int n=0; n<4; n++) {
              sum += boost_cmv(k,j,i,l,n)*boost_lab(k,j,i,n,m);
            }
            printf("%d %d %g\n",l,m,sum);
          }
          }*/
      } // loop over i
    } // loop over j
  } // loop over k
}

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::TransformToComoving(Photon *pphot, int ips, int ipe)
//! \brief convert momentum vector to comoving (or tetrad) frame

void MonteCarloBlock::TransformToComoving(Photon *pphot, int ips, int ipe) {

  if (GENERAL_RELATIVITY) {
    for(int ip=ips; ip<=ipe; ip++) {
      // Construct the tetrad
      Real gcov[4][4];
      Real x[4];
      x[IMC0] = pphot->x0p[ip];
      x[IMC1] = pphot->x1p[ip];
      x[IMC2] = pphot->x2p[ip];
      x[IMC3] = pphot->x3p[ip];
      pcoord->Metric(x, gcov);

      // Create tetrad basis 
      Real ucon[4];
      ucon[IMC0] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],0);
      ucon[IMC1] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],1);
      ucon[IMC2] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],2);
      ucon[IMC3] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],3);
      Real econ[4][4], ecov[4][4];
      ConstructTetrad(ucon, gcov, econ, ecov);

      Real k0init = pphot->k0p[ip];
      // Transform to comoving tetrad
      Real kcopy[4];
      kcopy[IMC0] = pphot->k0p[ip];
      kcopy[IMC1] = pphot->k1p[ip];
      kcopy[IMC2] = pphot->k2p[ip];
      kcopy[IMC3] = pphot->k3p[ip];
      Real k[4];
      CoordinateToTetrad(kcopy, k, ecov);
      pphot->k0p[ip] = 1.;
      pphot->k1p[ip] = k[IMC1]/k[IMC0];
      pphot->k2p[ip] = k[IMC2]/k[IMC0];
      pphot->k3p[ip] = k[IMC3]/k[IMC0];

      Real nufact = k[IMC0]/k0init;
      //printf("com: %g %g %g %g %g\n",nufact,pphot->k0p[ip],pphot->k1p[ip],pphot->k2p[ip],pphot->k3p[ip]);
      pphot->ep[ip] *= nufact;
      pphot->acp[ip] /= nufact;
      pphot->scp[ip] /= nufact;
    }
  } else {
    for(int ip=ips; ip<=ipe; ip++) {

      int i1 = pphot->i1p[ip];
      int i2 = pphot->i2p[ip];
      int i3 = pphot->i3p[ip];

      Real k0init = pphot->k0p[ip];
      Real ki[4], kf[4];
      if (pmy_mc->general_pusher_flag) {
        ki[0] = pphot->k0p[ip];
        ki[1] = pphot->k1p[ip];
        ki[2] = pphot->k2p[ip];
        ki[3] = pphot->k3p[ip];

        Real x[4], invtet[4][4];
        x[0] = pphot->x0p[ip];
        x[1] = pphot->x1p[ip];
        x[2] = pphot->x2p[ip];
        x[3] = pphot->x3p[ip];
        pcoord->InverseTetrad(x,invtet);
        for (int j=0; j<4; j++) {
          kf[j] = 0.;
          for (int i=0; i<4; i++) {
            kf[j] += invtet[j][i] * ki[i];
          }
        }
        pphot->k0p[ip] = kf[0];
        pphot->k1p[ip] = kf[1];
        pphot->k2p[ip] = kf[2];
        pphot->k3p[ip] = kf[3];
      }
      if (boosts) {
        ki[0] = pphot->k0p[ip];
        ki[1] = pphot->k1p[ip];
        ki[2] = pphot->k2p[ip];
        ki[3] = pphot->k3p[ip];
        for (int j=0; j<4; j++) {
          kf[j] = 0.;
          for (int i=0; i<4; i++) {
            kf[j] += boost_cmv(i3,i2,i1,j,i) * ki[i];
          }
        }
        pphot->k0p[ip] = kf[0];
        pphot->k1p[ip] = kf[1];
        pphot->k2p[ip] = kf[2];
        pphot->k3p[ip] = kf[3];
      }

      //pphot->PrintPhoton("to com",ip);
      Real nufact = pphot->k0p[ip]/k0init;

      pphot->k0p[ip] = 1.;
      // SWD: maybe better to renormalize
      if ((COORDINATE_SYSTEM == "spherical_polar") && pmy_mc->polarized) {
        Real k3[3];
        Real cth = cos(pphot->x2p[ip]);
        Real sth = sin(pphot->x2p[ip]);
        Real cph = cos(pphot->x3p[ip]);
        Real sph = sin(pphot->x3p[ip]);
        k3[0] = pphot->k1p[ip];
        k3[1] = pphot->k2p[ip];
        k3[2] = pphot->k3p[ip];
        pphot->k1p[ip] = k3[0]*sth*cph + k3[1]*cth*cph - k3[2]*sph;
        pphot->k2p[ip] = k3[0]*sth*sph + k3[1]*cth*sph + k3[2]*cph;
        pphot->k3p[ip] = k3[0]*cth     - k3[1]*sth;
        Real knorm = std::sqrt(SQR(pphot->k1p[ip])+SQR(pphot->k2p[ip])+SQR(pphot->k3p[ip]));
        pphot->k1p[ip] /= knorm;
        pphot->k2p[ip] /= knorm;
        pphot->k3p[ip] /= knorm;
      } else {
        pphot->k1p[ip] = kf[1]/kf[0];
        pphot->k2p[ip] = kf[2]/kf[0];
        pphot->k3p[ip] = kf[3]/kf[0];
      }
      //if (std::isinf(nufact) || std::isnan(nufact))
      //  pphot->PrintPhoton("to com: nan in boost",ip);
      pphot->ep[ip] *= nufact;
      pphot->acp[ip] /= nufact;
      pphot->scp[ip] /= nufact;

    } //end loop overip
  } // end if GENERAL_RELATIVITY else

}


//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::TransformToCoordinate(Photon *pphot, int ips, int ipe)
//! \brief convert momentum vector to coordinate frame

void MonteCarloBlock::TransformToCoordinate(Photon *pphot, int ips, int ipe) {

  if (GENERAL_RELATIVITY) {
    for(int ip=ips; ip<=ipe; ip++) {
      Real gcov[4][4];
      Real x[4];
      x[IMC0] = pphot->x0p[ip];
      x[IMC1] = pphot->x1p[ip];
      x[IMC2] = pphot->x2p[ip];
      x[IMC3] = pphot->x3p[ip];
      pcoord->Metric(x, gcov);

      // Create tetrad basis 
      Real ucon[4];
      ucon[IMC0] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],0);
      ucon[IMC1] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],1);
      ucon[IMC2] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],2);
      ucon[IMC3] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],3);
      Real econ[4][4], ecov[4][4];
      ConstructTetrad(ucon, gcov, econ, ecov);

   
      // Transform to comoving tetrad
      Real kcopy[4];
      kcopy[IMC0] = pphot->k0p[ip] * pphot->ep[ip];
      kcopy[IMC1] = pphot->k1p[ip] * pphot->ep[ip];
      kcopy[IMC2] = pphot->k2p[ip] * pphot->ep[ip];
      kcopy[IMC3] = pphot->k3p[ip] * pphot->ep[ip];
      Real k0init = kcopy[IMC0];

      Real k[4];
      TetradToCoordinate(kcopy, k, ecov);
      pphot->k0p[ip] = k[IMC0];
      pphot->k1p[ip] = k[IMC1];
      pphot->k2p[ip] = k[IMC2];
      pphot->k3p[ip] = k[IMC3];

      Real nufact = pphot->k0p[ip]/k0init;
      //printf("coord: %g %g %g %g %g\n",nufact,pphot->k0p[ip],pphot->k1p[ip],pphot->k2p[ip],pphot->k3p[ip]);
      pphot->ep[ip] *= nufact;
      pphot->acp[ip] /= nufact;
      pphot->scp[ip] /= nufact;
    }
  } else {
    for(int ip=ips; ip<=ipe; ip++) {
      int i1 = pphot->i1p[ip];
      int i2 = pphot->i2p[ip];
      int i3 = pphot->i3p[ip];

      if (pmy_mc->general_pusher_flag) {
        pphot->k0p[ip] *= pphot->ep[ip];
        pphot->k1p[ip] *= pphot->ep[ip];
        pphot->k2p[ip] *= pphot->ep[ip];
        pphot->k3p[ip] *= pphot->ep[ip];
      } else if ((pmy_mc->polarized) && (COORDINATE_SYSTEM == "spherical_polar")) {
        // rotate cartesian to to spherical polar
        Real k3[3];
        Real cth = cos(pphot->x2p[ip]);
        Real sth = sin(pphot->x2p[ip]);
        Real cph = cos(pphot->x3p[ip]);
        Real sph = sin(pphot->x3p[ip]);
        k3[0] = pphot->k1p[ip];
        k3[1] = pphot->k2p[ip];
        k3[2] = pphot->k3p[ip];
        pphot->k1p[ip] = k3[0]*sth*cph + k3[1]*sth*sph + k3[2]*cth;
        pphot->k2p[ip] = k3[0]*cth*cph + k3[1]*cth*sph - k3[2]*sth;
        pphot->k3p[ip] = -k3[0]*sph + k3[1]*cph;
        Real knorm = std::sqrt(SQR(pphot->k1p[ip])+SQR(pphot->k2p[ip])+SQR(pphot->k3p[ip]));
        pphot->k1p[ip] /= knorm;
        pphot->k2p[ip] /= knorm;
        pphot->k3p[ip] /= knorm;
      }

      Real k0init = pphot->k0p[ip];
      Real ki[4], kf[4];
      if (boosts) {
        ki[0] = pphot->k0p[ip];
        ki[1] = pphot->k1p[ip];
        ki[2] = pphot->k2p[ip];
        ki[3] = pphot->k3p[ip];
        for (int j=0; j<4; j++) {
          kf[j] = 0.;
          for (int i=0; i<4; i++) {
            kf[j] += boost_lab(i3,i2,i1,j,i) * ki[i];
          }
        }
        pphot->k0p[ip] = kf[0];
        pphot->k1p[ip] = kf[1];
        pphot->k2p[ip] = kf[2];
        pphot->k3p[ip] = kf[3];
      }
      Real nufact;
      if (pmy_mc->general_pusher_flag) {
        ki[0] = pphot->k0p[ip];
        ki[1] = pphot->k1p[ip];
        ki[2] = pphot->k2p[ip];
        ki[3] = pphot->k3p[ip];

        Real x[4], tetrad[4][4];
        x[0] = pphot->x0p[ip];
        x[1] = pphot->x1p[ip];
        x[2] = pphot->x2p[ip];
        x[3] = pphot->x3p[ip];
        pcoord->Tetrad(x,tetrad);
        for (int j=0; j<4; j++) {
          kf[j] = 0.;
          for (int i=0; i<4; i++) {
            kf[j] += tetrad[j][i] * ki[i];
          }
        }
        pphot->k0p[ip] = kf[0];
        pphot->k1p[ip] = kf[1];
        pphot->k2p[ip] = kf[2];
        pphot->k3p[ip] = kf[3];
        nufact = kf[0]/k0init;
      } else {
        // spatial components of k are unit vectors
        nufact = pphot->k0p[ip]/k0init;

        Real knorm = std::sqrt(SQR(pphot->k1p[ip])+SQR(pphot->k2p[ip])+SQR(pphot->k3p[ip]));
        pphot->k0p[ip] = 1.;
        pphot->k1p[ip] /= knorm;
        pphot->k2p[ip] /= knorm;
        pphot->k3p[ip] /= knorm;
      }
      //if (std::isinf(nufact) || std::isnan(nufact))
      //  pphot->PrintPhoton("to cord: nan in boost",ip);
      // update energy and opacities
      pphot->ep[ip] *= nufact;
      pphot->acp[ip] /= nufact;
      pphot->scp[ip] /= nufact;

    } // for ip loop
  } // end if GENERAL_RELATIVITY else
}

//----------------------------------------------------------------------------------------
//! \fn Real MonteCarloBlock::FrequencyShiftComoving(Photon *pphot, int ip)
//! \brief frequency shift factor from coordinate to comoving

Real  MonteCarloBlock::FrequencyShiftComoving(Photon *pphot, int ip) {

  if (GENERAL_RELATIVITY) {
    Real gcov[4][4];
    Real x[4];
    x[IMC0] = pphot->x0p[ip];
    x[IMC1] = pphot->x1p[ip];
    x[IMC2] = pphot->x2p[ip];
    x[IMC3] = pphot->x3p[ip];
    pcoord->Metric(x, gcov);

    // Create tetrad basis 
    Real ucon[4];
    ucon[IMC0] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],0);
    ucon[IMC1] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],1);
    ucon[IMC2] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],2);
    ucon[IMC3] = vel(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip],3);
    Real econ[4][4], ecov[4][4];
    ConstructTetrad(ucon, gcov, econ, ecov);

    Real k0init = pphot->k0p[ip];
    // Transform to comoving tetrad
    Real kcopy[4];
    kcopy[IMC0] = pphot->k0p[ip];
    kcopy[IMC1] = pphot->k1p[ip];
    kcopy[IMC2] = pphot->k2p[ip];
    kcopy[IMC3] = pphot->k3p[ip];
    Real k[4];
    CoordinateToTetrad(kcopy, k, ecov);

    Real nufact = k[IMC0]/k0init;
    return nufact;
  } else {
    int i1 = pphot->i1p[ip];
    int i2 = pphot->i2p[ip];
    int i3 = pphot->i3p[ip];

    Real k0init = pphot->k0p[ip];
    Real ki[4], kf[4];
    if (tetrads) {
      ki[0] = pphot->k0p[ip];
      ki[1] = pphot->k1p[ip];
      ki[2] = pphot->k2p[ip];
      ki[3] = pphot->k3p[ip];

      Real x[4], invtet[4][4];
      x[0] = pphot->x0p[ip];
      x[1] = pphot->x1p[ip];
      x[2] = pphot->x2p[ip];
      x[3] = pphot->x3p[ip];
      pcoord->InverseTetrad(x,invtet);
      for (int j=0; j<4; j++) {
        kf[j] = 0.;
        for (int i=0; i<4; i++) {
          kf[j] += invtet[j][i] * ki[i];
        }
      }
    } else {
      kf[0] = pphot->k0p[ip];
      kf[1] = pphot->k1p[ip];
      kf[2] = pphot->k2p[ip];
      kf[3] = pphot->k3p[ip];
    }

    Real k0f = 0.;
    for (int i=0; i<4; i++) {
      k0f += boost_cmv(i3,i2,i1,0,i) * kf[i];
    }
    Real nufact =k0f/k0init;
    /*if (std::isinf(nufact) || std::isnan(nufact)) {
      printf("%g %g %g %g\n",boost_cmv(i3,i2,i1,0,0),boost_cmv(i3,i2,i1,0,1),
            boost_cmv(i3,i2,i1,0,2),boost_cmv(i3,i2,i1,0,3));
      printf("%g %g %g %g %g %g\n",k0init,k0f,kf[0],kf[1],kf[2],kf[3]);
      pphot->PrintPhoton("shift: nan in boost",ip);
      }*/
    return k0f/k0init;
  }

 }

//----------------------------------------------------------------------------------------
//! \fn void MonteCarloBlock::FrequencyAngelShiftComoving(Photon *pphot, int ip,
//             Real &shift, Real &k1, Real &k2, Real &k3) {
//! \brief frequency shift factor angular aberation from coordinate to comoving

void  MonteCarloBlock::FrequencyAngleShiftComoving(Photon *pphot, int ip, Real &shift,
                                                   Real &k1, Real &k2, Real &k3) {

    int i1 = pphot->i1p[ip];
    int i2 = pphot->i2p[ip];
    int i3 = pphot->i3p[ip];

    Real k0 = pphot->k0p[ip];
    Real ki[4], kf[4];
    if (tetrads) {
      ki[0] = pphot->k0p[ip];
      ki[1] = pphot->k1p[ip];
      ki[2] = pphot->k2p[ip];
      ki[3] = pphot->k3p[ip];

      Real x[4], invtet[4][4];
      x[0] = pphot->x0p[ip];
      x[1] = pphot->x1p[ip];
      x[2] = pphot->x2p[ip];
      x[3] = pphot->x3p[ip];
      pcoord->InverseTetrad(x,invtet);
      for (int j=0; j<4; j++) {
        kf[j] = 0.;
        for (int i=0; i<4; i++) {
          kf[j] += invtet[j][i] * ki[i];
        }
      }
    } else {
      kf[0] = pphot->k0p[ip];
      kf[1] = pphot->k1p[ip];
      kf[2] = pphot->k2p[ip];
      kf[3] = pphot->k3p[ip];
    }

    Real ke[4];
    for (int j=0; j<4; j++) {
      ke[j] = 0.;
      for (int i=0; i<4; i++) {
        ke[j] += boost_cmv(i3,i2,i1,j,i) * kf[i];
      }
    }
    //shift = k1 = k2 = k3 =1.;
    shift = ke[0]/k0;
    k1 = ke[1]/ke[0];
    k2 = ke[2]/ke[0];
    k3 = ke[3]/ke[0];

}

//----------------------------------------------------------------------------------------
//! \fn Real MonteCarloBlock::FrequencyShiftCoordinate(Photon *pphot, int ip)
//! \brief frequency shift factor from coordinate to comoving

Real  MonteCarloBlock::FrequencyShiftCoordinate(Photon *pphot, int ip) {

    int i1 = pphot->i1p[ip];
    int i2 = pphot->i2p[ip];
    int i3 = pphot->i3p[ip];

    Real ki[4];
    ki[0] = pphot->k0p[ip];
    ki[1] = pphot->k1p[ip];
    ki[2] = pphot->k2p[ip];
    ki[3] = pphot->k3p[ip];
    Real k0 = 0.;
    for (int i=0; i<4; i++) {
      k0 += boost_lab(i3,i2,i1,0,i) * ki[i];
    }
    return k0/ki[0];

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
