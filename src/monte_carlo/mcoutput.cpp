//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mcoutput.cpp
//  \brief implementation of functions in class MCOutput

// C++ headers
#include <stdio.h>
#include <stdlib.h>

// Athena++ headers
#include "montecarlo.hpp"
#include "mcoutput.hpp"

// constructor
Spectrum::Spectrum(MomentumRange input_range, bool pol) {

  next = NULL;
  face = FACE_UNDEF;
  range = input_range;
  polarized = pol;
  
 
  // Allocate and intialize energy bins
  energies.NewAthenaArray(range.ne+1);
  BuildFrequencyGrid(range.emin,range.emax,range.ne);

  // Allocate arrays for intensities
  intensity.NewAthenaArray(range.nphi,range.ncth,range.ne);
  intensity_sq.NewAthenaArray(range.nphi,range.ncth,range.ne);
  if (polarized) {
    stokesq.NewAthenaArray(range.nphi,range.ncth,range.ne);
    stokesq_sq.NewAthenaArray(range.nphi,range.ncth,range.ne);
    stokesu.NewAthenaArray(range.nphi,range.ncth,range.ne);
    stokesu_sq.NewAthenaArray(range.nphi,range.ncth,range.ne);
  }
}

// constructor
Spectrum::Spectrum(Spectrum *pspec) {

  next = NULL;
  range = pspec->range;
  polarized = pspec->polarized;
  face = pspec->face;
 
  // Allocate and intialize energy bins
  energies.NewAthenaArray(range.ne+1);
  BuildFrequencyGrid(range.emin,range.emax,range.ne);

  // Allocate arrays for intensities
  intensity.NewAthenaArray(range.nphi,range.ncth,range.ne);
  intensity_sq.NewAthenaArray(range.nphi,range.ncth,range.ne);
  if (polarized) {
    stokesq.NewAthenaArray(range.nphi,range.ncth,range.ne);
    stokesq_sq.NewAthenaArray(range.nphi,range.ncth,range.ne);
    stokesu.NewAthenaArray(range.nphi,range.ncth,range.ne);
    stokesu_sq.NewAthenaArray(range.nphi,range.ncth,range.ne);
  }
}
// destructor
Spectrum::~Spectrum() {

  energies.DeleteAthenaArray();
  intensity.DeleteAthenaArray();
  intensity_sq.DeleteAthenaArray();
  if (polarized) {
    stokesq.DeleteAthenaArray();
    stokesq_sq.DeleteAthenaArray();
    stokesu.DeleteAthenaArray();
    stokesu_sq.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Spectrum::BuildFrequencyGrid(Real emin, Real emax, int nfreq)
//  \brief initialize frequency bins

void Spectrum::BuildFrequencyGrid(Real emin, Real emax, int nfreq) {

  Real de = log10(emax/emin) / static_cast<Real>(nfreq);
  energies(0) = log10(emin);
  for(int i=0; i<nfreq; ++i) {
    energies(i+1) = energies(i) + de;
  }
  for(int i=0; i<=nfreq; ++i){
    energies(i) = exp(2.30258509299*energies(i));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Spectrum::SetSurface(std::string input_face)
//  \brief set corresponding surface

void Spectrum::SetSurface(std::string input_face) {

  if (input_face == "inner_x1") {
    face = INNER_X1;
  } else if (input_face == "outer_x1") {
    face = OUTER_X1;
  } else if (input_face == "inner_x2") {
    face = INNER_X2;
  } else if (input_face == "outer_x2") {
    face = OUTER_X2;
  } else if (input_face == "inner_x3") {
    face = INNER_X3;
  } else if (input_face == "outer_x3") {
    face = OUTER_X3;
  } else {
    std::cout << "Warning: face not set correctly in output spectrum: " << input_face
              << ", leaving undefined." << std::endl;
  }
}


//----------------------------------------------------------------------------------------
//! \fn void Spectrum::UpdateSpectrum(Photon *pphot)
//  \brief add photon contribution to spectrum
//  *** should be more general to explitly account for non-cartesian possibilities

void Spectrum::UpdateSpectrum(Photon *pphot) {

  Real mu = fabs(pphot->k[2]); //CARTESIAN ONLY
  int ebin,mubin,phibin;
  Real phi;
  int nphi = range.nphi;
  int nmu = range.ncth;

  // Get ebin
  ebin = GetEbin(pphot->energy);
  if (ebin < 0) return;
  // Get phi bin
  if (mu == 1.0) {
    phibin = 0;
  } else {
    Real stheta = sqrt(SQR(pphot->k[0])+SQR(pphot->k[1]));
    phi = acos(pphot->k[0]/stheta);
    if(pphot->k[1] < 0.0)
      phi = 2 * PI - phi;
    phibin = static_cast<int>(phi * static_cast<Real>(nphi) / (2.*PI));
  }
  if(phibin >= nphi) {
    std::cout << "Warning: phibin > nphi (phibin, phi, kx, ky): " << phibin << ' '
              << phi << ' ' << pphot->k[0] << ' ' << pphot->k[1] << std::endl;
    phibin = nphi-1;
  }
  if(phibin < 0) {
    std::cout << "Warning: phibin < 0 (phibin, phi, kx, ky): " << phibin << ' '
              << phi << ' ' << pphot->k[0] << ' ' << pphot->k[1] << std::endl;
    phibin = 0;
  }
  // Get mu bin
  mubin = static_cast<int>(mu * static_cast<Real>(nmu) );
  if(mubin >= nmu) {
    std::cout << "Warning: mubin > nmu (mubin, mu, kz): " << mubin << ' '
              << mu << ' ' << pphot->k[2]  << std::endl;
    mubin = nmu-1;
  }

  Real& weight = pphot->weight;
  weight *= pphot->eweight;
  if ((isinf(weight)) || (isnan(weight))) {
    std::cout << "Warning: weight is Nan or Inf: " << weight << std::endl;
  } else {
    //std::cout << "final: " << weight << std::endl;
    intensity(phibin,mubin,ebin) += pphot->stokes[0] * weight;
    intensity_sq(phibin,mubin,ebin) += SQR(pphot->stokes[0] * weight);  
    if (polarized) {
      stokesq(phibin,mubin,ebin) += pphot->stokes[1] * weight;
      stokesq_sq(phibin,mubin,ebin) += SQR(pphot->stokes[1] * weight); 
      stokesu(phibin,mubin,ebin) += pphot->stokes[2] * weight;
      stokesu_sq(phibin,mubin,ebin) += SQR(pphot->stokes[2] * weight);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Spectrum::GetEbin(Real energy)
//  \brief return bin number corresponding to energy

int Spectrum::GetEbin(Real energy)
{

  if ( (energy < energies(0)) || (energy > energies(range.ne)) )
    return -1;

  // Perform binary search
  int low = 0;
  int high = range.ne+1;
  int mid;
  while(low <= high) {
    mid = (low + high) / 2;
    if(energies(mid) <= energy) {
      if(energies(mid+1) > energy)
        return mid;
      else
        low = mid+1;
    } else
      high = mid-1;
  }
  if(mid == high)
    return mid;
  else {
    return -1;
    std::cout << "Warning: binary search failed in ebin: " << energy << ' '
              << mid <<  std::endl;
  }

}

// constructor
MCOutput::MCOutput(MonteCarlo *pmc, ParameterInput *pin) {

  pmy_mc = pmc;
  std::stringstream msg;
  InputBlock *pib = pin->pfirst_block;

  moments = false;
  pspec = NULL;
  // loop over input block names.  Find those that start with "output", read parameters,
  // and construct linked list of spectra if present, set moments flag if moments output
  // present
  Spectrum *pfirst = NULL, *plast;
  while (pib != NULL) {
    if (pib->block_name.compare(0,6,"output") == 0) {
      // Look for spectra
      std::string type = pin->GetString(pib->block_name,"file_type");
      if (type.compare("spec") == 0) {
        MomentumRange range;
        range.ne = pin->GetInteger(pib->block_name,"ne");
        Real emin = pin->GetReal("montecarlo","emin");
        Real emax = pin->GetReal("montecarlo","emax");
        Real everg = 1.6021772e-12;
        range.emin = everg * pin->GetOrAddReal(pib->block_name,"emin",emin);
        range.emax = everg * pin->GetOrAddReal(pib->block_name,"emax",emax);
        range.nphi = pin->GetOrAddInteger(pib->block_name,"nphi",8);
        range.phimin = pin->GetOrAddReal(pib->block_name,"phimin",0.);
        range.phimax = pin->GetOrAddReal(pib->block_name,"phimax",2.*PI);
        range.ncth = pin->GetOrAddInteger(pib->block_name,"ncth",8);
        range.cthmin = pin->GetOrAddReal(pib->block_name,"cthmin",0.);
        range.cthmax = pin->GetOrAddReal(pib->block_name,"cthmax",1.); 
        bool polarized = pin->GetOrAddBoolean(pib->block_name,"polarized",pmc->polarized);
        pspec = new Spectrum(range,polarized);
        std::string face = pin->GetOrAddString(pib->block_name,"face","absent");
        if (face.compare("absent") != 0) {
          pspec->SetSurface(face);
        }
        if (pfirst == NULL)
          pfirst = pspec;
        else
          plast->next = pspec;
        plast = pspec;
      } else {
        // Look for moments
        std::string var = pin->GetOrAddString(pib->block_name,"variable","absent");
        if (var.compare("mcmom") == 0 || var.compare("Ermc") == 0 ||
            var.compare("Frmc") == 0 || var.compare("Prmc") == 0) {
          moments = true;
        }
      }
    }
    pib = pib->pnext;  // move to next input block name
  }
  // set pspec to head node of specrum list
  pspec = pfirst;
 
}

// destructor
MCOutput::~MCOutput() {

}

//----------------------------------------------------------------------------------------
//! \fn void MCOutput::OutputSpectrum(MonteCarlo *pmc
//  \brief output intensity spectrum

void MCOutput::OutputSpectrum(Spectrum *pspec, Real norm, std::string outfile) {

  FILE *of_ptr;
  
  Real everg = 1.6021772e-12;  
  Real emin = pspec->range.emin / everg; // output in eV
  Real emax = pspec->range.emax / everg; // output in eV
  if ((of_ptr=fopen("intens_sums.out","w")) != NULL) {
    fprintf(of_ptr,"%d %d %d %g\n",pspec->range.ne,pspec->range.ncth,pspec->range.nphi,norm);
    fprintf(of_ptr,"%lG %lG %lG\n",everg,emin,emax);

    // Output intensity data at top of domain
    for(int k=0; k<pspec->range.ne; ++k) {
      for(int j=0; j<pspec->range.ncth; ++j) {
        for(int i=0; i<pspec->range.nphi; ++i) {
          fprintf(of_ptr,"%G %G ",pspec->intensity(i,j,k),
                  pspec->intensity_sq(i,j,k));
          if (pspec->polarized) {
            fprintf(of_ptr,"%G %G %G %G\n",pspec->stokesq(i,j,k),
                    pspec->stokesq_sq(i,j,k),pspec->stokesu(i,j,k),
                    pspec->stokesu_sq(i,j,k));
          } else {
            fprintf(of_ptr,"\n");
          }

        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MCOutput::OutputSpectra(MonteCarlo *pmc)
//  \brief output all intensity spectra over all blocks

void MCOutput::OutputSpectra(MonteCarlo *pmc) {
  
  if (pspec == NULL)
    return;

  std::string outfile;
  
  MonteCarloBlock *pmcb = pmc->pblock;
  if (pmcb == NULL)
    return;
  do {
    Real norm = static_cast<Real>(pmc->nphot)/ static_cast<Real>(pmc->ncells);    
    //for (int i=0; i<pmcb->nspec; ++i) {
    //  OutputSpectrum(pmcb->pspec[i],nphot/ntot,outfile);
    //}
    OutputSpectrum(pmcb->pspec,norm,outfile);
    pmcb = pmcb->next;
  } while (pmcb != NULL);

}
