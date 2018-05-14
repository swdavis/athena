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
Spectrum::Spectrum(Real elow, Real ehigh, int n0, int n1, int n2, 
                   bool pol) {

  emin = elow;
  emax = ehigh;
  nfreq = n0;
  nmu = n1;
  nphi = n2;
  polarized = pol;
  
  // Allocate and intialize energy bins
  energies.NewAthenaArray(nfreq+1);
  BuildFrequencyGrid(emin,emax,nfreq);

  // Allocate arrays for intensities
  intensity.NewAthenaArray(nphi,nmu,nfreq);
  intensity_sq.NewAthenaArray(nphi,nmu,nfreq);
  if (polarized) {
    stokesq.NewAthenaArray(nphi,nmu,nfreq);
    stokesq_sq.NewAthenaArray(nphi,nmu,nfreq);
    stokesu.NewAthenaArray(nphi,nmu,nfreq);
    stokesu_sq.NewAthenaArray(nphi,nmu,nfreq);
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
//! \fn void Spectrum::UpdateSpectrum(Photon *pphot)
//  \brief add photon contribution to spectrum
//  *** should be more general to explitly account for non-cartesian possibilities

void Spectrum::UpdateSpectrum(Photon *pphot) {

  Real mu = fabs(pphot->k[2]); //CARTESIAN ONLY
  int ebin,mubin,phibin;
  Real phi;

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

  if ( (energy < energies(0)) || (energy > energies(nfreq)) )
    return -1;

  // Perform binary search
  int low = 0;
  int high = nfreq+1;
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
  Real emin = pspec->emin / everg; // output in eV
  Real emax = pspec->emax / everg; // output in eV
  if ((of_ptr=fopen("intens_sums.out","w")) != NULL) {
    fprintf(of_ptr,"%d %d %d %g\n",pspec->nfreq,pspec->nmu,pspec->nphi,norm);
    fprintf(of_ptr,"%lG %lG %lG\n",everg,emin,emax);

    // Output intensity data at top of domain
    for(int k=0; k<pspec->nfreq; ++k) {
      for(int j=0; j<pspec->nmu; ++j) {
        for(int i=0; i<pspec->nphi; ++i) {
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
  
  std::string outfile;

  MonteCarloBlock *pmcb = pmc->pblock;
  do {
    Real ntot = static_cast<Real>(pmcb->ncells);
    Real nphot = static_cast<Real>(pmcb->nphot);
    //for (int i=0; i<pmcb->nspec; ++i) {
    //  OutputSpectrum(pmcb->pspec[i],nphot/ntot,outfile);
    //}
    OutputSpectrum(pmcb->pspec,nphot/ntot,outfile);
    pmcb = pmcb->next;
  } while (pmcb != NULL);

}
