#ifndef POLARIZATION_HPP
#define POLARIZATION_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file polarization.hpp
//! \brief the polarization interface shared by the legacy and general pushers
//
// The polarized scattering routines in scattering.cpp do not work in an arbitrary basis.
// They require:
//
//   - the wavevector in cartesian components, and
//   - Q and U referenced to the meridian plane containing the global cartesian z axis
//     and k-hat.  ScatterThomsonPolarized reads that reference straight off the stored
//     components (mu = kz, phi = acos(kx/sin(theta))), so it is not a free choice.
//
// Both legacy pushers run only in flat spacetime (cartesian, spherical-polar).  A flat
// spacetime has a globally covariantly-constant frame -- the cartesian one -- and the
// convention above references Q and U to exactly that frame.  A photon between
// scatterings travels a straight line, so k-hat is constant in cartesian components, the
// meridian plane is constant, and Q and U do not change.  The stored components of k in
// spherical does change along the ray, which is why the wavevector has to be
// rotated at a scattering even though the Stokes parameters do not.
//
// In a curved spacetime no global frame exists, so the general pusher must parallel
// transport the coherency tensor step by step.
//
// Faraday rotation rotates Q into U continuously along the ray, so the Stokes parameters
// stop being constant between scatterings even in flat spacetime.

class MonteCarloBlock;
class Photon;

//! \brief rotate the photon into the basis the polarized scattering routines assume
void ToScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip);

//! \brief rotate the photon back out of the scattering basis
void FromScatteringBasis(MonteCarloBlock *pmcb, Photon *pphot, int ip);

#endif // POLARIZATION_HPP
