#ifndef PARTICLES_PARTICLE_BUFFER_HPP_
#define PARTICLES_PARTICLE_BUFFER_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particle_buffer.hpp
//! \brief defines ParticleBuffer class for communication of particles.
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include <complex>

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//--------------------------------------------------------------------------------------
//! \class ParticleBuffer
//! \brief defines the class for managing buffers for transporting particles.

class ParticleBuffer {
friend class Particles;
friend class Photon;
 public:
  // Class methods
  static void SetNumberOfProperties(int nint0, int nreal0, int ncplx0);

  // Constructors and destructor
  ParticleBuffer();
  explicit ParticleBuffer(int nparmax0);
  ~ParticleBuffer();

  // Instance method
  void Reallocate(int new_nparmax);

 protected:
  // Class variables
  static int nint;   //!> number of integer properties per particle
  static int nreal;  //!> number of real properties per particle
  static int ncplx;  //!> numbers of complex particle properties
  // Instance variables
  int* ibuf;   //!> ptr to integer buffer
  Real* rbuf;   //!> ptr to real buffer
  std::complex<Real>* cbuf; //!> ptr to complex buffer

  int nparmax;  //!> maximum number of particles
  int npar;     //!> actual number of particles in the buffer
#ifdef MPI_PARALLEL
  MPI_Request reqi, reqr, reqc;  //!> MPI request handles
  bool mpi_active;         //!> Whether or not receiving process is initiated
  int flagn;               //!> Flag indicating if the incoming number is known
  int flagi, flagr, flagc; //!> Flags indicating if the respective buffer is filled
  int tag;                 //!> MPI tag (allowing for from tag to tag + 2)
#endif
};
#endif  // PARTICLES_PARTICLE_BUFFER_HPP_
