
// Athena++ astrophysical MHD code
//=======================================================================================
//! \file particle.cpp
//  \brief defines functions from particle class
//=======================================================================================

#include <stdio.h>
#include <iostream>
#include <time.h>

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../field/field.hpp"
#include "../parameter_input.hpp"
#include "../globals.hpp"

// Primary header
#include "particle.hpp"
#include "../hybrid/hybrid.hpp"

class Hybrid;

void Particle::Move()
{
  MeshBlock *pmb = pmy_block;
  Real wei1[3], wei2[3], wei3[3];
  Real isgr, jsgr, ksgr;
  Real qom = 1.0;
  Real dt = pmy_block->pmy_mesh->dt;
  Real x1n, x2n, x3n;
  Real v1n, v2n, v3n;
  Real vp1, vp2, vp3;
  Real s1, s2, s3;
  Real a, d;
  int ig;
  int is, js, ks;
  int i0, j0, k0;
  Real w, totwei, totwei2;
  Real bfld1, bfld2, bfld3, efld1, efld2, efld3, eb;
  Real bsq, b, diff, edotb;

// start of Bei's modifications
  Real x1tmp,x2tmp,x3tmp,v1tmp,v2tmp,v3tmp;
  Real halfdt;
  Real bfld1v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real bfld2v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real bfld3v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real efld1v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real efld2v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real efld3v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real ebv[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
// end of Bei's modificatons

  Real x1s, x2s, x3s, dx1, dx2, dx3;
  x1s = pmy_block->block_size.x1min;
  x2s = pmy_block->block_size.x2min;
  x3s = pmy_block->block_size.x3min;
  dx1 = pmy_block->block_size.nx1/(pmy_block->block_size.x1max-pmy_block->block_size.x1min);
  dx2 = pmy_block->block_size.nx2/(pmy_block->block_size.x2max-pmy_block->block_size.x2min);
  dx3 = pmy_block->block_size.nx3/(pmy_block->block_size.x3max-pmy_block->block_size.x3min);

  AthenaArray<Real> fcoup;
  fcoup.InitWithShallowCopy(pmy_block->phybrid->fcoup);
  isgr = (Real)pmb->is;
  jsgr = (Real)pmb->js;
  ksgr = (Real)pmb->ks;
  halfdt = 0.5 * dt;
  totwei = qom * dt * 0.5;
  totwei2 = totwei * totwei;
  for (long p=0; p<nparticle; p+=SIMD_WIDTH)
  {
    int remains = nparticle - p;

    Real *__restrict__ x1p = x1 + p;
    Real *__restrict__ x2p = x2 + p;
    Real *__restrict__ x3p = x3 + p;
    Real *__restrict__ v1p = v1 + p;
    Real *__restrict__ v2p = v2 + p;
    Real *__restrict__ v3p = v3 + p;


// gather fields
#pragma omp simd aligned(x1p,x2p,x3p,v1p,v2p,v3p : CACHELINE_BYTES) simdlen(SIMD_WIDTH) private(x1tmp,x2tmp,x3tmp,v1tmp,v2tmp,v3tmp,x1n,x2n,x3n,a,ig,is,js,ks,d,wei1,wei2,wei3,k0,j0,i0)
    for (int pp=0; pp<std::min(SIMD_WIDTH,remains); pp++) {

      x1tmp = x1p[pp];
      x2tmp = x2p[pp];
      x3tmp = x3p[pp];
      v1tmp = v1p[pp];
      v2tmp = v2p[pp];
      v3tmp = v3p[pp];

      x1n = x1tmp + v1tmp * halfdt;
      x2n = x2tmp + v2tmp * halfdt;
      x3n = x3tmp + v3tmp * halfdt;

      a = (x1tmp - x1s) * dx1 + isgr;
      ig = (int)(a);
      is = ig - 1;
      d = a - ig;
      wei1[0] = 0.5 * (1.0 - d) * (1.0 - d);
      wei1[1] = 0.75 - (d - 0.5) * (d - 0.5);
      wei1[2] = 0.5 * d * d;

      a = (x2tmp - x2s) * dx2 + jsgr;
      ig = (int)(a);
      js = ig - 1;
      d = a - ig;
      wei2[0] = 0.5 * (1.0 - d) * (1.0 - d);
      wei2[1] = 0.75 - (d - 0.5) * (d - 0.5);
      wei2[2] = 0.5 * d * d;

      a = (x3tmp - x3s) * dx3 + ksgr;
      ig = (int)(a);
      ks = ig - 1;
      d = a - ig;
      wei3[0] = 0.5 * (1.0 - d) * (1.0 - d);
      wei3[1] = 0.75 - (d - 0.5) * (d - 0.5);
      wei3[2] = 0.5 * d * d;
      
      bfld1v[pp] = 0.0; bfld2v[pp] = 0.0; bfld3v[pp] = 0.0;
      efld1v[pp] = 0.0; efld2v[pp] = 0.0; efld3v[pp] = 0.0;
      ebv[pp]    = 0.0;
    
      for (k0=0; k0<=2; k0++){
        for (j0=0; j0<=2; j0++){
          for (i0=0; i0<=2; i0++){
            w = wei3[k0] * wei2[j0] * wei1[i0];
            bfld1v[pp] += w * fcoup(ks+k0,js+j0,is+i0,IB1);
            bfld2v[pp] += w * fcoup(ks+k0,js+j0,is+i0,IB2);
            bfld3v[pp] += w * fcoup(ks+k0,js+j0,is+i0,IB3);
            efld1v[pp] += w * fcoup(ks+k0,js+j0,is+i0,IE1);
            efld2v[pp] += w * fcoup(ks+k0,js+j0,is+i0,IE2);
            efld3v[pp] += w * fcoup(ks+k0,js+j0,is+i0,IE3);
            ebv[pp]    += w * fcoup(ks+k0,js+j0,is+i0,IEB);
          }
        }
      }
    }

// push
#pragma omp simd simdlen(SIMD_WIDTH) aligned(x1p,x2p,x3p,v1p,v2p,v3p : CACHELINE_BYTES) private(bsq,edotb,diff,x1n,x2n,x3n,v1n,v2n,v3n,x1tmp,x2tmp,x3tmp,v1tmp,v2tmp,v3tmp,vp1,vp2,vp3,b,s1,s2,s3)
    for (int pp=0; pp<std::min(SIMD_WIDTH,remains); pp++) { 
      x1tmp = x1p[pp];
      x2tmp = x2p[pp];
      x3tmp = x3p[pp];
      v1tmp = v1p[pp];
      v2tmp = v2p[pp];
      v3tmp = v3p[pp];

      x1n = x1tmp + v1tmp * halfdt;
      x2n = x2tmp + v2tmp * halfdt;
      x3n = x3tmp + v3tmp * halfdt;
    
      bsq = std::max(bfld1v[pp] * bfld1v[pp] + bfld2v[pp] * bfld2v[pp] 
                   + bfld3v[pp] * bfld3v[pp], TINY_NUMBER);
      edotb = efld1v[pp] * bfld1v[pp] + efld2v[pp] * bfld2v[pp] 
            + efld3v[pp] * bfld3v[pp];
      diff = (ebv[pp] - edotb)/bsq;
      efld1v[pp] += diff*bfld1v[pp];
      efld2v[pp] += diff*bfld2v[pp];
      efld3v[pp] += diff*bfld3v[pp];

      bfld1v[pp] *= totwei; bfld2v[pp] *= totwei; bfld3v[pp] *= totwei;
      efld1v[pp] *= totwei; efld2v[pp] *= totwei; efld3v[pp] *= totwei;
      b = 2.0 / (1 + bsq * totwei2);

      v1n = v1tmp + efld1v[pp];
      v2n = v2tmp + efld2v[pp];
      v3n = v3tmp + efld3v[pp];

      vp1 = v1n + v2n * bfld3v[pp] - v3n * bfld2v[pp];
      vp2 = v2n + v3n * bfld1v[pp] - v1n * bfld3v[pp];
      vp3 = v3n + v1n * bfld2v[pp] - v2n * bfld1v[pp];

      s1 = bfld1v[pp] * b;
      s2 = bfld2v[pp] * b;
      s3 = bfld3v[pp] * b;

      v1tmp = v1n + efld1v[pp] + vp2 * s3 - vp3 * s2;
      v2tmp = v2n + efld2v[pp] + vp3 * s1 - vp1 * s3;
      v3tmp = v3n + efld3v[pp] + vp1 * s2 - vp2 * s1;
      x1p[pp] = x1n + v1tmp * halfdt;
      x2p[pp] = x2n + v2tmp * halfdt;
      x3p[pp] = x3n + v3tmp * halfdt;
      v1p[pp] = v1tmp;
      v2p[pp] = v2tmp;
      v3p[pp] = v3tmp;
    }
  }
  return;
}
