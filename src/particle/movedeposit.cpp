
// Athena++ astrophysical MHD code
//=======================================================================================
//! \file particle.cpp
//  \brief defines functions from particle class
//=======================================================================================

#include <stdio.h>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <cstdlib>

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

void Particle::MoveDeposit()
{
  MeshBlock *pmb = pmy_block;
  Real wei1[3], wei2[3], wei3[3];
  int isg, jsg, ksg;
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

// for mover  
  Real x1tmp,x2tmp,x3tmp,v1tmp,v2tmp,v3tmp,wei3tmp;
  Real halfdt;
  Real bfld1v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real bfld2v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real bfld3v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real efld1v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real efld2v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real efld3v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real ebv[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));

// predicted positions
  Real x1pred[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real x2pred[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real x3pred[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real v1pred[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real v2pred[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real v3pred[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));

  Real x1s, x2s, x3s, dx1, dx2, dx3;
  x1s = pmy_block->block_size.x1min;
  x2s = pmy_block->block_size.x2min;
  x3s = pmy_block->block_size.x3min;
  dx1 = pmy_block->block_size.nx1/(pmy_block->block_size.x1max-pmy_block->block_size.x1min);
  dx2 = pmy_block->block_size.nx2/(pmy_block->block_size.x2max-pmy_block->block_size.x2min);
  dx3 = pmy_block->block_size.nx3/(pmy_block->block_size.x3max-pmy_block->block_size.x3min);

  AthenaArray<Real> mcoup, fcoup;
  mcoup.InitWithShallowCopy(pmy_block->phybrid->mcoup);
  fcoup.InitWithShallowCopy(pmy_block->phybrid->fcoup);

  isg = pmb->is;
  jsg = pmb->js;
  ksg = pmb->ks;
  isgr = (Real)pmb->is;
  jsgr = (Real)pmb->js;
  ksgr = (Real)pmb->ks;

  halfdt = 0.5 * dt;
  totwei = qom * dt * 0.5;
  totwei2 = totwei * totwei;

// for deposit
  int nx1=pmy_block->block_size.nx1;
  int nx2=pmy_block->block_size.nx2;
  int nx3=pmy_block->block_size.nx3;
  int nx1g=nx1+2*NGHOST;
  int nx2g=nx2+2*NGHOST;
  int nx3g=nx3+2*NGHOST;
  Real wei1v[3][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real wei2v[3][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real wei3v[3][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real weitotv[9][SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  int ind1[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  int ind2[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  int ind3[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  Real *mcoupv;
  int error = posix_memalign((void **)&mcoupv, CACHELINE_BYTES, sizeof(Real)*8*NMCOUP*nx1g*nx2g*nx3g);
  if (error!=0) {
    printf("failed to allocate memory for mcoupv");
    exit(1);
  }

#pragma omp simd aligned(mcoupv: CACHELINE_BYTES) simdlen(SIMD_WIDTH)
  for (int i=0; i<8*4*nx1g*nx2g*nx3g; i++){
    mcoupv[i] = 0.0;
  }

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
      x1pred[pp] = x1n + v1tmp * halfdt;
      x2pred[pp] = x2n + v2tmp * halfdt;
      x3pred[pp] = x3n + v3tmp * halfdt;
      v1pred[pp] = v1tmp;
      v2pred[pp] = v2tmp;
      v3pred[pp] = v3tmp;
    }

// gather weights
#pragma omp simd aligned(x1p,x2p,x3p,v1p,v2p,v3p : CACHELINE_BYTES) simdlen(SIMD_WIDTH) private(x1tmp,x2tmp,x3tmp,v1tmp,v2tmp,v3tmp,a,ig,is,js,ks,d,k0,j0,i0)
    for (int pp=0; pp<std::min(SIMD_WIDTH,remains); pp++) {
      x1tmp = x1pred[pp];
      x2tmp = x2pred[pp];
      x3tmp = x3pred[pp];
      v1tmp = v1pred[pp];
      v2tmp = v2pred[pp];
      v3tmp = v3pred[pp];
      a = (x1tmp - x1s) * dx1 + isgr;
      ig = (int)(a);
      is = ig - 1;
      d = a - ig;
      wei1v[0][pp] = 0.5 * (1.0 - d) * (1.0 - d);
      wei1v[1][pp] = 0.75 - (d - 0.5) * (d - 0.5);
      wei1v[2][pp] = 0.5 * d * d;
      ind1[pp] = is - isg + NGHOST;

      a = (x2tmp - x2s) * dx2 + jsgr;
      ig = (int)(a);
      js = ig - 1;
      d = a - ig;
      wei2v[0][pp] = 0.5 * (1.0 - d) * (1.0 - d);
      wei2v[1][pp] = 0.75 - (d - 0.5) * (d - 0.5);
      wei2v[2][pp] = 0.5 * d * d;
      ind2[pp] = js - jsg + NGHOST;

      a = (x3tmp - x3s) * dx3 + ksgr;
      ig = (int)(a);
      ks = ig - 1;
      d = a - ig;
      wei3v[0][pp] = 0.5 * (1.0 - d) * (1.0 - d);
      wei3v[1][pp] = 0.75 - (d - 0.5) * (d - 0.5);
      wei3v[2][pp] = 0.5 * d * d;
      ind3[pp] = ks - ksg + NGHOST;

      weitotv[0][pp] = wei1v[0][pp]*wei2v[0][pp];
      weitotv[1][pp] = wei1v[1][pp]*wei2v[0][pp];
      weitotv[2][pp] = wei1v[2][pp]*wei2v[0][pp];
      weitotv[3][pp] = wei1v[0][pp]*wei2v[1][pp];
      weitotv[4][pp] = wei1v[1][pp]*wei2v[1][pp];
      weitotv[5][pp] = wei1v[2][pp]*wei2v[1][pp];
      weitotv[6][pp] = wei1v[0][pp]*wei2v[2][pp];
      weitotv[7][pp] = wei1v[1][pp]*wei2v[2][pp];
      weitotv[8][pp] = wei1v[2][pp]*wei2v[2][pp];
    }

// deposit
    for (int pp=0; pp<std::min(SIMD_WIDTH,remains); pp++) {
      int index = ind1[pp] + ind2[pp]*nx1g + ind3[pp]*nx1g*nx2g;
      Real *mcoupvp = mcoupv+index*32;
      v1tmp = v1pred[pp];
      v2tmp = v2pred[pp];
      v3tmp = v3pred[pp];
      for (int k=0; k<=2; k++){
        wei3tmp = wei3v[k][pp];
#pragma omp simd aligned(mcoupvp) simdlen(8)
        for (int ij=0; ij<8; ij++){
          Real weight = wei3tmp*weitotv[ij][pp];
          mcoupvp[ij] += weight;
          mcoupvp[8+ij] += weight*v1tmp;
          mcoupvp[16+ij]+= weight*v2tmp;
          mcoupvp[24+ij]+= weight*v3tmp;
        }
        Real weight8 = wei3tmp*weitotv[8][pp];
        int index8 = index + NGHOST + NGHOST*nx1g + k*nx1g*nx2g;
        mcoupv[index8*32] += weight8;
        mcoupv[index8*32+8] += weight8*v1tmp;
        mcoupv[index8*32+16] += weight8*v2tmp;
        mcoupv[index8*32+24] += weight8*v3tmp;
      }
    }
  }

  for (int k=pmy_block->ks-NGHOST; k<=pmy_block->ke; k++) {
    for (int j=pmy_block->js-NGHOST; j<=pmy_block->je; j++) {
#pragma omp simd
      for (int i=pmy_block->is-NGHOST; i<=pmy_block->ie; i++) {
        int index = (i-isg+NGHOST)+(j-jsg+NGHOST)*nx1g+(k-ksg+NGHOST)*nx1g*nx2g;
        mcoup(k,j,i,IDN) += mcoupv[index*32];
        mcoup(k,j,i+1,IDN) += mcoupv[index*32+1];
        mcoup(k,j,i+2,IDN) += mcoupv[index*32+2];
        mcoup(k,j+1,i,IDN) += mcoupv[index*32+3];
        mcoup(k,j+1,i+1,IDN) += mcoupv[index*32+4];
        mcoup(k,j+1,i+2,IDN) += mcoupv[index*32+5];
        mcoup(k,j+2,i,IDN) += mcoupv[index*32+6];
        mcoup(k,j+2,i+1,IDN) += mcoupv[index*32+7];

        mcoup(k,j,i,IM1) += mcoupv[index*32+8];
        mcoup(k,j,i+1,IM1) += mcoupv[index*32+9];
        mcoup(k,j,i+2,IM1) += mcoupv[index*32+10];
        mcoup(k,j+1,i,IM1) += mcoupv[index*32+11];
        mcoup(k,j+1,i+1,IM1) += mcoupv[index*32+12];
        mcoup(k,j+1,i+2,IM1) += mcoupv[index*32+13];
        mcoup(k,j+2,i,IM1) += mcoupv[index*32+14];
        mcoup(k,j+2,i+1,IM1) += mcoupv[index*32+15];

        mcoup(k,j,i,IM2) += mcoupv[index*32+16];
        mcoup(k,j,i+1,IM2) += mcoupv[index*32+17];
        mcoup(k,j,i+2,IM2) += mcoupv[index*32+18];
        mcoup(k,j+1,i,IM2) += mcoupv[index*32+19];
        mcoup(k,j+1,i+1,IM2) += mcoupv[index*32+20];
        mcoup(k,j+1,i+2,IM2) += mcoupv[index*32+21];
        mcoup(k,j+2,i,IM2) += mcoupv[index*32+22];
        mcoup(k,j+2,i+1,IM2) += mcoupv[index*32+23];

        mcoup(k,j,i,IM3) += mcoupv[index*32+24];
        mcoup(k,j,i+1,IM3) += mcoupv[index*32+25];
        mcoup(k,j,i+2,IM3) += mcoupv[index*32+26];
        mcoup(k,j+1,i,IM3) += mcoupv[index*32+27];
        mcoup(k,j+1,i+1,IM3) += mcoupv[index*32+28];
        mcoup(k,j+1,i+2,IM3) += mcoupv[index*32+29];
        mcoup(k,j+2,i,IM3) += mcoupv[index*32+30];
        mcoup(k,j+2,i+1,IM3) += mcoupv[index*32+31];
      }
    }
  }

  free(mcoupv);

 return;
}
