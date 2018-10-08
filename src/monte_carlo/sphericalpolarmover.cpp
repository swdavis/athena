//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sphericalpolarmover.cpp
//  \brief implementation for moving photons through sphericalpolar grid

// Athena++ headers
#include "photon.hpp"
#include "photonmover.hpp"
#include "../mesh/mesh.hpp"
#include "debug.hpp"
#define MAXITER 1000000
//#define DEBUG
//#define OUTTEST
#define NBUFFER 100
#define TAU_TOL_COH  5.

// Implementation of Sphericalpolar Photon mover

SphericalPolarMover::SphericalPolarMover(MonteCarloBlock *pmcb) 
  : PhotonMover(pmcb) {

}

SphericalPolarMover::~SphericalPolarMover() {

}

//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::Move(Photon *pphot)
//  \brief Moves photon along straight line specified number of mean free paths or until
//         photon leave monte carlo block

void SphericalPolarMover::Move(Photon *pphot) {

  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  MCCoord *pco = pmy_mcb->pcoord;

  // get number of mean free paths photon will travel
  Real TauRemaining = GetOpticalDepth(pran);

  // References for momentum vectors
  Real& kx = pphot->kcart[0];
  Real& ky = pphot->kcart[1];
  Real& kz = pphot->kcart[2];
  Real& kr  = pphot->k[0];
  Real& kth = pphot->k[1];
  Real& kph = pphot->k[2];
  bool thface = false;

#ifdef DEBUG
  typedef struct {
    Real dl, dlr, dlt, dlp;
    Real cth, sth, cph, sph;
    Real kr, kth, kph;
    Real kx, ky, kz;
    Real x,y,z;
    int i,j,k;
    bool ascend[3];
    Real l_ext,det,lm,lp;
  } debug_t;
  debug_t db[NBUFFER];
#endif
#ifdef OUTTEST
  Real rf,thf,phf,dl0;
  FinalPositionSphericalPolar(pmcb,pco,pphot,rf,thf,phf,dl0);
#endif


  Real cth = cos(pphot->x[1]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[2]);
  Real sph = sin(pphot->x[2]);

  int iter = 0;
  // Move photon until requisite # of mean free paths or escape
  // MAXITER is present to account for (near) infinite trajectories in optically thin,
  // periodic domains.
  while( (TauRemaining > 0.) && (pphot->status == EVOLVING) && (iter < MAXITER)) {
    iter++;

    kr  = kx * sth * cph + ky * sth * sph + kz * cth;
    kth = kx * cth * cph + ky * cth * sph - kz * sth;
    kph = -kx * sph + ky * cph;
    // Compute cartesian positions
    Real r0 = pphot->x[0];
    Real x0 = r0 * sth * cph;
    Real y0 = r0 * sth * sph;
    Real z0 = r0 * cth;
    
    // Compute distance to all faces
    Real dl, dlr, dlt, dlp;
    bool ascend[3];

    // r face
    if (kr > 0.0) { // Can only intersect outer sphere
      Real ri = pco->x1f(pphot->i1+1) / r0;
      Real det = 1. + (SQR(ri) - 1.) / SQR(kr);
      dlr = kr * r0 * (sqrt(det) - 1.);
      ascend[0] = true;
    } else if (kr < 0.0) { // Can intersect either sphere
      Real ri = pco->x1f(pphot->i1) / r0;
      Real det = 1.0 + (SQR(ri)-1.)/SQR(kr);
      if (det < 0.) {
        // ray does not intersect inner sphere so must intersect outer
        ri = pco->x1f(pphot->i1+1) / r0;
        det = 1. + (SQR(ri) - 1.) / SQR(kr);
        dlr = -kr * r0 * (sqrt(det) + 1.); //positive solution
        ascend[0] = true;
      } else {
        // ray intersects inner sphere first
        // one or two solutions -- pick shorter dlr
        dlr = r0 * kr * (sqrt(det) - 1.);
        ascend[0] = false;   
      }
    } else { // kr == 0
      Real ri = pco->x1f(pphot->i1+1) / r0;
      if (ri > 1.) {
        dlr = r0 * sqrt(SQR(ri) - 1.);
	ascend[0] = true;
      } else {
	std::cout << "Warning: kr == 0 and ri < r0, absorbing photon" << std::endl;
        pphot->status = DESTROYED;
        return;
      }
    }
    
    // theta face
    Real thi;
    if (kth > 0) {
      thi = pco->x2f(pphot->i2+1);
      ascend[1] = true;
    } else if (kth < 0) {
      thi = pco->x2f(pphot->i2);
      ascend[1] = false;
    } else {
      if (pphot->x[1] < 0.5*PI) {
        thi = pco->x2f(pphot->i2+1);
	ascend[1] = true;
      } else {
        thi = pco->x2f(pphot->i2);
	ascend[1] = false;
      }
    }
    Real cthi = cos(thi);

    // impact parameters
    Real b2 = SQR(r0) * (1. - SQR(kr));
    Real bz = r0 * (cth - kr * kz);
    Real th_ninf = acos(-kz);
    Real th_pinf = acos(kz);
    Real l_ext = (kz * b2 - r0 * kr * bz) / bz;
    Real cth_ext2 = SQR(kz) + SQR(bz)/b2;
    Real th_ext;
    if (bz < 0) {
      th_ext = acos(-sqrt(cth_ext2));
    } else {
      th_ext = acos(sqrt(cth_ext2));
    }

    // check for turning point
    if ((l_ext > 0) && ((thi-th_ninf) * (thi-th_ext)) > 0.) {
      dlt = HUGE_NUMBER;
    } else if ((l_ext <= 0) && ((thi-th_ext) * (thi-th_pinf)) > 0.) {
      dlt = HUGE_NUMBER;
    } else { //No turning point
      Real det = b2*SQR(cthi)*(cth_ext2-SQR(cthi)); //det should be positive
      Real lm = (-bz*kz - sqrt(det)) / (SQR(kz)-SQR(cthi)) - r0*kr;
      Real lp = (-bz*kz + sqrt(det)) / (SQR(kz)-SQR(cthi)) - r0*kr;
      if (lm <= 0) {
        dlt = lp;
      } else if (lp <= 0) {
        dlt = lm;
      } else {
        if (lp >= lm) 
          dlt = lm;
        else
          dlt = lp;
      }
    }

    // phi face
    Real phii;
    if (kph > 0) {
      phii = pco->x3f(pphot->i3+1);
      ascend[2] = true;
    } else if (kph < 0.0) {
      phii = pco->x3f(pphot->i3);
      ascend[2] = false;
    }
    Real tphii;
    if (kph != 0.) {
      tphii = tan(phii);
      if (fabs(kx * tphii - ky) < TINY_NUMBER) {
        dlp = HUGE_NUMBER;
      } else {
        dlp = r0 * sth * (sph-tphii*cph) / (kx*tphii-ky);
        if (dlp < 0.) dlp = HUGE_NUMBER;
      }
    } else
      dlp = HUGE_NUMBER;
    //if (fabs(dlp) < 1.e-8)
    //if (dlp == 0.0)
    // printf("dlp=0: %g %g %g %g %g %g %g %d\n",phii,tphii,kph,(sph-tphii*cph),(kx*tphii-ky),sth,r0,pphot->i2);

    if ((dlr <= dlt) && (dlr <= dlp)) {
      dl = dlr;
    } else if ((dlp <= dlt) && (dlp <= dlr)) {
      dl = dlp;
    } else if ((dlt <= dlr) && (dlt <= dlp)) {
      dl = dlt;
    }

    //checking turning point for theta
    if ((dl > l_ext) && (l_ext > 0.)) {
      if (kth < 0.) {
        cthi = cos(pco->x2f(pphot->i2+1));
      } else if (kth > 0.) {
        cthi = cos(pco->x2f(pphot->i2));
      }
      b2 = SQR(r0) * (1. - SQR(kr));
      bz = r0 * (cth - kr * kz);
      Real det = b2 * SQR(cthi) * (cth_ext2 - SQR(cthi)); //det should be positive
      Real lm = (-bz*kz - sqrt(det)) / (SQR(kz) - SQR(cthi)) - kr*r0;
      Real lp = (-bz*kz + sqrt(det)) / (SQR(kz) - SQR(cthi)) - kr*r0; 
#ifdef DEBUG
      if (iter < NBUFFER) {
        db[iter-1].lm = lm;
        db[iter-1].lp = lp;
        db[iter-1].det = det;
      }
#endif
      if (lm <= 0) {
        if ( (lp > 0) && !(thface) )
          dlt = lp;
        else
          dlt = HUGE_NUMBER;
      } else {
        if (lp > 0)
	  dlt = (lm > lp) ? lm : lp;
        else
          dlt = lm;
      }
      if (dlt <= dl) {
	ascend[1] = !(ascend[1]);
      }
    }

    int face;
    NextFace(dlr,dlt,dlp,face,dl);
#ifdef DEBUG
    if (iter < NBUFFER) {
      db[iter-1].dl = dl; db[iter-1].dlr = dlr; db[iter-1].dlt = dlt; db[iter-1].dlp = dlp;
      db[iter-1].cth = cth; db[iter-1].sth = sth; db[iter-1].cph = cph; db[iter-1].sph = sph;
      db[iter-1].kr = kr; db[iter-1].kth = kth; db[iter-1].kph = kph;
      db[iter-1].kx = kx; db[iter-1].ky = ky; db[iter-1].kz = kz;
      db[iter-1].x = pphot->x[0]; db[iter-1].y = pphot->x[1]; db[iter-1].z = pphot->x[2];
      db[iter-1].i = pphot->i1; db[iter-1].j = pphot->i2; db[iter-1].k = pphot->i3;
      db[iter-1].ascend[0] = ascend[0]; db[iter-1].ascend[1] = ascend[1]; db[iter-1].ascend[2] = ascend[2];
      db[iter-1].l_ext = l_ext;
    }
#endif
    

    Real chi = pphot->sct_coef + pphot->abs_coef;
    chi = (chi > TINY_NUMBER) ? chi : TINY_NUMBER; // return max
    if (dl > TauRemaining / chi) { // Photon remains in zone
      if ((acceleration) && ((pphot->abs_coef+pphot->sct_coef) * dl > TAU_TOL_COH)) {
        printf("a ");
        // position packet on sphere of radius dl
        Real mu = 2.*pran->uniform()-1.0;
        Real stheta = sqrt(1.0-mu*mu);
        Real phi = 2.*PI*pran->uniform();
        
        
        x0 += stheta*cos(phi) * dl;
        y0 += stheta*sin(phi) * dl;
        z0 += mu * dl;
        pphot->x[0] = sqrt(SQR(x0)+SQR(y0)+SQR(z0)); 
        pphot->x[1] = acos((z0) / pphot->x[0]);
        pphot->x[2] = atan2(y0,x0);
        
        // current assume photon k parallel to rhat
        kr = 1.0;
        kth = 0.;
        kph = 0.;
        
        // draw from path length distribution and reduce weight accordingly
        Real mrw = MRWDist(pran);          
        while (mrw <= 0.)
          mrw = MRWDist(pran);          
        Real ct = -log(mrw)*SQR(dl)/SQR(PI)*(pphot->abs_coef+pphot->sct_coef)*3.;
        pphot->weight *= exp(-ct*pphot->abs_coef);
      } else {
        dl = TauRemaining/chi;
        // Update moments
        if (pmcb->moments_flag)
          pmcb->UpdateMoments(pphot,dl);
        // Update postions
        pphot->x[0] = sqrt(SQR(r0) + 2. * dl * kr * r0 + SQR(dl));
        pphot->x[1] = acos((z0 + kz * dl) / pphot->x[0]);
        pphot->x[2] = atan2(y0 + ky * dl,x0 + kx * dl);
        
        if (pphot->x[2] < 0.)
          pphot->x[2] += 2.*PI;
        return;
      }
    } else { // Photon moves to next zone and reduce TauRemaining
      // Update moments
      if (pmcb->moments_flag)
	pmcb->UpdateMoments(pphot,dl);
      // Update postions
      pphot->x[0] = sqrt(SQR(r0) + 2. * dl * kr * r0 + SQR(dl));
      pphot->x[1] = acos((z0 + kz * dl) / pphot->x[0]);
      pphot->x[2] = atan2(y0 + ky * dl,x0 + kx * dl);
      if (pphot->x[2] < 0.)
	pphot->x[2] += 2.*PI;
    
      TauRemaining -= chi * dl;
      // move photon to next zone and pdate angular positions
      MovePhotonToNextZone(pphot,pco,pmcb,face,ascend);
      if ((face == 1) || (face == 3) || (face == 4) || (face == 6))
        thface = true;
      cth = cos(pphot->x[1]);
      sth = sqrt(1. - SQR(cth));
      cph = cos(pphot->x[2]);
      sph = sin(pphot->x[2]);
      
    }
  }
  if (iter >= MAXITER) {
    std::cout << "Warning: iter exceeded ITERMAX in photon mover." << std::endl;
#ifdef DEBUG
    int nmax = (NBUFFER > iter) ? iter : NBUFFER;
    for (int i=0; i < nmax; ++i) {
      // printf("dl: %16.12e %16.12e %16.12e %16.12e\n");
      printf("--------------------------------\n %d\n",i);
      printf("dl: %16.12e %16.12e %16.12e %16.12e\n",db[i].dl,db[i].dlr,db[i].dlt,db[i].dlp);
      printf("ang: %16.12e %16.12e %16.12e %16.12e\n", db[i].cth,db[i].sth,db[i].cph,db[i].sph);
      printf("k: %16.12e %16.12e %16.12e\n",db[i].kr,db[i].kth,db[i].kph);
      printf("kc: %16.12e %16.12e %16.12e\n",db[i].kx,db[i].ky,db[i].kz);
      printf("x: %16.12e %16.12e %16.12e\n",db[i].x,db[i].y,db[i].z);
      printf("i: %d %d %d\n",db[i].i,db[i].j,db[i].k);
      printf("ascend: %d %d %d\n",db[i].ascend[0],db[i].ascend[1],db[i].ascend[2]);
      printf("l_ext: %g %g %g %g\n",db[i].l_ext,db[i].lm,db[i].lp,db[i].det);
      printf("xf: ");
      if ((db[i].i >= pmcb->is) && (db[i].i <= pmcb->ie))
	printf("%16.12e %16.12e ",pco->x1f(db[i].i),pco->x1f(db[i].i+1));
      if ((db[i].j >= pmcb->js) && (db[i].j <= pmcb->je))
	printf("%16.12e %16.12e ",pco->x2f(db[i].j),pco->x2f(db[i].j+1));
      if ((db[i].k >= pmcb->ks) && (db[i].k <= pmcb->ke))
	printf("%16.12e %16.12e ",pco->x3f(db[i].k),pco->x3f(db[i].k+1));
      printf("\n");
    }
#endif
    pphot->status = DESTROYED;
  }

#ifdef OUTTEST
  Real xf = rf*sin(thf)*cos(phf);
  Real yf = rf*sin(thf)*sin(phf);
  Real zf = rf*cos(thf);
  Real xp =  pphot->x[0]*sin(pphot->x[1])*cos(pphot->x[2]);
  Real yp =  pphot->x[0]*sin(pphot->x[1])*sin(pphot->x[2]);
  Real zp =  pphot->x[0]*cos(pphot->x[1]);
  Real delta = sqrt(SQR(xf-xp)+SQR(yf-yp)+SQR(zf-zp));
  Real dmax = 1.e-8*rf;
  //std::cout << delta << " " << xf << " " << xp << std::endl;
  if ((delta > dmax)&&(iter < MAXITER)) {
    std::cout << "-----------------------" << std::endl;
    std::cout << delta <<  ' ' << iter << std::endl;
    std::cout << "k: " << pphot->k[0] << ' ' << pphot->k[1] << ' ' << pphot->k[2] << std::endl;
    std::cout << "xf: " << rf << ' ' << thf << ' ' << phf << ' ' << dl0 << std::endl;
    std::cout << "xp: " << pphot->x[0] << ' ' <<  pphot->x[1] << ' ' <<  pphot->x[2] << std::endl;
  }
#endif

}

//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::CartesianToCurvalinear(Photon *pphot)
//  \brief convert k vector from cartesian to curvalinear

void SphericalPolarMover::CartesianToCurvalinear(Photon *pphot) {

  Real cth = cos(pphot->x[1]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[2]);
  Real sph = sin(pphot->x[2]);
  // Compute spherical-polar
  pphot->k[0] = pphot->kcart[0]*sth*cph + pphot->kcart[1]*sth*sph + pphot->kcart[2]*cth;
  pphot->k[1] = pphot->kcart[0]*cth*cph + pphot->kcart[1]*cth*sph - pphot->kcart[2]*sth;
  pphot->k[2] = -pphot->kcart[0]*sph + pphot->kcart[1]*cph;
  
}


//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::CurvalinearToCartesian(Photon *pphot)
//  \brief convert k vector from curvalinear to cartesian

void SphericalPolarMover::CurvalinearToCartesian(Photon *pphot) {

  Real cth = cos(pphot->x[1]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[2]);
  Real sph = sin(pphot->x[2]);
  // Compute cartesian
  pphot->kcart[0] = pphot->k[0]*sth*cph + pphot->k[1]*cth*cph - pphot->k[2]*sph;
  pphot->kcart[1] = pphot->k[0]*sth*sph + pphot->k[1]*cth*sph + pphot->k[2]*cph;
  pphot->kcart[2] = pphot->k[0]*cph - pphot->k[1]*sph;
  
}

