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


//#define DEBUG_SM
//#define NBUFFER 50

//----------------------------------------------------------------------------------------
//! SphericalPolarMover class constructor, derived from PhotonMover base class

SphericalPolarMover::SphericalPolarMover(MonteCarloBlock *pmcb)
  : PhotonMover(pmcb) {

}

//----------------------------------------------------------------------------------------
//! destructor

SphericalPolarMover::~SphericalPolarMover() {

}

//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::Move(Photon *pphot, int ips, int ipe)
//! \brief Moves photon using cell-by-cell approach through spherical polar grid

void SphericalPolarMover::Move(Photon *pphot, int ips, int ipe) {

  MonteCarloBlock *pmcb = pmy_mcb;
  MCRandom *pran = pmy_mcb->pran;
  MCCoord *pco = pmy_mcb->pcoord;
  PhotonTrajectoryList *ptraj = pmy_mcb->ptraj;

  for (int ip=ips; ip<=ipe; ip++) {
    // get number of mean free paths photon will travel
    Real tauremaining = GetOpticalDepth(pran);

    // References for momentum vectors
    //CurvalinearToCartesian(pphot);// SWD: Redundant calculation of cth,sth,cph,sph
    //Real& kx = pphot->kcart[0];
    //Real& ky = pphot->kcart[1];
    //Real& kz = pphot->kcart[2];
    Real& kr  = pphot->k1p[ip];
    Real& kth = pphot->k2p[ip];
    Real& kph = pphot->k3p[ip];
    bool thface = false;

#ifdef DEBUG_SM
    typedef struct {
      Real dlr, dlt, dlp;
      Real cth, sth, cph, sph;
      Real x,y,z;
      int i,j,k;
      bool ascend[3];
      Real l_ext,det,lm,lp;
    } debug_t;
    debug_t db[NBUFFER];
#endif

    Real cth = cos(pphot->x2p[ip]);
    Real sth = sin(pphot->x2p[ip]);
    Real cph = cos(pphot->x3p[ip]);
    Real sph = sin(pphot->x3p[ip]);

    // Make sure kcart is set
    Real kx = kr*sth*cph + kth*cth*cph - kph*sph;
    Real ky = kr*sth*sph + kth*cth*sph + kph*cph;
    Real kz = kr*cth - kth*sth;

    int iter = 0;
    Real c_cgs = 2.99792458e10;

    // Move photon until requisite # of mean free paths or escape
    while( (tauremaining > 0.) && (pphot->statp[ip] == EVOLVING) && (iter < checkmove) &&
           (pphot->dtp[ip] > 0.) ) {
      iter++;

      // Compute cartesian positions
      Real r0 = pphot->x1p[ip];
      Real x0 = r0 * sth * cph;
      Real y0 = r0 * sth * sph;
      Real z0 = r0 * cth;

      // Compute distance to all faces
      Real dlr, dlt, dlp;
      bool ascend[3];

      // r face
      if (kr > 0.0) { // Can only intersect outer sphere
        Real ri = pco->x1f(pphot->i1p[ip]+1) / r0;
        Real det = 1. + (SQR(ri) - 1.) / SQR(kr);
        dlr = kr * r0 * (sqrt(det) - 1.);
        ascend[0] = true;
      } else if (kr < 0.0) { // Can intersect either sphere
        Real ri = pco->x1f(pphot->i1p[ip]) / r0;
        Real det = 1.0 + (SQR(ri)-1.)/SQR(kr);
        if (det < 0.) {
          // ray does not intersect inner sphere so must intersect outer
          ri = pco->x1f(pphot->i1p[ip]+1) / r0;
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
        Real ri = pco->x1f(pphot->i1p[ip]+1) / r0;
        if (ri > 1.) {
          dlr = r0 * sqrt(SQR(ri) - 1.);
          ascend[0] = true;
        } else {
          std::cout << "Warning: kr == 0 and ri < r0, absorbing photon" << std::endl;
          pphot->PrintPhoton(ip);
          pphot->statp[ip] = DESTROYED;
          break; // break out of while loop
        }
      }

      // theta face
      Real thi;
      if (kth > 0) {
        thi = pco->x2f(pphot->i2p[ip]+1);
        ascend[1] = true;
      } else if (kth < 0) {
        thi = pco->x2f(pphot->i2p[ip]);
        ascend[1] = false;
      } else {
        if (pphot->x2p[ip] < 0.5*PI) {
          thi = pco->x2f(pphot->i2p[ip]+1);
          ascend[1] = true;
        } else {
          thi = pco->x2f(pphot->i2p[ip]);
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
        phii = pco->x3f(pphot->i3p[ip]+1);
        ascend[2] = true;
      } else if (kph < 0.0) {
        phii = pco->x3f(pphot->i3p[ip]);
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
      } else {
        dlp = HUGE_NUMBER;
      }


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
          cthi = cos(pco->x2f(pphot->i2p[ip]+1));
        } else if (kth > 0.) {
          cthi = cos(pco->x2f(pphot->i2p[ip]));
        }
        b2 = SQR(r0) * (1. - SQR(kr));
        bz = r0 * (cth - kr * kz);
        Real det = b2 * SQR(cthi) * (cth_ext2 - SQR(cthi)); //det should be positive
        Real lm = (-bz*kz - sqrt(det)) / (SQR(kz) - SQR(cthi)) - kr*r0;
        Real lp = (-bz*kz + sqrt(det)) / (SQR(kz) - SQR(cthi)) - kr*r0;
#ifdef DEBUG_SM
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
#ifdef DEBUG_SM
      if (iter < NBUFFER) {
        db[iter-1].dl = dl; db[iter-1].dlr = dlr;
        db[iter-1].dlt = dlt; db[iter-1].dlp = dlp;
        db[iter-1].cth = cth; db[iter-1].sth = sth;
        db[iter-1].cph = cph; db[iter-1].sph = sph;
        db[iter-1].kr = kr; db[iter-1].kth = kth; db[iter-1].kph = kph;
        db[iter-1].kx = kx; db[iter-1].ky = ky; db[iter-1].kz = kz;
        db[iter-1].x = pphot->x[IMC1];
        db[iter-1].y = pphot->x[IMC2]; db[iter-1].z = pphot->x[IMC3];
        db[iter-1].i = pphot->i1;
        db[iter-1].j = pphot->i2; db[iter-1].k = pphot->i3;
        db[iter-1].ascend[0] = ascend[0]; db[iter-1].ascend[1] = ascend[1];
        db[iter-1].ascend[2] = ascend[2];
        db[iter-1].l_ext = l_ext;
      }
#endif

      Real chi = GetExtinctionCoefficient(pphot->acp[ip],pphot->scp[ip]);
      chi = (chi > TINY_NUMBER) ? chi : TINY_NUMBER; // return max
      bool test = false;
      if (dl > tauremaining / chi) { // Photon remains in zone
        bool accel_success = false;
        if (acceleration) {
          Real dist = pco->dmin(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip]);
          // Try/perform MRW acceleration if optical depth is large enough
          if (pmcb->coherent_scattering) {
            Real tauacc = 10.;
            if ((pphot->acp[ip]+pphot->scp[ip]) * dist > tauacc)
              accel_success = MRWAcceleration(pphot,pran,dist,tauacc,ip);
          } else {
            Real tauacc = 10.;
            Real pio = pmcb->planck_inv_opacity(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip]);
            if (pio * dist > tauacc)
              accel_success = MRWAcceleration(pphot,pran,dist,tauacc,ip);
          }
        }

        if (!accel_success) {
          if (pphot->statp[ip] != EVOLVING)
            break; // break out of while loop
          // compute distance remaining in zone
          dl = tauremaining/chi;
          pphot->dtp[ip] -= dl/c_cgs;
          // Update moments
          if (pmcb->call_moments)
            pmcb->UpdateMoments(pphot,dl,1.,ip);
          // Update postions
          pphot->x1p[ip] = sqrt(SQR(r0) + 2. * dl * kr * r0 + SQR(dl));
          pphot->x2p[ip] = acos((z0 + kz * dl) / pphot->x1p[ip]);
          pphot->x3p[ip] = atan2(y0 + ky * dl,x0 + kx * dl);
          if (pphot->x3p[ip] < 0.)
            pphot->x3p[ip] += 2.*PI;
          pphot->x0p[ip] += pphot->k0p[ip] * dl;

          // Update k vector
          cth = cos(pphot->x2p[ip]);
          sth = sin(pphot->x2p[ip]);
          cph = cos(pphot->x3p[ip]);
          sph = sin(pphot->x3p[ip]);
          kr  = kx * sth * cph + ky * sth * sph + kz * cth;
          kth = kx * cth * cph + ky * cth * sph - kz * sth;
          kph = -kx * sph + ky * cph;
        }
        if (ptraj != NULL) ptraj->AddToTrajectory(pphot,ip);
        tauremaining = 0.; // will cause break out while loop

      } else { // Photon moves to next zone and reduce tauremaining
        // Update moments
        pphot->dtp[ip] -= dl/c_cgs;

        if (pmcb->call_moments)
          pmcb->UpdateMoments(pphot,dl,1.,ip);
        // Update positions
        pphot->x1p[ip] = sqrt(SQR(r0) + 2. * dl * kr * r0 + SQR(dl));
        pphot->x2p[ip] = acos((z0 + kz * dl) / pphot->x1p[ip]);
        pphot->x3p[ip] = atan2(y0 + ky * dl,x0 + kx * dl);
        if (pphot->x3p[ip] < 0.)
          pphot->x3p[ip] += 2.*PI;
        pphot->x0p[ip] += pphot->k0p[ip] * dl;

        tauremaining -= chi * dl;
        // move photon to next zone and update angular positions
        MovePhotonToNextZone(pphot,pco,pmcb,face,ascend,ip);
        if ((face == 1) || (face == 3) || (face == 4) || (face == 6))
          thface = true;
        cth = cos(pphot->x2p[ip]);
        sth = sin(pphot->x2p[ip]);
        cph = cos(pphot->x3p[ip]);
        sph = sin(pphot->x3p[ip]);
        kr  = kx * sth * cph + ky * sth * sph + kz * cth;
        kth = kx * cth * cph + ky * cth * sph - kz * sth;
        kph = -kx * sph + ky * cph;
        if (ptraj != NULL) ptraj->AddToTrajectory(pphot,ip);
      }

    }

    // -------------------------- Debugging ----------------------------------------------
    if (iter >= checkmove) {
      std::stringstream msg;
      msg << "Warning: iter exceeded " << checkmove << " in photon mover.";
      pphot->PrintPhoton(msg.str(),ip);
#ifdef DEBUG_SM
      int nmax = (NBUFFER > iter) ? iter : NBUFFER;
      for (int i=0; i < nmax; ++i) {
        // printf("dl: %16.12e %16.12e %16.12e %16.12e\n");
        printf("--------------------------------\n %d\n",i);
        printf("dl: %16.12e %16.12e %16.12e %16.12e\n",db[i].dl,db[i].dlr,db[i].dlt,
               db[i].dlp);
        printf("ang: %16.12e %16.12e %16.12e %16.12e\n", db[i].cth,db[i].sth,db[i].cph,
               db[i].sph);
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
      pphot->statp[ip] = DESTROYED;
    } // iter >= checkmove

  } // loop over ip

}

// SWD: Deprecated and slated for removal
//----------------------------------------------------------------------------------------
//! \fn void SphericalPolarMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4])
//! \brief convert k vector from curvalinear to cartesian

void SphericalPolarMover::CurvalinearToCartesian(Photon *pphot, Real kcart[4]) {

  // SWD: Broken by parallelization
  /*
  Real cth = cos(pphot->x[IMC2]);
  Real sth = sqrt(1. - SQR(cth));
  Real cph = cos(pphot->x[IMC3]);
  Real sph = sin(pphot->x[IMC3]);
  // Compute cartesian
  kcart[IMC1] = pphot->k[IMC1]*sth*cph + pphot->k[IMC2]*cth*cph - pphot->k[IMC3]*sph;
  kcart[IMC2] = pphot->k[IMC1]*sth*sph + pphot->k[IMC2]*cth*sph + pphot->k[IMC3]*cph;
  kcart[IMC3] = pphot->k[IMC1]*cth - pphot->k[IMC2]*sth;
  */
}
