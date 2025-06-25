//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file from_array.cpp
//! \brief Problem generator for initializing with preexisting array from HDF5 input

// C headers

// C++ headers
#include <algorithm>  // max()
#include <string>     // c_str(), string

// Athena++ headers
#include "../athena.hpp"              // Real
#include "../athena_arrays.hpp"       // AthenaArray
#include "../field/field.hpp"         // Field
#include "../globals.hpp"             // Globals
#include "../hydro/hydro.hpp"         // Hydro
#include "../eos/eos.hpp"                  // EquationOfState
#include "../inputs/hdf5_reader.hpp"  // HDF5ReadRealArray()
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"     // ParameterInput
#include "../monte_carlo/montecarlo.hpp"
#include "../monte_carlo/photon.hpp"
#include "../monte_carlo/mcutils.hpp"

namespace {
  // Global variables
  bool tnorm;
  Real logemin, logemax;
  std::string emission_type;
  // frequency table parameters
  int nfre, nrho, ntem;
  Real lmine, lmaxe, dle, lmint, lmaxt, dlt, lmind, lmaxd, dld;
  AthenaArray<Real> fre_grid;
  AthenaArray<Real> temp_grid;
  AthenaArray<Real> rho_grid;
  AthenaArray<Real> ross_tab;
  AthenaArray<Real> plan_tab;
  AthenaArray<Real> eta_cum_tab;
  AthenaArray<Real> eta_tab;

  //functions
  Real TableOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip);
  Real IntegrateEmission(Real temp, Real num, Real nup, Real am, Real ap);
  Real Planck(Real temp, Real nu);
  Real TableEmission(MonteCarloBlock *pmcb, int k, int j, int i);
  Real SampleEmissivity(MonteCarloBlock *pmcb, Photon *pphot, int ip);
}

std::vector<float> x1coord;
std::vector<float> x2coord;
std::vector<float> x3coord;

int getindex(std::vector<float> vec, float val){
  std::vector<float>::iterator it = std::find(vec.begin(), vec.end(), val);
  int index = std::distance(vec.begin(), it);
  return index;
}

void MonteCarlo::InitUserMonteCarloData(ParameterInput *pin) {

  nuser_var = 3;
  emission_type = pin->GetOrAddString("montecarlo","emission","none");
  if (emission_type == "freefree")
    return;

  // Read in opacity table
  FILE  *opac_file;
  if ( (opac_file=fopen("out_opacity_table_nfreq32.txt","r"))==NULL) {
    std::stringstream msg;
    msg << "FATAL ERROR: Could not open out_opacity_table_nfreq16.txt." << std::endl;
    ATHENA_ERROR(msg);
  }

  fscanf(opac_file,"%d",&(nfre));
  fscanf(opac_file,"%d",&(ntem));
  fscanf(opac_file,"%d",&(nrho));

  // Create arrays for opacity
  fre_grid.NewAthenaArray(nfre+1);
  temp_grid.NewAthenaArray(ntem);
  rho_grid.NewAthenaArray(nrho);
  ross_tab.NewAthenaArray(nfre,ntem,nrho);
  plan_tab.NewAthenaArray(nfre,ntem,nrho);

  for(int i=1; i<=nfre; ++i){
    fscanf(opac_file,"%lf",&(fre_grid(i)));
  }
  fre_grid(0) = fre_grid(1)*fre_grid(1)/fre_grid(2);
  // convert to erg
  Real keverg = 1.602176634e-9;
  for(int i=0; i<=nfre; ++i)
    fre_grid(i) *= keverg;
  lmine = std::log10(fre_grid(1));
  lmaxe = std::log10(fre_grid(nfre-1));
  dle = (lmaxe-lmine)/static_cast<Real>(nfre-2);
  // temperature grid (keV)
  for(int i=0; i<ntem; ++i){
    fscanf(opac_file,"%lf",&(temp_grid(i)));
  }
  // convert to kelvin
  Real kb = 1.380649e-16;
  for(int i=0; i<ntem; ++i)
    temp_grid(i) *= keverg/kb;
  lmint = std::log10(temp_grid(0));
  lmaxt = std::log10(temp_grid(ntem-1));
  dlt = (lmaxt-lmint)/static_cast<Real>(ntem-1);
  // note temperature grid not uniform in log

  // density grid (g/cm^3)
  for(int i=0; i<nrho; ++i) {
    fscanf(opac_file,"%lf",&(rho_grid(i)));
  }
  lmind = std::log10(rho_grid(0));
  lmaxd = std::log10(rho_grid(nrho-1));
  dld = (lmaxd-lmind)/static_cast<Real>(nrho-1);

  // frequency integrated rosseland mean
  // Read in but not used
  Real buf;
  for(int j=0; j<ntem; ++j) {
    for(int i=0; i<nrho; ++i) {
      fscanf(opac_file,"%lf",&buf);
    }
  }

  // frequency integrated planck mean
  // Read in but not used
  for(int j=0; j<ntem; ++j) {
    for(int i=0; i<nrho; ++i) {
      fscanf(opac_file,"%lf",&buf);
    }
  }

  // ross mean for each frequency group
  for(int k=0; k<nfre; ++k) {
    for(int j=0; j<ntem; ++j) {
      for(int i=0; i<nrho; ++i) {
        fscanf(opac_file,"%lf",&(ross_tab(k,j,i)));
        ross_tab(k,j,i) *= rho_grid(i);
      }
    }
  }

  // planck mean for each frequency group
  for(int k=0; k<nfre; ++k) {
    for(int j=0; j<ntem; ++j) {
      for(int i=0; i<nrho; ++i) {
        fscanf(opac_file,"%lf",&(plan_tab(k,j,i)));
        plan_tab(k,j,i) *= rho_grid(i);
      }
    }
  }

  // Replaces plan_tab with free-free values
  /*Real dummy;
  for(int k=0; k<nfre; ++k) {
    Real ffnrm = 3.692146e8;
    Real heabund = 0.09; //hardcode for now (should be parameter)
    Real mp = 1.67262192369e-24;
    Real h = 6.62607015e-27;
    Real kb = 1.380649e-16;
    Real nu = fre_grid(k) / h;
    for(int j=0; j<ntem; ++j) {
      Real tgas = temp_grid(j);
      Real ehnu = exp(-h*nu / (kb * tgas) );
      for(int i=0; i<nrho; ++i) {
        Real nh = rho_grid(i) / (mp*(1.+4.*heabund));
        Real nhe = nh*heabund;
        Real ne = nh + 2.*nhe;
        fscanf(opac_file,"%lf",&(dummy));
        Real aff = ffnrm/sqrt(tgas)/pow(nu,3);
        Real opac = ne * (nh + 4. * nhe) * aff * (1. - ehnu);
        plan_tab(k,j,i) = opac;
      }
    }
    }*/

  // planck mean for each frequency group
  for(int j=0; j<ntem; ++j) {
    for(int i=0; i<nrho; ++i) {
      Real min = 1.e40;
      Real max = 1.e-40;
      for(int k=0; k<nfre; ++k) {
        min = (min > plan_tab(k,j,i)) ? plan_tab(k,j,i) : min;
        max = (max < plan_tab(k,j,i)) ? plan_tab(k,j,i) : max;
      }
      if (max/min < 1.2) {
        Real ffnrm = 3.692146e8;
        Real heabund = 0.09; //hardcode for now (should be parameter)
        Real mp = 1.67262192369e-24;
        Real h = 6.62607015e-27;
        Real kb = 1.380649e-16;
        Real tgas = temp_grid(j);
        Real nh = rho_grid(i) / (mp*(1.+4.*heabund));
        Real nhe = nh*heabund;
        Real ne = nh + 2.*nhe;
        for(int k=0; k<nfre; ++k) {
          Real nu = fre_grid(k) / h;
          Real ehnu = exp(-h*nu / (kb * tgas) );
          Real aff = ffnrm/sqrt(tgas)/pow(nu,3);
          Real opac = ne * (nh + 4. * nhe) * aff * (1. - ehnu);
          //printf("%d %g %g %g %g\n",k,temp_grid(i),rho_grid(j),plan_tab(k,j,i),opac);
          plan_tab(k,j,i) = opac;
          // plan_tab(k,j,i) = 1.e-60;

        }

      }
    }
  }

  fclose(opac_file);

  EnrollUserEmissionFunction(TableEmission);
  EnrollUserOpacityFunction(TableOpacity,true);
}

void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin) {

  if (emission_type == "freefree") {
    // Set the energy boundaries for free-free emission
    tnorm = pin->GetOrAddBoolean("problem","tnorm",false);
    if (tnorm) {
      // interpret as xmin/xmax with x=E/(kb*T)
      const Real kb = 1.380649e-16;
      logemin = log(kb*pin->GetReal("problem", "emin"));
      logemax = log(kb*pin->GetReal("problem", "emax"));
    } else {
      // interpret as emin/emax in eV
      const Real everg = 1.6021772e-12;
      logemin = log(everg*pin->GetReal("problem", "emin"));
      logemax = log(everg*pin->GetReal("problem", "emax"));
    }
  } else {
    // Create array for emissivity
    AthenaArray<Real> eta_nu_tab;
    eta_nu_tab.NewAthenaArray(nfre,ntem,nrho);
    Real h = 6.62607015e-27;
    for(int k=0; k<nfre; ++k) {
      Real nup = fre_grid(k+1)/h;
      Real num = fre_grid(k)/h;
      if (num > nup)
        printf("%d %g %g\n",k,fre_grid(k),fre_grid(k+1));
      for(int j=0; j<ntem; ++j) {
        Real temp = temp_grid(j);
        for(int i=0; i<nrho; ++i) {
          Real opp = plan_tab(k,j,i);
          Real opm;
          if (k == 0)
            opm = 0.;
          else
            opm = plan_tab(k-1,j,i);
          eta_nu_tab(k,j,i) = IntegrateEmission(temp,num,nup,opm,opp);
        }
      }
    }

    eta_tab.NewAthenaArray(ntem,nrho);
    for(int j=0; j<ntem; ++j) {
      Real temp = temp_grid(j);
      for(int i=0; i<nrho; ++i) {
        eta_tab(j,i) = 0.;
        for(int k=0; k<nfre; ++k) {
          eta_tab(j,i) += 4.*PI*eta_nu_tab(k,j,i);
          //printf("%d %d %d %g %g\n",j,i,k,eta_nu_tab(k,j,i),eta_tab(j,i));
        }
      }
    }


    int ncells1 = nx1 + 2*(NGHOST);
    int ncells2 = 1, ncells3 = 1;
    if (nx2 > 1) ncells2 = nx2 + 2*(NGHOST);
    if (nx3 > 1) ncells3 = nx3 + 2*(NGHOST);
    eta_cum_tab.NewAthenaArray(ncells3,ncells2,ncells1,nfre+1);
    Real keverg = 1.602176634e-9;
    for(int k=ks; k<=ke; ++k) {
      for(int j=js; j<=je; ++j) {
        for(int i=is; i<=ie; ++i) {
          Real ld = log10(rho(k,j,i));
          ld = (ld < lmind) ? lmind : ld;
          ld = (ld > lmaxd) ? lmaxd : ld;
          Real temp = tgas(k,j,i);
          Real lt = log10(temp);
          //printf("%d %d %d %g %g %g\n",k,j,i,lt,tgas(k,j,i),rho(k,j,i));
          lt = (lt < lmint) ? lmint : lt;
          lt = (lt > lmaxt) ? lmaxt : lt;
          Real xi = (ld - lmind) / dld;
          int ii = std::floor(xi);
          xi -= static_cast<Real>(ii);
          Real xj = (lt - lmint) / dlt;
          int jj = std::floor(xj);
          while ((jj<ntem-2) && (temp_grid(jj+1) < temp)){
            jj++;
          }
          while ((jj>0) && (temp_grid(jj) > temp)){
            jj--;
          }
          if(jj > ntem-2) {
            jj = ntem-2;
          }
          xj = (temp-temp_grid(jj))/(temp_grid(jj+1)-temp_grid(jj));

          eta_cum_tab(k,j,i,0) = 0.;
          for(int l=0; l<nfre; ++l) {
            Real eta = (1.-xi) * ((1.-xj)*eta_nu_tab(l,jj,ii)+xj*eta_nu_tab(l,jj+1,ii))
              + xi * ((1.-xj)*eta_nu_tab(l,jj,ii+1)+xj*eta_nu_tab(l,jj+1,ii+1));
            eta_cum_tab(k,j,i,l+1) =  eta_cum_tab(k,j,i,l)+eta;
          }

          for(int l=0; l<nfre+1; ++l) {
            eta_cum_tab(k,j,i,l) /= eta_cum_tab(k,j,i,nfre);
          //if ((Globals::my_rank == 0) && (k==ks) && (j == js) && (i == is)) {
          //  printf("%d %g %g\n",l,fre_grid(l)/keverg,eta_cum_tab(k,j,i,l));
          //}
          }
        }
      }
    }
    eta_nu_tab.DeleteAthenaArray();
  }

}


void Mesh::InitUserMeshData(ParameterInput *pin) {

  // Read in hdf5 file to initialize pgen
  std::string input_filename = pin->GetString("problem", "input_filename");

  int mesh_nx1 = pin->GetInteger("mesh", "nx1");
  Real mesh_x1min = pin->GetReal("mesh", "x1min");
  Real mesh_x1max = pin->GetReal("mesh", "x1max");
  int mesh_nx2 = pin->GetInteger("mesh", "nx2");
  Real mesh_x2min = pin->GetReal("mesh", "x2min");
  Real mesh_x2max = pin->GetReal("mesh", "x2max");
  int mesh_nx3 = pin->GetInteger("mesh", "nx3");
  Real mesh_x3min = pin->GetReal("mesh", "x3min");
  Real mesh_x3max = pin->GetReal("mesh", "x3max");
  Real x1ratio = pin->GetReal("mesh", "x1rat");

  //load data file
  int start_file[3] = {0,0,0};
  int count_file[3] = {mesh_nx3, mesh_nx2, mesh_nx1};
  int start_mem[3] = {0,0,0};
  int count_mem[3] = {mesh_nx3, mesh_nx2, mesh_nx1};

  //load data to user mesh data for later use
  AllocateRealUserMeshDataField(5);
  ruser_mesh_data[0].NewAthenaArray(mesh_nx3, mesh_nx2, mesh_nx1);
  ruser_mesh_data[1].NewAthenaArray(mesh_nx3, mesh_nx2, mesh_nx1);
  ruser_mesh_data[2].NewAthenaArray(mesh_nx3, mesh_nx2, mesh_nx1);
  ruser_mesh_data[3].NewAthenaArray(mesh_nx3, mesh_nx2, mesh_nx1);
  ruser_mesh_data[4].NewAthenaArray(mesh_nx3, mesh_nx2, mesh_nx1);
  HDF5ReadRealArray(input_filename.c_str(), "prim/rho", 3, start_file, count_file,
                    3, start_mem, count_mem, ruser_mesh_data[0], true);
  HDF5ReadRealArray(input_filename.c_str(), "prim/vel1", 3, start_file, count_file,
                    3, start_mem, count_mem, ruser_mesh_data[1], true);
  HDF5ReadRealArray(input_filename.c_str(), "prim/vel2", 3, start_file, count_file,
                    3, start_mem, count_mem, ruser_mesh_data[2], true);
  HDF5ReadRealArray(input_filename.c_str(), "prim/vel3", 3, start_file, count_file,
                    3, start_mem, count_mem, ruser_mesh_data[3], true);
  HDF5ReadRealArray(input_filename.c_str(), "prim/press", 3, start_file, count_file,
                    3, start_mem, count_mem, ruser_mesh_data[4], true);

  //Real dx1 = (mesh_x1max - mesh_x1min)/mesh_nx1;
  Real dx2 = (mesh_x2max - mesh_x2min)/mesh_nx2;
  Real dx3 = (mesh_x3max - mesh_x3min)/mesh_nx3;

  //prepare three vectors for index finding of x1 x2 x3 coordinates
  //the vector are equivalent to pcoord->x1v, x2v, x3v
  for(int i=0; i<mesh_nx1; i++){
    Real x1coord_now = (pow(x1ratio, i)-1.0)/(pow(x1ratio, mesh_nx1)-1.0) *
                       (mesh_x1max - mesh_x1min) + mesh_x1min;
    x1coord.push_back(x1coord_now);
  }
  for(int j=0; j<mesh_nx2; j++){
    x2coord.push_back(mesh_x2min+j*dx2);
  }
  for(int k=0; k<mesh_nx3; k++){
    x3coord.push_back(mesh_x3min+k*dx3);
  }

}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief monte carlo test problem generator
//! Inputs:
//! - pin: parameters
//! Outputs: (none)
//! Notes:
//! - reads in array using a slightly modified verions of from_array.cpp
//!   - NHYDRO
//!   - total number of MeshBlocks
//!   - MeshBlock/nx3
//!   - MeshBlock/nx2
//!   - MeshBlock/nx1

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  //  Initialize primitive
  for (int k=ks; k<=ke; ++k) {

    //load data
    Real z_now = pcoord->x3f(k);
    int index_znow = getindex(x3coord, z_now);

    for (int j=js; j<=je; ++j) {

      //load data
      Real y_now = pcoord->x2f(j);
      int index_ynow = getindex(x2coord, y_now);


      for (int i=is; i<=ie; ++i) {
        //load data
        Real x_now = pcoord->x1f(i);
        int index_xnow = getindex(x1coord, x_now);

        //printf("x_now:%g, y_now:%g, z_now:%g, index_x:%d, index_y:%d, index_z:%d\n", x_now, y_now, z_now, index_xnow, index_ynow, index_znow);

        phydro->w(IDN,k,j,i) = pmy_mesh->ruser_mesh_data[0](index_znow, index_ynow,
                                                            index_xnow);
        phydro->w(IVX,k,j,i) = pmy_mesh->ruser_mesh_data[1](index_znow, index_ynow,
                                                            index_xnow);
        phydro->w(IVY,k,j,i) = pmy_mesh->ruser_mesh_data[2](index_znow, index_ynow,
                                                            index_xnow);
        phydro->w(IVZ,k,j,i) = pmy_mesh->ruser_mesh_data[3](index_znow, index_ynow,
                                                            index_xnow);
        phydro->w(IPR,k,j,i) = pmy_mesh->ruser_mesh_data[4](index_znow, index_ynow,
                                                            index_xnow);

      }// end i
    }//end j
  }// end k

  // Set index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }
  // Initialize conserved
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, il, iu, jl, ju,
                             kl, ku);


}

//========================================================================================
//! \fn void MonteCarloBlock::MonteCarloProblemGenerator(ParameterInput *pin)
//! \brief Analogous to problem generator but used in support of InitializePhoton
//========================================================================================


//========================================================================================
//! \fn void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe)
//! \brief Initializes Photon packets before integration
//========================================================================================

void MonteCarloBlock::InitializePhoton(Photon *pphot, int ips, int ipe) {


  // Set initial cells and emission weights for all photon samples
  SetEmissionCellWeight(pphot,ips,ipe);

  for (int ip=ips; ip<=ipe; ip++) {

    // Obtain initial position within zone
    GetZonePosition(pphot,pran,pcoord,ip);

    // Set maximum integration time
    pphot->dtp[ip] = pphot->pmy_mcb->pmy_mc->tmax;

    //xs: store gas temperature
    pphot->user[0][ip] = tgas(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip]);

    if (emission_type == "freefree") {
      // Obtain intitial energy, polarization, direction and weight
      // Utilize free-free emission function in emission.cpp
      if(tnorm) {
        Real logtg = log(tgas(pphot->i3p[ip],pphot->i2p[ip],pphot->i1p[ip]));
        PhotonEmitFreeFree(this,pphot,logemin+logtg,logemax+logtg,ip);
      } else{
        PhotonEmitFreeFree(this,pphot,logemin,logemax,ip);
      }
    } else {
      pphot->ep[ip] = SampleEmissivity(this,pphot,ip);
      if (pphot->IsNanPhoton(ip))
        pphot->PrintPhoton("initialization: ",ip);
      //printf("en: %g\n",pphot->ep[ip]);
      //pphot->PrintPhoton("initialization: ",ip)

      // Generate initial angle parameters
      Real phi = 2. * PI * pran->uniform();
      Real cphi = cos(phi);
      Real sphi = sin(phi);
      Real cth = 2. * pran->uniform() - 1.;
      Real sth = sqrt(1. - SQR(cth));


      // Initialize wave vector with isotropic distribution
      pphot->k0p[ip] = 1.;
      pphot->k1p[ip] = sth*cphi;
      pphot->k2p[ip] = sth*sphi;
      pphot->k3p[ip] = cth;
    }

    if (pmy_mc->polarized) {
      // Initialize Stokes vector
      pphot->sip[ip] = 1.0;
      pphot->sup[ip] = 0.0;
      pphot->sqp[ip] = 0.0;
      pphot->svp[ip] = 0.0;
    }

    //xs: store photon energy
    pphot->user[1][ip] = pphot->ep[ip];

    // Set status flag
    if (pphot->wp[ip] < 0.0)
      pphot->statp[ip] = DESTROYED;
    else
      pphot->statp[ip] = EVOLVING;

    // initialize scattering number
    pphot->nscp[ip] = 0;

    // Initialize the absorption and scattering extinction coefficients
    // to the values appropriate in the emitted zone
    pphot->acp[ip] = AbsorptionOpacity(this,pphot,ip);
    pphot->scp[ip] = ScatteringOpacity(this,pphot,ip);

  } // loop over ip

}

//========================================================================================
//! \fn void MonteCarloBlock::FinalizePhoton(Photon *pphot, int ip)
//! \brief Complete work at end of photon packets before integration
//========================================================================================

void MonteCarloBlock::FinalizePhoton(Photon *pphot, int ip) {

  //xs: store scatter number
  pphot->user[2][ip] = pphot->nscp[ip];

}

namespace {

Real TableOpacity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  // Sets energy, temp, dens to table minimum if outside bounds
  Real le = log10(pphot->ep[ip]);
  le = (le < lmine) ? lmine : le;
  le = (le > lmaxe) ? lmaxe : le;
  Real ld = log10(pmcb->rho(i3,i2,i1));
  ld = (ld < lmind) ? lmind : ld;
  ld = (ld > lmaxd) ? lmaxd : ld;
  Real temp = pmcb->tgas(i3,i2,i1);
  Real lt = log10(temp);
  lt = (lt < lmint) ? lmint : lt;
  lt = (lt > lmaxt) ? lmaxt : lt;
  Real xi = (ld - lmind) / dld;
  Real xj = (lt - lmint) / dlt;
  Real xk = (le - lmine) / dle;
  int i = std::floor(xi);
  int j = std::floor(xj);
  int k = std::floor(xk);
  xi -= static_cast<Real>(i);
  xk -= static_cast<Real>(k);
  while ((j<ntem-2) && (temp_grid(j+1) < temp)){
    j++;
  }
  while ((j>0) && (temp_grid(j) > temp)){
    j--;
  }
  if(j > ntem-2) {
    j = ntem-2;
    temp = temp_grid(j+1);
  }
  xj = (temp-temp_grid(j))/(temp_grid(j+1)-temp_grid(j));

  Real opacl = (1.-xi) * ((1.-xj) * plan_tab(k,j,i) + xj * plan_tab(k,j+1,i))
                 + xi  * ((1.-xj) * plan_tab(k,j,i+1) + xj * plan_tab(k,j+1,i+1));
  Real opach = (1.-xi) * ((1.-xj) * plan_tab(k+1,j,i) + xj * plan_tab(k+1,j+1,i))
                 + xi  * ((1.-xj) * plan_tab(k+1,j,i+1) + xj * plan_tab(k+1,j+1,i+1));

  return (1.-xk) * opacl + xk * opach;

}

Real IntegrateEmission(Real temp, Real num, Real nup, Real am, Real ap) {

  int n = 20;
  Real h = 6.62607015e-27;
  Real dlnu = std::log(nup/num)/static_cast<Real>(n);
  Real dadnu = (ap-am)/(nup-num);
  Real lnu = std::log(num);
  Real sum = Planck(temp,num)*am*dlnu/h/2.;
  for(int i=1; i<n-1; ++i) {
    lnu += dlnu;
    Real nu = std::exp(lnu);
    Real alpha = dadnu*(nu-num)+am;
    sum += Planck(temp,nu)*alpha*dlnu/h;
  }
  sum += Planck(temp,nup)*ap*dlnu/h/2.;
  //if (sum < 0)
  //  printf("%g %g %g %g\n",num,nup,am,ap);
  return sum;
}

Real Planck(Real temp, Real nu) {

  Real h = 6.62607015e-27;
  Real c_cgs = 2.99792458e10;
  Real kb = 1.380649e-16;

  return 2.*h/c_cgs/c_cgs*pow(nu,3)/(std::exp(h*nu/kb/temp)-1.);

}

Real TableEmission(MonteCarloBlock *pmcb, int i3, int i2, int i1) {

  Real ld = log10(pmcb->rho(i3,i2,i1));
  ld = (ld < lmind) ? lmind : ld;
  ld = (ld > lmaxd) ? lmaxd : ld;
  Real temp = pmcb->tgas(i3,i2,i1);
  Real lt = log10(temp);
  lt = (lt < lmint) ? lmint : lt;
  lt = (lt > lmaxt) ? lmaxt : lt;
  Real xi = (ld - lmind) / dld;
  Real xj = (lt - lmint) / dlt;

  int i = std::floor(xi);
  int j = std::floor(xj);

  Real xi0 = xi;
  xi -= static_cast<Real>(i);
  while ((j<ntem-2) && (temp_grid(j+1) < temp)){
    j++;
  }
  while ((j>0) && (temp_grid(j) > temp)){
    j--;
  }
  if(j > ntem-2) {
    j = ntem-2;
    temp = temp_grid(j+1);
  }
  xj = (temp-temp_grid(j))/(temp_grid(j+1)-temp_grid(j));

  Real eta = (1.-xi) * ((1.-xj) * eta_tab(j,i) + xj * eta_tab(j+1,i))
               + xi  * ((1.-xj) * eta_tab(j,i+1) + xj * eta_tab(j+1,i+1));
  if (eta <= 0.) {
    printf("%d %d %g %g %g %g\n",i,j,xi,xj,xi0,eta);
  }
  return eta;
}


Real SampleEmissivity(MonteCarloBlock *pmcb, Photon *pphot, int ip) {

  Real dev = pmcb->pran->uniform();
  int i1 = pphot->i1p[ip];
  int i2 = pphot->i2p[ip];
  int i3 = pphot->i3p[ip];

  Real *prob = &(eta_cum_tab(i3,i2,i1,0));
  int i = mcbisect(dev,prob,nfre+1);
  Real a = (dev-prob[i])/(prob[i+1]-prob[i]);
  Real a1 = 1.-a;
  printf("%d %g %g\n",i,a,a1);
  //if (std::isinf(a)) {
  /*printf("%d %d %d\n",i3,i2,i1);
    for (int j=0; j< nfre+1; ++j)
      printf("%d %e\n",j,1-prob[j]);
      printf("%d %g %g %g %g\n",i,dev,fre_grid(i),a,a1);*/
  //}
  return a*fre_grid(i+1)+a1*fre_grid(i);
}

}
