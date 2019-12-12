#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <fftw3.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>


#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/init_baryon.c"
#include "init_emu.c"

double C_shear_tomo_sys(double ell,int z1,int z2);
double C_gl_tomo_sys(double ell,int zl,int zs);
void set_data_shear(int Ncl, double *ell, double *data, int start);
void set_data_ggl(int Ncl, double *ell, double *data, int start);
void set_data_clustering(int Ncl, double *ell, double *data, int start);
void compute_data_vector(char *details, double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double Q1, double Q2, double Q3);
double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz,double Q1, double Q2, double Q3);
void write_vector_wrapper(char *details, input_cosmo_params ic, input_nuisance_params in);
double log_like_wrapper(input_cosmo_params ic, input_nuisance_params in);
int get_N_tomo_shear(void);
int get_N_tomo_clustering(void);
int get_N_ggl(void);
int get_N_ell(void);
double log_L_clphotoz();
double log_L_wlphotoz();
double log_L_shear_calib();

int get_N_tomo_shear(void){
  return tomo.shear_Nbin;
}
int get_N_tomo_clustering(void){
  return tomo.clustering_Nbin;
}
int get_N_ggl(void){
  return tomo.ggl_Npowerspectra;
}
int get_N_ell(void){
  return like.Ncl;
}


double C_shear_tomo_sys(double ell, int z1, int z2)
{
  double C;
  // C= C_shear_tomo_nointerp(ell,z1,z2);
  // if(like.IA==1) C+=C_II_nointerp(ell,z1,z2)+C_GI_nointerp(ell,z1,z2);
  
  if(like.IA!=1) C= C_shear_tomo_nointerp(ell,z1,z2);
  //if(like.IA==1) C= C_shear_shear_IA(ell,z1,z2);
  if(like.IA==1) C = C_shear_tomo_nointerp(ell,z1,z2)+C_II_nointerp(ell,z1,z2)+C_GI_nointerp(ell,z1,z2);
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[z1])*(1.0+nuisance.shear_calibration_m[z2]);
  //printf("%le %d %d %le\n",ell,z1,z2,C_shear_tomo_nointerp(ell,z1,z2)+C_II_JB_nointerp(ell,z1,z2)+C_GI_JB_nointerp(ell,z1,z2));
return C;
}

double C_gl_tomo_sys(double ell,int zl,int zs)
{
  double C;
  // C=C_gl_tomo_nointerp(ell,zl,zs); 
  // if(like.IA==1) C += C_gI_nointerp(ell,zl,zs);
  
  if(like.IA!=1) C=C_gl_tomo_nointerp(ell,zl,zs);
  if(like.IA==1) C = C_ggl_IA(ell,zl,zs);
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
return C;
}
 

void set_data_shear(int Ncl, double *ell, double *data, int start)
{
  int i,z1,z2,nz;
  double a;
  for (nz = 0; nz < tomo.shear_Npowerspectra; nz++){
    z1 = Z1(nz); z2 = Z2(nz);
    for (i = 0; i < Ncl; i++){
      if (ell[i] < like.lmax_shear){ data[Ncl*nz+i] = C_shear_tomo_sys(ell[i],z1,z2);}
      else {data[Ncl*nz+i] = 0.;}
    }
  }
}

void set_data_ggl(int Ncl, double *ell, double *data, int start)
{
  int i, zl,zs,nz;  
  for (nz = 0; nz < tomo.ggl_Npowerspectra; nz++){
    zl = ZL(nz); zs = ZS(nz);
    for (i = 0; i < Ncl; i++){
      if (test_kmax(ell[i],zl)){
        data[start+(Ncl*nz)+i] = C_gl_tomo_sys(ell[i],zl,zs);
      }
      else{
        data[start+(Ncl*nz)+i] = 0.;
      }
    } 
  }
}

void set_data_clustering(int Ncl, double *ell, double *data, int start){
  int i, nz;
  for (nz = 0; nz < tomo.clustering_Npowerspectra; nz++){
    //printf("%d %e %e\n",nz, gbias.b[nz][1],pf_photoz(gbias.b[nz][1],nz));
    for (i = 0; i < Ncl; i++){
      if (test_kmax(ell[i],nz)){data[start+(Ncl*nz)+i] = C_cl_tomo_nointerp(ell[i],nz,nz);}
      else{data[start+(Ncl*nz)+i] = 0.;}
      //printf("%d %d %le %le\n",nz,nz,ell[i],data[Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + nz)+i]);
    }
  }
}



int set_cosmology_params(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu)
{
  cosmology.Omega_m=OMM;
  cosmology.Omega_v= 1.0-cosmology.Omega_m;
  cosmology.sigma_8=S8;
  cosmology.n_spec= NS;
  cosmology.w0=W0;
  cosmology.wa=WA;
  cosmology.omb=OMB;
  cosmology.h0=H0;
  cosmology.MGSigma=MGSigma;
  cosmology.MGmu=MGmu;

  if (cosmology.Omega_m < 0.05 || cosmology.Omega_m > 0.6) return 0;
  if (cosmology.omb < 0.04 || cosmology.omb > 0.055) return 0;
  if (cosmology.sigma_8 < 0.5 || cosmology.sigma_8 > 1.1) return 0;
  if (cosmology.n_spec < 0.84 || cosmology.n_spec > 1.06) return 0;
  if (cosmology.w0 < -2.1 || cosmology.w0 > -0.0) return 0;
  if (cosmology.wa < -2.6 || cosmology.wa > 2.6) return 0;
  if (cosmology.h0 < 0.4 || cosmology.h0 > 0.9) return 0;
  //CH BEGINS 
  //CH: to use for running planck15_BA0_w0_wa prior alone) 
  //printf("like_fourier.c from WFIRST_forecasts: cosmology bounds set for running with planck15_BA0_w0_wa prior\n");
  //if (cosmology.Omega_m < 0.05 || cosmology.Omega_m > 0.6) return 0; 
  //if (cosmology.omb < 0.01 || cosmology.omb > 0.1) return 0; 
  //if (cosmology.sigma_8 < 0.5 || cosmology.sigma_8 > 1.1) return 0; 
  //if (cosmology.n_spec < 0.84 || cosmology.n_spec > 1.06) return 0; 
  //if (cosmology.w0 < -2.1 || cosmology.w0 > 1.5) return 0; 
  //if (cosmology.wa < -5.0 || cosmology.wa > 2.6) return 0; 
  //if (cosmology.h0 < 0.3 || cosmology.h0 > 0.9) return 0; 
  //CH ENDS
  return 1;
}

void set_nuisance_shear_calib(double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10)
{
  nuisance.shear_calibration_m[0] = M1;
  nuisance.shear_calibration_m[1] = M2;
  nuisance.shear_calibration_m[2] = M3;
  nuisance.shear_calibration_m[3] = M4;
  nuisance.shear_calibration_m[4] = M5;
  nuisance.shear_calibration_m[5] = M6;
  nuisance.shear_calibration_m[6] = M7;
  nuisance.shear_calibration_m[7] = M8;
  nuisance.shear_calibration_m[8] = M9;
  nuisance.shear_calibration_m[9] = M10;
}

int set_nuisance_shear_photoz(double SP1,double SP2,double SP3,double SP4,double SP5,double SP6,double SP7,double SP8,double SP9,double SP10,double SPS1)
{
  int i;
  nuisance.bias_zphot_shear[0]=SP1;
  nuisance.bias_zphot_shear[1]=SP2;
  nuisance.bias_zphot_shear[2]=SP3;
  nuisance.bias_zphot_shear[3]=SP4;
  nuisance.bias_zphot_shear[4]=SP5;
  nuisance.bias_zphot_shear[5]=SP6;
  nuisance.bias_zphot_shear[6]=SP7;
  nuisance.bias_zphot_shear[7]=SP8;
  nuisance.bias_zphot_shear[8]=SP9;
  nuisance.bias_zphot_shear[9]=SP10;
  
  for (i=0;i<tomo.shear_Nbin; i++){ 
    nuisance.sigma_zphot_shear[i]=SPS1;
    if (nuisance.sigma_zphot_shear[i]<0.0001) return 0;
  }
  return 1;
}

int set_nuisance_clustering_photoz(double CP1,double CP2,double CP3,double CP4,double CP5,double CP6,double CP7,double CP8,double CP9,double CP10,double CPS1)
{
  int i;
  nuisance.bias_zphot_clustering[0]=CP1;
  nuisance.bias_zphot_clustering[1]=CP2;
  nuisance.bias_zphot_clustering[2]=CP3;
  nuisance.bias_zphot_clustering[3]=CP4;
  nuisance.bias_zphot_clustering[4]=CP5;
  nuisance.bias_zphot_clustering[5]=CP6;
  nuisance.bias_zphot_clustering[6]=CP7;
  nuisance.bias_zphot_clustering[7]=CP8;
  nuisance.bias_zphot_clustering[8]=CP9;
  nuisance.bias_zphot_clustering[9]=CP10;
  
  for (i=0;i<tomo.clustering_Nbin; i++){ 
    nuisance.sigma_zphot_clustering[i]=CPS1;
    if (nuisance.sigma_zphot_clustering[i]<0.0001) return 0;
  }
  return 1;
}

int set_nuisance_ia(double A_ia, double beta_ia, double eta_ia, double eta_ia_highz)
{
  nuisance.A_ia=A_ia;  
  nuisance.beta_ia=beta_ia;
  nuisance.eta_ia=eta_ia;
  nuisance.eta_ia_highz=eta_ia_highz;
  if (nuisance.A_ia < 0.0 || nuisance.A_ia > 10.0) return 0;
  if (nuisance.beta_ia < -1.0 || nuisance.beta_ia > 3.0) return 0;
  if (nuisance.eta_ia < -3.0 || nuisance.eta_ia> 3.0) return 0;
  if (nuisance.eta_ia_highz < -1.0 || nuisance.eta_ia_highz> 1.0) return 0;
  // if(like.IA!=0){
  //  if (check_LF()) return 0;
  // }
return 1;
}

int set_nuisance_gbias(double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8,double B9, double B10)
{
  int i;
  gbias.b[0] = B1;
  gbias.b[1] = B2;
  gbias.b[2] = B3;
  gbias.b[3] = B4;
  gbias.b[4] = B5;
  gbias.b[5] = B6;
  gbias.b[6] = B7;
  gbias.b[7] = B8;
  gbias.b[8] = B9;
  gbias.b[9] = B10;
  if(like.bias==1){
    for (i = 0; i < 10; i++){
      if (gbias.b[i] < 0.8 || gbias.b[i] > 3.0) return 0;
    }
  }
  return 1;
} 

double log_L_clphotoz()
{
  int i;
  double log_L = 0.;
  for (i=0;i<tomo.clustering_Nbin; i++){
      if (prior.bias_zphot_clustering[i][1] == 0.){
      printf("external_prior.c: called log_L_clphotoz while prior.bias_zphot_clustering[%d][1] not set.\nEXIT\n",i);
      exit(1);
    }
    log_L -= pow((nuisance.bias_zphot_clustering[i] - prior.bias_zphot_clustering[i][0])/ prior.bias_zphot_clustering[i][1],2.0);
  }
  if(redshift.clustering_photoz!=4) {
    if (prior.sigma_zphot_clustering[0][1] == 0.){
      printf("external_prior.c: called log_L_clphotoz while prior.sigma_zphot_clustering[0][1] not set.\nEXIT\n");
      exit(1);
    }
    log_L -= pow((nuisance.sigma_zphot_clustering[0] - prior.sigma_zphot_clustering[0][0])/ prior.sigma_zphot_clustering[0][1],2.0);
  }
  return 0.5*log_L;
}


double log_L_wlphotoz()
{
  int i;
  double log_L = 0.;
  for (i=0;i<tomo.shear_Nbin; i++){
    if (prior.bias_zphot_shear[i][1] == 0.){
      printf("external_prior.c: called log_L_wlphotoz while prior.bias_zphot_shear[%d][1] not set.\nEXIT\n",i);
      exit(1);
    }
    log_L -=  pow((nuisance.bias_zphot_shear[i] - prior.bias_zphot_shear[i][0])/ prior.bias_zphot_shear[i][1],2.0);
  }
  if(redshift.shear_photoz!=4) {
    if (prior.sigma_zphot_shear[0][1] == 0.){
      printf("external_prior.c: called log_L_wlphotoz while prior.sigma_zphot_shear[0][1] not set.\nEXIT\n");
      exit(1);
    }
    log_L -= pow((nuisance.sigma_zphot_shear[0] - prior.sigma_zphot_shear[0][0])/ prior.sigma_zphot_shear[0][1],2.0);
  }
  return 0.5*log_L;
}

double log_L_shear_calib()
{
  int i;
  double log_L = 0.; 
  for (i=0;i<tomo.shear_Nbin; i++){
    if (prior.shear_calibration_m[i][1] == 0.){
      printf("external_prior.c: called log_L_clphotoz while prior.shear_calibration_m[%d][1] not set.\nEXIT\n",i);
      exit(1);
    }
    log_L -= pow((nuisance.shear_calibration_m[i] - prior.shear_calibration_m[i][0])/ prior.shear_calibration_m[i][1],2.0);
  }
  return 0.5*log_L;
}


double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double Q1, double Q2, double Q3)
{
  int i,j,k,m=0,l;
  static double *pred;
  static double *ell;
  static double *ell_Cluster;
  static double darg;
  double chisqr,a,log_L_prior=0.0, log_L=0.0;
  
  if(ell==0){
    pred= create_double_vector(0, like.Ndata-1);
    ell= create_double_vector(0, like.Ncl-1);
    darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
    for (l=0;l<like.Ncl;l++){
      ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
    }
  }
  if (set_cosmology_params(OMM,S8,NS,W0,WA,OMB,H0,MGSigma,MGmu)==0){
    printf("Cosmology out of bounds\n");
    return -1.0e12;
  }
  set_nuisance_shear_calib(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10);
  if (set_nuisance_shear_photoz(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9,SP10,SPS1)==0){
    printf("Shear photo-z sigma too small\n");
    return -1.0e12;
  }
  if (set_nuisance_clustering_photoz(CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,CPS1)==0){
    printf("Clustering photo-z sigma too small\n");
    return -1.0e12;
  }
  if (set_nuisance_ia(A_ia,beta_ia,eta_ia,eta_ia_highz)==0){
    printf("IA parameters out of bounds\n");
    return -1.0e12; 
  }
  if (set_nuisance_gbias(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10)==0){
    printf("Bias out of bounds\n");
    return -1.0e12;
  }
       
  //printf("like %le %le %le %le %le %le %le %le\n",cosmology.Omega_m, cosmology.Omega_v,cosmology.sigma_8,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,cosmology.h0); 
  // printf("like %le %le %le %le\n",gbias.b[0][0], gbias.b[1][0], gbias.b[2][0], gbias.b[3][0]);    
  // for (i=0; i<10; i++){
  //   printf("nuisance %le %le %le\n",nuisance.shear_calibration_m[i],nuisance.bias_zphot_shear[i],nuisance.sigma_zphot_shear[i]);
  // }

  log_L_prior=0.0;
  if(like.wlphotoz!=0) log_L_prior+=log_L_wlphotoz();
  if(like.clphotoz!=0) log_L_prior+=log_L_clphotoz();
  if(like.shearcalib==1) log_L_prior+=log_L_shear_calib();
  
  if(like.IA!=0) {
    log_L = 0.0;
    log_L -= pow((nuisance.A_ia - prior.A_ia[0])/prior.A_ia[1],2.0);
    log_L -= pow((nuisance.beta_ia - prior.beta_ia[0])/prior.beta_ia[1],2.0);
    log_L -= pow((nuisance.eta_ia - prior.eta_ia[0])/prior.eta_ia[1],2.0);
    log_L -= pow((nuisance.eta_ia_highz - prior.eta_ia_highz[0])/prior.eta_ia_highz[1],2.0);
    log_L_prior+=0.5*log_L;
  }
  if(like.baryons==1){
    log_L = 0.0;
    log_L -= pow((Q1 - prior.bary_Q1[0])/prior.bary_Q1[1],2.0);
    log_L -= pow((Q2 - prior.bary_Q2[0])/prior.bary_Q2[1],2.0);
    log_L -= pow((Q3 - prior.bary_Q3[0])/prior.bary_Q3[1],2.0);
    log_L_prior+=0.5*log_L;
  }
 
  // printf("%d %d %d %d\n",like.BAO,like.wlphotoz,like.clphotoz,like.shearcalib);
  // printf("logl %le %le %le %le\n",log_L_shear_calib(),log_L_wlphotoz(),log_L_clphotoz(),log_L_clusterMobs());
  int start=0;  
  
  if(like.shear_shear==1) {
    set_data_shear(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.shear_Npowerspectra;
  }
  if(like.shear_pos==1){
    set_data_ggl(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.ggl_Npowerspectra;
  } 
  if(like.pos_pos==1){
    set_data_clustering(like.Ncl,ell,pred, start);
    start=start+like.Ncl*tomo.clustering_Npowerspectra;
  }

  chisqr=0.0;
  for (i=0; i<like.Ndata; i++){
    for (j=0; j<like.Ndata; j++){
      a=(pred[i]-data_read(1,i)+Q1*bary_read(1,1,i)+Q2*bary_read(1,2,i))*invcov_read(1,i,j)*(pred[j]-data_read(1,j)+Q1*bary_read(1,1,j)+Q2*bary_read(1,2,j));
      chisqr=chisqr+a;
    }
    // if (fabs(data_read(1,i)) < 1.e-25){
    //    printf("%d %le %le %le\n",i,data_read(1,i),pred[i],invcov_read(1,i,i));
    // }
  }
  if (chisqr<0.0){
    printf("error: chisqr = %le\n",chisqr);
    //exit(EXIT_FAILURE);
  }
//  printf("%le\n",chisqr);
  return -0.5*chisqr+log_L_prior;
}

void compute_data_vector(char *details, double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5,double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double Q1, double Q2, double Q3)
{

  int i,j,k,m=0,l;
  static double *pred;
  static double *ell;
  static double *ell_Cluster;
  static double darg;
  double chisqr,a,log_L_prior=0.0;
  
  if(ell==0){
    pred= create_double_vector(0, like.Ndata-1);
    ell= create_double_vector(0, like.Ncl-1);
    darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
    for (l=0;l<like.Ncl;l++){
      ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
    }
  }
  set_cosmology_params(OMM,S8,NS,W0,WA,OMB,H0,MGSigma,MGmu);
  set_nuisance_shear_calib(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10);
  set_nuisance_shear_photoz(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9,SP10,SPS1);
  set_nuisance_clustering_photoz(CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,CPS1);
  set_nuisance_ia(A_ia,beta_ia,eta_ia,eta_ia_highz);
  set_nuisance_gbias(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10);
   
  int start=0;  
  if(like.shear_shear==1) {
    set_data_shear(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.shear_Npowerspectra;
  }
  if(like.shear_pos==1){
    //printf("ggl\n");
    set_data_ggl(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.ggl_Npowerspectra;
  } 
  if(like.pos_pos==1){
    //printf("clustering\n");
    set_data_clustering(like.Ncl,ell,pred, start);
    start=start+like.Ncl*tomo.clustering_Npowerspectra;
  }


  FILE *F;
  char filename[300];
  if (strstr(details,"FM") != NULL){
    sprintf(filename,"%s",details);
  }
  else {sprintf(filename,"datav/%s_%s_%s",survey.name,like.probes,details);}
  F=fopen(filename,"w");
  for (i=0;i<like.Ndata; i++){  
    fprintf(F,"%d %le\n",i,pred[i]);
    //printf("%d %le\n",i,pred[i]);
  }
  fclose(F);
}


void write_vector_wrapper(char *details, input_cosmo_params ic, input_nuisance_params in)
{
  compute_data_vector(details, ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.h0, ic.MGSigma, ic.MGmu,
    in.bias[0], in.bias[1], in.bias[2], in.bias[3],in.bias[4], in.bias[5], in.bias[6], in.bias[7],in.bias[8], in.bias[9], 
    in.source_z_bias[0], in.source_z_bias[1], in.source_z_bias[2], in.source_z_bias[3], in.source_z_bias[4], 
    in.source_z_bias[5], in.source_z_bias[6], in.source_z_bias[7], in.source_z_bias[8], in.source_z_bias[9], 
    in.source_z_s, 
    in.lens_z_bias[0], in.lens_z_bias[1], in.lens_z_bias[2], in.lens_z_bias[3], in.lens_z_bias[4], 
    in.lens_z_bias[5], in.lens_z_bias[6], in.lens_z_bias[7], in.lens_z_bias[8], in.lens_z_bias[9], 
    in.lens_z_s, 
    in.shear_m[0], in.shear_m[1], in.shear_m[2], in.shear_m[3], in.shear_m[4], 
    in.shear_m[5], in.shear_m[6], in.shear_m[7], in.shear_m[8], in.shear_m[9], 
    in.A_ia, in.beta_ia, in.eta_ia, in.eta_ia_highz,in.bary[0], in.bary[1], in.bary[2]);
}

double log_like_wrapper(input_cosmo_params ic, input_nuisance_params in)
{
  double like = log_multi_like(ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.h0, ic.MGSigma, ic.MGmu,
    in.bias[0], in.bias[1], in.bias[2], in.bias[3],in.bias[4], in.bias[5], in.bias[6], in.bias[7],in.bias[8], in.bias[9], 
    in.source_z_bias[0], in.source_z_bias[1], in.source_z_bias[2], in.source_z_bias[3], in.source_z_bias[4], 
    in.source_z_bias[5], in.source_z_bias[6], in.source_z_bias[7], in.source_z_bias[8], in.source_z_bias[9], 
    in.source_z_s, 
    in.lens_z_bias[0], in.lens_z_bias[1], in.lens_z_bias[2], in.lens_z_bias[3], in.lens_z_bias[4], 
    in.lens_z_bias[5], in.lens_z_bias[6], in.lens_z_bias[7], in.lens_z_bias[8], in.lens_z_bias[9], 
    in.lens_z_s, 
    in.shear_m[0], in.shear_m[1], in.shear_m[2], in.shear_m[3], in.shear_m[4], 
    in.shear_m[5], in.shear_m[6], in.shear_m[7], in.shear_m[8], in.shear_m[9], 
    in.A_ia, in.beta_ia, in.eta_ia, in.eta_ia_highz,in.bary[0], in.bary[1], in.bary[2]); 
  return like;
}



void save_zdistr_sources(int zs){
  double z,dz =(redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/300.0;
  printf("Printing redshift distribution n(z) for source redshift bin %d\n",zs);
  
   FILE *F1;
   char filename[300];
   sprintf(filename,"zdistris/zdist_sources_bin%d.txt",zs);
   F1 = fopen(filename,"w");
   for (z =redshift.shear_zdistrpar_zmin; z< redshift.shear_zdistrpar_zmax; z+= dz){
      fprintf(F1,"%e %e\n", z, zdistr_photoz(z,zs));
   }
}


void save_zdistr_lenses(int zl){
   double z,dz =(redshift.clustering_zdistrpar_zmax-redshift.clustering_zdistrpar_zmin)/300.0;
  printf("Printing redshift distribution n(z) and bias b(z) for lens redshift bin %d\n",zl);
   
   FILE *F1;
   char filename[300];
   sprintf(filename,"zdistris/zdist_lenses_bin%d.txt", zl);
   F1 = fopen(filename,"w");
   for (z =redshift.clustering_zdistrpar_zmin; z< redshift.clustering_zdistrpar_zmax; z+= dz){
      fprintf(F1,"%e %e\n", z, pf_photoz(z,zl));
   }
}


 int main(int argc, char** argv)
{
  clock_t begin, end;
  double time_spent;
  int i;
  char arg1[400],arg2[400],arg3[400];
/* here, do your time-consuming job */
  int sce=atoi(argv[1]);

  int N_scenarios=36;
  
  double area_table[36]={7623.22,14786.3,9931.47,8585.43,17681.8,15126.9,9747.99,8335.08,9533.42,18331.3,12867.8,17418.9,19783.1,12538.8,15260.0,16540.7,19636.8,11112.7,10385.5,16140.2,18920.1,17976.2,11352.0,9214.77,16910.7,11995.6,16199.8,14395.1,8133.86,13510.5,19122.3,15684.5,12014.8,14059.7,10919.3,13212.7};
  
  double nsource_table[36]={13.991,33.3975,17.069,28.4875,35.3643,10.3802,11.1105,29.5904,31.4608,15.7463,13.3066,9.07218,9.58719,16.3234,25.8327,10.7415,38.0412,32.8821,19.8893,27.0369,15.1711,14.2418,19.1851,26.926,22.0012,12.6553,18.6304,11.9787,36.8728,22.5265,17.4381,12.3424,10.0095,23.5979,20.8771,24.8956};
  
  double nlens_table[36]={22.9726 ,61.5948 ,28.7809 ,51.4354 ,65.7226 ,16.3777 ,17.6899 ,53.6984 ,57.562 ,26.266 ,21.703 ,14.0587 ,14.9667 ,27.36 ,46.0364 ,17.0253 ,71.39 ,60.5184 ,34.2283 ,48.4767 ,25.1812 ,23.44 ,32.8578 ,48.2513 ,38.3766 ,20.5029 ,31.783 ,19.2647 ,68.9094 ,39.4168 ,29.4874 ,19.9292 ,15.7162 ,41.5487 ,36.1616 ,44.148};

char survey_designation[1][200]={"LSST_Y10"}; 

char source_zfile[36][400]={"wl_redshift_model0_WLz01.880307e-01_WLalpha8.485694e-01.txt", "wl_redshift_model1_WLz01.731166e-01_WLalpha7.662434e-01.txt", "wl_redshift_model2_WLz01.846221e-01_WLalpha8.297540e-01.txt", "wl_redshift_model3_WLz01.758423e-01_WLalpha7.812893e-01.txt", "wl_redshift_model4_WLz01.721357e-01_WLalpha7.608290e-01.txt", "wl_redshift_model5_WLz01.931476e-01_WLalpha8.768147e-01.txt", "wl_redshift_model6_WLz01.919821e-01_WLalpha8.703814e-01.txt", "wl_redshift_model7_WLz01.751912e-01_WLalpha7.776952e-01.txt", "wl_redshift_model8_WLz01.741405e-01_WLalpha7.718958e-01.txt", "wl_redshift_model9_WLz01.860048e-01_WLalpha8.373865e-01.txt", "wl_redshift_model10_WLz01.888904e-01_WLalpha8.533151e-01.txt", "wl_redshift_model11_WLz01.954564e-01_WLalpha8.895592e-01.txt", "wl_redshift_model12_WLz01.945099e-01_WLalpha8.843347e-01.txt", "wl_redshift_model13_WLz01.853877e-01_WLalpha8.339802e-01.txt", "wl_redshift_model14_WLz01.775192e-01_WLalpha7.905458e-01.txt", "wl_redshift_model15_WLz01.925611e-01_WLalpha8.735775e-01.txt", "wl_redshift_model16_WLz01.708849e-01_WLalpha7.539247e-01.txt", "wl_redshift_model17_WLz01.733832e-01_WLalpha7.677150e-01.txt", "wl_redshift_model18_WLz01.820009e-01_WLalpha8.152851e-01.txt", "wl_redshift_model19_WLz01.767381e-01_WLalpha7.862345e-01.txt", "wl_redshift_model20_WLz01.866426e-01_WLalpha8.409073e-01.txt", "wl_redshift_model21_WLz01.877261e-01_WLalpha8.468882e-01.txt", "wl_redshift_model22_WLz01.826188e-01_WLalpha8.186960e-01.txt", "wl_redshift_model23_WLz01.768086e-01_WLalpha7.866234e-01.txt", "wl_redshift_model24_WLz01.802711e-01_WLalpha8.057362e-01.txt", "wl_redshift_model25_WLz01.897506e-01_WLalpha8.580632e-01.txt", "wl_redshift_model26_WLz01.831217e-01_WLalpha8.214720e-01.txt", "wl_redshift_model27_WLz01.906926e-01_WLalpha8.632630e-01.txt", "wl_redshift_model28_WLz01.714197e-01_WLalpha7.568766e-01.txt", "wl_redshift_model29_WLz01.798667e-01_WLalpha8.035040e-01.txt", "wl_redshift_model30_WLz01.842554e-01_WLalpha8.277299e-01.txt", "wl_redshift_model31_WLz01.901797e-01_WLalpha8.604322e-01.txt", "wl_redshift_model32_WLz01.937710e-01_WLalpha8.802561e-01.txt", "wl_redshift_model33_WLz01.790701e-01_WLalpha7.991072e-01.txt", "wl_redshift_model34_WLz01.811701e-01_WLalpha8.106987e-01.txt", "wl_redshift_model35_WLz01.781525e-01_WLalpha7.940419e-01.txt"};

char lens_zfile[36][400]={"LSS_redshift_model0_LSSz02.629496e-01_LSSalpha9.285983e-01.txt","LSS_redshift_model1_LSSz02.852923e-01_LSSalpha8.985943e-01.txt","LSS_redshift_model2_LSSz02.664822e-01_LSSalpha9.186036e-01.txt","LSS_redshift_model3_LSSz02.798758e-01_LSSalpha9.014201e-01.txt","LSS_redshift_model4_LSSz02.873873e-01_LSSalpha8.978683e-01.txt","LSS_redshift_model5_LSSz02.593970e-01_LSSalpha9.470921e-01.txt","LSS_redshift_model6_LSSz02.600213e-01_LSSalpha9.425114e-01.txt","LSS_redshift_model7_LSSz02.811155e-01_LSSalpha9.006370e-01.txt","LSS_redshift_model8_LSSz02.831875e-01_LSSalpha8.995165e-01.txt","LSS_redshift_model9_LSSz02.649368e-01_LSSalpha9.224338e-01.txt","LSS_redshift_model10_LSSz02.622058e-01_LSSalpha9.314128e-01.txt","LSS_redshift_model11_LSSz02.584820e-01_LSSalpha9.568082e-01.txt","LSS_redshift_model12_LSSz02.588053e-01_LSSalpha9.527220e-01.txt","LSS_redshift_model13_LSSz02.656075e-01_LSSalpha9.206866e-01.txt","LSS_redshift_model14_LSSz02.768397e-01_LSSalpha9.037492e-01.txt","LSS_redshift_model15_LSSz02.596975e-01_LSSalpha9.447600e-01.txt","LSS_redshift_model16_LSSz02.901709e-01_LSSalpha8.971658e-01.txt","LSS_redshift_model17_LSSz02.847362e-01_LSSalpha8.988183e-01.txt","LSS_redshift_model18_LSSz02.698330e-01_LSSalpha9.121821e-01.txt","LSS_redshift_model19_LSSz02.782257e-01_LSSalpha9.026084e-01.txt","LSS_redshift_model20_LSSz02.642756e-01_LSSalpha9.243038e-01.txt","LSS_redshift_model21_LSSz02.632273e-01_LSSalpha9.276296e-01.txt","LSS_redshift_model22_LSSz02.689934e-01_LSSalpha9.135968e-01.txt","LSS_redshift_model23_LSSz02.780987e-01_LSSalpha9.027073e-01.txt","LSS_redshift_model24_LSSz02.723464e-01_LSSalpha9.085463e-01.txt","LSS_redshift_model25_LSSz02.615210e-01_LSSalpha9.343470e-01.txt","LSS_redshift_model26_LSSz02.683327e-01_LSSalpha9.147934e-01.txt","LSS_redshift_model27_LSSz02.608392e-01_LSSalpha9.376962e-01.txt","LSS_redshift_model28_LSSz02.889655e-01_LSSalpha8.974355e-01.txt","LSS_redshift_model29_LSSz02.729686e-01_LSSalpha9.077654e-01.txt","LSS_redshift_model30_LSSz02.669178e-01_LSSalpha9.176391e-01.txt","LSS_redshift_model31_LSSz02.612016e-01_LSSalpha9.358553e-01.txt","LSS_redshift_model32_LSSz02.591077e-01_LSSalpha9.496317e-01.txt","LSS_redshift_model33_LSSz02.742326e-01_LSSalpha9.063038e-01.txt","LSS_redshift_model34_LSSz02.710103e-01_LSSalpha9.103760e-01.txt","LSS_redshift_model35_LSSz02.757517e-01_LSSalpha9.047459e-01.txt"};

double shear_prior[36]={0.00891915,0.0104498 ,0.0145972 ,0.0191916 ,0.00450246 ,0.00567828 ,0.00294841 ,0.00530922 ,0.0118632 ,0.0151849 ,0.00410151 ,0.0170622 ,0.0197331 ,0.0106615 ,0.0124445 ,0.00994507 ,0.0136251 ,0.0143491 ,0.0164314 ,0.016962 ,0.0186608 ,0.00945903 ,0.0113246 ,0.0155225 ,0.00800846 ,0.00732104 ,0.00649453 ,0.00243976 ,0.0125932 ,0.0182587 ,0.00335859 ,0.00682287 ,0.0177269 ,0.0035219 ,0.00773304 ,0.0134886};

double delta_z_prior[36]={0.0032537,0.00135316,0.00168787,0.00215043,0.00406031,0.00222358,0.00334993,0.00255186,0.00266499,0.00159226,0.00183664,0.00384965,0.00427765,0.00314377,0.00456113,0.00347868,0.00487938,0.00418152,0.00469911,0.00367598,0.0028009,0.00234161,0.00194964,0.00200982,0.00122739,0.00310886,0.00275168,0.00492736,0.00437241,0.00113931,0.00104864,0.00292328,0.00452082,0.00394114,0.00150756,0.003613};

double sigma_z[36]={0.0849973 ,0.0986032 ,0.0875521 ,0.0968222 ,0.0225239 ,0.0718278 ,0.0733675 ,0.0385274 ,0.0425549 ,0.0605867 ,0.0178555 ,0.0853407 ,0.0124119 ,0.0531027 ,0.0304032 ,0.0503145 ,0.0132213 ,0.0941765 ,0.0416444 ,0.0668198 ,0.063227 ,0.0291332 ,0.0481633 ,0.0595606 ,0.0818742 ,0.0472518 ,0.0270185 ,0.0767401 ,0.0219945 ,0.0902663 ,0.0779705 ,0.0337666 ,0.0362358 ,0.0692429 ,0.0558841 ,0.0150457};

double sigma_z_prior[36]={0.00331909,0.00529541,0.00478151,0.00437497,0.00443062,0.00486333,0.00467423,0.0036723,0.00426963,0.00515357,0.0054553,0.00310132,0.00305971,0.00406327,0.00594293,0.00348709,0.00562526,0.00396025,0.00540537,0.00500447,0.00318595,0.00460592,0.00412137,0.00336418,0.00524988,0.00390092,0.00498349,0.0056667,0.0036384,0.00455861,0.00554822,0.00381061,0.0057615,0.00357705,0.00590572,0.00422393};


  init_cosmo_runmode("halofit");
  init_bary("TNG100");
  init_binning_fourier(15,20.0,3000.0,3000.0,21.0,10,10);
  init_priors(shear_prior[sce],sigma_z[sce],delta_z_prior[sce],sigma_z_prior[sce]);
  init_survey("LSST",nsource_table[sce],nlens_table[sce],area_table[sce]);
  sprintf(arg1,"zdistris/%s",source_zfile[sce]);
  sprintf(arg2,"zdistris/%s",lens_zfile[sce]); 
  init_galaxies(arg1,arg2,"gaussian","gaussian","LSST_Y10");
  init_IA("NLA_HF","GAMA"); 
  init_probes("3x2pt");



sprintf(arg3,"TNG100_Y10_area%le",area_table[sce]);
compute_data_vector(arg3,0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,gbias.b[0],gbias.b[1],gbias.b[2],gbias.b[3],gbias.b[4],gbias.b[5],gbias.b[6],gbias.b[7],gbias.b[8],gbias.b[9],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,sigma_z[sce],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,sigma_z[sce]*0.6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0);

// compute_data_vector("mu1_Sigma0",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,1.,1.35,1.5,1.65,1.8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);
  // compute_data_vector("mu1_Sigma1",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,1.,1.,1.35,1.5,1.65,1.8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);

  // init_data_inv("cov/WFIRST_3x2pt_clusterN_clusterWL_inv","datav/WFIRST_all_2pt_clusterN_clusterWL_fid");
  

  // begin = clock();
  // log_multi_like(0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  // // printf("knonlin %le\n",nonlinear_scale_computation(1.0));
  // // printf("knonlin %le\n",nonlinear_scale_computation(0.5));
  // end = clock();
  // time_spent = (double)(end - begin) / CLOCKS_PER_SEC;      
  // printf("timespent %le\n",time_spent);
  
  //CH BEGINS
  //for testing Planck15_BAO_w0wa prior alone
  //compute_data_vector("fid",3.50989e-01,8.04675e-01,9.64061e-01,-5.05518e-01,-1.46884e+00,5.46245e-02,6.39839e-01,0.,0.,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);
  //expect 0.0 for the following
  //log_multi_like(3.50989e-01,8.04675e-01.1,9.64061e-01,-5.05518e-01,-1.46884e+00,5.46245e-02,6.39839e-01,0.,0.,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);
  //expect -13.195605 for the following
  //log_multi_like(0.35449914, 0.81272201,0.97370111, -0.51057289,-1.48353327,0.05517077,0.64623714,0.,0.,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9); 
  //CH ENDS
  
  return 0;
}


