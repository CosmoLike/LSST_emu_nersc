#include "IA.h"

void set_LF_GAMA(void);
void set_LF_DEEP2(void);
int check_LF (void);
double M_abs(double mag, double a); //absolute magnitude corresponding to aparent magnitude at a; h = 1 units, incl k-corrections
double n_all_LF(double mag, double a);
double f_red_LF(double mag, double a);
double A_LF(double mag, double a);

double W_source(double a, double nz);

double A_IA_Joachimi(double a);
double P_II (double k, double a);
double P_dI (double k, double a);
double P_II_JB(double k, double a);
double P_dI_JB(double k, double a);
/*** fast routines for total shear-shear + ggl signal for like.IA = 1,3****/
double C_ggl_IA(double s, int ni, int nj);

double C_gI_nointerp(double s, int ni, int nj); //gI term using NL-LA P_deltaI power spectrum

double C_II_nointerp(double s, int ni, int nj);
double C_GI_nointerp(double s, int ni, int nj);


/******** select LF ************/
void set_LF_GAMA(void){
  int i;
  for (i = 0; i <5; i++){
    LF_coefficients[0][i] =LF_coefficients_GAMA[0][i];
    LF_coefficients[1][i] =LF_coefficients_GAMA[1][i];
  }
}
void set_LF_DEEP2(void){
  int i;
  for (i = 0; i <5; i++){
    LF_coefficients[0][i] =LF_coefficients_DEEP2[0][i];
    LF_coefficients[1][i] =LF_coefficients_DEEP2[1][i];
  }
}
/************** normalization rouintes *********************/
double M_abs(double mag, double a){ //in h = 1 units, incl. Poggianti 1997 k+e-corrections
  static double *table;
  static double dz;
  double ke;
  if (table ==0){
    int i;
    //read in + tabulate k+e corrections for early types, restframe r band
    //interpolated from http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A%2BAS/122/399
    table = create_double_vector(0,30);
    dz = 0.1;
    /*FILE *ein;
    double d1,d2,d3;
    ein = fopen(survey.Kcorrect_File,"r");
    EXIT_MISSING_FILE(ein, "LoverL0",survey.Kcorrect_File)*/
    for (i = 0; i< 31; i++){
   /*   fscanf(ein,"%le %le %le\n",&d1,&d2,&d3);
       table[i] = d2+d3; */
      table[i] = KE[i];
    }
   /* fclose(ein);*/
  }
  double z=1./a-1.0;
  if(z >= 3.0) {z=2.99;}  // no acceptable k-korrection exists for k>3, also no meaningful IA model
  ke = interpol(table, 31, 0., 3.0, 0.1, z, 1.0, 1.0);
  
  return mag - 5.0*log10(f_K(chi(a))/a*cosmology.coverH0) -25.0 -ke;
}
double n_all_LF(double mag, double a){
   double Q, alpha, Mstar, Phistar,P,Mlim, LF_all[3];
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}
  
  //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
  //all galaxies
  alpha =LF_coefficients[1][1]+ nuisance.LF_alpha;
  Mstar = LF_coefficients[1][2];
  Q = LF_coefficients[1][3]+nuisance.LF_Q;
  P = LF_coefficients[1][4]+nuisance.LF_P;
  Phistar = LF_coefficients[1][0];
  
  LF_all[0] = Phistar*pow(10.0,0.4*P*(1./a-1));
  LF_all[1] = Mstar -Q*(1./a-1. - 0.1);
  LF_all[2]= alpha;
  Mlim = M_abs(mag,a); //also in h = 1 units
  
   
  return LF_all[0]*gsl_sf_gamma_inc(LF_all[2]+1,pow(10.0,-0.4*(Mlim-LF_all[1])));
}

double f_red_LF(double mag, double a){
  double Q, alpha, Mstar, Phistar,P,Q_red, alpha_red, Mstar_red, Phistar_red,P_red,Mlim, LF_all[3], LF_red[3];
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}

  //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
  //red galaxies
  alpha_red = LF_coefficients[0][1] + nuisance.LF_red_alpha;
  Mstar_red = LF_coefficients[0][2]; //in h = 1 units
  Q_red = LF_coefficients[0][3]+nuisance.LF_red_Q;
  P_red = LF_coefficients[0][4]+nuisance.LF_red_P;
  Phistar_red = LF_coefficients[0][0];
  
  LF_red[0] = Phistar_red*pow(10.0,0.4*P_red*(1./a-1));
  LF_red[1] = Mstar_red-Q_red*(1./a-1. - 0.1);
  LF_red[2]= alpha_red;
  
  //all galaxies
  alpha =LF_coefficients[1][1]+ nuisance.LF_alpha;
  Mstar = LF_coefficients[1][2];
  Q = LF_coefficients[1][3]+nuisance.LF_Q;
  P = LF_coefficients[1][4]+nuisance.LF_P;
  Phistar = LF_coefficients[1][0];
  
  LF_all[0] = Phistar*pow(10.0,0.4*P*(1./a-1));
  LF_all[1] = Mstar -Q*(1./a-1. - 0.1);
  LF_all[2]= alpha;
  Mlim = M_abs(mag,a); //also in h = 1 units
  
  return LF_red[0]/LF_all[0]*gsl_sf_gamma_inc(LF_red[2]+1,pow(10.0,-0.4*(Mlim-LF_red[1])))/gsl_sf_gamma_inc(LF_all[2]+1,pow(10.0,-0.4*(Mlim-LF_all[1])));
}

double A_LF(double mag, double a){// averaged (L/L_0)^beta over red galaxy LF
    double Q_red, alpha_red, Mstar_red,x,Mlim, LF_red[3], Lstar,L0;
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}
    
    //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
    //red galaxies
  alpha_red = LF_coefficients[0][1]+ nuisance.LF_red_alpha;
  Mstar_red = LF_coefficients[0][2]; //in h = 1 units
  Q_red = LF_coefficients[0][3]+nuisance.LF_red_Q;
  
    LF_red[1] = Mstar_red-Q_red*(1./a-1. -.1);
    LF_red[2]= alpha_red;
  
    Mlim = M_abs(mag,a);
    Lstar = pow(10.,-0.4*LF_red[1]);
    x = pow(10.,-0.4*(Mlim-LF_red[1])); //Llim/Lstar
    L0 =pow(10.,-0.4*(-22.)); //all in h = 1 units
    return pow(Lstar/L0, nuisance.beta_ia)*gsl_sf_gamma_inc(LF_red[2]+nuisance.beta_ia+1,x)/gsl_sf_gamma_inc(LF_red[2]+1,x);
}


double A_LF_all(double mag, double a){// averaged (L/L_0)^beta over red galaxy LF
    double Q_red, alpha_red, Mstar_red,x,Mlim, LF_red[3], Lstar,L0;
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}
    
    //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
    //all galaxies
  alpha_red = LF_coefficients[1][1]+ nuisance.LF_alpha;
  Mstar_red = LF_coefficients[1][2]; //in h = 1 units
  Q_red = LF_coefficients[1][3]+nuisance.LF_Q;
  
    LF_red[1] = Mstar_red-Q_red*(1./a-1. -.1);
    LF_red[2]= alpha_red;
  
    Mlim = M_abs(mag,a);
    Lstar = pow(10.,-0.4*LF_red[1]);
    x = pow(10.,-0.4*(Mlim-LF_red[1])); //Llim/Lstar
    L0 =pow(10.,-0.4*(-22.)); //all in h = 1 units
    return pow(Lstar/L0, nuisance.beta_ia)*gsl_sf_gamma_inc(LF_red[2]+nuisance.beta_ia+1,x)/gsl_sf_gamma_inc(LF_red[2]+1,x);
}


int check_LF(void){ //return 1 if combination of all + red galaxy LF parameters is unphysical, i.e. if f_red > 1 for some z < redshift.shear_zdistrpar_zmax
  double a=1./(1+redshift.shear_zdistrpar_zmax)+0.005;
  while (a < 1.){
    if( M_abs(survey.m_lim,a) < LF_coefficients[1][2]-(LF_coefficients[1][3]+nuisance.LF_Q)*(1./a-1. - 0.1) || M_abs(survey.m_lim,a) < LF_coefficients[0][2]-(LF_coefficients[0][3]+nuisance.LF_red_Q)*(1./a-1. - 0.1)){return 1;}
    if (f_red_LF(survey.m_lim,a) > 1.0){return 1;}
    a+=0.01;
  }
  return 0;
}

/*=========================================================*/
/*=============  Intrinsic Alignment models  ==============*/
/*=========================================================*/


double A_IA_Joachimi(double a){
  double z, A_red, highz = 0.75;
  z = 1./a-1;
  A_red = nuisance.A_ia*A_LF(survey.m_lim,a)*f_red_LF(survey.m_lim,a); //A_0*<(L/L_0)^beta>*f_red
  if (a < 1./(1.+highz)){ //z > highz, factor in uncertainty in extrapolation of redshift scaling
    return A_red*pow((1.+z)/nuisance.oneplusz0_ia,nuisance.eta_ia)*pow((1.+z)/(1.+highz),nuisance.eta_ia_highz);
  }
  else{//standard redshift scaling
    return A_red*pow((1.+z)/nuisance.oneplusz0_ia,nuisance.eta_ia);
  }
}

double P_II (double k, double a){
  return pow(A_IA_Joachimi(a),2.0)*pow(cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a),2.0)*Pdelta(k,a); //Joachimi+ 11 Eq.(23) + Eq. (B.6)
}
double P_dI (double k, double a){
  return A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*Pdelta(k,a); // Joachimi+ 11 Eq.(6)
}


/*=========================================================*/
/*============= Intrinsic Alignment II-Term ==============*/
/*=========================================================*/


double int_for_C_II(double a, void *params)
{
  double res, ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
   ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;

  res= W_source(a,ar[0])*W_source(a,ar[1])*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_II(k,a);
  return res;
}
     
double C_II_nointerp(double s, int ni, int nj) 
{
  if(ni!=nj && redshift.shear_photoz==0) return 0.0;
  if (abs(ni-nj) >1 && redshift.shear_photoz==3) return 0.0;

  double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return int_gsl_integrate_low_precision(int_for_C_II,(void*)array,fmax(amin_source(ni),amin_source(nj)),fmin(amax_source_IA(ni),amax_source_IA(nj)),NULL,1000); //restrict integral to redshift overlap of redshift bins
  //use low precision to avoid slow down for ni != nj with redshift.shear_photoz ==3
}



/*=========================================================*/
/*============= gI routines ==============*/
/*=========================================================*/

double int_for_C_gI(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell,fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  ell = ar[2]+0.5;
  fK = f_K(chi(a));
  k = ell/fK;
  res = W_gal(a,ar[0])*W_source(a,ar[1])*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_dI(k,a);
  return res;
}

double C_gI_nointerp(double s, int ni, int nj)
{
  if (amax_source_IA(nj)< amin_lens(nj)) return 0.; // source bin and lens bin don't overlap, no gI term
  double array[3] = {(double) ni, (double) nj,s};
  return -int_gsl_integrate_low_precision(int_for_C_gI,(void*)array,fmax(amin_source(nj),amin_lens(ni)),fmin(amax_source_IA(nj),amax_lens(ni)),NULL,1000);  //restrict integral to redshift overlap of redshift bins
  //use low precision to avoid slowdown for ni != nj with redshift.shear_photoz ==3
}


/*=========================================================*/
/*============= GI routines ==============*/
/*=========================================================*/

double int_for_C_GI(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell,fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  ell = ar[2]+0.5;
  fK = f_K(chi(a));
  k = ell/fK;
  res = (W_source(a,ar[0])*W_kappa(a,fK,ar[1])+W_source(a,ar[1])*W_kappa(a,fK,ar[0]))*dchi_da(a)/fK/fK;
  if (res > 0) res *= P_dI(k,a);
  return res;
}

double C_GI_nointerp(double s, int ni, int nj) 
{
  double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return -int_gsl_integrate_low_precision(int_for_C_GI,(void*)array,amin_source(j),amax_source_IA(k),NULL,1000);
}


/******** faster shear/ggl + IA routines *********/

double int_for_C_ggl_IA(double a, void *params)
{
  double res, ell, fK, k,norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a);
  res= W_gal(a,ar[0])*(W_kappa(a,fK,ar[1])-W_source(a,ar[1])*norm);
  
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}


double int_for_C_shear_shear_IA(double a, void *params)
{
  double res, ell, fK, k,ws1,ws2,wk1,wk2, norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  ws1 = W_source(a,ar[0]);
  ws2 = W_source(a,ar[1]);
  wk1 = W_kappa(a,fK,ar[0]);
  wk2 = W_kappa(a,fK,ar[1]);
  norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a);
  res= ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm+wk1*wk2;
  
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}


double C_shear_shear_IA(double s, int ni, int nj)
{
 double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  return int_gsl_integrate_medium_precision(int_for_C_shear_shear_IA,(void*)array,amin_source(j),amax_source(k),NULL,1000);
} 

double C_ggl_IA(double s, int nl, int ns)
{
  double array[3] = {(double) nl, (double) ns,s};
  return int_gsl_integrate_medium_precision(int_for_C_ggl_IA,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);
}
