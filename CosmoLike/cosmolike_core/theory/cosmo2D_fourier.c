//#include "pt.c"

double W_kappa(double a, double fK, double nz);//complete lens efficiency weight
double W_source(double a, double nz); //source redshift distribution (radial weight for IA,source clustering)
double W_gal(double a, double nz); //complete weight for galaxy statistics

double C_cl_tomo(double l, int ni, int nj);  //galaxy clustering power spectrum of galaxies in bins ni, nj
double C_cl_tomo_nointerp(double l, int ni, int nj);

double C_gl_tomo(double l, int ni, int nj);  //G-G lensing power spectrum, lens bin ni, source bin nj
double C_gl_tomo_nointerp(double l, int ni, int nj);

double C_shear_tomo(double l, int ni, int nj); //shear tomography power spectra
double C_shear_tomo_nointerp(double l, int ni, int nj);
/**********************************************************/

double MG_Sigma(double a)
{
//  double aa=a*a;
//  double omegam=cosmology.Omega_m/(aa*a);
  double omegav=omv_vareos(a);
  double hub=hoverh0(a);
  hub = hub*hub;
 
  return cosmology.MGSigma*omegav/hub/cosmology.Omega_v;
}

double dchi_da(double a){
  return 1./(a*a*hoverh0(a));
}
double W_kappa(double a, double fK, double nz){
  double wkappa = 1.5*cosmology.Omega_m*fK/a*g_tomo(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wkappa *= (1.+MG_Sigma(a));
  }
  return wkappa;
}


double W_gal(double a, double nz){
  double wgal = gbias.b1_function(1./a-1.,(int)nz)*pf_photoz(1./a-1.,(int)nz)*hoverh0(a);
  double wmag = gbias.b_mag[(int)nz]*1.5*cosmology.Omega_m*f_K(chi(a))/a*g_lens(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wmag *= (1.+MG_Sigma(a));
  }
  return wgal + wmag;
  // return wgal;
}

double W_source(double a, double nz){
  return zdistr_photoz(1./a-1.,(int)nz)*hoverh0(a);
}

double W_mag(double a, double fK, double nz){
  double wmag = 1.5*cosmology.Omega_m*fK/a*g_lens(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wmag *= (1.+MG_Sigma(a));
  }
  return wmag;
}
/*********** Limber approximation integrands for angular power spectra *************/
// use W_i W_j/fK^2 dchi_da(a) as Limber weight for power spectra


double int_for_C_cl_tomo(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res=W_gal(a,ar[0])*W_gal(a,ar[1])*dchi_da(a)/fK/fK;
  res= res*Pdelta(k,a);
  return res;
}


double int_for_C_gl_tomo(double a, void *params) // Add RSD
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_C_gl_tomo");

  double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if(ell_prefactor2<=0.) 
    ell_prefactor2=0.;
  else
    ell_prefactor2=sqrt(ell_prefactor2);

  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  double chi_0,chi_1,a_0,a_1;
  chi_0 = f_K(ell/k);
  chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;}
  a_0 = a_chi(chi_0);
  a_1 = a_chi(chi_1);

  double wgal = W_gal(a,ar[0]);
  wgal += W_mag(a,fK,ar[0])*(ell_prefactor1/ell/ell -1.) ;
  
  res= (wgal)*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK  * ell_prefactor2/ell/ell;
  res= res*Pdelta(k,a);
  return res;
}



double int_for_C_shear_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_C_shear_tomo");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= W_kappa(a,fK,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
  res= res*Pdelta(k,a); 
  return res;
}

/*********** angular power spectra - without look-up tables ******************/

double C_cl_tomo_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{ static int init =-1;
  double array[3] = {1.0*ni,1.0*nj,l};

  return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
}


double C_gl_tomo_nointerp(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  // l=1.0*2;
  if(l==1.) {return 0.;}
  double array[3] = {(double)ni,(double)nj,l};
  double res = int_gsl_integrate_medium_precision(int_for_C_gl_tomo,(void*)array,amin_lens(ni),0.99999,NULL,1000);
  return res;
}

double C_shear_tomo_nointerp(double l, int ni, int nj) //shear tomography power spectra of source galaxy bins ni, nj
{
  double res;
  double array[3] = {(double)ni, (double)nj,l};
  // printf("ni,nj:%d,%d\n", ni,nj );
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  res=int_gsl_integrate_medium_precision(int_for_C_shear_tomo,(void*)array,amin_source(j),amax_source(k),NULL,1000);
  // printf("res:%lg\n", res);
  return res;
  
}

