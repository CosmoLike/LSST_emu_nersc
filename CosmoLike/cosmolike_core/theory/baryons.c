#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>


#include "../../theory/basics.c"
#include "../../emu13/emu.c"
#include "../../theory/parameters.c"
#include "../../theory/cosmo3D.c"
#include "../../theory/redshift.c"

typedef struct {
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  gsl_spline2d *spline;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
  double *grid_values;
} spline_interpolator;
spline_interpolator global_baryon_spline;


void init_baryon_spline_interpolation_to_TNG100(){
  
  const size_t N = 100;             /* number of points to interpolate */
  const double lnk_values[] = { 0.0, 1.0 };
  const double z_values[] = { 3.71, 3.49, 3.28, 2.90, 2.44, 2.10, 1.74, 1.41, 1.04, 0.7, 0.35, 0.18, 0.0 }; /* define unit square */
  const size_t nk = sizeof(lnk_values) / sizeof(double); /* y grid points */
  const size_t nz = sizeof(z_values) / sizeof(double); /* x grid points */
  global_baryon_spline.grid_values = malloc(nk * nz * sizeof(double));
  global_baryon_spline.spline = gsl_spline2d_alloc(global_baryon_spline.T, nk, nz);
  global_baryon_spline.xacc = gsl_interp_accel_alloc();
  global_baryon_spline.yacc = gsl_interp_accel_alloc();
  size_t i, j;
  
  for(int i = 0; i < nk; i++){
    for(int j = 0; j < nz; j++){
      gsl_spline2d_set(spline, global_baryon_spline.grid_values, i, j, TNG100[i][j]);
    }
  }

  /* initialize interpolation */
  gsl_spline2d_init(global_baryon_spline.spline, lnk_values, z_values, Pk_ratio, nk, nz);
  

}


void free_baryon_spline_interpolation_to_TNG100(){
  gsl_spline2d_free(global_baryon_spline.spline);
  gsl_interp_accel_free(global_baryon_spline.xacc);
  gsl_interp_accel_free(global_baryon_spline.yacc);
  free(global_baryon_spline.grid_values);
}


double fraction_Pdelta_baryon_from_sims_DES(double k,double z)
{
  return gsl_spline2d_eval(global_baryon_spline.spline, k, z, global_baryon_spline.xacc, global_baryon_spline.yacc);
}


double fraction_Pdelta_baryon_from_spline_DES(double k,double z)
{
  int i;
  double val;
  static int READ_TABLE=0;
  static double **table;
  FILE *F;
  int baryon_zbin=13;
  int baryon_kbin=325;
  double logk_min = -3.3010299956639813;
  double logk_max = 3.1760912590556813;
  double dz=(redshift.shear_zdistrpar_zmax)/(10.0);
  static double dlogk=(logk_max-logk_min)/(double(baryon_kbin-1));
  val = interpol2d(AGN_DESdepth, baryon_kbin, logk_min, logk_max, dk, log(k), baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0);
  return val; 
}

//   if (READ_TABLE==0){
//     if (strcmp(pdeltaparams.baryons,"AGN_DESdepth")==0) P_type = 0;
//     if (strcmp(pdeltaparams.baryons,"NOSN_DESdepth")==0) P_type = 1;
//     if (strcmp(pdeltaparams.baryons,"NOSN_NOZCOOL_DESdepth")==0) P_type = 2;
//     if (strcmp(pdeltaparams.baryons,"NOZCOOL_DESdepth")==0) P_type = 3;
//     if (strcmp(pdeltaparams.baryons,"REF_DESdepth")==0) P_type = 4;
//     if (strcmp(pdeltaparams.baryons,"WDENS_DESdepth") ==0) P_type = 5;
//     if (strcmp(pdeltaparams.baryons,"DBLIMFV1618_DESdepth")  ==0) P_type = 6;
//     if (strcmp(pdeltaparams.baryons,"WML4_DESdepth")==0) P_type = 7;
//     if (strcmp(pdeltaparams.baryons,"WML1V848_DESdepth")==0) P_type = 8;
//     if (strcmp(pdeltaparams.baryons,"AD_DESdepth")==0) P_type = 9;
//     if (strcmp(pdeltaparams.baryons,"CX_DESdepth")==0) P_type = 10;
//     if (strcmp(pdeltaparams.baryons,"CW_DESdepth")==0) P_type = 11;
//     if (strcmp(pdeltaparams.baryons,"A_DESdepth")==0) P_type = 12;
//     if (strcmp(pdeltaparams.baryons,"CSF_DESdepth")==0) P_type = 13;
//   }
//   switch (P_type){
//     case 0:  val = interpol2d(AGN_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 1:  val = interpol2d(NOSN_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 2:  val = interpol2d(NOSN_NOZCOOL_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 3:  val = interpol2d(NOZCOOL_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 4:  val = interpol2d(REF_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 5:  val = interpol2d(WDENS_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 6:  val = interpol2d(DBLIMFV1618_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 7:  val = interpol2d(WML4_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 8:  val = interpol2d(WML1V848_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 9:  val = interpol2d(AD_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 10:  val = interpol2d(CX_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 11:  val = interpol2d(CW_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 12:  val = interpol2d(A_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;
//     case 13:  val = interpol2d(CSF_DESdepth, baryon_kbin, 0.3, 10.0, dk, kintern, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0); break;

//     default: 
//             printf("baryons.c: %s baryonic scenario not defined\n",pdeltaparams.baryons);

//             break;
//     }
//     kintern=k/cosmology.coverH0;
 
//  return val; 
// }







//    if (strcmp(pdeltaparams.baryons,"AGN_LSSTdepth")==0) P_type = 14;
//     if (strcmp(pdeltaparams.baryons,"NOSN_LSSTdepth")==0) P_type = 15;
//     if (strcmp(pdeltaparams.baryons,"NOSN_NOZCOOL_LSSTdepth")==0) P_type = 16;
//     if (strcmp(pdeltaparams.baryons,"NOZCOOL_LSSTdepth")==0) P_type = 17;
//     if (strcmp(pdeltaparams.baryons,"REF_LSSTdepth")==0) P_type = 18;
//     if (strcmp(pdeltaparams.baryons,"WDENS_LSSTdepth") ==0) P_type = 19;
//     if (strcmp(pdeltaparams.baryons,"DBLIMFV1618_LSSTdepth")  ==0) P_type = 20;
//     if (strcmp(pdeltaparams.baryons,"WML4_LSSTdepth")==0) P_type = 21;
//     if (strcmp(pdeltaparams.baryons,"WML1V848_LSSTdepth")==0) P_type = 22;
//     if (strcmp(pdeltaparams.baryons,"AD_LSSTdepth")==0) P_type = 23;
//     if (strcmp(pdeltaparams.baryons,"CX_LSSTdepth")==0) P_type = 24;
//     if (strcmp(pdeltaparams.baryons,"CW_LSSTdepth")==0) P_type = 25;
//   if (table!=0) free_double_matrix(table, 0, OWLS_kbin-1, 0, OWLS_zbin-1);
//     table   = create_double_matrix(0, OWLS_kbin-1, 0, OWLS_zbin-1);
//     printf("%s\n",file.DATA_FILE);  
//     F=fopen(file.DATA_FILE,"r");
//       for(i=0;i<OWLS_kbin;i++){
// 	fscanf(F,"%le %le %le %le %le %le %le %le %le %le %le\n",&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9,&a10,&a11);
// //	printf("%le %le %le %le %le %le %le %le %le %le %le\n",a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11);
// 	table[i][0]=a1;
// 	table[i][1]=a2;
// 	table[i][2]=a3;
// 	table[i][3]=a4;
// 	table[i][4]=a5;
// 	table[i][5]=a6;
// 	table[i][6]=a7;
// 	table[i][7]=a8;
// 	table[i][8]=a9;
// 	table[i][9]=a10;
// 	table[i][10]=a11;
//       }
//       READ_TABLE=1;
//   }
//  kintern=k/cosmology.coverH0;
//  val = interpol2d(table, OWLS_kbin, 0.3, 10.0, dk, kintern,OWLS_zbin, 0.0, 2.0, dz, z, 0.0, 0.0);
//  return val; 
// }



