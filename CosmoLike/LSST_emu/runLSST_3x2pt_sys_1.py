#!/usr/bin/env python

import sys
sys.path.append('/home/u17/timeifler/WFIRST_LSST')

from cosmolike_libs_WF_LSST import * 
from schwimmbad import MPIPool

inv=["inv_3x2pt_WFIRST_area=2.000000e+03_dmo","inv_3x2pt_LSST_Y10_area=1.800000e+04_dmo","inv_3x2pt_WFIRST_area=1.800000e+04_dmo"]

data=["3x2pt_WFIRST_area=2.000000e+03_dmo","3x2pt_LSST_Y10_area=1.800000e+04_dmo","3x2pt_WFIRST_area=1.800000e+04_dmo"]

bary=["3x2pt_WFIRST_area=2.000000e+03","3x2pt_LSST_Y10_area=1.800000e+04",
"3x2pt_WFIRST_area=1.800000e+04"]

source_z=["zdistri_WFIRST_LSST_lensing_fine_bin","SRD_zdistri_model_z0=1.100000e-01_beta=6.800000e-01_Y10_source.txt","zdistri_WFIRST_LSST_lensing_fine_bin"] 

lens_z=["zdistri_WFIRST_LSST_clustering_fine_bin","SRD_zdistri_model_z0=2.800000e-01_beta=9.000000e-01_Y10_lens.txt","zdistri_WFIRST_LSST_clustering_fine_bin"]

shear_prior=[0.002,0.003,0.002] 
delta_z_prior_shear=[0.001,0.001,0.001]
delta_z_prior_clustering=[0.001,0.001,0.001]
sigma_z_shear=[0.01,0.05,0.02]
sigma_z_clustering=[0.01,0.03,0.02]
sigma_z_prior_shear=[0.002,0.003,0.002]
sigma_z_prior_clustering=[0.002,0.003,0.002]

nsource_table=[51.0,27.0,48.0]  
nlens_table=[66.0,48.0,52.0]
area_table=[2000.0,18000.0,18000.0]

survey_designation=["WFIRST","LSST_Y10","WFIRST"]
tomo_binning_source=["source_std","source_std","source_std"]
tomo_binning_lens=["WF_SN10","LSST_gold","WF_SN10"]

model=1
file_source_z = os.path.join(dirname, "zdistris/",source_z[model])
file_lens_z = os.path.join(dirname, "zdistris/",lens_z[model])
data_file = os.path.join(dirname, "datav/",data[model])
cov_file = os.path.join(dirname, "inv/",inv[model])
#cov_file = os.path.join("/Users/timeifler/Dropbox/cosmolike_store/WFIRST_LSST/inv/",inv[model])
chain_file = os.path.join(dirname, "like/like_"+data[model])
bary_file=os.path.join(dirname, "baryons/",bary[model])

initcosmo("halofit")
initbins(15,20.0,3000.0,3000.0,21.0,10,10)
initpriors(shear_prior[model],sigma_z_shear[model],delta_z_prior_shear[model],sigma_z_prior_shear[model],sigma_z_clustering[model],delta_z_prior_clustering[model],sigma_z_prior_clustering[model],3.0,1.2,3.8,2.0,16.0,5.0,0.8);
  
initsurvey(survey_designation[model],nsource_table[model],nlens_table[model],area_table[model])
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian",tomo_binning_source[model],tomo_binning_lens[model])
initclusters()
initia("NLA_HF","GAMA")

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("3x2pt")
initdatainv(cov_file ,data_file, bary_file)

#sample_params=sample_bias_only(get_N_tomo_clustering())
#sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_shear_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_photo_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_bias(get_N_tomo_clustering())
sample_params = sample_cosmology_2pt_nuisance_IA_bary_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

sample_main(sample_params,sigma_z_shear[model],sigma_z_clustering[model],8000,1120,1,chain_file, blind=False, pool=MPIPool())

