import numpy as np
import os
from numpy.linalg import inv
from numpy.linalg import det

list1=["like_3x2pt_LSST_Y10_area=1.800000e+04_dmo","like_3x2pt_WFIRST_area=2.000000e+03_dmo","like_3x2pt_WFIRST_area=1.800000e+04_dmo"]

FOMDE=np.array([1.0,1.0,1.0])
i=0
FOMDE[i]
print FOMDE[i]


for name in list1:
    i=0
    filen=name
    print filen
    d1 = np.genfromtxt(filen,skip_header=7800000)
    #d1.shape
    cov=np.cov(d1[:,:],rowvar=False)
    covDE=cov[np.ix_([3,4],[3,4])]
    FM=inv(cov) 
    FOMDE[i]=(np.power(det(inv(covDE)),1./2))
    print FOMDE[i]
    i=i+1
