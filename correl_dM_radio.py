import sys
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from scipy import stats
from pylab import *
from main import *
from scipy.stats import pearsonr
import os
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

cDMDM = CorrelR(DarkMatter,DarkMatter)
cSTDM = CorrelR(Stars,DarkMatter)
cGXDM = CorrelR(Galaxias,DarkMatter)

plotR(cSTDM.R,cSTDM.Sr,r'$S$','S_stars_dm')
plotR(cSTDM.R,cSTDM.Tr,r'$T$','T_stars_dm')
plotR(cSTDM.R,cSTDM.t3D,r'$\theta_{3D}$','t3D_stars_dm')

plotR(cSTDM.Rp,cSTDM.qr,r'$q$','q_stars_dm')
plotR(cSTDM.Rp,cSTDM.t2D,r'$\theta_{2D}$','t2D_stars_dm')

plotR(cSTDM.R,cSTDM.t3D_dm200,r'$\theta_{3D}$','t3D_stars_dm200')
plotR(cSTDM.Rp,cSTDM.t2D_dm200,r'$\theta_{2D}$','t2D_stars_dm200')

plotR(cSTDM.R,cSTDM.t3D_dm1000,r'$\theta_{3D}$','t3D_stars_dm1000')
plotR(cSTDM.Rp,cSTDM.t2D_dm1000,r'$\theta_{2D}$','t2D_stars_dm1000')

plotR(cGXDM.R,cGXDM.t3D_dm200,r'$\theta_{3D}$','t3D_gxs_dm200',gxs=True)
plotR(cGXDM.Rp,cGXDM.t2D_dm200,r'$\theta_{2D}$','t2D_gxs_dm200',gxs=True)
       
plotR(cGXDM.R,cGXDM.t3D_dm1000,r'$\theta_{3D}$','t3D_gxs_dm1000',gxs=True)
plotR(cGXDM.Rp,cGXDM.t2D_dm1000,r'$\theta_{2D}$','t2D_gxs_dm1000',gxs=True)

mgx = np.array([0,0,0,1,1,1]).astype(bool)
Rlegend = np.array(['30kpc','50kpc','0.1R500','R1000','R500','R200'])
path_plots = '../plots/correl_dM_radio/'

######## S(R) ########
plt.figure()
plotR_ind(cDMDM.R,cDMDM.Sr,'k','DM')
plotR_ind(cSTDM.R,cSTDM.Sr,'C0','stars')
plotR_ind(cGXDM.R[:,mgx],cGXDM.Sr[:,mgx],'C1','gxs')

plt.legend()
plt.ylabel('S')
plt.xticks(np.median(np.log10(cDMDM.R),axis=0),Rlegend)
plt.axis([np.log10(30),np.log10(2000),0.4,1])
plt.savefig(path_plots+'S_R.png')

######## T(R) ########
plt.figure()
plotR_ind(cDMDM.R,cDMDM.Tr,'k','DM')
plotR_ind(cSTDM.R,cSTDM.Tr,'C0','stars')
plotR_ind(cGXDM.R[:,mgx],cGXDM.Tr[:,mgx],'C1','gxs')

plt.legend()
plt.ylabel('T')
plt.xticks(np.median(np.log10(cDMDM.R),axis=0),Rlegend)
plt.axis([np.log10(30),np.log10(2000),0.4,1])
plt.savefig(path_plots+'T_R.png')

######## q(R) ########
plt.figure()
plotR_ind(cDMDM.Rp,cDMDM.qr,'k','DM')
plotR_ind(cSTDM.Rp,cSTDM.qr,'C0','stars')
plotR_ind(cGXDM.Rp[:,mgx],cGXDM.qr[:,mgx],'C1','gxs')

plt.legend()
plt.ylabel('q')
plt.xticks(np.median(np.log10(cDMDM.Rp),axis=0),Rlegend)
plt.axis([np.log10(30),np.log10(2000),0.4,1])
plt.savefig(path_plots+'q_R.png')
