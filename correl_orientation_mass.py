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

gral  = np.loadtxt('../catalog/gral_091_2.dat').T
path_plots = '../plots/correl_orientation_mass/'

lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)

dm = DarkMatter(200)

# CON GALAXIAS

gx1000 = Galaxias(radio = 1000)
gx500 = Galaxias(radio = 500)
gx200 = Galaxias(radio = 200)
gx_c  = Galaxias('concentradas')
gx_e  = Galaxias('extendidas')

mN =  np.isfinite(gx1000.S)
mN2D = mN.tolist()*3

ct1000_3D, t1000_3D    = cosangle(dm.a3D[mN],gx1000.a3D[mN]) 
ct500_3D , t500_3D     = cosangle(dm.a3D,gx500.a3D) 
ct200_3D , t200_3D     = cosangle(dm.a3D,gx200.a3D) 
ctc_3D , tc_3D     = cosangle(dm.a3D,gx_c.a3D) 
cte_3D , te_3D     = cosangle(dm.a3D,gx_e.a3D) 

ct1000_2D, t1000_2D    = cosangle(dm.a2D[mN2D],gx1000.a2D[mN2D]) 
ct500_2D , t500_2D     = cosangle(dm.a2D,gx500.a2D) 
ct200_2D , t200_2D     = cosangle(dm.a2D,gx200.a2D) 
ctc_2D , tc_2D     = cosangle(dm.a2D,gx_c.a2D) 
cte_2D , te_2D     = cosangle(dm.a2D,gx_e.a2D) 


plt.figure()

plot_binned(lM[mN],t1000_3D,'R1000','C0',':',nbins=5)
plot_binned(lM,t500_3D,'R500','C0','--',nbins=5)
plot_binned(lM,t200_3D,'R200','C0','-',nbins=5)

plot_binned(lM,tc_3D,'concentradas','C3','-',nbins=5)
plot_binned(lM,te_3D,'extendidas','C1','-',nbins=5)

plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{3D}$')
plt.ylim([0,50])
plt.savefig(path_plots+'gxs_theta3D_200.png')


plt.figure()

plot_binned(lMp[mN2D],t1000_2D,'R1000','C0',':',nbins=5)
plot_binned(lMp,t500_2D,'R500','C0','--',nbins=5)
plot_binned(lMp,t200_2D,'R200','C0','-',nbins=5)

plot_binned(lMp,tc_2D,'concentradas','C3','-',nbins=5)
plot_binned(lMp,te_2D,'extendidas','C1','-',nbins=5)

plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{2D}$')
plt.ylim([0,50])
plt.savefig(path_plots+'gxs_theta2D_200.png')

# CON ESTRELLS

st1000  = Stars(radio = 1000)
st500   = Stars(radio = 500)
st200   = Stars(radio = 200)
st30    = Stars(radio = 30)
st50    = Stars(radio = 50)
st100   = Stars(radio = 100)

ct1000_3D, t1000_3D  = cosangle(dm.a3D,st1000.a3D) 
ct500_3D , t500_3D   = cosangle(dm.a3D,st500.a3D) 
ct200_3D , t200_3D   = cosangle(dm.a3D,st200.a3D) 
ct100_3D, t100_3D    = cosangle(dm.a3D,st100.a3D) 
ct50_3D , t50_3D     = cosangle(dm.a3D,st50.a3D) 
ct30_3D , t30_3D     = cosangle(dm.a3D,st30.a3D) 

ct1000_2D, t1000_2D  = cosangle(dm.a2D,st1000.a2D) 
ct500_2D , t500_2D   = cosangle(dm.a2D,st500.a2D) 
ct200_2D , t200_2D   = cosangle(dm.a2D,st200.a2D) 
ct100_2D, t100_2D    = cosangle(dm.a2D,st100.a2D)
ct50_2D , t50_2D     = cosangle(dm.a2D,st50.a2D) 
ct30_2D , t30_2D     = cosangle(dm.a2D,st30.a2D) 

plt.figure()

plot_binned(lM,t1000_3D,'R1000','C1',':',nbins=5)
plot_binned(lM,t500_3D,'R500','C1','--',nbins=5)
plot_binned(lM,t200_3D,'R200','C1','-',nbins=5)
plot_binned(lM,t100_3D,'0.1R500','C3',':',nbins=5)
plot_binned(lM,t50_3D,'R50kpc','C3','--',nbins=5)
plot_binned(lM,t30_3D,'R30kpc','C3','-',nbins=5)

plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{3D}$')
plt.ylim([0,50])
plt.savefig(path_plots+'st_theta3D_200.png')

plt.figure()

plot_binned(lMp,t1000_2D,'R1000','C1',':',nbins=5)
plot_binned(lMp,t500_2D,'R500','C1','--',nbins=5)
plot_binned(lMp,t200_2D,'R200','C1','-',nbins=5)
plot_binned(lMp,t100_2D,'0.1R500','C3',':',nbins=5)
plot_binned(lMp,t50_2D,'R50kpc','C3','--',nbins=5)
plot_binned(lMp,t30_2D,'R30kpc','C3','-',nbins=5)

plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{2D}$')
plt.ylim([0,50])
plt.savefig(path_plots+'st_theta2D_200.png')
