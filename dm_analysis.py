import sys
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from scipy import stats
from pylab import *
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

path = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/catalog/'
path_plots = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/plots/dm_plots/'


# LOAD CATALOGS

gral  = np.loadtxt(path+'gral_091.dat').T
dm    = np.loadtxt(path+'dm_091.dat').T

# COMPUTE PARAMETERS

def make_plots(mass):

    # COMPUTE PARAMETERS
    if mass == 1000:
        lMass = np.log10(gral[7])
        lab = '1000'
    elif mass == 500:
        lMass = np.log10(gral[8])
        lab = '500'
    elif mass == 200:
        lMass = np.log10(gral[9])
        lab = '200'
    elif mass == 30:
        lMass = np.log10(gral[10])
        lab = 'BCG-30kpc'
    elif mass == 50:
        lMass = np.log10(gral[11])
        lab = 'BCG-50kpc'
    elif mass == 100:
        lMass = np.log10(gral[12])
        lab = 'BCG-0.1R500'
    
    a_dm = dm[2]
    b_dm = dm[3]
    c_dm = dm[4]
    qx_dm = dm[5]
    qy_dm = dm[6]
    qz_dm = dm[7]
    
    T_dm = (a_dm**2 - b_dm**2)/(a_dm**2 - c_dm**2)
    S_dm = c_dm/a_dm
    
    q_dm = np.concatenate((qx_dm,qy_dm,qz_dm))
    
    a_dm0 = dm[17]
    b_dm0 = dm[18]
    c_dm0 = dm[19]
    
    qx_dm0 = dm[20]
    qy_dm0 = dm[21]
    qz_dm0 = dm[22]
    
    T_dm0 = (a_dm0**2 - b_dm0**2)/(a_dm0**2 - c_dm0**2)
    S_dm0 = c_dm0/a_dm0
    
    q_dm0 = np.concatenate((qx_dm0,qy_dm0,qz_dm0))
    
    # T bins plots
    
    bined = stats.binned_statistic(lMass,T_dm,statistic='median', bins=10)
    lMass_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    T_b     = bined.statistic
    
    plt.figure()
    plt.plot(lMass,T_dm,'C0o',alpha=0.5)
    plt.plot(lMass_b,T_b,'C0')
    
    bined = stats.binned_statistic(lMass,T_dm0,statistic='median', bins=10)
    lMass_b0 = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    T_b0     = bined.statistic
    
    plt.plot(lMass,T_dm0,'C0x',alpha=0.5)
    plt.plot(lMass_b0,T_b0,'C0--')
    plt.xlabel('$\log M_{'+lab+'}$')
    plt.ylabel('$T$')
    plt.ylim([0,1])
    plt.savefig(path_plots+'T_mass'+lab+'.png')
    
    # S bins plots
    
    bined = stats.binned_statistic(lMass,S_dm,statistic='median', bins=10)
    lMass_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    S_b     = bined.statistic
    
    plt.figure()
    plt.plot(lMass,S_dm,'C1o',alpha=0.5)
    plt.plot(lMass_b,S_b,'C1')
    
    bined = stats.binned_statistic(lMass,S_dm0,statistic='median', bins=10)
    lMass_b0 = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    S_b0     = bined.statistic
    
    plt.plot(lMass,S_dm0,'C1x',alpha=0.5)
    plt.plot(lMass_b0,S_b0,'C1--')
    plt.xlabel('$\log M_{'+lab+'}$')
    plt.ylabel('$S$')
    plt.ylim([0,1])
    plt.savefig(path_plots+'S_mass'+lab+'.png')
    
    # q bins plots
    
    lMass = np.array((lMass.tolist())*3)
    
    bined = stats.binned_statistic(lMass,q_dm,statistic='median', bins=10)
    lMass_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    q_b     = bined.statistic
    
    plt.figure()
    plt.plot(lMass,q_dm,'C2o',alpha=0.5)
    plt.plot(lMass_b,q_b,'C2')
    
    bined = stats.binned_statistic(lMass,q_dm0,statistic='median', bins=10)
    lMass_b0 = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    q_b0     = bined.statistic
    
    plt.plot(lMass,q_dm0,'C2x',alpha=0.5)
    plt.plot(lMass_b0,q_b0,'C2--')
    plt.xlabel('$\log M_{'+lab+'}$')
    plt.ylabel('$q$')
    plt.ylim([0,1])
    plt.savefig(path_plots+'q_mass'+lab+'.png')
    

for i in [30,50,100,200,500,1000]:
    make_plots(i)
