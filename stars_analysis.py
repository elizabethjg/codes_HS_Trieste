import sys
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from scipy import stats
import os
from pylab import *
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

path = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/catalog/'
path_plots = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/plots/stars_plots/'

os.system('mkdir '+path_plots)

def ninbin(binnumber):
    
    N = np.array([])
    for n in np.arange(binnumber.min(),binnumber.max()+1):
        N = np.append(N,sum(binnumber==n))
    
    return N

# LOAD CATALOGS

gral  = np.loadtxt(path+'gral_091.dat').T
stars = np.loadtxt(path+'stars_091.dat').T
dm    = np.loadtxt(path+'dm_091.dat').T

# COMPUTE PARAMETERS

def make_plots(radio,mass):

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

    if radio == 1000:
        lab2 = '1000'
    elif radio == 500:
        lab2 = '500'
    elif radio == 200:
        lab2 = '200'
    elif radio == 30:
        lab2 = 'BCG-30kpc'
    elif radio == 50:
        lab2 = 'BCG-50kpc'
    elif radio == 100:
        lab2 = 'BCG-0.1R500'


    #-----------------
    # DM-only
    
    a_dm = dm[2]
    b_dm = dm[3]
    c_dm = dm[4]
    qx_dm = dm[5]
    qy_dm = dm[6]
    qz_dm = dm[7]
    
    T_dm = (a_dm**2 - b_dm**2)/(a_dm**2 - c_dm**2)
    S_dm = c_dm/a_dm
    
    q_dm = np.concatenate((qx_dm,qy_dm,qz_dm))
    lMassp = np.array((lMass.tolist())*3)
    
    # T bins
    bined = stats.binned_statistic(lMass,T_dm,statistic='median', bins=10)
    lMass_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    T_b     = bined.statistic
    bined = stats.binned_statistic(lMass,T_dm,statistic='std', bins=10)
    sT_b     = bined.statistic/np.sqrt(ninbin(bined.binnumber))
    # S bins 
    bined = stats.binned_statistic(lMass,S_dm,statistic='median', bins=10)
    S_b     = bined.statistic
    bined = stats.binned_statistic(lMass,S_dm,statistic='std', bins=10)
    sS_b     = bined.statistic/np.sqrt(ninbin(bined.binnumber))
    # q bins 
    bined = stats.binned_statistic(lMassp,q_dm,statistic='median', bins=10)
    lMass_bp = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    q_b     = bined.statistic
    bined = stats.binned_statistic(lMassp,q_dm,statistic='std', bins=10)
    sq_b     = bined.statistic/np.sqrt(ninbin(bined.binnumber))
    
    #-----------------
    # GAL DIST
    
    def data_gal(radio):
    
        if radio == 30:
            ind = 2
        elif radio == 50:
            ind = 2+15
        elif radio == 100:
            ind = 2+15*2
        elif radio == 1000:
            ind = 2+15*3
        elif radio == 500:
            ind = 2+15*4
        elif radio == 200:
            ind = 2+15*5
    
        a  = stars[ind+0]
        b  = stars[ind+1]
        c  = stars[ind+2]
        qx = stars[ind+3]
        qy = stars[ind+4]
        qz = stars[ind+5]
                    
        T = (a**2 - b**2)/(a**2 - c**2)
        S = c/a
        q = np.concatenate((qx,qy,qz))
        
        
        return T,S,q



    T,S,q = data_gal(radio)

    # T bins
    bined = stats.binned_statistic(lMass,T,statistic='median', bins=10)
    T_bs     = bined.statistic
    bined = stats.binned_statistic(lMass,T,statistic='std', bins=10)
    sT_bs     = bined.statistic/np.sqrt(ninbin(bined.binnumber))

    
    plt.figure()
    plt.plot(lMass,T,'C0o',alpha=0.5)

    plt.plot(lMass_b,T_b,'k')
    plt.plot(lMass_b,T_b+sT_b,'k--')
    plt.plot(lMass_b,T_b-sT_b,'k--')
    
    plt.plot(lMass_b,T_bs,'C0')
    plt.plot(lMass_b,T_bs+sT_bs,'C0--')
    plt.plot(lMass_b,T_bs-sT_bs,'C0--')
    
    plt.xlabel('$\log M_{'+lab+'}$')
    plt.ylabel('$T_{'+lab2+'}$')
    plt.ylim([0,1])
    plt.savefig(path_plots+'T_s'+lab2+'_'+lab+'.png')

    # S bins
    bined = stats.binned_statistic(lMass,S,statistic='median', bins=10)
    S_bs     = bined.statistic
    bined = stats.binned_statistic(lMass,S,statistic='std', bins=10)
    sS_bs     = bined.statistic/np.sqrt(ninbin(bined.binnumber))
    
    plt.figure()
    plt.plot(lMass,S,'C1o',alpha=0.5)
    plt.plot(lMass_b,S_b,'k')
    plt.plot(lMass_b,S_b+sS_b,'k--')
    plt.plot(lMass_b,S_b-sS_b,'k--')
    plt.plot(lMass_b,S_bs,'C1')
    plt.plot(lMass_b,S_bs+sS_bs,'C1--')
    plt.plot(lMass_b,S_bs-sS_bs,'C1--')
    plt.xlabel('$\log M_{'+lab+'}$')
    plt.ylabel('$S_{'+lab2+'}$')
    plt.ylim([0,1])
    plt.savefig(path_plots+'S_s'+lab2+'_'+lab+'.png')

    # S bins
    bined = stats.binned_statistic(lMassp,q,statistic='median', bins=10)
    q_bs     = bined.statistic
    bined = stats.binned_statistic(lMassp,q,statistic='std', bins=10)
    sq_bs     = bined.statistic/np.sqrt(ninbin(bined.binnumber))

    
    plt.figure()
    plt.plot(lMassp,q,'C2o',alpha=0.5)
    plt.plot(lMass_bp,q_b,'k')
    plt.plot(lMass_bp,q_b+sq_b,'k--')
    plt.plot(lMass_bp,q_b-sq_b,'k--')
    plt.plot(lMass_bp,q_bs,'C2')
    plt.plot(lMass_bp,q_bs+sq_bs,'C2--')
    plt.plot(lMass_bp,q_bs-sq_bs,'C2--')
    plt.xlabel('$\log M_{'+lab+'}$')
    plt.ylabel('$q_{'+lab2+'}$')
    plt.ylim([0,1])
    plt.savefig(path_plots+'q_s'+lab2+'_M'+lab+'.png')


for j in [30,50,100,1000,500,200]:
    # for i in [30,50,100,200,500,1000]:
    for i in [200]:
        make_plots(j,i)
