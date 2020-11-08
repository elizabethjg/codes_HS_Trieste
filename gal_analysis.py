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
path_plots = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/plots/gal_plots/'

os.system('mkdir '+path_plots)

def ninbin(binnumber):
    
    N = np.array([])
    for n in np.arange(binnumber.min(),binnumber.max()+1):
        N = np.append(N,sum(binnumber==n))
    
    return N

# LOAD CATALOGS

gral  = np.loadtxt(path+'gral_091.dat').T
gal   = np.loadtxt(path+'glxs_091.dat').T
dm    = np.loadtxt(path+'dm_091.dat').T

def make_plots(radio,masa_cut,mass):

    # COMPUTE PARAMETERS
    if mass == 1000:
        lMass = np.log10(gral[7])
    elif mass == 500:
        lMass = np.log10(gral[8])
    elif mass == 200:
        lMass = np.log10(gral[9])
    elif mass == 30:
        lMass = np.log10(gral[10])
    elif mass == 50:
        lMass = np.log10(gral[11])
    elif mass == 100:
        lMass = np.log10(gral[12])
    
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
    
    def data_gal(radio,masa_cut=True):
    
        if radio == 1000:
            ind = 2
        elif radio == 500:
            ind = 2+19*2
        elif radio == 200:
            ind = 2+19*2*2
    
        if masa_cut:
            ind += 19
    
        N  = gal[ind]
        nx = gal[ind+1]
        ny = gal[ind+2]
        nz = gal[ind+3]
        a  = gal[ind+4]
        b  = gal[ind+5]
        c  = gal[ind+6]
        qx = gal[ind+7]
        qy = gal[ind+8]
        qz = gal[ind+9]
        
        T = (a**2 - b**2)/(a**2 - c**2)
        S = c/a
        q = np.concatenate((qx,qy,qz))
        n = np.concatenate((nx,ny,nz))
        
        return N,n,T,S,q


    N,n,T,S,q = data_gal(radio,masa_cut)
    
    
    plt.figure()
    plt.scatter(lMass,T,c=N,alpha=0.5)
    plt.plot(lMass_b,T_b,'C0')
    plt.plot(lMass_b,T_b+sT_b,'C0--')
    plt.plot(lMass_b,T_b-sT_b,'C0--')
    plt.xlabel('$\log M_{200}$')
    plt.ylabel('$T$')
    plt.ylim([0,1])
    plt.savefig(path_plots+'T_ratio_'+str(mass)+'_'+str(radio)+'_'+str(masa_cut)+'.png')
    
    plt.figure()
    plt.scatter(lMass,S,c=N,alpha=0.5)
    plt.plot(lMass_b,S_b,'C1')
    plt.plot(lMass_b,S_b+sS_b,'C1--')
    plt.plot(lMass_b,S_b-sS_b,'C1--')
    plt.xlabel('$\log M_{200}$')
    plt.ylabel('$S$')
    plt.ylim([0,1])
    plt.savefig(path_plots+'S_ratio_'+str(mass)+'_'+str(radio)+'_'+str(masa_cut)+'.png')
    
    plt.figure()
    plt.scatter(lMassp,q,c=n,alpha=0.5)
    plt.plot(lMass_bp,q_b,'C2')
    plt.plot(lMass_bp,q_b+sq_b,'C2--')
    plt.plot(lMass_bp,q_b-sq_b,'C2--')
    plt.xlabel('$\log Mass$')
    plt.ylabel('$q$')
    plt.ylim([0,1])
    plt.savefig(path_plots+'q_ratio_'+str(mass)+'_'+str(radio)+'_'+str(masa_cut)+'.png')


for j in [200,500,1000]:
    # for i in [30,50,100,200,500,1000]:
        make_plots(j,True,200)
        make_plots(j,False,200)
