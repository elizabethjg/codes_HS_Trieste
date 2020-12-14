import sys
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from scipy import stats
from pylab import *
from main import *
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)
from scipy.stats import pearsonr

# path = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/catalog/'
# path_plots = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/plots/stars_plots/'

path = '../catalog/'
path_plots = '../plots/stars_plots/'

gral  = np.loadtxt(path+'gral_091.dat').T
lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)


dm = DarkMatter()


# SHAPE PARAMETERS

s30   = Stars(30)
s50   = Stars(50)
s100  = Stars(100)
s200  = Stars(200)
s500  = Stars(500)
s1000 = Stars(1000)

S = np.vstack((s30.S,s50.S,s100.S,s1000.S,s500.S,s200.S)).T
T = np.vstack((s30.T,s50.T,s100.T,s1000.T,s500.T,s200.T)).T
q = np.vstack((s30.q,s50.q,s100.q,s1000.q,s500.q,s200.q)).T


# ANGLES ACCORDINT TO DM

ct30_3D  , t30_3D     = cosangle(dm.a3D,s30.a3D) 
ct50_3D  , t50_3D     = cosangle(dm.a3D,s50.a3D)  
ct100_3D , t100_3D    = cosangle(dm.a3D,s100.a3D) 
ct200_3D , t200_3D    = cosangle(dm.a3D,s200.a3D) 
ct500_3D , t500_3D    = cosangle(dm.a3D,s500.a3D)
ct1000_3D, t1000_3D   = cosangle(dm.a3D,s1000.a3D)

ct30_2D  , t30_2D     = cosangle(dm.a2D,s30.a2D) 
ct50_2D  , t50_2D     = cosangle(dm.a2D,s50.a2D)  
ct100_2D , t100_2D    = cosangle(dm.a2D,s100.a2D) 
ct200_2D , t200_2D    = cosangle(dm.a2D,s200.a2D) 
ct500_2D , t500_2D    = cosangle(dm.a2D,s500.a2D)
ct1000_2D, t1000_2D   = cosangle(dm.a2D,s1000.a2D)

t3D = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
t2D = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T


## VARIATION WITH RADIUS

R1000 = gral[4]
R500  = gral[5]
R200  = gral[6]
R30   = np.ones(len(R200))*30.
R50   = np.ones(len(R200))*50.
R1    = 0.1*R500
R = np.vstack((R30,R50,R1,R1000,R500,R200)).T
Rp = np.array((R.tolist())*3)


plt.figure()
for i in range(len(R1)):
    j = np.argsort(lM)[i]
    plt.plot(np.log10(R[j]),S[j],'C1',alpha=0.25+i*0.01)
plt.plot(np.median(np.log10(R),axis=0),np.median(S,axis=0),'k')
plt.xlabel('$\log(R/kpc)$')
plt.ylabel('$S$')
plt.axis([np.log10(30),np.log10(R200.max()),0,1])
plt.savefig(path_plots+'S_R.png')

plt.figure()
for i in range(len(R1)):
    j = np.argsort(lM)[i]
    plt.plot(np.log10(R[j]),T[j],'C0',alpha=0.25+i*0.01)
plt.plot(np.median(np.log10(R),axis=0),np.median(T,axis=0),'k')
plt.axis([np.log10(30),np.log10(R200.max()),0,1])
plt.xlabel('$\log(R/kpc)$')
plt.ylabel('$T$')
plt.savefig(path_plots+'T_R.png')


plt.figure()
for i in range(len(R1)*3):
    j = np.argsort(lMp)[i]
    plt.plot(np.log10(Rp[j]),q[j],'C2',alpha=0.25+i*0.003)
plt.plot(np.log10(Rp[0]),np.median(q,axis=0),'k')
plt.xlabel('$\log(R/kpc)$')
plt.ylabel('$q$')
plt.axis([np.log10(30),np.log10(R200.max()),0,1])
plt.savefig(path_plots+'q_R.png')

plt.figure()
for i in range(len(R1)):
    j = np.argsort(lM)[i]
    plt.plot(np.log10(R[j]),t3D[j],'C3',alpha=0.25+i*0.01)
plt.plot(np.median(np.log10(R),axis=0),np.median(t3D,axis=0),'k')
plt.axis([np.log10(30),np.log10(R200.max()),0,90])
plt.xlabel('$\log(R/kpc)$')
plt.ylabel(r'$\theta$')
plt.savefig(path_plots+'t3D_R.png')


plt.figure()
for i in range(len(R1)*3):
    j = np.argsort(lMp)[i]
    plt.plot(np.log10(Rp[j]),t2D[j],'C4',alpha=0.25+i*0.003)
plt.plot(np.median(np.log10(Rp),axis=0),np.median(t2D,axis=0),'k')
plt.xlabel('$\log(R/kpc)$')
plt.ylabel(r'$\theta_{2D}$')
plt.axis([np.log10(30),np.log10(R200.max()),0,90])
plt.savefig(path_plots+'t2D_R.png')

######################

m = lM < 15.4
mp = lMp < 15.4

# S plots

plt.figure()
plt.plot(dm.S,s30.S,'b^',label='30 kpc p='+np.str(np.round(pearsonr(dm.S,s30.S)[0],2)))
plt.plot(dm.S,s50.S,'C3o',label='50 kpc p='+np.str(np.round(pearsonr(dm.S,s50.S)[0],2)))
plt.plot(dm.S,s100.S,'C1v',label='0.1R500 p='+np.str(np.round(pearsonr(dm.S,s100.S)[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*}$')
plt.savefig(path_plots+'SS_BCG.png')

plt.figure()
plt.plot(dm.S,s1000.S,'b^',label='R1000 p='+np.str(np.round(pearsonr(dm.S,s1000.S)[0],2)))
plt.plot(dm.S,s500.S,'C3o',label='R500 p='+np.str(np.round(pearsonr(dm.S,s500.S)[0],2)))
plt.plot(dm.S,s200.S,'C1v',label='R200 p='+np.str(np.round(pearsonr(dm.S,s200.S)[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*}$')
plt.savefig(path_plots+'SS_ext.png')


    
plt.figure()
plot_binned(lM[m],dm.S[m],'Dark matter','k')
plot_binned(lM[m],s30.S[m],'30 kpc','C1','')
plot_binned(lM[m],s50.S[m],'50 kpc','C1','--')
plot_binned(lM[m],s100.S[m],'0.1R500','C1',':')
plt.plot(lM[m],s_dm(lM[m]),'C7--',label='Tenneti et al 2014')
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')
plt.savefig(path_plots+'S_BCG.png')

plt.figure()
plot_binned(lM[m],dm.S[m],'Dark matter','k')
plot_binned(lM[m],s1000.S[m],'R1000','C1','')
plot_binned(lM[m],s500.S[m],'R500','C1','--')
plot_binned(lM[m],s200.S[m],'R200','C1',':')
plt.plot(lM[m],s_dm(lM[m]),'C7--',label='Tenneti et al 2014')
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')
plt.savefig(path_plots+'S_ext.png')

# T plots

plt.figure()
plt.plot(dm.T,s30.T,'b^',label='30 kpc p='+np.str(np.round(pearsonr(dm.T,s30.T)[0],2)))
plt.plot(dm.T,s50.T,'C3o',label='50 kpc p='+np.str(np.round(pearsonr(dm.T,s50.T)[0],2)))
plt.plot(dm.T,s100.T,'C1v',label='0.1R500 p='+np.str(np.round(pearsonr(dm.T,s100.T)[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*}$')
plt.savefig(path_plots+'TT_BCG.png')

plt.figure()
plt.plot(dm.T,s1000.T,'b^',label='R1000 p='+np.str(np.round(pearsonr(dm.T,s1000.T)[0],2)))
plt.plot(dm.T,s500.T,'C3o',label='R500 p='+np.str(np.round(pearsonr(dm.T,s500.T)[0],2)))
plt.plot(dm.T,s200.T,'C1v',label='R200 p='+np.str(np.round(pearsonr(dm.T,s200.T)[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*}$')
plt.savefig(path_plots+'TT_ext.png')



plt.figure()
plot_binned(lM[m],dm.T[m],'Dark matter','k')
plot_binned(lM[m],s30.T[m],'30 kpc','C0','')
plot_binned(lM[m],s50.T[m],'50 kpc','C0','--')
plot_binned(lM[m],s100.T[m],'0.1R500','C0',':')
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$T$')
plt.savefig(path_plots+'T_BCG.png')

plt.figure()
plot_binned(lM[m],dm.T[m],'Dark matter','k')
plot_binned(lM[m],s1000.T[m],'R1000','C0','')
plot_binned(lM[m],s500.T[m],'R500','C0','--')
plot_binned(lM[m],s200.T[m],'R200','C0',':')
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$T$')
plt.savefig(path_plots+'T_ext.png')

# q plots

plt.figure()
plt.plot(dm.q,s30.q,'b^',label='30 kpc p='+np.str(np.round(pearsonr(dm.q,s30.q)[0],2)))
plt.plot(dm.q,s50.q,'C3o',label='50 kpc p='+np.str(np.round(pearsonr(dm.q,s50.q)[0],2)))
plt.plot(dm.q,s100.q,'C1v',label='0.1R500 p='+np.str(np.round(pearsonr(dm.q,s100.q)[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*}$')
plt.savefig(path_plots+'qq_BCG.png')

plt.figure()
plt.plot(dm.q,s1000.q,'b^',label='R1000 p='+np.str(np.round(pearsonr(dm.q,s1000.q)[0],2)))
plt.plot(dm.q,s500.q,'C3o',label='R500 p='+np.str(np.round(pearsonr(dm.q,s500.q)[0],2)))
plt.plot(dm.q,s200.q,'C1v',label='R200 p='+np.str(np.round(pearsonr(dm.q,s200.q)[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*}$')
plt.savefig(path_plots+'qq_ext.png')



plt.figure()
plot_binned(lMp[mp],dm.q[mp],'Dark matter','k')
plot_binned(lMp[mp],s30.q[mp],'30 kpc','C2','')
plot_binned(lMp[mp],s50.q[mp],'50 kpc','C2','--')
plot_binned(lMp[mp],s100.q[mp],'0.1R500','C2',':')
plt.plot(lM[m],q_dm(lM[m]),'C7--',label='Tenneti et al 2014')
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q$')
plt.savefig(path_plots+'q_BCG.png')

plt.figure()
plot_binned(lMp[mp],dm.q[mp],'Dark matter','k')
plot_binned(lMp[mp],s1000.q[mp],'R1000','C2','')
plot_binned(lMp[mp],s500.q[mp],'R500','C2','--')
plot_binned(lMp[mp],s200.q[mp],'R200','C2',':')
plt.plot(lM[m],q_dm(lM[m]),'C7--',label='Tenneti et al 2014')
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q$')
plt.savefig(path_plots+'q_ext.png')

# THETA PLOTS

plt.figure()
plot_binned(lM[m],t30_3D[m],'30 kpc','C3','')
plot_binned(lM[m],t50_3D[m],'50 kpc','C3','--')
plot_binned(lM[m],t100_3D[m],'0.1R500','C3',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_BCG.png')

plot_binned(lM[m*(s30.S < 0.6)] , t30_3D[m*(s30.S < 0.6)],'30 kpc','C1','')
plot_binned(lM[m*(s50.S < 0.6)] , t50_3D[m*(s50.S < 0.6)],'50 kpc','C1','--')
plot_binned(lM[m*(s100.S < 0.5)],t100_3D[m*(s100.S < 0.5)],'0.1R500','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_BCG_2.png')


plt.figure()
plot_binned(lM[m],t1000_3D[m],'R1000','C3','')
plot_binned(lM[m],t500_3D[m],'R500','C3','--')
plot_binned(lM[m],t200_3D[m],'R200','C3',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_ext.png')



plt.figure()
plot_binned(lMp[mp],t30_2D[mp],'30 kpc','C4','')
plot_binned(lMp[mp],t50_2D[mp],'50 kpc','C4','--')
plot_binned(lMp[mp],t100_2D[mp],'0.1R500','C4',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{2D}$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_BCG.png')

plot_binned(lMp[mp*(s30.q < 0.7)] , t30_2D[mp*(s30.q < 0.7)],'30 kpc','C1','')
plot_binned(lMp[mp*(s50.q < 0.7)] , t50_2D[mp*(s50.q < 0.7)],'50 kpc','C1','--')
plot_binned(lMp[mp*(s100.q < 0.7)],t100_2D[mp*(s100.q < 0.7)],'0.1R500','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{2D}$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_BCG_2.png')


plt.figure()
plot_binned(lMp[mp],t1000_2D[mp],'R1000','C4','')
plot_binned(lMp[mp],t500_2D[mp],'R500','C4','--')
plot_binned(lMp[mp],t200_2D[mp],'R200','C4',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{2D}$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_ext.png')
