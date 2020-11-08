import sys
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from scipy import stats
from pylab import *
from main import *
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

path = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/catalog/'
gral  = np.loadtxt(path+'gral_091.dat').T
lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)

dm = DarkMatter()

s30   = Stars(30)
s50   = Stars(50)
s100  = Stars(100)
s200  = Stars(200)
s500  = Stars(500)
s1000 = Stars(1000)

ct30_3D  , t30_3D     = cosangle(dm.a3D,s30.a3D) 
ct50_3D  , t50_3D     = cosangle(dm.a3D,s50.a3D)  
ct100_3D , t100_3D    = cosangle(dm.a3D,s100.a3D) 
ct200_3D , t200_3D    = cosangle(dm.a3D,s200.a3D) 
ct500_3D , t500_3D    = cosangle(dm.a3D,s500.a3D)
ct1000_3D, t1000_3D   = cosangle(dm.a3D,s1000.a3D)


m = lM < 15.4
mp = lMp < 15.4

def plot_binned(X,Y,label,color='C3',shade=1,nbins=10):
    x,y,s = binned(X,Y,nbins)
    plt.plot(x,y,color,label=label,alpha = shade)
    # plt.plot(x,y+s,color+'--',alpha = shade)
    # plt.plot(x,y-s,color+'--', alpha = shade)
    plt.errorbar(x,y,yerr=s,fmt='none',ecolor=color,alpha = shade)
    
plt.figure()
plot_binned(lM[m],dm.S[m],'Dark matter','k')
plot_binned(lM[m],s30.S[m],'30 kpc','C1',1)
plot_binned(lM[m],s50.S[m],'50 kpc','C1',0.7)
plot_binned(lM[m],s100.S[m],'0.1R500','C1',0.5)
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')

plt.figure()
plot_binned(lM[m],dm.S[m],'Dark matter','k')
plot_binned(lM[m],s1000.S[m],'R1000','C1',1)
plot_binned(lM[m],s500.S[m],'R500','C1',0.7)
plot_binned(lM[m],s200.S[m],'R200','C1',0.5)
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')

plt.figure()
plot_binned(lM[m],dm.T[m],'Dark matter','k')
plot_binned(lM[m],T30.T[m],'30 kpc','C0',1)
plot_binned(lM[m],T50.T[m],'50 kpc','C0',0.7)
plot_binned(lM[m],T100.T[m],'0.1R500','C0',0.5)
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$T$')

plt.figure()
plot_binned(lM[m],dm.T[m],'Dark matter','k')
plot_binned(lM[m],T1000.T[m],'R1000','C0',1)
plot_binned(lM[m],T500.T[m],'R500','C0',0.7)
plot_binned(lM[m],T200.T[m],'R200','C0',0.5)
plt.legend()
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$T$')

plq.figure()
ploq_binned(lMp[mp],dm.q[mp],'Dark maqqer','k')
ploq_binned(lMp[mp],q30.q[mp],'30 kpc','C2',1)
ploq_binned(lMp[mp],q50.q[mp],'50 kpc','C2',0.7)
ploq_binned(lMp[mp],q100.q[mp],'0.1R500','C2',0.5)
plq.legend()
plq.ylim([0,1])
plq.xlabel('$\log M_{200}$')
plq.ylabel('$q$')

plq.figure()
ploq_binned(lMp[mp],dm.q[mp],'Dark maqqer','k')
ploq_binned(lMp[mp],q1000.q[mp],'R1000','C2',1)
ploq_binned(lMp[mp],q500.q[mp],'R500','C2',0.7)
ploq_binned(lMp[mp],q200.q[mp],'R200','C2',0.5)
plq.legend()
plq.ylim([0,1])
plq.xlabel('$\log M_{200}$')
plq.ylabel('$q$')


plt.figure()
plot_binned(lM[m],t30_3D[m],'30 kpc','C1',1)
plot_binned(lM[m],t50_3D[m],'50 kpc','C1',0.7)
plot_binned(lM[m],t100_3D[m],'0.1R500','C1',0.5)
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta$')
plt.ylim([0,50])

plt.figure()
plot_binned(lM[m],t1000_3D[m],'R1000','C2',1)
plot_binned(lM[m],t500_3D[m],'R500','C2',0.7)
plot_binned(lM[m],t200_3D[m],'R200','C2',0.5)
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta$')
plt.ylim([0,50])

