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
path_plots = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/plots/gal_plots/'
gral  = np.loadtxt(path+'gral_091.dat').T
lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)

dm = DarkMatter()

g200  = Galaxias(200)
g500  = Galaxias(500)
g1000 = Galaxias(1000)

g200x  = Galaxias(200,True)
g500x  = Galaxias(500,True)
g1000x = Galaxias(1000,True)

mN200  = (g200.N > 9)
mN500  = (g500.N > 9)
mN1000 = (g1000.N > 9)

mN200x  = (g200x.N > 9)
mN500x  = (g500x.N > 9)
mN1000x = (g1000x.N > 9)

mn200  = (np.array(g200.N.tolist()*3) > 9)
mn500  = (np.array(g500.N.tolist()*3) > 9)
mn1000 = (np.array(g1000.N.tolist()*3) > 9)

mn200x  = (np.array(g200x.N.tolist()*3) > 9)
mn500x  = (np.array(g500x.N.tolist()*3) > 9)
mn1000x = (np.array(g1000x.N.tolist()*3) > 9)


ct200_3D , t200_3D    = cosangle(dm.a3D,g200.a3D) 
ct500_3D , t500_3D    = cosangle(dm.a3D,g500.a3D)
ct1000_3D, t1000_3D   = cosangle(dm.a3D,g1000.a3D)
ct200x_3D , t200x_3D    = cosangle(dm.a3D,g200x.a3D) 
ct500x_3D , t500x_3D    = cosangle(dm.a3D,g500x.a3D)
ct1000x_3D, t1000x_3D   = cosangle(dm.a3D,g1000x.a3D)

ct200_2D , t200_2D    = cosangle(dm.a2D,g200.a2D) 
ct500_2D , t500_2D    = cosangle(dm.a2D,g500.a2D)
ct1000_2D, t1000_2D   = cosangle(dm.a2D,g1000.a2D)
ct200x_2D , t200x_2D    = cosangle(dm.a2D,g200x.a2D) 
ct500x_2D , t500x_2D    = cosangle(dm.a2D,g500x.a2D)
ct1000x_2D, t1000x_2D   = cosangle(dm.a2D,g1000x.a2D)


m = (lM < 15.4)
mp = lMp < 15.4

# S plots

plt.figure()
plt.plot(dm.S,g1000.S,'b^',label='R1000',alpha=0.5)
plt.plot(dm.S,g500.S,'C3o',label='R500')
plt.plot(dm.S,g200.S,'C1v',label='R200')
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{GAL}$')
plt.savefig(path_plots+'SS_allgal.png')


plt.figure()
plt.plot(dm.S,g1000x.S,'b^',label='R1000',alpha=0.5)
plt.plot(dm.S,g500x.S,'C3o',label='R500')
plt.plot(dm.S,g200x.S,'C1v',label='R200')
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{GAL}$')
plt.savefig(path_plots+'SS_galcut.png')


plt.figure()
plot_binned(lM[m],dm.S[m],'Dark matter','k')
plot_binned(lM[m*mN1000],g1000.S[m*mN1000],'R1000','C1','-')
plot_binned(lM[m*mN500] ,g500.S[m*mN500],'R500','C1','--')
plot_binned(lM[m*mN200],g200.S[m*mN200],'R200','C1',':')
plt.legend(loc = 3)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')
plt.savefig(path_plots+'S_allgal.png')

plt.figure()
plot_binned(lM[m],dm.S[m],'Dark matter','k')
plot_binned(lM[m*mN1000x],g1000x.S[m*mN1000x],'R1000','C1','')
plot_binned(lM[m*mN500x] ,g500x.S[m*mN500x],'R500','C1','--')
plot_binned(lM[m*mN200x],g200x.S[m*mN200x],'R200','C1',':')
plt.legend(loc = 3)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S_x$')
plt.savefig(path_plots+'S_galcut.png')


plt.figure()
plot_binned(g1000.N[m*mN1000],g1000.S[m*mN1000],'R1000','C1','-')
plot_binned(g500.N[m*mN500]  ,g500.S[m*mN500]  ,'R500','C1','--')
plot_binned(g200.N[m*mN200]  ,g200.S[m*mN200]  ,'R200','C1',':')
plt.legend(loc = 3)
plt.axis([10,650,0,1])
plt.xlabel('$N_{GAL}$')
plt.ylabel('$S$')
plt.savefig(path_plots+'S_allgal_N.png')

plt.figure()
plot_binned(g1000x.N[m*mN1000x],g1000x.S[m*mN1000x],'R1000','C1','-')
plot_binned(g500x.N[m*mN500x]  ,g500x.S[m*mN500x]  ,'R500','C1','--')
plot_binned(g200x.N[m*mN200x]  ,g200x.S[m*mN200x]  ,'R200','C1',':')
plt.legend(loc = 3)
plt.axis([10,650,0,1])
plt.xlabel('$N_{GAL}$')
plt.ylabel('$S_x$')
plt.savefig(path_plots+'S_galcut_N.png')


# T plots

plt.figure()
plt.plot(dm.T,g1000.T,'b^',label='R1000',alpha=0.5)
plt.plot(dm.T,g500.T,'C3o',label='R500')
plt.plot(dm.T,g200.T,'C1v',label='R200')
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{GAL}$')
plt.savefig(path_plots+'TT_allgal.png')


plt.figure()
plt.plot(dm.T,g1000x.T,'b^',label='R1000',alpha=0.5)
plt.plot(dm.T,g500x.T,'C3o',label='R500')
plt.plot(dm.T,g200x.T,'C1v',label='R200')
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{GAL}$')
plt.savefig(path_plots+'TT_galcut.png')



plt.figure()
plot_binned(lM[m],dm.T[m],'Dark matter','k')
plot_binned(lM[m*mN1000],g1000.T[m*mN1000],'R1000','C0','')
plot_binned(lM[m*mN500] ,g500.T[m*mN500],'R500','C0','--')
plot_binned(lM[m*mN200],g200.T[m*mN200],'R200','C0',':')
plt.legend(loc = 3)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$T$')
plt.savefig(path_plots+'T_allgal.png')

plt.figure()
plot_binned(lM[m],dm.T[m],'Dark matter','k')
plot_binned(lM[m*mN1000x],g1000x.T[m*mN1000x],'R1000','C0','')
plot_binned(lM[m*mN500x] ,g500x.T[m*mN500x],'R500','C0','--')
plot_binned(lM[m*mN200x],g200x.T[m*mN200x],'R200','C0',':')
plt.legend(loc = 3)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$Tx$')
plt.savefig(path_plots+'T_galcut.png')


plt.figure()
plot_binned(g1000.N[m*mN1000],g1000.T[m*mN1000],'R1000','C0','-')
plot_binned(g500.N[m*mN500]  ,g500.T[m*mN500]  ,'R500','C0','--')
plot_binned(g200.N[m*mN200]  ,g200.T[m*mN200]  ,'R200','C0',':')
plt.legend(loc = 3)
plt.axis([10,650,0,1])
plt.xlabel('$N_{GAL}$')
plt.ylabel('$T$')
plt.savefig(path_plots+'T_allgal_N.png')

plt.figure()
plot_binned(g1000x.N[m*mN1000x],g1000x.T[m*mN1000x],'R1000','C0','-')
plot_binned(g500x.N[m*mN500x]  ,g500x.T[m*mN500x]  ,'R500','C0','--')
plot_binned(g200x.N[m*mN200x]  ,g200x.T[m*mN200x]  ,'R200','C0',':')
plt.legend(loc = 3)
plt.axis([10,650,0,1])
plt.xlabel('$N_{GAL}$')
plt.ylabel('$T_x$')
plt.savefig(path_plots+'T_galcut_N.png')


# q plots

plt.figure()
plt.plot(dm.q,g1000.q,'b^',label='R1000',alpha=0.5)
plt.plot(dm.q,g500.q,'C3o',label='R500')
plt.plot(dm.q,g200.q,'C1v',label='R200')
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{GAL}$')
plt.savefig(path_plots+'qq_allgal.png')


plt.figure()
plt.plot(dm.q,g1000x.q,'b^',label='R1000',alpha=0.5)
plt.plot(dm.q,g500x.q,'C3o',label='R500')
plt.plot(dm.q,g200x.q,'C1v',label='R200')
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{GAL}$')
plt.savefig(path_plots+'qq_galcut.png')



plt.figure()
plot_binned(lMp[mp],dm.q[mp],'Dark matter','k')
plot_binned(lMp[mp*mn1000],g1000.q[mp*mn1000],'R1000','C2','')
plot_binned(lMp[mp*mn500] ,g500.q[mp*mn500],'R500','C2','--')
plot_binned(lMp[mp*mn200],g200.q[mp*mn200],'R200','C2',':')
plt.legend(loc = 3)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q$')
plt.savefig(path_plots+'q_allgal.png')

plt.figure()
plot_binned(lMp[mp],dm.q[mp],'Dark matter','k')
plot_binned(lMp[mp*mn1000x],g1000x.q[mp*mn1000x],'R1000','C2','')
plot_binned(lMp[mp*mn500x] ,g500x.q[mp*mn500x],'R500','C2','--')
plot_binned(lMp[mp*mn200x],g200x.q[mp*mn200x],'R200','C2',':')
plt.legend(loc = 3)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q_x$')
plt.savefig(path_plots+'q_galcut.png')

plt.figure()
plot_binned(g1000.n[mp*mn1000],g1000.q[mp*mn1000],'R1000','C2','-')
plot_binned(g500.n[mp*mn500]  ,g500.q[mp*mn500]  ,'R500','C2','--')
plot_binned(g200.n[mp*mn200]  ,g200.q[mp*mn200]  ,'R200','C2',':')
plt.legend(loc = 3)
plt.axis([10,650,0,1])
plt.xlabel('$N_{GAL}$')
plt.ylabel('$q$')
plt.savefig(path_plots+'q_allgal_N.png')

plt.figure()
plot_binned(g1000x.n[mp*mn1000x],g1000x.q[mp*mn1000x],'R1000','C2','-')
plot_binned(g500x.n[mp*mn500x]  ,g500x.q[mp*mn500x]  ,'R500','C2','--')
plot_binned(g200x.n[mp*mn200x]  ,g200x.q[mp*mn200x]  ,'R200','C2',':')
plt.legend(loc = 3)
plt.axis([10,650,0,1])
plt.xlabel('$N_{GAL}$')
plt.ylabel('$q_x$')
plt.savefig(path_plots+'q_galcut_N.png')



# theta plots

plt.figure()
plot_binned(lM[m*mN1000],t1000_3D[m*mN1000],'R1000','C3','')
plot_binned(lM[m*mN500] ,t500_3D[m*mN500],'R500','C3','--')
plot_binned(lM[m*mN200],t200_3D[m*mN200],'R200','C3',':')
plt.legend(loc = 1)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta$')
plt.savefig(path_plots+'theta3D_allgal.png')

plt.figure()
plot_binned(lM[m*mN1000x],t1000x_3D[m*mN1000x],'R1000','C3','')
plot_binned(lM[m*mN500x] ,t500x_3D[m*mN500x],'R500','C3','--')
plot_binned(lM[m*mN200x],t200x_3D[m*mN200x],'R200','C3',':')
plt.legend(loc = 1)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_x$')
plt.savefig(path_plots+'theta3D_galcut.png')


plt.figure()
plot_binned(lMp[mp*mn1000],t1000_2D[mp*mn1000],'R1000','C4','')
plot_binned(lMp[mp*mn500] ,t500_2D[mp*mn500],'R500','C4','--')
plot_binned(lMp[mp*mn200],t200_2D[mp*mn200],'R200','C4',':')
plt.legend(loc = 1)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta^{2D}$')
plt.savefig(path_plots+'theta2D_allgal.png')

plt.figure()
plot_binned(lMp[mp*mn1000x],t1000x_2D[mp*mn1000x],'R1000','C4','')
plot_binned(lMp[mp*mn500x] ,t500x_2D[mp*mn500x],'R500','C4','--')
plot_binned(lMp[mp*mn200x],t200x_2D[mp*mn200x],'R200','C4',':')
plt.legend(loc = 1)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta^{2D}_x$')
plt.savefig(path_plots+'theta2D_galcut.png')


# with n

plt.figure()
plot_binned(g1000.N[m*mN1000],t1000_3D[m*mN1000],'R1000','C3','')
plot_binned(g500.N[m*mN500] ,t500_3D[m*mN500],'R500','C3','--')
plot_binned(g200.N[m*mN200],t200_3D[m*mN200],'R200','C3',':')
plt.legend(loc = 1)
plt.ylim([0,50])
plt.xlabel('$N_{GAL}$')
plt.ylabel(r'$\theta$')
plt.savefig(path_plots+'theta3D_allgal_N.png')

plt.figure()
plot_binned(g1000.N[m*mN1000x],t1000x_3D[m*mN1000x],'R1000','C3','')
plot_binned(g500.N[m*mN500x] ,t500x_3D[m*mN500x],'R500','C3','--')
plot_binned(g200.N[m*mN200x],t200x_3D[m*mN200x],'R200','C3',':')
plt.legend(loc = 1)
plt.ylim([0,50])
plt.xlabel('$N_{GAL}$')
plt.ylabel(r'$\theta_x$')
plt.savefig(path_plots+'theta3D_galcut_N.png')


plt.figure()
plot_binned(g1000.n[mp*mn1000],t1000_2D[mp*mn1000],'R1000','C4','')
plot_binned(g500.n[mp*mn500] ,t500_2D[mp*mn500],'R500','C4','--')
plot_binned(g200.n[mp*mn200],t200_2D[mp*mn200],'R200','C4',':')
plt.legend(loc = 1)
plt.ylim([0,50])
plt.xlabel('$N_{GAL}$')
plt.ylabel(r'$\theta^{2D}$')
plt.savefig(path_plots+'theta2D_allgal_N.png')

plt.figure()
plot_binned(g1000.n[mp*mn1000x],t1000x_2D[mp*mn1000x],'R1000','C4','')
plot_binned(g500.n[mp*mn500x] ,t500x_2D[mp*mn500x],'R500','C4','--')
plot_binned(g200.n[mp*mn200x],t200x_2D[mp*mn200x],'R200','C4',':')
plt.legend(loc = 1)
plt.ylim([0,50])
plt.xlabel('$N_{GAL}$')
plt.ylabel(r'$\theta^{2D}_x$')
plt.savefig(path_plots+'theta2D_galcut_N.png')
