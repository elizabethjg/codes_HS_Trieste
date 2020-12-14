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

path = '../catalog/'


indicator = 'off'

path_plots = '../plots/T_plots/stars/'
os.system('mkdir '+path_plots)

gral  = np.loadtxt(path+'gral_091_2.dat').T
lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)

dm = DarkMatter()

gral  = np.loadtxt(path+'gral_091_2.dat').T
lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)

off   = gral[13]
off2D = np.concatenate((gral[14],gral[15],gral[16]))
DV    = gral[17]
DV2D  = np.concatenate((gral[18],gral[19],gral[20]))
gap  = gral[21]
gap2D = np.concatenate((gral[22],gral[23],gral[24]))

if 'off' in indicator:
    gap = off
    gap2D = off2D
elif 'DV' in indicator:
    gap = DV
    gap2D = DV2D

T2D  = np.array((dm.T.tolist())*3)

mpro = dm.T > 0.5
mobl = ~mpro

mpro2D = T2D > 0.5
mobl2D = ~mpro2D



m = lM < 18
mp = lMp < 18

mgap = binned(lM,gap)[-1]
mgap2D = binned(lMp,gap2D)[-1]

mold = ~mgap
mold2D = ~mgap2D

mnew = mgap
mnew2D = mgap2D


# mold = gap < 0.2
# mnew = gap > 0.3

# mold2D = gap2D < 0.2
# mnew2D = gap2D > 0.3



# '''
# WITH STARS

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

ct30_2D  , t30_2D     = cosangle(dm.a2D,s30.a2D) 
ct50_2D  , t50_2D     = cosangle(dm.a2D,s50.a2D)  
ct100_2D , t100_2D    = cosangle(dm.a2D,s100.a2D) 
ct200_2D , t200_2D    = cosangle(dm.a2D,s200.a2D) 
ct500_2D , t500_2D    = cosangle(dm.a2D,s500.a2D)
ct1000_2D, t1000_2D   = cosangle(dm.a2D,s1000.a2D)

t3D = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
t2D = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T

S = np.vstack((s30.S,s50.S,s100.S,s1000.S,s500.S,s200.S)).T
T = np.vstack((s30.T,s50.T,s100.T,s1000.T,s500.T,s200.T)).T
q = np.vstack((s30.q,s50.q,s100.q,s1000.q,s500.q,s200.q)).T


# ANGLES AND RADIUS

R1000 = gral[4]
R500  = gral[5]
R200  = gral[6]
R30   = np.ones(len(R200))*30.
R50   = np.ones(len(R200))*50.
R1    = 0.1*R500
R = np.vstack((R30,R50,R1,R1000,R500,R200)).T
Rp = np.array((R.tolist())*3)


plt.figure()
for j in range(mpro.sum()):
        plt.plot(np.log10(R[mpro][j]),S[mpro][j],'C1',alpha=0.5)
for j in range(mobl.sum()):
        plt.plot(np.log10(R[mobl][j]),S[mobl][j],'C2',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'S_R_ind.png')

plt.figure()
for j in range(mpro.sum()):
        plt.plot(np.log10(R[mpro][j]),T[mpro][j],'C1',alpha=0.5)
for j in range(mobl.sum()):
        plt.plot(np.log10(R[mobl][j]),T[mobl][j],'C2',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'T_R_ind.png')

plt.figure()
for j in range(mpro2D.sum()):
        plt.plot(np.log10(Rp[mpro2D][j]),q[mpro2D][j],'C1',alpha=0.5)
for j in range(mobl2D.sum()):
        plt.plot(np.log10(Rp[mobl2D][j]),q[mobl2D][j],'C2',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.,0.9])
plt.savefig(path_plots+'q_R_ind.png')


plt.figure()
for j in range(mpro.sum()):
        plt.plot(np.log10(R[mpro][j]),t3D[mpro][j],'C1',alpha=0.5)
for j in range(mobl.sum()):
        plt.plot(np.log10(R[mobl][j]),t3D[mobl][j],'C2',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.,90])
plt.savefig(path_plots+'t3D_R_ind.png')


plt.figure()
for j in range(mpro2D.sum()):
        plt.plot(np.log10(Rp[mpro2D][j]),t2D[mpro2D][j],'C1',alpha=0.5)
for j in range(mobl2D.sum()):
        plt.plot(np.log10(Rp[mobl2D][j]),t2D[mobl2D][j],'C2',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.,90])
plt.savefig(path_plots+'t2D_R_ind.png')


plt.figure()    
plt.plot(np.median(np.log10(R),axis=0),np.median(S,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(S,axis=0)+np.std(S,axis=0),np.median(S,axis=0)-np.std(S,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mpro]),axis=0),np.median(S[mpro],axis=0),'C1')
plt.fill_between(np.median(np.log10(R[mpro]),axis=0),np.median(S[mpro],axis=0)+np.std(S[mpro],axis=0),np.median(S[mpro],axis=0)-np.std(S[mpro],axis=0),color = 'C1',alpha=0.1)
plt.plot(np.median(np.log10(R[mobl]),axis=0),np.median(S[mobl],axis=0),'C2')
plt.fill_between(np.median(np.log10(R[mobl]),axis=0),np.median(S[mobl],axis=0)+np.std(S[mobl],axis=0),np.median(S[mobl],axis=0)-np.std(S[mobl],axis=0),color = 'C2',alpha=0.1)
plt.ylabel('$S$')
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'S_R.png')

plt.figure()
plt.plot(np.median(np.log10(R),axis=0),np.median(T,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(T,axis=0)+np.std(T,axis=0),np.median(T,axis=0)-np.std(T,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mpro]),axis=0),np.median(T[mpro],axis=0),'C1')
plt.fill_between(np.median(np.log10(R[mpro]),axis=0),np.median(T[mpro],axis=0)+np.std(T[mpro],axis=0),np.median(T[mpro],axis=0)-np.std(T[mpro],axis=0),color = 'C1',alpha=0.1)
plt.plot(np.median(np.log10(R[mobl]),axis=0),np.median(T[mobl],axis=0),'C2')
plt.fill_between(np.median(np.log10(R[mobl]),axis=0),np.median(T[mobl],axis=0)+np.std(T[mobl],axis=0),np.median(T[mobl],axis=0)-np.std(T[mobl],axis=0),color = 'C2',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel('$T$')
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'T_R.png')

plt.figure()
plt.plot(np.median(np.log10(Rp),axis=0),np.median(q,axis=0),'k')
plt.fill_between(np.median(np.log10(Rp),axis=0),np.median(q,axis=0)+np.std(q,axis=0),np.median(q,axis=0)-np.std(q,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mpro2D]),axis=0),np.median(q[mpro2D],axis=0),'C1')
plt.fill_between(np.median(np.log10(Rp[mpro2D]),axis=0),np.median(q[mpro2D],axis=0)+np.std(q[mpro2D],axis=0),np.median(q[mpro2D],axis=0)-np.std(q[mpro2D],axis=0),color = 'C1',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mobl2D]),axis=0),np.median(q[mobl2D],axis=0),'C2')
plt.fill_between(np.median(np.log10(Rp[mobl2D]),axis=0),np.median(q[mobl2D],axis=0)+np.std(q[mobl2D],axis=0),np.median(q[mobl2D],axis=0)-np.std(q[mobl2D],axis=0),color = 'C2',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel('$q$')
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'q_Rp.png')



plt.figure()
plt.plot(np.median(np.log10(R),axis=0),np.median(t3D,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(t3D,axis=0)+np.std(t3D,axis=0),np.median(t3D,axis=0)-np.std(t3D,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mpro]),axis=0),np.median(t3D[mpro],axis=0),'C1')
plt.fill_between(np.median(np.log10(R[mpro]),axis=0),np.median(t3D[mpro],axis=0)+np.std(t3D[mpro],axis=0),np.median(t3D[mpro],axis=0)-np.std(t3D[mpro],axis=0),color = 'C1',alpha=0.1)
plt.plot(np.median(np.log10(R[mobl]),axis=0),np.median(t3D[mobl],axis=0),'C2')
plt.fill_between(np.median(np.log10(R[mobl]),axis=0),np.median(t3D[mobl],axis=0)+np.std(t3D[mobl],axis=0),np.median(t3D[mobl],axis=0)-np.std(t3D[mobl],axis=0),color = 'C2',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel(r'$\theta$')
plt.axis([np.log10(30),np.log10(2000),0.,50.])
plt.savefig(path_plots+'t3D_R.png')

plt.figure()
plt.plot(np.median(np.log10(Rp),axis=0),np.median(t2D,axis=0),'k')
plt.fill_between(np.median(np.log10(Rp),axis=0),np.median(t2D,axis=0)+np.std(t2D,axis=0),np.median(t2D,axis=0)-np.std(t2D,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mpro2D]),axis=0),np.median(t2D[mpro2D],axis=0),'C1')
plt.fill_between(np.median(np.log10(Rp[mpro2D]),axis=0),np.median(t2D[mpro2D],axis=0)+np.std(t2D[mpro2D],axis=0),np.median(t2D[mpro2D],axis=0)-np.std(t2D[mpro2D],axis=0),color = 'C1',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mobl2D]),axis=0),np.median(t2D[mobl2D],axis=0),'C2')
plt.fill_between(np.median(np.log10(Rp[mobl2D]),axis=0),np.median(t2D[mobl2D],axis=0)+np.std(t2D[mobl2D],axis=0),np.median(t2D[mobl2D],axis=0)-np.std(t2D[mobl2D],axis=0),color = 'C2',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel(r'$\theta_{2D}$')
plt.axis([np.log10(30),np.log10(2000),0.,50.])
plt.savefig(path_plots+'t2D_Rp.png')


############ S plots
plt.scatter(dm.S[mobl],s30.S[mobl],c=lM[mobl])

plt.figure()
plt.plot(dm.S[mobl],s30.S[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.S[mobl],s30.S[mobl])[0],2)))
plt.plot(dm.S[mpro],s30.S[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.S[mpro],s30.S[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < 30kpc)$')
plt.savefig(path_plots+'S30_.png')

plt.figure()
plt.plot(dm.S[mobl],s50.S[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.S[mobl],s50.S[mobl])[0],2)))
plt.plot(dm.S[mpro],s50.S[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.S[mpro],s50.S[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < 50kpc)$')
plt.savefig(path_plots+'S50_.png')

plt.figure()
plt.plot(dm.S[mobl],s100.S[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.S[mobl],s100.S[mobl])[0],2)))
plt.plot(dm.S[mpro],s100.S[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.S[mpro],s100.S[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < 0.1R500)$')
plt.savefig(path_plots+'S100_.png')

plt.figure()
plt.plot(dm.S[mobl],s1000.S[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.S[mobl],s1000.S[mobl])[0],2)))
plt.plot(dm.S[mpro],s1000.S[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.S[mpro],s1000.S[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < R1000)$')
plt.savefig(path_plots+'S1000_.png')

plt.figure()
plt.plot(dm.S[mobl],s500.S[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.S[mobl],s500.S[mobl])[0],2)))
plt.plot(dm.S[mpro],s500.S[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.S[mpro],s500.S[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < R500)$')
plt.savefig(path_plots+'S500_.png')

plt.figure()
plt.plot(dm.S[mobl],s200.S[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.S[mobl],s200.S[mobl])[0],2)))
plt.plot(dm.S[mpro],s200.S[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.S[mpro],s200.S[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < R200)$')
plt.savefig(path_plots+'S200_.png')


############ T plots

plt.figure()
plt.plot(dm.T[mobl],s30.T[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.T[mobl],s30.T[mobl])[0],2)))
plt.plot(dm.T[mpro],s30.T[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.T[mpro],s30.T[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < 30kpc)$')
plt.savefig(path_plots+'T30_.png')

plt.figure()
plt.plot(dm.T[mobl],s50.T[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.T[mobl],s50.T[mobl])[0],2)))
plt.plot(dm.T[mpro],s50.T[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.T[mpro],s50.T[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < 50kpc)$')
plt.savefig(path_plots+'T50_.png')

plt.figure()
plt.plot(dm.T[mobl],s100.T[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.T[mobl],s100.T[mobl])[0],2)))
plt.plot(dm.T[mpro],s100.T[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.T[mpro],s100.T[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < 0.1R500)$')
plt.savefig(path_plots+'T100_.png')

plt.figure()
plt.plot(dm.T[mobl],s1000.T[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.T[mobl],s1000.T[mobl])[0],2)))
plt.plot(dm.T[mpro],s1000.T[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.T[mpro],s1000.T[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < R1000)$')
plt.savefig(path_plots+'T1000_.png')

plt.figure()
plt.plot(dm.T[mobl],s500.T[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.T[mobl],s500.T[mobl])[0],2)))
plt.plot(dm.T[mpro],s500.T[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.T[mpro],s500.T[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < R500)$')
plt.savefig(path_plots+'T500_.png')

plt.figure()
plt.plot(dm.T[mobl],s200.T[mobl],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.T[mobl],s200.T[mobl])[0],2)))
plt.plot(dm.T[mpro],s200.T[mpro],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.T[mpro],s200.T[mpro])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < R200)$')
plt.savefig(path_plots+'T200_.png')

############ q plots

plt.figure()
plt.plot(dm.q[mobl2D],s30.q[mobl2D],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.q[mobl2D],s30.q[mobl2D])[0],2)))
plt.plot(dm.q[mpro2D],s30.q[mpro2D],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.q[mpro2D],s30.q[mpro2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < 30kpc)$')
plt.savefig(path_plots+'q30_.png')

plt.figure()
plt.plot(dm.q[mobl2D],s50.q[mobl2D],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.q[mobl2D],s50.q[mobl2D])[0],2)))
plt.plot(dm.q[mpro2D],s50.q[mpro2D],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.q[mpro2D],s50.q[mpro2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < 50kpc)$')
plt.savefig(path_plots+'q50_.png')

plt.figure()
plt.plot(dm.q[mobl2D],s100.q[mobl2D],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.q[mobl2D],s100.q[mobl2D])[0],2)))
plt.plot(dm.q[mpro2D],s100.q[mpro2D],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.q[mpro2D],s100.q[mpro2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < 0.1R500)$')
plt.savefig(path_plots+'q100_.png')

plt.figure()
plt.plot(dm.q[mobl2D],s1000.q[mobl2D],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.q[mobl2D],s1000.q[mobl2D])[0],2)))
plt.plot(dm.q[mpro2D],s1000.q[mpro2D],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.q[mpro2D],s1000.q[mpro2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < R1000)$')
plt.savefig(path_plots+'q1000_.png')

plt.figure()
plt.plot(dm.q[mobl2D],s500.q[mobl2D],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.q[mobl2D],s500.q[mobl2D])[0],2)))
plt.plot(dm.q[mpro2D],s500.q[mpro2D],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.q[mpro2D],s500.q[mpro2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < R500)$')
plt.savefig(path_plots+'q500_.png')

plt.figure()
plt.plot(dm.q[mobl2D],s200.q[mobl2D],'C2^',label='oblate p='+np.str(np.round(pearsonr(dm.q[mobl2D],s200.q[mobl2D])[0],2)))
plt.plot(dm.q[mpro2D],s200.q[mpro2D],'C1^',label='prolate p='+np.str(np.round(pearsonr(dm.q[mpro2D],s200.q[mpro2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < R200)$')
plt.savefig(path_plots+'q200_.png')


### theta plots

plt.figure()
plot_binned(lM,t30_3D,'all','k','')
plot_binned(lM[mobl],t30_3D[mobl],'old','C2','--')
plot_binned(lM[mpro],t30_3D[mpro],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 30kpc)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_30.png')

plt.figure()
plot_binned(lM,t50_3D,'all','k','')
plot_binned(lM[mobl],t50_3D[mobl],'old','C2','--')
plot_binned(lM[mpro],t50_3D[mpro],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 50kpc)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_50.png')

plt.figure()
plot_binned(lM,t100_3D,'all','k','')
plot_binned(lM[mobl],t100_3D[mobl],'old','C2','--')
plot_binned(lM[mpro],t100_3D[mpro],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 0.1R500)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_100.png')

plt.figure()
plot_binned(lM,t1000_3D,'all','k','')
plot_binned(lM[mobl],t1000_3D[mobl],'old','C2','--')
plot_binned(lM[mpro],t1000_3D[mpro],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R1000)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_1000.png')


plt.figure()
plot_binned(lM,t500_3D,'all','k','')
plot_binned(lM[mobl],t500_3D[mobl],'old','C2','--')
plot_binned(lM[mpro],t500_3D[mpro],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R500)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_500.png')

plt.figure()
plot_binned(lM,t200_3D,'all','k','')
plot_binned(lM[mobl],t200_3D[mobl],'old','C2','--')
plot_binned(lM[mpro],t200_3D[mpro],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R200)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_200.png')


plt.figure()
plot_binned(lMp,t30_2D,'all','k','')
plot_binned(lMp[mobl2D],t30_2D[mobl2D],'old','C2','--')
plot_binned(lMp[mpro2D],t30_2D[mpro2D],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 30kpc)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_30.png')

plt.figure()
plot_binned(lMp,t50_2D,'all','k','')
plot_binned(lMp[mobl2D],t50_2D[mobl2D],'old','C2','--')
plot_binned(lMp[mpro2D],t50_2D[mpro2D],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 50kpc)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_50.png')

plt.figure()
plot_binned(lMp,t100_2D,'all','k','')
plot_binned(lMp[mobl2D],t100_2D[mobl2D],'old','C2','--')
plot_binned(lMp[mpro2D],t100_2D[mpro2D],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 0.1R500)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_100.png')

plt.figure()
plot_binned(lMp,t1000_2D,'all','k','')
plot_binned(lMp[mobl2D],t1000_2D[mobl2D],'old','C2','--')
plot_binned(lMp[mpro2D],t1000_2D[mpro2D],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R1000)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_1000.png')


plt.figure()
plot_binned(lMp,t500_2D,'all','k','')
plot_binned(lMp[mobl2D],t500_2D[mobl2D],'old','C2','--')
plot_binned(lMp[mpro2D],t500_2D[mpro2D],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R500)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_500.png')

plt.figure()
plot_binned(lMp,t200_2D,'all','k','')
plot_binned(lMp[mobl2D],t200_2D[mobl2D],'old','C2','--')
plot_binned(lMp[mpro2D],t200_2D[mpro2D],'new','C1',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R200)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_200.png')
