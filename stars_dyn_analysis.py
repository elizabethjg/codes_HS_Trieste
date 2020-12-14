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

path_plots = '../plots/age_plots/stars_'+indicator+'/'
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
for j in range(mnew.sum()):
        plt.plot(np.log10(R[mnew][j]),S[mnew][j],'C0',alpha=0.5)
for j in range(mold.sum()):
        plt.plot(np.log10(R[mold][j]),S[mold][j],'C3',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'S_R'+indicator+'_ind.png')

plt.figure()
for j in range(mnew.sum()):
        plt.plot(np.log10(R[mnew][j]),T[mnew][j],'C0',alpha=0.5)
for j in range(mold.sum()):
        plt.plot(np.log10(R[mold][j]),T[mold][j],'C3',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'T_R'+indicator+'_ind.png')

plt.figure()
for j in range(mnew2D.sum()):
        plt.plot(np.log10(Rp[mnew2D][j]),q[mnew2D][j],'C0',alpha=0.5)
for j in range(mold2D.sum()):
        plt.plot(np.log10(Rp[mold2D][j]),q[mold2D][j],'C3',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.,0.9])
plt.savefig(path_plots+'q_R'+indicator+'_ind.png')


plt.figure()
for j in range(mnew.sum()):
        plt.plot(np.log10(R[mnew][j]),t3D[mnew][j],'C0',alpha=0.5)
for j in range(mold.sum()):
        plt.plot(np.log10(R[mold][j]),t3D[mold][j],'C3',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.,90])
plt.savefig(path_plots+'t3D_R'+indicator+'_ind.png')


plt.figure()
for j in range(mnew2D.sum()):
        plt.plot(np.log10(Rp[mnew2D][j]),t2D[mnew2D][j],'C0',alpha=0.5)
for j in range(mold2D.sum()):
        plt.plot(np.log10(Rp[mold2D][j]),t2D[mold2D][j],'C3',alpha=0.5)
plt.axis([np.log10(30),np.log10(2000),0.,90])
plt.savefig(path_plots+'t2D_R'+indicator+'_ind.png')


plt.figure()    
plt.plot(np.median(np.log10(R),axis=0),np.median(S,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(S,axis=0)+np.std(S,axis=0),np.median(S,axis=0)-np.std(S,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mnew]),axis=0),np.median(S[mnew],axis=0),'C0')
plt.fill_between(np.median(np.log10(R[mnew]),axis=0),np.median(S[mnew],axis=0)+np.std(S[mnew],axis=0),np.median(S[mnew],axis=0)-np.std(S[mnew],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(R[mold]),axis=0),np.median(S[mold],axis=0),'C3')
plt.fill_between(np.median(np.log10(R[mold]),axis=0),np.median(S[mold],axis=0)+np.std(S[mold],axis=0),np.median(S[mold],axis=0)-np.std(S[mold],axis=0),color = 'C3',alpha=0.1)
plt.ylabel('$S$')
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'S_R'+indicator+'.png')

plt.figure()
plt.plot(np.median(np.log10(R),axis=0),np.median(T,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(T,axis=0)+np.std(T,axis=0),np.median(T,axis=0)-np.std(T,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mnew]),axis=0),np.median(T[mnew],axis=0),'C0')
plt.fill_between(np.median(np.log10(R[mnew]),axis=0),np.median(T[mnew],axis=0)+np.std(T[mnew],axis=0),np.median(T[mnew],axis=0)-np.std(T[mnew],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(R[mold]),axis=0),np.median(T[mold],axis=0),'C3')
plt.fill_between(np.median(np.log10(R[mold]),axis=0),np.median(T[mold],axis=0)+np.std(T[mold],axis=0),np.median(T[mold],axis=0)-np.std(T[mold],axis=0),color = 'C3',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel('$T$')
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'T_R'+indicator+'.png')

plt.figure()
plt.plot(np.median(np.log10(Rp),axis=0),np.median(q,axis=0),'k')
plt.fill_between(np.median(np.log10(Rp),axis=0),np.median(q,axis=0)+np.std(q,axis=0),np.median(q,axis=0)-np.std(q,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mnew2D]),axis=0),np.median(q[mnew2D],axis=0),'C0')
plt.fill_between(np.median(np.log10(Rp[mnew2D]),axis=0),np.median(q[mnew2D],axis=0)+np.std(q[mnew2D],axis=0),np.median(q[mnew2D],axis=0)-np.std(q[mnew2D],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mold2D]),axis=0),np.median(q[mold2D],axis=0),'C3')
plt.fill_between(np.median(np.log10(Rp[mold2D]),axis=0),np.median(q[mold2D],axis=0)+np.std(q[mold2D],axis=0),np.median(q[mold2D],axis=0)-np.std(q[mold2D],axis=0),color = 'C3',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel('$q$')
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'q_Rp'+indicator+'.png')

plt.figure()
plt.plot(np.median(np.log10(R),axis=0),np.median(t3D,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(t3D,axis=0)+np.std(t3D,axis=0),np.median(t3D,axis=0)-np.std(t3D,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mnew]),axis=0),np.median(t3D[mnew],axis=0),'C0')
plt.fill_between(np.median(np.log10(R[mnew]),axis=0),np.median(t3D[mnew],axis=0)+np.std(t3D[mnew],axis=0),np.median(t3D[mnew],axis=0)-np.std(t3D[mnew],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(R[mold]),axis=0),np.median(t3D[mold],axis=0),'C3')
plt.fill_between(np.median(np.log10(R[mold]),axis=0),np.median(t3D[mold],axis=0)+np.std(t3D[mold],axis=0),np.median(t3D[mold],axis=0)-np.std(t3D[mold],axis=0),color = 'C3',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel(r'$\theta$')
plt.axis([np.log10(30),np.log10(2000),0.,50.])
plt.savefig(path_plots+'t3D_R'+indicator+'.png')

plt.figure()
plt.plot(np.median(np.log10(Rp),axis=0),np.median(t2D,axis=0),'k')
plt.fill_between(np.median(np.log10(Rp),axis=0),np.median(t2D,axis=0)+np.std(t2D,axis=0),np.median(t2D,axis=0)-np.std(t2D,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mnew2D]),axis=0),np.median(t2D[mnew2D],axis=0),'C0')
plt.fill_between(np.median(np.log10(Rp[mnew2D]),axis=0),np.median(t2D[mnew2D],axis=0)+np.std(t2D[mnew2D],axis=0),np.median(t2D[mnew2D],axis=0)-np.std(t2D[mnew2D],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mold2D]),axis=0),np.median(t2D[mold2D],axis=0),'C3')
plt.fill_between(np.median(np.log10(Rp[mold2D]),axis=0),np.median(t2D[mold2D],axis=0)+np.std(t2D[mold2D],axis=0),np.median(t2D[mold2D],axis=0)-np.std(t2D[mold2D],axis=0),color = 'C3',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel(r'$\theta_{2D}$')
plt.axis([np.log10(30),np.log10(2000),0.,50.])
plt.savefig(path_plots+'t2D_Rp'+indicator+'.png')


############ S plots
plt.scatter(dm.S[mold],s30.S[mold],c=lM[mold])

plt.figure()
plt.plot(dm.S[mold],s30.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold],s30.S[mold])[0],2)))
plt.plot(dm.S[mnew],s30.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew],s30.S[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < 30kpc)$')
plt.savefig(path_plots+'S30_'+indicator+'.png')

plt.figure()
plt.plot(dm.S[mold],s50.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold],s50.S[mold])[0],2)))
plt.plot(dm.S[mnew],s50.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew],s50.S[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < 50kpc)$')
plt.savefig(path_plots+'S50_'+indicator+'.png')

plt.figure()
plt.plot(dm.S[mold],s100.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold],s100.S[mold])[0],2)))
plt.plot(dm.S[mnew],s100.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew],s100.S[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < 0.1R500)$')
plt.savefig(path_plots+'S100_'+indicator+'.png')

plt.figure()
plt.plot(dm.S[mold],s1000.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold],s1000.S[mold])[0],2)))
plt.plot(dm.S[mnew],s1000.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew],s1000.S[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < R1000)$')
plt.savefig(path_plots+'S1000_'+indicator+'.png')

plt.figure()
plt.plot(dm.S[mold],s500.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold],s500.S[mold])[0],2)))
plt.plot(dm.S[mnew],s500.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew],s500.S[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < R500)$')
plt.savefig(path_plots+'S500_'+indicator+'.png')

plt.figure()
plt.plot(dm.S[mold],s200.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold],s200.S[mold])[0],2)))
plt.plot(dm.S[mnew],s200.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew],s200.S[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{*} (R < R200)$')
plt.savefig(path_plots+'S200_'+indicator+'.png')


############ T plots

plt.figure()
plt.plot(dm.T[mold],s30.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold],s30.T[mold])[0],2)))
plt.plot(dm.T[mnew],s30.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew],s30.T[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < 30kpc)$')
plt.savefig(path_plots+'T30_'+indicator+'.png')

plt.figure()
plt.plot(dm.T[mold],s50.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold],s50.T[mold])[0],2)))
plt.plot(dm.T[mnew],s50.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew],s50.T[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < 50kpc)$')
plt.savefig(path_plots+'T50_'+indicator+'.png')

plt.figure()
plt.plot(dm.T[mold],s100.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold],s100.T[mold])[0],2)))
plt.plot(dm.T[mnew],s100.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew],s100.T[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < 0.1R500)$')
plt.savefig(path_plots+'T100_'+indicator+'.png')

plt.figure()
plt.plot(dm.T[mold],s1000.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold],s1000.T[mold])[0],2)))
plt.plot(dm.T[mnew],s1000.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew],s1000.T[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < R1000)$')
plt.savefig(path_plots+'T1000_'+indicator+'.png')

plt.figure()
plt.plot(dm.T[mold],s500.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold],s500.T[mold])[0],2)))
plt.plot(dm.T[mnew],s500.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew],s500.T[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < R500)$')
plt.savefig(path_plots+'T500_'+indicator+'.png')

plt.figure()
plt.plot(dm.T[mold],s200.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold],s200.T[mold])[0],2)))
plt.plot(dm.T[mnew],s200.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew],s200.T[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{*} (R < R200)$')
plt.savefig(path_plots+'T200_'+indicator+'.png')

############ q plots

plt.figure()
plt.plot(dm.q[mold2D],s30.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D],s30.q[mold2D])[0],2)))
plt.plot(dm.q[mnew2D],s30.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D],s30.q[mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < 30kpc)$')
plt.savefig(path_plots+'q30_'+indicator+'.png')

plt.figure()
plt.plot(dm.q[mold2D],s50.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D],s50.q[mold2D])[0],2)))
plt.plot(dm.q[mnew2D],s50.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D],s50.q[mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < 50kpc)$')
plt.savefig(path_plots+'q50_'+indicator+'.png')

plt.figure()
plt.plot(dm.q[mold2D],s100.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D],s100.q[mold2D])[0],2)))
plt.plot(dm.q[mnew2D],s100.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D],s100.q[mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < 0.1R500)$')
plt.savefig(path_plots+'q100_'+indicator+'.png')

plt.figure()
plt.plot(dm.q[mold2D],s1000.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D],s1000.q[mold2D])[0],2)))
plt.plot(dm.q[mnew2D],s1000.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D],s1000.q[mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < R1000)$')
plt.savefig(path_plots+'q1000_'+indicator+'.png')

plt.figure()
plt.plot(dm.q[mold2D],s500.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D],s500.q[mold2D])[0],2)))
plt.plot(dm.q[mnew2D],s500.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D],s500.q[mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < R500)$')
plt.savefig(path_plots+'q500_'+indicator+'.png')

plt.figure()
plt.plot(dm.q[mold2D],s200.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D],s200.q[mold2D])[0],2)))
plt.plot(dm.q[mnew2D],s200.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D],s200.q[mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{*} (R < R200)$')
plt.savefig(path_plots+'q200_'+indicator+'.png')


### theta plots

plt.figure()
plot_binned(lM,t30_3D,'all','k','')
plot_binned(lM[mold],t30_3D[mold],'old','C3','--')
plot_binned(lM[mnew],t30_3D[mnew],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 30kpc)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_30.png')

plt.figure()
plot_binned(lM,t50_3D,'all','k','')
plot_binned(lM[mold],t50_3D[mold],'old','C3','--')
plot_binned(lM[mnew],t50_3D[mnew],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 50kpc)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_50.png')

plt.figure()
plot_binned(lM,t100_3D,'all','k','')
plot_binned(lM[mold],t100_3D[mold],'old','C3','--')
plot_binned(lM[mnew],t100_3D[mnew],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 0.1R500)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_100.png')

plt.figure()
plot_binned(lM,t1000_3D,'all','k','')
plot_binned(lM[mold],t1000_3D[mold],'old','C3','--')
plot_binned(lM[mnew],t1000_3D[mnew],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R1000)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_1000.png')


plt.figure()
plot_binned(lM,t500_3D,'all','k','')
plot_binned(lM[mold],t500_3D[mold],'old','C3','--')
plot_binned(lM[mnew],t500_3D[mnew],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R500)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_500.png')

plt.figure()
plot_binned(lM,t200_3D,'all','k','')
plot_binned(lM[mold],t200_3D[mold],'old','C3','--')
plot_binned(lM[mnew],t200_3D[mnew],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R200)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta3D_200.png')


plt.figure()
plot_binned(lMp,t30_2D,'all','k','')
plot_binned(lMp[mold2D],t30_2D[mold2D],'old','C3','--')
plot_binned(lMp[mnew2D],t30_2D[mnew2D],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 30kpc)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_30.png')

plt.figure()
plot_binned(lMp,t50_2D,'all','k','')
plot_binned(lMp[mold2D],t50_2D[mold2D],'old','C3','--')
plot_binned(lMp[mnew2D],t50_2D[mnew2D],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 50kpc)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_50.png')

plt.figure()
plot_binned(lMp,t100_2D,'all','k','')
plot_binned(lMp[mold2D],t100_2D[mold2D],'old','C3','--')
plot_binned(lMp[mnew2D],t100_2D[mnew2D],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < 0.1R500)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_100.png')

plt.figure()
plot_binned(lMp,t1000_2D,'all','k','')
plot_binned(lMp[mold2D],t1000_2D[mold2D],'old','C3','--')
plot_binned(lMp[mnew2D],t1000_2D[mnew2D],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R1000)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_1000.png')


plt.figure()
plot_binned(lMp,t500_2D,'all','k','')
plot_binned(lMp[mold2D],t500_2D[mold2D],'old','C3','--')
plot_binned(lMp[mnew2D],t500_2D[mnew2D],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R500)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_500.png')

plt.figure()
plot_binned(lMp,t200_2D,'all','k','')
plot_binned(lMp[mold2D],t200_2D[mold2D],'old','C3','--')
plot_binned(lMp[mnew2D],t200_2D[mnew2D],'new','C0',':')
plt.legend()
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R200)$')
plt.ylim([0,50])
plt.savefig(path_plots+'theta2D_200.png')
