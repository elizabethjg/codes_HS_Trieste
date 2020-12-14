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

path_plots = '../plots/age_plots/gx_'+indicator+'/'
os.system('mkdir '+path_plots)

gral  = np.loadtxt(path+'gral_091_2.dat').T
lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)

dm = DarkMatter()

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
# WITH GXs

# shape_parameters

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

# ANGLES WITH DARK MATTER

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

# STACKED SHAPES

R1000 = gral[4]
R500  = gral[5]
R200  = gral[6]

R = np.vstack((R1000,R500,R200)).T
Rp = np.array((R.tolist())*3)

R  = R[mN1000x]
Rp = Rp[mn1000x]


S = np.vstack((g1000.S,g500.S,g200.S)).T[mN1000x]
T = np.vstack((g1000.T,g500.T,g200.T)).T[mN1000x]
q = np.vstack((g1000.q,g500.q,g200.q)).T[mn1000x]

Sx = np.vstack((g1000x.S,g500x.S,g200x.S)).T[mN1000x]
Tx = np.vstack((g1000x.T,g500x.T,g200x.T)).T[mN1000x]
qx = np.vstack((g1000x.q,g500x.q,g200x.q)).T[mn1000x]


t3D = np.vstack((t1000_3D,t500_3D,t200_3D)).T[mN1000x]
t2D = np.vstack((t1000_2D,t500_2D,t200_2D)).T[mn1000x]
t3Dx = np.vstack((t1000x_3D,t500x_3D,t200x_3D)).T[mN1000x]
t2Dx = np.vstack((t1000x_2D,t500x_2D,t200x_2D)).T[mn1000x]


####### S plots

plt.figure()
plt.plot(dm.S[mold],g200.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold],g200.S[mold])[0],2)))
plt.plot(dm.S[mnew],g200.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew],g200.S[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{gx} (R < R200)$')
plt.savefig(path_plots+'S200_'+indicator+'.png')

plt.figure()
plt.plot(dm.S[mold],g500.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold],g500.S[mold])[0],2)))
plt.plot(dm.S[mnew],g500.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew],g500.S[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{gx} (R < R500)$')
plt.savefig(path_plots+'S500_'+indicator+'.png')

plt.figure()
plt.plot(dm.S[mold],g1000.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mN1000*mold],g1000.S[mN1000*mold])[0],2)))
plt.plot(dm.S[mnew],g1000.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mN1000*mnew],g1000.S[mN1000*mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S_{gx} (R < R1000)$')
plt.savefig(path_plots+'S1000_'+indicator+'.png')


plt.figure()
plt.plot(dm.S[mold],g200x.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold],g200x.S[mold])[0],2)))
plt.plot(dm.S[mnew],g200x.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew],g200x.S[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S^x_{gx} (R < R200)$')
plt.savefig(path_plots+'S200x_'+indicator+'.png')

plt.figure()
plt.plot(dm.S[mold],g500x.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold*mN500x],g500x.S[mold*mN500x])[0],2)))
plt.plot(dm.S[mnew],g500x.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew*mN500x],g500x.S[mnew*mN500x])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S^x_{gx} (R < R500)$')
plt.savefig(path_plots+'S500x_'+indicator+'.png')

plt.figure()
plt.plot(dm.S[mold],g1000x.S[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.S[mold*mN1000x],g1000x.S[mold*mN1000x])[0],2)))
plt.plot(dm.S[mnew],g1000x.S[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.S[mnew*mN1000x],g1000x.S[mnew*mN1000x])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$S_{DM}$')
plt.ylabel('$S^x_{gx} (R < R1000)$')
plt.savefig(path_plots+'S1000x_'+indicator+'.png')


plt.figure()    
plt.plot(np.median(np.log10(R),axis=0),np.median(S,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(S,axis=0)+np.std(S,axis=0),np.median(S,axis=0)-np.std(S,axis=0),'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mnew[mN1000x]]),axis=0),np.median(S[mnew[mN1000x]],axis=0),'C0')
plt.fill_between(np.median(np.log10(R[mnew[mN1000x]]),axis=0),np.median(S[mnew[mN1000x]],axis=0)+np.std(S[mnew[mN1000x]],axis=0),np.median(S[mnew[mN1000x]],axis=0)-np.std(S[mnew[mN1000x]],axis=0),'C0',alpha=0.1)
plt.plot(np.median(np.log10(R[mold[mN1000x]]),axis=0),np.median(S[mold[mN1000x]],axis=0),'C3')
plt.fill_between(np.median(np.log10(R[mold[mN1000x]]),axis=0),np.median(S[mold[mN1000x]],axis=0)+np.std(S[mold[mN1000x]],axis=0),np.median(S[mold[mN1000x]],axis=0)-np.std(S[mnew[mN1000x]],axis=0),'C3',alpha=0.1)
plt.ylabel('$S$')
plt.xticks(np.median(np.log10(R),axis=0),['R1000','R500','R200'])
plt.axis([2.8,np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'S_R'+indicator+'.png')

####### T plots

plt.figure()
plt.plot(dm.T[mold],g200.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold],g200.T[mold])[0],2)))
plt.plot(dm.T[mnew],g200.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew],g200.T[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{gx} (R < R200)$')
plt.savefig(path_plots+'T200_'+indicator+'.png')

plt.figure()
plt.plot(dm.T[mold],g500.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold],g500.T[mold])[0],2)))
plt.plot(dm.T[mnew],g500.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew],g500.T[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{gx} (R < R500)$')
plt.savefig(path_plots+'T500_'+indicator+'.png')

plt.figure()
plt.plot(dm.T[mold],g1000.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mN1000*mold],g1000.T[mN1000*mold])[0],2)))
plt.plot(dm.T[mnew],g1000.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mN1000*mnew],g1000.T[mN1000*mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T_{gx} (R < R1000)$')
plt.savefig(path_plots+'T1000_'+indicator+'.png')


plt.figure()
plt.plot(dm.T[mold],g200x.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold],g200x.T[mold])[0],2)))
plt.plot(dm.T[mnew],g200x.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew],g200x.T[mnew])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T^x_{gx} (R < R200)$')
plt.savefig(path_plots+'T200x_'+indicator+'.png')

plt.figure()
plt.plot(dm.T[mold],g500x.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold*mN500x],g500x.T[mold*mN500x])[0],2)))
plt.plot(dm.T[mnew],g500x.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew*mN500x],g500x.T[mnew*mN500x])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T^x_{gx} (R < R500)$')
plt.savefig(path_plots+'T500x_'+indicator+'.png')

plt.figure()
plt.plot(dm.T[mold],g1000x.T[mold],'C3^',label='older p='+np.str(np.round(pearsonr(dm.T[mold*mN1000x],g1000x.T[mold*mN1000x])[0],2)))
plt.plot(dm.T[mnew],g1000x.T[mnew],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.T[mnew*mN1000x],g1000x.T[mnew*mN1000x])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$T_{DM}$')
plt.ylabel('$T^x_{gx} (R < R1000)$')
plt.savefig(path_plots+'T1000x_'+indicator+'.png')


####### q plots

plt.figure()
plt.plot(dm.q[mold2D],g200.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D],g200.q[mold2D])[0],2)))
plt.plot(dm.q[mnew2D],g200.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D],g200.q[mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{gx} (R < R200)$')
plt.savefig(path_plots+'q200_'+indicator+'.png')

plt.figure()
plt.plot(dm.q[mold2D],g500.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D],g500.q[mold2D])[0],2)))
plt.plot(dm.q[mnew2D],g500.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D],g500.q[mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{gx} (R < R500)$')
plt.savefig(path_plots+'q500_'+indicator+'.png')

plt.figure()
plt.plot(dm.q[mold2D],g1000.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mn1000*mold2D],g1000.q[mn1000*mold2D])[0],2)))
plt.plot(dm.q[mnew2D],g1000.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mn1000*mnew2D],g1000.q[mn1000*mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q_{gx} (R < R1000)$')
plt.savefig(path_plots+'q1000_'+indicator+'.png')


plt.figure()
plt.plot(dm.q[mold2D],g200x.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D],g200x.q[mold2D])[0],2)))
plt.plot(dm.q[mnew2D],g200x.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D],g200x.q[mnew2D])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q^x_{gx} (R < R200)$')
plt.savefig(path_plots+'q200x_'+indicator+'.png')

plt.figure()
plt.plot(dm.q[mold2D],g500x.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D*mn500x],g500x.q[mold2D*mn500x])[0],2)))
plt.plot(dm.q[mnew2D],g500x.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D*mn500x],g500x.q[mnew2D*mn500x])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q^x_{gx} (R < R500)$')
plt.savefig(path_plots+'q500x_'+indicator+'.png')

plt.figure()
plt.plot(dm.q[mold2D],g1000x.q[mold2D],'C3^',label='older p='+np.str(np.round(pearsonr(dm.q[mold2D*mn1000x],g1000x.q[mold2D*mn1000x])[0],2)))
plt.plot(dm.q[mnew2D],g1000x.q[mnew2D],'C0^',label='recent p='+np.str(np.round(pearsonr(dm.q[mnew2D*mn1000x],g1000x.q[mnew2D*mn1000x])[0],2)))
plt.plot([0.3,0.85],[0.3,0.85],'C7--')
plt.legend()
plt.xlabel('$q_{DM}$')
plt.ylabel('$q^x_{gx} (R < R1000)$')
plt.savefig(path_plots+'q1000x_'+indicator+'.png')


# theta plots

plt.figure()
plot_binned(lM[mN1000],t1000_3D[mN1000],'all','k','',nbins=5)
plot_binned(lM[mN1000*mold],t1000_3D[mN1000*mold],'old','C3','',nbins=5)
plot_binned(lM[mN1000*mnew],t1000_3D[mN1000*mnew],'new','C0','',nbins=5)
plot_binned(lM[mN1000x],t1000_3D[mN1000x],'gal_cut','k','--',nbins=5)
plt.legend(loc = 1)
plot_binned(lM[mN1000x*mold],t1000x_3D[mN1000x*mold],'old','C3','--',nbins=5)
plot_binned(lM[mN1000x*mnew],t1000x_3D[mN1000x*mnew],'new','C0','--',nbins=5)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R1000)$')
plt.savefig(path_plots+'theta3D_1000.png')

plt.figure()
plot_binned(lM[mN500],t500_3D[mN500],'all','k','',nbins=5)
plot_binned(lM[mN500*mold],t500_3D[mN500*mold],'old','C3','',nbins=5)
plot_binned(lM[mN500*mnew],t500_3D[mN500*mnew],'new','C0','',nbins=5)
plot_binned(lM[mN500x],t500_3D[mN500x],'gal_cut','k','--',nbins=5)
plt.legend(loc = 1)
plot_binned(lM[mN500x*mold],t500x_3D[mN500x*mold],'old','C3','--',nbins=5)
plot_binned(lM[mN500x*mnew],t500x_3D[mN500x*mnew],'new','C0','--',nbins=5)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R500)$')
plt.savefig(path_plots+'theta3D_500.png')

plt.figure()
plot_binned(lM[mN200],t200_3D[mN200],'all','k','',nbins=5)
plot_binned(lM[mN200*mold],t200_3D[mN200*mold],'old','C3','',nbins=5)
plot_binned(lM[mN200*mnew],t200_3D[mN200*mnew],'new','C0','',nbins=5)
plot_binned(lM[mN200x],t200_3D[mN200x],'gal_cut','k','--',nbins=5)
plt.legend(loc = 1)
plot_binned(lM[mN200x*mold],t200x_3D[mN200x*mold],'old','C3','--',nbins=5)
plot_binned(lM[mN200x*mnew],t200x_3D[mN200x*mnew],'new','C0','--',nbins=5)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta (R < R200)$')
plt.savefig(path_plots+'theta3D_200.png')

plt.figure()
plot_binned(lMp[mn1000],t1000_2D[mn1000],'all','k','',nbins=5)
plot_binned(lMp[mn1000*mold2D],t1000_2D[mn1000*mold2D],'old','C3','',nbins=5)
plot_binned(lMp[mn1000*mnew2D],t1000_2D[mn1000*mnew2D],'new','C0','',nbins=5)
plot_binned(lMp[mn1000x],t1000_2D[mn1000x],'gal_cut','k','--',nbins=5)
plt.legend(loc = 1)
plot_binned(lMp[mn1000x*mold2D],t1000x_2D[mn1000x*mold2D],'old','C3','--',nbins=5)
plot_binned(lMp[mn1000x*mnew2D],t1000x_2D[mn1000x*mnew2D],'new','C0','--',nbins=5)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{2D} (R < R1000)$')
plt.savefig(path_plots+'theta2D_1000.png')

plt.figure()
plot_binned(lMp[mn500],t500_2D[mn500],'all','k','',nbins=5)
plot_binned(lMp[mn500*mold2D],t500_2D[mn500*mold2D],'old','C3','',nbins=5)
plot_binned(lMp[mn500*mnew2D],t500_2D[mn500*mnew2D],'new','C0','',nbins=5)
plot_binned(lMp[mn500x],t500_2D[mn500x],'gal_cut','k','--',nbins=5)
plt.legend(loc = 1)
plot_binned(lMp[mn500x*mold2D],t500x_2D[mn500x*mold2D],'old','C3','--',nbins=5)
plot_binned(lMp[mn500x*mnew2D],t500x_2D[mn500x*mnew2D],'new','C0','--',nbins=5)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{2D} (R < R500)$')
plt.savefig(path_plots+'theta2D_500.png')

plt.figure()
plot_binned(lMp[mn200],t200_2D[mn200],'all','k','',nbins=5)
plot_binned(lMp[mn200*mold2D],t200_2D[mn200*mold2D],'old','C3','',nbins=5)
plot_binned(lMp[mn200*mnew2D],t200_2D[mn200*mnew2D],'new','C0','',nbins=5)
plot_binned(lMp[mn200x],t200_2D[mn200x],'gal_cut','k','--',nbins=5)
plt.legend(loc = 1)
plot_binned(lMp[mn200x*mold2D],t200x_2D[mn200x*mold2D],'old','C3','--',nbins=5)
plot_binned(lMp[mn200x*mnew2D],t200x_2D[mn200x*mnew2D],'new','C0','--',nbins=5)
plt.ylim([0,50])
plt.xlabel('$\log M_{200}$')
plt.ylabel(r'$\theta_{2D} (R < R200)$')
plt.savefig(path_plots+'theta2D_200.png')

