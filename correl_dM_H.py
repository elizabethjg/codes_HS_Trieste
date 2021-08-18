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
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 14})
plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'





path = '../catalog/nuevosdats/'

gral  = np.loadtxt(path+'gral_nounb_091.dat').T


# gral  = np.loadtxt('../catalog/gral_091_2.dat').T

plotspath = '../final_plots/'

lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)

lM1000 = np.log10(gral[7])
lM1000p = np.array((lM1000.tolist())*3)


off   = gral[13]
off2D = np.concatenate((gral[14],gral[15],gral[16]))
DV    = gral[17]
DV2D  = np.concatenate((gral[18],gral[19],gral[20]))
gap  = gral[21]
gap2D = np.concatenate((gral[22],gral[23],gral[24]))

mn_off,mo_off,mn2d_off,mo2d_off = newold('off')
mn_gap,mo_gap,mn2d_gap,mo2d_gap = newold('gap')
mn_dv ,mo_dv, mn2d_dv ,mo2d_dv  = newold('DV')

N200 = Galaxias(radio=200).N
N500 = Galaxias(radio=500).N
N1000 = Galaxias(radio=1000).N

n200 = Galaxias(radio=200).n
n500 = Galaxias(radio=500).n
n1000 = Galaxias(radio=1000).n


H30  = DarkMatter(30,False)
H50  = DarkMatter(50,False)
H100  = DarkMatter(100,False)
H1000  = DarkMatter(1000,False)
H500  = DarkMatter(500,False)
H200  = DarkMatter(200,False)

DM30 = DarkMatter(30)
DM50 = DarkMatter(50)
DM100 = DarkMatter(100)
DM1000 = DarkMatter(1000)
DM500 = DarkMatter(500)
DM200 = DarkMatter(200)

R1000 = gral[4]
R500  = gral[5]
R200  = gral[6]
R30   = np.ones(len(R200))*30.
R50   = np.ones(len(R200))*50.
R1    = 0.1*R500
R = np.vstack((R30/R200,R50/R200,R1/R200,R1000/R200,R500/R200,R200/R200)).T
Rp = np.array((R.tolist())*3)



ct30_3D  , t30_3D     = cosangle(H30.a3D  ,DM30.a3D) 
ct50_3D  , t50_3D     = cosangle(H50.a3D  ,DM50.a3D)
ct100_3D , t100_3D    = cosangle(H100.a3D ,DM100.a3D)
ct200_3D , t200_3D    = cosangle(H200.a3D ,DM200.a3D)
ct500_3D , t500_3D    = cosangle(H500.a3D ,DM500.a3D)
ct1000_3D, t1000_3D   = cosangle(H1000.a3D,DM1000.a3D)

ct30_2D  , t30_2D     = cosangle(H30.a2D  ,DM30.a2D) 
ct50_2D  , t50_2D     = cosangle(H50.a2D  ,DM50.a2D)
ct100_2D , t100_2D    = cosangle(H100.a2D ,DM100.a2D)
ct200_2D , t200_2D    = cosangle(H200.a2D ,DM200.a2D)
ct500_2D , t500_2D    = cosangle(H500.a2D ,DM500.a2D)
ct1000_2D, t1000_2D   = cosangle(H1000.a2D,DM1000.a2D)

t3D = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
t2D = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T

SH_r = np.vstack((H30.S,H50.S,H100.S,H1000.S,H500.S,H200.S)).T
SDM_r = np.vstack((DM30.S,DM50.S,DM100.S,DM1000.S,DM500.S,DM200.S)).T

TH_r = np.vstack((H30.T,H50.T,H100.T,H1000.T,H500.T,H200.T)).T
TDM_r = np.vstack((DM30.T,DM50.T,DM100.T,DM1000.T,DM500.T,DM200.T)).T

qH_r = np.vstack((H30.q,H50.q,H100.q,H1000.q,H500.q,H200.q)).T
qDM_r = np.vstack((DM30.q,DM50.q,DM100.q,DM1000.q,DM500.q,DM200.q)).T

limites = [0.02,1.,0.88,1.12]

Rlegend = np.array(['30kpc','50kpc','0.1R$_{500}$','R$_{1000}$','R$_{500}$','R$_{200}$'])

f, ax = plt.subplots(3,1, figsize=(5,12), sharex=True,sharey=True)
f.subplots_adjust(hspace=0,wspace=0)
for ax2 in ax.flatten():
    [ax2.axvline(x, ls='--', color='k',lw=0.5,alpha=0.5) for x in np.median(R,axis=0)]

# S plot


plotR_ind(R,SDM_r/SH_r,'k','all',ax = ax[1])
plotR_ind(R[mn_gap],(SDM_r/SH_r)[mn_gap],'C0','non-relaxed',style='',ax = ax[1])
plotR_ind(R[mo_gap],(SDM_r/SH_r)[mo_gap],'sienna','relaxed',style='',ax = ax[1])

ax[1].legend(frameon=False,loc=3)

# plotR_ind(R[mn_off],(SDM_r/SH_r)[mn_off],'C0','all',style='--',ax = ax[1])
# plotR_ind(R[mn_dv] ,(SDM_r/SH_r)[mn_dv] ,'C0','all',style=':',ax = ax[1])


# plotR_ind(R[mo_off],(SDM_r/SH_r)[mo_off],'sienna','all',style='--',ax = ax[1])
# plotR_ind(R[mo_dv] ,(SDM_r/SH_r)[mo_dv] ,'sienna','all',style=':',ax = ax[1])

ax[0].plot([R.min(),1],[1,1],'C7')
ax[0].set_xlim([0.02,1.])
ax[0].set_ylim([0.88,1.12])

ax[1].set_ylabel('$S_{DM}/S_H$')
ax[2].set_xlabel('$R/R_{200}$')

ax[0].set_xscale('log')

# T plot

plotR_ind(R,TDM_r/TH_r,'k','all',ax = ax[0])
plotR_ind(R[mn_gap],(TDM_r/TH_r)[mn_gap],'C0','non-relaxed',style='',ax = ax[0])
# plotR_ind(R[mn_off],(TDM_r/TH_r)[mn_off],'C0','all',style='--',ax = ax[0])
# plotR_ind(R[mn_dv] ,(TDM_r/TH_r)[mn_dv] ,'C0','all',style=':',ax = ax[0])
                           
plotR_ind(R[mo_gap],(TDM_r/TH_r)[mo_gap],'sienna','relaxed',style='',ax = ax[0])
# plotR_ind(R[mo_off],(TDM_r/TH_r)[mo_off],'sienna','all',style='--',ax = ax[0])
# plotR_ind(R[mo_dv] ,(TDM_r/TH_r)[mo_dv] ,'sienna','all',style=':',ax = ax[0])

# plt.xticks(np.median(np.log10(R),axis=0),Rlegend,fontsize=10.5)
ax[1].plot([0,5],[1,1],'C7')


ax[0].set_ylabel('$T_{DM}/T_H$')


# q plot

plotR_ind(Rp,qDM_r/qH_r,'k','all',ax = ax[2])
plotR_ind(Rp[mn2d_gap],(qDM_r/qH_r)[mn2d_gap],'C0','non-relaxed',style='',ax = ax[2])
# plotR_ind(Rp[mn2d_off],(qDM_r/qH_r)[mn2d_off],'C0','all',style='--',ax = ax[2])
# plotR_ind(Rp[mn2d_dv] ,(qDM_r/qH_r)[mn2d_dv] ,'C0','all',style=':',ax = ax[2])      
plotR_ind(Rp[mo2d_gap],(qDM_r/qH_r)[mo2d_gap],'sienna','relaxed',style='',ax = ax[2])
# plotR_ind(Rp[mo2d_off],(qDM_r/qH_r)[mo2d_off],'sienna','all',style='--',ax = ax[2])
# plotR_ind(Rp[mo2d_dv] ,(qDM_r/qH_r)[mo2d_dv] ,'sienna','all',style=':',ax = ax[2])


ax[2].plot([0,5],[1,1],'C7')



ax[2].set_ylabel('$q_{DM}/q_H$')

plt.savefig(plotspath+'shape_DM_H.pdf',bbox_inches='tight')


# t3d plot

limites = [0.02,1.,0.,10.]

f, ax = plt.subplots(2,1, figsize=(5,6), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

for ax2 in ax.flatten():
    [ax2.axvline(x, ls='--', color='k',lw=0.5,alpha=0.5) for x in np.median(R,axis=0)]


plotR_ind(R,t3D,'k','all',ax = ax[0])
plotR_ind(R[mn_gap],t3D[mn_gap],'C0','non-relaxed',style='',ax = ax[0])
plotR_ind(R[mo_gap],t3D[mo_gap],'sienna','relaxed',style='',ax = ax[0])

ax[0].legend(frameon=False,loc=2)

# plotR_ind(R[mn_off],t3D[mn_off],'C0','all',style='--',ax = ax[0])
# plotR_ind(R[mn_dv] ,t3D[mn_dv] ,'C0','all',style=':',ax = ax[0])

# plotR_ind(R[mo_off],t3D[mo_off],'sienna','all',style='--',ax = ax[0])
# plotR_ind(R[mo_dv] ,t3D[mo_dv] ,'sienna','all',style=':',ax = ax[0])


ax[0].axis(limites)


ax[0].set_ylabel(r'$\theta^{3D} [\circ]$')

# t2d plot


plotR_ind(Rp,t2D,'k','all',ax = ax[1])
plotR_ind(Rp[mn2d_gap],t2D[mn2d_gap],'C0','non-relaxed',style='',ax = ax[1])
# plotR_ind(Rp[mn2d_off],t2D[mn2d_off],'C0','all',style='--',ax = ax[1])
# plotR_ind(Rp[mn2d_dv] ,t2D[mn2d_dv] ,'C0','all',style=':',ax = ax[1])
plotR_ind(Rp[mo2d_gap],t2D[mo2d_gap],'sienna','relaxed',style='',ax = ax[1])
# plotR_ind(Rp[mo2d_off],t2D[mo2d_off],'sienna','all',style='--',ax = ax[1])
# plotR_ind(Rp[mo2d_dv] ,t2D[mo2d_dv] ,'sienna','all',style=':',ax = ax[1])


# ax[1].set_xticks(np.median(np.log10(R),axis=0))
# ax[1].set_xticklabels(Rlegend,fontsize=10.5)


ax[1].axis(limites)

ax[1].set_ylabel(r'$\theta [\circ]$')
ax[1].set_xlabel('$R/R_{200}$')

ax[1].set_xscale('log')

plt.savefig(plotspath+'theta_DM_H.pdf',bbox_inches='tight')
