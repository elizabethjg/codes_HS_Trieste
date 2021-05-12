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
import pandas as pd
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



plotspath = '../final_plots/'
mice = pd.read_pickle('../catalog/data_eli_01.pkl')

C = Clusters()

DM30   = DarkMatter(30)
DM50   = DarkMatter(50)
DM100  = DarkMatter(100)
DM1000 = DarkMatter(1000)
DM500  = DarkMatter(500)
DM200  = DarkMatter(200)


SDM_r = np.vstack((DM30.S,DM50.S,DM100.S,DM1000.S,DM500.S,DM200.S)).T
TDM_r = np.vstack((DM30.T,DM50.T,DM100.T,DM1000.T,DM500.T,DM200.T)).T
qDM_r = np.vstack((DM30.q,DM50.q,DM100.q,DM1000.q,DM500.q,DM200.q)).T




mmas = C.lM200 > 14.6
mmasp = C.lM200p > 14.6



mlow = ~mmas
mlowp = ~mmasp



cRL_gap = CorrelR(Stars,DarkMatter,mN = mlow*C.mo_gap, mN2D = mlowp*C.mo2d_gap)
cRM_gap = CorrelR(Stars,DarkMatter,mN = mmas*C.mo_gap, mN2D = mmasp*C.mo2d_gap)
cUL_gap = CorrelR(Stars,DarkMatter,mN = mlow*C.mn_gap, mN2D = mlowp*C.mn2d_gap)
cUM_gap = CorrelR(Stars,DarkMatter,mN = mmas*C.mn_gap, mN2D = mmasp*C.mn2d_gap)
cRL_off = CorrelR(Stars,DarkMatter,mN = mlow*C.mo_off, mN2D = mlowp*C.mo2d_off)
cRM_off = CorrelR(Stars,DarkMatter,mN = mmas*C.mo_off, mN2D = mmasp*C.mo2d_off)
cUL_off = CorrelR(Stars,DarkMatter,mN = mlow*C.mn_off, mN2D = mlowp*C.mn2d_off)
cUM_off = CorrelR(Stars,DarkMatter,mN = mmas*C.mn_off, mN2D = mmasp*C.mn2d_off)
cRL_dv  = CorrelR(Stars,DarkMatter,mN = mlow*C.mo_dv, mN2D = mlowp*C.mo2d_dv)
cRM_dv  = CorrelR(Stars,DarkMatter,mN = mmas*C.mo_dv, mN2D = mmasp*C.mo2d_dv)
cUL_dv  = CorrelR(Stars,DarkMatter,mN = mlow*C.mn_dv, mN2D = mlowp*C.mn2d_dv)
cUM_dv  = CorrelR(Stars,DarkMatter,mN = mmas*C.mn_dv, mN2D = mmasp*C.mn2d_dv)


mo = DM200.T < 1/3.
mt = (DM200.T > 1/3.)*(DM200.T < 2/3.)
mp = DM200.T > 2/3.




limites = [0.02,1.,0.88,1.12]



f, ax = plt.subplots(2,1, figsize=(5,6), sharex=True,sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

for ax2 in ax.flatten():
    [ax2.axvline(x, ls='--', color='k',lw=0.5,alpha=0.5) for x in np.median(C.Rsp,axis=0)]


mmas = C.lM200 > 14.6
mmasp = C.lM200p > 14.6

ax[1].plot(200,200,'k-',label = '$\log(M_{200}/M_\odot) > 14.6$')
ax[1].plot(200,200,'k--',label = '$\log(M_{200}/M_\odot) < 14.6$')


mlow = ~mmas
mlowp = ~mmasp


ax[1].legend(loc=2,frameon=False,fontsize=13)

plotR_ind(C.Rs[mmas*C.mo_gap],cRM_gap.t3D,'sienna','relaxed',ax = ax[0],style='-')
plotR_ind(C.Rs[mmas*C.mn_gap],cUM_gap.t3D,'C0','non-relaxed',ax = ax[0],style='-')
ax[0].legend(loc=2,frameon=False,fontsize=13)
plotR_ind(C.Rs[mlow*C.mo_gap],cRL_gap.t3D,'sienna','all',ax = ax[0],style='--')
plotR_ind(C.Rs[mlow*C.mn_gap],cUL_gap.t3D,'C0','all',ax = ax[0],style='--')

plotR_ind(C.Rsp[mlowp*C.mo2d_gap],cRL_gap.t2D, 'sienna','all',ax = ax[1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_gap],cRM_gap.t2D, 'sienna','all',ax = ax[1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],cUL_gap.t2D, 'C0','all',ax = ax[1],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],cUM_gap.t2D, 'C0','all',ax = ax[1],style='-')                                                                     

ax[0].set_xlim([0.012,1.01])
ax[0].set_xscale('log')

ax[1].set_xlabel('$R/R_{200}$')


ax[0].set_ylabel(r'$\theta^{3D}$')
ax[1].set_ylabel(r'$\theta$')

ax[0].set_ylim([0,30])

plt.savefig(plotspath+'theta_stars.pdf',bbox_inches='tight')


f, ax = plt.subplots(3,2, figsize=(10,10), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

for ax2 in ax.flatten():
    [ax2.axvline(x, ls='--', color='k',lw=0.5,alpha=0.5) for x in np.median(C.Rsp,axis=0)]

plotR_ind(C.Rs[mmas*C.mo_gap],TDM_r[mmas*C.mo_gap],'k','$\log(M_{200}/M_\odot) > 14.6$',ax=ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_gap],TDM_r[mlow*C.mo_gap],'k','$\log(M_{200}/M_\odot) < 14.6$',ax=ax[0,0],style='--')
ax[0,0].legend(loc=3,frameon=False)
plotR_ind(C.Rs[mmas*C.mo_gap],cRM_gap.Tr,'teal','Stellar',ax=ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_gap],cRL_gap.Tr,'teal','all',ax=ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],TDM_r[mmas*C.mn_gap],'k','Dark matter',ax=ax[0,1],style='-')
plotR_ind(C.Rs[mmas*C.mn_gap],cUM_gap.Tr,'teal','Stellar',ax=ax[0,1],style='-')
ax[0,1].legend(loc=3,frameon=False)
plotR_ind(C.Rs[mlow*C.mn_gap],TDM_r[mlow*C.mn_gap],'k','all',ax=ax[0,1],style='--')
plotR_ind(C.Rs[mlow*C.mn_gap],cUL_gap.Tr,'teal','all',ax=ax[0,1],style='--')

plotR_ind(C.Rs[mmas*C.mo_gap],SDM_r[mmas*C.mo_gap],'k','Dark matter',ax=ax[1,0],style='-')
plotR_ind(C.Rs[mmas*C.mo_gap],cRM_gap.Sr,'teal','Stellar',ax=ax[1,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_gap],SDM_r[mlow*C.mo_gap],'k','all',ax=ax[1,0],style='--')
plotR_ind(C.Rs[mlow*C.mo_gap],cRL_gap.Sr,'teal','all',ax=ax[1,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],SDM_r[mmas*C.mn_gap],'k','Dark matter',ax=ax[1,1],style='-')
plotR_ind(C.Rs[mmas*C.mn_gap],cUM_gap.Sr,'teal','Stellar',ax=ax[1,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap],SDM_r[mlow*C.mn_gap],'k','all',ax=ax[1,1],style='--')
plotR_ind(C.Rs[mlow*C.mn_gap],cUL_gap.Sr,'teal','all',ax=ax[1,1],style='--')

plotR_ind(C.Rsp[mmasp*C.mo2d_gap],qDM_r[mmasp*C.mo2d_gap],'k','Dark matter',ax=ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_gap],qDM_r[mlowp*C.mo2d_gap],'k','all',ax=ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],qDM_r[mmasp*C.mn2d_gap],'k','Dark matter',ax=ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],qDM_r[mlowp*C.mn2d_gap],'k','all',ax=ax[2,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_gap],cRM_gap.qr,'teal','Stellar',ax=ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_gap],cRL_gap.qr,'teal','all',ax=ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],cUM_gap.qr,'teal','Stellar',ax=ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],cUL_gap.qr,'teal','all',ax=ax[2,1],style='--')



ax[0,0].set_xlim([0.012,1.01])
ax[0,0].set_xscale('log')

ax[2,0].set_xlabel('$R/R_{200}$')
ax[2,1].set_xlabel('$R/R_{200}$')

ax[0,1].set_yticklabels([])
ax[1,1].set_yticklabels([])
ax[2,1].set_yticklabels([])


ax[0,0].set_ylabel('$T$')
ax[1,0].set_ylabel('$S$')
ax[2,0].set_ylabel('$q$')


ax[0,0].set_ylim([0.2,1.])
ax[0,1].set_ylim([0.2,1.])

ax[1,0].set_ylim([0.4,0.85])
ax[1,1].set_ylim([0.4,0.85])

ax[2,0].set_ylim([0.5,0.95])
ax[2,1].set_ylim([0.5,0.95])

ax[1,0].set_yticks([0.4,0.5,0.6,0.7,0.8])
ax[1,1].set_yticks([0.4,0.5,0.6,0.7,0.8])

ax[2,0].set_yticks([0.5,0.6,0.7,0.8,0.9])
ax[2,1].set_yticks([0.5,0.6,0.7,0.8,0.9])


ax[0,0].text(0.07,0.9,'Relaxed')
ax[0,1].text(0.07,0.9,'Non-relaxed')

plt.savefig(plotspath+'stars_gap.pdf',bbox_inches='tight')

