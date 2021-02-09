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



plotspath = '../final_plots/'
mice = pd.read_pickle('../catalog/data_eli_01.pkl')

C = Clusters()

DM30   = DarkMatter(30,False)
DM50   = DarkMatter(50,False)
DM100  = DarkMatter(100,False)
DM1000 = DarkMatter(1000,False)
DM500  = DarkMatter(500,False)
DM200  = DarkMatter(200,False)


SDM_r = np.vstack((DM30.S,DM50.S,DM100.S,DM1000.S,DM500.S,DM200.S)).T
TDM_r = np.vstack((DM30.T,DM50.T,DM100.T,DM1000.T,DM500.T,DM200.T)).T
qDM_r = np.vstack((DM30.q,DM50.q,DM100.q,DM1000.q,DM500.q,DM200.q)).T


f, ax = plt.subplots(3,3, figsize=(14,14), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

# mmas = C.lM30 > 11.5
# mmasp = C.lM30p > 11.5

# ax[1,2].plot(200,200,'k-',label = '$\log(M_{< 30kpc}) > 11.5$')
# ax[1,2].plot(200,200,'k--',label = '$\log(M_{< 30kpc}) < 11.5$')


mmas = C.lM200 > 14.6
mmasp = C.lM200p > 14.6

ax[1,2].plot(200,200,'k-',label = '$\log(M_{200}) > 14.6$')
ax[1,2].plot(200,200,'k--',label = '$\log(M_{200}) < 14.6$')


mlow = ~mmas
mlowp = ~mmasp


ax[1,2].legend(loc=1,frameon=False)

cRL_gap = CorrelR(Stars,DarkMatter,m = mlow*C.mo_gap, m2d = mlowp*C.mo2d_gap)
cRM_gap = CorrelR(Stars,DarkMatter,m = mmas*C.mo_gap, m2d = mmasp*C.mo2d_gap)
cUL_gap = CorrelR(Stars,DarkMatter,m = mlow*C.mn_gap, m2d = mlowp*C.mn2d_gap)
cUM_gap = CorrelR(Stars,DarkMatter,m = mmas*C.mn_gap, m2d = mmasp*C.mn2d_gap)
cRL_off = CorrelR(Stars,DarkMatter,m = mlow*C.mo_off, m2d = mlowp*C.mo2d_off)
cRM_off = CorrelR(Stars,DarkMatter,m = mmas*C.mo_off, m2d = mmasp*C.mo2d_off)
cUL_off = CorrelR(Stars,DarkMatter,m = mlow*C.mn_off, m2d = mlowp*C.mn2d_off)
cUM_off = CorrelR(Stars,DarkMatter,m = mmas*C.mn_off, m2d = mmasp*C.mn2d_off)
cRL_dv  = CorrelR(Stars,DarkMatter,m = mlow*C.mo_dv, m2d = mlowp*C.mo2d_dv)
cRM_dv  = CorrelR(Stars,DarkMatter,m = mmas*C.mo_dv, m2d = mmasp*C.mo2d_dv)
cUL_dv  = CorrelR(Stars,DarkMatter,m = mlow*C.mn_dv, m2d = mlowp*C.mn2d_dv)
cUM_dv  = CorrelR(Stars,DarkMatter,m = mmas*C.mn_dv, m2d = mmasp*C.mn2d_dv)


mo = DM200.T < 1/3.
mt = (DM200.T > 1/3.)*(DM200.T < 2/3.)
mp = DM200.T > 2/3.




limites = [0.02,1.,0.88,1.12]





plotR_ind(C.Rs[mlow*C.mo_gap], cRL_gap.Tr/TDM_r[mlow*C.mo_gap],'C3','all',ax = ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mo_gap], cRM_gap.Tr/TDM_r[mmas*C.mo_gap],'C3','all',ax = ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap], cUL_gap.Tr/TDM_r[mlow*C.mn_gap],'C0','all',ax = ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap], cUM_gap.Tr/TDM_r[mmas*C.mn_gap],'C0','all',ax = ax[0,0],style='-')                                                                     
plotR_ind(C.Rs[mlow*C.mo_off], cRL_off.Tr/TDM_r[mlow*C.mo_off],'C3','all',ax = ax[0,1],style='--')
plotR_ind(C.Rs[mmas*C.mo_off], cRM_off.Tr/TDM_r[mmas*C.mo_off],'C3','all',ax = ax[0,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_off], cUL_off.Tr/TDM_r[mlow*C.mn_off],'C0','all',ax = ax[0,1],style='--')
plotR_ind(C.Rs[mmas*C.mn_off], cUM_off.Tr/TDM_r[mmas*C.mn_off],'C0','all',ax = ax[0,1],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv] , cRL_dv.Tr/TDM_r[mlow*C.mo_dv],'C3','all',ax = ax[0,2],style='--')
plotR_ind(C.Rs[mmas*C.mo_dv] , cRM_dv.Tr/TDM_r[mmas*C.mo_dv],'C3','all',ax = ax[0,2],style='-')
plotR_ind(C.Rs[mlow*C.mn_dv] , cUL_dv.Tr/TDM_r[mlow*C.mn_dv],'C0','all',ax = ax[0,2],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv] , cUM_dv.Tr/TDM_r[mmas*C.mn_dv],'C0','all',ax = ax[0,2],style='-')

plotR_ind(C.Rs[mlow*C.mo_gap],cRL_gap.Sr/SDM_r[mlow*C.mo_gap],'C3','all',ax = ax[1,0],style='--')
plotR_ind(C.Rs[mmas*C.mo_gap],cRM_gap.Sr/SDM_r[mmas*C.mo_gap],'C3','all',ax = ax[1,0],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap],cUL_gap.Sr/SDM_r[mlow*C.mn_gap],'C0','all',ax = ax[1,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],cUM_gap.Sr/SDM_r[mmas*C.mn_gap],'C0','all',ax = ax[1,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_off],cRL_off.Sr/SDM_r[mlow*C.mo_off],'C3','all',ax = ax[1,1],style='--')
plotR_ind(C.Rs[mmas*C.mo_off],cRM_off.Sr/SDM_r[mmas*C.mo_off],'C3','all',ax = ax[1,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_off],cUL_off.Sr/SDM_r[mlow*C.mn_off],'C0','all',ax = ax[1,1],style='--')
plotR_ind(C.Rs[mmas*C.mn_off],cUM_off.Sr/SDM_r[mmas*C.mn_off],'C0','all',ax = ax[1,1],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv] ,cRL_dv.Sr/SDM_r[mlow*C.mo_dv] ,'C3','all',ax = ax[1,2],style='--')
plotR_ind(C.Rs[mmas*C.mo_dv] ,cRM_dv.Sr/SDM_r[mmas*C.mo_dv] ,'C3','all',ax = ax[1,2],style='-')
plotR_ind(C.Rs[mlow*C.mn_dv] ,cUL_dv.Sr/SDM_r[mlow*C.mn_dv] ,'C0','all',ax = ax[1,2],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv] ,cUM_dv.Sr/SDM_r[mmas*C.mn_dv] ,'C0','all',ax = ax[1,2],style='-')

plotR_ind(C.Rsp[mlowp*C.mo2d_gap],cRL_gap.qr/qDM_r[mlowp*C.mo2d_gap],'C3','all',ax = ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_gap],cRM_gap.qr/qDM_r[mmasp*C.mo2d_gap],'C3','all',ax = ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],cUL_gap.qr/qDM_r[mlowp*C.mn2d_gap],'C0','all',ax = ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],cUM_gap.qr/qDM_r[mmasp*C.mn2d_gap],'C0','all',ax = ax[2,0],style='-')                                                                     
plotR_ind(C.Rsp[mlowp*C.mo2d_off],cRL_off.qr/qDM_r[mlowp*C.mo2d_off],'C3','all',ax = ax[2,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_off],cRM_off.qr/qDM_r[mmasp*C.mo2d_off],'C3','all',ax = ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_off],cUL_off.qr/qDM_r[mlowp*C.mn2d_off],'C0','all',ax = ax[2,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_off],cUM_off.qr/qDM_r[mmasp*C.mn2d_off],'C0','all',ax = ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_dv] ,cRL_dv.qr/qDM_r[mlowp*C.mo2d_dv],'C3','all',ax = ax[2,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_dv] ,cRM_dv.qr/qDM_r[mmasp*C.mo2d_dv],'C3','all',ax = ax[2,2],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_dv] ,cUL_dv.qr/qDM_r[mlowp*C.mn2d_dv],'C0','all',ax = ax[2,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_dv] ,cUM_dv.qr/qDM_r[mmasp*C.mn2d_dv],'C0','all',ax = ax[2,2],style='-')

ax[0,0].set_xlim([0.012,1.01])
ax[0,0].set_xscale('log')

ax[2,0].set_xlabel('$R/R_{200}$')
ax[2,1].set_xlabel('$R/R_{200}$')
ax[2,2].set_xlabel('$R/R_{200}$')

ax[0,1].set_yticklabels([])
ax[0,2].set_yticklabels([])
ax[1,1].set_yticklabels([])
ax[1,2].set_yticklabels([])
ax[2,1].set_yticklabels([])
ax[2,2].set_yticklabels([])

ax[0,0].set_ylabel('$T_{\star}/T_{DM}$')
ax[1,0].set_ylabel('$S_{\star}/S_{DM}$')
ax[2,0].set_ylabel('$q_{\star}/q_{DM}$')

ax[0,0].text(0.07,1.8,'$M_{sat}/M_{BCG}$')
ax[0,1].text(0.07,1.8,'$D_{offset}$')
ax[0,2].text(0.07,1.8,'$\Delta V$')


ax[0,0].set_ylim([1.0,2.])
ax[0,1].set_ylim([1.0,2.])
ax[0,2].set_ylim([1.0,2.])

ax[1,0].set_ylim([0.65,0.93])
ax[1,1].set_ylim([0.65,0.93])
ax[1,2].set_ylim([0.65,0.93])

ax[2,0].set_ylim([0.65,0.93])
ax[2,1].set_ylim([0.65,0.93])
ax[2,2].set_ylim([0.65,0.93])



plt.savefig(plotspath+'StarsH.pdf',bbox_inches='tight')

f, ax = plt.subplots(4,3, figsize=(14,14), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

# mmas = C.lM30 > 11.5
# mmasp = C.lM30p > 11.5

# ax[1,2].plot(200,200,'k-',label = '$\log(M_{< 30kpc}) > 11.5$')
# ax[1,2].plot(200,200,'k--',label = '$\log(M_{< 30kpc}) < 11.5$')


mmas = C.lM200 > 14.6
mmasp = C.lM200p > 14.6

ax[1,2].plot(200,200,'k-',label = '$\log(M_{200}) > 14.6$')
ax[1,2].plot(200,200,'k--',label = '$\log(M_{200}) < 14.6$')


mlow = ~mmas
mlowp = ~mmasp


ax[1,2].legend(loc=2,frameon=False)

plotR_ind(C.Rs[mlow*C.mo_gap],cRL_gap.t3D,'C3','all',ax = ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mo_gap],cRM_gap.t3D,'C3','all',ax = ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap],cUL_gap.t3D,'C0','all',ax = ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],cUM_gap.t3D,'C0','all',ax = ax[0,0],style='-')                                                                     
plotR_ind(C.Rs[mlow*C.mo_off],cRL_off.t3D,'C3','all',ax = ax[0,1],style='--')
plotR_ind(C.Rs[mmas*C.mo_off],cRM_off.t3D,'C3','all',ax = ax[0,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_off],cUL_off.t3D,'C0','all',ax = ax[0,1],style='--')
plotR_ind(C.Rs[mmas*C.mn_off],cUM_off.t3D,'C0','all',ax = ax[0,1],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv] ,cRL_dv.t3D,'C3','all',ax = ax[0,2],style='--')
plotR_ind(C.Rs[mmas*C.mo_dv] ,cRM_dv.t3D,'C3','all',ax = ax[0,2],style='-')
plotR_ind(C.Rs[mlow*C.mn_dv] ,cUL_dv.t3D,'C0','all',ax = ax[0,2],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv] ,cUM_dv.t3D,'C0','all',ax = ax[0,2],style='-')

plotR_ind(C.Rsp[mlowp*C.mo2d_gap],cRL_gap.t2D, 'C3','all',ax = ax[1,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_gap],cRM_gap.t2D, 'C3','all',ax = ax[1,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],cUL_gap.t2D, 'C0','all',ax = ax[1,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],cUM_gap.t2D, 'C0','all',ax = ax[1,0],style='-')                                                                     
plotR_ind(C.Rsp[mlowp*C.mo2d_off],cRL_off.t2D, 'C3','all',ax = ax[1,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_off],cRM_off.t2D, 'C3','all',ax = ax[1,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_off],cUL_off.t2D, 'C0','all',ax = ax[1,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_off],cUM_off.t2D, 'C0','all',ax = ax[1,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_dv] ,cRL_dv.t2D, 'C3','all',ax = ax[1,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_dv] ,cRM_dv.t2D, 'C3','all',ax = ax[1,2],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_dv] ,cUL_dv.t2D, 'C0','all',ax = ax[1,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_dv] ,cUM_dv.t2D, 'C0','all',ax = ax[1,2],style='-')

plotR_ind(C.Rs[mlow*C.mo_gap],cRL_gap.t3D_dm200,'C3','all',ax = ax[2,0],style='--')
plotR_ind(C.Rs[mmas*C.mo_gap],cRM_gap.t3D_dm200,'C3','all',ax = ax[2,0],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap],cUL_gap.t3D_dm200,'C0','all',ax = ax[2,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],cUM_gap.t3D_dm200,'C0','all',ax = ax[2,0],style='-')                                                                     
plotR_ind(C.Rs[mlow*C.mo_off],cRL_off.t3D_dm200,'C3','all',ax = ax[2,1],style='--')
plotR_ind(C.Rs[mmas*C.mo_off],cRM_off.t3D_dm200,'C3','all',ax = ax[2,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_off],cUL_off.t3D_dm200,'C0','all',ax = ax[2,1],style='--')
plotR_ind(C.Rs[mmas*C.mn_off],cUM_off.t3D_dm200,'C0','all',ax = ax[2,1],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv] ,cRL_dv.t3D_dm200,'C3','all',ax = ax[2,2],style='--')
plotR_ind(C.Rs[mmas*C.mo_dv] ,cRM_dv.t3D_dm200,'C3','all',ax = ax[2,2],style='-')
plotR_ind(C.Rs[mlow*C.mn_dv] ,cUL_dv.t3D_dm200,'C0','all',ax = ax[2,2],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv] ,cUM_dv.t3D_dm200,'C0','all',ax = ax[2,2],style='-')

plotR_ind(C.Rsp[mlowp*C.mo2d_gap],cRL_gap.t2D_dm200, 'C3','all',ax = ax[3,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_gap],cRM_gap.t2D_dm200, 'C3','all',ax = ax[3,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],cUL_gap.t2D_dm200, 'C0','all',ax = ax[3,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],cUM_gap.t2D_dm200, 'C0','all',ax = ax[3,0],style='-')                                                                     
plotR_ind(C.Rsp[mlowp*C.mo2d_off],cRL_off.t2D_dm200, 'C3','all',ax = ax[3,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_off],cRM_off.t2D_dm200, 'C3','all',ax = ax[3,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_off],cUL_off.t2D_dm200, 'C0','all',ax = ax[3,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_off],cUM_off.t2D_dm200, 'C0','all',ax = ax[3,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_dv] ,cRL_dv.t2D_dm200, 'C3','all',ax = ax[3,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_dv] ,cRM_dv.t2D_dm200, 'C3','all',ax = ax[3,2],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_dv] ,cUL_dv.t2D_dm200, 'C0','all',ax = ax[3,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_dv] ,cUM_dv.t2D_dm200, 'C0','all',ax = ax[3,2],style='-')


ax[0,0].set_xlim([0.012,1.01])
ax[0,0].set_xscale('log')

ax[3,0].set_xlabel('$R/R_{200}$')
ax[3,1].set_xlabel('$R/R_{200}$')
ax[3,2].set_xlabel('$R/R_{200}$')

ax[0,1].set_yticklabels([])
ax[0,2].set_yticklabels([])
ax[1,1].set_yticklabels([])
ax[1,2].set_yticklabels([])
ax[2,1].set_yticklabels([])
ax[2,2].set_yticklabels([])
ax[3,1].set_yticklabels([])
ax[3,2].set_yticklabels([])


# '''
# DM


ax[0,0].set_ylabel(r'$\theta$')
ax[1,0].set_ylabel(r'$\theta^*$')
ax[2,0].set_ylabel(r'$\theta_{R200}$')
ax[3,0].set_ylabel(r'$\theta^*_{R200}$')

ax[0,0].text(0.07,25,'$M_{sat}/M_{BCG}$')
ax[0,1].text(0.07,25,'$D_{offset}$')
ax[0,2].text(0.07,25,'$\Delta V$')

ax[0,0].set_ylim([0,30])
ax[0,1].set_ylim([0,30])
ax[0,2].set_ylim([0,30])

ax[1,0].set_ylim([0,30])
ax[1,1].set_ylim([0,30])
ax[1,2].set_ylim([0,30])

ax[2,0].set_ylim([0,65])
ax[2,1].set_ylim([0,65])
ax[2,2].set_ylim([0,65])

ax[3,0].set_ylim([0,65])
ax[3,1].set_ylim([0,65])
ax[3,2].set_ylim([0,65])


plt.savefig(plotspath+'theta_starsH.pdf',bbox_inches='tight')

f, ax = plt.subplots(3,2, figsize=(12,14), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)


plotR_ind(C.Rs[mmas*C.mo_gap],TDM_r[mmas*C.mo_gap],'k','$\log(M_{200}) > 14.6$',ax=ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_gap],TDM_r[mlow*C.mo_gap],'k','$\log(M_{200}) < 14.6$',ax=ax[0,0],style='--')
ax[0,0].legend(loc=3,frameon=False)
plotR_ind(C.Rs[mmas*C.mo_gap],cRM_gap.Tr,'C8','Stellar',ax=ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_gap],cRL_gap.Tr,'C8','all',ax=ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],TDM_r[mmas*C.mn_gap],'k','Dark matter',ax=ax[0,1],style='-')
plotR_ind(C.Rs[mmas*C.mn_gap],cUM_gap.Tr,'C8','Stellar',ax=ax[0,1],style='-')
ax[0,1].legend(loc=3,frameon=False)
plotR_ind(C.Rs[mlow*C.mn_gap],TDM_r[mlow*C.mn_gap],'k','all',ax=ax[0,1],style='--')
plotR_ind(C.Rs[mlow*C.mn_gap],cUL_gap.Tr,'C8','all',ax=ax[0,1],style='--')

plotR_ind(C.Rs[mmas*C.mo_gap],SDM_r[mmas*C.mo_gap],'k','Dark matter',ax=ax[1,0],style='-')
plotR_ind(C.Rs[mmas*C.mo_gap],cRM_gap.Sr,'C8','Stellar',ax=ax[1,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_gap],SDM_r[mlow*C.mo_gap],'k','all',ax=ax[1,0],style='--')
plotR_ind(C.Rs[mlow*C.mo_gap],cRL_gap.Sr,'C8','all',ax=ax[1,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],SDM_r[mmas*C.mn_gap],'k','Dark matter',ax=ax[1,1],style='-')
plotR_ind(C.Rs[mmas*C.mn_gap],cUM_gap.Sr,'C8','Stellar',ax=ax[1,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap],SDM_r[mlow*C.mn_gap],'k','all',ax=ax[1,1],style='--')
plotR_ind(C.Rs[mlow*C.mn_gap],cUL_gap.Sr,'C8','all',ax=ax[1,1],style='--')

plotR_ind(C.Rsp[mmasp*C.mo2d_gap],qDM_r[mmasp*C.mo2d_gap],'k','Dark matter',ax=ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_gap],qDM_r[mlowp*C.mo2d_gap],'k','all',ax=ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],qDM_r[mmasp*C.mn2d_gap],'k','Dark matter',ax=ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],qDM_r[mlowp*C.mn2d_gap],'k','all',ax=ax[2,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_gap],cRM_gap.qr,'C8','Stellar',ax=ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_gap],cRL_gap.qr,'C8','all',ax=ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],cUM_gap.qr,'C8','Stellar',ax=ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],cUL_gap.qr,'C8','all',ax=ax[2,1],style='--')



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


ax[0,0].text(0.1,0.9,'Relaxed')
ax[0,1].text(0.1,0.9,'Non-relaxed')

plt.savefig(plotspath+'Hstars_gap.pdf',bbox_inches='tight')

f, ax = plt.subplots(3,2, figsize=(12,14), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)


plotR_ind(C.Rs[mmas*C.mo_dv],TDM_r[mmas*C.mo_dv],'k','Dark matter',ax=ax[0,0],style='-')
plotR_ind(C.Rs[mmas*C.mo_dv],cRM_dv.Tr,'C8','Stellar',ax=ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv],TDM_r[mlow*C.mo_dv],'k','all',ax=ax[0,0],style='--')
plotR_ind(C.Rs[mlow*C.mo_dv],cRL_dv.Tr,'C8','all',ax=ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv],TDM_r[mmas*C.mn_dv],'k','Dark matter',ax=ax[0,1],style='-')
plotR_ind(C.Rs[mmas*C.mn_dv],cUM_dv.Tr,'C8','Stellar',ax=ax[0,1],style='-')
ax[0,1].legend(loc=3,frameon=False)
plotR_ind(C.Rs[mlow*C.mn_dv],TDM_r[mlow*C.mn_dv],'k','all',ax=ax[0,1],style='--')
plotR_ind(C.Rs[mlow*C.mn_dv],cUL_dv.Tr,'C8','all',ax=ax[0,1],style='--')

plotR_ind(C.Rs[mmas*C.mo_dv],SDM_r[mmas*C.mo_dv],'k','Dark matter',ax=ax[1,0],style='-')
plotR_ind(C.Rs[mmas*C.mo_dv],cRM_dv.Sr,'C8','Stellar',ax=ax[1,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv],SDM_r[mlow*C.mo_dv],'k','all',ax=ax[1,0],style='--')
plotR_ind(C.Rs[mlow*C.mo_dv],cRL_dv.Sr,'C8','all',ax=ax[1,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv],SDM_r[mmas*C.mn_dv],'k','Dark matter',ax=ax[1,1],style='-')
plotR_ind(C.Rs[mmas*C.mn_dv],cUM_dv.Sr,'C8','Stellar',ax=ax[1,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_dv],SDM_r[mlow*C.mn_dv],'k','all',ax=ax[1,1],style='--')
plotR_ind(C.Rs[mlow*C.mn_dv],cUL_dv.Sr,'C8','all',ax=ax[1,1],style='--')

plotR_ind(C.Rsp[mmasp*C.mo2d_dv],qDM_r[mmasp*C.mo2d_dv],'k','Dark matter',ax=ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_dv],qDM_r[mlowp*C.mo2d_dv],'k','all',ax=ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_dv],qDM_r[mmasp*C.mn2d_dv],'k','Dark matter',ax=ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_dv],qDM_r[mlowp*C.mn2d_dv],'k','all',ax=ax[2,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_dv],cRM_dv.qr,'C8','Stellar',ax=ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_dv],cRL_dv.qr,'C8','all',ax=ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_dv],cUM_dv.qr,'C8','Stellar',ax=ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_dv],cUL_dv.qr,'C8','all',ax=ax[2,1],style='--')



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


ax[0,0].set_ylim([0.2,0.9])
ax[0,1].set_ylim([0.2,0.9])

ax[1,0].set_ylim([0.4,0.85])
ax[1,1].set_ylim([0.4,0.85])

ax[2,0].set_ylim([0.5,0.95])
ax[2,1].set_ylim([0.5,0.95])

ax[1,0].set_yticks([0.4,0.5,0.6,0.7,0.8])
ax[1,1].set_yticks([0.4,0.5,0.6,0.7,0.8])

ax[2,0].set_yticks([0.5,0.6,0.7,0.8,0.9])
ax[2,1].set_yticks([0.5,0.6,0.7,0.8,0.9])


ax[1,0].text(0.3,0.8,'Relaxed')
ax[1,1].text(0.3,0.8,'Non-relaxed')

plt.savefig(plotspath+'Hstars_dv.pdf',bbox_inches='tight')

f, ax = plt.subplots(3,2, figsize=(12,14), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)


plotR_ind(C.Rs[mmas*C.mo_off],TDM_r[mmas*C.mo_off],'k','Dark matter',ax=ax[0,0],style='-')
plotR_ind(C.Rs[mmas*C.mo_off],cRM_off.Tr,'C8','Stellar',ax=ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_off],TDM_r[mlow*C.mo_off],'k','all',ax=ax[0,0],style='--')
plotR_ind(C.Rs[mlow*C.mo_off],cRL_off.Tr,'C8','all',ax=ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_off],TDM_r[mmas*C.mn_off],'k','Dark matter',ax=ax[0,1],style='-')
plotR_ind(C.Rs[mmas*C.mn_off],cUM_off.Tr,'C8','Stellar',ax=ax[0,1],style='-')
ax[0,1].legend(loc=3,frameon=False)
plotR_ind(C.Rs[mlow*C.mn_off],TDM_r[mlow*C.mn_off],'k','all',ax=ax[0,1],style='--')
plotR_ind(C.Rs[mlow*C.mn_off],cUL_off.Tr,'C8','all',ax=ax[0,1],style='--')

plotR_ind(C.Rs[mmas*C.mo_off],SDM_r[mmas*C.mo_off],'k','Dark matter',ax=ax[1,0],style='-')
plotR_ind(C.Rs[mmas*C.mo_off],cRM_off.Sr,'C8','Stellar',ax=ax[1,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_off],SDM_r[mlow*C.mo_off],'k','all',ax=ax[1,0],style='--')
plotR_ind(C.Rs[mlow*C.mo_off],cRL_off.Sr,'C8','all',ax=ax[1,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_off],SDM_r[mmas*C.mn_off],'k','Dark matter',ax=ax[1,1],style='-')
plotR_ind(C.Rs[mmas*C.mn_off],cUM_off.Sr,'C8','Stellar',ax=ax[1,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_off],SDM_r[mlow*C.mn_off],'k','all',ax=ax[1,1],style='--')
plotR_ind(C.Rs[mlow*C.mn_off],cUL_off.Sr,'C8','all',ax=ax[1,1],style='--')

plotR_ind(C.Rsp[mmasp*C.mo2d_off],qDM_r[mmasp*C.mo2d_off],'k','Dark matter',ax=ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_off],qDM_r[mlowp*C.mo2d_off],'k','all',ax=ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_off],qDM_r[mmasp*C.mn2d_off],'k','Dark matter',ax=ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_off],qDM_r[mlowp*C.mn2d_off],'k','all',ax=ax[2,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_off],cRM_off.qr,'C8','Stellar',ax=ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_off],cRL_off.qr,'C8','all',ax=ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_off],cUM_off.qr,'C8','Stellar',ax=ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_off],cUL_off.qr,'C8','all',ax=ax[2,1],style='--')



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


ax[0,0].set_ylim([0.2,0.9])
ax[0,1].set_ylim([0.2,0.9])

ax[1,0].set_ylim([0.4,0.85])
ax[1,1].set_ylim([0.4,0.85])

ax[2,0].set_ylim([0.5,0.95])
ax[2,1].set_ylim([0.5,0.95])


ax[1,0].set_yticks([0.4,0.5,0.6,0.7,0.8])
ax[1,1].set_yticks([0.4,0.5,0.6,0.7,0.8])

ax[2,0].set_yticks([0.5,0.6,0.7,0.8,0.9])
ax[2,1].set_yticks([0.5,0.6,0.7,0.8,0.9])



ax[1,0].text(0.3,0.8,'Relaxed')
ax[1,1].text(0.3,0.8,'Non-relaxed')

plt.savefig(plotspath+'Hstars_off.pdf',bbox_inches='tight')
