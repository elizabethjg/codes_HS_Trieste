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

DM30   = DarkMatter(30)
DM50   = DarkMatter(50)
DM100  = DarkMatter(100)
DM1000 = DarkMatter(1000)
DM500  = DarkMatter(500)
DM200  = DarkMatter(200)


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


ax[1,2].legend(loc=3,frameon=False)

cDMRL_gap = CorrelR(DarkMatter,DarkMatter,m = mlow*C.mo_gap, m2d = mlowp*C.mo2d_gap)
cDMRM_gap = CorrelR(DarkMatter,DarkMatter,m = mmas*C.mo_gap, m2d = mmasp*C.mo2d_gap)
cDMUL_gap = CorrelR(DarkMatter,DarkMatter,m = mlow*C.mn_gap, m2d = mlowp*C.mn2d_gap)
cDMUM_gap = CorrelR(DarkMatter,DarkMatter,m = mmas*C.mn_gap, m2d = mmasp*C.mn2d_gap)

cDMRL_off = CorrelR(DarkMatter,DarkMatter,m = mlow*C.mo_off, m2d = mlowp*C.mo2d_off)
cDMRM_off = CorrelR(DarkMatter,DarkMatter,m = mmas*C.mo_off, m2d = mmasp*C.mo2d_off)
cDMUL_off = CorrelR(DarkMatter,DarkMatter,m = mlow*C.mn_off, m2d = mlowp*C.mn2d_off)
cDMUM_off = CorrelR(DarkMatter,DarkMatter,m = mmas*C.mn_off, m2d = mmasp*C.mn2d_off)

cDMRL_dv  = CorrelR(DarkMatter,DarkMatter,m = mlow*C.mo_dv, m2d = mlowp*C.mo2d_dv)
cDMRM_dv  = CorrelR(DarkMatter,DarkMatter,m = mmas*C.mo_dv, m2d = mmasp*C.mo2d_dv)
cDMUL_dv  = CorrelR(DarkMatter,DarkMatter,m = mlow*C.mn_dv, m2d = mlowp*C.mn2d_dv)
cDMUM_dv  = CorrelR(DarkMatter,DarkMatter,m = mmas*C.mn_dv, m2d = mmasp*C.mn2d_dv)


mo = DM200.T < 1/3.
mt = (DM200.T > 1/3.)*(DM200.T < 2/3.)
mp = DM200.T > 2/3.




limites = [0.02,1.,0.88,1.12]





plotR_ind(C.Rs[mlow*C.mo_gap],TDM_r[mlow*C.mo_gap],'C3','all',ax = ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mo_gap],TDM_r[mmas*C.mo_gap],'C3','all',ax = ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap],TDM_r[mlow*C.mn_gap],'C0','all',ax = ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],TDM_r[mmas*C.mn_gap],'C0','all',ax = ax[0,0],style='-')                                                                     
plotR_ind(C.Rs[mlow*C.mo_off],TDM_r[mlow*C.mo_off],'C3','all',ax = ax[0,1],style='--')
plotR_ind(C.Rs[mmas*C.mo_off],TDM_r[mmas*C.mo_off],'C3','all',ax = ax[0,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_off],TDM_r[mlow*C.mn_off],'C0','all',ax = ax[0,1],style='--')
plotR_ind(C.Rs[mmas*C.mn_off],TDM_r[mmas*C.mn_off],'C0','all',ax = ax[0,1],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv] ,TDM_r[mlow*C.mo_dv] ,'C3','all',ax = ax[0,2],style='--')
plotR_ind(C.Rs[mmas*C.mo_dv] ,TDM_r[mmas*C.mo_dv] ,'C3','all',ax = ax[0,2],style='-')
plotR_ind(C.Rs[mlow*C.mn_dv] ,TDM_r[mlow*C.mn_dv] ,'C0','all',ax = ax[0,2],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv] ,TDM_r[mmas*C.mn_dv] ,'C0','all',ax = ax[0,2],style='-')

plotR_ind(C.Rs[mlow*C.mo_gap],SDM_r[mlow*C.mo_gap],'C3','all',ax = ax[1,0],style='--')
plotR_ind(C.Rs[mmas*C.mo_gap],SDM_r[mmas*C.mo_gap],'C3','all',ax = ax[1,0],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap],SDM_r[mlow*C.mn_gap],'C0','all',ax = ax[1,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],SDM_r[mmas*C.mn_gap],'C0','all',ax = ax[1,0],style='-')
plotR_ind(C.Rs[mlow*C.mo_off],SDM_r[mlow*C.mo_off],'C3','all',ax = ax[1,1],style='--')
plotR_ind(C.Rs[mmas*C.mo_off],SDM_r[mmas*C.mo_off],'C3','all',ax = ax[1,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_off],SDM_r[mlow*C.mn_off],'C0','all',ax = ax[1,1],style='--')
plotR_ind(C.Rs[mmas*C.mn_off],SDM_r[mmas*C.mn_off],'C0','all',ax = ax[1,1],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv] ,SDM_r[mlow*C.mo_dv] ,'C3','all',ax = ax[1,2],style='--')
plotR_ind(C.Rs[mmas*C.mo_dv] ,SDM_r[mmas*C.mo_dv] ,'C3','all',ax = ax[1,2],style='-')
plotR_ind(C.Rs[mlow*C.mn_dv] ,SDM_r[mlow*C.mn_dv] ,'C0','all',ax = ax[1,2],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv] ,SDM_r[mmas*C.mn_dv] ,'C0','all',ax = ax[1,2],style='-')

plotR_ind(C.Rsp[mlowp*C.mo2d_gap],qDM_r[mlowp*C.mo2d_gap],'C3','all',ax = ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_gap],qDM_r[mmasp*C.mo2d_gap],'C3','all',ax = ax[2,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],qDM_r[mlowp*C.mn2d_gap],'C0','all',ax = ax[2,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],qDM_r[mmasp*C.mn2d_gap],'C0','all',ax = ax[2,0],style='-')                                                                     
plotR_ind(C.Rsp[mlowp*C.mo2d_off],qDM_r[mlowp*C.mo2d_off],'C3','all',ax = ax[2,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_off],qDM_r[mmasp*C.mo2d_off],'C3','all',ax = ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_off],qDM_r[mlowp*C.mn2d_off],'C0','all',ax = ax[2,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_off],qDM_r[mmasp*C.mn2d_off],'C0','all',ax = ax[2,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_dv] ,qDM_r[mlowp*C.mo2d_dv] ,'C3','all',ax = ax[2,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_dv] ,qDM_r[mmasp*C.mo2d_dv] ,'C3','all',ax = ax[2,2],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_dv] ,qDM_r[mlowp*C.mn2d_dv] ,'C0','all',ax = ax[2,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_dv] ,qDM_r[mmasp*C.mn2d_dv] ,'C0','all',ax = ax[2,2],style='-')

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

ax[0,0].set_ylabel('$T_{DM}$')
ax[1,0].set_ylabel('$S_{DM}$')
ax[2,0].set_ylabel('$q_{DM}$')

ax[0,0].text(0.07,0.8,'$M_{sat}/M_{BCG}$')
ax[0,1].text(0.07,0.8,'$D_{offset}$')
ax[0,2].text(0.07,0.8,'$\Delta V$')

ax[0,0].set_ylim([0.2,0.9])
ax[0,1].set_ylim([0.2,0.9])
ax[0,2].set_ylim([0.2,0.9])

ax[1,0].set_ylim([0.55,0.85])
ax[1,1].set_ylim([0.55,0.85])
ax[1,2].set_ylim([0.55,0.85])

ax[2,0].set_ylim([0.7,0.95])
ax[2,1].set_ylim([0.7,0.95])
ax[2,2].set_ylim([0.7,0.95])

ax[0,0].set_yticks([0.3,0.5,0.7,0.9])
ax[0,1].set_yticks([0.3,0.5,0.7,0.9])
ax[0,2].set_yticks([0.3,0.5,0.7,0.9])

ax[1,0].set_yticks([0.6,0.7,0.8])
ax[1,1].set_yticks([0.6,0.7,0.8])
ax[1,2].set_yticks([0.6,0.7,0.8])

ax[2,0].set_yticks([0.7,0.8,0.9])
ax[2,1].set_yticks([0.7,0.8,0.9])
ax[2,2].set_yticks([0.7,0.8,0.9])

plt.savefig(plotspath+'DM.pdf',bbox_inches='tight')

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


ax[1,2].legend(loc=1,frameon=False)

plotR_ind(C.Rs[mlow*C.mo_gap],cDMRL_gap.t3D_2,'C3','all',ax = ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mo_gap],cDMRM_gap.t3D_2,'C3','all',ax = ax[0,0],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap],cDMUL_gap.t3D_2,'C0','all',ax = ax[0,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],cDMUM_gap.t3D_2,'C0','all',ax = ax[0,0],style='-')                                                                     
plotR_ind(C.Rs[mlow*C.mo_off],cDMRL_off.t3D_2,'C3','all',ax = ax[0,1],style='--')
plotR_ind(C.Rs[mmas*C.mo_off],cDMRM_off.t3D_2,'C3','all',ax = ax[0,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_off],cDMUL_off.t3D_2,'C0','all',ax = ax[0,1],style='--')
plotR_ind(C.Rs[mmas*C.mn_off],cDMUM_off.t3D_2,'C0','all',ax = ax[0,1],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv] , cDMRL_dv.t3D_2,'C3','all',ax = ax[0,2],style='--')
plotR_ind(C.Rs[mmas*C.mo_dv] , cDMRM_dv.t3D_2,'C3','all',ax = ax[0,2],style='-')
plotR_ind(C.Rs[mlow*C.mn_dv] , cDMUL_dv.t3D_2,'C0','all',ax = ax[0,2],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv] , cDMUM_dv.t3D_2,'C0','all',ax = ax[0,2],style='-')

plotR_ind(C.Rsp[mlowp*C.mo2d_gap],cDMRL_gap.t2D_2, 'C3','all',ax = ax[1,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_gap],cDMRM_gap.t2D_2, 'C3','all',ax = ax[1,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],cDMUL_gap.t2D_2, 'C0','all',ax = ax[1,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],cDMUM_gap.t2D_2, 'C0','all',ax = ax[1,0],style='-')                                                                     
plotR_ind(C.Rsp[mlowp*C.mo2d_off],cDMRL_off.t2D_2, 'C3','all',ax = ax[1,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_off],cDMRM_off.t2D_2, 'C3','all',ax = ax[1,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_off],cDMUL_off.t2D_2, 'C0','all',ax = ax[1,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_off],cDMUM_off.t2D_2, 'C0','all',ax = ax[1,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_dv] , cDMRL_dv.t2D_2, 'C3','all',ax = ax[1,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_dv] , cDMRM_dv.t2D_2, 'C3','all',ax = ax[1,2],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_dv] , cDMUL_dv.t2D_2, 'C0','all',ax = ax[1,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_dv] , cDMUM_dv.t2D_2, 'C0','all',ax = ax[1,2],style='-')

plotR_ind(C.Rs[mlow*C.mo_gap],cDMRL_gap.t3D_dm200,'C3','all',ax = ax[2,0],style='--')
plotR_ind(C.Rs[mmas*C.mo_gap],cDMRM_gap.t3D_dm200,'C3','all',ax = ax[2,0],style='-')
plotR_ind(C.Rs[mlow*C.mn_gap],cDMUL_gap.t3D_dm200,'C0','all',ax = ax[2,0],style='--')
plotR_ind(C.Rs[mmas*C.mn_gap],cDMUM_gap.t3D_dm200,'C0','all',ax = ax[2,0],style='-')                                                                     
plotR_ind(C.Rs[mlow*C.mo_off],cDMRL_off.t3D_dm200,'C3','all',ax = ax[2,1],style='--')
plotR_ind(C.Rs[mmas*C.mo_off],cDMRM_off.t3D_dm200,'C3','all',ax = ax[2,1],style='-')
plotR_ind(C.Rs[mlow*C.mn_off],cDMUL_off.t3D_dm200,'C0','all',ax = ax[2,1],style='--')
plotR_ind(C.Rs[mmas*C.mn_off],cDMUM_off.t3D_dm200,'C0','all',ax = ax[2,1],style='-')
plotR_ind(C.Rs[mlow*C.mo_dv] , cDMRL_dv.t3D_dm200,'C3','all',ax = ax[2,2],style='--')
plotR_ind(C.Rs[mmas*C.mo_dv] , cDMRM_dv.t3D_dm200,'C3','all',ax = ax[2,2],style='-')
plotR_ind(C.Rs[mlow*C.mn_dv] , cDMUL_dv.t3D_dm200,'C0','all',ax = ax[2,2],style='--')
plotR_ind(C.Rs[mmas*C.mn_dv] , cDMUM_dv.t3D_dm200,'C0','all',ax = ax[2,2],style='-')

plotR_ind(C.Rsp[mlowp*C.mo2d_gap],cDMRL_gap.t2D_dm200, 'C3','all',ax = ax[3,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_gap],cDMRM_gap.t2D_dm200, 'C3','all',ax = ax[3,0],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_gap],cDMUL_gap.t2D_dm200, 'C0','all',ax = ax[3,0],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_gap],cDMUM_gap.t2D_dm200, 'C0','all',ax = ax[3,0],style='-')                                                                     
plotR_ind(C.Rsp[mlowp*C.mo2d_off],cDMRL_off.t2D_dm200, 'C3','all',ax = ax[3,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_off],cDMRM_off.t2D_dm200, 'C3','all',ax = ax[3,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_off],cDMUL_off.t2D_dm200, 'C0','all',ax = ax[3,1],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_off],cDMUM_off.t2D_dm200, 'C0','all',ax = ax[3,1],style='-')
plotR_ind(C.Rsp[mlowp*C.mo2d_dv] , cDMRL_dv.t2D_dm200, 'C3','all',ax = ax[3,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mo2d_dv] , cDMRM_dv.t2D_dm200, 'C3','all',ax = ax[3,2],style='-')
plotR_ind(C.Rsp[mlowp*C.mn2d_dv] , cDMUL_dv.t2D_dm200, 'C0','all',ax = ax[3,2],style='--')
plotR_ind(C.Rsp[mmasp*C.mn2d_dv] , cDMUM_dv.t2D_dm200, 'C0','all',ax = ax[3,2],style='-')


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


ax[0,0].set_ylabel(r'$\theta_{R_{i+1}}$')
ax[1,0].set_ylabel(r'$\theta^*_{R_{i+1}}$')
ax[2,0].set_ylabel(r'$\theta_{R200}$')
ax[3,0].set_ylabel(r'$\theta^*_{R200}$')

ax[0,0].text(0.07,50,'$M_{sat}/M_{BCG}$')
ax[0,1].text(0.07,50,'$D_{offset}$')
ax[0,2].text(0.07,50,'$\Delta V$')

ax[0,0].set_ylim([0,65])
ax[0,1].set_ylim([0,65])
ax[0,2].set_ylim([0,65])

ax[1,0].set_ylim([0,65])
ax[1,1].set_ylim([0,65])
ax[1,2].set_ylim([0,65])

ax[2,0].set_ylim([0,65])
ax[2,1].set_ylim([0,65])
ax[2,2].set_ylim([0,65])

ax[3,0].set_ylim([0,65])
ax[3,1].set_ylim([0,65])
ax[3,2].set_ylim([0,65])


plt.savefig(plotspath+'theta_DM.pdf',bbox_inches='tight')


'''

# Stars

ax[0,0].set_ylabel('$T\star$')
ax[1,0].set_ylabel('$S\star$')
ax[2,0].set_ylabel('$q\star$')

ax[0,0].set_ylim([0.5,1.0])
ax[0,1].set_ylim([0.5,1.0])
ax[0,2].set_ylim([0.5,1.0])

ax[1,0].set_ylim([0.4,0.7])
ax[1,1].set_ylim([0.4,0.7])
ax[1,2].set_ylim([0.4,0.7])

ax[2,0].set_ylim([0.5,0.8])
ax[2,1].set_ylim([0.5,0.8])
ax[2,2].set_ylim([0.5,0.8])

ax[0,0].set_yticks([0.5,0.7,0.9])
ax[0,1].set_yticks([0.5,0.7,0.9])
ax[0,2].set_yticks([0.5,0.7,0.9])

ax[1,0].set_yticks([0.4,0.5,0.6])
ax[1,1].set_yticks([0.4,0.5,0.6])
ax[1,2].set_yticks([0.4,0.5,0.6])

ax[2,0].set_yticks([0.5,0.6,0.7])
ax[2,1].set_yticks([0.5,0.6,0.7])
ax[2,2].set_yticks([0.5,0.6,0.7])

ax[0,0].text(0.07,0.95,'$M_{sat}/M_{BCG}$')
ax[0,1].text(0.07,0.95,'$D_{offset}$')
ax[0,2].text(0.07,0.95,'$\Delta V$')

plt.savefig(plotspath+'Stars.pdf',bbox_inches='tight')
'''
