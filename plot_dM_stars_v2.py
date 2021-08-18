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


st1000 = Stars(1000)
st500  = Stars(500)
st200  = Stars(200)
st100  = Stars(100)
st50   = Stars(50)
st30   = Stars(30)

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

##################

colors = ['navy','indigo','brown','C3','C1','y']
radio = [30,50,100,1000,500,200]
labels = ['30kpc','50kpc','0.1$R_{500}$','$R_{1000}$','$R_{500}$','$R_{200}$']

mold = (C.ltime > 0.)*(C.ltime > 4.9)
mold_p = (C.ltimep > 0.)*(C.ltimep > 4.9)
mnew = (C.ltime > 0.)*(C.ltime < 4.9)
mnew_p = (C.ltimep > 0.)*(C.ltimep < 4.9)

mold   = C.mo_gap
mold_p = C.mo2d_gap
mnew   = C.mn_gap
mnew_p = C.mn2d_gap

f, ax = plt.subplots(3,1, figsize=(5,12))
f2, ax2 = plt.subplots(3,1, figsize=(5,12))
f3, ax3 = plt.subplots(3,1, figsize=(5,12))
f4, ax4 = plt.subplots(3,1, figsize=(5,12))
f5, ax5 = plt.subplots(3,1, figsize=(5,12))

for j in range(6):
    
    plot_fig(Stars(radio[j]).T,Stars(radio[j]).T/DarkMatter(radio[j]).T,5,color=colors[j],ax=ax[0])
    plot_fig(Stars(radio[j]).S,Stars(radio[j]).S/DarkMatter(radio[j]).S,5,color=colors[j],ax=ax[1],label=labels[j])
    plot_fig(Stars(radio[j]).q,Stars(radio[j]).q/DarkMatter(radio[j]).q,5,color=colors[j],ax=ax[2],label=labels[j])

    plot_fig(DarkMatter(radio[j]).T,Stars(radio[j]).T/DarkMatter(radio[j]).T,5,color=colors[j],ax=ax3[0])
    plot_fig(DarkMatter(radio[j]).S,Stars(radio[j]).S/DarkMatter(radio[j]).S,5,color=colors[j],ax=ax3[1],label=labels[j])
    plot_fig(DarkMatter(radio[j]).q,Stars(radio[j]).q/DarkMatter(radio[j]).q,5,color=colors[j],ax=ax3[2],label=labels[j])

    plot_fig(Stars(radio[j]).T,DarkMatter(radio[j]).T/Stars(radio[j]).T,5,color=colors[j],ax=ax4[0])
    plot_fig(Stars(radio[j]).S,DarkMatter(radio[j]).S/Stars(radio[j]).S,5,color=colors[j],ax=ax4[1],label=labels[j])
    plot_fig(Stars(radio[j]).q,DarkMatter(radio[j]).q/Stars(radio[j]).q,5,color=colors[j],ax=ax4[2],label=labels[j])

    plot_fig(DarkMatter(radio[j]).T,DarkMatter(radio[j]).T/Stars(radio[j]).T,5,color=colors[j],ax=ax5[0])
    plot_fig(DarkMatter(radio[j]).S,DarkMatter(radio[j]).S/Stars(radio[j]).S,5,color=colors[j],ax=ax5[1],label=labels[j])
    plot_fig(DarkMatter(radio[j]).q,DarkMatter(radio[j]).q/Stars(radio[j]).q,5,color=colors[j],ax=ax5[2],label=labels[j])
    
    if j == 2 or j == 5:
        
        plot_fig(Stars(radio[j]).T,Stars(radio[j]).T/DarkMatter(radio[j]).T,5,color=colors[j],ax=ax2[0])
        plot_fig(Stars(radio[j]).S,Stars(radio[j]).S/DarkMatter(radio[j]).S,5,color=colors[j],ax=ax2[1],label=labels[j])
        plot_fig(Stars(radio[j]).q,Stars(radio[j]).q/DarkMatter(radio[j]).q,5,color=colors[j],ax=ax2[2],label=labels[j])

        plot_fig(Stars(radio[j]).T[mold],(Stars(radio[j]).T/DarkMatter(radio[j]).T)[mold],5,color=colors[j],ax=ax2[0],style='--')
        plot_fig(Stars(radio[j]).S[mold],(Stars(radio[j]).S/DarkMatter(radio[j]).S)[mold],5,color=colors[j],ax=ax2[1],label=labels[j],style='--')
        plot_fig(Stars(radio[j]).q[mold_p],(Stars(radio[j]).q/DarkMatter(radio[j]).q)[mold_p],5,color=colors[j],ax=ax2[2],label=labels[j],style='--')
        plot_fig(Stars(radio[j]).T[mnew],(Stars(radio[j]).T/DarkMatter(radio[j]).T)[mnew],5,color=colors[j],ax=ax2[0],style=':')
        plot_fig(Stars(radio[j]).S[mnew],(Stars(radio[j]).S/DarkMatter(radio[j]).S)[mnew],5,color=colors[j],ax=ax2[1],label=labels[j],style=':')
        plot_fig(Stars(radio[j]).q[mnew_p],(Stars(radio[j]).q/DarkMatter(radio[j]).q)[mnew_p],5,color=colors[j],ax=ax2[2],label=labels[j],style=':')


    
ax[2].legend(loc=2,frameon=False,fontsize = 12,ncol=2)
ax[0].set_ylabel('$T\star / T_{DM}$')
ax[1].set_ylabel('$S\star / S_{DM}$')
ax[2].set_ylabel('$q\star / q_{DM}$')
ax[2].set_xlabel('$q\star$')
ax[1].set_xlabel('$S\star$')
ax[0].set_xlabel('$T\star$')


ax2[0].set_ylabel('$T\star / T_{DM}$')
ax2[1].set_ylabel('$S\star / S_{DM}$')
ax2[2].set_ylabel('$q\star / q_{DM}$')
ax2[2].set_xlabel('$q\star$')
ax2[1].set_xlabel('$S\star$')
ax2[0].set_xlabel('$T\star$')

ax3[0].set_ylabel('$T\star / T_{DM}$')
ax3[1].set_ylabel('$S\star / S_{DM}$')
ax3[2].set_ylabel('$q\star / q_{DM}$')
ax3[2].set_xlabel('$q_{DM}$')
ax3[1].set_xlabel('$S_{DM}$')
ax3[0].set_xlabel('$T_{DM}$')

ax4[0].set_ylabel('$T_{DM} / T\star$')
ax4[1].set_ylabel('$S_{DM} / S\star$')
ax4[2].set_ylabel('$q_{DM} / q\star$')
ax4[2].set_xlabel('$q\star$')
ax4[1].set_xlabel('$S\star$')
ax4[0].set_xlabel('$T\star$')

ax5[0].set_ylabel('$T_{DM} / T\star$')
ax5[1].set_ylabel('$S_{DM} / S\star$')
ax5[2].set_ylabel('$q_{DM} / q\star$')
ax5[2].set_xlabel('$q_{DM}$')
ax5[1].set_xlabel('$S_{DM}$')
ax5[0].set_xlabel('$T_{DM}$')

ax[0].set_xlim([0.15,0.9])
ax[1].set_xlim([0.32,0.75])
ax[2].set_xlim([0.4,0.95])
ax[0].set_ylim([0.0,5.0])
ax[1].set_ylim([0.5,1.0])
ax[2].set_ylim([0.5,1.2])

ax2[0].set_xlim([0.15,0.9])
ax2[1].set_xlim([0.32,0.75])
ax2[2].set_xlim([0.4,0.95])
ax2[0].set_ylim([0.0,5.0])
ax2[1].set_ylim([0.5,1.0])
ax2[2].set_ylim([0.5,1.2])

ax3[0].set_xlim([0.1,0.9])
ax3[1].set_xlim([0.5,0.85])
ax3[2].set_xlim([0.5,1.0])
ax3[0].set_ylim([0.0,5.0])
ax3[1].set_ylim([0.5,1.0])
ax3[2].set_ylim([0.5,1.2])

ax4[0].set_xlim([0.15,0.9])
ax4[1].set_xlim([0.32,0.75])
ax4[2].set_xlim([0.4,0.95])
ax4[0].set_ylim([0.0,1.2])
ax4[1].set_ylim([0.7,2.0])
ax4[2].set_ylim([0.7,2.0])

ax5[0].set_xlim([0.1,0.9])
ax5[1].set_xlim([0.5,0.85])
ax5[2].set_xlim([0.5,1.0])
ax5[0].set_ylim([0.0,1.2])
ax5[1].set_ylim([0.7,2.0])
ax5[2].set_ylim([0.7,2.0])



f.savefig(plotspath+'ratio_stars.pdf',bbox_inches='tight')
f4.savefig(plotspath+'ratio_stars_inv.pdf',bbox_inches='tight')
f2.savefig(plotspath+'ratio_stars_relax_nonrelax.pdf',bbox_inches='tight')
f3.savefig(plotspath+'ratio_stars_dm.pdf',bbox_inches='tight')
f5.savefig(plotspath+'ratio_stars_dm_inv.pdf',bbox_inches='tight')

f, ax = plt.subplots(3,1, figsize=(5,12))
f2, ax2 = plt.subplots(3,1, figsize=(5,12))
# f.subplots_adjust(hspace=0,wspace=0)

for j in range(6):
    plot_fig(Stars(radio[j]).T,DarkMatter(radio[j]).T,5,color=colors[j],ax=ax[0])
    plot_fig(Stars(radio[j]).S,DarkMatter(radio[j]).S,5,color=colors[j],ax=ax[1])
    plot_fig(Stars(radio[j]).q,DarkMatter(radio[j]).q,5,color=colors[j],ax=ax[2],label=labels[j])

    plot_fig(DarkMatter(radio[j]).T,Stars(radio[j]).T,5,color=colors[j],ax=ax2[0])
    plot_fig(DarkMatter(radio[j]).S,Stars(radio[j]).S,5,color=colors[j],ax=ax2[1])
    plot_fig(DarkMatter(radio[j]).q,Stars(radio[j]).q,5,color=colors[j],ax=ax2[2],label=labels[j])
    
ax[2].legend(loc=4,frameon=False,fontsize = 12,ncol=2)
ax[0].set_ylabel('$T_{DM}$')
ax[1].set_ylabel('$S_{DM}$')
ax[2].set_ylabel('$q_{DM}$')
ax[0].set_xlabel('$T\star$')
ax[2].set_xlabel('$q\star$')
ax[1].set_xlabel('$S\star$')

ax2[0].set_xlabel('$T_{DM}$')
ax2[1].set_xlabel('$S_{DM}$')
ax2[2].set_xlabel('$q_{DM}$')
ax2[0].set_ylabel('$T\star$')
ax2[2].set_ylabel('$q\star$')
ax2[1].set_ylabel('$S\star$')

ax[0].plot([0,1],[0,1],'C7--')
ax[1].plot([0,1],[0,1],'C7--')
ax[2].plot([0,1],[0,1],'C7--')

ax2[0].plot([0,1],[0,1],'C7--')
ax2[1].plot([0,1],[0,1],'C7--')
ax2[2].plot([0,1],[0,1],'C7--')

ax[0].set_xlim([0.15,0.9])
ax[1].set_xlim([0.32,0.75])
ax[2].set_xlim([0.4,0.95])
ax[0].set_ylim([0.1,0.9])
ax[1].set_ylim([0.4,0.9])
ax[2].set_ylim([0.5,1.0])

ax2[0].set_ylim([0.15,0.9])
ax2[1].set_ylim([0.32,0.75])
ax2[2].set_ylim([0.4,0.95])
ax2[0].set_xlim([0.1,0.9])
ax2[1].set_xlim([0.4,0.9])
ax2[2].set_xlim([0.5,1.0])

f.savefig(plotspath+'compare_stars.pdf',bbox_inches='tight')
f2.savefig(plotspath+'compare_stars_dm.pdf',bbox_inches='tight')

'''
f, ax = plt.subplots(3,1, figsize=(5,12))

for j in range(6):
    plot_fig(DarkMatter(radio[j]).T,Stars(radio[j]).T/DarkMatter(radio[j]).T,5,color=colors[j],ax=ax[0])
    plot_fig(DarkMatter(radio[j]).S,Stars(radio[j]).S/DarkMatter(radio[j]).S,5,color=colors[j],ax=ax[1],label=labels[j])
    plot_fig(DarkMatter(radio[j]).q,Stars(radio[j]).q/DarkMatter(radio[j]).q,5,color=colors[j],ax=ax[2],label=labels[j])
    
ax[2].legend(loc=2,frameon=False,fontsize = 12,ncol=2)
ax[0].set_ylabel('$T\star / T_{DM}$')
ax[1].set_ylabel('$S\star / S_{DM}$')
ax[2].set_ylabel('$q\star / q_{DM}$')
ax[2].set_xlabel('$q_{DM}$')
ax[1].set_xlabel('$S_{DM}$')
ax[0].set_xlabel('$T_{DM}$')

plt.savefig(plotspath+'ratio_stars_dm.pdf',bbox_inches='tight')

f, ax = plt.subplots(3,1, figsize=(5,12))
# f.subplots_adjust(hspace=0,wspace=0)

for j in range(6):
    plot_fig(DarkMatter(radio[j]).T,DarkMatter(radio[j]).T,5,color=colors[j],ax=ax[0])
    plot_fig(DarkMatter(radio[j]).S,DarkMatter(radio[j]).S,5,color=colors[j],ax=ax[1])
    plot_fig(DarkMatter(radio[j]).q,DarkMatter(radio[j]).q,5,color=colors[j],ax=ax[2],label=labels[j])
    
ax[2].legend(loc=4,frameon=False,fontsize = 12,ncol=2)
ax[0].set_ylabel('$T_{DM}$')
ax[1].set_ylabel('$S_{DM}$')
ax[2].set_ylabel('$q_{DM}$')
ax[0].set_xlabel('$T_{DM}$')
ax[2].set_xlabel('$q_{DM}$')
ax[1].set_xlabel('$S_{DM}$')

plt.savefig(plotspath+'compare_stars_dm.pdf',bbox_inches='tight')
'''
