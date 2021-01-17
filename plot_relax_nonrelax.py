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
matplotlib.rcParams.update({'font.size': 12})


path = '../catalog/nuevosdats/'

gral  = np.loadtxt(path+'gral_nounb_091.dat').T

plotspath = '../final_plots/'

lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)

off   = gral[13]
off2D = np.concatenate((gral[14],gral[15],gral[16]))
DV    = gral[17]
DV2D  = np.concatenate((gral[18],gral[19],gral[20]))
gap  = gral[21]
gap2D = np.concatenate((gral[22],gral[23],gral[24]))

mn_off,mo_off,mn2d_off,mo2d_off = newold('off')
mn_gap,mo_gap,mn2d_gap,mo2d_gap = newold('gap')
mn_dv ,mo_dv, mn2d_dv ,mo2d_dv  = newold('DV')

f, ax = plt.subplots(2,1, figsize=(4,6), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].plot(lM[mo_gap],gap[mo_gap],'C3o',label='relaxed',alpha=0.7)
ax[0].plot(lM[mn_gap],gap[mn_gap],'C0s',label='non-relaxed',alpha=0.7)
ax[0].set_ylabel(r'$M_{sat}/M_{BCG}$')
ax[0].set_xlabel(u'$\log(M_{200})$')
ax[0].legend(loc=2,frameon=False)

ax[1].plot(lMp[mo2d_gap],gap2D[mo2d_gap],'C3o',label='relaxed',alpha=0.7)
ax[1].plot(lMp[mn2d_gap],gap2D[mn2d_gap],'C0s',label='non-relaxed',alpha=0.7)
ax[1].set_ylabel(r'$M_{sat}/M_{BCG}$')
ax[1].set_xlabel(u'$\log(M_{200})$')

f.savefig(plotspath+'ratio.pdf',bbox_inches='tight')

f, ax = plt.subplots(2,1, figsize=(4,6), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].plot(lM[mo_dv],DV[mo_dv],'C3o',label='relaxed',alpha=0.7)
ax[0].plot(lM[mn_dv],DV[mn_dv],'C0s',label='non-relaxed',alpha=0.7)
ax[0].set_ylabel(r'$\Delta V$ [km/s]')
ax[0].set_xlabel(u'$\log(M_{200})$')


ax[1].plot(lMp[mo2d_dv],DV2D[mo2d_dv],'C3o',label='relaxed',alpha=0.7)
ax[1].plot(lMp[mn2d_dv],DV2D[mn2d_dv],'C0s',label='non-relaxed',alpha=0.7)
ax[1].set_ylabel(r'$\Delta V$ [km/s]')
ax[1].set_xlabel(u'$\log(M_{200})$')


f.savefig(plotspath+'DV.pdf',bbox_inches='tight')

f, ax = plt.subplots(2,1, figsize=(4,6), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].plot(lM[mo_off],off[mo_off],'C3o',label='relaxed',alpha=0.7)
ax[0].plot(lM[mn_off],off[mn_off],'C0s',label='non-relaxed',alpha=0.7)
ax[0].text(14,300,'3D')
ax[0].set_ylabel(r'$D_{offset}$ [kpc]')
ax[0].set_xlabel(u'$\log(M_{200})$')


ax[1].plot(lMp[mo2d_off],off2D[mo2d_off],'C3o',label='relaxed',alpha=0.7)
ax[1].plot(lMp[mn2d_off],off2D[mn2d_off],'C0s',label='non-relaxed',alpha=0.7)
ax[1].text(14,270,'2D')
ax[1].set_ylabel(r'$D_{offset}$ [kpc]')
ax[1].set_xlabel(u'$\log(M_{200})$')


f.savefig(plotspath+'Doff.pdf',bbox_inches='tight')
