import sys
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from scipy import stats
from pylab import *
from main import *

cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
plt.rcParams['axes.grid'] =True
plt.rcParams['grid.color'] = '0.8'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

matplotlib.rcParams.update({'font.size': 14})

plotspath = '../final_plots/'

C = Clusters()

mtp = C.ltimep > 0
mt  = C.ltime > 0

f, ax = plt.subplots(2,1, figsize=(5,6.3), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)


ax[0].text(1,1.35,'3D')
ax[0].plot(C.ltime[mt*C.mo_gap],C.gap[mt*C.mo_gap],'.',color='sienna',label='relaxed')
ax[0].plot(C.ltime[mt*C.mn_gap],C.gap[mt*C.mn_gap],'C0.',label='non-relaxed')
ax[0].legend(frameon=False,loc=1)
plot_fig(C.ltime[mt],C.gap[mt],5,ax=ax[0],color='C7',label='R200')

ax[1].text(1,1.35,'2D')
ax[1].plot(C.ltimep[mtp*C.mn2d_gap],C.gapp[mtp*C.mn2d_gap],'C0.')
ax[1].plot(C.ltimep[mtp*C.mo2d_gap],C.gapp[mtp*C.mo2d_gap],'.',color='sienna')
ax[1].legend(frameon=False,loc=1)
plot_fig(C.ltimep[mtp],C.gapp[mtp],5,ax=ax[1],color='C7',label='R200')

ax[0].set_ylabel(r'$M_{sat}/M_{BCG}$')
ax[1].set_ylabel(r'$M_{sat}/M_{BCG}$')
ax[1].set_xlabel(r'look-back time [Gyr]')

f.savefig(plotspath+'relaxed_class.pdf',bbox_inches='tight')

DM  = DarkMatter(200)

f, ax = plt.subplots(2,1, figsize=(5,6.3), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)



ax[0].plot(C.ltime[mt*C.mo_gap],DM.T[mt*C.mo_gap],'.',color='sienna',label='relaxed')
ax[0].plot(C.ltime[mt*C.mn_gap],DM.T[mt*C.mn_gap],'C0.',label='non-relaxed')
plot_fig(C.ltime[mt],DM.T[mt],5,ax=ax[0],color='C7',label='R200')


ax[1].plot(C.ltime[mt*C.mo_gap],DM.S[mt*C.mo_gap],'.',color='sienna',label='relaxed')
ax[1].plot(C.ltime[mt*C.mn_gap],DM.S[mt*C.mn_gap],'C0.',label='non-relaxed')
ax[1].legend(frameon=False,loc=1)
plot_fig(C.ltime[mt],DM.S[mt],5,ax=ax[1],color='C7',label='R200')

ax[0].set_ylabel(r'$T$')
ax[1].set_ylabel(r'$S$')
ax[1].set_xlabel(r'look-back time [Gyr]')

f.savefig(plotspath+'3D_wtime.pdf',bbox_inches='tight')
