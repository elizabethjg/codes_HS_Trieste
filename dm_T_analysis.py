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

path_plots = '../plots/T_plots/DM/'

os.system('mkdir '+path_plots)

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

dm = DarkMatter()
dm0 = DarkMatter(False)

T2D  = np.array((dm.T.tolist())*3)
T2D0 = np.array((dm0.T.tolist())*3)

Nh  = Galaxias(200).N
nh  = Galaxias(200).n
Nh_x = Galaxias(200,True).N
nh_x = Galaxias(200,True).n

mgap = binned(lM,gap)[-1]
mgap2D = binned(lMp,gap2D)[-1]

mold = ~mgap
mold2D = ~mgap2D

mnew = mgap
mnew2D = mgap2D

mpro = dm.T > 0.5
mobl = ~mpro

mpro2D = T2D > 0.5
mobl2D = ~mpro2D



# '''
plt.figure()
plt.hist(dm.S[mpro],histtype='step',color='C1')
plt.hist(dm.S[mobl],histtype='step',color='C2')
plt.xlabel('$S$')
plt.ylabel('$N$')
plt.savefig(path_plots+'Sdist.png')



plt.figure()
plt.plot(lM,dm.S,'C7.')
plot_binned(lM,dm.S,'Dark matter','k',nbins=5)
plt.plot(lM[mpro],dm.S[mpro],'C1.')
plot_binned(lM[mpro],dm.S[mpro],'Prolate','C1',nbins=5)
plot_binned(lM[mobl],dm.S[mobl],'Oblate','C2',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')
plt.legend()
plt.savefig(path_plots+'S_DM_'+indicator+'.png')


plt.figure()
plt.hist(dm.q[mpro2D],histtype='step',color='C1')
plt.hist(dm.q[mobl2D],histtype='step',color='C2')
plt.xlabel('$q$')
plt.ylabel('$N$')
plt.savefig(path_plots+'qdist.png')


plt.figure()
plt.plot(lMp,dm.q,'C7.')
plot_binned(lMp,dm.q,'Dark matter','k',nbins=5)
plt.plot(lMp[mpro2D],dm.q[mpro2D],'C1.')
plot_binned(lMp[mold2D],dm.q[mold2D],'Prolate','C1',nbins=5)
plot_binned(lMp[mobl2D],dm.q[mobl2D],'Oblate','C2',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q$')
plt.savefig(path_plots+'q_DM.png')




# ''' DM whitout subhalos

plt.figure()
plt.hist(dm0.S[mpro],histtype='step',color='C1')
plt.hist(dm0.S[mobl],histtype='step',color='C2')
plt.xlabel('$S$')
plt.ylabel('$N$')
plt.savefig(path_plots+'Sdist_halo.png')



plt.figure()
plt.plot(lM,dm0.S,'C7.')
plot_binned(lM,dm0.S,'Dark matter','k',nbins=5)
plt.plot(lM[mpro],dm0.S[mpro],'C1.')
plot_binned(lM[mpro],dm0.S[mpro],'Prolate','C1',nbins=5)
plot_binned(lM[mobl],dm0.S[mobl],'Oblate','C2',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')
plt.legend()
plt.savefig(path_plots+'S_halo.png')


plt.figure()
plt.hist(dm0.q[mpro2D],histtype='step',color='C1')
plt.hist(dm0.q[mobl2D],histtype='step',color='C2')
plt.xlabel('$q$')
plt.ylabel('$N$')
plt.savefig(path_plots+'qdist_halo.png')


plt.figure()
plt.plot(lMp,dm0.q,'C7.')
plot_binned(lMp,dm0.q,'Dark matter','k',nbins=5)
plt.plot(lMp[mpro2D],dm0.q[mpro2D],'C1.')
plot_binned(lMp[mold2D],dm0.q[mold2D],'Prolate','C1',nbins=5)
plot_binned(lMp[mobl2D],dm0.q[mobl2D],'Oblate','C2',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q$')
plt.savefig(path_plots+'q_halo.png')



