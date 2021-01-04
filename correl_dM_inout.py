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

gral  = np.loadtxt(path+'gral_091_2.dat').T

mnew,mold,mnew2D,mold2D = newold()

dm30 = DarkMatter(30)
dm200 = DarkMatter(200)

path_plots = '../plots/correl_dM_inout/'

plt.figure()
plt.plot(dm30.T[mnew],dm200.T[mnew],'C0.')
plt.plot(dm30.T[mold],dm200.T[mold],'C3.')
plot_binned(dm30.T[mnew],dm200.T[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm30.T[mold],dm200.T[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.xlabel('T_dm(30kpc)')
plt.ylabel('T_dm(R200)')
plt.savefig(path_plots+'T_inout.png')

plt.figure()
plt.plot(dm30.S[mnew],dm200.S[mnew],'C0.')
plt.plot(dm30.S[mold],dm200.S[mold],'C3.')
plot_binned(dm30.S[mnew],dm200.S[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm30.S[mold],dm200.S[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.6,1,0.4,1])
plt.xlabel('S_dm(30kpc)')
plt.ylabel('S_dm(R200)')
plt.savefig(path_plots+'S_inout.png')

plt.figure()
plt.plot(dm30.q[mnew2D],dm200.q[mnew2D],'C0.')
plt.plot(dm30.q[mold2D],dm200.q[mold2D],'C3.')
plot_binned(dm30.q[mnew2D],dm200.q[mnew2D],'Dark matter','C0',nbins=5)
plot_binned(dm30.q[mold2D],dm200.q[mold2D],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.7,1,0.4,1])
plt.xlabel('q_dm(30kpc)')
plt.ylabel('q_dm(R200)')
plt.savefig(path_plots+'q_inout.png')

lM1000 = np.log10(gral[7])
lM500 = np.log10(gral[8])
lM200 = np.log10(gral[9])
lM30 = np.log10(gral[10])
lM50 = np.log10(gral[11])
lM100 = np.log10(gral[12])

plt.figure()
plt.plot(lM30[mnew],lM200[mnew],'C0.')
plt.plot(lM30[mold],lM200[mold],'C3.')
plot_binned(lM30[mnew],lM200[mnew],'Dark matter','C0',nbins=5)
plot_binned(lM30[mold],lM200[mold],'Dark matter','C3',nbins=5)
plt.xlabel('logM(30kpc)')
plt.ylabel('logM(R200)')
plt.plot([11,16],[11,16],'C7--')
plt.axis([11.2,11.85,13.9,15.7])
plt.savefig(path_plots+'logM_inout.png')
