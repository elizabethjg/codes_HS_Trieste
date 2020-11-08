import sys
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from scipy import stats
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

path = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/Trieste/catalog/'


# LOAD CATALOGS

gral  = np.loadtxt(path+'gral_091.dat',skiprows=10).T
gal   = np.loadtxt(path+'glxs_091.dat',skiprows=14).T
stars = np.loadtxt(path+'stars_091.dat',skiprows=14).T
dm    = np.loadtxt(path+'dm_091.dat',skiprows=13).T

# COMPUTE PARAMETERS

M200 = gral[9]
lM200 = np.log10(gral[9])

a_dm = dm[2]
b_dm = dm[3]
c_dm = dm[4]
qx_dm = dm[5]
qy_dm = dm[6]
qz_dm = dm[7]

T_dm = (a_dm**2 - b_dm**2)/(a_dm**2 - c_dm**2)
S_dm = c_dm/a_dm

q_dm = np.concatenate((qx_dm,qy_dm,qz_dm))

a_dm0 = dm[17]
b_dm0 = dm[18]
c_dm0 = dm[19]

qx_dm0 = dm[20]
qy_dm0 = dm[21]
qz_dm0 = dm[22]

T_dm0 = (a_dm0**2 - b_dm0**2)/(a_dm0**2 - c_dm0**2)
S_dm0 = c_dm0/a_dm0

q_dm0 = np.concatenate((qx_dm0,qy_dm0,qz_dm0))

# T bins plots

bined = stats.binned_statistic(lM200,T_dm,statistic='median', bins=10)
lM200_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
T_b     = bined.statistic

plt.figure()
plt.plot(lM200,T_dm,'C0o',alpha=0.5)
plt.plot(lM200_b,T_b,'C0')

bined = stats.binned_statistic(lM200,T_dm0,statistic='median', bins=10)
lM200_b0 = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
T_b0     = bined.statistic

plt.plot(lM200,T_dm0,'C0x',alpha=0.5)
plt.plot(lM200_b0,T_b0,'C0--')
plt.xlabel('$\log M_{200}$')
plt.ylabel('$T$')

plt.figure()
plt.plot(lM200,T_dm/T_dm0,'o',alpha=0.5)
plt.plot([14,15.65],[1,1],'C7--')
plt.xlabel('$\log M_{200}$')
plt.ylabel('$T/T_0$')

# S bins plots

bined = stats.binned_statistic(lM200,S_dm,statistic='median', bins=10)
lM200_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
S_b     = bined.statistic

plt.figure()
plt.plot(lM200,S_dm,'C1o',alpha=0.5)
plt.plot(lM200_b,S_b,'C1')

bined = stats.binned_statistic(lM200,S_dm0,statistic='median', bins=10)
lM200_b0 = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
S_b0     = bined.statistic

plt.plot(lM200,S_dm0,'C1x',alpha=0.5)
plt.plot(lM200_b0,S_b0,'C1--')
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')

plt.figure()
plt.plot(lM200,S_dm/S_dm0,'C1o',alpha=0.5)
plt.plot([14,15.65],[1,1],'C7--')
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S/S_0$')

# q bins plots

lM200 = np.array((lM200.tolist())*3)

bined = stats.binned_statistic(lM200,q_dm,statistic='median', bins=10)
lM200_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
q_b     = bined.statistic

plt.figure()
plt.plot(lM200,q_dm,'C2o',alpha=0.5)
plt.plot(lM200_b,q_b,'C2')

bined = stats.binned_statistic(lM200,q_dm0,statistic='median', bins=10)
lM200_b0 = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
q_b0     = bined.statistic

plt.plot(lM200,q_dm0,'C2x',alpha=0.5)
plt.plot(lM200_b0,q_b0,'C2--')
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q$')

plt.figure()
plt.plot(lM200,q_dm/q_dm0,'C2o',alpha=0.5)
plt.plot([14,15.65],[1,1],'C2--')
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q/q_0$')

