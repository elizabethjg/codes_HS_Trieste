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

path_plots = '../plots/age_plots/DM_'+indicator+'/'

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



plt.figure()
# plt.hist(dm.T-dm0.T,np.linspace(-0.2,0.2,25),histtype='step',color='C1')
plt.hist((dm.T-dm0.T)[mnew],np.linspace(-0.2,0.2,25),histtype='step',color='C0')
plt.hist((dm.T-dm0.T)[mold],np.linspace(-0.2,0.2,25),histtype='step',color='C3')
plt.xlabel('$T - T0$')
plt.ylabel('$N$')
plt.savefig(path_plots+'difT_dist.png')

plt.figure()
# plt.hist(dm.S-dm0.S,np.linspace(-0.2,0.2,25),histtype='step',color='C1')
plt.hist((dm.S-dm0.S)[mnew],np.linspace(-0.2,0.2,25),histtype='step',color='C0')
plt.hist((dm.S-dm0.S)[mold],np.linspace(-0.2,0.2,25),histtype='step',color='C3')
plt.xlabel('$S - S0$')
plt.ylabel('$N$')
plt.savefig(path_plots+'difS_dist.png')

plt.figure()
# plt.hist(dm.q-dm0.q,np.linspace(-0.2,0.2,25),histtype='step',color='C2')
plt.hist((dm.q-dm0.q)[mnew2D],np.linspace(-0.2,0.2,25),histtype='step',color='C0')
plt.hist((dm.q-dm0.q)[mold2D],np.linspace(-0.2,0.2,25),histtype='step',color='C3')
plt.xlabel('$q - q0$')
plt.ylabel('$N$')
plt.savefig(path_plots+'difq_dist.png')


plt.figure()
plt.plot(Nh,dm.T-dm0.T,'.')
plt.plot(Nh_x,dm.T-dm0.T,'.')
plt.plot([1,Nh.max()+1],[0,0],'C7--')
plt.ylabel('$T - T0$')
plt.xlabel('$N_{GAL}$')
plt.savefig(path_plots+'NdifT.png')

plt.figure()
plt.plot(Nh,dm.S-dm0.S,'.')
plt.plot(Nh_x,dm.S-dm0.S,'.')
plt.plot([1,Nh.max()+1],[0,0],'C7--')
plt.ylabel('$S - S0$')
plt.xlabel('$N_{GAL}$')
plt.savefig(path_plots+'NdifS.png')

plt.figure()
plt.plot(nh,dm.q-dm0.q,'.')
plt.plot(nh_x,dm.q-dm0.q,'.')
plt.plot([1,nh.max()+1],[0,0],'C7--')
plt.ylabel('$q - q0$')
plt.xlabel('$N_{GAL}$')
plt.savefig(path_plots+'Ndifq.png')

plt.figure()
plt.plot(gap,dm.T-dm0.T,'.')
plt.plot([gap.min(),gap.max()],[0,0],'C7--')
plt.ylabel('$T - T0$')
plt.xlabel(indicator)
plt.savefig(path_plots+indicator+'NdifT.png')

plt.figure()
plt.plot(gap,dm.S-dm0.S,'.')
plt.plot([gap.min(),gap.max()],[0,0],'C7--')
plt.ylabel('$S - S0$')
plt.xlabel(indicator)
plt.savefig(path_plots+indicator+'NdifS.png')

plt.figure()
plt.plot(gap2D,dm.q-dm0.q,'.')
plt.plot([gap2D.min(),gap2D.max()],[0,0],'C7--')
plt.ylabel('$q - q0$')
plt.xlabel(indicator)
plt.savefig(path_plots+indicator+'Ndifq.png')


m = lM < 18
mp = lMp < 18



# mold = gap < 0.2
# mnew = gap > 0.3

# mold2D = gap2D < 0.2
# mnew2D = gap2D > 0.3


# '''
plt.figure()
plt.hist(dm.S[mnew],histtype='step',color='C0')
plt.hist(dm.S[mold],histtype='step',color='C3')
plt.xlabel('$S$')
plt.ylabel('$N$')
plt.savefig(path_plots+'Sdist_'+indicator+'.png')



plt.figure()
plt.plot(lM,dm.S,'C7.')
plot_binned(lM,dm.S,'Dark matter','k',nbins=5)
plt.plot(lM[mold],dm.S[mold],'C3.')
plot_binned(lM[mold],dm.S[mold],'Dark matter - old','C3',nbins=5)
plot_binned(lM[mnew],dm.S[mnew],'Dark matter - new','C0',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')
plt.savefig(path_plots+'S_DM_'+indicator+'.png')

plt.figure()
plt.hist(dm.T[mnew],histtype='step',color='C0')
plt.hist(dm.T[mold],histtype='step',color='C3')
plt.xlabel('$T$')
plt.ylabel('$N$')
plt.savefig(path_plots+'Tdist_'+indicator+'.png')


plt.figure()
plt.plot(lM,dm.T,'C7.')
plot_binned(lM,dm.T,'Dark matter','k',nbins=5)
plt.plot(lM[mold],dm.T[mold],'C3.')
plot_binned(lM[mold],dm.T[mold],'Dark matter - old','C3',nbins=5)
plot_binned(lM[mnew],dm.T[mnew],'Dark matter - new','C0',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$T$')
plt.savefig(path_plots+'T_DM_'+indicator+'.png')

plt.figure()
plt.plot(lM,gap,'C7.')
plt.plot(lM[mold],gap[mold],'C3.')
plt.xlabel('$\log M_{200}$')
plt.ylabel(indicator)
plt.savefig(path_plots+indicator+'.png')

plt.figure()
plt.hist(dm.q[mnew2D],histtype='step',color='C0')
plt.hist(dm.q[mold2D],histtype='step',color='C3')
plt.xlabel('$q$')
plt.ylabel('$N$')
plt.savefig(path_plots+'qdist_'+indicator+'.png')


plt.figure()
plt.plot(lMp,dm.q,'C7.')
plot_binned(lMp,dm.q,'Dark matter','k',nbins=5)
plt.plot(lMp[mold2D],dm.q[mold2D],'C3.')
plot_binned(lMp[mold2D],dm.q[mold2D],'Dark matter - old','C3',nbins=5)
plot_binned(lMp[mnew2D],dm.q[mnew2D],'Dark matter - new','C0',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q$')
plt.savefig(path_plots+'q_DM_'+indicator+'.png')

plt.figure()
plt.plot(lMp,gap2D,'C7.')
plt.plot(lMp[mold2D],gap2D[mold2D],'C3.')
plt.xlabel('$\log M_{200}$')
plt.ylabel(indicator+'_2D')
plt.savefig(path_plots+indicator+'_2D.png')



# ''' DM whitout subhalos

plt.figure()
plt.plot(lM,dm0.S,'C7.')
plot_binned(lM,dm0.S,'Dark matter','k')
plt.plot(lM[mold],dm0.S[mold],'C3.')
plot_binned(lM[mold],dm0.S[mold],'Dark matter - old','C3',nbins=5)
plot_binned(lM[mnew],dm0.S[mnew],'Dark matter - new','C0',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$S$')
plt.savefig(path_plots+'S_DM_'+indicator+'_ssub.png')


plt.figure()
plt.plot(lM,dm0.T,'C7.')
plot_binned(lM,dm0.T,'Dark matter','k')
plt.plot(lM[mold],dm0.T[mold],'C3.')
plot_binned(lM[mold],dm0.T[mold],'Dark matter - old','C3',nbins=5)
plot_binned(lM[mnew],dm0.T[mnew],'Dark matter - new','C0',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$T$')
plt.savefig(path_plots+'T_DM_'+indicator+'_ssub.png')

plt.figure()
plt.plot(lMp,dm0.q,'C7.')
plot_binned(lMp,dm0.q,'Dark matter','k')
plt.plot(lMp[mold2D],dm0.q[mold2D],'C3.')
plot_binned(lMp[mold2D],dm0.q[mold2D],'Dark matter - old','C3',nbins=5)
plot_binned(lMp[mnew2D],dm0.q[mnew2D],'Dark matter - new','C0',nbins=5)
plt.ylim([0,1])
plt.xlabel('$\log M_{200}$')
plt.ylabel('$q$')
plt.savefig(path_plots+'q_DM_'+indicator+'_ssub.png')

plt.figure()
plt.plot(lMp,gap2D,'C7.')
plt.plot(lMp[mold2D],gap2D[mold2D],'C3.')
plt.xlabel('$\log M_{200}$')
plt.ylabel(indicator+'_2D')
plt.savefig(path_plots+indicator+'_2D.png')
