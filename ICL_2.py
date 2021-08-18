import sys
import time
import numpy as np
from pylab import *
from main import *
import os


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


matplotlib.rcParams.update({'font.size': 13})

plotspath = '../final_plots/'


C = Clusters()
ICL_ = ICL()
DM = ICL('dm')

DM1000 = DarkMatter(1000)
DM500 = DarkMatter(500)
DM200 = DarkMatter(200)
DM100 = DarkMatter(100)
DM50 = DarkMatter(50)
DM30 = DarkMatter(30)

H1000 = DarkMatter(1000,False)


st1000 = Stars(1000)
st500  = Stars(500)
st200  = Stars(200)
st100  = Stars(100)
st50   = Stars(50)
st30   = Stars(30)

lt =  C.ltimep

m = (C.sub == 0)
micl = np.array(m.tolist()*3)

mo = C.mo2d_gap[micl]
mn = C.mn2d_gap[micl]

qicl = ICL_.qicl
qbcg = ICL_.qbcg

qdm = DM.qicl
qbcg_dm = DM.qbcg

mbdm = DM.Rbcg > 0.
mb   = ICL_.Rbcg > 0.

r200  = C.Rp[micl,-1]
r1000  = C.Rp[micl,-3]


t_st = cosangle(st1000.a2D[micl],ICL_.a2Dicl)[1]
t_dmp = cosangle(DM.a2Dicl,ICL_.a2Dicl)[1]
t_dmp_bcg = cosangle(DM.a2Dbcg,ICL_.a2Dicl)[1]

t_dmdm = cosangle(DM.a2Dicl,DM1000.a2D[micl])[1]
t_dm = cosangle(ICL_.a2Dicl,DM1000.a2D[micl])[1]
t_dm_bcg = cosangle(ICL_.a2Dicl,DM100.a2D[micl])[1]
t_dmh = cosangle(DM.a2Dicl,H1000.a2D[micl])[1]

t_bcg = cosangle(ICL_.a2Dicl,st100.a2D[micl])[1]
t_bcg_dm = cosangle(DM.a2Dicl,DM100.a2D[micl])[1]

rlim = np.zeros(micl.sum())
rlim_dm = np.zeros(micl.sum())

for j in range(micl.sum()):
    rlim[j]    = ICL_.a[j][-1]/r1000[j]
    rlim_dm[j] = DM.a[j][-1]/r1000[j]

f, ax = plt.subplots(2,1, figsize=(5.5,8), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(lt[micl],qdm/DM1000.q[micl],4,color='k',ax=ax[0],label = 'DM Isocontour - DM Inertial')
# plot_fig(lt[micl],qdm/H1000.q[micl],5,color='teal',ax=ax[0],label = 'DM Iscontour - H Inertial')
plot_fig(lt[micl],qicl/st1000.q[micl],4,color='C1',ax=ax[0],label = 'ICL - Stars Inertial')
plot_fig(lt[micl][mbdm],(qbcg_dm/DM100.q[micl])[mbdm],4,color='k',ax=ax[0],style='--',label = 'DM Isocontour - DM Inertial')
plot_fig(lt[micl][mb],(qbcg/st100.q[micl])[mb],4,color='C1',style='--',ax=ax[0],label = 'ICL - Stars Inertial')
ax[0].set_ylabel(r'$q_{iso}/q_{ine}$')


plot_fig(lt[micl],t_dmdm,4,color='k',ax=ax[1],label = 'DM')
# plot_fig(lt[micl],t_dmh,5,color='teal',ax=ax[1],label = 'DM Iscontour - H Inertial')
plot_fig(lt[micl],t_st,4,color='C1',ax=ax[1],label = 'Stellar')
plot_fig(lt[micl],t_dmdm,4,color='k',ax=ax[1],label = 'ICM')
plot_fig(lt[micl][mbdm],t_bcg_dm[mbdm],4,color='k',ax=ax[1],style='--',label = 'BCG')
plot_fig(lt[micl][mb],t_bcg[mb],4,color='C1',ax=ax[1],style='--')
ax[1].set_xlabel(r'Formation time [Gyr]')
ax[1].set_ylabel(r'$\theta$')
ax[1].legend(loc=1,frameon=False,fontsize = 12,ncol=2)

plt.savefig(plotspath+'comparison_iso_iner_time.pdf',bbox_inches='tight')

f, ax = plt.subplots(2,1, figsize=(5.5,8), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(qdm,qdm/DM1000.q[micl],4,color='k',ax=ax[0],label = 'DM Isocontour - DM Inertial')
plot_fig(qicl[rlim>0.7],(qicl/st1000.q[micl])[rlim>0.7],4,color='C1',ax=ax[0],label = 'ICL - Stars Inertial')
plot_fig(qbcg_dm[mbdm],(qbcg_dm/DM100.q[micl])[mbdm],4,color='k',ax=ax[0],style='--',label = 'DM Isocontour - DM Inertial')
plot_fig(qbcg[mb],(qbcg/st100.q[micl])[mb],4,color='C1',style='--',ax=ax[0],label = 'ICL - Stars Inertial')
ax[0].set_ylabel(r'$q_{iso}/q_{ine}$')


plot_fig(qdm,t_dmdm,4,color='k',ax=ax[1],label = 'DM')
plot_fig(qicl[rlim>0.7],t_st[rlim>0.7],4,color='C1',ax=ax[1],label = 'Stellar')
plot_fig(qdm,t_dmdm,4,color='k',ax=ax[1],label = 'ICM')
plot_fig(qbcg_dm[mbdm],t_bcg_dm[mbdm],4,color='k',ax=ax[1],style='--',label = 'BCG')
plot_fig(qbcg[mb],t_bcg[mb],4,color='C1',ax=ax[1],style='--')
ax[1].set_xlabel(r'$q_{iso}$')
ax[1].set_ylabel(r'$\theta^{iso}_{ine}$')
ax[1].legend(loc=2,frameon=False,fontsize = 12,ncol=2)

plt.savefig(plotspath+'comparison_iso_iner.pdf',bbox_inches='tight')



# plt.figure(figsize=(5.5,4))
# plot_fig(lt[micl],rlim,5,color='C1',label = 'ICL')
# plot_fig(lt[micl],rlim_dm,5,color='teal',label = 'DM isodensity')
# plt.plot(lt[micl],rlim,'.',color='C1')
# plt.plot(lt[micl],rlim_dm,'.',color='teal')
# plt.xlabel(r'Formation time [Gyr]')
# plt.hist(rlim,10,histtype='step',color='C1',lw=2,label='ICL')
# plt.hist(rlim_dm,10,histtype='step',color='teal',lw=2,label = 'DM isodensity')
# plt.axvline(np.median(rlim),color='C1')
# plt.axvline(np.median(rlim_dm),color='teal')
# plt.legend(loc=1,frameon=False,fontsize = 12)
# plt.ylabel('$N$')
# plt.xlabel('$R_{MAX}/R_{1000}$')
# plt.savefig(plotspath+'RMAX_ICL.pdf',bbox_inches='tight')


