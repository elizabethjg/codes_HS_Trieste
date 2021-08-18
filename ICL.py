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
ICL_ = ICL(cutmax=True)
DM = ICL('dm',cutmax=True)

DM1000 = DarkMatter(1000)
DM100 = DarkMatter(100)
H1000 = DarkMatter(1000,False)
st1000 = Stars(1000)
st100 = Stars(100)

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

rlim = np.zeros(micl.sum())
rlim_dm = np.zeros(micl.sum())

for j in range(micl.sum()):
    rlim[j]    = ICL_.a[j][-1]/r1000[j]
    rlim_dm[j] = DM.a[j][-1]/r1000[j]


#####################

mold = lt[micl] > 4.9
mnew = lt[micl] <= 4.9


f, ax = plt.subplots(2,1, figsize=(5.5,8), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(qicl[mnew],qicl[mnew]/qdm[mnew],4,ax=ax[0],color='teal',label='late-formed')
plot_fig(qicl[mold],qicl[mold]/qdm[mold],4,ax=ax[0],color='C1',label='early-formed')

# ax[0].plot([0.3,0.9],[0.3,0.9],'C7--')
ax[0].set_ylim([0.76,1.12])
# ax[0].set_xlim([0.38,0.85])
ax[1].set_ylim([0,35])


plot_fig(qicl[mnew],t_dmp[mnew],4,ax=ax[1],color='teal',label='DM Isocountours shape')
plot_fig(qicl[mold],t_dmp[mold],4,ax=ax[1],color='C1',label='DM Isocountours shape')
    
ax[1].set_xlabel('')
ax[0].legend(loc=4,frameon=False,fontsize = 12)

ax[1].set_xlabel('$q_{ICL}$')
ax[0].set_ylabel('$q_{ICL}/q_{DM}$')
ax[1].set_ylabel(r'$\theta$')

plt.savefig(plotspath+'ICL.pdf',bbox_inches='tight')


#####################



