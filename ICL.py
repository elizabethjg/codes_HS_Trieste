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
from scipy.optimize import curve_fit
from scipy import integrate

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

def plot_fig(x,y,nbins,ax=plt,color = 'sienna', style = '',label=''):
                
        X,q50,q25,q75,mz = binned(x,y,nbins)
        ax.plot(X,q50,color+style,label=label)
        ax.plot(X,q75,color+style,alpha=0.2)
        ax.plot(X,q25,color+style,alpha=0.2)
        ax.fill_between(X,q75,q25,color = color,alpha=0.1)



C = Clusters()
ICL_ = ICL()
DM = ICL('dm')

DM1000 = DarkMatter(1000)
DM500  = DarkMatter(500)
DM200  = DarkMatter(200)
DM30  = DarkMatter(30)
DM50  = DarkMatter(50)
DM100  = DarkMatter(100)

st200  = Stars(200)
st1000 = Stars(1000)
st500  = Stars(500)
st30   = Stars(30)
st50   = Stars(50)
st100  = Stars(100)

lt =  np.array((C.ltime.tolist())*3) 

m = (C.sub == 0)
micl = (m.tolist()*3)

mo = C.mo2d_gap[micl]
mn = C.mn2d_gap[micl]

qicl = ICL_.qicl
qbcg = ICL_.qbcg

qdm = DM.qicl
qbcg_dm = DM.qbcg

mbdm = DM.Rbcg > 0.
mb   = ICL_.Rbcg > 0.

r200  = C.Rp[micl,-1]


t_dmp = cosangle(DM.a2Dicl,ICL_.a2Dicl)[1]
t_dmp_bcg = cosangle(DM.a2Dbcg,ICL_.a2Dicl)[1]

t_dm = cosangle(ICL_.a2Dicl,DM1000.a2D[micl])[1]
t_dm_bcg = cosangle(ICL_.a2Dicl,DM100.a2D[micl])[1]



#####################

f, ax = plt.subplots(2,1, figsize=(5,7), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)


plot_fig(lt[micl],qicl/qdm,5,ax=ax[0],color='teal')
plot_fig(lt[micl],qicl/DM1000.q[micl],5,ax=ax[0],color='C2',style='--')

plot_fig(lt[micl],t_dmp,5,ax=ax[1],color='teal')
plot_fig(lt[micl],t_dm,5,ax=ax[1],color='C2',style='--')

    
ax[1].set_xlabel('look-back time [Gyr]')


ax[0].set_ylabel('$q_{ICL}/q_{DM}$')
ax[1].set_ylabel(r'$\theta_{ICL}^{DM}$')

plt.savefig(plotspath+'ICL_time.pdf',bbox_inches='tight')

#####################

f, ax = plt.subplots(1,2, figsize=(7,5),sharey=True)
f.subplots_adjust(hspace=0,wspace=0)


ax[0].hist((qicl/qdm)[C.mo2d_gap[micl]],np.linspace(0.5,1,6,10),color='sienna',histtype='step')
ax[0].hist((qicl/qdm)[C.mn2d_gap[micl]],np.linspace(0.5,1,6,10),color='C0',histtype='step')

ax[1].hist(t_dmp[C.mo2d_gap[micl]],np.linspace(0,70,10),color='sienna',histtype='step')
ax[1].hist(t_dmp[C.mn2d_gap[micl]],np.linspace(0,70,10),color='C0',histtype='step')
    
ax[1].set_ylabel('$N$')


ax[0].set_xlabel('$q_{ICL}/q_{DM}$')
ax[1].set_xlabel(r'$\theta_{ICL}^{DM}$')

plt.savefig(plotspath+'ICL_dist.pdf',bbox_inches='tight')


'''
header = 'D \n proj: 0 xy, 1 xz, 2 yz \n (5) BCG properties: R [kpc], q=b/a, a_11, a_12, PA \n (5) same for ICL'


out = np.array([ICL_.D,ICL_.proj,
                ICL_.Rbcg*r200,ICL_.qbcg,
                ICL_.a2Dbcg[:,0],ICL_.a2Dbcg[:,1],ICL_.PAbcg,
                ICL_.Ricl*r200,ICL_.qicl,
                ICL_.a2Dicl[:,0],ICL_.a2Dicl[:,1],ICL_.PAicl])

out[np.isnan(out)] = -99.
np.savetxt('../ISOPHOTES_stars.list', out.T,fmt=['%5i']*2+['%7.1f']+['%12.5f']*4+['%7.1f']+['%12.5f']*4,header = header)

out_DM = np.array([DM.D,DM.proj,
                DM.Rbcg*r200,DM.qbcg,
                DM.a2Dbcg[:,0],DM.a2Dbcg[:,1],DM.PAbcg,
                DM.Ricl*r200,DM.qicl,
                DM.a2Dicl[:,0],DM.a2Dicl[:,1],DM.PAicl])

out_DM[np.isnan(out_DM)] = -99.

np.savetxt('../ISOPHOTES_DM.list', out_DM.T,fmt=['%5i']*2+['%7.1f']+['%12.5f']*4+['%7.1f']+['%12.5f']*4,header = header)
'''
