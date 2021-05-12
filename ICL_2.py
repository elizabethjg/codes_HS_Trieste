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



plotspath = '../plots_ICL/'

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

t_dm = cosangle(DM.a2Dicl,DM1000.a2D[micl])[1]
t_dm_bcg = cosangle(DM.a2Dbcg,DM1000.a2D[micl])[1]

t_st = cosangle(ICL_.a2Dicl,st1000.a2D[micl])[1]
t_st_bcg = cosangle(ICL_.a2Dbcg,st100.a2D[micl])[1]


f, ax = plt.subplots(2,2, figsize=(10,10), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax2 = ax.flatten()

ax[0,0].plot(DM100.q[micl][mn*mbdm],qbcg_dm[mn*mbdm],'C0.')
ax[0,0].plot(DM100.q[micl][mo*mbdm],qbcg_dm[mo*mbdm],'.',c='sienna')
ax[0,1].plot(DM1000.q[micl][mn],qdm[mn],'C0.')
ax[0,1].plot(DM1000.q[micl][mo],qdm[mo],'.',c='sienna')

ax[1,0].plot(st100.q[micl][mn*mb],qbcg[mn*mb],'C0*')
ax[1,0].plot(st100.q[micl][mo*mb],qbcg[mo*mb],'*',c='sienna')
ax[1,1].plot(st1000.q[micl][mn],qicl[mn],'C0*')
ax[1,1].plot(st1000.q[micl][mo],qicl[mo],'*',c='sienna')

for j in range(len(ax2)):
    ax2[j].plot([0.4,1.0],[0.4,1.0],'C7--')
    
ax[1,0].set_xlabel('$q (R < 0.1R_{500})$')
ax[1,1].set_xlabel('$q (R < R_{1000})$')

ax[0,0].set_ylabel('$q^*_{DM} $')
ax[1,0].set_ylabel('$q^*\star$')

plt.savefig(plotspath+'q_icl_scatter_gap.pdf',bbox_inches='tight')

f, ax = plt.subplots(2,2, figsize=(10,10), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax2 = ax.flatten()

ax[0,0].hist((DM100.q[micl][mn*mbdm]/qbcg_dm[mn*mbdm]),histtype='step',color='C0',density=True)
ax[0,0].hist((DM100.q[micl][mo*mbdm]/qbcg_dm[mo*mbdm]),histtype='step',color='peru',density=True)
ax[0,1].hist((DM1000.q[micl][mn]/qdm[mn]),histtype='step',color='C0',density=True)
ax[0,1].hist((DM1000.q[micl][mo]/qdm[mo]),histtype='step',color='peru',density=True)

ax[0,0].plot(np.mean(DM100.q[micl][mn*mbdm]/qbcg_dm[mn*mbdm]),4.1,'C0o')
ax[0,0].plot(np.mean(DM100.q[micl][mo*mbdm]/qbcg_dm[mo*mbdm]),4,'o',c='sienna')
ax[0,1].plot(np.mean(DM1000.q[micl][mn]/qdm[mn]),4.1,'C0o')
ax[0,1].plot(np.mean(DM1000.q[micl][mo]/qdm[mo]),4,'o',c='sienna')

ax[0,0].errorbar(np.mean(DM100.q[micl][mn*mbdm]/qbcg_dm[mn*mbdm]),4.1,xerr=np.std(DM100.q[micl][mn*mbdm]/qbcg_dm[mn*mbdm]),ecolor='C0',fmt='none')
ax[0,0].errorbar(np.mean(DM100.q[micl][mo*mbdm]/qbcg_dm[mo*mbdm]),4,xerr=np.std(DM100.q[micl][mo*mbdm]/qbcg_dm[mo*mbdm]),ecolor='peru',fmt='none')
ax[0,1].errorbar(np.mean(DM1000.q[micl][mn]/qdm[mn]),4.1,xerr=np.std(DM1000.q[micl][mn]/qdm[mn]),ecolor='C0',fmt='none')
ax[0,1].errorbar(np.mean(DM1000.q[micl][mo]/qdm[mo]),4,xerr=np.std(DM1000.q[micl][mo]/qdm[mo]),ecolor='peru',fmt='none')


ax[1,0].hist((st100.q[micl][mn*mb]/qbcg[mn*mb]),histtype='step',color='C0',density=True)
ax[1,0].hist((st100.q[micl][mo*mb]/qbcg[mo*mb]),histtype='step',color='peru',density=True)
ax[1,1].hist((st1000.q[micl][mn]/qicl[mn]),histtype='step',color='C0',density=True)
ax[1,1].hist((st1000.q[micl][mo]/qicl[mo]),histtype='step',color='peru',density=True)

ax[1,0].plot(np.mean(st100.q[micl][mn*mb]/qbcg[mn*mb]),4.1,'C0o')
ax[1,0].plot(np.mean(st100.q[micl][mo*mb]/qbcg[mo*mb]),4,'o',c='sienna')
ax[1,1].plot(np.mean(st1000.q[micl][mn]/qicl[mn]),4.1,'C0o')
ax[1,1].plot(np.mean(st1000.q[micl][mo]/qicl[mo]),4,'o',c='sienna')

ax[1,0].errorbar(np.mean(st100.q[micl][mn*mb]/qbcg[mn*mb]),4.1,xerr=np.std(st100.q[micl][mn*mb]/qbcg[mn*mb]),ecolor='C0',fmt='none')
ax[1,0].errorbar(np.mean(st100.q[micl][mo*mb]/qbcg[mo*mb]),4,xerr=np.std(st100.q[micl][mo*mb]/qbcg[mo*mb]),ecolor='peru',fmt='none')
ax[1,1].errorbar(np.mean(st1000.q[micl][mn]/qicl[mn]),4.1,xerr=np.std(st1000.q[micl][mn]/qicl[mn]),ecolor='C0',fmt='none')
ax[1,1].errorbar(np.mean(st1000.q[micl][mo]/qicl[mo]),4,xerr=np.std(st1000.q[micl][mo]/qicl[mo]),ecolor='peru',fmt='none')

ax[0,0].text(0.8,3,'DM')
ax[0,1].text(0.8,3,'DM')
ax[1,0].text(0.8,3,'Stars')
ax[1,1].text(0.8,3,'Stars')
    
ax[0,0].set_ylabel('$n$')
ax[1,0].set_ylabel('$n$')

ax[1,0].set_xlabel('$q/q^* (R < 0.1R_{500})$')
ax[1,1].set_xlabel('$q/q^*(R < R_{1000})$')

plt.savefig(plotspath+'q_dist_scatter_gap.pdf',bbox_inches='tight')


##### Alignment

f, ax = plt.subplots(2,2, figsize=(10,10), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax2 = ax.flatten()

ax[0,0].hist(t_dm_bcg[mbdm*mn],histtype='step',color='C0',density=True)
ax[0,0].hist(t_dm_bcg[mbdm*mo],histtype='step',color='peru',density=True)
ax[0,1].hist(t_dm[mn],histtype='step',color='C0',density=True)
ax[0,1].hist(t_dm[mo],histtype='step',color='peru',density=True)

ax[0,0].plot(np.mean(t_dm_bcg[mbdm*mn]),0.13,'C0o')
ax[0,0].plot(np.mean(t_dm_bcg[mbdm*mo]),0.12,'o',c='sienna')
ax[0,1].plot(np.mean(t_dm[mn]),0.13,'C0o')
ax[0,1].plot(np.mean(t_dm[mo]),0.12,'o',c='sienna')

ax[0,0].errorbar(np.mean(t_dm_bcg[mbdm*mn]),0.13,xerr=np.std(t_dm_bcg[mbdm*mn]),ecolor='C0',fmt='none')
ax[0,0].errorbar(np.mean(t_dm_bcg[mbdm*mo]),0.12,xerr=np.std(t_dm_bcg[mbdm*mo]),ecolor='peru',fmt='none')
ax[0,1].errorbar(np.mean(t_dm[mn]),0.13,xerr=np.std(t_dm[mn]),ecolor='C0',fmt='none')
ax[0,1].errorbar(np.mean(t_dm[mo]),0.12,xerr=np.std(t_dm[mo]),ecolor='peru',fmt='none')


ax[1,0].hist((t_st_bcg[mb*mn]),histtype='step',color='C0',density=True)
ax[1,0].hist((t_st_bcg[mb*mo]),histtype='step',color='peru',density=True)
ax[1,1].hist((t_st[mn]),histtype='step',color='C0',density=True)
ax[1,1].hist((t_st[mo]),histtype='step',color='peru',density=True)

ax[1,0].plot(np.mean(t_st_bcg[mb*mn]),0.13,'C0o')
ax[1,0].plot(np.mean(t_st_bcg[mb*mo]),0.12,'o',c='sienna')
ax[1,1].plot(np.mean(t_st[mn]),0.13,'C0o')
ax[1,1].plot(np.mean(t_st[mo]),0.12,'o',c='sienna')

ax[1,0].errorbar(np.mean(t_st_bcg[mb*mn]),0.13,xerr=np.std(t_st_bcg[mb*mn]),ecolor='C0',fmt='none')
ax[1,0].errorbar(np.mean(t_st_bcg[mb*mo]),0.12,xerr=np.std(t_st_bcg[mb*mo]),ecolor='peru',fmt='none')
ax[1,1].errorbar(np.mean(t_st[mn]),0.13,xerr=np.std(t_st[mn]),ecolor='C0',fmt='none')
ax[1,1].errorbar(np.mean(t_st[mo]),0.12,xerr=np.std(t_st[mo]),ecolor='peru',fmt='none')

ax[0,0].text(60,0.1,'DM')
ax[0,1].text(60,0.1,'DM')
ax[1,0].text(60,0.1,'Stars')
ax[1,1].text(60,0.1,'Stars')
    
ax[0,0].set_ylabel('$n$')
ax[1,0].set_ylabel('$n$')

ax[1,0].set_xlabel(r'$\theta (R < 0.1R_{500})$')
ax[1,1].set_xlabel(r'$\theta (R < R_{1000})$')

plt.savefig(plotspath+'theta_dist_scatter_gap.pdf',bbox_inches='tight')


#####################

f, ax = plt.subplots(2,1, figsize=(7,8), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(lt[micl][mbdm],DM100.q[micl][mbdm]/qbcg_dm[mbdm],5,ax=ax[0],color='k')
plot_fig(lt[micl][mb],st100.q[micl][mb]/qbcg[mb],5,ax=ax[0],color='C9')

plot_fig(lt[micl],DM1000.q[micl]/qdm,5,ax=ax[1],color='k')
plot_fig(lt[micl],st1000.q[micl]/qicl,5,ax=ax[1],color='C9')
    
ax[1].set_xlabel('look-back time [Gyr]')


ax[0].set_ylabel('$q/q^* (R < 0.1R_{500})$')
ax[1].set_ylabel('$q/q^*(R < R_{1000})$')

plt.savefig(plotspath+'q_icl_time.pdf',bbox_inches='tight')

#####################

f, ax = plt.subplots(2,1, figsize=(7,8), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(lt[micl][mbdm],t_dm_bcg[mbdm],5,ax=ax[0],color='k')
plot_fig(lt[micl][mb],t_st_bcg[mb],5,ax=ax[0],color='C9')

plot_fig(lt[micl],t_dm,5,ax=ax[1],color='k')
plot_fig(lt[micl],t_st,5,ax=ax[1],color='C9')
    
ax[1].set_xlabel('look-back time [Gyr]')


ax[0].set_ylabel(r'$\theta (R < 0.1R_{500})$')
ax[1].set_ylabel(r'$\theta (R < R_{1000})$')

plt.savefig(plotspath+'align_icl_time.pdf',bbox_inches='tight')

#####################

f, ax = plt.subplots(2,1, figsize=(7,8), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(lt[micl][mbdm],DM100.q[micl][mbdm]/qbcg_dm[mbdm],5,ax=ax[0],color='k')
plot_fig(lt[micl][mb],st100.q[micl][mb]/qbcg[mb],5,ax=ax[0],color='C9')

plot_fig(lt[micl],DM1000.q[micl]/qdm,5,ax=ax[1],color='k')
plot_fig(lt[micl],st1000.q[micl]/qicl,5,ax=ax[1],color='C9')
    
ax[1].set_xlabel('look-back time [Gyr]')


ax[0].set_ylabel('$q/q^* (R < 0.1R_{500})$')
ax[1].set_ylabel('$q/q^*(R < R_{1000})$')

plt.savefig(plotspath+'q_icl_time.pdf',bbox_inches='tight')

#####################

f, ax = plt.subplots(2,1, figsize=(7,8), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(lt[micl][mbdm],t_dm_bcg[mbdm],5,ax=ax[0],color='k')
plot_fig(lt[micl][mb],t_st_bcg[mb],5,ax=ax[0],color='C9')

plot_fig(lt[micl],t_dm,5,ax=ax[1],color='k')
plot_fig(lt[micl],t_st,5,ax=ax[1],color='C9')
    
ax[1].set_xlabel('look-back time [Gyr]')


ax[0].set_ylabel(r'$\theta (R < 0.1R_{500})$')
ax[1].set_ylabel(r'$\theta (R < R_{1000})$')

plt.savefig(plotspath+'align_icl_time.pdf',bbox_inches='tight')
