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
import pandas as pd
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




def plot_fig(x,y,xc,yc,ycM,ycm,nbins,ax,mask,style,ly='',lyc=''):
        
        x,y,xc,yc,ycM,ycm = x[mask],y[mask],xc[mask],yc[mask],ycM[mask],ycm[mask]
        
        X,Y,q25,q75,m = binned(xc,yc,nbins)
        ax.plot(X,Y,'k'+style,label=lyc)
        X,YM,q25,q75,m = binned(xc,ycM,nbins)
        ax.plot(X,YM,'k'+style,alpha=0.2)
        X,Ym,q25,q75,m = binned(xc,ycm,nbins)
        ax.plot(X,Ym,'k'+style,alpha=0.2)
        ax.fill_between(X,YM,Ym,color = 'k',alpha=0.05)
        X,q50,q25,q75,mz = binned(x,y,nbins)
        ax.plot(X,q50,'C3'+style,label=ly)
        ax.plot(X,q75,'C3'+style,alpha=0.2)
        ax.plot(X,q25,'C3'+style,alpha=0.2)
        ax.fill_between(X,q75,q25,color = 'C1',alpha=0.1)


plotspath = '../final_plots/'

C = Clusters()

DM1000 = DarkMatter(1000)
DM500  = DarkMatter(500)
DM200  = DarkMatter(200)

r200 = Random()
r1000 = Random(1000)
r500 = Random(500)

gx200 = Galaxias(radio=200)
gx1000 = Galaxias(radio=1000)
gx500 = Galaxias(radio=500)

gx_c = Galaxias(tipo='concentradas')
gx_e = Galaxias(tipo='extendidas')

mall    = np.ones(72).astype(bool)
mall2D  = np.ones(72*3).astype(bool)

m1000  = gx1000.N > 9
m1000p = np.array((m1000.tolist())*3)
m500  = gx500.N > 9
m500p = np.array((m500.tolist())*3)
mc     = gx_c.N > 9
mcp    = np.array((mc.tolist())*3)
me     = gx_e.N > 9
mep    = np.array((me.tolist())*3)




# PLOT formas


f, ax = plt.subplots(3,1, figsize=(6,12), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(np.log10(gx1000.N),gx1000.T/DM1000.T,
         np.log10(r1000.N),r1000.T50/DM1000.T,r1000.T75/DM1000.T,
         r1000.T25/DM1000.T,6,ax[0],m1000,'-')

plot_fig(np.log10(gx500.N),gx500.T/DM500.T,
         np.log10(r500.N),r500.T50/DM500.T,r500.T75/DM500.T,
         r500.T25/DM500.T,6,ax[0],m500,'--')

plot_fig(np.log10(gx200.N),gx200.T/DM200.T,
         np.log10(r200.N),r200.T50/DM200.T,r200.T75/DM200.T,
         r200.T25/DM200.T,6,ax[0],mall,'-.')

         
plot_fig(np.log10(gx1000.N),gx1000.S/DM1000.S,
         np.log10(r1000.N),r1000.S50/DM1000.S,r1000.S75/DM1000.S,
         r1000.S25/DM1000.S,6,ax[1],m1000,'-',lyc='$R < R_{1000}$')

plot_fig(np.log10(gx500.N),gx500.S/DM500.S,
         np.log10(r500.N),r500.S50/DM500.S,r500.S75/DM500.S,
         r500.S25/DM500.S,6,ax[1],m500,'--',lyc='$R < R_{500}$')

plot_fig(np.log10(gx200.N),gx200.S/DM200.S,
         np.log10(r200.N),r200.S50/DM200.S,r200.S75/DM200.S,
         r200.S25/DM200.S,6,ax[1],mall,'-.',lyc='$R < R_{200}$')

ax[1].legend(frameon=False)

plot_fig(np.log10(gx1000.n),gx1000.q/DM1000.q,
         np.log10(r1000.n),r1000.q50/DM1000.q,r1000.q75/DM1000.q,
         r1000.q25/DM1000.q,6,ax[2],m1000p,'-')

plot_fig(np.log10(gx500.n),gx500.q/DM500.q,
         np.log10(r500.n),r500.q50/DM500.q,r500.q75/DM500.q,
         r500.q25/DM500.q,6,ax[2],m500p,'--')

plot_fig(np.log10(gx200.n),gx200.q/DM200.q,
         np.log10(r200.n),r200.q50/DM200.q,r200.q75/DM200.q,
         r200.q25/DM200.q,6,ax[2],mall2D,'-.')


ax[0].set_xlim([1.0,2.8])

ax[0].set_ylim([0.6,1.86])
ax[1].set_ylim([0.4,1.15])
ax[2].set_ylim([0.6,1.15])

ax[2].set_xlabel('$\log(N)$')
ax[0].set_ylabel('$T/T_{DM}$')
ax[1].set_ylabel('$S/S_{DM}$')
ax[2].set_ylabel('$q/q_{DM}$')

ax[0].plot([0.5,3],[1,1],'C7')
ax[1].plot([0.5,3],[1,1],'C7')
ax[2].plot([0.5,3],[1,1],'C7')

f.savefig(plotspath+'shape_galaxies.pdf',bbox_inches='tight')

# PLOT concentradas / extendidas

f, ax = plt.subplots(3,1, figsize=(6,12), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(np.log10(gx_c.N),gx_c.T/DM200.T,
         np.log10(r200.N),r200.T50/DM200.T,r200.T75/DM200.T,
         r200.T25/DM200.T,6,ax[0],mc,'-')

plot_fig(np.log10(gx_e.N),gx_e.T/DM200.T,
         np.log10(r200.N),r200.T50/DM200.T,r200.T75/DM200.T,
         r200.T25/DM200.T,6,ax[0],me,'--')

         
plot_fig(np.log10(gx_c.N),gx_c.S/DM200.S,
         np.log10(r200.N),r200.S50/DM200.S,r200.S75/DM200.S,
         r200.S25/DM200.S,6,ax[1],mc,'-',ly='Concentrated')

plot_fig(np.log10(gx_e.N),gx_e.S/DM200.S,
         np.log10(r200.N),r200.S50/DM200.S,r200.S75/DM200.S,
         r200.S25/DM200.S,6,ax[1],me,'--',ly='Extended')

ax[1].legend(frameon=False)

plot_fig(np.log10(gx_c.n),gx_c.q/DM200.q,
         np.log10(r200.n),r200.q50/DM200.q,r200.q75/DM200.q,
         r200.q25/DM200.q,6,ax[2],mcp,'-')

plot_fig(np.log10(gx_e.n),gx_e.q/DM200.q,
         np.log10(r200.n),r200.q50/DM200.q,r200.q75/DM200.q,
         r200.q25/DM200.q,6,ax[2],mep,'--')



ax[0].set_xlim([1.0,2.8])

ax[0].set_ylim([0.6,1.86])
ax[1].set_ylim([0.4,1.15])
ax[2].set_ylim([0.6,1.15])

ax[2].set_xlabel('$\log(N)$')
ax[0].set_ylabel('$T/T_{DM}$')
ax[1].set_ylabel('$S/S_{DM}$')
ax[2].set_ylabel('$q/q_{DM}$')

ax[0].plot([0.5,3],[1,1],'C7')
ax[1].plot([0.5,3],[1,1],'C7')
ax[2].plot([0.5,3],[1,1],'C7')


f.savefig(plotspath+'shape_galaxies_ce_R200.pdf',bbox_inches='tight')

# '''

# PLOT relajados / no-relajados

f, ax = plt.subplots(3,1, figsize=(6,14), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(np.log10(gx200.N),gx200.T/DM200.T,
         np.log10(r200.N),r200.T50/DM200.T,r200.T75/DM200.T,
         r200.T25/DM200.T,6,ax[0],C.mo_gap,'-')

plot_fig(np.log10(gx200.N),gx200.T/DM200.T,
         np.log10(r200.N),r200.T50/DM200.T,r200.T75/DM200.T,
         r200.T25/DM200.T,6,ax[0],C.mn_gap,'--')

         
plot_fig(np.log10(gx200.N),gx200.S/DM200.S,
         np.log10(r200.N),r200.S50/DM200.S,r200.S75/DM200.S,
         r200.S25/DM200.S,6,ax[1],C.mo_gap,'-',ly='relaxed')

plot_fig(np.log10(gx200.N),gx200.S/DM200.S,
         np.log10(r200.N),r200.S50/DM200.S,r200.S75/DM200.S,
         r200.S25/DM200.S,6,ax[1],C.mn_gap,'--',ly='non-relaxed')

ax[1].legend(frameon=False)

plot_fig(np.log10(gx200.n),gx200.q/DM200.q,
         np.log10(r200.n),r200.q50/DM200.q,r200.q75/DM200.q,
         r200.q25/DM200.q,6,ax[2],C.mo2d_gap,'-')

plot_fig(np.log10(gx200.n),gx200.q/DM200.q,
         np.log10(r200.n),r200.q50/DM200.q,r200.q75/DM200.q,
         r200.q25/DM200.q,6,ax[2],C.mn2d_gap,'--')



ax[0].set_xlim([1.0,2.8])

ax[0].set_ylim([0.6,1.86])
ax[1].set_ylim([0.4,1.15])
ax[2].set_ylim([0.6,1.15])

ax[2].set_xlabel('$\log(N)$')
ax[0].set_ylabel('$T/T_{DM}$')
ax[1].set_ylabel('$S/S_{DM}$')
ax[2].set_ylabel('$q/q_{DM}$')

ax[0].plot([0.5,3],[1,1],'C7')
ax[1].plot([0.5,3],[1,1],'C7')
ax[2].plot([0.5,3],[1,1],'C7')


f.savefig(plotspath+'shape_galaxies_r_R200.pdf',bbox_inches='tight')


# '''


cgall = CorrelR(Galaxias,DarkMatter)
cgce = CorrelR(DarkMatter,DarkMatter)

tgx  = cgall.t3D
tgxp = cgall.t2D

tgxc  = cgce.t3D_gxc
tgxe  = cgce.t3D_gxe

tgxcp  = cgce.t2D_gxc
tgxep  = cgce.t2D_gxe


f, ax = plt.subplots(2,1, figsize=(6,8), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(np.log10(gx1000.N),tgx.T[-3],
         np.log10(r1000.N),r1000.tdm50,r1000.tdm75,
         r1000.tdm25,4,ax[0],m1000,'-',lyc='$R < R_{1000}$')

plot_fig(np.log10(gx500.N),tgx.T[-2],
         np.log10(r500.N),r500.tdm50,r500.tdm75,
         r500.tdm25,6,ax[0],m500,'--',lyc='$R < R_{500}$')

plot_fig(np.log10(gx200.N),tgx.T[-1],
         np.log10(r200.N),r200.tdm50,r200.tdm75,
         r200.tdm25,6,ax[0],mall,'-.',lyc='$R < R_{200}$')

ax[0].legend(frameon=False)

plot_fig(np.log10(gx1000.n),tgxp.T[-3],
         np.log10(r1000.n),r1000.t2Ddm50,r1000.t2Ddm75,
         r1000.t2Ddm25,4,ax[1],m1000p,'-')

plot_fig(np.log10(gx500.n),tgxp.T[-2],
         np.log10(r500.n),r500.t2Ddm50,r500.t2Ddm75,
         r500.t2Ddm25,6,ax[1],m500p,'--')

plot_fig(np.log10(gx200.n),tgxp.T[-1],
         np.log10(r200.n),r200.t2Ddm50,r200.t2Ddm75,
         r200.t2Ddm25,6,ax[1],mall2D,'-.')


ax[0].set_xlim([1.0,2.8])

ax[0].set_ylim([0.,60])
ax[1].set_ylim([0.,60])


ax[1].set_xlabel('$\log(N)$')
ax[0].set_ylabel(r'$\theta$')
ax[1].set_ylabel(r'$\theta^*$')


f.savefig(plotspath+'align_galaxies.pdf',bbox_inches='tight')

# ALIGN CON Y EXT

f, ax = plt.subplots(2,1, figsize=(6,8), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(np.log10(gx_c.N),tgxc.T[-1],
         np.log10(r200.N),r200.tdm50,r200.tdm75,
         r200.tdm25,5,ax[0],mc,'-',ly='concentrated')

plot_fig(np.log10(gx_e.N),tgxe.T[-1],
         np.log10(r200.N),r200.tdm50,r200.tdm75,
         r200.tdm25,5,ax[0],me,'--',ly='extended')


ax[0].legend(frameon=False)

plot_fig(np.log10(gx_c.n),tgxcp.T[-1],
         np.log10(r200.n),r200.t2Ddm50,r200.t2Ddm75,
         r200.t2Ddm25,5,ax[1],mcp,'-',ly='concentrated')

plot_fig(np.log10(gx_e.n),tgxep.T[-1],
         np.log10(r200.n),r200.t2Ddm50,r200.t2Ddm75,
         r200.t2Ddm25,5,ax[1],mep,'--',ly='extended')


ax[0].set_xlim([1.0,2.8])

ax[0].set_ylim([0.,60])
ax[1].set_ylim([0.,60])


ax[1].set_xlabel('$\log(N)$')
ax[0].set_ylabel(r'$\theta$')
ax[1].set_ylabel(r'$\theta^*$')


f.savefig(plotspath+'align_gx_re.pdf',bbox_inches='tight')
# '''

# ALIGN CON Y EXT

f, ax = plt.subplots(2,1, figsize=(6,9.5), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(np.log10(gx_c.N),tgxc.T[-3],
         np.log10(r1000.N),r1000.tdm50,r1000.tdm75,
         r200.tdm25,5,ax[0],mc*m1000,'-',ly='concentrated')

plot_fig(np.log10(gx_e.N),tgxe.T[-3],
         np.log10(r1000.N),r1000.tdm50,r1000.tdm75,
         r1000.tdm25,5,ax[0],me*m1000,'--',ly='extended')


ax[0].legend(frameon=False)

plot_fig(np.log10(gx_c.n),tgxcp.T[-3],
         np.log10(r1000.n),r1000.t2Ddm50,r1000.t2Ddm75,
         r1000.t2Ddm25,5,ax[1],mcp*m1000p,'-',ly='concentrated')

plot_fig(np.log10(gx_e.n),tgxep.T[-3],
         np.log10(r1000.n),r1000.t2Ddm50,r1000.t2Ddm75,
         r1000.t2Ddm25,5,ax[1],mep*m1000p,'--',ly='extended')


ax[0].set_xlim([1.0,2.8])

ax[0].set_ylim([0.,60])
ax[1].set_ylim([0.,60])


ax[1].set_xlabel('$\log(N)$')
ax[0].set_ylabel(r'$\theta$')
ax[1].set_ylabel(r'$\theta^*$')


f.savefig(plotspath+'align_gx_re_R1000.pdf',bbox_inches='tight')
# '''

f, ax = plt.subplots(2,1, figsize=(6,9.5), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(np.log10(gx_c.N),tgxc.T[-2],
         np.log10(r500.N),r500.tdm50,r500.tdm75,
         r200.tdm25,5,ax[0],mc*m500,'-',ly='concentrated')

plot_fig(np.log10(gx_e.N),tgxe.T[-2],
         np.log10(r500.N),r500.tdm50,r500.tdm75,
         r500.tdm25,5,ax[0],me*m500,'--',ly='extended')


ax[0].legend(frameon=False)

plot_fig(np.log10(gx_c.n),tgxcp.T[-2],
         np.log10(r500.n),r500.t2Ddm50,r500.t2Ddm75,
         r1000.t2Ddm25,5,ax[1],mcp*m500p,'-',ly='concentrated')

plot_fig(np.log10(gx_e.n),tgxep.T[-2],
         np.log10(r500.n),r500.t2Ddm50,r500.t2Ddm75,
         r500.t2Ddm25,5,ax[1],mep*m500p,'--',ly='extended')


ax[0].set_xlim([1.0,2.8])


ax[0].set_ylim([0.,60])
ax[1].set_ylim([0.,60])


ax[1].set_xlabel('$\log(N)$')
ax[0].set_ylabel(r'$\theta$')
ax[1].set_ylabel(r'$\theta^*$')


f.savefig(plotspath+'align_gx_re_R500.pdf',bbox_inches='tight')



# ALIGN RELAX

f, ax = plt.subplots(2,1, figsize=(6,9.5), sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

plot_fig(np.log10(gx200.N),tgxc.T[-1],
         np.log10(r200.N),r200.tdm50,r200.tdm75,
         r200.tdm25,5,ax[0],C.mo_gap,'-',ly='relaxed')

plot_fig(np.log10(gx_e.N),tgxe.T[-1],
         np.log10(r200.N),r200.tdm50,r200.tdm75,
         r200.tdm25,5,ax[0],C.mn_gap,'--',ly='non-relaxed')


ax[0].legend(frameon=False)

plot_fig(np.log10(gx_c.n),tgxcp.T[-1],
         np.log10(r200.n),r200.t2Ddm50,r200.t2Ddm75,
         r200.t2Ddm25,5,ax[1],C.mo2d_gap,'-',ly='concentrated')

plot_fig(np.log10(gx_e.n),tgxep.T[-1],
         np.log10(r200.n),r200.t2Ddm50,r200.t2Ddm75,
         r200.t2Ddm25,5,ax[1],C.mn2d_gap,'--',ly='extended')


ax[0].set_xlim([1.0,2.8])

ax[0].set_ylim([0.,60])
ax[1].set_ylim([0.,60])


ax[1].set_xlabel('$\log(N)$')
ax[0].set_ylabel(r'$\theta$')
ax[1].set_ylabel(r'$\theta^*$')


f.savefig(plotspath+'align_gx_r.pdf',bbox_inches='tight')
