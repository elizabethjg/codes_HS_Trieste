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


def plot_dm(radio,radio_dm,mask,mask2D,plot = False,sample = 'all'):    
    

    dm = DarkMatter(radio_dm)
    st = Stars(radio)

    pT = np.round(pearsonr(st.T[mask],dm.T[mask])[0],3)
    pS = np.round(pearsonr(st.S[mask],dm.S[mask])[0],3)
    pq = np.round(pearsonr(st.q[mask2D],dm.q[mask2D])[0],3)
    

    if plot: 

        path_plots = '../plots/correl_dM_stars/'
        
       
        plt.figure()
        
        plt.plot(dm.S[mask],st.S[mask],'C7.')
        plot_binned(dm.S[mask],st.S[mask],'Dark matter','k',nbins=5)
        plt.plot([0,1],[0,1],'C7--')
        plt.axis([0.,1,0.,1])
        plt.xlabel('$S_{DM}(R'+str(radio_dm)+')$')
        plt.ylabel('$S\star(R'+str(radio)+')$')
        plt.savefig(path_plots+'S_DM_R'+str(radio_dm)+'_stars_R'+str(radio)+'_'+sample+'.png')
        
        plt.figure()
        
        plt.plot(dm.T[mask],st.T[mask],'C7.')
        plot_binned(dm.T[mask],st.T[mask],'Dark matter','k',nbins=5)
        plt.plot([0,1],[0,1],'C7--')
        plt.axis([0.,1,0.,1])
        plt.xlabel('$T_{DM}(R'+str(radio_dm)+')$')
        plt.ylabel('$T\star(R'+str(radio)+')$')
        plt.savefig(path_plots+'T_DM_R'+str(radio_dm)+'_stars_R'+str(radio)+'_'+sample+'.png')

        plt.figure()
        
        plt.plot(dm.q[mask2D],st.q[mask2D],'C7.')
        plot_binned(dm.q[mask2D],st.q[mask2D],'Dark matter','k',nbins=5)
        plt.plot([0,1],[0,1],'C7--')
        plt.axis([0.,1,0.,1])
        plt.xlabel('$q_{DM}(R'+str(radio_dm)+')$')
        plt.ylabel('$q\star(R'+str(radio)+')$')
        plt.savefig(path_plots+'q_DM_R'+str(radio_dm)+'_stars_R'+str(radio)+'_'+sample+'.png')

        
    return pT, pS, pq


options = [30,50,100,1000,500,200]
order   = np.arange(len(options))
Rlabels   = ['30kpc','50kpc','0.1R500','R1000','R500','R200']

path_plots = '../plots/correl_dM_stars/'

lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)



indicator = 'gap'


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


mgap = binned(lM,gap)[-1]
mgap2D = binned(lMp,gap2D)[-1]

mold = ~mgap
mold2D = ~mgap2D

mnew = mgap
mnew2D = mgap2D

mall   = np.ones(len(lM)).astype(bool)
mall2D = np.ones(len(lMp)).astype(bool)

mask = [mall,mold,mnew]
mask2D = [mall2D,mold2D,mnew2D]
samples = ['all','old','new']

for m in range(3):
    
    sample    = samples[m]
    
    RDM = np.array([])
    R   = np.array([])
    
    PT = np.array([])
    PS = np.array([])
    Pq = np.array([])

    
    
    for j in order:
        for i in order:
            pT, pS, pq = plot_dm(options[order[j]],options[order[i]],mask[m],mask2D[m],sample=samples[m])
        
            R = np.append(R,order[j])
            RDM = np.append(RDM,order[i])
        
            PT = np.append(PT,pT)
            PS = np.append(PS,pS)
            Pq = np.append(Pq,pq)
    
    
    
    
            
    plt.figure()
    for x in order:
        plt.plot(R[RDM==x],PT[RDM==x],label=Rlabels[x],lw = (x+3)/4.)
    plt.xticks(order,Rlabels)
    plt.xlabel('R')
    plt.ylabel('p_T(R)')
    plt.axis([-1,6,-1,1])
    plt.legend(loc=3,ncol=4)
    plt.savefig(path_plots+'pT_DM_stars_'+sample+'.png')
    
    plt.figure()
    for x in order:
        plt.plot(R[RDM==x],PS[RDM==x],label=Rlabels[x],lw = (x+3)/4.)
    plt.xticks(order,Rlabels)
    plt.xlabel('R')
    plt.ylabel('p_S(R)')
    plt.axis([-1,6,-1,1])
    plt.legend(loc=3,ncol=4)
    plt.savefig(path_plots+'pS_DM_stars_'+sample+'.png')
    
    plt.figure()
    for x in order:
        plt.plot(R[RDM==x],Pq[RDM==x],label=Rlabels[x],lw = (x+3)/4.)
    plt.xticks(order,Rlabels)
    plt.xlabel('R')
    plt.ylabel('p_q(R)')
    plt.axis([-1,6,-1,1])
    plt.legend(loc=3,ncol=4)
    plt.savefig(path_plots+'pq_DM_stars_'+sample+'.png')


dm = DarkMatter(200)
gx = Stars(radio=30)

mnew,mold,mnew2D,mold2D = newold()

plt.figure()

plt.plot(dm.S[mnew],gx.S[mnew],'C0.')
plt.plot(dm.S[mold],gx.S[mold],'C3.')
plot_binned(dm.S[mnew],gx.S[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm.S[mold],gx.S[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,0.9,0.2,0.8])
plt.xlabel('$S_{DM}(R200)$')
plt.ylabel('$S\star(30kpc)$')
plt.savefig(path_plots+'S_DM_R200_stars30.png')

plt.figure()

plt.plot(dm.T[mnew],gx.T[mnew],'C0.')
plt.plot(dm.T[mold],gx.T[mold],'C3.')
plot_binned(dm.T[mnew],gx.T[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm.T[mold],gx.T[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.,1,0.,1])
plt.xlabel('$T_{DM}(R200)$')
plt.ylabel('$T\star(30kpc)$')
plt.savefig(path_plots+'T_DM_R200_stars30.png')


plt.figure()

plt.plot(dm.q[mnew2D],gx.q[mnew2D],'C0.')
plt.plot(dm.q[mold2D],gx.q[mold2D],'C3.')
plot_binned(dm.q[mnew2D],gx.q[mnew2D],'Dark matter','C0',nbins=5)
plot_binned(dm.q[mold2D],gx.q[mold2D],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,1,0.3,1])
plt.xlabel('$q_{DM}(R200)$')
plt.ylabel('$q\star(30kpc)$')
plt.savefig(path_plots+'q_DM_R200_stars30.png')

dm = DarkMatter(200)
gx = Stars(radio=200)

mnew,mold,mnew2D,mold2D = newold()

plt.figure()

plt.plot(dm.S[mnew],gx.S[mnew],'C0.')
plt.plot(dm.S[mold],gx.S[mold],'C3.')
plot_binned(dm.S[mnew],gx.S[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm.S[mold],gx.S[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,0.9,0.2,0.8])
plt.xlabel('$S_{DM}(R200)$')
plt.ylabel('$S\star(R200)$')
plt.savefig(path_plots+'S_DM_R200_starsR200.png')

plt.figure()

plt.plot(dm.T[mnew],gx.T[mnew],'C0.')
plt.plot(dm.T[mold],gx.T[mold],'C3.')
plot_binned(dm.T[mnew],gx.T[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm.T[mold],gx.T[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.,1,0.,1])
plt.xlabel('$T_{DM}(R200)$')
plt.ylabel('$T\star(R200)$')
plt.savefig(path_plots+'T_DM_R200_starsR200.png')


plt.figure()

plt.plot(dm.q[mnew2D],gx.q[mnew2D],'C0.')
plt.plot(dm.q[mold2D],gx.q[mold2D],'C3.')
plot_binned(dm.q[mnew2D],gx.q[mnew2D],'Dark matter','C0',nbins=5)
plot_binned(dm.q[mold2D],gx.q[mold2D],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,1,0.3,1])
plt.xlabel('$q_{DM}(R200)$')
plt.ylabel('$q\star(R200)$')
plt.savefig(path_plots+'q_DM_R200_starsR200.png')


dm = DarkMatter(30)
gx = Stars(radio=30)

mnew,mold,mnew2D,mold2D = newold()

plt.figure()

plt.plot(dm.S[mnew],gx.S[mnew],'C0.')
plt.plot(dm.S[mold],gx.S[mold],'C3.')
plot_binned(dm.S[mnew],gx.S[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm.S[mold],gx.S[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,0.9,0.2,0.8])
plt.xlabel('$S_{DM}(30kpc)$')
plt.ylabel('$S\star(30kpc)$')
plt.savefig(path_plots+'S_DM_30_stars30.png')

plt.figure()

plt.plot(dm.T[mnew],gx.T[mnew],'C0.')
plt.plot(dm.T[mold],gx.T[mold],'C3.')
plot_binned(dm.T[mnew],gx.T[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm.T[mold],gx.T[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.,1,0.,1])
plt.xlabel('$T_{DM}(30kpc)$')
plt.ylabel('$T\star(30kpc)$')
plt.savefig(path_plots+'T_DM_30_stars30.png')


plt.figure()

plt.plot(dm.q[mnew2D],gx.q[mnew2D],'C0.')
plt.plot(dm.q[mold2D],gx.q[mold2D],'C3.')
plot_binned(dm.q[mnew2D],gx.q[mnew2D],'Dark matter','C0',nbins=5)
plot_binned(dm.q[mold2D],gx.q[mold2D],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,1,0.3,1])
plt.xlabel('$q_{DM}(30kpc)$')
plt.ylabel('$q\star(30kpc)$')
plt.savefig(path_plots+'q_DM_30_stars30.png')
