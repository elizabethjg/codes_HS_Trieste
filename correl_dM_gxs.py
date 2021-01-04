import sys
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from scipy import stats
from pylab import *
from main_newdata import *
from scipy.stats import pearsonr
import os
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

path = '../catalog/nuevosdats/'

gral  = np.loadtxt(path+'gral_nounb_091.dat').T


def plot_dm(tipo,radio,radio_dm,mask,mask2D,plot = False,sample = 'all'):    
    

    dm = DarkMatter(radio_dm)
    st = Galaxias(tipo = tipo, radio=radio)

    mN = st.N > 9
    mN2D = np.array((st.N.tolist())*3) > 9
    
    mask = mask*mN
    mask2D = mask2D*mN2D
    
    pT = np.round(pearsonr(st.T[mask],dm.T[mask])[0],3)
    pS = np.round(pearsonr(st.S[mask],dm.S[mask])[0],3)
    pq = np.round(pearsonr(st.q[mask2D],dm.q[mask2D])[0],3)
    

    if plot: 

        path_plots = '../plots_newdata/correl_dM_gxs/'
        
       
        plt.figure()
        plt.title(sample)
        plt.plot(dm.S[mask],st.S[mask],'C7.')
        plot_binned(dm.S[mask],st.S[mask],'Dark matter','k',nbins=5)
        plt.plot([0,1],[0,1],'C7--')
        plt.axis([0.,1,0.,1])
        plt.xlabel('$S_{DM}(R'+str(radio_dm)+')$')
        plt.ylabel('$S_{gxs}(R'+str(radio)+')$')
        plt.savefig(path_plots+'S_DM_R'+str(radio_dm)+'_gxs_'+tipo+'_'+sample+'_new.png')
        
        plt.figure()
        plt.title(sample)
        plt.plot(dm.T[mask],st.T[mask],'C7.')
        plot_binned(dm.T[mask],st.T[mask],'Dark matter','k',nbins=5)
        plt.plot([0,1],[0,1],'C7--')
        plt.axis([0.,1,0.,1])
        plt.xlabel('$T_{DM}(R'+str(radio_dm)+')$')
        plt.ylabel('$T_{gxs}(R'+str(radio)+')$')
        plt.savefig(path_plots+'T_DM_R'+str(radio_dm)+'_gxs_'+tipo+'_'+sample+'_new.png')

        plt.figure()
        plt.title(sample)
        plt.plot(dm.q[mask2D],st.q[mask2D],'C7.')
        plot_binned(dm.q[mask2D],st.q[mask2D],'Dark matter','k',nbins=5)
        plt.plot([0,1],[0,1],'C7--')
        plt.axis([0.,1,0.,1])
        plt.xlabel('$q_{DM}(R'+str(radio_dm)+')$')
        plt.ylabel('$q_{gxs}(R'+str(radio)+')$')
        plt.savefig(path_plots+'q_DM_R'+str(radio_dm)+'_gxs_'+tipo+'_'+sample+'_new.png')

        
    return pT, pS, pq


options = [30,50,100,1000,500,200]

tipos = ['all','concentradas','extendidas']
order   = np.arange(len(options))
Rlabels   = ['30kpc','50kpc','0.1R500','R1000','R500','R200']

path_plots = '../plots_newdata/correl_dM_gxs/'

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
    T   = np.array([])
    
    PT = np.array([])
    PS = np.array([])
    Pq = np.array([])

    
    
    for j in order[:3]:
        for i in order:
            pT, pS, pq = plot_dm(tipos[j],200,options[order[i]],mask[m],mask2D[m],sample=samples[m])
        
            T   = np.append(T,j)
            RDM = np.append(RDM,order[i])
        
            PT = np.append(PT,pT)
            PS = np.append(PS,pS)
            Pq = np.append(Pq,pq)
    
    
    
    
            
    plt.figure()
    for x in order:
        plt.plot(T[RDM==x],PT[RDM==x],label=Rlabels[x],lw = (x+3)/4.)
    plt.xticks(order[:3],tipos)
    plt.xlabel('Type')
    plt.ylabel('p_T(R)')
    plt.axis([-1,3,-1,1])
    plt.legend(loc=3,ncol=4)
    plt.savefig(path_plots+'pT_DM_gxs_'+sample+'_type_new.png')
    
    plt.figure()
    for x in order:
        plt.plot(T[RDM==x],PS[RDM==x],label=Rlabels[x],lw = (x+3)/4.)
    plt.xticks(order[:3],tipos)
    plt.xlabel('R')
    plt.ylabel('p_S(R)')
    plt.axis([-1,3,-1,1])
    plt.legend(loc=3,ncol=4)
    plt.savefig(path_plots+'pS_DM_gxs_'+sample+'_type_new.png')
    
    plt.figure()
    for x in order:
        plt.plot(T[RDM==x],Pq[RDM==x],label=Rlabels[x],lw = (x+3)/4.)
    plt.xticks(order[:3],tipos)
    plt.xlabel('R')
    plt.ylabel('p_q(R)')
    plt.axis([-1,3,-1,1])
    plt.legend(loc=3,ncol=4)
    plt.savefig(path_plots+'pq_DM_gxs_'+sample+'_type_new.png')

dm = DarkMatter(200)
gx = Galaxias(radio=200)

mnew,mold,mnew2D,mold2D = newold()

plt.figure()
plt.title(sample)
plt.plot(dm.S[mnew],gx.S[mnew],'C0.')
plt.plot(dm.S[mold],gx.S[mold],'C3.')
plot_binned(dm.S[mnew],gx.S[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm.S[mold],gx.S[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,0.9,0.2,0.8])
plt.xlabel('$S_{DM}(R200)$')
plt.ylabel('$S_{gxs}(R200)$')
plt.savefig(path_plots+'S_DM_R200_gxs_all_new.png')

plt.figure()
plt.title(sample)
plt.plot(dm.T[mnew],gx.T[mnew],'C0.')
plt.plot(dm.T[mold],gx.T[mold],'C3.')
plot_binned(dm.T[mnew],gx.T[mnew],'Dark matter','C0',nbins=5)
plot_binned(dm.T[mold],gx.T[mold],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.,1,0.,1])
plt.xlabel('$T_{DM}(R200)$')
plt.ylabel('$T_{gxs}(R200)$')
plt.savefig(path_plots+'T_DM_R200_gxs_all_new.png')


plt.figure()
plt.title(sample)
plt.plot(dm.q[mnew2D],gx.q[mnew2D],'C0.')
plt.plot(dm.q[mold2D],gx.q[mold2D],'C3.')
plot_binned(dm.q[mnew2D],gx.q[mnew2D],'Dark matter','C0',nbins=5)
plot_binned(dm.q[mold2D],gx.q[mold2D],'Dark matter','C3',nbins=5)
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,1,0.3,1])
plt.xlabel('$q_{DM}(R200)$')
plt.ylabel('$q_{gxs}(R200)$')
plt.savefig(path_plots+'q_DM_R200_gxs_all_new.png')


dm = DarkMatter(200)
gx = Galaxias(radio=200)
gx_c = Galaxias('concentradas')
gx_e = Galaxias('extendidas')

mnew,mold,mnew2D,mold2D = newold()

plt.figure()
plot_binned(dm.S[mnew],gx_c.S[mnew],'concentradas','C0',style='--',nbins=5)
plot_binned(dm.S[mnew],gx_e.S[mnew],'extendidas','C0',style=':',nbins=5)
plt.legend()
plot_binned(dm.S[mold],gx_e.S[mold],'extendidas','C3',style=':',nbins=5)
plot_binned(dm.S[mold],gx_c.S[mold],'concentradas','C3',style='--',nbins=5)
plt.xlabel('$S_{DM}(R200)$')
plt.ylabel('$S_{gxs}(R200)$')
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,0.9,0.2,0.8])
plt.savefig(path_plots+'S_concentradas_extendidas_new.png')

plt.figure()
plot_binned(dm.T[mnew],gx_c.T[mnew],'concentradas','C0',style='--',nbins=5)
plot_binned(dm.T[mnew],gx_e.T[mnew],'extendidas','C0',style=':',nbins=5)
plt.legend()
plot_binned(dm.T[mold],gx_e.T[mold],'extendidas','C3',style=':',nbins=5)
plot_binned(dm.T[mold],gx_c.T[mold],'concentradas','C3',style='--',nbins=5)
plt.xlabel('$T_{DM}(R200)$')
plt.ylabel('$T_{gxs}(R200)$')
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.,1.,0.,1.])
plt.savefig(path_plots+'T_concentradas_extendidas_new.png')

plt.figure()
plot_binned(dm.q[mnew2D],gx_c.q[mnew2D],'concentradas','C0',style='--',nbins=5)
plot_binned(dm.q[mnew2D],gx_e.q[mnew2D],'extendidas','C0',style=':',nbins=5)
plt.legend()
plot_binned(dm.q[mold2D],gx_e.q[mold2D],'extendidas','C3',style=':',nbins=5)
plot_binned(dm.q[mold2D],gx_c.q[mold2D],'concentradas','C3',style='--',nbins=5)
plt.xlabel('$q_{DM}(R200)$')
plt.ylabel('$q_{gxs}(R200)$')
plt.plot([0,1],[0,1],'C7--')
plt.axis([0.4,1,0.2,1])
plt.savefig(path_plots+'q_concentradas_extendidas_new.png')


plt.figure()
plot_binned(lM[mnew],gx_c.S[mnew],'concentradas','C0',style='--',nbins=5)
plot_binned(lM[mnew],gx_e.S[mnew],'extendidas','C0',style=':',nbins=5)
plot_binned(lM,dm.S,'dark Matter','k',nbins=5)
plot_binned(lM,gx.S,'galaxias','C2',nbins=5)
plt.legend()
plot_binned(lM[mold],gx_e.S[mold],'extendidas','C3',style=':',nbins=5)
plot_binned(lM[mold],gx_c.S[mold],'concentradas','C3',style='--',nbins=5)
plt.xlabel('$lM(R200)$')
plt.ylabel('$S_{gxs}(R200)$')
plt.axis([14,15.6,0.,1])
plt.savefig(path_plots+'SlM_concentradas_extendidas_new.png')

plt.figure()
plot_binned(lM[mnew],gx_c.T[mnew],'concentradas','C0',style='--',nbins=5)
plot_binned(lM[mnew],gx_e.T[mnew],'extendidas','C0',style=':',nbins=5)
plot_binned(lM[mnew],dm.T[mnew],'dark Matter','C0',nbins=5)
plot_binned(lM,dm.T,'dark Matter','k',nbins=5)
plot_binned(lM,gx.T,'galaxias','C2',nbins=5)
plt.legend()
plot_binned(lM[mold],gx_e.T[mold],'extendidas','C3',style=':',nbins=5)
plot_binned(lM[mold],gx_c.T[mold],'concentradas','C3',style='--',nbins=5)
plt.xlabel('$lM(R200)$')
plt.ylabel('$T_{gxs}(R200)$')
plt.axis([14,15.6,0.,1])
plt.savefig(path_plots+'TlM_concentradas_extendidas_new.png')

plt.figure()
plot_binned(lMp[mnew2D],gx_c.q[mnew2D],'concentradas','C0',style='--',nbins=5)
plot_binned(lMp[mnew2D],gx_e.q[mnew2D],'extendidas','C0',style=':',nbins=5)
plot_binned(lMp,dm.q,'dark Matter','k',nbins=5)
plot_binned(lMp,gx.q,'galaxias','C2',nbins=5)
plt.legend()
plot_binned(lMp[mold2D],gx_e.q[mold2D],'extendidas','C3',style=':',nbins=5)
plot_binned(lMp[mold2D],gx_c.q[mold2D],'concentradas','C3',style='--',nbins=5)
plt.xlabel('$lM(R200)$')
plt.ylabel('$q_{gxs}(R200)$')
plt.axis([14,15.6,0.,1])
plt.savefig(path_plots+'qlM_concentradas_extendidas_new.png')
