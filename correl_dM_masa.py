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


def plot_dm(radio,mass,mask,mask2D,plot = False,sample = 'all'):
    
    
    if mass == 1000:
        lM = np.log10(gral[7])
        lab = '1000'
    elif mass == 500:
        lM = np.log10(gral[8])
        lab = '500'
    elif mass == 200:
        lM = np.log10(gral[9])
        lab = '200'
    elif mass == 30:
        lM = np.log10(gral[10])
        lab = 'BCG-30kpc'
    elif mass == 50:
        lM = np.log10(gral[11])
        lab = 'BCG-50kpc'
    elif mass == 100:
        lM = np.log10(gral[12])
        lab = 'BCG-0.1R500'
    
    lMp = np.array((lM.tolist())*3)
    
    lM = lM[mask]
    lMp = lMp[mask2D]
    
    
    dm = DarkMatter(radio)
    dm0 = DarkMatter(radio,False)

    
    pT = np.round(pearsonr(lM,dm.T[mask])[0],3)
    pS = np.round(pearsonr(lM,dm.S[mask])[0],3)
    pq = np.round(pearsonr(lMp,dm.q[mask2D])[0],3)
    
    
    
    # print('DM within R'+str(radio)+' M'+lab)
    # print('Difference with non subhalos')
    # print('T diff: mean = '+np.str(mT)+' std = '+np.str(sT))
    # print('S diff: mean = '+np.str(mS)+' std = '+np.str(sS))
    # print('q diff: mean = '+np.str(mq)+' std = '+np.str(sq))
    # print('Correlation with mass')
    # print('T : p = '+np.str(pT))
    # print('S : p = '+np.str(pS))
    # print('q : p = '+np.str(pq))
    
    if plot: 

        path_plots = '../plots/correl_dM_masa/'
        
       
        plt.figure()
        plt.title(sample)
        plt.plot(lM,dm.S[mask],'C7.')
        plot_binned(lM,dm.S[mask],'Dark matter','k',nbins=5)
        plt.plot(lM,dm0.S[mask],'C7x')
        plot_binned(lM,dm0.S[mask],'Dark matter','C7',nbins=5)
        plt.ylim([0,1])
        plt.xlabel('$\log M_{'+lab+'}$')
        plt.ylabel('$S(R'+str(radio)+'$')
        plt.savefig(path_plots+'S_DM_R'+str(radio)+'_M'+lab+'_'+sample+'.png')
        
        plt.figure()
        plt.title(sample)
        plt.plot(lM,dm.T[mask],'C7.')
        plot_binned(lM,dm.T[mask],'Dark matter','k',nbins=5)
        plt.plot(lM,dm0.T[mask],'C7x')
        plot_binned(lM,dm0.T[mask],'Dark matter','C7',nbins=5)
        plt.ylim([0,1])
        plt.xlabel('$\log M_{'+lab+'}$')
        plt.ylabel('$T(R'+str(radio)+'$')
        plt.savefig(path_plots+'T_DM_R'+str(radio)+'_M'+lab+'_'+sample+'.png')
        
        plt.figure()
        plt.title(sample)
        plt.plot(lMp,dm.q[mask2D],'C7.')
        plot_binned(lMp,dm.q[mask2D],'Dark matter','k',nbins=5)
        plt.plot(lMp,dm0.q[mask2D],'C7x')
        plot_binned(lMp,dm0.q[mask2D],'Dark matter','C7',nbins=5)
        plt.ylim([0,1])
        plt.xlabel('$\log M_{'+lab+'}$')
        plt.ylabel('$q(R'+str(radio)+'$')
        plt.savefig(path_plots+'q_DM_R'+str(radio)+'_M'+lab+'_'+sample+'.png')
        
    return pT, pS, pq


options = [30,50,100,1000,500,200]
order   = np.arange(len(options))
Rlabels   = ['30kpc','50kpc','0.1R500','R1000','R500','R200']
Mlabels   = ['M(30kpc)','M(50kpc)','M(0.1R500)','M(R1000)','M(R500)','M(R200)']


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
    
    M = np.array([])
    R   = np.array([])
    
    PT = np.array([])
    PS = np.array([])
    Pq = np.array([])

    for j in order:
        for i in order:
            pT, pS, pq = plot_dm(options[order[j]],options[order[i]],mask[m],mask2D[m],sample=samples[m])
        
            R = np.append(R,order[j])
            M = np.append(M,order[i])
        
            PT = np.append(PT,pT)
            PS = np.append(PS,pS)
            Pq = np.append(Pq,pq)
    
    
    
    path_plots = '../plots/correl_dM_masa/'
            
    plt.figure()
    for x in order:
        plt.plot(M[R==x],PT[R==x],label=Rlabels[x],lw = (x+3)/4.)
    plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
    plt.xlabel('M')
    plt.ylabel('p_T')
    plt.axis([-1,6,-0.7,0.5])
    plt.legend(loc=3,ncol=4)
    plt.savefig(path_plots+'pT_DM_mass_'+sample+'.png')
    
    plt.figure()
    for x in order:
        plt.plot(M[R==x],PS[R==x],label=Rlabels[x],lw = (x+3)/4.)
    plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
    plt.xlabel('M')
    plt.ylabel('p_S')
    plt.axis([-1,6,-0.7,0.5])
    plt.legend(loc=3,ncol=4)
    plt.savefig(path_plots+'pS_DM_mass_'+sample+'.png')
    
    plt.figure()
    for x in order:
        plt.plot(M[R==x],Pq[R==x],label=Rlabels[x],lw = (x+3)/4.)
    plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
    plt.xlabel('M')
    plt.ylabel('p_q')
    plt.axis([-1,6,-0.7,0.5])
    plt.legend(loc=3,ncol=4)
    plt.savefig(path_plots+'pq_DM_mass_'+sample+'.png')


lM200 = np.log10(gral[9])
lM30 = np.log10(gral[10])

lM200p = np.array((lM200.tolist())*3)
lM30p = np.array((lM30.tolist())*3)


dm200 = DarkMatter(200)
dm1000 = DarkMatter(1000)
dm30 = DarkMatter(30)

mnew,mold,mnew2D,mold2D = newold()

plt.figure()
plot_binned(lM30[mnew],dm30.S[mnew],'R < 30kpc','C0',style='',nbins=5)
plot_binned(lM30[mnew],dm1000.S[mnew],'R < R1000','C0',style='--',nbins=5)
plot_binned(lM30[mnew],dm200.S[mnew],'R < R200','C0',style=':',nbins=5)
plt.legend()
plot_binned(lM30[mold],dm30.S[mold],'dark Matter','C3',nbins=5)
plot_binned(lM30[mold],dm1000.S[mold],'R < R1000','C3',style='--',nbins=5)
plot_binned(lM30[mold],dm200.S[mold],'R < R200','C3',style=':',nbins=5)
plt.xlabel('$lM(30kpc)$')
plt.ylabel('$S_{DM}(R)$')
plt.axis([11.25,11.8,0.5,0.9])
plt.savefig(path_plots+'SlM30_dm.png')

plt.figure()
plot_binned(lM30[mnew],dm30.T[mnew],'R < 30kpc','C0',style='',nbins=5)
plot_binned(lM30[mnew],dm1000.T[mnew],'R < R1000','C0',style='--',nbins=5)
plot_binned(lM30[mnew],dm200.T[mnew],'R < R200','C0',style=':',nbins=5)
plt.legend()
plot_binned(lM30[mold],dm30.T[mold],'dark Matter','C3',nbins=5)
plot_binned(lM30[mold],dm1000.T[mold],'R < R1000','C3',style='--',nbins=5)
plot_binned(lM30[mold],dm200.T[mold],'R < R200','C3',style=':',nbins=5)
plt.xlabel('$lM(30kpc)$')
plt.ylabel('$T_{DM}(R)$')
plt.axis([11.25,11.8,0.,1.])
plt.savefig(path_plots+'TlM30_dm.png')

plt.figure()
plot_binned(lM30[mnew],dm30.T[mnew],'new','C0',style='',nbins=5)
plot_binned(lM30[mold],dm30.T[mold],'old','C3',nbins=5)
plt.plot(lM30[mnew],dm30.T[mnew],'C0.')
plt.plot(lM30[mold],dm30.T[mold],'C3.')
plt.legend()
plt.xlabel('$lM(30kpc)$')
plt.ylabel('$T_{DM}(R < 30kpc)$')
plt.axis([11.25,11.8,0.,1.])
plt.savefig(path_plots+'TlM30_dm30.png')

plt.figure()
plot_binned(lM30[mnew],dm1000.T[mnew],'new','C0',style='--',nbins=5)
plot_binned(lM30[mold],dm1000.T[mold],'old','C3',style='--',nbins=5)
plt.plot(lM30[mnew],dm1000.T[mnew],'C0.')
plt.plot(lM30[mold],dm1000.T[mold],'C3.')
plt.legend()
plt.xlabel('$lM(30kpc)$')
plt.ylabel('$T_{DM}(R < R1000)$')
plt.axis([11.25,11.8,0.,1.])
plt.savefig(path_plots+'TlM30_dm1000.png')

plt.figure()
plot_binned(lM30[mnew],dm200.T[mnew],'new','C0',style=':',nbins=5)
plot_binned(lM30[mold],dm200.T[mold],'old','C3',style=':',nbins=5)
plt.plot(lM30[mnew],dm200.T[mnew],'C0.')
plt.plot(lM30[mold],dm200.T[mold],'C3.')
plt.legend()
plt.xlabel('$lM(30kpc)$')
plt.ylabel('$T_{DM}(R < R200)$')
plt.axis([11.25,11.8,0.,1.])
plt.savefig(path_plots+'TlM30_dm200.png')

plt.figure()
plot_binned(lM200[mnew],dm30.T[mnew],'new','C0',style='',nbins=5)
plot_binned(lM200[mold],dm30.T[mold],'old','C3',nbins=5)
plt.plot(lM200[mnew],dm30.T[mnew],'C0.')
plt.plot(lM200[mold],dm30.T[mold],'C3.')
plt.legend()
plt.xlabel('$lM(R200)$')
plt.ylabel('$T_{DM}(R < 30kpc)$')
plt.axis([14.,15.6,0.,1.0])
plt.savefig(path_plots+'TlM200_dm30.png')

plt.figure()
plot_binned(lM200[mnew],dm1000.T[mnew],'new','C0',style='--',nbins=5)
plot_binned(lM200[mold],dm1000.T[mold],'old','C3',style='--',nbins=5)
plt.plot(lM200[mnew],dm1000.T[mnew],'C0.')
plt.plot(lM200[mold],dm1000.T[mold],'C3.')
plt.legend()
plt.xlabel('$lM(R200)$')
plt.ylabel('$T_{DM}(R < R1000)$')
plt.axis([14.,15.6,0.,1.0])
plt.savefig(path_plots+'TlM200_dm1000.png')

plt.figure()
plot_binned(lM200,dm200.T,'all','k',style=':',nbins=3)
plot_binned(lM200[mnew],dm200.T[mnew],'new','C0',style=':',nbins=3)
plot_binned(lM200[mold],dm200.T[mold],'old','C3',style=':',nbins=3)
# plt.plot(lM200[mnew],dm200.T[mnew],'C0.')
# plt.plot(lM200[mold],dm200.T[mold],'C3.')
plt.legend()
plt.xlabel('$lM(R200)$')
plt.ylabel('$T_{DM}(R < R200)$')
plt.axis([14.,15.6,0.,1.0])
plt.savefig(path_plots+'TlM200_dm200.png')


plt.figure()
plot_binned(lM30p[mnew2D],dm30.q[mnew2D],'R < 30kpc','C0',style='',nbins=5)
plot_binned(lM30p[mnew2D],dm1000.q[mnew2D],'R < R1000','C0',style='--',nbins=5)
plot_binned(lM30p[mnew2D],dm200.q[mnew2D],'R < R200','C0',style=':',nbins=5)
plt.legend()
plot_binned(lM30p[mold2D],dm30.q[mold2D],'dark Matter','C3',nbins=5)
plot_binned(lM30p[mold2D],dm1000.q[mold2D],'R < R1000','C3',style='--',nbins=5)
plot_binned(lM30p[mold2D],dm200.q[mold2D],'R < R200','C3',style=':',nbins=5)
plt.xlabel('$lM(30kpc)$')
plt.ylabel('$q_{DM}(R)$')
plt.axis([11.25,11.8,0.6,1.0])
plt.savefig(path_plots+'qlM30p_dm.png')


plt.figure()
plot_binned(lM200[mnew],dm30.S[mnew],'R < 30kpc','C0',style='',nbins=5)
plot_binned(lM200[mnew],dm1000.S[mnew],'R < R1000','C0',style='--',nbins=5)
plot_binned(lM200[mnew],dm200.S[mnew],'R < R200','C0',style=':',nbins=5)
plt.legend()
plot_binned(lM200[mold],dm30.S[mold],'dark Matter','C3',nbins=5)
plot_binned(lM200[mold],dm1000.S[mold],'R < R1000','C3',style='--',nbins=5)
plot_binned(lM200[mold],dm200.S[mold],'R < R200','C3',style=':',nbins=5)
plt.xlabel('$lM(R200)$')
plt.ylabel('$S_{DM}(R)$')
plt.axis([14.,15.6,0.5,0.9])
plt.savefig(path_plots+'SlM200_dm.png')

plt.figure()
plot_binned(lM200[mnew],dm30.T[mnew],'R < 30kpc','C0',style='',nbins=5)
plot_binned(lM200[mnew],dm1000.T[mnew],'R < R1000','C0',style='--',nbins=5)
plot_binned(lM200[mnew],dm200.T[mnew],'R < R200','C0',style=':',nbins=5)
plt.legend()
plot_binned(lM200[mold],dm30.T[mold],'dark Matter','C3',nbins=5)
plot_binned(lM200[mold],dm1000.T[mold],'R < R1000','C3',style='--',nbins=5)
plot_binned(lM200[mold],dm200.T[mold],'R < R200','C3',style=':',nbins=5)
plt.xlabel('$lM(R200)$')
plt.ylabel('$T_{DM}(R)$')
plt.axis([14.,15.6,0.,1.0])
plt.savefig(path_plots+'TlM200_dm.png')

plt.figure()
plot_binned(lM200p[mnew2D],dm30.q[mnew2D],'R < 30kpc','C0',style='',nbins=5)
plot_binned(lM200p[mnew2D],dm1000.q[mnew2D],'R < R1000','C0',style='--',nbins=5)
plot_binned(lM200p[mnew2D],dm200.q[mnew2D],'R < R200','C0',style=':',nbins=5)
plt.legend()
plot_binned(lM200p[mold2D],dm30.q[mold2D],'dark Matter','C3',nbins=5)
plot_binned(lM200p[mold2D],dm1000.q[mold2D],'R < R1000','C3',style='--',nbins=5)
plot_binned(lM200p[mold2D],dm200.q[mold2D],'R < R200','C3',style=':',nbins=5)
plt.xlabel('$lM(R200)$')
plt.ylabel('$q_{DM}(R)$')
plt.axis([14.,15.6,0.6,1.0])
plt.savefig(path_plots+'qlM200p_dm.png')
