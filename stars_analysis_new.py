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


def plot_dm(radio,mass,plot = False):
    
    
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
    
    
    

    dm = DarkMatter(200)
    dmr = DarkMatter(radio)
    
    stars = Stars(radio)
    
    pT = np.round(pearsonr(lM,stars.T)[0],2)
    pS = np.round(pearsonr(lM,stars.S)[0],2)
    pq = np.round(pearsonr(lMp,stars.q)[0],2)
    
    pT_dM = np.round(pearsonr(dm.T,stars.T)[0],2)
    pS_dM = np.round(pearsonr(dm.S,stars.S)[0],2)
    pq_dM = np.round(pearsonr(dm.q,stars.q)[0],2)

    pT_dMr = np.round(pearsonr(dmr.T,stars.T)[0],2)
    pS_dMr = np.round(pearsonr(dmr.S,stars.S)[0],2)
    pq_dMr = np.round(pearsonr(dmr.q,stars.q)[0],2)
    
    
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

        path_plots = '../plots/stars_plots/R'+str(radio)+'/'
        
        os.system('mkdir '+path_plots)


        plt.figure()
        plt.plot(lM,stars.S,'C7.')
        plot_binned(lM,dm.S,'Dark matter','k',nbins=5)
        plt.plot(lM,dm.S,'C0x')
        plot_binned(lM,stars.S,'stars','C0',nbins=5)
        plt.ylim([0,1])
        plt.xlabel('$\log M_{'+lab+'}$')
        plt.ylabel('$S$')
        plt.savefig(path_plots+'S_stars_'+str(radio)+'_'+lab+'.png')
        
        plt.figure()
        plt.plot(lM,stars.T,'C7.')
        plot_binned(lM,dm.T,'Dark matter','k',nbins=5)
        plt.plot(lM,dm.T,'C0x')
        plot_binned(lM,stars.T,'stars','C0',nbins=5)
        plt.ylim([0,1])
        plt.xlabel('$\log M_{'+lab+'}$')
        plt.ylabel('$T$')
        plt.savefig(path_plots+'T_stars_'+str(radio)+'_'+lab+'.png')
        
        plt.figure()
        plt.plot(lMp,stars.q,'C7.')
        plot_binned(lMp,dm.q,'Dark matter','k',nbins=5)
        plt.plot(lMp,dm.q,'C0x')
        plot_binned(lMp,stars.q,'stars','C0',nbins=5)
        plt.ylim([0,1])
        plt.xlabel('$\log M_{'+lab+'}$')
        plt.ylabel('$q$')
        plt.savefig(path_plots+'q_stars_'+str(radio)+'_'+lab+'.png')
        
    return pT, pS, pq, pT_dM, pS_dM, pq_dM, pT_dMr, pS_dMr, pq_dMr


options = [30,50,100,1000,500,200]
order   = np.arange(len(options))

lM = np.log10(gral[9])
lMp = np.array((lM.tolist())*3)

M = np.array([])
R = np.array([])

PT = np.array([])
PS = np.array([])
Pq = np.array([])

PT_dM = np.array([])
PS_dM = np.array([])
Pq_dM = np.array([])

PT_dMr = np.array([])
PS_dMr = np.array([])
Pq_dMr = np.array([])


indicator = 'DV'

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

for j in order:
    for i in order:
        pT, pS, pq, pT_dM, pS_dM, pq_dM, pT_dMr, pS_dMr, pq_dMr = plot_dm(options[order[j]],options[order[i]],True)
    
        R = np.append(R,order[j])
        M = np.append(M,order[i])
    
        PT = np.append(PT,pT)
        PS = np.append(PS,pS)
        Pq = np.append(Pq,pq)

        PT_dM = np.append(PT_dM,pT_dM)
        PS_dM = np.append(PS_dM,pS_dM)
        Pq_dM = np.append(Pq_dM,pq_dM)

        PT_dMr = np.append(PT_dMr,pT_dMr)
        PS_dMr = np.append(PS_dMr,pS_dMr)
        Pq_dMr = np.append(Pq_dMr,pq_dMr)


path_plots = '../plots/news/stars/'

# CORRELATION WITH MASS
        
plt.figure()
plt.scatter(R,PT,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_T')
plt.savefig(path_plots+'pT_mass_radio.png')

plt.figure()
plt.scatter(M,PT,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_T')
plt.savefig(path_plots+'pT_mass_mass.png')

plt.figure()
plt.scatter(R,PS,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_S')
plt.savefig(path_plots+'pS_mass_radio.png')

plt.figure()
plt.scatter(M,PS,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_S')
plt.savefig(path_plots+'pS_mass_mass.png')


plt.figure()
plt.scatter(R,Pq,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_q')
plt.savefig(path_plots+'pS_mass_radio.png')

plt.figure()
plt.scatter(M,Pq,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_q')
plt.savefig(path_plots+'pq_mass_mass.png')


# CORRELATION WITH DM
        
plt.figure()
plt.scatter(R,PT_dMr,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_T')
plt.savefig(path_plots+'pT_dM_radio.png')

plt.figure()
plt.scatter(M,PT_dMr,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_T')
plt.savefig(path_plots+'pT_dM_mass.png')

plt.figure()
plt.scatter(R,PS_dMr,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_S')
plt.savefig(path_plots+'pS_dM_radio.png')

plt.figure()
plt.scatter(M,PS_dMr,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_S')
plt.savefig(path_plots+'pS_dM_mass.png')

plt.figure()
plt.scatter(R,Pq_dMr,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_q')
plt.savefig(path_plots+'pS_dM_radio.png')

plt.figure()
plt.scatter(M,Pq_dMr,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_q')
plt.savefig(path_plots+'pq_dM_mass.png')

# CORRELATION WITH DM200
        
plt.figure()
plt.scatter(R,PT_dM,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_T')
plt.savefig(path_plots+'pT_dM_radio.png')

plt.figure()
plt.scatter(M,PT_dMr,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_T')
plt.savefig(path_plots+'pT_dM_mass.png')

plt.figure()
plt.scatter(R,PS_dMr,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_S')
plt.savefig(path_plots+'pS_dM_radio.png')

plt.figure()
plt.scatter(M,PS_dMr,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_S')
plt.savefig(path_plots+'pS_dM_mass.png')

plt.figure()
plt.scatter(R,Pq_dMr,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_q')
plt.savefig(path_plots+'pS_dM_radio.png')

plt.figure()
plt.scatter(M,Pq_dMr,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_q')
plt.savefig(path_plots+'pq_dM_mass.png')


# WITH RADIUS

dm200  = DarkMatter(200)

s30   = Stars(30)
s50   = Stars(50)
s100  = Stars(100)
s500  = Stars(500)
s200  = Stars(200)
s1000 = Stars(1000)

ct30_3D  , t30_3D     = cosangle(dm200.a3D,s30.a3D) 
ct50_3D  , t50_3D     = cosangle(dm200.a3D,s50.a3D)  
ct100_3D , t100_3D    = cosangle(dm200.a3D,s100.a3D) 
ct200_3D , t200_3D    = cosangle(dm200.a3D,s200.a3D) 
ct500_3D , t500_3D    = cosangle(dm200.a3D,s500.a3D)
ct1000_3D, t1000_3D   = cosangle(dm200.a3D,s1000.a3D)

ct30_2D  , t30_2D     = cosangle(dm200.a2D,s30.a2D) 
ct50_2D  , t50_2D     = cosangle(dm200.a2D,s50.a2D)  
ct100_2D , t100_2D    = cosangle(dm200.a2D,s100.a2D) 
ct200_2D , t200_2D    = cosangle(dm200.a2D,s200.a2D) 
ct500_2D , t500_2D    = cosangle(dm200.a2D,s500.a2D)
ct1000_2D, t1000_2D   = cosangle(dm200.a2D,s1000.a2D)

t3D = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
t2D = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T

S = np.vstack((dm30.S,dm50.S,dm100.S,dm1000.S,dm500.S,dm200.S)).T
T = np.vstack((dm30.T,dm50.T,dm100.T,dm1000.T,dm500.T,dm200.T)).T
q = np.vstack((dm30.q,dm50.q,dm100.q,dm1000.q,dm500.q,dm200.q)).T


# ANGLES AND RADIUS

R1000 = gral[4]
R500  = gral[5]
R200  = gral[6]
R30   = np.ones(len(R200))*30.
R50   = np.ones(len(R200))*50.
R1    = 0.1*R500
R = np.vstack((R30,R50,R1,R1000,R500,R200)).T
Rp = np.array((R.tolist())*3)


plt.figure()    
plt.plot(np.median(np.log10(R),axis=0),np.median(S,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(S,axis=0)+np.std(S,axis=0),np.median(S,axis=0)-np.std(S,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mnew]),axis=0),np.median(S[mnew],axis=0),'C0')
plt.fill_between(np.median(np.log10(R[mnew]),axis=0),np.median(S[mnew],axis=0)+np.std(S[mnew],axis=0),np.median(S[mnew],axis=0)-np.std(S[mnew],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(R[mold]),axis=0),np.median(S[mold],axis=0),'C3')
plt.fill_between(np.median(np.log10(R[mold]),axis=0),np.median(S[mold],axis=0)+np.std(S[mold],axis=0),np.median(S[mold],axis=0)-np.std(S[mold],axis=0),color = 'C3',alpha=0.1)
plt.ylabel('$S$')
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'S_R'+indicator+'.png')

plt.figure()
plt.plot(np.median(np.log10(R),axis=0),np.median(T,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(T,axis=0)+np.std(T,axis=0),np.median(T,axis=0)-np.std(T,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mnew]),axis=0),np.median(T[mnew],axis=0),'C0')
plt.fill_between(np.median(np.log10(R[mnew]),axis=0),np.median(T[mnew],axis=0)+np.std(T[mnew],axis=0),np.median(T[mnew],axis=0)-np.std(T[mnew],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(R[mold]),axis=0),np.median(T[mold],axis=0),'C3')
plt.fill_between(np.median(np.log10(R[mold]),axis=0),np.median(T[mold],axis=0)+np.std(T[mold],axis=0),np.median(T[mold],axis=0)-np.std(T[mold],axis=0),color = 'C3',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel('$T$')
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'T_R'+indicator+'.png')

plt.figure()
plt.plot(np.median(np.log10(Rp),axis=0),np.median(q,axis=0),'k')
plt.fill_between(np.median(np.log10(Rp),axis=0),np.median(q,axis=0)+np.std(q,axis=0),np.median(q,axis=0)-np.std(q,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mnew2D]),axis=0),np.median(q[mnew2D],axis=0),'C0')
plt.fill_between(np.median(np.log10(Rp[mnew2D]),axis=0),np.median(q[mnew2D],axis=0)+np.std(q[mnew2D],axis=0),np.median(q[mnew2D],axis=0)-np.std(q[mnew2D],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mold2D]),axis=0),np.median(q[mold2D],axis=0),'C3')
plt.fill_between(np.median(np.log10(Rp[mold2D]),axis=0),np.median(q[mold2D],axis=0)+np.std(q[mold2D],axis=0),np.median(q[mold2D],axis=0)-np.std(q[mold2D],axis=0),color = 'C3',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel('$q$')
plt.axis([np.log10(30),np.log10(2000),0.4,0.9])
plt.savefig(path_plots+'q_Rp'+indicator+'.png')

plt.figure()
plt.plot(np.median(np.log10(R),axis=0),np.median(t3D,axis=0),'k')
plt.fill_between(np.median(np.log10(R),axis=0),np.median(t3D,axis=0)+np.std(t3D,axis=0),np.median(t3D,axis=0)-np.std(t3D,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(R[mnew]),axis=0),np.median(t3D[mnew],axis=0),'C0')
plt.fill_between(np.median(np.log10(R[mnew]),axis=0),np.median(t3D[mnew],axis=0)+np.std(t3D[mnew],axis=0),np.median(t3D[mnew],axis=0)-np.std(t3D[mnew],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(R[mold]),axis=0),np.median(t3D[mold],axis=0),'C3')
plt.fill_between(np.median(np.log10(R[mold]),axis=0),np.median(t3D[mold],axis=0)+np.std(t3D[mold],axis=0),np.median(t3D[mold],axis=0)-np.std(t3D[mold],axis=0),color = 'C3',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel(r'$\theta$')
plt.axis([np.log10(30),np.log10(2000),0.,50.])
plt.savefig(path_plots+'t3D_R'+indicator+'.png')

plt.figure()
plt.plot(np.median(np.log10(Rp),axis=0),np.median(t2D,axis=0),'k')
plt.fill_between(np.median(np.log10(Rp),axis=0),np.median(t2D,axis=0)+np.std(t2D,axis=0),np.median(t2D,axis=0)-np.std(t2D,axis=0),color = 'k',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mnew2D]),axis=0),np.median(t2D[mnew2D],axis=0),'C0')
plt.fill_between(np.median(np.log10(Rp[mnew2D]),axis=0),np.median(t2D[mnew2D],axis=0)+np.std(t2D[mnew2D],axis=0),np.median(t2D[mnew2D],axis=0)-np.std(t2D[mnew2D],axis=0),color = 'C0',alpha=0.1)
plt.plot(np.median(np.log10(Rp[mold2D]),axis=0),np.median(t2D[mold2D],axis=0),'C3')
plt.fill_between(np.median(np.log10(Rp[mold2D]),axis=0),np.median(t2D[mold2D],axis=0)+np.std(t2D[mold2D],axis=0),np.median(t2D[mold2D],axis=0)-np.std(t2D[mold2D],axis=0),color = 'C3',alpha=0.1)
plt.xticks(np.median(np.log10(R),axis=0),['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.ylabel(r'$\theta_{2D}$')
plt.axis([np.log10(30),np.log10(2000),0.,50.])
plt.savefig(path_plots+'t2D_Rp'+indicator+'.png')

