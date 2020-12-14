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
    
    
    
    dm = DarkMatter(radio)
    dm0 = DarkMatter(radio,False)
    
    Nh  = Galaxias(200).N
    nh  = Galaxias(200).n
    Nh_x = Galaxias(200,True).N
    nh_x = Galaxias(200,True).n
    
    
    pT = np.round(pearsonr(lM,dm.T)[0],2)
    pS = np.round(pearsonr(lM,dm.S)[0],2)
    pq = np.round(pearsonr(lMp,dm.q)[0],2)
    
    mT = np.round(np.mean(dm.T-dm0.T),2)
    mS = np.round(np.mean(dm.S-dm0.S),2)
    mq = np.round(np.mean(dm.q-dm0.q),2)
    
    sT = np.round(np.std(dm.T-dm0.T),2)
    sS = np.round(np.std(dm.S-dm0.S),2)
    sq = np.round(np.std(dm.q-dm0.q),2)
    
    
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

        path_plots = '../plots/dm_plots/R'+str(radio)+'/'
        
        os.system('mkdir '+path_plots)

    
        plt.figure()
        # plt.hist(dm.T-dm0.T,np.linspace(-0.2,0.2,25),histtype='step',color='C1')
        plt.hist((dm.T-dm0.T),np.linspace(-0.2,0.2,25),histtype='step',color='C0')
        plt.xlabel('$T - T0$')
        plt.ylabel('$N$')
        plt.savefig(path_plots+'difT_dist_'+str(radio)+'.png')
        
        plt.figure()
        # plt.hist(dm.S-dm0.S,np.linspace(-0.2,0.2,25),histtype='step',color='C1')
        plt.hist((dm.S-dm0.S),np.linspace(-0.2,0.2,25),histtype='step',color='C0')
        plt.xlabel('$S - S0$')
        plt.ylabel('$N$')
        plt.savefig(path_plots+'difS_dist_'+str(radio)+'.png')
        
        plt.figure()
        plt.hist((dm.q-dm0.q),np.linspace(-0.2,0.2,25),histtype='step',color='C0')
        plt.xlabel('$q - q0$')
        plt.ylabel('$N$')
        plt.savefig(path_plots+'difq_dist_'+str(radio)+'.png')
        
        
        plt.figure()
        plt.plot(Nh,dm.T-dm0.T,'.')
        plt.plot(Nh_x,dm.T-dm0.T,'.')
        plt.plot([1,Nh.max()+1],[0,0],'C7--')
        plt.ylabel('$T - T0$')
        plt.xlabel('$N_{GAL}$')
        plt.savefig(path_plots+'NdifT_'+str(radio)+'.png')
        
        plt.figure()
        plt.plot(Nh,dm.S-dm0.S,'.')
        plt.plot(Nh_x,dm.S-dm0.S,'.')
        plt.plot([1,Nh.max()+1],[0,0],'C7--')
        plt.ylabel('$S - S0$')
        plt.xlabel('$N_{GAL}$')
        plt.savefig(path_plots+'NdifS_'+str(radio)+'.png')
        
        plt.figure()
        plt.plot(nh,dm.q-dm0.q,'.')
        plt.plot(nh_x,dm.q-dm0.q,'.')
        plt.plot([1,nh.max()+1],[0,0],'C7--')
        plt.ylabel('$q - q0$')
        plt.xlabel('$N_{GAL}$')
        plt.savefig(path_plots+'Ndifq_'+str(radio)+'.png')
        
        
        plt.figure()
        plt.plot(lM,dm.S,'C7.')
        plot_binned(lM,dm.S,'Dark matter','k',nbins=5)
        plt.plot(lM,dm0.S,'C7x')
        plot_binned(lM,dm0.S,'Dark matter','C7',nbins=5)
        plt.ylim([0,1])
        plt.xlabel('$\log M_{'+lab+'}$')
        plt.ylabel('$S$')
        plt.savefig(path_plots+'S_DM_'+str(radio)+'_'+lab+'.png')
        
        plt.figure()
        plt.plot(lM,dm.T,'C7.')
        plot_binned(lM,dm.T,'Dark matter','k',nbins=5)
        plt.plot(lM,dm0.T,'C7x')
        plot_binned(lM,dm0.T,'Dark matter','C7',nbins=5)
        plt.ylim([0,1])
        plt.xlabel('$\log M_{'+lab+'}$')
        plt.ylabel('$T$')
        plt.savefig(path_plots+'T_DM_'+str(radio)+'_'+lab+'.png')
        
        plt.figure()
        plt.plot(lMp,dm.q,'C7.')
        plot_binned(lMp,dm.q,'Dark matter','k',nbins=5)
        plt.plot(lMp,dm0.q,'C7x')
        plot_binned(lMp,dm0.q,'Dark matter','C7',nbins=5)
        plt.ylim([0,1])
        plt.xlabel('$\log M_{'+lab+'}$')
        plt.ylabel('$q$')
        plt.savefig(path_plots+'q_DM_'+str(radio)+'_'+lab+'.png')
        
    return pT, pS, pq, mT, mS, mq, sT, sS, sq


options = [30,50,100,1000,500,200]
order   = np.arange(len(options))

lM = np.log10(gral[10])
lMp = np.array((lM.tolist())*3)

M = np.array([])
R = np.array([])

PT = np.array([])
PS = np.array([])
Pq = np.array([])

MT = np.array([])
MS = np.array([])
Mq = np.array([])

ST = np.array([])
SS = np.array([])
Sq = np.array([])

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
        pT, pS, pq, mT, mS, mq, sT, sS, sq = plot_dm(options[order[j]],options[order[i]])
    
        R = np.append(R,order[j])
        M = np.append(M,order[i])
    
        PT = np.append(PT,pT)
        PS = np.append(PS,pS)
        Pq = np.append(Pq,pq)

        MT = np.append(MT,mT)
        MS = np.append(MS,mS)
        Mq = np.append(Mq,mq)

        ST = np.append(ST,sT)
        SS = np.append(SS,sS)
        Sq = np.append(Sq,sq)


path_plots = '../plots/news/DM/'
        
plt.figure()
plt.scatter(R,PT,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_T')
plt.savefig(path_plots+'pT_DM_radio.png')

plt.figure()
plt.scatter(M,PT,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_T')
plt.savefig(path_plots+'pT_DM_mass.png')

plt.figure()
plt.scatter(R,PS,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_S')
plt.savefig(path_plots+'pS_DM_radio.png')

plt.figure()
plt.scatter(M,PS,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_S')
plt.savefig(path_plots+'pS_DM_mass.png')


plt.figure()
plt.scatter(R,Pq,c=M)
plt.xticks(order,['30kpc','50kpc','0.1R500','R1000','R500','R200'])
plt.xlabel('R')
plt.ylabel('p_q')
plt.savefig(path_plots+'pS_DM_radio.png')

plt.figure()
plt.scatter(M,Pq,c=R)
plt.xticks(order,['M30','M50','M0.1R500','M1000','M500','M200'])
plt.xlabel('M')
plt.ylabel('p_q')
plt.savefig(path_plots+'pq_DM_mass.png')


# WITH SUBHALOS

dm30   = DarkMatter(30)
dm50   = DarkMatter(50)
dm100  = DarkMatter(100)
dm200  = DarkMatter(200)
dm500  = DarkMatter(500)
dm1000 = DarkMatter(1000)

ct30_3D  , t30_3D     = cosangle(dm200.a3D,dm30.a3D) 
ct50_3D  , t50_3D     = cosangle(dm200.a3D,dm50.a3D)  
ct100_3D , t100_3D    = cosangle(dm200.a3D,dm100.a3D) 
ct200_3D , t200_3D    = cosangle(dm200.a3D,dm200.a3D) 
ct500_3D , t500_3D    = cosangle(dm200.a3D,dm500.a3D)
ct1000_3D, t1000_3D   = cosangle(dm200.a3D,dm1000.a3D)

ct30_2D  , t30_2D     = cosangle(dm200.a2D,dm30.a2D) 
ct50_2D  , t50_2D     = cosangle(dm200.a2D,dm50.a2D)  
ct100_2D , t100_2D    = cosangle(dm200.a2D,dm100.a2D) 
ct200_2D , t200_2D    = cosangle(dm200.a2D,dm200.a2D) 
ct500_2D , t500_2D    = cosangle(dm200.a2D,dm500.a2D)
ct1000_2D, t1000_2D   = cosangle(dm200.a2D,dm1000.a2D)

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

# WITHOUT SUBHALOS

dm30   = DarkMatter(30,False)
dm50   = DarkMatter(50,False)
dm100  = DarkMatter(100,False)
dm200  = DarkMatter(200,False)
dm500  = DarkMatter(500,False)
dm1000 = DarkMatter(1000,False)

ct30_3D  , t30_3D     = cosangle(dm200.a3D,dm30.a3D) 
ct50_3D  , t50_3D     = cosangle(dm200.a3D,dm50.a3D)  
ct100_3D , t100_3D    = cosangle(dm200.a3D,dm100.a3D) 
ct200_3D , t200_3D    = cosangle(dm200.a3D,dm200.a3D) 
ct500_3D , t500_3D    = cosangle(dm200.a3D,dm500.a3D)
ct1000_3D, t1000_3D   = cosangle(dm200.a3D,dm1000.a3D)

ct30_2D  , t30_2D     = cosangle(dm200.a2D,dm30.a2D) 
ct50_2D  , t50_2D     = cosangle(dm200.a2D,dm50.a2D)  
ct100_2D , t100_2D    = cosangle(dm200.a2D,dm100.a2D) 
ct200_2D , t200_2D    = cosangle(dm200.a2D,dm200.a2D) 
ct500_2D , t500_2D    = cosangle(dm200.a2D,dm500.a2D)
ct1000_2D, t1000_2D   = cosangle(dm200.a2D,dm1000.a2D)

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
plt.savefig(path_plots+'S_R'+indicator+'_woSH.png')

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
plt.savefig(path_plots+'T_R'+indicator+'_woSH.png')

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
plt.savefig(path_plots+'q_Rp'+indicator+'_woSH.png')

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
plt.savefig(path_plots+'t3D_R'+indicator+'_woSH.png')

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
plt.savefig(path_plots+'t2D_Rp'+indicator+'_woSH.png')
