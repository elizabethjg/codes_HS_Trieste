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


def plot_fig(x,y,nbins,ax=plt,color = 'C3', style = '',label=''):
                
        X,q50,q25,q75,mz = binned(x,y,nbins)
        ax.plot(X,q50,color+style,label=label)
        ax.plot(X,q75,color+style,alpha=0.2)
        ax.plot(X,q25,color+style,alpha=0.2)
        ax.fill_between(X,q75,q25,color = color,alpha=0.1)


plotspath = '../plots_ICL/'

C = Clusters()
ICL_ = ICL()

DM1000 = DarkMatter(1000)
DM500  = DarkMatter(500)
DM200  = DarkMatter(200)

gx200 = Galaxias(radio=200)
gx1000 = Galaxias(radio=1000)
gx500 = Galaxias(radio=500)

st200  = Stars(200)
st1000 = Stars(1000)
st500  = Stars(500)
st30   = Stars(30)
st50   = Stars(50)
st100  = Stars(100)

gx_c = Galaxias(tipo='concentradas')
gx_e = Galaxias(tipo='extendidas')

mall = np.ones(72).astype(bool)
mall2D = np.ones(72*3).astype(bool)

m = (C.sub == 0)

mp = (m.tolist()*3)

yaxis = np.zeros((len(DM200.a2D),2))
yaxis[:,1] = 1.

PA200 = cosangle2(DM200.a2D[mp],yaxis[mp])[1]
PA500 = cosangle2(DM500.a2D[mp],yaxis[mp])[1]
PA1000 = cosangle2(DM1000.a2D[mp],yaxis[mp])[1]

PA200st = cosangle2(st200.a2D[mp],yaxis[mp])[1]
PA500st = cosangle2(st500.a2D[mp],yaxis[mp])[1]
PA1000st = cosangle2(st1000.a2D[mp],yaxis[mp])[1]

r200  = C.Rp[mp,-1]
r500  = C.Rp[mp,-2]
r1000 = C.Rp[mp,-3]

mo = C.mo2d_gap[mp]


p    = ['xy','xz','yz']



for d in np.arange(29):

    plt.figure()

    D = [d,d+29,d+58]
    
    qdm = DM200.q[mp][D]

    if not all(np.diff(ICL_.D[D]) == 0):
        print('something wrong')

    for i in range(3): 
    
        j = D[i]
    
        
        dif1000 = abs(ICL_.PA - PA1000)[j]
        dif1000[dif1000 > 90.] = 180. - dif1000[dif1000 > 90.]    
    
        dif500 = abs(ICL_.PA - PA500)[j]
        dif500[dif500 > 90.] = 180. - dif500[dif500 > 90.]    
    
        dif200 = abs(ICL_.PA - PA200)[j]
        dif200[dif200 > 90.] = 180. - dif200[dif200 > 90.]    
        
        plt.plot(ICL_.a[j]/r200[j],dif1000,'C'+str(i+3)+'-',label = p[i]+' '+str(qdm[i]))
        plt.plot(ICL_.a[j]/r200[j],dif500,'C'+str(i+3)+'--')
        plt.plot(ICL_.a[j]/r200[j],dif200,'C'+str(i+3)+':')

        plt.legend(loc=1)

        # if mo[j]:
            # plt.plot(ICL_.a[j]/r200[j],dif1000,'C3o')
            # plt.plot(ICL_.a[j]/r200[j],dif500,'C'+str(i+5)+'o')
            # plt.plot(ICL_.a[j]/r200[j],dif200,'C'+str(i+5)+'o')

    
    
    plt.xlabel('$R/R_{200}$')
    plt.ylabel(r'$\theta_{DM}$')
    
    plt.axis([0,0.8,0,90])
    
    plt.savefig(plotspath+'D'+str(int(ICL_.D[d]))+'.png')



plt.figure()
for j in np.arange(87): 
    dif = abs(ICL_.PA - PA1000)[j]
    dif[dif > 90.] = 180. - dif[dif > 90.]    
    if mo[j]:
        plt.plot(ICL_.a[j]/r200[j],dif,'r',alpha = 0.5)
    else:
        plt.plot(ICL_.a[j]/r200[j],dif,'C0',alpha = 0.5)

plt.xlabel('$R/R_{200}$')
plt.ylabel(r'$\theta_{R1000}$')
plt.savefig(plotspath+'align_R1000.png')


plt.figure()
for j in np.arange(87): 
    dif = abs(ICL_.PA - PA500)[j]
    dif[dif > 90.] = 180. - dif[dif > 90.]    
    if mo[j]:
        plt.plot(ICL_.a[j]/r200[j],dif,'r',alpha = 0.5)
    else:
        plt.plot(ICL_.a[j]/r200[j],dif,'C0',alpha = 0.5)

plt.xlabel('$R/R_{200}$')
plt.ylabel(r'$\theta_{R500}$')
plt.savefig(plotspath+'align_R500.png')

plt.figure()
for j in np.arange(87): 
    dif = abs(ICL_.PA - PA200)[j]
    dif[dif > 90.] = 180. - dif[dif > 90.]    
    if mo[j]:
        plt.plot(ICL_.a[j]/r200[j],dif,'r',alpha = 0.5)
    else:
        plt.plot(ICL_.a[j]/r200[j],dif,'C0',alpha = 0.5)

plt.xlabel('$R/R_{200}$')
plt.ylabel(r'$\theta_{R200}$')
plt.savefig(plotspath+'align_R200.png')


rs_o    = np.array([])
t200_o  = np.array([])
t500_o  = np.array([])
t1000_o = np.array([])

rs_n    = np.array([])
t200_n  = np.array([])
t500_n  = np.array([])
t1000_n = np.array([])

for j in np.arange(87):
    
    dif1000 = abs(ICL_.PA - PA1000)[j]
    dif1000[dif1000 > 90.] = 180. - dif1000[dif1000 > 90.]    

    dif500 = abs(ICL_.PA - PA500)[j]
    dif500[dif500 > 90.] = 180. - dif500[dif500 > 90.]    

    dif200 = abs(ICL_.PA - PA200)[j]
    dif200[dif200 > 90.] = 180. - dif200[dif200 > 90.]    
    
    
    if mo[j]:
        rs_o = np.append(rs_o,(ICL_.a/r200)[j])
        t200_o = np.append(t200_o,dif200)
        t500_o = np.append(t500_o,dif500)
        t1000_o = np.append(t1000_o,dif1000)
    else:
        rs_n = np.append(rs_n,(ICL_.a/r200)[j])
        t200_n = np.append(t200_n,dif200)
        t500_n = np.append(t500_n,dif500)
        t1000_n = np.append(t1000_n,dif1000)


f, ax = plt.subplots(1,2, figsize=(12,5), sharex=True,sharey=True)
f.subplots_adjust(hspace=0,wspace=0)
    
plot_fig(np.log10(rs_o),t200_o,6,ax=ax[0],style=':')
plot_fig(np.log10(rs_o),t500_o,6,ax=ax[0],style='--')
plot_fig(np.log10(rs_o),t1000_o,6,ax=ax[0],style='-')

plot_fig(np.log10(rs_n),t200_n,6,color='C0',ax=ax[1],style=':')
plot_fig(np.log10(rs_n),t500_n,6,color='C0',ax=ax[1],style='--')
plot_fig(np.log10(rs_n),t1000_n,6,color='C0',ax=ax[1],style='-')

ax[0].set_xlabel('$\log(R/R200)$')
ax[1].set_xlabel('$\log(R/R200)$')
ax[0].set_ylabel(r'$\theta_{DM}$')

plt.savefig(plotspath+'align_medians_gap.png')


rs_o    = np.array([])
t200_o  = np.array([])
t500_o  = np.array([])
t1000_o = np.array([])

rs_n    = np.array([])
t200_n  = np.array([])
t500_n  = np.array([])
t1000_n = np.array([])

for j in np.arange(87):
    
    dif1000 = abs(ICL_.PA - PA1000st)[j]
    dif1000[dif1000 > 90.] = 180. - dif1000[dif1000 > 90.]    

    dif500 = abs(ICL_.PA - PA500st)[j]
    dif500[dif500 > 90.] = 180. - dif500[dif500 > 90.]    

    dif200 = abs(ICL_.PA - PA200st)[j]
    dif200[dif200 > 90.] = 180. - dif200[dif200 > 90.]    
    
    
    if mo[j]:
        rs_o = np.append(rs_o,(ICL_.a/r200)[j])
        t200_o = np.append(t200_o,dif200)
        t500_o = np.append(t500_o,dif500)
        t1000_o = np.append(t1000_o,dif1000)
    else:
        rs_n = np.append(rs_n,(ICL_.a/r200)[j])
        t200_n = np.append(t200_n,dif200)
        t500_n = np.append(t500_n,dif500)
        t1000_n = np.append(t1000_n,dif1000)


f, ax = plt.subplots(1,2, figsize=(12,5), sharex=True,sharey=True)
f.subplots_adjust(hspace=0,wspace=0)
    
plot_fig(np.log10(rs_o),t200_o,5,ax=ax[0],style=':')
plot_fig(np.log10(rs_o),t500_o,5,ax=ax[0],style='--')
plot_fig(np.log10(rs_o),t1000_o,5,ax=ax[0],style='-')

plot_fig(np.log10(rs_n),t200_n,5,color='C0',ax=ax[1],style=':')
plot_fig(np.log10(rs_n),t500_n,5,color='C0',ax=ax[1],style='--')
plot_fig(np.log10(rs_n),t1000_n,5,color='C0',ax=ax[1],style='-')

ax[0].set_xlabel('$\log(R/R200)$')
ax[1].set_xlabel('$\log(R/R200)$')
ax[0].set_ylabel(r'$\theta \star}$')

plt.savefig(plotspath+'align_medians_stars_gap.png')

#######################

for d in np.arange(29):

    plt.figure()

    D = [d,d+29,d+58]
    
    qdm = DM200.q[mp][D]

    if not all(np.diff(ICL_.D[D]) == 0):
        print('something wrong')

    for i in range(3): 
    
        j = D[i]
    
        
        dif1000 =  ICL_.q[j]/DM1000.q[mp][j]
        dif500  =  ICL_.q[j]/DM500.q[mp][j]
        dif200  =  ICL_.q[j]/DM200.q[mp][j]

        
        plt.plot(ICL_.a[j]/r200[j],dif1000,'C'+str(i+3)+'-',label = p[i]+' '+str(qdm[i]))
        plt.plot(ICL_.a[j]/r200[j],dif500,'C'+str(i+3)+'--')
        plt.plot(ICL_.a[j]/r200[j],dif200,'C'+str(i+3)+':')

        plt.legend(loc=1)

        # if mo[j]:
            # plt.plot(ICL_.a[j]/r200[j],dif1000,'C3o')
            # plt.plot(ICL_.a[j]/r200[j],dif500,'C'+str(i+5)+'o')
            # plt.plot(ICL_.a[j]/r200[j],dif200,'C'+str(i+5)+'o')

    
    
    plt.xlabel('$R/R_{200}$')
    plt.ylabel(r'$q_{ICL}/q_{DM}$')
    
    plt.axis([0,0.8,0,2])
    
    plt.savefig(plotspath+'D'+str(int(ICL_.D[d]))+'_q.png')



plt.figure()
for j in np.arange(87): 
    dif = ICL_.q[j]/DM1000.q[mp][j]

    if mo[j]:
        plt.plot(ICL_.a[j]/r200[j],dif,'r',alpha = 0.5)
    else:
        plt.plot(ICL_.a[j]/r200[j],dif,'C0',alpha = 0.5)

plt.xlabel('$R/R_{200}$')
plt.ylabel(r'$q_{ICL}/q_{DM}$')
plt.savefig(plotspath+'ratio_q_R1000.png')


plt.figure()
for j in np.arange(87): 

    dif = ICL_.q[j]/DM500.q[mp][j]

    if mo[j]:
        plt.plot(ICL_.a[j]/r200[j],dif,'r',alpha = 0.5)
    else:
        plt.plot(ICL_.a[j]/r200[j],dif,'C0',alpha = 0.5)

plt.xlabel('$R/R_{200}$')
plt.ylabel(r'$q_{ICL}/q_{DM}$')
plt.savefig(plotspath+'ratio_q_R500.png')

plt.figure()
for j in np.arange(87): 
    
    dif = ICL_.q[j]/DM200.q[mp][j]
        
    if mo[j]:
        plt.plot(ICL_.a[j]/r200[j],dif,'r',alpha = 0.5)
    else:
        plt.plot(ICL_.a[j]/r200[j],dif,'C0',alpha = 0.5)

plt.xlabel('$R/R_{200}$')
plt.ylabel(r'$q_{ICL}/q_{DM}$')
plt.savefig(plotspath+'ratio_q_R200.png')


rs_o    = np.array([])
qicl    = np.array([])
eqicl    = np.array([])
t200_o  = np.array([])
t500_o  = np.array([])
t1000_o = np.array([])

rs_n    = np.array([])
t200_n  = np.array([])
t500_n  = np.array([])
t1000_n = np.array([])

for j in np.arange(87):
    
    dif200  = ICL_.q[j]/DM200.q[mp][j]    
    dif500  = ICL_.q[j]/DM500.q[mp][j]    
    dif1000 = ICL_.q[j]/DM1000.q[mp][j]    
    
    qicl    = np.append(qicl,ICL_.q[j][-1])
    eqicl    = np.append(qicl,ICL_.eq[j][-1])
    
    if mo[j]:
        rs_o = np.append(rs_o,(ICL_.a/r200)[j])
        t200_o = np.append(t200_o,dif200)
        t500_o = np.append(t500_o,dif500)
        t1000_o = np.append(t1000_o,dif1000)
    else:
        rs_n = np.append(rs_n,(ICL_.a/r200)[j])
        t200_n = np.append(t200_n,dif200)
        t500_n = np.append(t500_n,dif500)
        t1000_n = np.append(t1000_n,dif1000)


f, ax = plt.subplots(1,2, figsize=(12,5), sharex=True,sharey=True)
f.subplots_adjust(hspace=0,wspace=0)
    
plot_fig(np.log10(rs_o),t200_o,6,ax=ax[0],style=':')
plot_fig(np.log10(rs_o),t500_o,6,ax=ax[0],style='--')
plot_fig(np.log10(rs_o),t1000_o,6,ax=ax[0],style='-')

plot_fig(np.log10(rs_n),t200_n,6,color='C0',ax=ax[1],style=':')
plot_fig(np.log10(rs_n),t500_n,6,color='C0',ax=ax[1],style='--')
plot_fig(np.log10(rs_n),t1000_n,6,color='C0',ax=ax[1],style='-')

ax[0].set_xlabel('$\log(R/R200)$')
ax[1].set_xlabel('$\log(R/R200)$')
ax[0].set_ylabel(r'$q_{ICL}/q_{DM}$')

ax[0].set_ylim([0.5,1.2])

plt.savefig(plotspath+'ratio_q_gap.png')

# with stars

maxa  = np.array([])
maxas200 = np.array([])
maxas500 = np.array([])
maxas1000 = np.array([])

rs_o    = np.array([])
t200_o  = np.array([])
t500_o  = np.array([])
t1000_o = np.array([])

rs_n    = np.array([])
t200_n  = np.array([])
t500_n  = np.array([])
t1000_n = np.array([])

for j in np.arange(87):
    
    maxas200  = np.append(maxas200,(ICL_.a/r200)[j][-1])
    maxas500  = np.append(maxas500,(ICL_.a/r500)[j][-1])
    maxas1000 = np.append(maxas1000,(ICL_.a/r1000)[j][-1])
    maxa  = np.append(maxa,ICL_.a[j][-1])
    
    dif200  = ICL_.q[j]/st200.q[mp][j]    
    dif500  = ICL_.q[j]/st500.q[mp][j]    
    dif1000 = ICL_.q[j]/st1000.q[mp][j]    
    
    if mo[j]:
        rs_o = np.append(rs_o,(ICL_.a/r200)[j])
        t200_o = np.append(t200_o,dif200)
        t500_o = np.append(t500_o,dif500)
        t1000_o = np.append(t1000_o,dif1000)
    else:
        rs_n = np.append(rs_n,(ICL_.a/r200)[j])
        t200_n = np.append(t200_n,dif200)
        t500_n = np.append(t500_n,dif500)
        t1000_n = np.append(t1000_n,dif1000)


f, ax = plt.subplots(1,2, figsize=(12,5), sharex=True,sharey=True)
f.subplots_adjust(hspace=0,wspace=0)
    
plot_fig(np.log10(rs_o),t200_o,6,ax=ax[0],style=':')
plot_fig(np.log10(rs_o),t500_o,6,ax=ax[0],style='--')
plot_fig(np.log10(rs_o),t1000_o,6,ax=ax[0],style='-')

plot_fig(np.log10(rs_n),t200_n,6,color='C0',ax=ax[1],style=':')
plot_fig(np.log10(rs_n),t500_n,6,color='C0',ax=ax[1],style='--')
plot_fig(np.log10(rs_n),t1000_n,6,color='C0',ax=ax[1],style='-')

ax[0].set_ylim([0.5,1.2])

ax[0].set_xlabel('$\log(R/R200)$')
ax[1].set_xlabel('$\log(R/R200)$')
ax[0].set_ylabel(r'$q_{ICL}/q_\star$')
plt.savefig(plotspath+'ratio_q_stars_gap.png')


plt.figure()
plt.hist(maxa,histtype='step')
plt.xlabel('R [kpc]')
plt.ylabel('N')
plt.savefig(plotspath+'Rmax_dist.png')

plt.figure()
plt.hist(maxas200,np.linspace(0.1,1.2,15),histtype='step')
plt.xlabel('R/R200')
plt.ylabel('N')
plt.savefig(plotspath+'Rmax_dist_scaled200.png')

plt.figure()
plt.hist(maxas500,np.linspace(0.1,1.2,15),histtype='step')
plt.xlabel('R/R500')
plt.ylabel('N')
plt.savefig(plotspath+'Rmax_dist_scaled500.png')

plt.figure()
plt.hist(maxas1000,np.linspace(0.1,1.2,15),histtype='step')
plt.xlabel('R/R1000')
plt.ylabel('N')
plt.savefig(plotspath+'Rmax_dist_scaled1000.png')


f, ax = plt.subplots(3,2, figsize=(12,12), sharey=True,sharex=True)
f.subplots_adjust(hspace=0,wspace=0)


ax[0,1].text(0.85,0.3,r'$R < R_{200}$')
ax[1,1].text(0.85,0.3,r'$R < R_{500}$')
ax[2,1].text(0.85,0.3,r'$R < R_{1000}$')

ax[0,0].plot(0.76*DM200.q[mp][mo],qicl[mo],'C3.')
ax[1,0].plot(0.76*DM500.q[mp][mo],qicl[mo],'C3.')
ax[2,0].plot(0.76*DM1000.q[mp][mo],qicl[mo],'C3.')
ax[0,0].plot(0.76*DM200.q[mp][~mo],qicl[~mo],'C0.')
ax[1,0].plot(0.76*DM500.q[mp][~mo],qicl[~mo],'C0.')
ax[2,0].plot(0.76*DM1000.q[mp][~mo],qicl[~mo],'C0.')


# ax[0,0].errorbar(DM200.q[mp][mo],qicl[mo],yerr=qicl[mo],ecolor='C3',fmt='none')
# ax[1,0].errorbar(DM500.q[mp][mo],qicl[mo],yerr=qicl[mo],ecolor='C3',fmt='none')
# ax[2,0].errorbar(DM1000.q[mp][mo],qicl[mo],yerr=qicl[mo],ecolor='C3',fmt='none')
# ax[0,0].errorbar(DM200.q[mp][~mo],qicl[~mo],yerr=qicl[~mo],ecolor='C0',fmt='none')
# ax[1,0].errorbar(DM500.q[mp][~mo],qicl[~mo],yerr=qicl[~mo],ecolor='C0',fmt='none')
# ax[2,0].errorbar(DM1000.q[mp][~mo],qicl[~mo],yerr=qicl[~mo],ecolor='C0',fmt='none')


ax[0,1].plot(0.86*st200.q[mp][mo],qicl[mo],'C3.')
ax[1,1].plot(0.86*st500.q[mp][mo],qicl[mo],'C3.')
ax[2,1].plot(0.86*st1000.q[mp][mo],qicl[mo],'C3.')
ax[0,1].plot(0.86*st200.q[mp][~mo],qicl[~mo],'C0.')
ax[1,1].plot(0.86*st500.q[mp][~mo],qicl[~mo],'C0.')
ax[2,1].plot(0.86*st1000.q[mp][~mo],qicl[~mo],'C0.')

# ax[0,1].errorbar(st200.q[mp][mo],qicl[mo],yerr=qicl[mo],ecolor='C3',fmt='none')
# ax[1,1].errorbar(st500.q[mp][mo],qicl[mo],yerr=qicl[mo],ecolor='C3',fmt='none')
# ax[2,1].errorbar(st1000.q[mp][mo],qicl[mo],yerr=qicl[mo],ecolor='C3',fmt='none')
# ax[0,1].errorbar(st200.q[mp][~mo],qicl[~mo],yerr=qicl[~mo],ecolor='C0',fmt='none')
# ax[1,1].errorbar(st500.q[mp][~mo],qicl[~mo],yerr=qicl[~mo],ecolor='C0',fmt='none')
# ax[2,1].errorbar(st1000.q[mp][~mo],qicl[~mo],yerr=qicl[~mo],ecolor='C0',fmt='none')

ax[0,0].plot([0,1],[0,1],'C7--')
ax[1,0].plot([0,1],[0,1],'C7--')
ax[2,0].plot([0,1],[0,1],'C7--')
ax[0,0].plot([0,1],[0,1],'C7--')
ax[1,0].plot([0,1],[0,1],'C7--')
ax[2,0].plot([0,1],[0,1],'C7--')
ax[0,1].plot([0,1],[0,1],'C7--')
ax[1,1].plot([0,1],[0,1],'C7--')
ax[2,1].plot([0,1],[0,1],'C7--')
ax[0,1].plot([0,1],[0,1],'C7--')
ax[1,1].plot([0,1],[0,1],'C7--')
ax[2,1].plot([0,1],[0,1],'C7--')


ax[0,1].set_xlim([0.26,1.0])
ax[0,1].set_ylim([0.26,1.0])

ax[2,0].set_xlabel('$q_{DM}$')
ax[2,1].set_xlabel('$q_{\star}$')

ax[0,0].set_ylabel('$q_{ICL}$')
ax[1,0].set_ylabel('$q_{ICL}$')
ax[2,0].set_ylabel('$q_{ICL}$')

f.savefig(plotspath+'qicl_q_scaled.png')
