import sys
import numpy as np
from scipy import stats
from pylab import *
from scipy.optimize import curve_fit
from scipy import integrate

# path = '../catalog/nuevosdats/'
path = '../catalog/'

def g(x,disp): 
    return (1./(np.sqrt(2*np.pi)*disp))*np.exp(-0.5*(x/disp)**2)

def q_75(y):
    return np.quantile(y, 0.75)

def q_25(y):
    return np.quantile(y, 0.25)


class Dcompute():
    
    
    def __init__(self, theta):
    
        m = (theta < 90.)
        
        theta = theta[m]
    
        theta = np.append(theta,-1.*theta)
        
        n,c = np.histogram(theta,50,density=True)
        c   = (c+(c[1]-c[0])*0.5)[:-1]
        
        D    = np.mean(np.cos(2.*(np.deg2rad(theta))))
        Dstd = np.std(np.cos(2.*(np.deg2rad(theta))))
        
        fitg = curve_fit(g,c,n,sigma=np.ones(len(c)),absolute_sigma=True)
        std_fit = fitg[0][0]
        
        argumento = lambda x: g(x,std_fit)*np.cos(2*np.deg2rad(x))
        Dfit = integrate.quad(argumento, -90., 90.)[0]
        
        argumento = lambda x: g(x,np.std(theta))*np.cos(2*np.deg2rad(x))
        Dg  = integrate.quad(argumento, -90., 90.)[0]
    
        self.theta   = theta
        self.std_fit = std_fit
        self.std     = np.std(theta)
        self.D       = D
        self.Dfit    = Dfit
        self.Dstd    = Dstd
        self.Dg      = Dg


def ninbin(binnumber):
    
    N = np.array([])
    for n in np.arange(binnumber.min(),binnumber.max()+1):
        N = np.append(N,sum(binnumber==n))
    
    return N
    
def cosangle(a,b):
    costheta = np.array([])

    for j in range(len(a)):
        cos = np.abs(np.dot(a[j],b[j]))
        costheta = np.append(costheta,cos)
    
    costheta[abs(costheta) > 1.] = 1.
    theta = np.rad2deg(np.arccos(costheta))
    
    return costheta,theta

def cosangle2(a,b):
    costheta = np.array([])

    for j in range(len(a)):
        cos = np.dot(a[j],b[j])
        costheta = np.append(costheta,cos)
    
    costheta[abs(costheta) > 1.] = 1.
    theta = np.rad2deg(np.arccos(costheta))
    
    return costheta,theta


def binned(x,y,nbins=10):
    
    
    bined = stats.binned_statistic(x,y,statistic='median', bins=nbins)
    x_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    q50     = bined.statistic
    
    bined = stats.binned_statistic(x,y,statistic=q_25, bins=nbins)
    q25     = bined.statistic

    bined = stats.binned_statistic(x,y,statistic=q_75, bins=nbins)
    q75     = bined.statistic
    
    dig   = np.digitize(x,bined.bin_edges)
    mz    = np.ones(len(x))
    for j in range(nbins):
        mbin = dig == (j+1)
        mz[mbin] = y[mbin] >= q50[j]   
    mz = mz.astype(bool)
    return x_b,q50,q25,q75,mz
            

def plot_binned(X,Y,label,color='C3',style='',nbins=10):
    x,y,s,m = binned(X,Y,nbins)
    plt.errorbar(x,y,yerr=s,fmt=style,ecolor=color,color=color,label=label)
    

def newold(indicator = 'gap',gxs = False):
    
    gral  = np.loadtxt(path+'gral_091.dat').T
    
    lM = np.log10(gral[9])
    lMp = np.array((lM.tolist())*3)
    
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
    
    if 'gap' in indicator:
        mold = gap <= 0.22
        mold2D = gap2D <= 0.22
        
        mnew = ~mold
        mnew2D = ~mold2D
    
    if gxs:
            mN = Galaxias(radio=1000).N > 9
            mN2D = np.array((Galaxias(radio=1000).N.tolist())*3) > 9
            mnew = mnew[mN]
            mnew2D = mnew2D[mN2D]
            mold = mold[mN]
            mold2D = mold2D[mN2D] 
    
    return mnew,mold,mnew2D,mold2D
    
def plot_fig(x,y,nbins,ax=plt,color = 'sienna', style = '',label=''):
                
        X,q50,q25,q75,mz = binned(x,y,nbins)
        ax.plot(X,q50,style,color=color,label=label)
        ax.plot(X,q75,style,color=color,alpha=0.2)
        ax.plot(X,q25,style,color=color,alpha=0.2)
        ax.fill_between(X,q75,q25,color = color,alpha=0.1)


class Shape:
    
    def __init__(self, ind=2, name_cat=path+'dm_091.dat'):
        self.ind = ind
        self.name_cat = name_cat
        cat = np.loadtxt(name_cat).T
    
        a  = cat[ind+0]
        b  = cat[ind+1]
        c  = cat[ind+2]
        qx = cat[ind+3]
        qy = cat[ind+4]
        qz = cat[ind+5]
        
        self.T = (a**2 - b**2)/(a**2 - c**2)
        self.S = c/a
        self.Q = b/a
        self.q = np.concatenate((qx,qy,qz))
        
        self.a3D = np.array([cat[ind+6],cat[ind+7],cat[ind+8]]).T
        
        a2Dx = np.array([cat[ind+9] ,cat[ind+10]]).T
        a2Dy = np.array([cat[ind+11],cat[ind+12]]).T
        a2Dz = np.array([cat[ind+13],cat[ind+14]]).T
        
        self.a2D  = np.concatenate((a2Dx,a2Dy,a2Dz))
    
class DarkMatter(Shape):
    
    def __init__(self, radio,subhalos=True):
        
        self.name_cat = path+'dm_091.dat'      

        if radio == 30:
            ind = 2
        elif radio == 50:
            ind = 2 + (15*2)
        elif radio == 100:
            ind = 2 + (15*4)
        elif radio == 1000:
            ind = 2 + (15*6)
        elif radio == 500:
            ind = 2 + (15*8)
        elif radio == 200:
            ind = 2 + (15*10)
        
        if not subhalos:
            ind = ind + 15
        Shape.__init__(self, ind=ind, name_cat=self.name_cat)

class Stars(Shape):
    
    def __init__(self, radio):
        
        if radio == 30:
            ind = 2
        elif radio == 50:
            ind = 2+15
        elif radio == 100:
            ind = 2+15*2
        elif radio == 1000:
            ind = 2+15*3
        elif radio == 500:
            ind = 2+15*4
        elif radio == 200:
            ind = 2+15*5
            
        self.name_cat = path+'stars_091.dat'      
        Shape.__init__(self, ind=ind, name_cat=self.name_cat)

class Galaxias(Shape):

    def __init__(self, tipo = 'all', radio = 200, masa_cut= True):
        
        if tipo == 'all':
        
            self.name_cat = path+'glxs_091.dat'
    
            if radio == 1000:
                ind = 2
            elif radio == 500:
                ind = 2+19*2
            elif radio == 200:
                ind = 2+19*2*2
            if masa_cut:
                ind += 19
            
            gal = np.loadtxt(self.name_cat).T
            self.N  = gal[ind]
            nx = gal[ind+1]
            ny = gal[ind+2]
            nz = gal[ind+3]
            self.n = np.concatenate((nx,ny,nz))
            
            ind = ind + 4
        
        else:
            self.name_cat = path+'glxs_hmr_091.dat'
            gal = np.loadtxt(self.name_cat).T
            ind = 2
            self.N  = gal[ind]
            nx = gal[ind+1]
            ny = gal[ind+2]
            nz = gal[ind+3]
            self.n = np.concatenate((nx,ny,nz))
            
            
            if tipo == 'concentradas':
                ind = ind + 4
            elif tipo == 'extendidas':
                ind = ind + 4 + 15
            else:
                print('non correct type selection')
                ind = 999
        
        Shape.__init__(self, ind=ind, name_cat=self.name_cat)

class Random():

    def __init__(self, radio = 200, half = False):
        
        if half:
            gal = np.loadtxt(path+'dmrand_halfnglxs_091.dat').T
        else:
            gal = np.loadtxt(path+'dmrand_091.dat').T
        
        if radio == 1000:
            ind = 2
        elif radio == 500:
            ind = 2+69
        elif radio == 200:
            ind = 2+69*2
            
        self.N  = gal[ind]
        nx = gal[ind+1]
        ny = gal[ind+2]
        nz = gal[ind+3]
        self.n = np.concatenate((nx,ny,nz))
        
        self.T25 = gal[ind+4]
        self.T50 = gal[ind+5]
        self.T75 = gal[ind+6]

        self.S25 = gal[ind+7]
        self.S50 = gal[ind+8]
        self.S75 = gal[ind+9]

        self.tgx25 = gal[ind+10]
        self.tgx50 = gal[ind+11]
        self.tgx75 = gal[ind+12]

        self.tdm25 = gal[ind+13]
        self.tdm50 = gal[ind+14]
        self.tdm75 = gal[ind+15]
        
        self.T_m = gal[ind+16]
        self.S_m = gal[ind+17]
        self.tgx_m = gal[ind+18]
        self.tdm_m = gal[ind+19]
        
        self.T_std = gal[ind+20]
        self.S_std = gal[ind+21]
        self.tgx_std = gal[ind+22]
        self.tdm_std = gal[ind+23]

        self.q25 = np.concatenate((gal[ind+24],gal[ind+24+15],gal[ind+24+15*2]))
        self.q50 = np.concatenate((gal[ind+25],gal[ind+25+15],gal[ind+25+15*2]))
        self.q75 = np.concatenate((gal[ind+26],gal[ind+26+15],gal[ind+26+15*2]))

        self.t2Dgx25 = np.concatenate((gal[ind+27],gal[ind+27+15],gal[ind+27+15*2]))
        self.t2Dgx50 = np.concatenate((gal[ind+28],gal[ind+28+15],gal[ind+28+15*2]))
        self.t2Dgx75 = np.concatenate((gal[ind+29],gal[ind+29+15],gal[ind+29+15*2]))

        self.t2Ddm25 = np.concatenate((gal[ind+30],gal[ind+30+15],gal[ind+30+15*2]))
        self.t2Ddm50 = np.concatenate((gal[ind+31],gal[ind+31+15],gal[ind+31+15*2]))
        self.t2Ddm75 = np.concatenate((gal[ind+32],gal[ind+32+15],gal[ind+32+15*2]))

        self.q_m       = np.concatenate((gal[ind+33],gal[ind+33+15],gal[ind+33+15*2]))
        self.t2Dgx_m   = np.concatenate((gal[ind+34],gal[ind+34+15],gal[ind+34+15*2]))
        self.t2Ddm_m   = np.concatenate((gal[ind+35],gal[ind+35+15],gal[ind+35+15*2]))

        self.q_std       = np.concatenate((gal[ind+36],gal[ind+36+15],gal[ind+36+15*2]))
        self.t2Dgx_std   = np.concatenate((gal[ind+37],gal[ind+37+15],gal[ind+37+15*2]))
        self.t2Ddm_std   = np.concatenate((gal[ind+38],gal[ind+38+15],gal[ind+38+15*2]))
        
        
        

def Tenneti(lM,a):
    lMpiv = np.log10(1e12)
    param = 0
    for i in range(len(a)):
        param += a[i]*((lM - lMpiv)**i)
    return param
    

class Clusters:
    
    def __init__(self):

        gral  = np.loadtxt(path+'gral_091.dat').T
        times  = np.loadtxt('../catalog/relevant_times.dat').T
        
        ltime = np.ones(len(gral.T))*-999.

        for j in range(len(gral.T)):
            m = (gral[0,j] == times[0])*(gral[2,j] == times[1])
            if m.sum() == 1:
                ltime[j] = times[3,m]
        
        self.D     = gral[0]
        self.sub   = gral[1]
        self.ltime   = ltime
        
        self.R1000 = gral[4]
        self.R500  = gral[5]
        self.R200  = gral[6]

        self.lM1000 = np.log10(gral[7])
        self.lM500  = np.log10(gral[8])
        self.lM200  = np.log10(gral[9])
        self.lM30   = np.log10(gral[10])
        self.lM50   = np.log10(gral[11])
        self.lM100  = np.log10(gral[12])

        self.lM1000p = np.array((np.log10(gral[7]).tolist())*3)
        self.lM500p  = np.array((np.log10(gral[8]).tolist())*3)
        self.lM200p  = np.array((np.log10(gral[9]).tolist())*3)
        self.lM30p   = np.array((np.log10(gral[10]).tolist())*3)
        self.lM50p   = np.array((np.log10(gral[11]).tolist())*3)
        self.lM100p  = np.array((np.log10(gral[12]).tolist())*3)


        self.off  = gral[13]
        self.offp = np.concatenate((gral[14],gral[15],gral[16]))
        self.DV   = gral[17]
        self.DVp  = np.concatenate((gral[18],gral[19],gral[20]))
        self.gap  = gral[21]
        self.gapp = np.concatenate((gral[22],gral[23],gral[24]))
        
        m_gap = newold('gap')
        m_off = newold('off')
        m_dv = newold('DV')
        
        self.mn_gap   = m_gap[0]
        self.mo_gap   = m_gap[1]
        self.mn2d_gap = m_gap[2]
        self.mo2d_gap = m_gap[3]

        self.mn_off   = m_off[0]
        self.mo_off   = m_off[1]
        self.mn2d_off = m_off[2]
        self.mo2d_off = m_off[3]

        self.mn_dv   = m_dv[0]
        self.mo_dv   = m_dv[1]
        self.mn2d_dv = m_dv[2]
        self.mo2d_dv = m_dv[3]
        
        R30   = np.ones(gral.shape[1])*30.
        R50   = np.ones(gral.shape[1])*50.
        R100    = 0.1*gral[5]

        
        
        self.Rs = np.vstack((R30/gral[6],R50/gral[6],R100/gral[6],gral[4]/gral[6],gral[5]/gral[6],np.ones(gral.shape[1]))).T
        self.R  = np.vstack((R30,R50,R100,gral[4],gral[5],gral[6])).T
        self.Rp = np.array((self.R.tolist())*3)
        self.Rsp = np.array((self.Rs.tolist())*3)
        
        self.ltimep = np.array((ltime.tolist())*3)
        
class ICL:
    
    def __init__(self,trazer='icl',cutmax=False):
        
        
        patht = '../catalog/tablas_'+trazer+'/'
        pathicl    = '../catalog/tablas_icl/'
               
        D    = np.array([])
        proj = np.array([])
        A    = []
        q    = []
        eq   = []
        PA   = []
        ePA  = []

        pbcg  = []
        A0    = []
        q0    = []
        eq0   = []
        PA0   = []
        ePA0  = []
        
        qicl      = np.array([])
        eqicl     = np.array([])
        PAicl     = np.array([])
        ePAicl    = np.array([])
        
        Ricl = np.array([])
        Rbcg = np.array([])
        
        xicl = np.array([])
        yicl = np.array([])
        
        xbcg = np.array([])
        ybcg = np.array([])
        
        qbcg   = np.array([])
        eqbcg  = np.array([])
        PAbcg  = np.array([])
        ePAbcg = np.array([])
        
        C = Clusters()
        m = (C.sub == 0)
        mp = (m.tolist()*3)

        r200  = C.Rp[mp,-1]
        r100  = C.Rp[mp,2]
        r1000  = C.Rp[mp,3]

        
        p    = ['xy','xz','yz']
        
        csample = np.arange(1,30)
        
        for i in range(3):
        
            for j in csample:
                
                ## icl data
                
                tab = np.loadtxt(patht+'d'+str(int(j))+'_'+p[i]+'.dat',dtype='str').T
                tab_icl = np.loadtxt(pathicl+'d'+str(int(j))+'_'+p[i]+'.dat',dtype='str').T
            
                mindef = ~(tab[-5] == 'INDEF')
                tab[-5][~mindef] = '99.'
                musar  = (tab[-5].astype(float) < 90.)

                mindef_icl = ~(tab_icl[-5] == 'INDEF')
                tab_icl[-5][~mindef_icl] = '99.'
                musar_icl  = (tab_icl[-5].astype(float) < 90.)
            
                D    = np.append(D,j)
                proj = np.append(proj,i)
                
                try:
                    a, i0, ei0, e, ee, pa, epa, x0, ex0, y0, ey0 = tab[:,musar].astype(float)                
                except:
                    a, i0, e, ee, pa, epa, x0, ex0, y0, ey0 = tab[:,musar].astype(float)                
                                    
                
                aicl = tab_icl[0,musar_icl].astype(float)                

                
                r = (1+a)*20.
                aicl = (1+aicl[-1])*20.
                A   = A + [r]
                PA  = PA + [pa]
                ePA = ePA + [epa]
                q   = q + [1.-e]
                eq = eq + [ee]

                A0   = A0 + [(1+a[-1])*20.]
                PA0  = PA0 + [pa[-1]]
                ePA0 = ePA0 + [epa[-1]]
                q0   = q0 + [1.-e[-1]]
                eq0 = eq0 + [ee[-1]]
                
                if cutmax:
                    micl = (r <= aicl)#(r > r100[j-1])
                else:
                    micl = (r <= r1000[j-1])#(r > r100[j-1])
                mbcg = r < r100[j-1]
                # try:
                    # pbcg = pbcg + [max(a[mbcg])]
                    # print(i,j,max(a[mbcg]))
                # except:
                    # print('no bcg')
                
                if any(pa > 80.):
                    pa[pa < 0.] = 180. + pa[pa < 0.]
                
                
                paicl = np.mean(pa[micl])
                pabcg = np.mean(pa[mbcg])

                ricl = np.max(r[micl]/r200[j-1])
                
                if mbcg.sum() > 2.:
                    rbcg = np.max(r[mbcg]/r200[j-1])
                else:
                    rbcg = 0.

                if paicl < 0.:
                    paicl += 360.
                if paicl > 180.:
                    paicl -= 180.

                if pabcg < 0.:
                    pabcg += 360.
                if pabcg > 180.:
                    pabcg -= 180.

                
                qicl     = np.append(qicl,np.mean((1.-e)[micl]))
                eqicl    = np.append(eqicl,np.std((1.-e)[micl]))

                Ricl    = np.append(Ricl,ricl)
                Rbcg    = np.append(Rbcg,rbcg)
                
                PAicl    = np.append(PAicl,paicl)
                ePAicl   = np.append(ePAicl,np.std(epa[micl]))
            
                qbcg = np.append(qbcg,np.mean((1.-e)[mbcg]))
                eqbcg = np.append(eqbcg,np.std((1.-e)[mbcg]))
                PAbcg    = np.append(PAbcg,np.mean(pa[mbcg]))
                ePAbcg   = np.append(ePAbcg,np.std(epa[mbcg]))
                
                xbcg = np.append(xbcg,-1.*np.sin(np.deg2rad(PAbcg[-1])))
                ybcg = np.append(ybcg,np.cos(np.deg2rad(PAbcg[-1])))

                xicl = np.append(xicl,-1.*np.sin(np.deg2rad(PAicl[-1])))
                yicl = np.append(yicl,np.cos(np.deg2rad(PAicl[-1])))
                
                
        a2Dicl = np.zeros((len(xicl),2))
        a2Dbcg = np.zeros((len(xicl),2))
        
        a2Dicl[:,0] = xicl
        a2Dicl[:,1] = yicl

        a2Dbcg[:,0] = xbcg
        a2Dbcg[:,1] = ybcg

        self.a2Dicl = a2Dicl
        self.a2Dbcg = a2Dbcg

        self.Ricl = Ricl
        self.Rbcg = Rbcg

        self.qicl  = qicl
        self.eqicl = eqicl
        self.PAicl = PAicl
        self.ePAicl = ePAicl

        self.qbcg  = qbcg
        self.eqbcg = eqbcg
        self.PAbcg = PAbcg
        self.ePAbcg = ePAbcg
        
        
        self.D    = D
        self.proj = proj
        self.a    = A
        self.q    = q
        self.eq   = eq
        self.PA   = PA
        self.ePA  = ePA

        self.a0   = np.array(A0)
        self.q0   = np.array(q0)
        self.eq0  = np.array(eq0)
        self.PA0  = np.array(PA0)
        self.ePA0 = np.array(ePA0)
        
        
                
            

class CorrelR():
    
    def __init__(self,trazer1,trazer2,mN = [],mN2D = []):
    
        gral  = np.loadtxt(path+'gral_091.dat').T
        
        
        t1_200  = trazer1(radio=200)
        t1_500  = trazer1(radio=500)
        t1_1000 = trazer1(radio=1000)
    
        t2_30   = trazer2(radio=30)
        t2_50   = trazer2(radio=50)
        t2_100  = trazer2(radio=100)
        t2_200  = trazer2(radio=200)
        t2_500  = trazer2(radio=500)
        t2_1000 = trazer2(radio=1000)

        if not len(mN):
            mN   = (np.ones(len(t1_200.S)).astype(bool))
            mN2D = (np.ones(len(t1_200.q)).astype(bool))
        
        try:
            t1_30   = trazer1(radio=30)
            t1_50   = trazer1(radio=50)
            t1_100  = trazer1(radio=100)
        except:
            t1_30   = t2_30
            t1_50   = t2_50
            t1_100  = t2_100

            
        dm   = DarkMatter(1000)
        dm200   = DarkMatter(200)
    
        gx_c   = Galaxias(tipo='concentradas')
        gx_e   = Galaxias(tipo='extendidas')


        ct30_3D  , t30_3D     = cosangle(t1_30.a3D[mN],t2_30.a3D[mN]) 
        ct50_3D  , t50_3D     = cosangle(t1_50.a3D[mN],t2_50.a3D[mN])
        ct100_3D , t100_3D    = cosangle(t1_100.a3D[mN],t2_100.a3D[mN])
        ct200_3D , t200_3D    = cosangle(t1_200.a3D[mN],t2_200.a3D[mN])
        ct500_3D , t500_3D    = cosangle(t1_500.a3D[mN],t2_500.a3D[mN])
        ct1000_3D, t1000_3D   = cosangle(t1_1000.a3D[mN],t2_1000.a3D[mN])
        
        ct30_2D  , t30_2D     = cosangle(t1_30.a2D[mN2D],t2_30.a2D[mN2D]) 
        ct50_2D  , t50_2D     = cosangle(t1_50.a2D[mN2D],t2_50.a2D[mN2D])
        ct100_2D , t100_2D    = cosangle(t1_100.a2D[mN2D],t2_100.a2D[mN2D])
        ct200_2D , t200_2D    = cosangle(t1_200.a2D[mN2D],t2_200.a2D[mN2D])
        ct500_2D , t500_2D    = cosangle(t1_500.a2D[mN2D],t2_500.a2D[mN2D])
        ct1000_2D, t1000_2D   = cosangle(t1_1000.a2D[mN2D],t2_1000.a2D[mN2D])
        
        t3D = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
        t2D = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T


        ct30_3D  , t30_3D     = cosangle(t1_30.a3D[mN],t2_50.a3D[mN]) 
        ct50_3D  , t50_3D     = cosangle(t1_50.a3D[mN],t2_100.a3D[mN])
        ct100_3D , t100_3D    = cosangle(t1_100.a3D[mN],t2_1000.a3D[mN])
        ct1000_3D, t1000_3D   = cosangle(t1_1000.a3D[mN],t2_500.a3D[mN])
        ct500_3D , t500_3D    = cosangle(t1_500.a3D[mN],t2_200.a3D[mN])
        ct200_3D , t200_3D    = cosangle(t1_200.a3D[mN],t2_200.a3D[mN])

        ct30_2D  , t30_2D     = cosangle(t1_30.a2D[mN2D],t2_50.a2D[mN2D]) 
        ct50_2D  , t50_2D     = cosangle(t1_50.a2D[mN2D],t2_100.a2D[mN2D])
        ct100_2D , t100_2D    = cosangle(t1_100.a2D[mN2D],t2_1000.a2D[mN2D])
        ct1000_2D, t1000_2D   = cosangle(t1_1000.a2D[mN2D],t2_500.a2D[mN2D])
        ct500_2D , t500_2D    = cosangle(t1_500.a2D[mN2D],t2_200.a2D[mN2D])
        ct200_2D , t200_2D    = cosangle(t1_200.a2D[mN2D],t2_200.a2D[mN2D])
                
        t3D_2 = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
        t2D_2 = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T

    
        ct30_3D  , t30_3D     = cosangle(dm.a3D[mN],t1_30.a3D[mN]) 
        ct50_3D  , t50_3D     = cosangle(dm.a3D[mN],t1_50.a3D[mN])
        ct100_3D , t100_3D    = cosangle(dm.a3D[mN],t1_100.a3D[mN])
        ct200_3D , t200_3D    = cosangle(dm.a3D[mN],t1_200.a3D[mN])
        ct500_3D , t500_3D    = cosangle(dm.a3D[mN],t1_500.a3D[mN])
        ct1000_3D, t1000_3D   = cosangle(dm.a3D[mN],t1_1000.a3D[mN])
        
        ct30_2D  , t30_2D     = cosangle(dm.a2D[mN2D],t1_30.a2D[mN2D]) 
        ct50_2D  , t50_2D     = cosangle(dm.a2D[mN2D],t1_50.a2D[mN2D])
        ct100_2D , t100_2D    = cosangle(dm.a2D[mN2D],t1_100.a2D[mN2D])
        ct200_2D , t200_2D    = cosangle(dm.a2D[mN2D],t1_200.a2D[mN2D])
        ct500_2D , t500_2D    = cosangle(dm.a2D[mN2D],t1_500.a2D[mN2D])
        ct1000_2D, t1000_2D   = cosangle(dm.a2D[mN2D],t1_1000.a2D[mN2D])
        
        t3D_dm = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
        t2D_dm = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T

        ct30_3D  , t30_3D     = cosangle(dm200.a3D[mN],t1_30.a3D[mN]) 
        ct50_3D  , t50_3D     = cosangle(dm200.a3D[mN],t1_50.a3D[mN])
        ct100_3D , t100_3D    = cosangle(dm200.a3D[mN],t1_100.a3D[mN])
        ct200_3D , t200_3D    = cosangle(dm200.a3D[mN],t1_200.a3D[mN])
        ct500_3D , t500_3D    = cosangle(dm200.a3D[mN],t1_500.a3D[mN])
        ct1000_3D, t1000_3D   = cosangle(dm200.a3D[mN],t1_1000.a3D[mN])
        
        ct30_2D  , t30_2D     = cosangle(dm200.a2D[mN2D],t1_30.a2D[mN2D]) 
        ct50_2D  , t50_2D     = cosangle(dm200.a2D[mN2D],t1_50.a2D[mN2D])
        ct100_2D , t100_2D    = cosangle(dm200.a2D[mN2D],t1_100.a2D[mN2D])
        ct200_2D , t200_2D    = cosangle(dm200.a2D[mN2D],t1_200.a2D[mN2D])
        ct500_2D , t500_2D    = cosangle(dm200.a2D[mN2D],t1_500.a2D[mN2D])
        ct1000_2D, t1000_2D   = cosangle(dm200.a2D[mN2D],t1_1000.a2D[mN2D])
        
        t3D_dm200 = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
        t2D_dm200 = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T

        ct30_3D  , t30_3D     = cosangle(gx_c.a3D[mN],t1_30.a3D[mN]) 
        ct50_3D  , t50_3D     = cosangle(gx_c.a3D[mN],t1_50.a3D[mN])
        ct100_3D , t100_3D    = cosangle(gx_c.a3D[mN],t1_100.a3D[mN])
        ct200_3D , t200_3D    = cosangle(gx_c.a3D[mN],t1_200.a3D[mN])
        ct500_3D , t500_3D    = cosangle(gx_c.a3D[mN],t1_500.a3D[mN])
        ct1000_3D, t1000_3D   = cosangle(gx_c.a3D[mN],t1_1000.a3D[mN])
                                         
        ct30_2D  , t30_2D     = cosangle(gx_c.a2D[mN2D],t1_30.a2D[mN2D]) 
        ct50_2D  , t50_2D     = cosangle(gx_c.a2D[mN2D],t1_50.a2D[mN2D])
        ct100_2D , t100_2D    = cosangle(gx_c.a2D[mN2D],t1_100.a2D[mN2D])
        ct200_2D , t200_2D    = cosangle(gx_c.a2D[mN2D],t1_200.a2D[mN2D])
        ct500_2D , t500_2D    = cosangle(gx_c.a2D[mN2D],t1_500.a2D[mN2D])
        ct1000_2D, t1000_2D   = cosangle(gx_c.a2D[mN2D],t1_1000.a2D[mN2D])
        
        t3D_gxc = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
        t2D_gxc = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T

        ct30_3D  , t30_3D     = cosangle(gx_e.a3D[mN],t1_30.a3D[mN]) 
        ct50_3D  , t50_3D     = cosangle(gx_e.a3D[mN],t1_50.a3D[mN])
        ct100_3D , t100_3D    = cosangle(gx_e.a3D[mN],t1_100.a3D[mN])
        ct200_3D , t200_3D    = cosangle(gx_e.a3D[mN],t1_200.a3D[mN])
        ct500_3D , t500_3D    = cosangle(gx_e.a3D[mN],t1_500.a3D[mN])
        ct1000_3D, t1000_3D   = cosangle(gx_e.a3D[mN],t1_1000.a3D[mN])
                                         
        ct30_2D  , t30_2D     = cosangle(gx_e.a2D[mN2D],t1_30.a2D[mN2D]) 
        ct50_2D  , t50_2D     = cosangle(gx_e.a2D[mN2D],t1_50.a2D[mN2D])
        ct100_2D , t100_2D    = cosangle(gx_e.a2D[mN2D],t1_100.a2D[mN2D])
        ct200_2D , t200_2D    = cosangle(gx_e.a2D[mN2D],t1_200.a2D[mN2D])
        ct500_2D , t500_2D    = cosangle(gx_e.a2D[mN2D],t1_500.a2D[mN2D])
        ct1000_2D, t1000_2D   = cosangle(gx_e.a2D[mN2D],t1_1000.a2D[mN2D])
        
        t3D_gxe = np.vstack((t30_3D,t50_3D,t100_3D,t1000_3D,t500_3D,t200_3D)).T
        t2D_gxe = np.vstack((t30_2D,t50_2D,t100_2D,t1000_2D,t500_2D,t200_2D)).T

    
        S = np.vstack((t1_30.S[mN],t1_50.S[mN],t1_100.S[mN],t1_1000.S[mN],t1_500.S[mN],t1_200.S[mN])).T
        T = np.vstack((t1_30.T[mN],t1_50.T[mN],t1_100.T[mN],t1_1000.T[mN],t1_500.T[mN],t1_200.T[mN])).T
        q = np.vstack((t1_30.q[mN2D],t1_50.q[mN2D],t1_100.q[mN2D],t1_1000.q[mN2D],t1_500.q[mN2D],t1_200.q[mN2D])).T
        
        self.Sr = S
        self.Tr = T
        self.qr = q
        
        self.t3D_dm200 = t3D_dm200
        self.t2D_dm200 = t2D_dm200
        
        self.t3D_dm1000 = t3D_dm
        self.t2D_dm1000 = t2D_dm

        self.t3D_gxe = t3D_gxe
        self.t2D_gxe = t2D_gxe

        self.t3D_gxc = t3D_gxc
        self.t2D_gxc = t2D_gxc
        
        self.t3D = t3D
        self.t2D = t2D
        self.t3D_2 = t3D_2
        self.t2D_2 = t2D_2
        


def plotR_ind(R,S,color,label,style='',ax=plt):
    ax.plot(np.median(R,axis=0),np.median(S,axis=0),style,color=color,label = label)
    # ax.fill_between(np.median(np.log10(R),axis=0),np.median(S,axis=0)+np.std(S,axis=0),np.median(S,axis=0)-np.std(S,axis=0),color = color,alpha=0.1)
    ax.plot(np.median(R,axis=0),np.quantile(S, 0.75, axis=0),style,color=color,alpha=0.2)
    ax.plot(np.median(R,axis=0),np.quantile(S, 0.25, axis=0),style,color=color,alpha=0.2)
    ax.fill_between(np.median(R,axis=0),np.quantile(S, 0.75, axis=0),np.quantile(S, 0.25, axis=0),color = color,alpha=0.1)


def plotR(R,S,param_name,infile,
          path_plots = '../plots/correl_dM_radio/' ,
          split_age=True,indicator='gap',gxs = False,
          label0 = 'all',label1 = 'non-relaxed',
          label2 = 'relaxed',limitesy = None):
    
    plt.figure()    
    
    if gxs:
        mask = np.array([0,0,0,1,1,1]).astype(bool)
    else:
        mask = np.ones(6).astype(bool)
    
    plotR_ind(R[:,mask],S[:,mask],'k',label0)
    
    
    if not limitesy:
        if S.max() < 2.:
            limitesy = [0.4,1.0]
        else:
            limitesy = [0,50]
        
    if gxs:
        limites = [2.5,np.log10(2000)] + limitesy
    else:
        limites = [np.log10(30),np.log10(2000)] + limitesy
    
    if split_age:
    
        mnew,mold,mnew2D,mold2D = newold(indicator,gxs)
        
        if len(R) == len(mnew2D):
            mnew = mnew2D
            mold = mold2D
        
        plotR_ind(R[:,mask][mnew],S[:,mask][mnew],'C0',label1)
        plotR_ind(R[:,mask][mold],S[:,mask][mold],'C3',label2)
    
    Rlegend = np.array(['30kpc','50kpc','0.1R500','R1000','R500','R200'])
    
    plt.plot(limites[:-2],[1,1],'C7--')
    plt.legend(loc=2,frameon=False)
    plt.ylabel(param_name)
    plt.xticks(np.median(np.log10(R),axis=0)[mask],Rlegend[mask])
    plt.axis(limites)
    plt.savefig(path_plots+infile+'_R'+indicator+'.png')



def q_dm(lM):
    aq_dm = [0.797,-0.049]
    return Tenneti(lM,aq_dm)
    
def s_dm(lM):
    as_dm = [0.663,-0.059]
    return Tenneti(lM,as_dm)
    
def q_stars(lM):
    aq_s = [0.771,-0.004,-0.068,-0.017,-0.061,-0.003,-0.015]
    return Tenneti(lM,aq_s)
    
def s_stars(lM):
    as_s = [-0.585, 0.031, -0.089, -0.034, 0.075, -0.001, -0.016]
    return Tenneti(lM,as_s)
