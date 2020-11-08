import sys
import numpy as np
from scipy import stats
from pylab import *

def ninbin(binnumber):
    
    N = np.array([])
    for n in np.arange(binnumber.min(),binnumber.max()+1):
        N = np.append(N,sum(binnumber==n))
    
    return N
    
def cosangle(a,b):
    costheta = np.array([])

    for j in range(len(a)):
        cos = np.max([np.dot(-1.*a[j],b[j]),np.dot(a[j],b[j])])
        costheta = np.append(costheta,cos)
    
    costheta = np.round(costheta,3)
    theta = np.rad2deg(np.arccos(costheta))
    
    return costheta,theta


def binned(x,y,nbins=10):
    
    bined = stats.binned_statistic(x,y,statistic='median', bins=nbins)
    x_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    y_b     = bined.statistic
    bined = stats.binned_statistic(x,y,statistic='std', bins=nbins)
    sy_b     = bined.statistic/np.sqrt(ninbin(bined.binnumber))
    return x_b,y_b,sy_b

def plot_binned(X,Y,label,color='C3',style='',nbins=10):
    x,y,s = binned(X,Y,nbins)
    plt.errorbar(x,y,yerr=s,fmt=style,ecolor=color,color=color,label=label)

class Shape:
    
    def __init__(self, ind=2, name_cat='../catalog/dm_091.dat'):
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
        self.q = np.concatenate((qx,qy,qz))
        
        self.a3D = np.array([cat[ind+6],cat[ind+7],cat[ind+8]]).T
        
        a2Dx = np.array([cat[ind+9] ,cat[ind+10]]).T
        a2Dy = np.array([cat[ind+11],cat[ind+12]]).T
        a2Dz = np.array([cat[ind+13],cat[ind+14]]).T
        
        self.a2D  = np.concatenate((a2Dx,a2Dy,a2Dz))
    
class DarkMatter(Shape):
    
    def __init__(self,nbins=10):
        
        self.name_cat = '../catalog/dm_091.dat'      
        ind = 2
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
            
        self.name_cat = '../catalog/stars_091.dat'      
        Shape.__init__(self, ind=ind, name_cat=self.name_cat)

class Galaxias(Shape):

    def __init__(self, radio, masa_cut= False):
        
        self.name_cat = '../catalog/glxs_091.dat'

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
        Shape.__init__(self, ind=ind, name_cat=self.name_cat)

def Tenneti(lM,a):
    lMpiv = np.log10(1e12)
    param = 0
    for i in range(len(a)):
        param += a[i]*((lM - lMpiv)**i)
    return param


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
