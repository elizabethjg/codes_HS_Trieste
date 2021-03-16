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
from scipy.optimize import curve_fit

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 14})

def f(erN,a):
   er,N = erN
   
   fraction = (1./N)*(1./er**3.7)*a + 1.
   et = er/fraction
   
   return et

def fN(er,a,b):
   
   fraction = (1./er**b)*a + 1.
   et = er/fraction
   
   return et


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


mqL = DM200.q <= 0.7
mN = gx200.n > 200.

XqL,YqL,q25,q75,m = binned(np.log10(r200.n)[mqL],(r200.q50/DM200.q)[mqL],6)
XN,YN,q25,q75,m = binned(((1.-r200.q50)/(1.+r200.q50))[mN],(r200.q50/DM200.q)[mN],6)

er = ((1.-r200.q50)/(1.+r200.q50))
et = ((1.-DM200.q)/(1.+DM200.q))
N  = r200.n

fita = curve_fit(f,(er,N),et,sigma=np.ones(len(er)),absolute_sigma=True)
a    = fita[0][0]

fitab = curve_fit(fN,er[mN],et[mN],sigma=np.ones(sum(mN)),absolute_sigma=True)
a  = fitab[0][0]
b  = fitab[0][1]
