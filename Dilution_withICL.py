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
from matplotlib.ticker import ScalarFormatter
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)
from matplotlib import rc
from matplotlib.ticker import NullFormatter
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


matplotlib.rcParams.update({'font.size': 13})

cbcg = 'teal'

def D(DM200,ax1,ax2,ax3):
    
    ax1.axvspan(0, 0.07, alpha=0.3, color=cbcg)
    ax2.axvspan(0, 0.07, alpha=0.3, color=cbcg)
    ax3.axvspan(0, 0.07, alpha=0.3, color=cbcg)


    D_gx200  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1])
    D_gx500 = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1])
    D_gx1000 = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1])
    
    Dgx     = np.array([D_gx1000.D,D_gx500.D,D_gx200.D])
    Dgx_std = np.array([D_gx1000.Dstd,D_gx500.D,D_gx200.Dstd])
        
    #################
    
    D_gxeRM_gap  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_gxcRM_gap  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_gx200RM_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_gx1000RM_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_gx500RM_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    
    Dgx_RMgap = np.array([D_gx1000RM_gap.D,D_gx500RM_gap.D,D_gx200RM_gap.D])
    Dgx_RMgap_std = np.array([D_gx1000RM_gap.Dstd,D_gx500RM_gap.Dstd,D_gx200RM_gap.Dstd])
    
    D_gxeUM_gap  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_gxcUM_gap  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_gx200UM_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_gx1000UM_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_gx500UM_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    
    Dgx_UMgap = np.array([D_gx1000UM_gap.D,D_gx500UM_gap.D,D_gx200UM_gap.D])
    Dgx_UMgap_std = np.array([D_gx1000UM_gap.Dstd,D_gx500UM_gap.Dstd,D_gx200UM_gap.Dstd])
    
    D_gxeRM_off  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_gxcRM_off  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_gx200RM_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_gx1000RM_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_gx500RM_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    
    Dgx_RMoff = np.array([D_gx1000RM_off.D,D_gx500RM_off.D,D_gx200RM_off.D])
    Dgx_RMoff_std = np.array([D_gx1000RM_off.Dstd,D_gx500RM_off.Dstd,D_gx200RM_off.Dstd])
    
    D_gxeUM_off  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_gxcUM_off  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_gx200UM_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_gx1000UM_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_gx500UM_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    
    Dgx_UMoff = np.array([D_gx1000UM_off.D,D_gx500UM_off.D,D_gx200UM_off.D])
    Dgx_UMoff_std = np.array([D_gx1000UM_off.Dstd,D_gx500UM_off.Dstd,D_gx200UM_off.Dstd])
    
    D_gxeRM_dv  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_gxcRM_dv  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_gx200RM_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_gx1000RM_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_gx500RM_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    
    Dgx_RMdv = np.array([D_gx1000RM_dv.D,D_gx500RM_dv.D,D_gx200RM_dv.D])
    Dgx_RMdv_std = np.array([D_gx1000RM_dv.Dstd,D_gx500RM_dv.Dstd,D_gx200RM_dv.Dstd])
    
    D_gxeUM_dv  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_gxcUM_dv  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_gx200UM_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_gx1000UM_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_gx500UM_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    
    Dgx_UMdv = np.array([D_gx1000UM_dv.D,D_gx500UM_dv.D,D_gx200UM_dv.D])
    Dgx_UMdv_std = np.array([D_gx1000UM_dv.Dstd,D_gx500UM_dv.Dstd,D_gx200UM_dv.Dstd])
    
    #######################
    
    D_gxeRL_gap  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_gxcRL_gap  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_gx200RL_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_gx1000RL_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_gx500RL_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    
    Dgx_RLgap = np.array([D_gx1000RL_gap.D,D_gx500RL_gap.D,D_gx200RL_gap.D])
    Dgx_RLgap_std = np.array([D_gx1000RL_gap.Dstd,D_gx500RL_gap.Dstd,D_gx200RL_gap.Dstd])
    
    D_gxeUL_gap  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_gxcUL_gap  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_gx200UL_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_gx1000UL_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_gx500UL_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    
    Dgx_ULgap = np.array([D_gx1000UL_gap.D,D_gx500UL_gap.D,D_gx200UL_gap.D])
    Dgx_ULgap_std = np.array([D_gx1000UL_gap.Dstd,D_gx500UL_gap.Dstd,D_gx200UL_gap.Dstd])
    
    D_gxeRL_off  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_gxcRL_off  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_gx200RL_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_gx1000RL_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_gx500RL_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    
    Dgx_RLoff = np.array([D_gx1000RL_off.D,D_gx500RL_off.D,D_gx200RL_off.D])
    Dgx_RLoff_std = np.array([D_gx1000RL_off.Dstd,D_gx500RL_off.Dstd,D_gx200RL_off.Dstd])
    
    D_gxeUL_off  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_gxcUL_off  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_gx200UL_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_gx1000UL_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_gx500UL_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    
    Dgx_ULoff = np.array([D_gx1000UL_off.D,D_gx500UL_off.D,D_gx200UL_off.D])
    Dgx_ULoff_std = np.array([D_gx1000UL_off.Dstd,D_gx500UL_off.Dstd,D_gx200UL_off.Dstd])
    
    D_gxeRL_dv  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_gxcRL_dv  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_gx200RL_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_gx1000RL_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_gx500RL_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    
    Dgx_RLdv = np.array([D_gx1000RL_dv.D,D_gx500RL_dv.D,D_gx200RL_dv.D])
    Dgx_RLdv_std = np.array([D_gx1000RL_dv.Dstd,D_gx500RL_dv.Dstd,D_gx200RL_dv.Dstd])
    
    D_gxeUL_dv  = Dcompute(cosangle(gxe.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_gxcUL_dv  = Dcompute(cosangle(gxc.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_gx200UL_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_gx1000UL_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_gx500UL_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    
    Dgx_ULdv = np.array([D_gx1000UL_dv.D,D_gx500UL_dv.D,D_gx200UL_dv.D])
    Dgx_ULdv_std = np.array([D_gx1000UL_dv.Dstd,D_gx500UL_dv.Dstd,D_gx200UL_dv.Dstd])
    
    
    
    ########################
    
    ax1.plot(Rgx,Dgx,'k')
    
    ax1.axvline(Rgx_RMgap[-1]+0.03,c='k',lw=0.3)
    ax1.axvline(Rgx_RMgap[-1]+0.16,c='k',lw=0.3)
    ax1.plot(Rgx_RMgap[-1]+0.23,D_gxeRM_gap.D,'^',color='sienna',markersize=8)
    ax1.plot(Rgx_UMgap[-1]+0.23,D_gxeUM_gap.D,'C0^',markersize=8)
    ax1.plot(Rgx_RMgap[-1]+0.1,D_gxcRM_gap.D,'^',color='sienna',markersize=8)
    ax1.plot(Rgx_UMgap[-1]+0.1,D_gxcUM_gap.D,'C0^',markersize=8)

    ax1.plot(Rgx_RLgap[-1]+0.23,D_gxeRL_gap.D,'v',color='sienna',markersize=8)
    ax1.plot(Rgx_ULgap[-1]+0.23,D_gxeUL_gap.D,'C0v',markersize=8)
    ax1.plot(Rgx_RLgap[-1]+0.1,D_gxcRL_gap.D,'v',color='sienna',markersize=8)
    ax1.plot(Rgx_ULgap[-1]+0.1,D_gxcUL_gap.D,'C0v',markersize=8)


    
    ax1.plot(Rgx_RMgap,Dgx_RMgap,'^',color='sienna',markersize=8)
    # ax1.plot(Rgx_RMoff,Dgx_RMoff,'s',color='sienna',markersize=8)
    # ax1.plot(Rgx_RMdv ,Dgx_RMdv ,'^',color='sienna',markersize=8)                        
    ax1.plot(Rgx_UMgap,Dgx_UMgap,'C0^',markersize=8)
    # ax1.plot(Rgx_UMoff,Dgx_UMoff,'sC0')
    # ax1.plot(Rgx_UMdv ,Dgx_UMdv ,'C0^',markersize=8)
    
    ax1.plot(Rgx_RMgap,Dgx_RMgap,'sienna',alpha=1)
    # ax1.plot(Rgx_RMoff,Dgx_RMoff,'sienna',alpha=0.75)
    # ax1.plot(Rgx_RMdv ,Dgx_RMdv ,'sienna',alpha=0.5)                        
    ax1.plot(Rgx_UMgap,Dgx_UMgap,'C0',alpha=1)
    # ax1.plot(Rgx_UMoff,Dgx_UMoff,'C0',alpha=0.75)
    # ax1.plot(Rgx_UMdv ,Dgx_UMdv ,'C0',alpha=0.5) 
    
    ax1.plot(Rgx_RLgap,Dgx_RLgap,'v',color='sienna',markersize=8)
    # ax1.plot(Rgx_RLoff,Dgx_RLoff,'s',color='sienna',markersize=8)
    # ax1.plot(Rgx_RLdv ,Dgx_RLdv ,'^',color='sienna',markersize=8)
    ax1.plot(Rgx_ULgap,Dgx_ULgap,'C0v',markersize=8)
    # ax1.plot(Rgx_ULoff,Dgx_ULoff,'sC0')
    # ax1.plot(Rgx_ULdv ,Dgx_ULdv ,'C0^',markersize=8)
    
    ax1.plot(Rgx_RLgap,Dgx_RLgap,'--',color='sienna',markersize=8)
    # ax1.plot(Rgx_RLoff,Dgx_RLoff,'--',color='sienna',alpha=0.75)
    # ax1.plot(Rgx_RLdv ,Dgx_RLdv ,'--',color='sienna',alpha=0.5)  
    ax1.plot(Rgx_ULgap,Dgx_ULgap,'C0--')
    # ax1.plot(Rgx_ULoff,Dgx_ULoff,'C0--',alpha=0.75)
    # ax1.plot(Rgx_ULdv ,Dgx_ULdv ,'C0--',alpha=0.5) 
    
    
    ########### ESTRELLAS
    
    D_st30  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1])
    D_st50  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1])
    D_st100  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1])
    D_st200  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1])
    D_st500  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1])
    D_st1000 = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1])
    
    Dst = np.array([D_st30.D,D_st50.D,D_st100.D,D_st1000.D,D_st500.D,D_st200.D])
    
    
    #################
    
    D_st30RM_gap   = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_st50RM_gap   = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_st100RM_gap  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_st200RM_gap  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_st1000RM_gap = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_st500RM_gap  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    
    Dst_RMgap = np.array([D_st30RM_gap.D,D_st50RM_gap.D,D_st100RM_gap.D,D_st1000RM_gap.D,D_st500RM_gap.D,D_st200RM_gap.D])
    
    D_st30UM_gap  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_st50UM_gap = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_st100UM_gap  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_st200UM_gap  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_st1000UM_gap = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_st500UM_gap  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    
    Dst_UMgap = np.array([D_st30UM_gap.D,D_st50UM_gap.D,D_st100UM_gap.D,D_st1000UM_gap.D,D_st500UM_gap.D,D_st200UM_gap.D])
    
    D_st30RM_off   = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_st50RM_off   = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_st100RM_off  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_st200RM_off  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_st1000RM_off = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_st500RM_off  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    
    Dst_RMoff = np.array([D_st30RM_off.D,D_st50RM_off.D,D_st100RM_off.D,D_st1000RM_off.D,D_st500RM_off.D,D_st200RM_off.D])
    
    D_st30UM_off  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_st50UM_off  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_st100UM_off  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_st200UM_off  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_st1000UM_off = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_st500UM_off  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    
    Dst_UMoff = np.array([D_st30UM_off.D,D_st50UM_off.D,D_st100UM_off.D,D_st1000UM_off.D,D_st500UM_off.D,D_st200UM_off.D])
    
    D_st30RM_dv  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_st50RM_dv  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_st100RM_dv  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_st200RM_dv  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_st1000RM_dv = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_st500RM_dv  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    
    Dst_RMdv = np.array([D_st30RM_dv.D,D_st50RM_dv.D,D_st100RM_dv.D,D_st1000RM_dv.D,D_st500RM_dv.D,D_st200RM_dv.D])
    
    D_st30UM_dv  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_st50UM_dv  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_st100UM_dv  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_st200UM_dv  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_st1000UM_dv = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_st500UM_dv  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    
    Dst_UMdv = np.array([D_st30UM_dv.D,D_st50UM_dv.D,D_st100UM_dv.D,D_st1000UM_dv.D,D_st500UM_dv.D,D_st200UM_dv.D])
    
    #######################
    
    D_st30RL_gap   = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_st50RL_gap   = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_st100RL_gap  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_st200RL_gap  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_st1000RL_gap = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_st500RL_gap  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    
    Dst_RLgap = np.array([D_st30RL_gap.D,D_st50RL_gap.D,D_st100RL_gap.D,D_st1000RL_gap.D,D_st500RL_gap.D,D_st200RL_gap.D])
    
    D_st30UL_gap  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_st50UL_gap = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_st100UL_gap  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_st200UL_gap  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_st1000UL_gap = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_st500UL_gap  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    
    Dst_ULgap = np.array([D_st30UL_gap.D,D_st50UL_gap.D,D_st100UL_gap.D,D_st1000UL_gap.D,D_st500UL_gap.D,D_st200UL_gap.D])
    
    D_st30RL_off   = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_st50RL_off   = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_st100RL_off  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_st200RL_off  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_st1000RL_off = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_st500RL_off  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    
    Dst_RLoff = np.array([D_st30RL_off.D,D_st50RL_off.D,D_st100RL_off.D,D_st1000RL_off.D,D_st500RL_off.D,D_st200RL_off.D])
    
    D_st30UL_off  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_st50UL_off  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_st100UL_off  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_st200UL_off  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_st1000UL_off = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_st500UL_off  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    
    Dst_ULoff = np.array([D_st30UL_off.D,D_st50UL_off.D,D_st100UL_off.D,D_st1000UL_off.D,D_st500UL_off.D,D_st200UL_off.D])
    
    D_st30RL_dv  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_st50RL_dv  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_st100RL_dv  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_st200RL_dv  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_st1000RL_dv = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_st500RL_dv  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    
    Dst_RLdv = np.array([D_st30RL_dv.D,D_st50RL_dv.D,D_st100RL_dv.D,D_st1000RL_dv.D,D_st500RL_dv.D,D_st200RL_dv.D])
    
    D_st30UL_dv  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_st50UL_dv  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_st100UL_dv  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_st200UL_dv  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_st1000UL_dv = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_st500UL_dv  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    
    Dst_ULdv = np.array([D_st30UL_dv.D,D_st50UL_dv.D,D_st100UL_dv.D,D_st1000UL_dv.D,D_st500UL_dv.D,D_st200UL_dv.D])
    
    ############
    
    ax2.plot(Rst,Dst,'k')
    
    
    ax2.plot(Rst_RMgap,Dst_RMgap,'^',color='sienna',markersize=8)
    # ax2.plot(Rst_RMoff,Dst_RMoff,'s',color='sienna',markersize=8)
    # ax2.plot(Rst_RMdv ,Dst_RMdv ,'^',color='sienna',markersize=8)
    ax2.plot(Rst_UMgap,Dst_UMgap,'C0^',markersize=8)
    # ax2.plot(Rst_UMoff,Dst_UMoff,'sC0')
    # ax2.plot(Rst_UMdv ,Dst_UMdv ,'C0^',markersize=8)
    
    ax2.plot(Rst_RMgap,Dst_RMgap,'sienna',alpha=1)
    # ax2.plot(Rst_RMoff,Dst_RMoff,'sienna',alpha=0.75)
    # ax2.plot(Rst_RMdv ,Dst_RMdv ,'sienna',alpha=0.5)
    ax2.plot(Rst_UMgap,Dst_UMgap,'C0',alpha=1)
    # ax2.plot(Rst_UMoff,Dst_UMoff,'C0',alpha=0.75)
    # ax2.plot(Rst_UMdv ,Dst_UMdv ,'C0',alpha=0.5)
    
    ax2.plot(Rst_RLgap,Dst_RLgap,'v',color='sienna',markersize=8)
    # ax2.plot(Rst_RLoff,Dst_RLoff,'s',color='sienna',markersize=8)
    # ax2.plot(Rst_RLdv ,Dst_RLdv ,'^',color='sienna',markersize=8)
    ax2.plot(Rst_ULgap,Dst_ULgap,'C0v',markersize=8)
    # ax2.plot(Rst_ULoff,Dst_ULoff,'sC0')
    # ax2.plot(Rst_ULdv ,Dst_ULdv ,'C0^',markersize=8)
    
    ax2.plot(Rst_RLgap,Dst_RLgap,'--',color='sienna',alpha=1)
    # ax2.plot(Rst_RLoff,Dst_RLoff,'--',color='sienna',alpha=0.75)
    # ax2.plot(Rst_RLdv ,Dst_RLdv ,'--',color='sienna',alpha=0.5)
    ax2.plot(Rst_ULgap,Dst_ULgap,'C0--',alpha=1)
    # ax2.plot(Rst_ULoff,Dst_ULoff,'C0--',alpha=0.75)
    # ax2.plot(Rst_ULdv ,Dst_ULdv ,'C0--',alpha=0.5)

    ax3.plot(100,100,'k^',label='$\log(M_{200}/M_\odot) > 14.6$',markersize=8)
    ax3.plot(100,100,'kv',label='$\log(M_{200}/M_\odot) < 14.6$',markersize=8)
    
    ##############################
    ## ICL
    ##############################
    
    tbcg = cosangle(ICL_.a2Dbcg[mmas[micl2]][mbcg],DM200.a2D[micl][mbcg])[1]
    ticl = cosangle(ICL_.a2Dicl[mmas[micl2]],DM200.a2D[micl])[1]
    # tbcg = cosangle(ICL_.a2Dbcg,DM200.a2D[micl])[1]
    # ticl = cosangle(ICL_.a2Dicl,DM200.a2D[micl])[1]
    
    D_bcg   = Dcompute(tbcg).D
    D_icl   = Dcompute(ticl).D    
    
    D_bcgR_gap   = Dcompute(tbcg[(C.mo2d_gap)[micl][mbcg]]).D
    D_iclR_gap   = Dcompute(ticl[(C.mo2d_gap)[micl]]).D
    D_bcgU_gap   = Dcompute(tbcg[(C.mn2d_gap)[micl][mbcg]]).D
    D_iclU_gap   = Dcompute(ticl[(C.mn2d_gap)[micl]]).D

    D_bcgR_dv   = Dcompute(tbcg[(C.mo2d_dv)[micl][mbcg]]).D
    D_iclR_dv   = Dcompute(ticl[(C.mo2d_dv)[micl]]).D
    D_bcgU_dv   = Dcompute(tbcg[(C.mn2d_dv)[micl][mbcg]]).D
    D_iclU_dv   = Dcompute(ticl[(C.mn2d_dv)[micl]]).D

    D_bcgR_off   = Dcompute(tbcg[(C.mo2d_off)[micl][mbcg]]).D
    D_iclR_off   = Dcompute(ticl[(C.mo2d_off)[micl]]).D
    D_bcgU_off   = Dcompute(tbcg[(C.mn2d_off)[micl][mbcg]]).D
    D_iclU_off   = Dcompute(ticl[(C.mn2d_off)[micl]]).D
    
    ax3.plot([Rbcg,Ricl],[D_bcg,D_icl],'k')
    
    
    # ax3.plot(Ricl_Rgap,D_iclR_gap,'o',color='sienna',markersize=8)
    # ax3.plot(Ricl_Roff,D_iclR_off,'s',color='sienna',markersize=8)
    # ax3.plot(Ricl_Rdv ,D_iclR_dv ,'^',color='sienna',markersize=8)
    # ax3.plot(Ricl_Ugap,D_iclU_gap,'oC0')
    # ax3.plot(Ricl_Uoff,D_iclU_off,'sC0')
    # ax3.plot(Ricl_Udv ,D_iclU_dv ,'C0^',markersize=8)

    # ax3.plot(Rbcg_Rgap,D_bcgR_gap,'o',color='sienna',markersize=8)
    # ax3.plot(Rbcg_Roff,D_bcgR_off,'s',color='sienna',markersize=8)
    # ax3.plot(Rbcg_Rdv ,D_bcgR_dv ,'^',color='sienna',markersize=8)
    # ax3.plot(Rbcg_Ugap,D_bcgU_gap,'oC0')
    # ax3.plot(Rbcg_Uoff,D_bcgU_off,'sC0')
    # ax3.plot(Rbcg_Udv ,D_bcgU_dv ,'C0^',markersize=8)
    

    ax3.plot([Rbcg_Rgap,Ricl_Rgap],[D_bcgR_gap,D_iclR_gap],'^',color='sienna',markersize=8)
    # ax3.plot([Rbcg_Roff,Ricl_Roff],[D_bcgR_off,D_iclR_off],'s',color='sienna',markersize=8)
    # ax3.plot([Rbcg_Rdv ,Ricl_Rdv ],[D_bcgR_dv ,D_iclR_dv ],'^',color='sienna',markersize=8)
    ax3.plot([Rbcg_Ugap,Ricl_Ugap],[D_bcgU_gap,D_iclU_gap],'C0^',markersize=8)
    # ax3.plot([Rbcg_Uoff,Ricl_Uoff],[D_bcgU_off,D_iclU_off],'sC0')
    # ax3.plot([Rbcg_Udv ,Ricl_Udv ],[D_bcgU_dv ,D_iclU_dv ],'C0^',markersize=8)

    ax3.plot([Rbcg_Rgap,Ricl_Rgap],[D_bcgR_gap,D_iclR_gap],'',color='sienna',markersize=8)
    # ax3.plot([Rbcg_Roff,Ricl_Roff],[D_bcgR_off,D_iclR_off],'',color='sienna',markersize=8)
    # ax3.plot([Rbcg_Rdv ,Ricl_Rdv ],[D_bcgR_dv ,D_iclR_dv ],'',color='sienna',markersize=8)
    ax3.plot([Rbcg_Ugap,Ricl_Ugap],[D_bcgU_gap,D_iclU_gap],'C0')
    # ax3.plot([Rbcg_Uoff,Ricl_Uoff],[D_bcgU_off,D_iclU_off],'C0')
    # ax3.plot([Rbcg_Udv ,Ricl_Udv ],[D_bcgU_dv ,D_iclU_dv ],'C0')

    

    # ax2.plot(100,100,'k-',label='$\log(M_{200}/M_\odot) > 14.6$')
    # ax2.plot(100,100,'k--',label='$\log(M_{200}/M_\odot) < 14.6$')





plotspath = '../final_plots/'

C = Clusters()
ICL_ = ICL()

DM1000 = DarkMatter(1000)
DM500  = DarkMatter(500)
DM200  = DarkMatter(200)

r200 = Random()
r1000 = Random(1000)
r500 = Random(500)

gx200 = Galaxias(radio=200)
gx1000 = Galaxias(radio=1000)
gx500 = Galaxias(radio=500)

gxe = Galaxias('extendidas')
gxc = Galaxias('concentradas')


m = (C.sub == 0)

mc     = gxc.N > 9
mcp    = np.array((mc.tolist())*3)


# ltimep = C.ltimep
# mt     = ltimep > 0.
# C.mo2d_gap[mt] = ltimep[mt] > 5.5
# C.mn2d_gap[mt] = ltimep[mt] <= 5.5

mmas = (C.lM200p > 14.6)#*mt
mlow = (C.lM200p < 14.6)#*mt

micl2 = (m.tolist()*3)
micl = (m.tolist()*3)*mmas
mbcg = (ICL_.Rbcg[mmas[micl2]] > 0)

yaxis = np.zeros((len(DM200.a2D),2))
yaxis[:,1] = 1.



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

gx = gx200
r  =  r200
DM = DM200



Rst = np.median(C.Rsp,axis=0)
Rgx = np.median(C.Rsp,axis=0)[3:]

Rgx_RMgap = np.median(C.Rsp[mmas*C.mo2d_gap],axis=0)[3:]
Rgx_UMgap = np.median(C.Rsp[mmas*C.mn2d_gap],axis=0)[3:]
Rgx_RMoff = np.median(C.Rsp[mmas*C.mo2d_off],axis=0)[3:]
Rgx_UMoff = np.median(C.Rsp[mmas*C.mn2d_off],axis=0)[3:]
Rgx_RMdv  = np.median(C.Rsp[mmas*C.mo2d_dv],axis=0)[3:]
Rgx_UMdv  = np.median(C.Rsp[mmas*C.mn2d_dv],axis=0)[3:]
Rgx_RLgap = np.median(C.Rsp[mlow*C.mo2d_gap],axis=0)[3:]
Rgx_ULgap = np.median(C.Rsp[mlow*C.mn2d_gap],axis=0)[3:]
Rgx_RLoff = np.median(C.Rsp[mlow*C.mo2d_off],axis=0)[3:]
Rgx_ULoff = np.median(C.Rsp[mlow*C.mn2d_off],axis=0)[3:]
Rgx_RLdv  = np.median(C.Rsp[mlow*C.mo2d_dv],axis=0)[3:]
Rgx_ULdv  = np.median(C.Rsp[mlow*C.mn2d_dv],axis=0)[3:]

Rst_RMgap = np.median(C.Rsp[mmas*C.mo2d_gap],axis=0)
Rst_UMgap = np.median(C.Rsp[mmas*C.mn2d_gap],axis=0)
Rst_RMoff = np.median(C.Rsp[mmas*C.mo2d_off],axis=0)
Rst_UMoff = np.median(C.Rsp[mmas*C.mn2d_off],axis=0)
Rst_RMdv  = np.median(C.Rsp[mmas*C.mo2d_dv],axis=0)
Rst_UMdv  = np.median(C.Rsp[mmas*C.mn2d_dv],axis=0)
Rst_RLgap = np.median(C.Rsp[mlow*C.mo2d_gap],axis=0)
Rst_ULgap = np.median(C.Rsp[mlow*C.mn2d_gap],axis=0)
Rst_RLoff = np.median(C.Rsp[mlow*C.mo2d_off],axis=0)
Rst_ULoff = np.median(C.Rsp[mlow*C.mn2d_off],axis=0)
Rst_RLdv  = np.median(C.Rsp[mlow*C.mo2d_dv],axis=0)
Rst_ULdv  = np.median(C.Rsp[mlow*C.mn2d_dv],axis=0)


Rbcg      = np.median(ICL_.Rbcg[mmas[micl2]][mbcg])
Rbcg_Rgap = np.median(ICL_.Rbcg[mmas[micl2]][(C.mo2d_gap)[micl]*mbcg])
Rbcg_Ugap = np.median(ICL_.Rbcg[mmas[micl2]][(C.mn2d_gap)[micl]*mbcg])
Rbcg_Roff = np.median(ICL_.Rbcg[mmas[micl2]][(C.mo2d_off)[micl]*mbcg])
Rbcg_Uoff = np.median(ICL_.Rbcg[mmas[micl2]][(C.mn2d_off)[micl]*mbcg])
Rbcg_Rdv  = np.median(ICL_.Rbcg[mmas[micl2]][(C.mo2d_dv)[micl]*mbcg])
Rbcg_Udv  = np.median(ICL_.Rbcg[mmas[micl2]][(C.mn2d_dv)[micl]*mbcg])

Ricl      = np.median(ICL_.Ricl[mmas[micl2]])
Ricl_Rgap = np.median(ICL_.Ricl[mmas[micl2]][(C.mo2d_gap)[micl]])
Ricl_Ugap = np.median(ICL_.Ricl[mmas[micl2]][(C.mn2d_gap)[micl]])
Ricl_Roff = np.median(ICL_.Ricl[mmas[micl2]][(C.mo2d_off)[micl]])
Ricl_Uoff = np.median(ICL_.Ricl[mmas[micl2]][(C.mn2d_off)[micl]])
Ricl_Rdv  = np.median(ICL_.Ricl[mmas[micl2]][(C.mo2d_dv)[micl]])
Ricl_Udv  = np.median(ICL_.Ricl[mmas[micl2]][(C.mn2d_dv)[micl]])

# Ricl       = 0.1*np.median(C.Rsp[micl],axis=0)[4:]
# Ricl_Rgap = 0.1*np.median(C.Rsp[micl*C.mo2d_gap],axis=0)[4:]
# Ricl_Ugap = 0.1*np.median(C.Rsp[micl*C.mn2d_gap],axis=0)[4:]
# Ricl_Roff = 0.1*np.median(C.Rsp[micl*C.mo2d_off],axis=0)[4:]
# Ricl_Uoff = 0.1*np.median(C.Rsp[micl*C.mn2d_off],axis=0)[4:]
# Ricl_Rdv  = 0.1*np.median(C.Rsp[micl*C.mo2d_dv],axis=0)[4:]
# Ricl_Udv  = 0.1*np.median(C.Rsp[micl*C.mn2d_dv],axis=0)[4:]



####### PLOT
f, ax = plt.subplots(3,3, figsize=(12,10), sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0,0].text(1.08,1.,'c')
ax[0,0].text(1.21,1.,'e')
ax[0,0].text(0.6,1.,'Galaxies')
ax[0,1].text(0.1,1.,'Stars')
ax[0,2].text(0.13,1.,'Isophotes')

ax[0,1].text(0.02,1.,'BCG',c=cbcg)
ax[0,2].text(0.052,1.,'BCG',c=cbcg)

ax[0,0].text(0.45,0.3,r'$\theta^*_{DM}(R < R_{200})$')
ax[1,0].text(0.45,0.3,r'$\theta^*_{DM}(R < R_{500})$')
ax[2,0].text(0.45,0.3,r'$\theta^*_{DM}(R < R_{1000})$')



ax[0,0].set_xscale('log')
ax[1,0].set_xscale('log')
ax[2,0].set_xscale('log')

ax[0,0].set_xlim([0.4,1.3])
ax[1,0].set_xlim([0.4,1.3])
ax[2,0].set_xlim([0.4,1.3])

ax[0,1].set_xscale('log')
ax[1,1].set_xscale('log')
ax[2,1].set_xscale('log')

ax[0,1].set_xlim([0.0105,1.5])
ax[1,1].set_xlim([0.0105,1.5])
ax[2,1].set_xlim([0.0105,1.5])

ax[0,2].set_xscale('log')
ax[1,2].set_xscale('log')
ax[2,2].set_xscale('log')

ax[0,2].set_xlim([0.05,0.5])
ax[1,2].set_xlim([0.05,0.5])
ax[2,2].set_xlim([0.05,0.5])


for j in np.arange(3):

    ax[j,2].xaxis.set_major_formatter(NullFormatter())
    ax[j,2].xaxis.set_minor_formatter(NullFormatter())
    ax[j,0].xaxis.set_major_formatter(NullFormatter())
    ax[j,0].xaxis.set_minor_formatter(NullFormatter())

ax[0,0].set_xticks([0.5,0.7,1.0])
ax[1,0].set_xticks([0.5,0.7,1.0])
ax[2,0].set_xticks([0.5,0.7,1.0])
ax[0,0].set_xticklabels([])
ax[1,0].set_xticklabels([])
ax[2,0].set_xticklabels(['0.5','0.7','1.0'])

ax[0,1].set_xticks([0.1,0.3,0.5,1.0])
ax[1,1].set_xticks([0.1,0.3,0.5,1.0])
ax[2,1].set_xticks([0.1,0.3,0.5,1.0])
ax[0,1].set_xticklabels([])
ax[1,1].set_xticklabels([])
ax[2,1].set_xticklabels(['0.1','0.3','0.5','1.0'])

ax[0,2].set_xticks([0.1,0.2,0.3,0.5])
ax[1,2].set_xticks([0.1,0.2,0.3,0.5])
ax[2,2].set_xticks([0.1,0.2,0.3,0.5])
ax[0,2].set_xticklabels([])
ax[1,2].set_xticklabels([])
ax[2,2].set_xticklabels(['0.1','0.2','0.3','0.5'])


ax[2,0].set_xlabel('$R/R_{200}$')
ax[2,1].set_xlabel('$R/R_{200}$')
ax[2,2].set_xlabel('$R/R_{200}$')

ax[0,0].set_ylabel('$D$')
ax[1,0].set_ylabel('$D$')
ax[2,0].set_ylabel('$D$')

ax[0,0].set_ylim([0.1,1.1])


# ax[2,2].plot(100,100,'ko',label=r'$M_{sat}/M_{BCG}$')
# ax[2,2].plot(100,100,'ks',label=r'$D_{offset}$')
# ax[2,2].plot(100,100,'k^',label=r'$\Delta V$')


'''
ax2 = ax[0,0].twiny()
ax2.set_xlim([0.4,1.1])
ax2.set_xscale('log')
ax2.xaxis.set_major_formatter(NullFormatter())
ax2.xaxis.set_minor_formatter(NullFormatter())
ax2.set_xticks(Rgx)
ax2.set_xticklabels(['$R_{1000}$','$R_{500}$','$R_{200}$'])

ax2 = ax[0,1].twiny()
ax2.set_xlim([0.0105,1.5])
ax2.set_xscale('log')
ax2.xaxis.set_major_formatter(NullFormatter())
ax2.xaxis.set_minor_formatter(NullFormatter())
ax2.set_xticks(Rst)
ax2.set_xticklabels(['30kpc','50kpc','$0.1R_{500}$','$R_{1000}$','$R_{500}$','$R_{200}$'])
'''


D(DM200,ax[0,0],ax[0,1],ax[0,2])
D(DM500,ax[1,0],ax[1,1],ax[1,2])
D(DM1000,ax[2,0],ax[2,1],ax[2,2])


ax[2,2].plot(100,100,'sienna',label=r'relaxed')
ax[2,2].plot(100,100,'C0',label=r'non-relaxed')
ax[2,2].legend(frameon=False,loc=4)


plt.savefig(plotspath+'D_R_withicl.pdf',bbox_inches='tight')





'''
mlist = [mmas*C.mo2d_gap,mmas*C.mn2d_gap,
         mlow*C.mo2d_gap,mlow*C.mn2d_gap]

# mlist = [mmas*C.mo2d_dv,mmas*C.mn2d_dv,
         # mlow*C.mo2d_dv,mlow*C.mn2d_dv]

# mlist = [mmas*C.mo2d_off,mmas*C.mn2d_off,
         # mlow*C.mo2d_off,mlow*C.mn2d_off]

samples = ['$\log(M_{200}/M_\odot) > 14.6$','$\log(M_{200}/M_\odot) > 14.6$',
            '$\log(M_{200}/M_\odot) < 14.6$','$\log(M_{200}/M_\odot) < 14.6$']

samples2 = ['relaxed','non-relaxed',
            'relaxed','non-relaxed']

f, ax = plt.subplots(3,4, figsize=(12,8), sharey=True,sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0,1].set_xlim([-0.15,1.05])

ax[0,0].set_ylabel('$n$')
ax[1,0].set_ylabel('$n$')
ax[2,0].set_ylabel('$n$')

ax[2,0].set_xlabel(r'$\cos(2 \Delta \theta)$')
ax[2,1].set_xlabel(r'$\cos(2 \Delta \theta)$')
ax[2,2].set_xlabel(r'$\cos(2 \Delta \theta)$')
ax[2,3].set_xlabel(r'$\cos(2 \Delta \theta)$')

ax[0,0].text(0.05,8,r'$\theta^*_{DM}(R < R_{200})$',fontsize=12)
ax[1,0].text(0.05,8,r'$\theta^*_{DM}(R < R_{500})$',fontsize=12)
ax[2,0].text(0.05,8,r'$\theta^*_{DM}(R < R_{1000})$',fontsize=12)

ax2 = ax.flatten()

for j in np.arange(4):
    
    c2te = np.cos(2.*(np.deg2rad(cosangle(gxe.a2D,DM200.a2D)[1][mlist[j]])))
    c2tc = np.cos(2.*(np.deg2rad(cosangle(gxc.a2D,DM200.a2D)[1][mlist[j]])))
    
    c2te = c2te[c2te > 0.]
    c2tc = c2tc[c2tc > 0.]
    ax2[j].hist(c2te,np.linspace(0,1,15),histtype='step',label=samples[j],density=True,color=cbcg)
    ax2[j].hist(c2tc,np.linspace(0,1,15),histtype='step',density=True,color=cbcg)
    
    ax2[j].text(0.05,6,samples[j],fontsize=11)
    ax2[j].text(0.05,5,samples2[j],fontsize=11)

    ax2[j].plot(np.mean(c2tc),10,'C1o')    
    ax2[j].plot(np.mean(c2te),10,'C4o')

    
    ax2[j].errorbar(np.mean(c2te),10,xerr=np.std(c2te),ecolor=cbcg,fmt='none')
    ax2[j].errorbar(np.mean(c2tc),10,xerr=np.std(c2tc),ecolor=cbcg,fmt='none')

for j in np.arange(4):
    
    c2te = np.cos(2.*(np.deg2rad(cosangle(gxe.a2D,DM500.a2D)[1][mlist[j]])))
    c2tc = np.cos(2.*(np.deg2rad(cosangle(gxc.a2D,DM500.a2D)[1][mlist[j]])))
    
    c2te = c2te[c2te > 0.]
    c2tc = c2tc[c2tc > 0.]
    ax2[j+4].hist(c2te,np.linspace(0,1,15),histtype='step',label=samples[j],density=True,color=cbcg)
    ax2[j+4].hist(c2tc,np.linspace(0,1,15),histtype='step',density=True,color=cbcg)
    
    # ax2[j].text(0.1,6,samples[j])
    # ax2[j].text(0.1,5,samples2[j])

    ax2[j+4].plot(np.mean(c2tc),10,'C1o')    
    ax2[j+4].plot(np.mean(c2te),10,'C4o')

    
    ax2[j+4].errorbar(np.mean(c2te),10,xerr=np.std(c2te),ecolor=cbcg,fmt='none')
    ax2[j+4].errorbar(np.mean(c2tc),10,xerr=np.std(c2tc),ecolor=cbcg,fmt='none')


for j in np.arange(4):
    
    c2te = np.cos(2.*(np.deg2rad(cosangle(gxe.a2D,DM1000.a2D)[1][mlist[j]])))
    c2tc = np.cos(2.*(np.deg2rad(cosangle(gxc.a2D,DM1000.a2D)[1][mlist[j]])))
    
    c2te = c2te[c2te > 0.]
    c2tc = c2tc[c2tc > 0.]
    ax2[j+8].hist(c2te,np.linspace(0,1,15),histtype='step',label=samples[j],density=True,color=cbcg)
    ax2[j+8].hist(c2tc,np.linspace(0,1,15),histtype='step',density=True,color=cbcg)
    
    # ax2[j].text(0.1,6,samples[j])
    # ax2[j].text(0.1,5,samples2[j])

    ax2[j+8].plot(np.mean(c2tc),10,'C1o')    
    ax2[j+8].plot(np.mean(c2te),10,'C4o')

    
    ax2[j+8].errorbar(np.mean(c2te),10,xerr=np.std(c2te),ecolor=cbcg,fmt='none')
    ax2[j+8].errorbar(np.mean(c2tc),10,xerr=np.std(c2tc),ecolor=cbcg,fmt='none')

    
    # ax2[j].legend()

plt.savefig(plotspath+'D_ec_gap.pdf',bbox_inches='tight')
'''
