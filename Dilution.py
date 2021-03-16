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


def D(DM200,ax1,ax2):
    

    D_gx200  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1])
    D_gx500 = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1])
    D_gx1000 = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1])
    
    Dgx     = np.array([D_gx1000.D,D_gx500.D,D_gx200.D])
    Dgx_std = np.array([D_gx1000.Dstd,D_gx500.D,D_gx200.Dstd])
    
    D_gx200M  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas])
    D_gx1000M = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas])
    D_gx500M  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas])
    
    Dgx_M = np.array([D_gx1000M.D,D_gx500M.D,D_gx200M.D])
    Dgx_M_std = np.array([D_gx1000M.Dstd,D_gx500M.Dstd,D_gx200M.Dstd])
    
    D_gx200L  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow])
    D_gx1000L = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow])
    D_gx500L  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow])
    
    Dgx_L = np.array([D_gx1000L.D,D_gx500L.D,D_gx200L.D])
    Dgx_L_std = np.array([D_gx1000L.Dstd,D_gx500L.Dstd,D_gx200L.Dstd])
    
    D_gx200R_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][C.mo2d_gap])
    D_gx1000R_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][C.mo2d_gap])
    D_gx500R_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][C.mo2d_gap])
    
    Dgx_Rgap = np.array([D_gx1000R_gap.D,D_gx500R_gap.D,D_gx200R_gap.D])
    Dgx_Rgap_std = np.array([D_gx1000R_gap.Dstd,D_gx500R_gap.Dstd,D_gx200R_gap.Dstd])
    
    D_gx200U_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][C.mn2d_gap])
    D_gx1000U_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][C.mn2d_gap])
    D_gx500U_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][C.mn2d_gap])
    
    Dgx_Ugap = np.array([D_gx1000U_gap.D,D_gx500U_gap.D,D_gx200U_gap.D])
    Dgx_Ugap_std = np.array([D_gx1000U_gap.Dstd,D_gx500U_gap.Dstd,D_gx200U_gap.Dstd])
    
    D_gx200R_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][C.mo2d_off])
    D_gx1000R_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][C.mo2d_off])
    D_gx500R_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][C.mo2d_off])
    
    Dgx_Roff = np.array([D_gx1000R_off.D,D_gx500R_off.D,D_gx200R_off.D])
    Dgx_Roff_std = np.array([D_gx1000R_off.Dstd,D_gx500R_off.Dstd,D_gx200R_off.Dstd])
    
    D_gx200U_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][C.mn2d_off])
    D_gx1000U_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][C.mn2d_off])
    D_gx500U_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][C.mn2d_off])
    
    Dgx_Uoff = np.array([D_gx1000U_off.D,D_gx500U_off.D,D_gx200U_off.D])
    Dgx_Uoff_std = np.array([D_gx1000U_off.Dstd,D_gx500U_off.Dstd,D_gx200U_off.Dstd])
    
    D_gx200R_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][C.mo2d_dv])
    D_gx1000R_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][C.mo2d_dv])
    D_gx500R_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][C.mo2d_dv])
    
    Dgx_Rdv = np.array([D_gx1000R_dv.D,D_gx500R_dv.D,D_gx200R_dv.D])
    Dgx_Rdv_std = np.array([D_gx1000R_dv.Dstd,D_gx500R_dv.Dstd,D_gx200R_dv.Dstd])
    
    D_gx200U_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][C.mn2d_dv])
    D_gx1000U_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][C.mn2d_dv])
    D_gx500U_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][C.mn2d_dv])
    
    Dgx_Udv = np.array([D_gx1000U_dv.D,D_gx500U_dv.D,D_gx200U_dv.D])
    Dgx_Udv_std = np.array([D_gx1000U_dv.Dstd,D_gx500U_dv.Dstd,D_gx200U_dv.Dstd])
    
    #################
    
    D_gx200RM_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_gx1000RM_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    D_gx500RM_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
    
    Dgx_RMgap = np.array([D_gx1000RM_gap.D,D_gx500RM_gap.D,D_gx200RM_gap.D])
    Dgx_RMgap_std = np.array([D_gx1000RM_gap.Dstd,D_gx500RM_gap.Dstd,D_gx200RM_gap.Dstd])
    
    D_gx200UM_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_gx1000UM_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    D_gx500UM_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
    
    Dgx_UMgap = np.array([D_gx1000UM_gap.D,D_gx500UM_gap.D,D_gx200UM_gap.D])
    Dgx_UMgap_std = np.array([D_gx1000UM_gap.Dstd,D_gx500UM_gap.Dstd,D_gx200UM_gap.Dstd])
    
    D_gx200RM_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_gx1000RM_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    D_gx500RM_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
    
    Dgx_RMoff = np.array([D_gx1000RM_off.D,D_gx500RM_off.D,D_gx200RM_off.D])
    Dgx_RMoff_std = np.array([D_gx1000RM_off.Dstd,D_gx500RM_off.Dstd,D_gx200RM_off.Dstd])
    
    D_gx200UM_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_gx1000UM_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    D_gx500UM_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
    
    Dgx_UMoff = np.array([D_gx1000UM_off.D,D_gx500UM_off.D,D_gx200UM_off.D])
    Dgx_UMoff_std = np.array([D_gx1000UM_off.Dstd,D_gx500UM_off.Dstd,D_gx200UM_off.Dstd])
    
    D_gx200RM_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_gx1000RM_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    D_gx500RM_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
    
    Dgx_RMdv = np.array([D_gx1000RM_dv.D,D_gx500RM_dv.D,D_gx200RM_dv.D])
    Dgx_RMdv_std = np.array([D_gx1000RM_dv.Dstd,D_gx500RM_dv.Dstd,D_gx200RM_dv.Dstd])
    
    D_gx200UM_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_gx1000UM_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    D_gx500UM_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
    
    Dgx_UMdv = np.array([D_gx1000UM_dv.D,D_gx500UM_dv.D,D_gx200UM_dv.D])
    Dgx_UMdv_std = np.array([D_gx1000UM_dv.Dstd,D_gx500UM_dv.Dstd,D_gx200UM_dv.Dstd])
    
    #######################
    
    D_gx200RL_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_gx1000RL_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    D_gx500RL_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
    
    Dgx_RLgap = np.array([D_gx1000RL_gap.D,D_gx500RL_gap.D,D_gx200RL_gap.D])
    Dgx_RLgap_std = np.array([D_gx1000RL_gap.Dstd,D_gx500RL_gap.Dstd,D_gx200RL_gap.Dstd])
    
    D_gx200UL_gap  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_gx1000UL_gap = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    D_gx500UL_gap  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
    
    Dgx_ULgap = np.array([D_gx1000UL_gap.D,D_gx500UL_gap.D,D_gx200UL_gap.D])
    Dgx_ULgap_std = np.array([D_gx1000UL_gap.Dstd,D_gx500UL_gap.Dstd,D_gx200UL_gap.Dstd])
    
    D_gx200RL_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_gx1000RL_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    D_gx500RL_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
    
    Dgx_RLoff = np.array([D_gx1000RL_off.D,D_gx500RL_off.D,D_gx200RL_off.D])
    Dgx_RLoff_std = np.array([D_gx1000RL_off.Dstd,D_gx500RL_off.Dstd,D_gx200RL_off.Dstd])
    
    D_gx200UL_off  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_gx1000UL_off = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    D_gx500UL_off  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
    
    Dgx_ULoff = np.array([D_gx1000UL_off.D,D_gx500UL_off.D,D_gx200UL_off.D])
    Dgx_ULoff_std = np.array([D_gx1000UL_off.Dstd,D_gx500UL_off.Dstd,D_gx200UL_off.Dstd])
    
    D_gx200RL_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_gx1000RL_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    D_gx500RL_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
    
    Dgx_RLdv = np.array([D_gx1000RL_dv.D,D_gx500RL_dv.D,D_gx200RL_dv.D])
    Dgx_RLdv_std = np.array([D_gx1000RL_dv.Dstd,D_gx500RL_dv.Dstd,D_gx200RL_dv.Dstd])
    
    D_gx200UL_dv  = Dcompute(cosangle(gx200.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_gx1000UL_dv = Dcompute(cosangle(gx1000.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    D_gx500UL_dv  = Dcompute(cosangle(gx500.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
    
    Dgx_ULdv = np.array([D_gx1000UL_dv.D,D_gx500UL_dv.D,D_gx200UL_dv.D])
    Dgx_ULdv_std = np.array([D_gx1000UL_dv.Dstd,D_gx500UL_dv.Dstd,D_gx200UL_dv.Dstd])
    
    ########################
    
    ax1.plot(Rgx,Dgx,'k')
    
    ax1.plot(Rgx_RMgap,Dgx_RMgap,'oC3')
    ax1.plot(Rgx_RMoff,Dgx_RMoff,'sC3')
    ax1.plot(Rgx_RMdv ,Dgx_RMdv ,'^C3')                        
    ax1.plot(Rgx_UMgap,Dgx_UMgap,'oC0')
    ax1.plot(Rgx_UMoff,Dgx_UMoff,'sC0')
    ax1.plot(Rgx_UMdv ,Dgx_UMdv ,'^C0')
    
    ax1.plot(Rgx_RMgap,Dgx_RMgap,'C3',alpha=1)
    ax1.plot(Rgx_RMoff,Dgx_RMoff,'C3',alpha=0.75)
    ax1.plot(Rgx_RMdv ,Dgx_RMdv ,'C3',alpha=0.5)                        
    ax1.plot(Rgx_UMgap,Dgx_UMgap,'C0',alpha=1)
    ax1.plot(Rgx_UMoff,Dgx_UMoff,'C0',alpha=0.75)
    ax1.plot(Rgx_UMdv ,Dgx_UMdv ,'C0',alpha=0.5) 
    # ax1.errorbar(Rgx_RMgap,Dgx_RMgap,Dgx_RMgap_std,ecolor='C3',alpha=1)
    # ax1.errorbar(Rgx_RMoff,Dgx_RMoff,Dgx_RMoff_std,ecolor='C3',alpha=0.75)
    # ax1.errorbar(Rgx_RMdv ,Dgx_RMdv ,Dgx_RMdv_std,ecolor= 'C3',alpha=0.5)                        
    # ax1.errorbar(Rgx_UMgap,Dgx_UMgap,Dgx_UMgap_std,ecolor='C0',alpha=1)
    # ax1.errorbar(Rgx_UMoff,Dgx_UMoff,Dgx_UMoff_std,ecolor='C0',alpha=0.75)
    # ax1.errorbar(Rgx_UMdv ,Dgx_UMdv ,Dgx_UMdv_std,ecolor='C0',alpha=0.5) 
    
    ax1.plot(Rgx_RLgap,Dgx_RLgap,'oC3')
    ax1.plot(Rgx_RLoff,Dgx_RLoff,'sC3')
    ax1.plot(Rgx_RLdv ,Dgx_RLdv ,'^C3')
    ax1.plot(Rgx_ULgap,Dgx_ULgap,'oC0')
    ax1.plot(Rgx_ULoff,Dgx_ULoff,'sC0')
    ax1.plot(Rgx_ULdv ,Dgx_ULdv ,'^C0')
    
    ax1.plot(Rgx_RLgap,Dgx_RLgap,'C3--',alpha=1)
    ax1.plot(Rgx_RLoff,Dgx_RLoff,'C3--',alpha=0.75)
    ax1.plot(Rgx_RLdv ,Dgx_RLdv ,'C3--',alpha=0.5)  
    ax1.plot(Rgx_ULgap,Dgx_ULgap,'C0--',alpha=1)
    ax1.plot(Rgx_ULoff,Dgx_ULoff,'C0--',alpha=0.75)
    ax1.plot(Rgx_ULdv ,Dgx_ULdv ,'C0--',alpha=0.5) 
    
    
    ########### ESTRELLAS
    
    D_st30  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1])
    D_st50  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1])
    D_st100  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1])
    D_st200  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1])
    D_st500  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1])
    D_st1000 = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1])
    
    Dst = np.array([D_st30.D,D_st50.D,D_st100.D,D_st1000.D,D_st500.D,D_st200.D])
    
    D_st30M   = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mmas])
    D_st50M   = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mmas])
    D_st100M  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mmas])
    D_st200M  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mmas])
    D_st1000M = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mmas])
    D_st500M  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mmas])
    
    Dst_M = np.array([D_st30M.D,D_st50M.D,D_st100M.D,D_st1000M.D,D_st500M.D,D_st200M.D])
    
    D_st30L   = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][mlow])
    D_st50L   = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][mlow])
    D_st100L  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][mlow])
    D_st200L  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][mlow])
    D_st1000L = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][mlow])
    D_st500L  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][mlow])
    
    Dst_L = np.array([D_st30L.D,D_st50L.D,D_st100L.D,D_st1000L.D,D_st500L.D,D_st200L.D])
    
    D_st30R_gap   = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][C.mo2d_gap])
    D_st50R_gap   = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][C.mo2d_gap])
    D_st100R_gap  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][C.mo2d_gap])
    D_st200R_gap  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][C.mo2d_gap])
    D_st1000R_gap = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][C.mo2d_gap])
    D_st500R_gap  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][C.mo2d_gap])
    
    Dst_Rgap = np.array([D_st30R_gap.D,D_st50R_gap.D,D_st100R_gap.D,D_st1000R_gap.D,D_st500R_gap.D,D_st200R_gap.D])
    
    D_st30U_gap  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][C.mn2d_gap])
    D_st50U_gap = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][C.mn2d_gap])
    D_st100U_gap  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][C.mn2d_gap])
    D_st200U_gap  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][C.mn2d_gap])
    D_st1000U_gap = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][C.mn2d_gap])
    D_st500U_gap  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][C.mn2d_gap])
    
    Dst_Ugap = np.array([D_st30U_gap.D,D_st50U_gap.D,D_st100U_gap.D,D_st1000U_gap.D,D_st500U_gap.D,D_st200U_gap.D])
    
    D_st30R_off   = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][C.mo2d_off])
    D_st50R_off   = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][C.mo2d_off])
    D_st100R_off  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][C.mo2d_off])
    D_st200R_off  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][C.mo2d_off])
    D_st1000R_off = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][C.mo2d_off])
    D_st500R_off  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][C.mo2d_off])
    
    Dst_Roff = np.array([D_st30R_off.D,D_st50R_off.D,D_st100R_off.D,D_st1000R_off.D,D_st500R_off.D,D_st200R_off.D])
    
    D_st30U_off  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][C.mn2d_off])
    D_st50U_off  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][C.mn2d_off])
    D_st100U_off  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][C.mn2d_off])
    D_st200U_off  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][C.mn2d_off])
    D_st1000U_off = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][C.mn2d_off])
    D_st500U_off  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][C.mn2d_off])
    
    Dst_Uoff = np.array([D_st30U_off.D,D_st50U_off.D,D_st100U_off.D,D_st1000U_off.D,D_st500U_off.D,D_st200U_off.D])
    
    D_st30R_dv  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][C.mo2d_dv])
    D_st50R_dv  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][C.mo2d_dv])
    D_st100R_dv  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][C.mo2d_dv])
    D_st200R_dv  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][C.mo2d_dv])
    D_st1000R_dv = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][C.mo2d_dv])
    D_st500R_dv  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][C.mo2d_dv])
    
    Dst_Rdv = np.array([D_st30R_dv.D,D_st50R_dv.D,D_st100R_dv.D,D_st1000R_dv.D,D_st500R_dv.D,D_st200R_dv.D])
    
    D_st30U_dv  = Dcompute(cosangle(st30.a2D,DM200.a2D)[1][C.mn2d_dv])
    D_st50U_dv  = Dcompute(cosangle(st50.a2D,DM200.a2D)[1][C.mn2d_dv])
    D_st100U_dv  = Dcompute(cosangle(st100.a2D,DM200.a2D)[1][C.mn2d_dv])
    D_st200U_dv  = Dcompute(cosangle(st200.a2D,DM200.a2D)[1][C.mn2d_dv])
    D_st1000U_dv = Dcompute(cosangle(st1000.a2D,DM200.a2D)[1][C.mn2d_dv])
    D_st500U_dv  = Dcompute(cosangle(st500.a2D,DM200.a2D)[1][C.mn2d_dv])
    
    Dst_Udv = np.array([D_st30U_dv.D,D_st50U_dv.D,D_st100U_dv.D,D_st1000U_dv.D,D_st500U_dv.D,D_st200U_dv.D])
    
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
    
    
    ax2.plot(Rst_RMgap,Dst_RMgap,'oC3')
    ax2.plot(Rst_RMoff,Dst_RMoff,'sC3')
    ax2.plot(Rst_RMdv ,Dst_RMdv ,'^C3')
    ax2.plot(Rst_UMgap,Dst_UMgap,'oC0')
    ax2.plot(Rst_UMoff,Dst_UMoff,'sC0')
    ax2.plot(Rst_UMdv ,Dst_UMdv ,'^C0')
    
    ax2.plot(Rst_RMgap,Dst_RMgap,'C3',alpha=1,label='relaxed')
    ax2.plot(Rst_RMoff,Dst_RMoff,'C3',alpha=0.75)
    ax2.plot(Rst_RMdv ,Dst_RMdv ,'C3',alpha=0.5)
    ax2.plot(Rst_UMgap,Dst_UMgap,'C0',alpha=1,label='non-relaxed')
    ax2.plot(Rst_UMoff,Dst_UMoff,'C0',alpha=0.75)
    ax2.plot(Rst_UMdv ,Dst_UMdv ,'C0',alpha=0.5)
    
    ax2.plot(Rst_RLgap,Dst_RLgap,'oC3')
    ax2.plot(Rst_RLoff,Dst_RLoff,'sC3')
    ax2.plot(Rst_RLdv ,Dst_RLdv ,'^C3')
    ax2.plot(Rst_ULgap,Dst_ULgap,'oC0')
    ax2.plot(Rst_ULoff,Dst_ULoff,'sC0')
    ax2.plot(Rst_ULdv ,Dst_ULdv ,'^C0')
    
    ax2.plot(Rst_RLgap,Dst_RLgap,'C3--',alpha=1)
    ax2.plot(Rst_RLoff,Dst_RLoff,'C3--',alpha=0.75)
    ax2.plot(Rst_RLdv ,Dst_RLdv ,'C3--',alpha=0.5)
    ax2.plot(Rst_ULgap,Dst_ULgap,'C0--',alpha=1)
    ax2.plot(Rst_ULoff,Dst_ULoff,'C0--',alpha=0.75)
    ax2.plot(Rst_ULdv ,Dst_ULdv ,'C0--',alpha=0.5)

    ax2.plot(100,100,'k-',label='$\log(M_{200}/M_\odot) > 14.6$')
    ax2.plot(100,100,'k--',label='$\log(M_{200}/M_\odot) < 14.6$')



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


mmas = C.lM200p > 14.6
mlow = C.lM200p < 14.6

Rst = np.median(C.Rsp,axis=0)
Rgx = np.median(C.Rsp,axis=0)[3:]

Rgx_M = np.median(C.Rsp[mmas],axis=0)[3:]
Rgx_L = np.median(C.Rsp[mlow],axis=0)[3:]
Rgx_Rgap = np.median(C.Rsp[C.mo2d_gap],axis=0)[3:]
Rgx_Ugap = np.median(C.Rsp[C.mn2d_gap],axis=0)[3:]
Rgx_Roff = np.median(C.Rsp[C.mo2d_off],axis=0)[3:]
Rgx_Uoff = np.median(C.Rsp[C.mn2d_off],axis=0)[3:]
Rgx_Rdv = np.median(C.Rsp[C.mo2d_dv],axis=0)[3:]
Rgx_Udv = np.median(C.Rsp[C.mn2d_dv],axis=0)[3:]
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

Rst_M = np.median(C.Rsp[mmas],axis=0)
Rst_L = np.median(C.Rsp[mlow],axis=0)
Rst_Rgap = np.median(C.Rsp[C.mo2d_gap],axis=0)
Rst_Ugap = np.median(C.Rsp[C.mn2d_gap],axis=0)
Rst_Roff = np.median(C.Rsp[C.mo2d_off],axis=0)
Rst_Uoff = np.median(C.Rsp[C.mn2d_off],axis=0)
Rst_Rdv = np.median(C.Rsp[C.mo2d_dv],axis=0)
Rst_Udv = np.median(C.Rsp[C.mn2d_dv],axis=0)
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

####### PLOT
f, ax = plt.subplots(3,2, figsize=(12,12), sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0,0].text(0.45,0.9,'Galaxies')
ax[0,1].text(0.02,0.9,'Stars')

ax[0,0].text(0.8,0.37,r'$\theta^*_{DM}(R < R_{200})$')
ax[1,0].text(0.8,0.37,r'$\theta^*_{DM}(R < R_{500})$')
ax[2,0].text(0.8,0.37,r'$\theta^*_{DM}(R < R_{1000})$')


ax[0,1].set_xlim([0.0105,1.5])
ax[1,1].set_xlim([0.0105,1.5])
ax[2,1].set_xlim([0.0105,1.5])

ax[0,0].set_xscale('log')
ax[1,0].set_xscale('log')
ax[2,0].set_xscale('log')

ax[0,0].set_xlim([0.4,1.1])
ax[1,0].set_xlim([0.4,1.1])
ax[2,0].set_xlim([0.4,1.1])

ax[0,1].set_xscale('log')
ax[1,1].set_xscale('log')
ax[2,1].set_xscale('log')

ax[0,0].set_xticks([0.4,0.6,1.0])
ax[1,0].set_xticks([0.4,0.6,1.0])
ax[2,0].set_xticks([0.4,0.6,1.0])

ax[0,0].set_xticklabels([])
ax[1,0].set_xticklabels([])
ax[2,0].set_xticklabels(['0.4','0.6','1.0'])

ax[2,0].set_xlabel('$R/R_{200}$')
ax[2,1].set_xlabel('$R/R_{200}$')

ax[0,0].set_ylabel('$D$')
ax[1,0].set_ylabel('$D$')
ax[2,0].set_ylabel('$D$')

ax[0,0].set_ylim([0.25,1.0])

ax[1,1].plot(100,100,'ko',label=r'$M_{sat}/M_{BCG}$')
ax[1,1].plot(100,100,'ks',label=r'$D_{offset}$')
ax[1,1].plot(100,100,'k^',label=r'$\Delta V$')
ax[1,1].legend(frameon=False,loc=4)


D(DM200,ax[0,0],ax[0,1])
D(DM500,ax[1,0],ax[1,1])
D(DM1000,ax[2,0],ax[2,1])

ax[2,1].legend(frameon=False,loc=4)


plt.savefig(plotspath+'D_R.pdf',bbox_inches='tight')

###################

############

mc     = gx_c.N > 9
mcp    = np.array((mc.tolist())*3)
me     = gx_e.N > 9
mep    = np.array((me.tolist())*3)





Dgxc_RMgap  = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
Dgxc_RMoff  = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
Dgxc_RMdv   = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])
Dgxe_RMgap  = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mmas*C.mo2d_gap])
Dgxe_RMoff  = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mmas*C.mo2d_off])
Dgxe_RMdv   = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mmas*C.mo2d_dv])

Dgxc_UMgap  = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
Dgxc_UMoff  = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
Dgxc_UMdv   = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])
Dgxe_UMgap  = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mmas*C.mn2d_gap])
Dgxe_UMoff  = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mmas*C.mn2d_off])
Dgxe_UMdv   = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mmas*C.mn2d_dv])

Dgxc_RLgap  = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
Dgxc_RLoff  = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
Dgxc_RLdv   = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])
Dgxe_RLgap  = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mlow*C.mo2d_gap])
Dgxe_RLoff  = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mlow*C.mo2d_off])
Dgxe_RLdv   = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mlow*C.mo2d_dv])

Dgxc_ULgap  = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
Dgxc_ULoff  = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
Dgxc_ULdv   = Dcompute(cosangle(gx_c.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])
Dgxe_ULgap  = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mlow*C.mn2d_gap])
Dgxe_ULoff  = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mlow*C.mn2d_off])
Dgxe_ULdv   = Dcompute(cosangle(gx_e.a2D,DM200.a2D)[1][mlow*C.mn2d_dv])

'''
######

f, ax = plt.subplots(2,2, figsize=(10,10), sharey=True, sharex = True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0,0].hist(np.cos(2.*(np.deg2rad(Dgxe_RMdv.theta[Dgxe_RMdv.theta > 0]))),20,histtype='step',color='C7')
ax[0,0].hist(np.cos(2.*(np.deg2rad(Dgxc_RMdv.theta[Dgxc_RMdv.theta > 0]))),20,histtype='step',color='C1')
ax[0,0].axvline(Dgxe_RMdv.D,color='C7')
ax[0,0].axvline(Dgxc_RMdv.D,color='C1')
ax[0,0].axvline(Dgx_RMdv[-1],color='k',linestyle='--')
ax[0,0].text(-1,40,'$\log(M_{200}/M_\odot) > 14.6$')

ax[1,0].hist(np.cos(2.*(np.deg2rad(Dgxe_RLdv.theta[Dgxe_RLdv.theta > 0]))),20,histtype='step',color='C7')
ax[1,0].hist(np.cos(2.*(np.deg2rad(Dgxc_RLdv.theta[Dgxc_RLdv.theta > 0]))),20,histtype='step',color='C1')
ax[1,0].axvline(Dgx_RLdv[-1],color='k',linestyle='--')
ax[1,0].axvline(Dgxc_RLdv.D,color='C1',linestyle=':')
ax[1,0].axvline(Dgxe_RLdv.D,color='C7')
ax[1,0].text(-1,40,'$\log(M_{200}/M_\odot) < 14.6$')

ax[0,1].hist(np.cos(2.*(np.deg2rad(Dgxe_UMdv.theta[Dgxe_UMdv.theta > 0]))),20,histtype='step',color='C7')
ax[0,1].hist(np.cos(2.*(np.deg2rad(Dgxc_UMdv.theta[Dgxc_UMdv.theta > 0]))),20,histtype='step',color='C1')
ax[0,1].axvline(Dgx_UMdv[-1],color='k',linestyle='--',label='all galaxies')
ax[0,1].axvline(Dgxc_UMdv.D,color='C1',linestyle=':',label='concentrated')
ax[0,1].axvline(Dgxe_UMdv.D,color='C7',label='extended')
ax[0,1].legend(frameon=False)

ax[1,1].hist(np.cos(2.*(np.deg2rad(Dgxe_ULdv.theta[Dgxe_ULdv.theta > 0]))),20,histtype='step',color='C7')
ax[1,1].hist(np.cos(2.*(np.deg2rad(Dgxc_ULdv.theta[Dgxc_ULdv.theta > 0]))),20,histtype='step',color='C1')
ax[1,1].axvline(Dgxe_ULdv.D,color='C7')
ax[1,1].axvline(Dgxc_ULdv.D,color='C1',linestyle=':')
ax[1,1].axvline(Dgx_ULdv[-1],color='k',linestyle='--')

ax[0,0].set_ylabel('$N$')
ax[1,0].set_ylabel('$N$')
ax[1,0].set_xlabel(r'$\cos(2\theta)$')
ax[1,1].set_xlabel(r'$\cos(2\theta)$')

plt.savefig(plotspath+'cos2t_DM1000.pdf',bbox_inches='tight')
'''
