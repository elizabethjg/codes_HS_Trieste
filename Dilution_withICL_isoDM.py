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

def D(DM,ax1,ax2,ax3):
    
    ax1.axvspan(0, 0.07, alpha=0.3, color=cbcg)
    ax2.axvspan(0, 0.07, alpha=0.3, color=cbcg)
    ax3.axvspan(0, 0.07, alpha=0.3, color=cbcg)


    D_gx200  = Dcompute(cosangle(gx200.a2D,DM.a2Dicl)[1])
    D_gx500 = Dcompute(cosangle(gx500.a2D,DM.a2Dicl)[1])
    D_gx1000 = Dcompute(cosangle(gx1000.a2D,DM.a2Dicl)[1])
    
    Dgx     = np.array([D_gx1000.D,D_gx500.D,D_gx200.D])
    Dgx_std = np.array([D_gx1000.Dstd,D_gx500.D,D_gx200.Dstd])
        
    #################
    
    D_gxeRM_gap  = Dcompute(cosangle(gxe.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    D_gxcRM_gap  = Dcompute(cosangle(gxc.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    D_gx200RM_gap  = Dcompute(cosangle(gx200.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    D_gx1000RM_gap = Dcompute(cosangle(gx1000.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    D_gx500RM_gap  = Dcompute(cosangle(gx500.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    
    Dgx_RMgap = np.array([D_gx1000RM_gap.D,D_gx500RM_gap.D,D_gx200RM_gap.D])
    Dgx_RMgap_std = np.array([D_gx1000RM_gap.Dstd,D_gx500RM_gap.Dstd,D_gx200RM_gap.Dstd])
    
    D_gxeUM_gap  = Dcompute(cosangle(gxe.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    D_gxcUM_gap  = Dcompute(cosangle(gxc.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    D_gx200UM_gap  = Dcompute(cosangle(gx200.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    D_gx1000UM_gap = Dcompute(cosangle(gx1000.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    D_gx500UM_gap  = Dcompute(cosangle(gx500.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    
    Dgx_UMgap = np.array([D_gx1000UM_gap.D,D_gx500UM_gap.D,D_gx200UM_gap.D])
    Dgx_UMgap_std = np.array([D_gx1000UM_gap.Dstd,D_gx500UM_gap.Dstd,D_gx200UM_gap.Dstd])
    
    #######################
    
    D_gxeRL_gap  = Dcompute(cosangle(gxe.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    D_gxcRL_gap  = Dcompute(cosangle(gxc.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    D_gx200RL_gap  = Dcompute(cosangle(gx200.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    D_gx1000RL_gap = Dcompute(cosangle(gx1000.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    D_gx500RL_gap  = Dcompute(cosangle(gx500.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    
    Dgx_RLgap = np.array([D_gx1000RL_gap.D,D_gx500RL_gap.D,D_gx200RL_gap.D])
    Dgx_RLgap_std = np.array([D_gx1000RL_gap.Dstd,D_gx500RL_gap.Dstd,D_gx200RL_gap.Dstd])
    
    D_gxeUL_gap  = Dcompute(cosangle(gxe.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    D_gxcUL_gap  = Dcompute(cosangle(gxc.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    D_gx200UL_gap  = Dcompute(cosangle(gx200.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    D_gx1000UL_gap = Dcompute(cosangle(gx1000.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    D_gx500UL_gap  = Dcompute(cosangle(gx500.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    
    Dgx_ULgap = np.array([D_gx1000UL_gap.D,D_gx500UL_gap.D,D_gx200UL_gap.D])
    Dgx_ULgap_std = np.array([D_gx1000UL_gap.Dstd,D_gx500UL_gap.Dstd,D_gx200UL_gap.Dstd])
    
    
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
    ax1.plot(Rgx_UMgap,Dgx_UMgap,'C0^',markersize=8)    
    ax1.plot(Rgx_RMgap,Dgx_RMgap,'sienna',alpha=1)
    ax1.plot(Rgx_UMgap,Dgx_UMgap,'C0',alpha=1)
    
    ax1.plot(Rgx_RLgap,Dgx_RLgap,'v',color='sienna',markersize=8)
    ax1.plot(Rgx_ULgap,Dgx_ULgap,'C0v',markersize=8)
    ax1.plot(Rgx_RLgap,Dgx_RLgap,'--',color='sienna',markersize=8)
    ax1.plot(Rgx_ULgap,Dgx_ULgap,'C0--')   
    
    ########### ESTRELLAS
    
    D_st30  = Dcompute(cosangle(st30.a2D,DM.a2Dicl)[1])
    D_st50  = Dcompute(cosangle(st50.a2D,DM.a2Dicl)[1])
    D_st100  = Dcompute(cosangle(st100.a2D,DM.a2Dicl)[1])
    D_st200  = Dcompute(cosangle(st200.a2D,DM.a2Dicl)[1])
    D_st500  = Dcompute(cosangle(st500.a2D,DM.a2Dicl)[1])
    D_st1000 = Dcompute(cosangle(st1000.a2D,DM.a2Dicl)[1])
    
    Dst = np.array([D_st30.D,D_st50.D,D_st100.D,D_st1000.D,D_st500.D,D_st200.D])
    
    
    #################
    
    D_st30RM_gap   = Dcompute(cosangle(st30.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    D_st50RM_gap   = Dcompute(cosangle(st50.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    D_st100RM_gap  = Dcompute(cosangle(st100.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    D_st200RM_gap  = Dcompute(cosangle(st200.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    D_st1000RM_gap = Dcompute(cosangle(st1000.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    D_st500RM_gap  = Dcompute(cosangle(st500.a2D,DM.a2Dicl)[1][mmas*C.mo2d_gap])
    
    Dst_RMgap = np.array([D_st30RM_gap.D,D_st50RM_gap.D,D_st100RM_gap.D,D_st1000RM_gap.D,D_st500RM_gap.D,D_st200RM_gap.D])
    
    D_st30UM_gap  = Dcompute(cosangle(st30.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    D_st50UM_gap = Dcompute(cosangle(st50.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    D_st100UM_gap  = Dcompute(cosangle(st100.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    D_st200UM_gap  = Dcompute(cosangle(st200.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    D_st1000UM_gap = Dcompute(cosangle(st1000.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    D_st500UM_gap  = Dcompute(cosangle(st500.a2D,DM.a2Dicl)[1][mmas*C.mn2d_gap])
    
    Dst_UMgap = np.array([D_st30UM_gap.D,D_st50UM_gap.D,D_st100UM_gap.D,D_st1000UM_gap.D,D_st500UM_gap.D,D_st200UM_gap.D])
    
    
    #######################
    
    D_st30RL_gap   = Dcompute(cosangle(st30.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    D_st50RL_gap   = Dcompute(cosangle(st50.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    D_st100RL_gap  = Dcompute(cosangle(st100.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    D_st200RL_gap  = Dcompute(cosangle(st200.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    D_st1000RL_gap = Dcompute(cosangle(st1000.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    D_st500RL_gap  = Dcompute(cosangle(st500.a2D,DM.a2Dicl)[1][mlow*C.mo2d_gap])
    
    Dst_RLgap = np.array([D_st30RL_gap.D,D_st50RL_gap.D,D_st100RL_gap.D,D_st1000RL_gap.D,D_st500RL_gap.D,D_st200RL_gap.D])
    
    D_st30UL_gap  = Dcompute(cosangle(st30.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    D_st50UL_gap = Dcompute(cosangle(st50.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    D_st100UL_gap  = Dcompute(cosangle(st100.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    D_st200UL_gap  = Dcompute(cosangle(st200.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    D_st1000UL_gap = Dcompute(cosangle(st1000.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    D_st500UL_gap  = Dcompute(cosangle(st500.a2D,DM.a2Dicl)[1][mlow*C.mn2d_gap])
    
    Dst_ULgap = np.array([D_st30UL_gap.D,D_st50UL_gap.D,D_st100UL_gap.D,D_st1000UL_gap.D,D_st500UL_gap.D,D_st200UL_gap.D])
        
    ############
    
    ax2.plot(Rst,Dst,'k')
    
    
    ax2.plot(Rst_RMgap,Dst_RMgap,'^',color='sienna',markersize=8)
    ax2.plot(Rst_UMgap,Dst_UMgap,'C0^',markersize=8)
    
    ax2.plot(Rst_RMgap,Dst_RMgap,'sienna',alpha=1)
    ax2.plot(Rst_UMgap,Dst_UMgap,'C0',alpha=1)
    
    ax2.plot(Rst_RLgap,Dst_RLgap,'v',color='sienna',markersize=8)
    ax2.plot(Rst_ULgap,Dst_ULgap,'C0v',markersize=8)
    
    ax2.plot(Rst_RLgap,Dst_RLgap,'--',color='sienna',alpha=1)
    ax2.plot(Rst_ULgap,Dst_ULgap,'C0--',alpha=1)

    ax3.plot(100,100,'k^',label='$\log(M_{200}/M_\odot) > 14.8$',markersize=8)
    ax3.plot(100,100,'kv',label='$\log(M_{200}/M_\odot) < 14.8$',markersize=8)
    
    ##############################
    ## ICL
    ##############################
    
    tbcg = cosangle(ICL_.a2Dbcg,DM.a2Dicl)[1]
    ticl = cosangle(ICL_.a2Dicl,DM.a2Dicl)[1]
    
    D_bcg   = Dcompute(tbcg).D
    D_icl   = Dcompute(ticl).D    
    
    D_bcgRM_gap   = Dcompute(tbcg[C.mo2d_gap*mmas*mbcg]).D
    D_iclRM_gap   = Dcompute(ticl[C.mo2d_gap*mmas]).D
    D_bcgUM_gap   = Dcompute(tbcg[C.mn2d_gap*mmas*mbcg]).D
    D_iclUM_gap   = Dcompute(ticl[C.mn2d_gap*mmas]).D

    D_bcgRL_gap   = Dcompute(tbcg[C.mo2d_gap*mlow*mbcg]).D
    D_iclRL_gap   = Dcompute(ticl[C.mo2d_gap*mlow]).D
    D_bcgUL_gap   = Dcompute(tbcg[C.mn2d_gap*mlow*mbcg]).D
    D_iclUL_gap   = Dcompute(ticl[C.mn2d_gap*mlow]).D

    
    ax3.plot([Rbcg,Ricl],[D_bcg,D_icl],'k')    
        

    ax3.plot([Rbcg_RMgap,Ricl_RMgap],[D_bcgRM_gap,D_iclRM_gap],'^',color='sienna',markersize=8)
    ax3.plot([Rbcg_UMgap,Ricl_UMgap],[D_bcgUM_gap,D_iclUM_gap],'C0^',markersize=8)
    ax3.plot([Rbcg_RMgap,Ricl_RMgap],[D_bcgRM_gap,D_iclRM_gap],'',color='sienna',markersize=8)
    ax3.plot([Rbcg_UMgap,Ricl_UMgap],[D_bcgUM_gap,D_iclUM_gap],'C0')

    ax3.plot([Rbcg_RLgap,Ricl_RLgap],[D_bcgRL_gap,D_iclRL_gap],'v',color='sienna',markersize=8)
    ax3.plot([Rbcg_ULgap,Ricl_ULgap],[D_bcgUL_gap,D_iclUL_gap],'C0v',markersize=8)
    ax3.plot([Rbcg_RLgap,Ricl_RLgap],[D_bcgRL_gap,D_iclRL_gap],'--',color='sienna',markersize=8)
    ax3.plot([Rbcg_ULgap,Ricl_ULgap],[D_bcgUL_gap,D_iclUL_gap],'C0--')
    
def D_tensor(DM,ax1,ax2,ax3):

    
    ax1.axvspan(0, 0.07, alpha=0.3, color=cbcg)
    ax2.axvspan(0, 0.07, alpha=0.3, color=cbcg)
    ax3.axvspan(0, 0.07, alpha=0.3, color=cbcg)


    D_gx200  = Dcompute(cosangle(gx200.a2D,DM.a2D)[1])
    D_gx500 = Dcompute(cosangle(gx500.a2D,DM.a2D)[1])
    D_gx1000 = Dcompute(cosangle(gx1000.a2D,DM.a2D)[1])
    
    Dgx     = np.array([D_gx1000.D,D_gx500.D,D_gx200.D])
    Dgx_std = np.array([D_gx1000.Dstd,D_gx500.D,D_gx200.Dstd])
        
    #################
    
    D_gxeRM_gap  = Dcompute(cosangle(gxe.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    D_gxcRM_gap  = Dcompute(cosangle(gxc.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    D_gx200RM_gap  = Dcompute(cosangle(gx200.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    D_gx1000RM_gap = Dcompute(cosangle(gx1000.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    D_gx500RM_gap  = Dcompute(cosangle(gx500.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    
    Dgx_RMgap = np.array([D_gx1000RM_gap.D,D_gx500RM_gap.D,D_gx200RM_gap.D])
    Dgx_RMgap_std = np.array([D_gx1000RM_gap.Dstd,D_gx500RM_gap.Dstd,D_gx200RM_gap.Dstd])
    
    D_gxeUM_gap  = Dcompute(cosangle(gxe.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    D_gxcUM_gap  = Dcompute(cosangle(gxc.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    D_gx200UM_gap  = Dcompute(cosangle(gx200.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    D_gx1000UM_gap = Dcompute(cosangle(gx1000.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    D_gx500UM_gap  = Dcompute(cosangle(gx500.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    
    Dgx_UMgap = np.array([D_gx1000UM_gap.D,D_gx500UM_gap.D,D_gx200UM_gap.D])
    Dgx_UMgap_std = np.array([D_gx1000UM_gap.Dstd,D_gx500UM_gap.Dstd,D_gx200UM_gap.Dstd])
    
    #######################
    
    D_gxeRL_gap  = Dcompute(cosangle(gxe.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    D_gxcRL_gap  = Dcompute(cosangle(gxc.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    D_gx200RL_gap  = Dcompute(cosangle(gx200.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    D_gx1000RL_gap = Dcompute(cosangle(gx1000.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    D_gx500RL_gap  = Dcompute(cosangle(gx500.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    
    Dgx_RLgap = np.array([D_gx1000RL_gap.D,D_gx500RL_gap.D,D_gx200RL_gap.D])
    Dgx_RLgap_std = np.array([D_gx1000RL_gap.Dstd,D_gx500RL_gap.Dstd,D_gx200RL_gap.Dstd])
    
    D_gxeUL_gap  = Dcompute(cosangle(gxe.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    D_gxcUL_gap  = Dcompute(cosangle(gxc.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    D_gx200UL_gap  = Dcompute(cosangle(gx200.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    D_gx1000UL_gap = Dcompute(cosangle(gx1000.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    D_gx500UL_gap  = Dcompute(cosangle(gx500.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    
    Dgx_ULgap = np.array([D_gx1000UL_gap.D,D_gx500UL_gap.D,D_gx200UL_gap.D])
    Dgx_ULgap_std = np.array([D_gx1000UL_gap.Dstd,D_gx500UL_gap.Dstd,D_gx200UL_gap.Dstd])
    
    
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
    ax1.plot(Rgx_UMgap,Dgx_UMgap,'C0^',markersize=8)    
    ax1.plot(Rgx_RMgap,Dgx_RMgap,'sienna',alpha=1)
    ax1.plot(Rgx_UMgap,Dgx_UMgap,'C0',alpha=1)
    
    ax1.plot(Rgx_RLgap,Dgx_RLgap,'v',color='sienna',markersize=8)
    ax1.plot(Rgx_ULgap,Dgx_ULgap,'C0v',markersize=8)
    ax1.plot(Rgx_RLgap,Dgx_RLgap,'--',color='sienna',markersize=8)
    ax1.plot(Rgx_ULgap,Dgx_ULgap,'C0--')   
    
    ########### ESTRELLAS
    
    D_st30  = Dcompute(cosangle(st30.a2D,DM.a2D)[1])
    D_st50  = Dcompute(cosangle(st50.a2D,DM.a2D)[1])
    D_st100  = Dcompute(cosangle(st100.a2D,DM.a2D)[1])
    D_st200  = Dcompute(cosangle(st200.a2D,DM.a2D)[1])
    D_st500  = Dcompute(cosangle(st500.a2D,DM.a2D)[1])
    D_st1000 = Dcompute(cosangle(st1000.a2D,DM.a2D)[1])
    
    Dst = np.array([D_st30.D,D_st50.D,D_st100.D,D_st1000.D,D_st500.D,D_st200.D])
    
    
    #################
    
    D_st30RM_gap   = Dcompute(cosangle(st30.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    D_st50RM_gap   = Dcompute(cosangle(st50.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    D_st100RM_gap  = Dcompute(cosangle(st100.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    D_st200RM_gap  = Dcompute(cosangle(st200.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    D_st1000RM_gap = Dcompute(cosangle(st1000.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    D_st500RM_gap  = Dcompute(cosangle(st500.a2D,DM.a2D)[1][mmas*C.mo2d_gap])
    
    Dst_RMgap = np.array([D_st30RM_gap.D,D_st50RM_gap.D,D_st100RM_gap.D,D_st1000RM_gap.D,D_st500RM_gap.D,D_st200RM_gap.D])
    
    D_st30UM_gap  = Dcompute(cosangle(st30.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    D_st50UM_gap = Dcompute(cosangle(st50.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    D_st100UM_gap  = Dcompute(cosangle(st100.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    D_st200UM_gap  = Dcompute(cosangle(st200.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    D_st1000UM_gap = Dcompute(cosangle(st1000.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    D_st500UM_gap  = Dcompute(cosangle(st500.a2D,DM.a2D)[1][mmas*C.mn2d_gap])
    
    Dst_UMgap = np.array([D_st30UM_gap.D,D_st50UM_gap.D,D_st100UM_gap.D,D_st1000UM_gap.D,D_st500UM_gap.D,D_st200UM_gap.D])
    
    
    #######################
    
    D_st30RL_gap   = Dcompute(cosangle(st30.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    D_st50RL_gap   = Dcompute(cosangle(st50.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    D_st100RL_gap  = Dcompute(cosangle(st100.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    D_st200RL_gap  = Dcompute(cosangle(st200.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    D_st1000RL_gap = Dcompute(cosangle(st1000.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    D_st500RL_gap  = Dcompute(cosangle(st500.a2D,DM.a2D)[1][mlow*C.mo2d_gap])
    
    Dst_RLgap = np.array([D_st30RL_gap.D,D_st50RL_gap.D,D_st100RL_gap.D,D_st1000RL_gap.D,D_st500RL_gap.D,D_st200RL_gap.D])
    
    D_st30UL_gap  = Dcompute(cosangle(st30.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    D_st50UL_gap = Dcompute(cosangle(st50.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    D_st100UL_gap  = Dcompute(cosangle(st100.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    D_st200UL_gap  = Dcompute(cosangle(st200.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    D_st1000UL_gap = Dcompute(cosangle(st1000.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    D_st500UL_gap  = Dcompute(cosangle(st500.a2D,DM.a2D)[1][mlow*C.mn2d_gap])
    
    Dst_ULgap = np.array([D_st30UL_gap.D,D_st50UL_gap.D,D_st100UL_gap.D,D_st1000UL_gap.D,D_st500UL_gap.D,D_st200UL_gap.D])
        
    ############
    
    ax2.plot(Rst,Dst,'k')
    
    
    ax2.plot(Rst_RMgap,Dst_RMgap,'^',color='sienna',markersize=8)
    ax2.plot(Rst_UMgap,Dst_UMgap,'C0^',markersize=8)
    
    ax2.plot(Rst_RMgap,Dst_RMgap,'sienna',alpha=1)
    ax2.plot(Rst_UMgap,Dst_UMgap,'C0',alpha=1)
    
    ax2.plot(Rst_RLgap,Dst_RLgap,'v',color='sienna',markersize=8)
    ax2.plot(Rst_ULgap,Dst_ULgap,'C0v',markersize=8)
    
    ax2.plot(Rst_RLgap,Dst_RLgap,'--',color='sienna',alpha=1)
    ax2.plot(Rst_ULgap,Dst_ULgap,'C0--',alpha=1)

    ax3.plot(100,100,'k^',label='$\log(M_{200}/M_\odot) > 14.8$',markersize=8)
    ax3.plot(100,100,'kv',label='$\log(M_{200}/M_\odot) < 14.8$',markersize=8)
    
    ##############################
    ## ICL
    ##############################
    
    tbcg = cosangle(ICL_.a2Dbcg,DM.a2D)[1]
    ticl = cosangle(ICL_.a2Dicl,DM.a2D)[1]
    
    D_bcg   = Dcompute(tbcg).D
    D_icl   = Dcompute(ticl).D    
    
    D_bcgRM_gap   = Dcompute(tbcg[C.mo2d_gap*mmas*mbcg]).D
    D_iclRM_gap   = Dcompute(ticl[C.mo2d_gap*mmas]).D
    D_bcgUM_gap   = Dcompute(tbcg[C.mn2d_gap*mmas*mbcg]).D
    D_iclUM_gap   = Dcompute(ticl[C.mn2d_gap*mmas]).D

    D_bcgRL_gap   = Dcompute(tbcg[C.mo2d_gap*mlow*mbcg]).D
    D_iclRL_gap   = Dcompute(ticl[C.mo2d_gap*mlow]).D
    D_bcgUL_gap   = Dcompute(tbcg[C.mn2d_gap*mlow*mbcg]).D
    D_iclUL_gap   = Dcompute(ticl[C.mn2d_gap*mlow]).D

    
    ax3.plot([Rbcg,Ricl],[D_bcg,D_icl],'k')    
        

    ax3.plot([Rbcg_RMgap,Ricl_RMgap],[D_bcgRM_gap,D_iclRM_gap],'^',color='sienna',markersize=8)
    ax3.plot([Rbcg_UMgap,Ricl_UMgap],[D_bcgUM_gap,D_iclUM_gap],'C0^',markersize=8)
    ax3.plot([Rbcg_RMgap,Ricl_RMgap],[D_bcgRM_gap,D_iclRM_gap],'',color='sienna',markersize=8)
    ax3.plot([Rbcg_UMgap,Ricl_UMgap],[D_bcgUM_gap,D_iclUM_gap],'C0')

    ax3.plot([Rbcg_RLgap,Ricl_RLgap],[D_bcgRL_gap,D_iclRL_gap],'v',color='sienna',markersize=8)
    ax3.plot([Rbcg_ULgap,Ricl_ULgap],[D_bcgUL_gap,D_iclUL_gap],'C0v',markersize=8)
    ax3.plot([Rbcg_RLgap,Ricl_RLgap],[D_bcgRL_gap,D_iclRL_gap],'--',color='sienna',markersize=8)
    ax3.plot([Rbcg_ULgap,Ricl_ULgap],[D_bcgUL_gap,D_iclUL_gap],'C0--')

plotspath = '../final_plots/'

C = Clusters()
ICL_ = ICL()
DM   = ICL(trazer='dm',cutmax=False)

m = (C.sub == 0)
micl2 = (m.tolist()*3)

C.mo2d_gap = C.mo2d_gap[micl2]
C.mn2d_gap = C.mn2d_gap[micl2]
C.Rsp      = C.Rsp[micl2]


ltimep = C.ltimep[micl2]
C.mo2d_gap = ltimep > 4.5
C.mn2d_gap = ltimep <= 4.5


gx200 = Galaxias(radio=200)
gx1000 = Galaxias(radio=1000)
gx500 = Galaxias(radio=500)

st200  = Stars(200)
st1000 = Stars(1000)
st500  = Stars(500)
st30   = Stars(30)
st50   = Stars(50)
st100  = Stars(100)

DM1000 = DarkMatter(1000)
DM500  = DarkMatter(500)
DM200  = DarkMatter(200)


gxe = Galaxias('extendidas')
gxc = Galaxias('concentradas')

st30.a2D   = st30.a2D[micl2]  
st50.a2D   = st50.a2D[micl2]  
st100.a2D  = st100.a2D[micl2] 
st200.a2D  = st200.a2D[micl2] 
st1000.a2D = st1000.a2D[micl2]
st500.a2D  = st500.a2D[micl2] 

gx200.a2D = gx200.a2D[micl2]  
gx500.a2D = gx500.a2D[micl2]  
gx1000.a2D = gx1000.a2D[micl2]  

DM200.a2D = DM200.a2D[micl2]  
DM500.a2D = DM500.a2D[micl2]  
DM1000.a2D = DM1000.a2D[micl2]  

gxe.a2D = gxe.a2D[micl2]
gxc.a2D = gxc.a2D[micl2]

mmas = (C.lM200p > 14)[micl2]
mlow = (C.lM200p < 16)[micl2]

mbcg = (ICL_.Rbcg > 0)

yaxis = np.zeros((len(DM.a2Dicl),2))
yaxis[:,1] = 1.

mall = np.ones(72).astype(bool)
mall2D = np.ones(72*3).astype(bool)



Rst = np.median(C.Rsp,axis=0)
Rgx = np.median(C.Rsp,axis=0)[3:]

Rgx_RMgap = np.median(C.Rsp[mmas*C.mo2d_gap],axis=0)[3:]
Rgx_UMgap = np.median(C.Rsp[mmas*C.mn2d_gap],axis=0)[3:]
Rgx_RLgap = np.median(C.Rsp[mlow*C.mo2d_gap],axis=0)[3:]
Rgx_ULgap = np.median(C.Rsp[mlow*C.mn2d_gap],axis=0)[3:]

Rst_RMgap = np.median(C.Rsp[mmas*C.mo2d_gap],axis=0)
Rst_UMgap = np.median(C.Rsp[mmas*C.mn2d_gap],axis=0)
Rst_RLgap = np.median(C.Rsp[mlow*C.mo2d_gap],axis=0)
Rst_ULgap = np.median(C.Rsp[mlow*C.mn2d_gap],axis=0)

Rbcg      = np.median(ICL_.Rbcg[mbcg])

Rbcg_RMgap = np.median(ICL_.Rbcg[mmas*(C.mo2d_gap)*mbcg])
Rbcg_UMgap = np.median(ICL_.Rbcg[mmas*(C.mn2d_gap)*mbcg])
Rbcg_RLgap = np.median(ICL_.Rbcg[mlow*(C.mo2d_gap)*mbcg])
Rbcg_ULgap = np.median(ICL_.Rbcg[mlow*(C.mn2d_gap)*mbcg])


Ricl      = np.median(ICL_.Ricl)

Ricl_RMgap = np.median(ICL_.Ricl[mmas*(C.mo2d_gap)])
Ricl_UMgap = np.median(ICL_.Ricl[mmas*(C.mn2d_gap)])
Ricl_RLgap = np.median(ICL_.Ricl[mlow*(C.mo2d_gap)])
Ricl_ULgap = np.median(ICL_.Ricl[mlow*(C.mn2d_gap)])


####### PLOT
f, ax = plt.subplots(1,3, figsize=(12,4), sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].text(1.08,1.,'c')
ax[0].text(1.21,1.,'e')
ax[0].text(0.6,1.,'Galaxies')
ax[1].text(0.1,1.,'Stars')
ax[2].text(0.13,1.,'Isophotes')

ax[1].text(0.02,1.,'BCG',c=cbcg)
ax[2].text(0.052,1.,'BCG',c=cbcg)




ax[0].set_xscale('log')
ax[0].set_xlim([0.4,1.3])

ax[1].set_xscale('log')
ax[1].set_xlim([0.0105,1.5])

ax[2].set_xscale('log')
ax[2].set_xlim([0.05,0.8])




ax[2].xaxis.set_major_formatter(NullFormatter())
ax[2].xaxis.set_minor_formatter(NullFormatter())
ax[0].xaxis.set_major_formatter(NullFormatter())
ax[0].xaxis.set_minor_formatter(NullFormatter())

ax[0].set_xticks([0.5,0.7,1.0])
ax[0].set_xticklabels(['0.5','0.7','1.0'])

ax[1].set_xticks([0.1,0.3,0.5,1.0])
ax[1].set_xticklabels(['0.1','0.3','0.5','1.0'])

ax[2].set_xticks([0.1,0.2,0.3,0.5])
ax[2].set_xticklabels(['0.1','0.2','0.3','0.5'])


ax[0].set_xlabel('$R/R_{200}$')
ax[1].set_xlabel('$R/R_{200}$')
ax[2].set_xlabel('$R/R_{200}$')

ax[0].set_ylabel('$D$')

ax[0].set_ylim([0.5,1.1])


D(DM,ax[0],ax[1],ax[2])


ax[2].plot(100,100,'sienna',label=r'relaxed')
ax[2].plot(100,100,'C0',label=r'non-relaxed')
ax[2].legend(frameon=False,loc=4)


plt.savefig(plotspath+'D_R_withicl_isoDM.pdf',bbox_inches='tight')


#########################

f, ax = plt.subplots(1,3, figsize=(12,4), sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].text(1.08,1.,'c')
ax[0].text(1.21,1.,'e')
ax[0].text(0.6,1.,'Galaxies')
ax[1].text(0.1,1.,'Stars')
ax[2].text(0.13,1.,'Isophotes')

ax[1].text(0.02,1.,'BCG',c=cbcg)
ax[2].text(0.052,1.,'BCG',c=cbcg)


ax[0].set_xscale('log')
ax[0].set_xlim([0.4,1.3])

ax[1].set_xscale('log')
ax[1].set_xlim([0.0105,1.5])

ax[2].set_xscale('log')
ax[2].set_xlim([0.05,0.8])


ax[2].xaxis.set_major_formatter(NullFormatter())
ax[2].xaxis.set_minor_formatter(NullFormatter())
ax[0].xaxis.set_major_formatter(NullFormatter())
ax[0].xaxis.set_minor_formatter(NullFormatter())

ax[0].set_xticks([0.5,0.7,1.0])
ax[0].set_xticklabels(['0.5','0.7','1.0'])

ax[1].set_xticks([0.1,0.3,0.5,1.0])
ax[1].set_xticklabels(['0.1','0.3','0.5','1.0'])

ax[2].set_xticks([0.1,0.2,0.3,0.5])
ax[2].set_xticklabels(['0.1','0.2','0.3','0.5'])


ax[0].set_xlabel('$R/R_{200}$')
ax[1].set_xlabel('$R/R_{200}$')
ax[2].set_xlabel('$R/R_{200}$')

ax[0].set_ylabel('$D$')

ax[0].set_ylim([0.5,1.1])


D_tensor(DM1000,ax[0],ax[1],ax[2])


ax[2].plot(100,100,'sienna',label=r'relaxed')
ax[2].plot(100,100,'C0',label=r'non-relaxed')
ax[2].legend(frameon=False,loc=4)


plt.savefig(plotspath+'D_R_withicl_ineDM.pdf',bbox_inches='tight')
