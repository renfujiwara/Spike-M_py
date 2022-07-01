import numpy as np
import argparse
import yaml
import sys

import tool as tl
import mhpspike
import fit

def exogenous(t, nc, Sc):
    if(t==nc):
        return Sc
    else:
        return 0
  
def decay_pl(n, beta, slope):
    if(n<1):
        return(0)
    else:
        return( beta* np.power(n, slope))
    
def period(t, Pp, Ps, Pa):
    if(Pa == 0):
        value = 1
    else:
        value = 1 - 0.5 * Pa * (np.sin(2*np.pi*(t + Ps)/Pp)+1)
    return value
    
# --- input ------------------------------
# T:   duration of the sequence
# params: SpikeM parameters 
#   params(1): N
#   params(2): beta*N
#   params(3): slope
#   params(4): nb  (starting time of breaking news)
#   params(5): Sb  (strength of external shock)
#   params(6): bgn (background noise)
#   params(7): Pp  (period)
#   params(8): Pa  (strength of periodicity
#   params(9): Ps  (phase Ps of periodicity
#   params(10):B0  (default, B0=0;)
# ----------------------------------------

#--- output ------------------------------
# idx:	index (size: Tx1)
# dB:    count of informed bloggers (size: Tx1)
# U:     count of un-informed bloggers (size: Tx1)
#----------------------------------------

def spikeM(params0,T):

    #--- params (please ignore)
    params = np.abs(params0)
    
    N = params[0]
    betaN = params[1]
    beta0 =  betaN/N
    slope = -abs(params[2])
    
    nc = round(params[3])
    Sc = params[4]
    bgnoise = params[5]
    
    Pp = params[6]
    Pa = params[7]
    Ps = params[8]
    
    wd = 1
    T=T*wd
    nc=nc*wd
    bgnoise=bgnoise/wd
    Sc=Sc/wd
    #--- for multi-wdsize (please ignore)

    U  = np.full(T,N)  # uninfected
    dB = np.zeros(T)  # infected (total)
    B  = np.zeros(T)  # blogged/infected

    # init
    B0=0
    U[0]  = N-B0
    dB[0] = B0
    B[0]  = B0

    for n in range(0,T-1):
        dsum = 0    # ~ number of sneezes
        for i in range(nc,n+1): # i.e, from zero to n
            Si = exogenous(i, nc, Sc)
            dsum += ( dB[i] + Si ) * decay_pl(n+1-i, beta0, slope)
        P_n1 = period(n+1, Pp, Ps, Pa)
        dB[n+1] = P_n1 * (U[n] * dsum + bgnoise)    

        #upper-bound b[n+1]:
        if( dB[n+1] > U[n]):
            dB[n+1] = U[n]
  
        U[n+1] = U[n] - dB[n+1] 
        B[n+1] = B[n] + dB[n+1]
        
        assert( abs(B[n+1] + U[n+1] - N) < 0.001)
        # if(abs(B[n+1] + U[n+1] - N) > 0.001 or U[n+1]==0):
        #     dB=np.full(T,np.nan)
        #     B=np.full(T,np.nan) 
        #     U=np.full(T,np.nan)

    #--- for multi-wdsize (please ignore)
    if(wd != 1):
        b_tmp = np.zeros(T)
        for j in range(0,np.floor(T/wd).astype(int)):
            for jj in range(0,wd):
                b_tmp[j] += dB[j*wd+jj]
        dB = b_tmp.copy()
        T=np.floor(T/wd).astype(int)
        
        #--- for multi-wdsize (please ignore)
        dB=dB[0:T]
        B=B[0:T]
        U=U[0:T]
    idx=T
    
    return (idx,dB, U)