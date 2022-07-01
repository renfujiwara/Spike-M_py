from warnings import catch_warnings
import numpy as np
import time
import sys
import lmfit

import mhpspike
import tool as tl

#-----------------------------#
# lmfit (default)
XTL=1.e-8
FTL=1.e-8
MAXFEV=100

INF = 1.e+10
TH = 1.e-4

#----------------------------#
#    model fitting           #
#----------------------------#
def LMFit(data, hparams):
    T = hparams['T']
    pfreq = hparams['pfreq']
    ITER = hparams['ITER']
    params = init_params(data, pfreq)
    RSEP=INF
    (LBs,UBs) = init_const(data, pfreq, T)
    order=np.array([3,4,1,0])
    if(pfreq != -1):
        order = np.append(order,8)
        order = np.append(order,7)
    order = np.append(order,5)
    for i in range(ITER):
        params0 = params.copy()
        if(i<(ITER/2)):
            type='log' 
        else:
            type='lin'
        for lo in range(len(order)):
            loc = order[lo]
            if(loc==5 and (type=='lin')):
                type='log'
            if(loc==3):
                params[loc] = FD_search(data, params, loc, T, type)
            else:
                try:
                    params = _nl_fit(data, params, loc, T, type)
                except:
                    params[loc] = params0[loc]
                params = const(params, LBs, UBs, loc)
        RSE = printRNF(data, params, T)
        if( (abs(RSEP-RSE)<TH)  or (RSE < 0.8*(10**-3)) ):
            break
        RSEP = RSE
    dic = dict_params(params)
    return (RSE,params,dic)

def init_params(data, Pp):
    # init params
    # RNF-base
    N_max = np.sum(data); 
    betaN = 1.0
    slope = -1.5
    # RNF-X
    nc=0
    Sc=0.1
    bgnoise=0.01 
    # RNF-P
    Pa=0.1
    Ps=0
    
    params=np.zeros(9)
    params[0] = N_max
    params[1] = betaN
    params[2] = slope
    params[3] = nc
    params[4] = Sc
    params[5] = bgnoise
    params[6] = Pp
    params[7] = Pa
    params[8] = Ps
    return params

def dict_params(params):
    dic = {
        'N_max':str(params[0]),
        'betaN':str(params[1]),
        'slope':str(params[2]),
        'nc':str(params[3]),
        'Sc':str(params[4]),
        'bgnoise':str(params[5]),
        'Pp':str(params[6]),
        'Pa':str(params[7]),
        'Ps':str(params[8])
    }
    return dic

def init_const(data, pfreq, T):
    LBs = np.array([np.sum(data), 0.01, -1.5, 0, 0.0, 0, pfreq, 0.05, 0])
    UBs = np.array([INF, 2.0, -1.5, T/2, INF, INF, pfreq, 1, pfreq])
    return (LBs,UBs)

def const(params, LBs, UBs, loc):
    params[8] = params[8] % params[6]
    params[loc] = max(LBs[loc],params[loc])
    params[loc] = min(UBs[loc],params[loc])
    return params

def FD_search(data, params0, loc, T, type):
    th=1.0
    params = params0.copy()
    # if starting point is too sparce, then, ignore the point
    spWD=4
    dat=removeSparse(data, spWD, th)
    loclist=np.where(dat>th)
    
    if(len(loclist) == 0):
        st=0
    else:
        st=loclist[0][0]-1
        
    if(st<0):
        st=0
    th = np.amax(dat)
    loclist = np.where(dat==th)
    ed=loclist[0][0]
    idxlen = ed - st + 1
    sselist = np.zeros(idxlen)
    
    for i in range(idxlen):
        params[loc] = i + st
        sselist[i]=F_RNF(dat, params, T, type)
        
    minlst=np.where(sselist==np.min(sselist))
    
    if(len(minlst) == 0):
        estimate = 0
    else:
        estimate = minlst[0][0] + st 
        
    return int(estimate)


def removeSparse(X, wd, th):
    wd=np.ceil(wd/2).astype(int)
    n=len(X)
    for t in range(n):
        st=t-wd
        ed=t+wd
        if(st<0):
            st=0
        if(ed>n):
            ed=n
        tmp = X[st:ed]
        counts=np.sum(tmp[tmp<th])
        length=ed-st
        if(counts>length/2):
            X[t]=0
    return X

def _nl_fit(data, params, loc, T, type):
    #---------------------------------------#
    # (1) create param set
    P=_createP(params,loc)
    #---------------------------------------#
    # (2) start lmfit
    lmsol = lmfit.Minimizer(_distfunc_RNF, P, fcn_args=(data, params, loc, T, type))
    res=lmsol.leastsq(xtol=XTL, ftol=FTL, max_nfev=MAXFEV)
    #---------------------------------------#
    # (3) update param set
    params=_updateP(res.params, params, loc)
    #---------------------------------------#
    return params


#----------------------------#
def _createP(params, loc):
    P = lmfit.Parameters()
    #PARAM_MX=0.1
    #pm=PARAM_MX
    #PARAM_INI=1.e-4 #6
    V=True
    P.add('param', value=params[loc], vary=V)
    return P
#----------------------------#
def _updateP(P, params, loc):
    params[loc]=P['param'].value
    return params
#----------------------------#

#--------------------------------------#
#    objective functions  #
#    return the array to be minimized
#--------------------------------------#
def _distfunc_RNF(P, data, params, loc, T, type):
    params=_updateP(P, params, loc)
    diff = F_RNF(data, params, T, type)
    return diff

def F_RNF(data, params, T, type):
    (idx,b,u)=mhpspike.spikeM(params, T)
    if(type =='log'):
        b=np.log(b+1)
        data=np.log(data+1)
    elif(type == 'R5'):
        b=np.power(b,1/5)
        data=np.power(data,1/5)
        
    diff = np.sqrt(np.mean(np.power(b-data,2)))
    return diff

def printRNF(dat, params, T):
    RSE_LIN=F_RNF(dat, params, T, 'lin')
    print(f'N = {str(params[0])}')
    print(f'beta*N = {str(params[1])}')
    print(f'slope = {str(params[2])}')
    print(f'nc = {str(params[3])}')
    print(f'Sc = {str(params[4])}')
    print(f'bgnoise = {str(params[5])}')
    print(f'pcycle (Pp, Pa, Ps) = ({str(params[6])},{str(params[7])},{str(params[8])})')
    print(f'err = {str(RSE_LIN)}')    
    return RSE_LIN