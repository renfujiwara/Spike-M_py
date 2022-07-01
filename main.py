import numpy as np
import argparse
import yaml
import sys
import matplotlib.pyplot as plt

import tool as tl
import mhpspike
import fit

isPlot=True

def plotsRNF(data, params, outfn):
    T=hparams['T']
    (idx,b,u)=mhpspike.spikeM(params, T)
    rmse=F_RNF(data, params, T, 'lin')
    
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

if __name__ == "__main__":
    #--- arguments ---#
    parser = argparse.ArgumentParser()
    parser.add_argument("-s",  "--seqfn", type=str, help="seqs filename")
    parser.add_argument("-o",  "--outdir",type=str, help="output dir")
    parser.add_argument("-p",  "--paramsfn",type=str, help="params file name")
    args = parser.parse_args()
    if(len(sys.argv)<2):
        parser.print_help(); tl.error("parser")
    #--- check sequencefn ---#
    if(args.seqfn!=None):
        seqfn = args.seqfn
    else: parser.print_help(); tl.error("parser")
    #--- check output dir ---#
    if(args.outdir!=None):
        outdir = args.outdir
    else: parser.print_help(); tl.error("parser")
    #--- check modelfn ---#
    if(args.paramsfn!=None):
        paramsfn= args.paramsfn
    else: parser.print_help(); tl.error("parser")
    
    with open(paramsfn) as file:
        hparams = yaml.safe_load(file)
        
    T = hparams['T']
    
    data = tl.loadseq(seqfn).T
    dat = data.flatten()
    dat = dat[:T]
    
    (rse,params,dic) = fit.LMFit(dat,hparams)
    with open("%s/best_params.yaml"%outdir, 'w') as file:
        yaml.dump(dic, file, default_flow_style=False)
    
    # plot result
    (idx,b,u)=mhpspike.spikeM(params, T)
    
    fig, ax = plt.subplots()
    ax.set_xlabel('t')  # x軸ラベル
    ax.set_ylabel('y')  # y軸ラベル
    # ax.set_title(r'$\sin(x)$ and $\cos(x)$') # グラフタイトル
    ax.plot(b,label = 'b')
    #ax.plot(u,label = 'u')
    ax.plot(dat,label = 'dat')
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(loc = 'upper right')
    ax.set_xlim([10,96])
    fig.savefig("%s/result.pdf"%outdir)
    
    