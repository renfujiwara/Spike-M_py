import sys
import pylab as pl

def error(msg):
    print("====================")
    print(" error              ")
    print("--------------------")
    print("%s"%msg)
    print("====================")
    #raise
    sys.exit(0);
    
def loadseq(fn):
    tmp=pl.loadtxt(fn, dtype = float, ndmin=2).T
    return tmp