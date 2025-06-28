import sys
import time
import random
import numpy as np
import re
def make_MboI(genome = 'hg19',add_chr=True):
    if genome=='hg19_dm3':
        file="/home/chenqw/lizj/fly_analysis/hg19_dm3/hg19_dm3_MboI.txt"
        chrs = [str(i) for i in range(1,23)]+['X','Y']+['2L','2R','2LHet','2RHet','3L','3R','3LHet','3RHet','4dm','U','Uextra','Xdm']
    if genome=='hg19':
        file="/home/lizj/software/juicer/restriction_sites/hg19_MboI.txt"
        chrs = [str(i) for i in range(1,23)]+['X','Y']
    if genome=='dm3':
        file="/home/chenqw/lizj/fly_analysis/hg19_dm3/dm3_MboI.txt"
        chrs = ['2L','2R','2LHet','2RHet','3L','3R','3LHet','3RHet','4dm','U','Uextra','Xdm']
    if add_chr:
        chrs = ['chr'+str(i) for i in chrs]
    chrindex={}
    for i in range(len(chrs)):
        chrindex[chrs[i]]=i+1
    fi=open(file,'r')
    MboI={}
    for i in fi:
        try:
            MboI[i.strip().split(' ')[0]]=np.array([1]+[int(x) for x in i.strip().split(' ')[1:]])
        except:
            continue
    nchr=[]        
    for k in MboI.keys():
        if k not in chrindex.keys():
            nchr.append(k)        
    for k in nchr:
        MboI.pop(k)
    fi.close()
    return chrindex,MboI