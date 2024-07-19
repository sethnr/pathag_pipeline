#!/usr/bin/env python

import allel
import numpy as np

import re
import sys
import argparse 
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--bed', '-b', help='primer / amplicon positions (bed)')
parser.add_argument('--depth', '-d', help='depth file (samtools / txt)')
parser.add_argument('--sample', '-s', help='sample name')
parser.add_argument('--factor', '-F', help='subsampling factor')
parser.add_argument('--out', '-o', help='outfile prefix')

args = parser.parse_args()
bedfile = args.bed
depthfile = args.depth
sample = args.sample
outfile = args.out


depth = pd.read_csv(depthfile, sep='\t')
ranges=pd.read_csv(bedfile, sep='\t')

nranges=len(ranges.iloc[:,1])

meandepth=pd.DataFrame({"sample":[sample]*nranges,
              "name":ranges.iloc[:,4],
              "start":ranges.iloc[:,1],
              "end":ranges.iloc[:,2],
              "depth":[-1]*nranges
              })

# meandepth=pd.concat([
#                     [sample]*nranges,
#                     ranges.iloc[4,:],
#                     ranges.iloc[1:2,:],
#                     [-1]*nranges
#                     ],axis=1)

for i in range(0,len(ranges)):
    st = int(ranges.iloc[i,1])
    en = int(ranges.iloc[i,2])
    name = ranges.iloc[i,3]
    indepth = np.logical_and((depth['POS'] >= st),(depth['POS'] <= en))

    if sum(indepth)>0:
        locdepth = round(depth[indepth].iloc[:,2].mean().tolist())
    else:
        locdepth = -1
    meandepth.iloc[i,:] = [sample,name,st,en,locdepth]


meandepth.to_csv(outfile,sep="\t",index=False,)