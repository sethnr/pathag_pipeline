import os
import sys
import math as m

def get_value_col(search,filename,fromcol,tocol):
        for line in open(filename):
                t=line.split()
                if t[fromcol]==search:
                        return t['tocol']
        return None

def get_amplicon_panel_id(wildcards):
        ampset=get_value_col(wildcards.sample, "sample","ampset")
	return ampset

def get_subsample_factor(fqsize,maxfqsize=1024):
        print("getting subsamp factor, "+str(fqsize)+" / "+str(maxfqsize), file=sys.stderr) 
        factor=fqsize/maxfqsize
        if factor <= 1: return 1
        else: return m.floor(factor)

