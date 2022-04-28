#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
import numpy as np
from scipy import stats

class Case:
    def __init__(self,alt1,ref1,alt2,ref2,cat):
        self.alt1=alt1
        self.alt2=alt2
        self.ref1=ref1
        self.ref2=ref2
        self.cat=cat

def binom(N,p):
    alt=np.random.binomial(N,p)
    ref=N-alt
    return (alt,ref)

def simNulls(numCases,N1,N2,altFreq):
    cases=[]
    for i in range(numCases):
        (alt1,ref1)=binom(N1,altFreq)
        (alt2,ref2)=binom(N2,altFreq)
        cases.append(Case(alt1,ref1,alt2,ref2,0))
    return cases

def simAlts(numCases,N1,N2,altFreq):
    cases=[]
    for i in range(numCases):
        (alt1,ref1)=binom(N1,altFreq)
        alt2=0; ref2=N2 # genotyping error!
        cases.append(Case(alt1,ref1,alt2,ref2,1))
    return cases

def fisher(data):
    for case in data:
        table=[[case.alt1,case.ref1],[case.alt2,case.ref2]]
        (odds,p)=stats.fisher_exact(table, alternative='two-sided')
        # alternative can be two-sided, greater, or less
        case.fisherP=p

def betabin(data):
    for case in data:
        p=stats.betabinom.cdf(case.alt2,case.alt2+case.ref2,
                              case.alt1+1,case.alt2+1)
        case.betabinomP=p
        
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=7):
    exit(ProgramName.get()+" <#cases-per-category> <N1> <N2> <alt-freq> <fisher-out> <betabin-out>\n")
(numCases,N1,N2,altFreq,fisherFile,betabinFile)=sys.argv[1:]
numCases=int(numCases)
N1=int(N1); N2=int(N2)
altFreq=float(altFreq)

# Simulate
data=simNulls(numCases,N1,N2,altFreq)
data.extend(simAlts(numCases,N1,N2,altFreq))

# Run tests
fisher(data)
betabin(data)

# Output results
FISHER=open(fisherFile,"wt")
BETABIN=open(betabinFile,"wt")
for case in data:
    print(1-case.fisherP,case.cat,sep="\t",file=FISHER)
    print(1-case.betabinomP,case.cat,sep="\t",file=BETABIN)
    #print(case.cat,case.fisherP,case.betabinomP,sep="\t")
    


