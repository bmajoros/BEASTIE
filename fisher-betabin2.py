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
import math
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

def simNulls_OLD(numCases,N1,N2,altFreq):
    cases=[]
    for i in range(numCases):
        while(True):
            (alt1,ref1)=binom(N1,altFreq)
            (alt2,ref2)=binom(N2,altFreq)
            if(alt2>0 and ref2>0): continue
            cases.append(Case(alt1,ref1,alt2,ref2,0))
            break
    return cases

def simNulls(numCases,N1,N2,altFreq):
    cases=[]
    for i in range(numCases):
        (alt1,ref1)=binom(N1,altFreq)
        while(True):
            n=stats.nbinom.rvs(2,0.06,size=1)[0]
            (alt2,ref2)=binom(n,altFreq)
            if(alt2>0 and ref2>0): continue
            cases.append(Case(alt1,ref1,alt2,ref2,0))
            break
    return cases


def fisher(data):
    P=[]
    for case in data:
        table=[[case.alt1,case.ref1],[case.alt2,case.ref2]]
        (odds,p)=stats.fisher_exact(table, alternative='two-sided')
        # alternative can be two-sided, greater, or less
        P.append(p)
    return P

def betabin(data):
    P=[]
    for case in data:
        p=stats.betabinom.cdf(case.alt2,case.alt2+case.ref2,
                              case.alt1+1,case.ref1+1)
        P.append(p)
    return P

def runNewModel(data):
    P=[]
    for case in data:
        #print(case.alt1,case.ref1,case.alt2,case.ref2)
        p=newModel(case.alt2+case.ref2,case.alt1,case.ref1)
        P.append(p)
    return P

def newModel(n,alpha,beta):
    s=math.log(alpha) if alpha>0 else 0
    for x in range(beta+2,beta+alpha+2):
        s+=math.log(x)
    for x in range(beta+n+1,beta+n+alpha+2):
        s-=math.log(x)
    return math.exp(s)

def getType1(P):
    errors=0
    for p in P:
        if(p<0.05): errors+=1
    rate=float(errors)/float(len(P))
    return rate
        
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

# Run tests
fisherP=fisher(data)
betabinP=betabin(data)
newModelP=runNewModel(data)
print(newModelP)

# Output results
fisherType1=getType1(fisherP)
betabinType1=getType1(betabinP)
newModelType1=getType1(newModelP)
print("Fisher Type I:",fisherType1)
print("Betabin Type I:",betabinType1)
print("New model Type I:",newModelType1)

#FISHER=open(fisherFile,"wt")
#BETABIN=open(betabinFile,"wt")
#for case in data:
    #print(1-case.fisherP,case.cat,sep="\t",file=FISHER)
    #print(1-case.betabinomP,case.cat,sep="\t",file=BETABIN)
    #print(case.cat,case.fisherP,case.betabinomP,sep="\t")
    


