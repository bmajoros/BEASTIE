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
from Pipe import Pipe
import numpy as np
from scipy import stats
from SummaryStats import SummaryStats

NB_P = 0.0294 # 0.06
NB_R = 0.496  # 2

MAX_N=1000
TOLERANCE=0.000001 #0.0001

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

def simNulls(numCases,N1,altFreq):
    cases=[]
    for i in range(numCases):
        alt1=None; ref1=None
        while(True):
            (alt1,ref1)=binom(N1,altFreq)
            if(alt1>0 and ref1>0): break
        while(True):
            n=stats.nbinom.rvs(NB_R,NB_P,size=1)[0]
            if(n<1): continue
            (alt2,ref2)=binom(n,altFreq)
            #if(alt2>0 and ref2>0): continue
            if(alt2>0): continue
            cases.append(Case(alt1,ref1,alt2,ref2,0))
            break
    return cases

def simAlts(numCases,N1,altFreq):
    cases=[]
    for i in range(numCases):
        alt1=None; ref1=None
        while(True):
            (alt1,ref1)=binom(N1,altFreq)
            if(alt1>0 and ref1>0): break
        while(True):
            n=stats.nbinom.rvs(NB_R,NB_P,size=1)[0]
            if(n<1): continue
            alt2=0; ref2=n
            cases.append(Case(alt1,ref1,alt2,ref2,0))
            break
    return cases

def nullPredictor(n,alpha,beta):
    Ns=[]
    i=0
    while(i<1000):
        p=stats.beta.rvs(alpha+1,beta+1,size=1)[0]
        n=stats.nbinom.rvs(NB_R,NB_P,size=1)[0]
        if(n<1): continue
        (alt2,ref2)=binom(n,p)
        if(alt2>0 and ref2>0): continue
        Ns.append(n)
        i+=1
    num=len(Ns)
    numGreater=0
    for x in Ns:
        if(x>=n): numGreater+=1
    return float(numGreater)/float(num)

def runNullPredictor(data):
    P=[]
    for case in data:
        p=nullPredictor(case.alt2+case.ref2,case.alt1,case.ref1)
        P.append(p)
    return P

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
        p=pvalue(case.alt2+case.ref2,case.alt1,case.ref1)
        P.append(p)
    return P

def pvalue(n,alpha,beta):
    s=0
    for i in range(1,n):
        dx=newModel(i,alpha,beta)
        s+=dx
    return 1-s

def pvalue_OLD(n,alpha,beta):
    s=0
    for i in range(n,MAX_N):
        dx=newModel(i,alpha,beta)
        s+=dx
        if(dx<TOLERANCE): break
    return s

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

def getPower(P):
    found=0
    for p in P:
        if(p<0.05): found+=1
    power=float(found)/float(len(P))
    return power

def runJags(data,nMu,nVar):
    P=[]
    for case in data:
        cmd="/hpc/group/majoroslab/BEASTIE/git/run-jags.py /hpc/group/majoroslab/BEASTIE/git/genotype.bugs "+str(case.alt1)+" "+\
            str(case.ref1)+" "+str(case.alt2+case.ref2)+\
            " "+str(nMu)+" "+str(nVar)
        p=float(Pipe.run(cmd))
        P.append(p)
    return P

def empiricalN(data):
    #TEMP=open("raw-counts.txt","wt") ###
    array=[]
    for case in data:
        n=case.alt2+case.ref2
        #print(n,file=TEMP) ###
        array.append(n)
    (mean,SD,Min,Max)=SummaryStats.summaryStats(array)
    #TEMP.close() ###
    return (mean,SD*SD)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <#cases-per-category> <N1> <alt-freq>\n")
(numCases,N1,altFreq)=sys.argv[1:]
numCases=int(numCases)
N1=int(N1)
altFreq=float(altFreq)

# Simulate
nulls=simNulls(numCases,N1,altFreq)
alts=simAlts(numCases,N1,altFreq)
train=simAlts(numCases,N1,altFreq)
(mu,Var)=empiricalN(train)

# Run on nulls
fisherP=fisher(nulls)
betabinP=betabin(nulls)
newModelP=runNewModel(nulls)
jagsP=runJags(nulls,mu,Var)
#nullPredP=runNullPredictor(nulls)

# Run on alts
fisherAltP=fisher(alts)
betabinAltP=betabin(alts)
newModelAltP=runNewModel(alts)
jagsAltP=runJags(alts,mu,Var)

# Compute Type 1 error rate
fisherType1=getType1(fisherP)
betabinType1=getType1(betabinP)
newModelType1=getType1(newModelP)
jagsType1=getType1(jagsP)
#nullPredType1=getType1(nullPredP)

# Compute power
fisherPower=getPower(fisherAltP)
betabinPower=getPower(betabinAltP)
newModelPower=getPower(newModelAltP)
jagsPower=getPower(jagsAltP)

print("Fisher Type I:",fisherType1)
print("Betabin Type I:",betabinType1)
print("New model Type I:",newModelType1)
print("JAGS Type I:",jagsType1)
#print("Null predictor Type I:",nullPredType1)

print("Fisher power:",fisherPower)
print("Betabin power:",betabinPower)
print("New model power:",newModelPower)
print("JAGS power:",jagsPower)




