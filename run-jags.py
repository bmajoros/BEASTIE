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
import TempFilename
from Pipe import Pipe

WARMUP=1000
KEEPERS=10000
SCRIPT_FILE=TempFilename.generate("jags.script")
DATA_FILE=TempFilename.generate("jags.data")
INIT_FILE=TempFilename.generate("jags.init")

def writeScript(scriptFile,modelFile,dataFile,initFile):
    OUT=open(scriptFile,"wt")
    text="model in "+modelFile+"\n"+\
        "data in "+dataFile+"\n"+\
        "compile, nchains(1)\n"+\
        "parameters in "+initFile+", chain(1)\n"+\
        "initialize\n"+\
        "update "+str(WARMUP)+"\n"+\
        "monitor n\n"+\
        "update "+str(KEEPERS)+"\n"+\
        "coda *\n"
    print(text,file=OUT)
    OUT.close()

def writeInitFile(filename):
    OUT=open(filename,"wt")
    text="p <- 0.5\n"+\
        "n <- 5\n"
    print(text,file=OUT)
    OUT.close()

def writeData(filename,alpha,beta):
    OUT=open(filename,"wt")
    text="alpha <- "+str(alpha)+"\n"+\
        "beta <- "+str(beta)+"\n"+\
        "x <- 0"
    print(text,file=OUT)
    OUT.close()

def readSamples(filename):
    samples=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=2): raise Exception("Cannot parse JAGS output")
            (index,value)=fields
            samples.append(value)
    return samples

def getP(samples,n):
    N=len(samples)
    numGreater=0
    for x in samples:
        if(int(x)>=n): numGreater+=1
    return float(numGreater)/float(N)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <model> <alpha> <beta> <n>\n")
(model,alpha,beta,n)=sys.argv[1:]
n=int(n)

writeScript(SCRIPT_FILE,model,DATA_FILE,INIT_FILE)
writeInitFile(INIT_FILE)
writeData(DATA_FILE,alpha,beta)
cmd="jags "+SCRIPT_FILE
Pipe.run(cmd)
samples=readSamples("CODAchain1.txt")
pvalue=getP(samples,n)
print(pvalue)


