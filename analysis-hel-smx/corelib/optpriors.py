#!/usr/bin/env python
import os,sys
import subprocess
import numpy as np
import scipy as sp
import pandas as pd
import copy

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config
from tools.randomstr import id_generator

#--from local
import core,classifier


def gen_priors(wdir,kc,nsamples):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    #cluster,colors,nc,cluster_order = core.get_clusters(wdir,istep) 
    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
  
    #--get data
    replicas=core.get_replicas(wdir)
    step=conf['steps'][istep]
    samples=[]
    for i in range(len(replicas)):

        if cluster[i]!=best_cluster: continue

        replica=replicas[i]
        order=replica['order'][istep]
        params=replica['params'][istep]
        samples.append(params)

    samples=np.array(samples)
    pmin=[np.amin(p) for p in samples.T]
    pmax=[np.amax(p) for p in samples.T]

    new=[]
    for i in range(len(pmin)):
        new.append(np.random.uniform(low=pmin[i], high=pmax[i], size=nsamples))
    new=np.transpose(new)

    checkdir('%s/msr-opt-priors'%wdir)
    replica=replicas[0] #--template replica
    for i in range(nsamples):
        lprint('%d/%d'%(i+1,nsamples))
        replica['params'][istep]=new[i]
        fname='%s/msr-opt-priors/%s.msr'%(wdir,id_generator(12))
        save(replica,fname)
    print 


