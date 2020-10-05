#!/usr/bin/env python
import os,sys
import subprocess
import numpy as np
import scipy as sp
import pandas as pd
import copy

#--matplotlib
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#matplotlib.rc('text',usetex=True)
import pylab as py

#--for biclustering
from sklearn.datasets import make_biclusters
from sklearn.datasets import samples_generator as sg
from sklearn.cluster.bicluster import SpectralBiclustering
from sklearn.cluster import SpectralCoclustering
from sklearn.metrics import consensus_score

#--for cluster finding
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import Birch
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs


#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config
from tools.inputmod  import INPUTMOD
from tools.randomstr import id_generator

#--from local  
import core

def _gen_labels_by_chi2(wdir,istep):

    print('\t by chi2 ...')

    #--load data where residuals are computed
    data=load('%s/data/predictions-%d.dat'%(wdir,istep))
    res=np.array(data['res'])
    #--compute chi2/npts
    chi2dof=[]
    for row in res:
        chi2dof.append(np.sum(row**2)/len(row))
    chi2dof=np.array(chi2dof)
    return chi2dof

def _gen_labels_by_echi2(wdir,istep):

    print('\t by echi2 ...')

    #--load data where residuals are computed
    data=load('%s/data/predictions-%d.dat'%(wdir,istep))
    res=np.array(data['res'])
    #--compute chi2/npts
    echi2=[]
    for row in res:
        echi2.append(0)#np.sum(row**2)/len(row))
    echi2=np.array(echi2)

    nexp=0
    for reaction in data['reactions']:
        for idx in data['reactions'][reaction]:
            nexp+=1
            res=np.array(data['reactions'][reaction][idx]['residuals-rep'])
            for i in range(len(res)):
                echi2[i]+=(np.sum(res[i]**2)/len(res[i]))
    return echi2/nexp

def hook(params,order):
    sample=[]
    for i in range(len(order)):
        if order[i][0]!=1: continue
        #if order[i][1]!='pdf': continue
        #if 'g1'  in order[i][2]: continue 
        #if 'uv1' in order[i][2]: continue 
        #if 'dv1' in order[i][2]: continue 
        #if 'db1' in order[i][2]: continue 
        #if 'ub1' in order[i][2]: continue 
        #if 's1' in order[i][2]: continue 
        #if 'sb1' in order[i][2]: continue 
        sample.append(params[i])
    return sample

def _gen_labels_by_cluster(wdir,istep,nc,hook=None):

    print('\t by kmeans ...')

    replicas=core.get_replicas(wdir)
    samples=[]
    for replica in replicas:
        params = replica['params'][istep]
        order  = replica['order'][istep]

        if hook==None:
            samples.append(params)
        else:
            samples.append(hook(params,order))

    #--affinity propagation
    #af = AffinityPropagation(preference=-50).fit(samples)
    #cluster_centers_indices = af.cluster_centers_indices_
    #return af.labels_

    #--agglomerative
    #clustering = AgglomerativeClustering(nc).fit(samples)
    #return clustering.labels_

    #brc = Birch(branching_factor=3, n_clusters=None, threshold=0.5,compute_labels=True)
    #brc.fit(samples) 
    #return brc.labels_

    #--kmean    
    kmeans = KMeans(n_clusters=nc, random_state=0).fit(samples)
    return kmeans.labels_
   
def gen_labels(wdir,kc):

    print('\ngenerating labels for %s\n'%wdir)

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    if istep not in kc.nc:
        print("istep not in kmeanconf.py")
        sys.exit()

    nc=kc.nc[istep]
    hook=kc.hooks[istep]

    istep=core.get_istep()
    labels={}
    labels['chi2dof']  = _gen_labels_by_chi2(wdir,istep)
    labels['echi2dof'] = _gen_labels_by_echi2(wdir,istep)
    labels['cluster']  = _gen_labels_by_cluster(wdir,istep,nc,hook)

    clusters = sorted(list(set(labels['cluster'])))
    chi_2_clusters = [0.0 for _ in clusters]
    length_clusters = [0.0 for _ in clusters]
    for i in range(len(labels['cluster'])):
        for j in range(len(clusters)):
            if labels['cluster'][i] == clusters[j]:
                chi_2_clusters[j] += labels['chi2dof'][i]
                length_clusters[j] += 1.0
    labels['cluster_average_chi2'] = np.array(chi_2_clusters) / np.array(length_clusters)

    checkdir('%s/data'%wdir)
    save(labels,'%s/data/labels-%d.dat'%(wdir,istep))

def get_clusters(wdir,istep,kc): 

    nc=kc.nc[istep]
    hook=kc.hooks[istep]

    labels   = load('%s/data/labels-%d.dat'%(wdir,istep))
    cluster  = labels['cluster']
    chi2dof  = labels['chi2dof']
    echi2dof = labels['echi2dof']

    #--get clusters idx ordering
    echi2_means = [np.mean([echi2dof[i] for i in range(len(echi2dof))
           if cluster[i]==j ])  for j in range(nc)]

    order=np.argsort(echi2_means)

    clist=['r','c','g','y']
    colors={order[i]:clist[i] for i in range(nc)}
    return cluster,colors,nc,order

#--plots

def plot_chi2_dist(wdir,istep):

    #--needs revisions NS (09/02/19)

    labels   = load('%s/data/labels-%d.dat'%(wdir,istep))
    chi2     = labels['chi2dof']
    cluster,colors,nc,order = get_clusters(wdir,istep) 

    nrows=1
    ncols=1
    #--plot labeled residuals
    ax=py.subplot(nrows,ncols,1)

    chi2=[_ for _ in chi2 if _<2]
    Mean=np.average(chi2)
    dchi=np.std(chi2)
    chi2min=Mean-2*dchi
    chi2max=Mean+2*dchi
    R=(chi2min,chi2max)
    ax.hist(chi2,bins=60,range=R,histtype='step',label='size=%d'%len(chi2))
    for j in range(nc):
        chi2_=[chi2[i] for i in range(len(chi2)) if cluster[i]==j]
        c=colors[j]
        print c,np.mean(chi2_),np.mean(chi2)
        label='cluster id=%d  size=%d'%(j,len(chi2_))
        ax.hist(chi2_,bins=30,range=R,histtype='step',label=label,color=c)
        mean=np.average(chi2_)
        ax.axvline(mean,color=c,ls=':')
        

    ax.set_xlabel('chi2/npts')
    ax.set_ylabel('yiled')
    ax.legend()
    ax.text(0.1,0.8,'step %d'%istep,transform=ax.transAxes)
    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/chi2-dist-%d.jpg'%(wdir,istep))
    py.close()

def plot_chi2_dist_per_exp(wdir,istep):

    #--needs revisions NS (09/02/19)

    labels   = load('%s/data/labels-%d.dat'%(wdir,istep))
    cluster,colors,nc,order = get_clusters(wdir,istep) 
    chi2     = labels['chi2dof']
    data=load('%s/data/predictions-%d.dat'%(wdir,istep))

    for reaction in data['reactions']:
        for idx in data['reactions'][reaction]:
            if reaction!='sidis': continue
            if idx!=2006:       continue
            
            res=np.array(data['reactions'][reaction][idx]['residuals-rep'])
            chi2=np.array([np.sum(r**2)/len(r) for r in res])
            #chi2=np.array([_ for _ in chi2 if _<5])

            Mean=np.average(chi2)
            dchi=np.std(chi2)
            chi2min=Mean-2*dchi
            chi2max=Mean+2*dchi
            R=(chi2min,chi2max)
            nrows=1
            ncols=1
            #--plot labeled residuals
            ax=py.subplot(nrows,ncols,1)

            ax.hist(chi2,bins=30,range=R,histtype='step',label='size=%d'%len(chi2))
            for j in range(nc):
                _chi2=[chi2[i] for i in range(len(chi2)) if cluster[i]==j]
                c=colors[j]
                print c, np.mean(_chi2),reaction,idx
                label='cluster id=%d  size=%d'%(j,len(_chi2))
                ax.hist(_chi2,bins=30,range=R,histtype='step',label=label,color=c)
            ax.set_xlabel('chi2/npts')
            ax.set_ylabel('yiled')
            ax.legend()
            ax.text(0.1,0.8,'step %d'%istep,transform=ax.transAxes)
            py.tight_layout()
            checkdir('%s/gallery'%wdir)
            py.savefig('%s/gallery/chi2-dist-%d-%s-%d.jpg'%(wdir,istep,reaction,idx))
            py.close()





