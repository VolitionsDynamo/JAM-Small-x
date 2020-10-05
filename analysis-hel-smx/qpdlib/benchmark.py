#!/usr/bin/env python
import sys,os
import numpy as np
import copy

#--matplotlib
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('text',usetex=True)
from matplotlib.lines import Line2D
import pylab as py

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core

import tolhapdf
try:
  import lhapdf
except:
  pass

def gen_grid(dist):
    if    dist.startswith('pdf') : X,Q2=tolhapdf.gen_cj_grid()
    elif  dist.startswith('ppdf'): X,Q2=tolhapdf.gen_cj_grid()
    #elif  dist.startswith('ff')  : X,Q2=tolhapdf.gen_dss_grid()
    elif  dist.startswith('ff')  : X,Q2=tolhapdf.gen_cj_grid()

    xmin,xmax=X[0],X[-1]
    Q2min,Q2max=Q2[0],Q2[-1]
    nQ2=len(Q2)

    X=10**np.linspace(np.log10(xmin),-1,10)
    X=np.append(X,np.linspace(0.11,xmax,10))
    #Q2=10**np.linspace(np.log10(Q2min),np.log10(Q2max),nQ2)
    Q2=[Q2min,10,100,1000,10000]

    return X,Q2

def gen_dist_from_fitpack(wdir,dist):

    print('\ngenerating  %s tables for benchmark using %s'%(dist,wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    core.mod_conf(istep) #--set conf as specified in istep   

    resman=RESMAN(nworkers=1,parallel=False,datasets=False)
    parman=resman.parman

    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas=jar['replicas']

    X,Q2=gen_grid(dist)
    qpd=conf[dist]


    nx=len(X)
    nQ2=len(Q2)

    flavs=[-5,-4,-3,-2,-1,1,2,3,4,5,21]
    fmap={-5:'bb',-4:'cb',-3:'sb',-2:'ub',-1:'db'}
    fmap.update({5:'b',4:'c',3:'s',2:'u',1:'d'})
    fmap.update({21:'g'})


    #--gen qpd data per replica
    checkdir('%s/benchmark'%wdir)
    cnt=0
    for par in replicas:
        lprint('progress: %d/%d'%(cnt+1,len(replicas)))
        parman.set_new_params(par)

        data={}
        for _ in Q2:
            data[_]={}
            for flav in flavs:
                data[_][flav]=np.array([qpd.get_xF(x,_,fmap[flav]) for  x in X])
        save(data,'%s/benchmark/%s-%d.dat'%(wdir,dist,cnt))
        cnt+=1
    print

def gen_dist_from_lhapdf(wdir,dist,dirname):

    print('\ngenerating  %s tables for benchmark using %s from lhapdf'%(dist,wdir))

    replicas=lhapdf.mkPDFs(dirname)
    X,Q2=gen_grid(dist)
    flavs=[-5,-4,-3,-2,-1,1,2,3,4,5,21]

    #--gen qpd data per replica
    checkdir('%s/benchmark'%wdir)
    cnt=0
    for i in range(len(replicas)):
        lprint('progress: %d/%d'%(cnt+1,len(replicas)))
        replica=replicas[i]
        data={}
        for _ in Q2:
            data[_]={}
            for flav in flavs:
                data[_][flav]=np.array([replica.xfxQ2(flav,x,_) for  x in X])
        save(data,'%s/benchmark/%s-%d-lha.dat'%(wdir,dist,cnt))
        cnt+=1
    print

def compare(wdir,dist,title,xscale='log'):

    F=os.listdir('%s/benchmark/'%(wdir))
    true = [_ for _ in F if '-lha' not in _ if dist in _]
    lha  = [_ for _ in F if '-lha' in _ if dist in _]
    itrue=np.argsort([int(_.split('-')[-1].split('.')[0]) for _ in true])
    true=[true[i] for i in itrue]
    ilha=np.argsort([int(_.replace('-lha','').split('-')[-1].split('.')[0]) for _ in lha])
    lha=[lha[i] for i in ilha]

    X,Q2=gen_grid(dist)
    flavs=[21,2,1,3,-2,-1,-3,4,5]
    fmap={}
    fmap[21] =r'$g$'
    fmap[ 2] =r'$u$'
    fmap[ 1] =r'$d$'
    fmap[ 3] =r'$s$'
    fmap[ 4] =r'$c$'
    fmap[ 5] =r'$b$'
    fmap[-2] =r'$\bar{u}$'
    fmap[-1] =r'$\bar{d}$'
    fmap[-3] =r'$\bar{s}$'

    nrows,ncols=3,3
    fig = py.figure(figsize=(ncols*5,nrows*5))

    AX={}
    for i in range(9): AX[flavs[i]]=py.subplot(nrows,ncols,i+1)

    measure=lambda A,B: np.abs((A-B)/A)*100

    nrep=50#len(true)
    handler={}
    flag=False
    for i in range(nrep):
        lprint('progress %d/%d'%(i+1,nrep))
        t=load('%s/benchmark/%s'%(wdir,true[i]))
        h=load('%s/benchmark/%s'%(wdir,lha[i]))
        if flag==False:
            for i in range(9):
                for Q2 in t:
                    flav=flavs[i]
                    handler[Q2],=AX[flav].plot(X,measure(h[Q2][flav],t[Q2][flav]) )
            flag=True
        else:
            for i in range(9):
                for Q2 in t:
                    flav=flavs[i]
                    c=handler[Q2].get_color()
                    AX[flav].plot(X,measure(h[Q2][flav],t[Q2][flav]),color=c)
    print 
    AX[4].legend([Line2D([0,1], [0,0], linewidth=10
                 , linestyle='-', color=handler[Q2].get_color()) 
                   for Q2 in sorted(handler)]\
                 ,[r'$Q^2=%0.1f$'%Q2 for Q2 in sorted(handler)]\
                 ,fontsize=30,loc=2,bbox_to_anchor=(-1.1,-0.25),ncol=3)
    AX[3].set_ylabel(r'$\rm rel.~err.(\%)$',size=40)
    AX[4].set_xlabel(r'$x$',size=40)

    AX[21].set_title(r'$\rm %s$'%(title).replace(' ',r'~'),size=40)

    for _ in AX:
        if xscale=='log':
            AX[_].semilogx()
            AX[_].set_xticks([1e-5,1e-3,1e-1])
        #else:
        #    AX[_].set_xticks([0.1])

        AX[_].semilogy()
        AX[_].text(0.1,0.8,r'$%s$'%fmap[_],transform=AX[_].transAxes,size=30)
        AX[_].set_ylim(1e-4,1)
        AX[_].set_xlim(np.amin(X),1)
        AX[_].tick_params(axis='both', which='major', labelsize=20) 
        AX[_].set_yticks([1e-3,1e-2,1e-1])
        AX[_].axhline(1e-1,ls=':',color='k')


    checkdir('%s/gallery'%wdir)
    #py.tight_layout()
    py.subplots_adjust(left=0.1
                      ,bottom=0.2
                      ,right=0.95
                      ,top=0.92
                      ,wspace=None
                      ,hspace=None)
    py.savefig('%s/gallery/benchmark-%s.jpg'%(wdir,dist))
    #py.savefig('%s/gallery/benchmark-%s.pdf'%(wdir,dist))
    py.close()

def compare2(wdir,dist,title,xscale='log'):

    F=os.listdir('%s/benchmark/'%(wdir))
    true = [_ for _ in F if '-lha' not in _ if dist in _]
    lha  = [_ for _ in F if '-lha' in _ if dist in _]
    itrue=np.argsort([int(_.split('-')[-1].split('.')[0]) for _ in true])
    true=[true[i] for i in itrue]
    ilha=np.argsort([int(_.replace('-lha','').split('-')[-1].split('.')[0]) for _ in lha])
    lha=[lha[i] for i in ilha]

    X,Q2=gen_grid(dist)
    flavs=[21,2,1,3,-2,-1,-3,4,5]
    fmap={}
    fmap[21] =r'$g$'
    fmap[ 2] =r'$u$'
    fmap[ 1] =r'$d$'
    fmap[ 3] =r'$s$'
    fmap[ 4] =r'$c$'
    fmap[ 5] =r'$b$'
    fmap[-2] =r'$\bar{u}$'
    fmap[-1] =r'$\bar{d}$'
    fmap[-3] =r'$\bar{s}$'

    nrows,ncols=3,3
    fig = py.figure(figsize=(ncols*5,nrows*5))

    AX={}
    for i in range(9): AX[flavs[i]]=py.subplot(nrows,ncols,i+1)

    measure=lambda A,B: np.abs((A-B)/A)*100
    #measure=lambda A,B: (B/A)

    nrep=50#len(true)
    handler={}
    flag=False
    for i in range(nrep):
        lprint('progress %d/%d'%(i+1,nrep))
        t=load('%s/benchmark/%s'%(wdir,true[i]))
        h=load('%s/benchmark/%s'%(wdir,lha[i]))
        if flag==False:
            for i in range(9):
                for Q2 in sorted(t):
                    flav=flavs[i]
                    handler[Q2],=AX[flav].plot(X,t[Q2][flav])#,alpha=0.1)
                    #AX[flav].plot(X,h[Q2][flav],color=handler[Q2].get_color(),ls='--')
                    break
            flag=True
        else:
            for i in range(9):
                for Q2 in sorted(t):
                    flav=flavs[i]
                    c=handler[Q2].get_color()
                    AX[flav].plot(X,t[Q2][flav],color=handler[Q2].get_color())#,alpha=0.1)
                    #AX[flav].plot(X,h[Q2][flav],color=handler[Q2].get_color(),ls='--')
                    break
    print 
    #AX[4].legend([Line2D([0,1], [0,0], linewidth=10
    #             , linestyle='-', color=handler[Q2].get_color()) 
    #               for Q2 in sorted(handler)]\
    #             ,[r'$Q^2=%0.1f$'%Q2 for Q2 in sorted(handler)]\
    #             ,fontsize=30,loc=2,bbox_to_anchor=(-1.1,-0.25),ncol=3)
    AX[3].set_ylabel(r'$\rm rel.~err.(\%)$',size=40)
    AX[4].set_xlabel(r'$x$',size=40)

    AX[21].set_title(r'$\rm %s$'%(title).replace(' ',r'~'),size=40)

    for _ in AX:
        if xscale=='log':
            AX[_].semilogx()
            AX[_].set_xticks([1e-5,1e-3,1e-1])
        #else:
        #    AX[_].set_xticks([0.1])

        #AX[_].semilogy()
        AX[_].text(0.1,0.8,r'$%s$'%fmap[_],transform=AX[_].transAxes,size=30)
        #AX[_].set_ylim(1e-4,1)
        AX[_].set_ylim(0,1)
        #AX[_].set_ylim(-2,2)
        AX[_].set_xlim(np.amin(X),1)
        AX[_].tick_params(axis='both', which='major', labelsize=20) 
        #AX[_].set_yticks([1e-3,1e-2,1e-1])
        AX[_].axhline(1e-1,ls=':',color='k')


    checkdir('%s/gallery'%wdir)
    #py.tight_layout()
    py.subplots_adjust(left=0.1
                      ,bottom=0.2
                      ,right=0.95
                      ,top=0.92
                      ,wspace=None
                      ,hspace=None)
    py.savefig('%s/gallery/benchmark2-%s.jpg'%(wdir,dist))
    py.close()






