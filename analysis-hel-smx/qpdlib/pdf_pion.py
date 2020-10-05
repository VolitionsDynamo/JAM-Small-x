import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT

#--matplotlib
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#matplotlib.rc('text',usetex=True)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import pylab as py
import matplotlib.gridspec as gridspec

#--from scipy stack 
from scipy.integrate import quad

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core
from analysis.corelib import classifier

FLAV=[]
FLAV.append('g')
FLAV.append('valence')
FLAV.append('sea')

def gen_xf(wdir,Q2=None):
    
    print('\ngenerating pdf-pions from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    if 'pdf-pion' not in conf['steps'][istep]['active distributions']:
        print('pdf-pion not in active distribution')
        return 

    resman=RESMAN(nworkers=1,parallel=False,datasets=False)
    parman=resman.parman

    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    jar_replicas=jar['replicas']
    parman.order=jar['order']

    pdf=conf['pdf-pion']

    #--setup kinematics
    X=10**np.linspace(-3,-1,100)
    X=np.append(X,np.linspace(0.1,0.99,100))
    if Q2==None: Q2=conf['Q20']

    #--compute XF for all replicas        
    XF={}
    cnt=0
    for par in jar_replicas:
        cnt+=1
        lprint('%d/%d'%(cnt,len(replicas)))

        ##--filter
        #flag=False
        #params=replica['params'][istep]
        #order=replica['order'][istep]
        #for i in range(len(order)):
        #    if order[i][0]!=1:continue
        #    if order[i][1]!='pdf':continue
        #    #if order[i][2]=='s1 a':
        #    #   if params[i]<-0.9: flag=True
        #if flag: continue
        core.mod_conf(istep,replicas[cnt-1])

        parman.set_new_params(par,initial=True)

        for flav in FLAV:
            if flav not in XF:  XF[flav]=[]

            if   flav=='valence':
                 func=lambda x: pdf.get_xF(x,Q2,'ub') - pdf.get_xF(x,Q2,'u')
            elif flav=='sea': 
                 func=lambda x: pdf.get_xF(x,Q2,'u') 
            else:
                 func=lambda x: pdf.get_xF(x,Q2,flav) 

            XF[flav].append([func(x) for x in X])
    print     
    checkdir('%s/data'%wdir)
    if Q2==conf['Q20']:
        save({'X':X,'Q2':Q2,'XF':XF},'%s/data/pdf-%d.dat'%(wdir,istep))
    else:
        save({'X':X,'Q2': Q2,'XF':XF},'%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
  

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax=py.subplot(nrows,ncols,1)

    X=data['X']
    for i in range(len(data['XF']['g'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['g'][i]
        ax.plot(X,np.array(f)/10,'r-',alpha=0.1)

    for i in range(len(data['XF']['valence'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['valence'][i]
        ax.plot(X,np.array(f),'g-',alpha=0.1)

    for i in range(len(data['XF']['sea'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['sea'][i]
        ax.plot(X,np.array(f),'b-',alpha=0.1)

    ax.set_xlim(1e-3,1)
    ax.semilogx()
    ax.set_ylim(0,0.6)

    ax.tick_params(axis='both', which='major', labelsize=20)
    #ax.legend(loc=2,bbox_to_anchor=(0.5,1))
    ax.set_xlabel(r'$x_{\pi}$',size=25)
    ax.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=25)
    ax.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6])
    ax.set_xticks([0.001,0.01,0.1,1])
    ax.set_xticklabels([r'$0.001$',r'$0.01$',r'$0.1$',r'$1$'])

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pdf-pion-1-%d.png'%(wdir,istep))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
  
    nrows,ncols=3,1
    fig = py.figure(figsize=(ncols*7,nrows*4))

    X=data['X']
    ax1=py.subplot(nrows,ncols,1)
    for i in range(len(data['XF']['g'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['g'][i]
        ax1.plot(X,np.array(f)/10,'r-',alpha=0.1)
    ax1.plot(X,np.mean(data['XF']['g'],axis=0)/10,'k-')

    ax2=py.subplot(nrows,ncols,2)
    for i in range(len(data['XF']['valence'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['valence'][i]
        ax2.plot(X,np.array(f),'g-',alpha=0.1)
    ax2.plot(X,np.mean(data['XF']['valence'],axis=0),'k-')

    ax3=py.subplot(nrows,ncols,3)
    for i in range(len(data['XF']['sea'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['sea'][i]
        ax3.plot(X,np.array(f),'b-',alpha=0.1)
    ax3.plot(X,np.mean(data['XF']['sea'],axis=0),'k-')

    for ax in [ax1,ax2,ax3]:
        ax.set_xlim(1e-3,1)
        ax.semilogx()
        ax.set_ylim(-0.5,0.5)

        ax.tick_params(axis='both', which='major', labelsize=20)
        #ax.legend(loc=2,bbox_to_anchor=(0.5,1))
        #ax.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6])
        ax.set_xticks([0.001,0.01,0.1,1])
        ax.set_yticks([-0.5,0.0,0.5])
        ax.set_xticklabels([r'$0.001$',r'$0.01$',r'$0.1$',r'$1$'])

    ax1.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=25)
    ax3.set_xlabel(r'$x_{\pi}$',size=25)

    ax1.text(0.1,0.1,r'\boldmath{$\rm glue/10$}', transform=ax1.transAxes,size=25)
    ax2.text(0.1,0.1,r'\boldmath{$\rm valence$}', transform=ax2.transAxes,size=25)
    ax3.text(0.1,0.1,r'\boldmath{$\rm sea$}', transform=ax3.transAxes,size=25)

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pdf-pion-2-%d.png'%(wdir,istep))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
  

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax=py.subplot(nrows,ncols,1)

    X=data['X']
    for i in range(len(data['XF']['g'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['g'][i]
        ax.plot(X,np.array(f)/10,'r-',alpha=0.1)

    for i in range(len(data['XF']['valence'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['valence'][i]
        ax.plot(X,np.array(f),'g-',alpha=0.1)

    for i in range(len(data['XF']['sea'])):
        if cluster[i]!=best_cluster: continue
        f=data['XF']['sea'][i]
        ax.plot(X,np.array(f),'b-',alpha=0.1)

    ax.set_xlim(1e-3,1)

    ax.tick_params(axis='both', which='major', labelsize=20)
    #ax.legend(loc=2,bbox_to_anchor=(0.5,1))
    ax.set_xlabel(r'$x_{\pi}$',size=25)
    ax.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=25)
    ax.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6])

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pdf-pion-2-%d.png'%(wdir,istep))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
 
    nrows,ncols=3,1
    fig = py.figure(figsize=(ncols*7.2,nrows*4))
    #gs = fig.add_gridspec(nrows=3, ncols=2, left=0.05, right=0.48, wspace=0.05)
    gs = gridspec.GridSpec(ncols=2, nrows=3, figure=fig
         ,left=0.17, right=0.95,top=0.99,bottom=0.05 ,wspace=0,hspace=0.05)
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax21 = fig.add_subplot(gs[1,0])
    ax22 = fig.add_subplot(gs[1,1])
    ax31 = fig.add_subplot(gs[2,0])
    ax32 = fig.add_subplot(gs[2,1])

    def plot_pdf(axl,axr,X,Y,ymin,ymax,factor=1,color='r'):
        for i in range(len(Y)):
            #if cluster[i]!=best_cluster: continue
            axl.plot(X,np.array(Y[i])*factor,color=color,ls='-',alpha=0.1)
            axr.plot(X,np.array(Y[i])*factor,color=color,ls='-',alpha=0.1)
        axl.plot(X,np.mean(Y,axis=0)*factor,'k-',lw=3)
        axr.plot(X,np.mean(Y,axis=0)*factor,'k-',lw=3)
        axl.semilogx()
        axl.set_xlim(1e-3,0.1)
        axr.set_xlim(0.1,1)
        axl.set_ylim(ymin,ymax)
        axr.set_ylim(ymin,ymax)
        axr.set_yticklabels([])
        axr.set_yticks([])
        axr.spines['left'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axr.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axl.set_xticks([1e-2,1e-1])
        axr.set_xticks([0.4,0.6,0.8])

    X=data['X']
    ymin,ymax=-0.2,0.5

    Y=data['XF']['g']
    plot_pdf(ax11,ax12,X,Y,ymin,ymax,factor=0.1,color='r')
    Y=data['XF']['valence']
    plot_pdf(ax21,ax22,X,Y,ymin,ymax,factor=1.0,color='g')
    Y=data['XF']['sea']
    plot_pdf(ax31,ax32,X,Y,ymin,ymax,factor=1.0,color='b')


    ax22.plot(X,(1-X),'k:',lw=3,label=r'$(1-x)$')
    ax22.plot(X,(1-X)**2,'k--',lw=3,label=r'$(1-x)^2$')
    ax21.plot(X,(1-X),'k:',lw=3,label=r'$(1-x)$')
    ax21.plot(X,(1-X)**2,'k--',lw=3,label=r'$(1-x)^2$')
    ax21.legend(loc=3,fontsize=20)

    for ax in [ax11,ax12,ax21,ax22]:
        ax.set_xticklabels([])

    ax11.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=30)
    ax11.yaxis.set_label_coords(-0.2, 0.5)
    ax32.set_xlabel(r'$x_{\pi}$',size=30)
    ax32.xaxis.set_label_coords(0.95, -0.05)

    ax31.set_xticklabels([r'$0.01$',r'$0.1$'])

    props = dict(facecolor='white', alpha=1.0)

    ax11.text(0.1,0.8,r'\boldmath{$\rm glue/10$}'
             ,transform=ax11.transAxes,size=25,bbox=props)
    ax21.text(0.1,0.8,r'\boldmath{$\rm valence$}'
             ,transform=ax21.transAxes,size=25,bbox=props)
    ax31.text(0.1,0.8,r'\boldmath{$\rm sea$}'
             ,transform=ax31.transAxes,size=25,bbox=props)

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    if Q2!=None: 
        py.savefig('%s/gallery/pdf-pion-%d.png'%(wdir,istep))
    else: 
        py.savefig('%s/gallery/pdf-pion-Q20-%d.png'%(wdir,istep))

def plot_xf(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
 
    nrows,ncols=3,1
    fig = py.figure(figsize=(ncols*7.2,nrows*4))
    #gs = fig.add_gridspec(nrows=3, ncols=2, left=0.05, right=0.48, wspace=0.05)
    gs = gridspec.GridSpec(ncols=2, nrows=3, figure=fig
         ,left=0.17, right=0.95,top=0.99,bottom=0.05 ,wspace=0,hspace=0.05)
    ax11 = fig.add_subplot(gs[0,0])
    ax12 = fig.add_subplot(gs[0,1])
    ax21 = fig.add_subplot(gs[1,0])
    ax22 = fig.add_subplot(gs[1,1])
    ax31 = fig.add_subplot(gs[2,0])
    ax32 = fig.add_subplot(gs[2,1])

    def plot_pdf(axl,axr,X,Y,ymin,ymax,factor=1,color='r'):

        for i in range(len(Y)):
            #if cluster[i]!=best_cluster: continue
            axl.plot(X,np.array(Y[i])*factor,color=color,ls='-',alpha=0.1)
            axr.plot(X,np.array(Y[i])*factor,color=color,ls='-',alpha=0.1)

        mean=np.mean(Y,axis=0)*factor
        std=np.std(Y,axis=0)*factor

        #axl.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='w'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')
        #axr.fill_between(X,mean-std,mean+std
        #      ,facecolor=None
        #      #,edgecolor='Y'
        #      ,color='Y'
        #      ,alpha=1,zorder=0,hatch='...')

        axl.plot(X,mean+std,'k-',lw=3)
        axl.plot(X,mean-std,'k-',lw=3)
        axr.plot(X,mean+std,'k-',lw=3)
        axr.plot(X,mean-std,'k-',lw=3)

        axl.plot(X,mean,'k--',lw=3)
        axr.plot(X,mean,'k--',lw=3)

        axl.semilogx()
        axl.set_xlim(1e-3,0.1)
        axr.set_xlim(0.1,1)
        axl.set_ylim(ymin,ymax)
        axr.set_ylim(ymin,ymax)
        axr.set_yticklabels([])
        axr.set_yticks([])
        axr.spines['left'].set_visible(False)
        axl.spines['right'].set_visible(False)
        axl.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axr.tick_params(axis='both', which='both',labelsize=20,direction='inout',length=5)
        axl.set_xticks([1e-2,1e-1])
        axr.set_xticks([0.4,0.6,0.8])

    X=data['X']
    ymin,ymax=-0.2,0.5

    Y=data['XF']['g']
    plot_pdf(ax11,ax12,X,Y,ymin,ymax,factor=0.1,color='r')
    Y=data['XF']['valence']
    plot_pdf(ax21,ax22,X,Y,ymin,ymax,factor=1.0,color='g')
    Y=data['XF']['sea']
    plot_pdf(ax31,ax32,X,Y,ymin,ymax,factor=1.0,color='b')


    ax22.plot(X,(1-X),'m:',lw=3,label=r'$(1-x)$')
    ax22.plot(X,(1-X)**2,'m--',lw=3,label=r'$(1-x)^2$')
    ax21.plot(X,(1-X),'m:',lw=3,label=r'$(1-x)$')
    ax21.plot(X,(1-X)**2,'m--',lw=3,label=r'$(1-x)^2$')
    ax22.legend(loc=3,fontsize=20)

    for ax in [ax11,ax12,ax21,ax22]:
        ax.set_xticklabels([])

    ax11.set_ylabel(r'$x_{\pi}f(x_{\pi})$',size=30)
    ax11.yaxis.set_label_coords(-0.2, 0.5)
    ax32.set_xlabel(r'$x_{\pi}$',size=30)
    ax32.xaxis.set_label_coords(0.95, -0.05)

    ax31.set_xticklabels([r'$0.01$',r'$0.1$'])

    props = dict(facecolor='white', alpha=1.0)

    ax11.text(0.1,0.8,r'\boldmath{$\rm glue/10$}'
             ,transform=ax11.transAxes,size=25,bbox=props)
    ax21.text(0.1,0.8,r'\boldmath{$\rm valence$}'
             ,transform=ax21.transAxes,size=25,bbox=props)
    ax31.text(0.1,0.8,r'\boldmath{$\rm sea$}'
             ,transform=ax31.transAxes,size=25,bbox=props)

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    if Q2!=conf['Q20']: 
        py.savefig('%s/gallery/pdf-pion-%d.png'%(wdir,istep))
    else: 
        py.savefig('%s/gallery/pdf-pion-Q20-%d.png'%(wdir,istep))
    py.close()

def gen_moms(wdir,Q2=None):

    print('\ngenerating pion pdf moments from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    if 'pdf-pion' not in conf['steps'][istep]['active distributions']:
        print('pdf-pion not in active distribution')
        return 

    resman=RESMAN(nworkers=1,parallel=False,datasets=False)
    parman=resman.parman

    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    jar_replicas=jar['replicas']
    parman.order=jar['order']

    pdf=conf['pdf-pion']

    #--setup kinematics
    if Q2==None: Q2=conf['Q20']

    #--compute mom for all replicas        
    mom1={}
    mom2={}
    cnt=0
    for par in jar_replicas:
        cnt+=1
        lprint('%d/%d'%(cnt,len(replicas)))
        core.mod_conf(istep,replicas[cnt-1]) #--set conf as specified in istep   
        parman.set_new_params(par,initial=True)

        #print conf['pdf-pion'].params['g1'][0],conf['pdf-pion'].params['g1'][1]

        for flav in FLAV:
            if flav not in mom1:  mom1[flav]=[]
            if flav not in mom2:  mom2[flav]=[]

            if   flav=='valence':
                 func=lambda x: pdf.get_xF(x,Q2,'ub') - pdf.get_xF(x,Q2,'u')
            elif flav=='sea': 
                 func=lambda x: pdf.get_xF(x,Q2,'u') 
            else:
                 func=lambda x: pdf.get_xF(x,Q2,flav) 

            mom1[flav].append(quad(lambda x: func(x),0,1)[0])
            mom2[flav].append(quad(lambda x: x*func(x),0,1)[0])
    print     
    checkdir('%s/data'%wdir)
    
    if Q2==conf['Q20']: fname='%s/data/pdf-moms-%d.dat'%(wdir,istep)
    else: fname='%s/data/pdf-moms-%d-%d.dat'%(wdir,istep,int(Q2))
    save({'Q2':Q2,'mom1':mom1,'mom2':mom2},fname)

def plot_mom1(wdir,kc,Q2=None):

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-moms-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-moms-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]
  

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax=py.subplot(nrows,ncols,1)

    #R=(0,3)
    R=(0,1)
    nbins=50
    mom1=data['mom1']['g']
    mom1=[mom1[i] for i in range(len(mom1)) if cluster[i]==best_cluster]
    ax.hist(mom1,color='r',histtype='step',
            range=R,bins=nbins,label=r'$\rm glue$',density=True)

    mom1=data['mom1']['sea']
    mom1=[6*mom1[i] for i in range(len(mom1)) if cluster[i]==best_cluster]
    ax.hist(mom1,color='b',histtype='step',
            range=R,bins=nbins,label=r'$\rm sea~(\rm tot)$',density=True)

    mom1=data['mom1']['valence']
    mom1=[2*mom1[i] for i in range(len(mom1)) if cluster[i]==best_cluster]
    ax.hist(mom1,color='g',histtype='step',
            range=R,bins=nbins,label=r'$\rm valence~(\rm tot)$',density=True)


    #ax.set_xlim(0,1)
    #ax.set_xlim(0,0.7)
    #ax.set_ylim(0,0.4)

    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.legend(loc=1,fontsize=15)
    ax.set_xlabel(r'$\left<x_{\pi}\right>$',size=25)
    ax.set_ylabel(r'$\rm normalized~yield$',size=25)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    #ax.set_yticks([0.1,0.2,0.3,0.4])

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pdf-mom-pion-%d.png'%(wdir,istep))

    py.close()

def print_mom1(wdir,kc,Q2=None):
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    if Q2==None: Q2=conf['Q20']

    if Q2==conf['Q20']: data=load('%s/data/pdf-moms-%d.dat'%(wdir,istep))
    else: data=load('%s/data/pdf-moms-%d-%d.dat'%(wdir,istep,int(Q2)))

    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]

    best_cluster=cluster_order[0]

    tab={}

    for flav in ['g','sea','valence']:
        
        if    flav=='sea'     : factor=6
        elif  flav=='valence' : factor=2
        else                  : factor=1

        mom1=data['mom1'][flav]
        mom1=[factor*mom1[i] for i in range(len(mom1)) if cluster[i]==best_cluster]
        tab[flav]={}
        tab[flav]['mean']=np.mean(mom1)
        tab[flav]['std']=np.std(mom1)

    print('\n moments from %s\n'%wdir)

    for flav in ['g','sea','valence']:
        msg='%10s = %10.2e +/- %10.2e'
        msg=msg%(flav,tab[flav]['mean'],tab[flav]['std'])
        print(msg)

def plot_pi2n(wdir,kc):

    print('\n plotting  pi2n lambda from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    if 'pdf-pion' not in conf['steps'][istep]['active distributions']:
        print('pdf-pion not in active distribution')
        return 

    Lambda=[]
    for replica in replicas:
        order  = replica['order'][istep]
        params = replica['params'][istep]
        for i in range(len(order)):
            if order[i][1]=='p->pi,n':
                Lambda.append(params[i])

    print np.mean(Lambda)
    print np.std(Lambda)

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax=py.subplot(nrows,ncols,1)

    ax.hist(Lambda,color='r',histtype='step',density=True,bins=100)

    ax.set_xlim(1,2)

    ax.tick_params(axis='both', which='major', labelsize=20)
    #ax.legend(loc=1,fontsize=15)
    ax.set_xlabel(r'$\lambda$',size=25)
    ax.set_ylabel(r'$\rm normalized~yield$',size=25)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
    #ax.set_yticks([0.1,0.2,0.3,0.4])

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/pi2n-%d.png'%(wdir,istep))








