import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd

#--matplotlib
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#matplotlib.rc('text',usetex=True)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import pylab as py

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core
from analysis.corelib import classifier


def plot_H1(wdir,istep,data):

    tab={}
    tab['xpi']=data['xpi']
    tab['Q2']=data['Q2']
    tab['value']=data['value']
    tab['alpha']=data['alpha']
    tab['rQ2']=np.around(data['Q2'],decimals=1)
    tab['ry']=np.around(data['y'],decimals=1)
    tab['thy']=np.mean(data['prediction-rep'],axis=0)
    tab['dthy']=np.std(data['prediction-rep'],axis=0)
    for k in range(len(data['prediction-rep'])):
        tab['thy%d'%k]=data['prediction-rep'][k]
    tab=pd.DataFrame(tab)


    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*7))
    ax=py.subplot(nrows,ncols,1)            

    Q2=sorted(set(tab['rQ2']))
    for i in range(len(Q2)):
        d=tab.query('rQ2==%f'%Q2[i])
        y=sorted(set(tab['ry']))
        for ii in range(len(y)):
            dd=d.query('ry==%f'%y[ii])
            if ii==0: c='r'
            if ii==1: c='b'
            ax.errorbar(dd['xpi'],dd['value']*3**i,dd['alpha']*3**i\
                        ,fmt='%s.'%c,label='y=%0.1f'%y[ii])
            iii=np.argsort(dd['xpi'].values)
            X=dd['xpi'].values[iii]
            Y=dd['thy'].values[iii]
            YP=dd['thy'].values[iii]+dd['dthy'].values[iii]
            YM=dd['thy'].values[iii]-dd['dthy'].values[iii]
            ax.fill_between(X,YM*3**i,YP*3**i,color=c,alpha=0.3)

    ax.semilogx()
    ax.semilogy()
    ax.tick_params(axis='both',which='both',direction='in',pad=5,labelsize=18)
    ax.set_xlim(1e-3,2.8e-1)
    ax.set_xlabel(r'$x_{\pi}$',size=22)
    ax.xaxis.set_label_coords(0.97,-0.03)
    #ax.set_ylim(1e-3,1e2)
    ax.text(0.06,0.83,r'$F_2^{\rm LN(3)}$',transform=ax.transAxes,size=22)
    ax.text(0.08,0.75,r'$(\times 3^i)$',transform=ax.transAxes,size=16)
    ax.text(0.47,0.10,r'$i=0$',transform=ax.transAxes,size=14)
    ax.text(0.54,0.855,r'$i=6$',transform=ax.transAxes,size=14)
    ax.text(0.72,0.350,r'$x_L=0.91$',color='r',transform=ax.transAxes,size=14.5)
    ax.text(0.06,0.625,r'$x_L=0.82$',color='b',transform=ax.transAxes,size=14.5)
    ax.text(0.85,0.05,r'$\bf H1$',transform=ax.transAxes,size=17)

    py.tight_layout()
    py.savefig('%s/gallery/ln-H1-%d.png'%(wdir,istep))

def plot_ZEUS(wdir,istep,data):

    tab={}
    tab['xpi']=data['xpi']
    tab['Q2']=data['Q2']
    tab['value']=data['value']
    tab['alpha']=data['alpha']
    tab['rQ2']=np.around(data['Q2'],decimals=1)
    tab['ry']=np.around(data['y'],decimals=1)
    tab['thy']=np.mean(data['prediction-rep'],axis=0)
    tab['dthy']=np.std(data['prediction-rep'],axis=0)
    for k in range(len(data['prediction-rep'])):
        tab['thy%d'%k]=data['prediction-rep'][k]
    tab=pd.DataFrame(tab)

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*7))
    ax=py.subplot(nrows,ncols,1)            
            
    Q2=sorted(set(tab['rQ2']))
    for i in range(len(Q2)):
        d=tab.query('rQ2==%f'%Q2[i])
        y=sorted(set(tab['ry']))
        for ii in range(len(y)):
            dd=d.query('ry==%f'%y[ii])
            if len(dd.xpi)==1: continue
            if ii==0: c='r'
            if ii==1: c='b'
            ax.errorbar(dd['xpi'],dd['value']*3**i,dd['alpha']*3**i\
                        ,fmt='%s.'%c,label='y=%0.1f'%y[ii])
            iii=np.argsort(dd['xpi'].values)
            X=dd['xpi'].values[iii]
            Y=dd['thy'].values[iii]
            YP=dd['thy'].values[iii]+dd['dthy'].values[iii]
            YM=dd['thy'].values[iii]-dd['dthy'].values[iii]
            ax.fill_between(X,YM*3**i,YP*3**i,color=c,alpha=0.3)
    ax.semilogx()
    ax.semilogy()

    ax.tick_params(axis='both',which='both',direction='in',pad=5,labelsize=18)
    ax.set_xlim(5e-4,0.9)
    ax.set_xlabel(r'$x_{\pi}$',size=22)
    ax.xaxis.set_label_coords(0.97,-0.03)
    ax.text(0.1,0.85,r'$r$',transform=ax.transAxes,size=28)
    ax.text(0.08,0.75,r'$(\times 3^i)$',transform=ax.transAxes,size=16)
    ax.text(0.47,0.11,r'$i=0$',transform=ax.transAxes,size=14)
    ax.text(0.84,0.88,r'$i=5$',transform=ax.transAxes,size=14)
    ax.text(0.71,0.44,r'$x_L=0.94$',color='r',transform=ax.transAxes,size=14.5)
    ax.text(0.06,0.51,r'$x_L=0.85$',color='b',transform=ax.transAxes,size=14.5)
    ax.text(0.74,0.05,r'$\bf ZEUS$',transform=ax.transAxes,size=17)

    py.tight_layout()
    py.savefig('%s/gallery/ln-ZEUS-%d.png'%(wdir,istep))

def plot_obs(wdir,kc): 

    print('\nplotting ln data from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    predictions=load('%s/data/predictions-%d.dat'%(wdir,istep))
    data=predictions['reactions']['ln']
    plot_H1(wdir,istep,data[1000])
    plot_ZEUS(wdir,istep,data[2000])





