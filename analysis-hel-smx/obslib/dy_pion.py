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

def plot_E615(wdir,istep,data):

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*7))
    ax=py.subplot(nrows,ncols,1)            

    Q2bins=[]
    Q2bins.append([20,23])
    Q2bins.append([23,28])
    Q2bins.append([28,35])
    Q2bins.append([35,40])
    Q2bins.append([40,45])
    Q2bins.append([45,50])
    Q2bins.append([50,60])
    Q2bins.append([60,65])

    tab={}
    tab['xF']=data['xF']
    tab['Q2']=data['Q2']
    tab['value']=data['value']
    tab['alpha']=data['alpha']
    tab['thy']=np.mean(data['prediction-rep'],axis=0)
    tab['dthy']=np.std(data['prediction-rep'],axis=0)
    for k in range(len(data['prediction-rep'])):
        tab['thy%d'%k]=data['prediction-rep'][k]

    tab=pd.DataFrame(tab)
    
    for i in range(len(Q2bins)):
        d=tab.query('Q2>%f and Q2<%f'%(Q2bins[i][0],Q2bins[i][1]))
        if len(d['Q2'])!=0:
            ax.errorbar(d['xF'][:-1],(d['value']*3**i)[:-1],(d['alpha']*3**i)[:-1],fmt='b.')
            ii=np.argsort(d['xF'].values)
            X=d['xF'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['xF'].values)>1:
                ax.fill_between(X[:-1],(YM*3**i)[:-1],(YP*3**i)[:-1],color='b',alpha=0.3)
                #for k in range(len(data['prediction-rep'])):
                #    Yk=d['thy%d'%k].values[ii]
                #    #print k, np.sum(((d['value']-Yk)/d['alpha'])**2)/X.size
                #    ax.plot(X[:-1],(Yk*3**i)[:-1],'k-',alpha=0.3) 

    ax.semilogy() 
    ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
    ax.set_xlim(0,0.8)
    ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r'',r''])
    ax.set_xlabel(r'$x_F$',size=20)
    ax.xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(2.5e-1,8e2)
    ax.text(0.06,0.06,r'$d^2\sigma/d\sqrt{\tau}dx_F$',transform=ax.transAxes,size=20)
    ax.text(0.56,0.08,r'$(\times 3^i)$',transform=ax.transAxes,size=16)
    ax.text(0.1,0.23,r'$i=0$',transform=ax.transAxes,size=14)
    ax.text(0.5,0.92,r'$i=7$',transform=ax.transAxes,size=14)
    ax.text(0.78,0.91,r'$\bf E615$',transform=ax.transAxes,size=17)
    py.tight_layout()
    py.savefig('%s/gallery/dy-pion-E615-%d.png'%(wdir,istep))

def plot_NA10(wdir,istep,data1,data2):

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*7))
    ax=py.subplot(nrows,ncols,1)            

    Q2bins=[]
    Q2bins.append([15,18])
    Q2bins.append([18,20])
    Q2bins.append([20,23])
    Q2bins.append([23,26])
    Q2bins.append([26,30])
    Q2bins.append([30,33])
    Q2bins.append([33,35])
    Q2bins.append([35,40.8])#
    Q2bins.append([40.8,43])#
    Q2bins.append([43,49])
    Q2bins.append([49,52])
    Q2bins.append([52,58])
    Q2bins.append([58,65])

    tab1={}
    tab1['xF']=data1['xF']
    tab1['Q2']=data1['Q2']
    tab1['value']=data1['value']
    tab1['alpha']=data1['alpha']
    tab1['thy']=np.mean(data1['prediction-rep'],axis=0)
    tab1['dthy']=np.std(data1['prediction-rep'],axis=0)
    for k in range(len(data1['prediction-rep'])):
        tab1['thy%d'%k]=data1['prediction-rep'][k]
    tab1=pd.DataFrame(tab1)

    tab2={}
    tab2['xF']=data2['xF']
    tab2['Q2']=data2['Q2']
    tab2['value']=data2['value']
    tab2['alpha']=data2['alpha']
    tab2['thy']=np.mean(data2['prediction-rep'],axis=0)
    tab2['dthy']=np.std(data2['prediction-rep'],axis=0)
    for k in range(len(data2['prediction-rep'])):
        tab2['thy%d'%k]=data2['prediction-rep'][k]
    tab2=pd.DataFrame(tab2)


    for i in range(len(Q2bins)):
        d=tab1.query('Q2>%f and Q2<%f'%(Q2bins[i][0],Q2bins[i][1]))

        if len(d['Q2'])!=0:
            ax.errorbar(d['xF'],d['value'],d['alpha'],fmt='g.')
            ii=np.argsort(d['xF'].values)
            X=d['xF'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['xF'].values)>1:ax.fill_between(X,YM,YP,color='g',alpha=0.3)
            #for k in range(len(data1['prediction-rep'][:100])):
            #    Yk=d['thy%d'%k].values[ii]
            #    #print k, np.sum(((d['value']-Yk)/d['alpha'])**2)/X.size
            #    ax.plot(X,Yk,'k-',alpha=0.3) 


    for i in range(len(Q2bins)):
        d=tab2.query('Q2>%f and Q2<%f'%(Q2bins[i][0],Q2bins[i][1]))
        if len(d['Q2'])!=0:
            ax.errorbar(d['xF'],d['value'],d['alpha'],fmt='b.')
            ii=np.argsort(d['xF'].values)
            X=d['xF'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['xF'].values)>1:ax.fill_between(X,YM,YP,color='b',alpha=0.3)
            #for k in range(len(data2['prediction-rep'][:100])):
            #    Yk=d['thy%d'%k].values[ii]
            #    #print k, np.sum(((d['value']-Yk)/d['alpha'])**2)/X.size
            #    ax.plot(X,Yk,'k-',alpha=0.3) 



    ax.semilogy() 
    ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
    ax.set_xlim(0,0.7)
    ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])
    ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r''])
    ax.set_xlabel(r'$x_F$',size=20)
    ax.xaxis.set_label_coords(0.97,-0.03)
    ax.set_ylim(0.06,None)
    ax.text(0.1,0.06,r'$d^2\sigma/d\sqrt{\tau}dx_F$',transform=ax.transAxes,size=20)
    ax.text(0.76,0.91,r'$\bf NA10$',transform=ax.transAxes,size=17)
    ax.text(0.38,0.89,r'$286~\rm GeV$',color='b',transform=ax.transAxes,size=14)
    ax.text(0.67,0.76,r'$194~\rm GeV$',color='g',transform=ax.transAxes,size=14)
 
    py.tight_layout()
    py.savefig('%s/gallery/dy-pion-NA10-%d.png'%(wdir,istep))
    py.close()

def plot_obs(wdir,kc): 

    print('\nplotting dy-pion data from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    predictions=load('%s/data/predictions-%d.dat'%(wdir,istep))
    if 'dy-pion' not in predictions['reactions']: return 


    data=predictions['reactions']['dy-pion']
    plot_E615(wdir,istep,data[10001])
    plot_NA10(wdir,istep,data[10002],data[10003])

    return 




