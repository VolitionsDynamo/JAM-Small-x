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

def plot_E615_Q(wdir,istep,data):

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*7))
    ax=py.subplot(nrows,ncols,1)            

    Qbins=[]
    #adjust bins
    Qbins.append([4.05,4.5])
    Qbins.append([4.5,4.95])
    Qbins.append([4.95,5.4])
    Qbins.append([5.4,5.85])
    Qbins.append([5.85,6.75])
    Qbins.append([6.75,7.65])
    Qbins.append([7.65,9])
    Qbins.append([9,10.35])
    Qbins.append([10.35,11.7])
    Qbins.append([11.7,13.05])   
    Qbins=Qbins[::-1]

    tab={}
    tab['pT']=data['pT']
    tab['Q']=data['Q']
    tab['value']=data['value']
    tab['alpha']=data['alpha']
    tab['thy']=np.mean(data['prediction-rep'],axis=0)
    tab['dthy']=np.std(data['prediction-rep'],axis=0)
    for k in range(len(data['prediction-rep'])):
        tab['thy%d'%k]=data['prediction-rep'][k]

    tab=pd.DataFrame(tab)
   
    f=40 
    for i in range(len(Qbins)):
        d=tab.query('Q>%f and Q<%f'%(Qbins[i][0],Qbins[i][1]))
        if len(d['Q'])!=0:
            ii=np.argsort(d['pT'].values)
            X=d['pT'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['pT'].values)>1:
                msg=r'$Q\in[%0.1f,%0.1f]$'%(Qbins[i][0],Qbins[i][1])
                p,=ax.plot(X[:-1],(Y*f**i)[:-1],alpha=0.3,label=None)
                #ax.fill_between(X[:-1],(YM*f**i)[:-1],(YP*f**i)[:-1]
                #  ,color=p.get_color(),alpha=0.3,label=None)
                for k in range(len(data['prediction-rep'])):
                #for k in range(70):
                    Yk=d['thy%d'%k].values[ii]
                    #print k, np.sum(((d['value']-Yk)/d['alpha'])**2)/X.size
                    ax.plot(X[:-1],(Yk*f**i)[:-1],ls='-'
                        ,zorder=0
                        ,color=p.get_color(),alpha=0.3) 

                ax.errorbar(d['pT'][:-1],(d['value']*f**i)[:-1],
                       (d['alpha']*f**i)[:-1],zorder=10
                        ,marker='o'
                        ,ls='none'
                        ,markeredgecolor='k'
                        ,ecolor='k'
                        ,color=p.get_color(),label=msg)
            #ax.set_title('Q bin: '+str(Qbins[i][0])+'-'+str(Qbins[i][1]))
            #ax.semilogy() 
            #ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
            #ax.set_xlabel(r'$pT$',size=20)
            #ax.xaxis.set_label_coords(0.95,-0.03)
            #ax.text(0.06,0.06,r'$d^2\sigma/dQdpT$',transform=ax.transAxes,size=20)
            #ax.text(0.78,0.91,r'$\bf E615 pT$',transform=ax.transAxes,size=17)
            #py.tight_layout()
    ax.legend(loc=4)
    ax.semilogy() 
    ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
    ax.set_xlim(None,4.5)
    ax.set_ylim(1e-38,1e-20)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    #ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r'',r''])
    #ax.set_xlabel(r'$pT$',size=20)
    ax.xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(2.5e-1,8e2)
    ax.set_ylabel(r'$d^2\sigma/dQdp_{\rm T}~\times %d^i$'%f,size=20)
    ax.set_xlabel(r'$p_{\rm T}$',size=20)
    ax.text(0.1,0.23,r'$i=0$',transform=ax.transAxes,size=14)
    ax.text(0.5,0.92,r'$i=10$',transform=ax.transAxes,size=14)
    ax.text(0.7,0.9,r'\boldmath{${\rm E615}$}',transform=ax.transAxes,size=17)
    py.tight_layout()
    py.savefig('%s/gallery/dy-pion-E615-Q-%d.png'%(wdir,istep))

def plot_E615_Q(wdir,istep,data):

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*7))
    ax=py.subplot(nrows,ncols,1)            

    Qbins=[]
    #adjust bins
    Qbins.append([4.05,4.5])
    Qbins.append([4.5,4.95])
    Qbins.append([4.95,5.4])
    Qbins.append([5.4,5.85])
    Qbins.append([5.85,6.75])
    Qbins.append([6.75,7.65])
    Qbins.append([7.65,9])
    Qbins.append([9,10.35])
    Qbins.append([10.35,11.7])
    Qbins.append([11.7,13.05])   
    Qbins=Qbins[::-1]

    tab={}
    tab['pT']=data['pT']
    tab['Q']=data['Q']
    tab['value']=data['value']
    tab['alpha']=data['alpha']
    tab['thy']=np.mean(data['prediction-rep'],axis=0)
    tab['dthy']=np.std(data['prediction-rep'],axis=0)
    for k in range(len(data['prediction-rep'])):
        tab['thy%d'%k]=data['prediction-rep'][k]

    tab=pd.DataFrame(tab)
   
    f=1 
    for i in range(len(Qbins)):
        d=tab.query('Q>%f and Q<%f'%(Qbins[i][0],Qbins[i][1]))
        if len(d['Q'])!=0:
            ii=np.argsort(d['pT'].values)
            X=d['pT'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['pT'].values)>1:
                msg=r'$Q\in[%0.1f,%0.1f]$'%(Qbins[i][0],Qbins[i][1])
                p,=ax.plot(X[:-1],(Y*f**i)[:-1],alpha=0.3,label=None)
                ax.fill_between(X[:-1],(YM*f**i)[:-1],(YP*f**i)[:-1]
                  ,color=p.get_color(),alpha=0.3,label=None)
                #for k in range(len(data['prediction-rep'])):
                ##for k in range(70):
                #    Yk=d['thy%d'%k].values[ii]
                #    #print k, np.sum(((d['value']-Yk)/d['alpha'])**2)/X.size
                #    ax.plot(X[:-1],(Yk*f**i)[:-1],ls='-'
                #        ,zorder=0
                #        ,color=p.get_color(),alpha=0.3) 

                ax.errorbar(d['pT'][:-1],(d['value']*f**i)[:-1],
                       (d['alpha']*f**i)[:-1],zorder=10
                        ,marker='o'
                        ,ls='none'
                        ,markeredgecolor='k'
                        ,ecolor='k'
                        ,color=p.get_color(),label=msg)
            #ax.set_title('Q bin: '+str(Qbins[i][0])+'-'+str(Qbins[i][1]))
            #ax.semilogy() 
            #ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
            #ax.set_xlabel(r'$pT$',size=20)
            #ax.xaxis.set_label_coords(0.95,-0.03)
            #ax.text(0.06,0.06,r'$d^2\sigma/dQdpT$',transform=ax.transAxes,size=20)
            #ax.text(0.78,0.91,r'$\bf E615 pT$',transform=ax.transAxes,size=17)
            #py.tight_layout()
    ax.legend(loc=4)
    ax.semilogy() 
    ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
    ax.set_xlim(None,4.5)
    #ax.set_ylim(1e-38,1e-20)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    #ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r'',r''])
    #ax.set_xlabel(r'$pT$',size=20)
    ax.xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(2.5e-1,8e2)
    ax.set_ylabel(r'$d^2\sigma/dQdp_{\rm T}~\times %d^i$'%f,size=20)
    ax.set_xlabel(r'$p_{\rm T}$',size=20)
    ax.text(0.1,0.23,r'$i=0$',transform=ax.transAxes,size=14)
    ax.text(0.5,0.92,r'$i=10$',transform=ax.transAxes,size=14)
    ax.text(0.7,0.9,r'\boldmath{${\rm E615}$}',transform=ax.transAxes,size=17)
    py.tight_layout()
    py.savefig('%s/gallery/dy-pion-E615-Q-%d.png'%(wdir,istep))

def plot_E615_Q(wdir,istep,data):

    Qbins=[]
    #adjust bins
    Qbins.append([4.05,4.5])
    Qbins.append([4.5,4.95])
    Qbins.append([4.95,5.4])
    Qbins.append([5.4,5.85])
    Qbins.append([5.85,6.75])
    Qbins.append([6.75,7.65])
    Qbins.append([7.65,9])
    Qbins.append([9,10.35])
    Qbins.append([10.35,11.7])
    Qbins.append([11.7,13.05])   

    nrows,ncols=6,1
    fig = py.figure(figsize=(ncols*5,nrows*1))
    AX={cnt: py.subplot(nrows,ncols,cnt) for cnt in range(1,7)}


    tab={}
    tab['pT']=data['pT']
    tab['Q']=data['Q']
    tab['value']=data['value']
    tab['alpha']=data['alpha']
    tab['thy']=np.mean(data['prediction-rep'],axis=0)
    tab['dthy']=np.std(data['prediction-rep'],axis=0)
    #for k in range(len(data['prediction-rep'])):
    #    tab['thy%d'%k]=data['prediction-rep'][k]

    tab=pd.DataFrame(tab)
   
    cnt=0
    for i in range(len(Qbins)):
        d=tab.query('Q>%f and Q<%f'%(Qbins[i][0],Qbins[i][1]))
        if len(d['Q'])!=0:
            ii=np.argsort(d['pT'].values)
            X=d['pT'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['pT'].values)>1:
                cnt+=1
                if cnt>sorted(AX.keys())[-1]: continue
                ax=AX[cnt]
                msg=r'$Q\in[%0.1f,%0.1f]$'%(Qbins[i][0],Qbins[i][1])
                ax.text(0.75,0.8,msg,transform=ax.transAxes)
                #p,=ax.plot(X[:-1],(Y*f**i)[:-1],alpha=0.3,label=None)
                ax.fill_between(X[:-1],(YM/Y)[:-1],(YP/Y)[:-1]
                  ,color='y',alpha=0.3,label=None)
                #for k in range(len(data['prediction-rep'])):
                ##for k in range(70):
                #    Yk=d['thy%d'%k].values[ii]
                #    #print k, np.sum(((d['value']-Yk)/d['alpha'])**2)/X.size
                #    ax.plot(X[:-1],(Yk*f**i)[:-1],ls='-'
                #        ,zorder=0
                #        ,color=p.get_color(),alpha=0.3) 

                ax.errorbar(d['pT'][:-1],(d['value']/Y)[:-1],
                       (d['alpha']/Y)[:-1],zorder=10
                        ,fmt='r.')
            #ax.set_title('Q bin: '+str(Qbins[i][0])+'-'+str(Qbins[i][1]))
            #ax.semilogy() 
            #ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
            #ax.set_xlabel(r'$pT$',size=20)
            #ax.xaxis.set_label_coords(0.95,-0.03)
            #ax.text(0.06,0.06,r'$d^2\sigma/dQdpT$',transform=ax.transAxes,size=20)
            #ax.text(0.78,0.91,r'$\bf E615 pT$',transform=ax.transAxes,size=17)
            #py.tight_layout()

    #ax.legend(loc=4)
    #ax.semilogy() 
    #ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
    for _ in AX:
        AX[_].set_xlim(2.0,4.5)
        AX[_].set_ylim(0,3)
        AX[_].axhline(1,color='k',alpha=0.3)
        if _!=6: AX[_].set_xticklabels([])
        AX[_].tick_params(axis='both',which='both',direction='in',labelsize=18)
        AX[_].set_yticks([1,2])

    AX[3].set_ylabel(r'$\rm Ratio~to~theory$',size=30)
    AX[3].yaxis.set_label_coords(-.1,-0.03)
    AX[6].set_xlabel(r'$p_{\rm T}$',size=20)
    AX[6].xaxis.set_label_coords(0.95,-0.03)
    AX[6].set_xticks([2.5,3,3.5,4])


    #AX[3].xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(1e-38,1e-20)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    #ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r'',r''])
    #ax.set_xlabel(r'$pT$',size=20)
    #ax.xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(2.5e-1,8e2)
    #ax.set_ylabel(r'$d^2\sigma/dQdp_{\rm T}~\times %d^i$'%f,size=20)
    #ax.text(0.1,0.23,r'$i=0$',transform=ax.transAxes,size=14)
    #ax.text(0.5,0.92,r'$i=10$',transform=ax.transAxes,size=14)
    AX[6].text(0.7,0.4,r'\boldmath{${\rm E615}$}',transform=AX[6].transAxes,size=17)
    py.tight_layout()
    py.subplots_adjust(wspace=0.1, hspace=0.0)
    py.savefig('%s/gallery/dy-pion-E615-Q-%d.png'%(wdir,istep))

def plot_E615_xF(wdir,istep,data):

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*7))
    ax=py.subplot(nrows,ncols,1)            

    Xbins=[]
    #adjust bins
    Xbins.append([0,0.1])
    Xbins.append([0.1,0.2])
    Xbins.append([0.2,0.3])
    Xbins.append([0.3,0.4])
    Xbins.append([0.4,0.5])
    Xbins.append([0.5,0.6])
    Xbins.append([0.7,0.8])
    Xbins.append([0.8,0.9])
    Xbins.append([0.9,1])
    Xbins=Xbins[::-1]

    tab={}
    tab['pT']=data['pT']
    tab['xF']=data['xF']
    tab['value']=data['value']
    tab['alpha']=data['alpha']
    tab['thy']=np.mean(data['prediction-rep'],axis=0)
    tab['dthy']=np.std(data['prediction-rep'],axis=0)
    for k in range(len(data['prediction-rep'])):
        tab['thy%d'%k]=data['prediction-rep'][k]

    tab=pd.DataFrame(tab)
   
    f=40 
    for i in range(len(Xbins)):
        d=tab.query('xF>%f and xF<%f'%(Xbins[i][0],Xbins[i][1]))
        if len(d['xF'])!=0:
            ii=np.argsort(d['pT'].values)
            X=d['pT'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['pT'].values)>1:
                msg=r'$xF\in[%0.1f,%0.1f]$'%(Xbins[i][0],Xbins[i][1])
                p,=ax.plot(X[:-1],(Y*f**i)[:-1],alpha=0.3,label=None)
                ax.fill_between(X[:-1],(YM*f**i)[:-1],(YP*f**i)[:-1]
                  ,color=p.get_color(),alpha=0.3,label=None)
                #for k in range(len(data['prediction-rep'])):
                ##for k in range(100):
                #    Yk=d['thy%d'%k].values[ii]
                #    #print k, np.sum(((d['value']-Yk)/d['alpha'])**2)/X.size
                #    ax.plot(X[:-1],(Yk*f**i)[:-1],ls='-'
                #        ,zorder=0
                #        ,color=p.get_color(),alpha=0.3) 

                ax.errorbar(d['pT'][:-1],(d['value']*f**i)[:-1],
                       (d['alpha']*f**i)[:-1],zorder=10
                        ,marker='o'
                        ,ls='none'
                        ,markeredgecolor='k'
                        ,ecolor='k'
                        ,color=p.get_color(),label=msg)
            #ax.set_title('Q bin: '+str(Qbins[i][0])+'-'+str(Qbins[i][1]))
            #ax.semilogy() 
            #ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
            #ax.set_xlabel(r'$pT$',size=20)
            #ax.xaxis.set_label_coords(0.95,-0.03)
            #ax.text(0.06,0.06,r'$d^2\sigma/dQdpT$',transform=ax.transAxes,size=20)
            #ax.text(0.78,0.91,r'$\bf E615 pT$',transform=ax.transAxes,size=17)
            #py.tight_layout()
    ax.legend(loc=4)
    ax.semilogy() 
    ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
    ax.set_xlim(None,4.5)
    ax.set_ylim(1e-38,1e-20)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    #ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r'',r''])
    #ax.set_xlabel(r'$pT$',size=20)
    ax.xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(2.5e-1,8e2)
    ax.set_ylabel(r'$d^2\sigma/dxFdp_{\rm T}~\times %d^i$'%f,size=20)
    ax.set_xlabel(r'$p_{\rm T}$',size=20)
    ax.text(0.1,0.23,r'$i=0$',transform=ax.transAxes,size=14)
    ax.text(0.5,0.92,r'$i=10$',transform=ax.transAxes,size=14)
    ax.text(0.7,0.9,r'\boldmath{${\rm E615}$}',transform=ax.transAxes,size=17)
    py.tight_layout()
    py.savefig('%s/gallery/dy-pion-E615-xF-%d.png'%(wdir,istep))

def plot_E615_xF(wdir,istep,data):

    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*5,nrows*7))
    ax=py.subplot(nrows,ncols,1)            

    Xbins=[]
    #adjust bins
    Xbins.append([0,0.1])
    Xbins.append([0.1,0.2])
    Xbins.append([0.2,0.3])
    Xbins.append([0.3,0.4])
    Xbins.append([0.4,0.5])
    Xbins.append([0.5,0.6])
    Xbins.append([0.7,0.8])
    Xbins.append([0.8,0.9])
    Xbins.append([0.9,1])
    Xbins=Xbins[::-1]

    tab={}
    tab['pT']=data['pT']
    tab['xF']=data['xF']
    tab['value']=data['value']
    tab['alpha']=data['alpha']
    tab['thy']=np.mean(data['prediction-rep'],axis=0)
    tab['dthy']=np.std(data['prediction-rep'],axis=0)
    for k in range(len(data['prediction-rep'])):
        tab['thy%d'%k]=data['prediction-rep'][k]

    tab=pd.DataFrame(tab)
   
    f=1#40 
    for i in range(len(Xbins)):
        d=tab.query('xF>%f and xF<%f'%(Xbins[i][0],Xbins[i][1]))
        if len(d['xF'])!=0:
            ii=np.argsort(d['pT'].values)
            X=d['pT'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['pT'].values)>1:
                msg=r'$xF\in[%0.1f,%0.1f]$'%(Xbins[i][0],Xbins[i][1])
                p,=ax.plot(X[:-1],(Y*f**i)[:-1],alpha=0.3,label=None)
                ax.fill_between(X[:-1],(YM*f**i)[:-1],(YP*f**i)[:-1]
                  ,color=p.get_color(),alpha=0.3,label=None)

                ax.errorbar(d['pT'][:-1],(d['value']*f**i)[:-1],
                       (d['alpha']*f**i)[:-1],zorder=10
                        ,marker='o'
                        ,ls='none'
                        ,markeredgecolor='k'
                        ,ecolor='k'
                        ,color=p.get_color(),label=msg)
            #ax.set_title('Q bin: '+str(Qbins[i][0])+'-'+str(Qbins[i][1]))
            #ax.semilogy() 
            #ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
            #ax.set_xlabel(r'$pT$',size=20)
            #ax.xaxis.set_label_coords(0.95,-0.03)
            #ax.text(0.06,0.06,r'$d^2\sigma/dQdpT$',transform=ax.transAxes,size=20)
            #ax.text(0.78,0.91,r'$\bf E615 pT$',transform=ax.transAxes,size=17)
            #py.tight_layout()
    ax.legend(loc=4)
    ax.semilogy() 
    ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
    ax.set_xlim(None,4.5)
    #ax.set_ylim(1e-38,1e-20)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    #ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r'',r''])
    #ax.set_xlabel(r'$pT$',size=20)
    ax.xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(2.5e-1,8e2)
    ax.set_ylabel(r'$d^2\sigma/dxFdp_{\rm T}~\times %d^i$'%f,size=20)
    ax.set_xlabel(r'$p_{\rm T}$',size=20)
    ax.text(0.1,0.23,r'$i=0$',transform=ax.transAxes,size=14)
    ax.text(0.5,0.92,r'$i=10$',transform=ax.transAxes,size=14)
    ax.text(0.7,0.9,r'\boldmath{${\rm E615}$}',transform=ax.transAxes,size=17)
    py.tight_layout()
    py.savefig('%s/gallery/dy-pion-E615-xF-%d.png'%(wdir,istep))

def plot_E615_xF(wdir,istep,data):


    Xbins=[]
    #adjust bins
    Xbins.append([0,0.1])
    Xbins.append([0.1,0.2])
    Xbins.append([0.2,0.3])
    Xbins.append([0.3,0.4])
    Xbins.append([0.4,0.5])
    Xbins.append([0.5,0.6])
    Xbins.append([0.7,0.8])
    Xbins.append([0.8,0.9])
    Xbins.append([0.9,1])
    Xbins=Xbins[::-1]

    tab={}
    tab['pT']=data['pT']
    tab['xF']=data['xF']
    tab['value']=data['value']
    tab['alpha']=data['alpha']
    tab['thy']=np.mean(data['prediction-rep'],axis=0)
    tab['dthy']=np.std(data['prediction-rep'],axis=0)
    #for k in range(len(data['prediction-rep'])):
    #    tab['thy%d'%k]=data['prediction-rep'][k]
    tab=pd.DataFrame(tab)

    nrows,ncols=8,1
    fig = py.figure(figsize=(ncols*5,nrows*1))
    AX={cnt: py.subplot(nrows,ncols,cnt) for cnt in range(1,9)}

    cnt=0 
    for i in range(len(Xbins)):
        d=tab.query('xF>%f and xF<%f'%(Xbins[i][0],Xbins[i][1]))
        if len(d['xF'])!=0:
            ii=np.argsort(d['pT'].values)
            X=d['pT'].values[ii]
            Y=d['thy'].values[ii]
            YP=Y+d['dthy'].values[ii]
            YM=Y-d['dthy'].values[ii]
            if len(d['pT'].values)>1:
                cnt+=1
                ax=AX[cnt]
                ax.fill_between(X[:-1],(YM/Y)[:-1],(YP/Y)[:-1]
                  ,color='y',alpha=0.3,label=None)

                ax.errorbar(d['pT'][:-1],(d['value']/Y)[:-1],
                       (d['alpha']/Y)[:-1],zorder=10
                        ,fmt='r.')
                msg=r'$x_F\in[%0.1f,%0.1f]$'%(Xbins[i][0],Xbins[i][1])
                ax.text(0.75,0.8,msg,transform=ax.transAxes)
            #ax.set_title('Q bin: '+str(Qbins[i][0])+'-'+str(Qbins[i][1]))
            #ax.semilogy() 
            #ax.tick_params(axis='both',which='both',direction='in',pad=4,labelsize=18)
            #ax.set_xlabel(r'$pT$',size=20)
            #ax.xaxis.set_label_coords(0.95,-0.03)
            #ax.text(0.06,0.06,r'$d^2\sigma/dQdpT$',transform=ax.transAxes,size=20)
            #ax.text(0.78,0.91,r'$\bf E615 pT$',transform=ax.transAxes,size=17)
            #py.tight_layout()

    for _ in AX:
        AX[_].set_xlim(2.0,4.5)
        AX[_].axhline(1,color='k',alpha=0.3)
        if _!=8: AX[_].set_xticklabels([])
        AX[_].tick_params(axis='both',which='both',direction='in',labelsize=18)
        #AX[_].set_ylim(0,10)
        #AX[_].set_yticks([1,5])

    AX[4].set_ylabel(r'$\rm Ratio~to~theory$',size=30)
    AX[4].yaxis.set_label_coords(-.1,0.03)
    AX[8].set_xlabel(r'$p_{\rm T}$',size=20)
    AX[8].xaxis.set_label_coords(0.95,-0.03)
    AX[8].set_xticks([2.5,3,3.5,4])

    AX[1].set_ylim(0,11)
    AX[2].set_ylim(0,3)
    AX[3].set_ylim(0,3)
    AX[4].set_ylim(0,3)
    AX[5].set_ylim(0,4)
    AX[6].set_ylim(0,4)
    AX[7].set_ylim(0,6)
    AX[8].set_ylim(0,11)

    AX[1].set_yticks([1,5,10])
    AX[2].set_yticks([1,2])
    AX[3].set_yticks([1,2])
    AX[4].set_yticks([1,2])
    AX[5].set_yticks([1,3])
    AX[6].set_yticks([1,3])
    AX[7].set_yticks([1,5])
    AX[8].set_yticks([1,5,10])

    AX[8].text(0.7,0.2,r'\boldmath{${\rm E615}$}',transform=AX[8].transAxes,size=17)

    #AX[3].xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(1e-38,1e-20)
    #ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    #ax.set_xticklabels([r'$0$',r'',r'$0.2$',r'',r'$0.4$',r'',r'$0.6$',r'',r''])
    #ax.set_xlabel(r'$pT$',size=20)
    #ax.xaxis.set_label_coords(0.95,-0.03)
    #ax.set_ylim(2.5e-1,8e2)
    #ax.set_ylabel(r'$d^2\sigma/dQdp_{\rm T}~\times %d^i$'%f,size=20)
    #ax.text(0.1,0.23,r'$i=0$',transform=ax.transAxes,size=14)
    #ax.text(0.5,0.92,r'$i=10$',transform=ax.transAxes,size=14)
    #ax.text(0.7,0.9,r'\boldmath{${\rm E615}$}',transform=ax.transAxes,size=17)
    py.tight_layout()
    py.subplots_adjust(wspace=0.1, hspace=0.0)
    py.savefig('%s/gallery/dy-pion-E615-xF-%d.png'%(wdir,istep))

def plot_obs(wdir,kc): 

    print('\nplotting DY-pion (qT) from %s'%(wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    predictions=load('%s/data/predictions-%d.dat'%(wdir,istep))
    #print(predictions['reactions'].keys())
    if 'pion_qT' not in predictions['reactions']: return 

    data=predictions['reactions']['pion_qT']
    plot_E615_Q(wdir,istep,data[1001])
    plot_E615_xF(wdir,istep,data[1002])

    return 




