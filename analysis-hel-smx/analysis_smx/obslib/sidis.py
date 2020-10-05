#!/usr/bin/env python
import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
import scipy as sp

## matplotlib
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
from matplotlib.ticker import MultipleLocator, FormatStrFormatter ## for minor ticks in x label
matplotlib.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
matplotlib.rc('text', usetex = True)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import matplotlib.pyplot as py

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
from analysis.corelib import classifier

def plot_mult(wdir, kc, istep, hadron):
    
    print 'Generating %s multiplicity plot'%hadron

    cluster,colors,nc,order = classifier.get_clusters(wdir,istep,kc)
    predictions=load('%s/data/predictions-%d.dat'%(wdir,istep))
    if 'sidis' not in predictions['reactions']: return
    data=predictions['reactions']['sidis']
     
    d_query={}
    d_plot={}
    
    for _ in data:
       del data[_]['prediction-rep']
       del data[_]['residuals-rep']
       d=pd.DataFrame(data[_])
       if d.hadron.values[0]!=hadron: continue
       
       d_query[_]={}
       d_query[_]['x1']={}
       d_query[_]['x1']['y1']=d.query('X>0.02 and X<0.03 and Y>0.10 and Y<0.15')
       d_query[_]['x1']['y2']=d.query('X>0.02 and X<0.03 and Y>0.15 and Y<0.20')
       d_query[_]['x1']['y3']=d.query('X>0.02 and X<0.03 and Y>0.20 and Y<0.30')
       d_query[_]['x1']['y4']=d.query('X>0.02 and X<0.03 and Y>0.30 and Y<0.50')
       d_query[_]['x1']['y5']=d.query('X>0.02 and X<0.03 and Y>0.50 and Y<0.70')
       d_query[_]['x2']={}
       d_query[_]['x2']['y1']=d.query('X>0.03 and X<0.04 and Y>0.10 and Y<0.15')
       d_query[_]['x2']['y2']=d.query('X>0.03 and X<0.04 and Y>0.15 and Y<0.20')
       d_query[_]['x2']['y3']=d.query('X>0.03 and X<0.04 and Y>0.20 and Y<0.30')
       d_query[_]['x2']['y4']=d.query('X>0.03 and X<0.04 and Y>0.30 and Y<0.50')
       d_query[_]['x2']['y5']=d.query('X>0.03 and X<0.04 and Y>0.50 and Y<0.70')
       d_query[_]['x3']={}
       d_query[_]['x3']['y1']=d.query('X>0.04 and X<0.06 and Y>0.10 and Y<0.15')
       d_query[_]['x3']['y2']=d.query('X>0.04 and X<0.06 and Y>0.15 and Y<0.20')
       d_query[_]['x3']['y3']=d.query('X>0.04 and X<0.06 and Y>0.20 and Y<0.30')
       d_query[_]['x3']['y4']=d.query('X>0.04 and X<0.06 and Y>0.30 and Y<0.50')
       d_query[_]['x3']['y5']=d.query('X>0.04 and X<0.06 and Y>0.50 and Y<0.70')
       d_query[_]['x4']={}
       d_query[_]['x4']['y1']=d.query('X>0.06 and X<0.1 and Y>0.10 and Y<0.15')
       d_query[_]['x4']['y2']=d.query('X>0.06 and X<0.1 and Y>0.15 and Y<0.20')
       d_query[_]['x4']['y3']=d.query('X>0.06 and X<0.1 and Y>0.20 and Y<0.30')
       d_query[_]['x4']['y4']=d.query('X>0.06 and X<0.1 and Y>0.30 and Y<0.50')
       d_query[_]['x4']['y5']=d.query('X>0.06 and X<0.1 and Y>0.50 and Y<0.70')
       d_query[_]['x5']={}
       d_query[_]['x5']['y1']=d.query('X>0.1 and X<0.14 and Y>0.10 and Y<0.15')
       d_query[_]['x5']['y2']=d.query('X>0.1 and X<0.14 and Y>0.15 and Y<0.20')
       d_query[_]['x5']['y3']=d.query('X>0.1 and X<0.14 and Y>0.20 and Y<0.30')
       d_query[_]['x5']['y4']=d.query('X>0.1 and X<0.14 and Y>0.30 and Y<0.50')
       d_query[_]['x5']['y5']=d.query('X>0.1 and X<0.14 and Y>0.50 and Y<0.70')
       d_query[_]['x6']={}
       d_query[_]['x6']['y1']=d.query('X>0.14 and X<0.18 and Y>0.10 and Y<0.15')
       d_query[_]['x6']['y2']=d.query('X>0.14 and X<0.18 and Y>0.15 and Y<0.20')
       d_query[_]['x6']['y3']=d.query('X>0.14 and X<0.18 and Y>0.20 and Y<0.30')
       d_query[_]['x6']['y4']=d.query('X>0.14 and X<0.18 and Y>0.30 and Y<0.50')
       d_query[_]['x6']['y5']=d.query('X>0.14 and X<0.18 and Y>0.50 and Y<0.70')
       d_query[_]['x7']={}
       d_query[_]['x7']['y1']=d.query('X>0.18 and X<0.4 and Y>0.10 and Y<0.15')
       d_query[_]['x7']['y2']=d.query('X>0.18 and X<0.4 and Y>0.15 and Y<0.20')
       d_query[_]['x7']['y3']=d.query('X>0.18 and X<0.4 and Y>0.20 and Y<0.30')
       d_query[_]['x7']['y4']=d.query('X>0.18 and X<0.4 and Y>0.30 and Y<0.50')
       d_query[_]['x7']['y5']=d.query('X>0.18 and X<0.4 and Y>0.50 and Y<0.70')

       d_plot[_]={}
       d_plot[_]['theory'] = {}
       d_plot[_]['Z'] = {}
       d_plot[_]['value'] = {}
       d_plot[_]['alpha'] = {}
       
       for x in d_query[_]:
           for i in d_plot[_]: d_plot[_][i][x]={}
           for y in d_query[_][x]:
               d_plot[_]['theory'][x][y] = d_query[_][x][y]['thy']
               d_plot[_]['Z'][x][y] = d_query[_][x][y]['Z']
               d_plot[_]['value'][x][y] = d_query[_][x][y]['value']
               d_plot[_]['alpha'][x][y] = d_query[_][x][y]['alpha']

       nrows, ncols = 2, 4
       fig,axs=py.subplots(nrows,ncols,sharex='all',sharey='all',figsize=(ncols*4.0,nrows*5.0))
       fig.subplots_adjust(hspace=0,wspace=0)       
       axs[1,3].remove()

       if hadron=='pi+': hadsym='\pi^+'
       elif hadron=='pi-': hadsym='\pi^-'
       elif hadron=='K+': hadsym='K^+'
       elif hadron=='K-': hadsym='K^-'
       elif hadron=='h+': hadsym='h^+'
       elif hadron=='h-': hadsym='h^-'
       alpha=[0.0,0.25,0.5,0.75,1.0]
       colors=['red','orange','green','cyan','magenta']
       xbins=[r'$0.02<x<0.03$',r'$0.03<x<0.04$',r'$0.04<x<0.06$',r'$0.06<x<0.1$',r'$0.1<x<0.14$',r'$0.14<x<0.18$',r'$0.18<x<0.4$']
       for i in range(nrows):
           for j in range(ncols):
               x='x%d'%(4*i+j+1)
               if x in d_query[_]:
                   ycount=0
                   for y in sorted(d_query[_][x]):
                       axs[i,j].plot(d_plot[_]['Z'][x][y],d_plot[_]['theory'][x][y]+alpha[ycount],color=colors[ycount],linestyle='solid')
                       axs[i,j].errorbar(d_plot[_]['Z'][x][y],d_plot[_]['value'][x][y]+alpha[ycount],d_plot[_]['alpha'][x][y],color=colors[ycount],marker='o',linestyle='none')
                       ycount+=1           
                   if hadron=='pi+' or hadron=='pi-': ymax=3.8
                   elif hadron=='K+' or hadron=='K-': ymax=1.8
                   elif hadron=='h+' or hadron=='h-': ymax=4.6
                   if i==0 and j==0: 
                       axs[i,j].set_ylim(-0.3,ymax)
                       axs[i,j].text(0.75,0.88*ymax,r'$%s$'%hadsym,fontsize=30)
                       axs[i,j].set_xlim(0.15,0.9)
                       axs[i,j].set_ylabel(r'$\frac{dM^{%s}}{dz}+\alpha$'%hadsym,fontsize=30,fontweight='bold')
                   if i==1 and j==2:
                       axs[i,j].set_xlabel(r'$z$',fontsize=30)
                       axs[i,j].text(1.0,0.5*ymax,r'$0.10<y<0.15,\alpha=0.00$',fontsize=20,color='red')
                       axs[i,j].text(1.0,0.6*ymax,r'$0.15<y<0.20,\alpha=0.25$',fontsize=20,color='orange')
                       axs[i,j].text(1.0,0.7*ymax,r'$0.20<y<0.30,\alpha=0.50$',fontsize=20,color='green')
                       axs[i,j].text(1.0,0.8*ymax,r'$0.30<y<0.50,\alpha=0.75$',fontsize=20,color='cyan')
                       axs[i,j].text(1.0,0.9*ymax,r'$0.50<y<0.70,\alpha=1.00$',fontsize=20,color='magenta')
                   if j==0:
                       if hadron=='pi+' or hadron=='pi-': 
                           axs[i,j].set_yticks([0,1,2,3],minor=False)
                           axs[i,j].set_yticks([-0.2,0.2,0.4,0.6,0.8,1.2,1.4,1.6,1.8,2.2,2.4,2.6,2.8,3.2,3.4,3.6],minor=True)
                           axs[i,j].set_yticklabels(['0','1','2','3'],fontsize=20,minor=False)
                       elif hadron=='K+' or hadron=='K-':
                           axs[i,j].set_yticks([0,1],minor=False)
                           axs[i,j].set_yticks([-0.2,0.2,0.4,0.6,0.8,1.2,1.4,1.6],minor=True)
                           axs[i,j].set_yticklabels(['0','1'],fontsize=20,minor=False)
                       elif hadron=='h+' or hadron=='h-':
                           axs[i,j].set_yticks([0,1,2,3,4],minor=False)
                           axs[i,j].set_yticks([-0.2,0.2,0.4,0.6,0.8,1.2,1.4,1.6,1.8,2.2,2.4,2.6,2.8,3.2,3.4,3.6,3.8,4.2,4.4],minor=True)
                           axs[i,j].set_yticklabels(['0','1','2','3','4'],fontsize=20,minor=False)
                   if (i==1 and j!=3) or (i==0 and j==3):
                       axs[i,j].tick_params(axis='x',which='major',labelbottom=True)
                       axs[i,j].set_xticks([0.2,0.4,0.6,0.8],minor=False)
                       axs[i,j].set_xticks([0.25,0.3,0.35,0.45,0.5,0.55,0.65,0.7,0.75,0.85],minor=True)
                       axs[i,j].set_xticklabels(['0.2','0.4','0.6','0.8'],fontsize=20,minor=False)
                   axs[i,j].tick_params(which='both',direction='in')
                   axs[i,j].axhline(y=0,color='k',linestyle='--')
                   axs[i,j].text(0.2,0.92*ymax,xbins[4*i+j],fontsize=20)
       
       py.tight_layout(h_pad=0,w_pad=0)
       py.savefig('%s/gallery/SIDIS-%d-%s-%s.png' % (wdir, istep, hadron, _))
       py.close()

