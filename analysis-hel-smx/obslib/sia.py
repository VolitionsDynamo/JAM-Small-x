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
import pylab as py

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
from analysis.corelib import classifier

def data_over_theory(wdir,istep,kc,had):

        print 'Generating data-over-theory plot'

        cluster,colors,nc,order = classifier.get_clusters(wdir,istep,kc)
        predictions=load('%s/data/predictions-%d.dat'%(wdir,istep))
        if 'sia' not in predictions['reactions']: return
        data=predictions['reactions']['sia']

        for _ in data:

            if _==4003:
                print data[4003]
  
            predictions=copy.copy(data[_]['prediction-rep'])
            del data[_]['prediction-rep']
            del data[_]['residuals-rep']
            for ic in range(nc):
                predictions_ic=[predictions[i] for i in range(len(predictions))
                               if cluster[i]==ic]
                data[_]['thy-%d'%ic]=np.mean(predictions_ic, axis=0)
                data[_]['dthy-%d'%ic]=np.std(predictions_ic, axis=0)**0.5

        npannels=0
        for _ in data:
            d=pd.DataFrame(data[_])
            if d.hadron.values[0]==had: npannels+=1

        ncols=4
        nrows=npannels/ncols
        if npannels%ncols>0: nrows+=1

        fig = py.figure(figsize=(ncols*3,nrows*1.5))
        cnt=0
        AX={}
        for _ in data:
            d=pd.DataFrame(data[_])
            col=d.col.values[0]
            if d.hadron.values[0]!=had: continue

            key=(col,had)
            if key not in AX:
                cnt+=1
                AX[key]=py.subplot(nrows,ncols,cnt)

            for ic in range(nc):
                c=colors[ic]
                if c=='r': zorder=10
                else:      zorder=0
                thy=d['thy-%d'%ic]
                Z=d.z
                value=d.value
                alpha=d.alpha
                hp=AX[key].errorbar(Z,value/thy,alpha/thy,fmt='%s.'%c,zorder=zorder)

            if cnt== 1: key0=(col,had)
            if cnt==17: key1=(col,had)
            AX[key].set_title(col)
        #AX[key0].legend()
        #AX[key1].legend()

        for _ in AX:
           AX[_].set_xticks([0.2,0.4,0.6,0.8])
           AX[_].axhline(1,color='b',ls=':')
           AX[_].set_ylim(0.7,1.3)
           AX[_].set_xlim(0,1)
           AX[_].set_xlabel(r'$z$',size=15)
           AX[_].xaxis.set_label_coords(0.95, -0.05)

        py.tight_layout()
        #py.subplots_adjust(left=0.08, bottom=0.08, right=0.99, top=0.97, wspace=None, hspace=0)
        py.savefig('%s/gallery/data-over-thy-sia-%d-%s.png'%(wdir,istep,had))
