#!/usr/bin/env python
import os,sys
import subprocess
import numpy as np
import scipy as sp
import pandas as pd
import copy

## matplotlib
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
matplotlib.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
matplotlib.rc('text', usetex = True)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import pylab as py
# from fflib.dss.dss   import DSS
# from fflib.akk.akk   import AKK
# from fflib.hkns.hkns import HKNS
# from fflib.kkp.kkp   import KKP

## from fitpack
from tools.tools import load, save, checkdir, lprint
from tools.config import conf, load_config
from analysis.corelib import core
# from pdfcalc  import PDFCALC

class PLOTS:

    def __init__(self, wdir, had, plot):

        self.msg = 'mplots.line_figures'
        if plot==1: self.line_figure1(wdir,had)
        elif plot==2: self.line_figure2(wdir,had)

    def plot_lines(self, ax, flav, c, n_replicas, lab = ''): ## plot all the FF replicas as separate lines
        if n_replicas<=10:
          line_width = 1.0 / float(n_replicas)
        else: line_width = 10.0 / float(n_replicas)
        X = self.data['X']
        nrep = len(self.data['XF'][flav])
        for i in range(nrep):
            if self.cluster[i] != self.best_cluster: continue
            ax.plot(X, self.data['XF'][flav][i], color = c, label = lab, linewidth = line_width)

    def line_figure1(self, wdir, had):
        # self.cluster, colors, nc, order = self.get_clusters(wdir, istep)
        
        print 'Generating ff plot'
        
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        nrows, ncols = 2, 3
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        # AX = [py.subplot(nrows, ncols, cnt) for cnt in range(1, 9)]
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(6)]

        ## FF
        for i in had:
          self.data = load('%s/data/ff%s-%d.dat' % (wdir, i, istep))

          if i=='pion':
            c='r'
            lab=r'$\pi$'
          elif i=='kaon': 
            c='b'
            lab=r'$K$'
          elif i=='hadron':
            c='k'
            lab=r'$h$'

          self.plot_lines(AX[0], 'u+ub', c, n_replicas, lab)
          self.plot_lines(AX[1], 'd+db' , c, n_replicas)
          self.plot_lines(AX[2], 's+sb', c, n_replicas)
          self.plot_lines(AX[3], 'c+cb', c, n_replicas)
          self.plot_lines(AX[4], 'b+bb', c, n_replicas)
          self.plot_lines(AX[5], 'g', c, n_replicas)

        AX[0].text(0.9,0.55,r'$u^+$', color = 'k', size = 20)
        AX[0].set_ylabel(r'$zD(z)$', size = 20)
        AX[0].legend(loc='upper center')
        AX[1].text(0.9,0.55,r'$d^+$', color = 'k', size = 20)
        AX[2].text(0.9,0.55,r'$s^+$', color = 'k', size = 20)
        AX[3].text(0.9,0.55,r'$c^+$', color = 'k', size = 20)
        AX[3].set_ylabel(r'$zD(z)$', size = 20)
        AX[3].set_xlabel(r'$z$', size = 20, x=0.9, va = 'baseline')
        AX[4].text(0.9,0.55,r'$b^+$', color = 'k', size = 20)
        AX[4].set_xlabel(r'$z$', size = 20, x=0.9, va = 'baseline')
        AX[5].text(0.9,0.55,r'$g$', color = 'k', size = 20)
        AX[5].set_xlabel(r'$z$', size = 20, x=0.9, va = 'baseline')

        for ax in AX:
            ax.set_ylim(0.0,0.6)
           # ax.set_yticks([0.2,0.6,1.0,1.4])
            ax.set_xlim(0.1,1.0)
            ax.set_xticks([0.2,0.4,0.6,0.8])

        py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        py.savefig('%s/gallery/ff-line-plots1-%d.png' % (wdir, istep), dpi = 1000)

    def line_figure2(self, wdir, had):
        # self.cluster, colors, nc, order = self.get_clusters(wdir, istep)
        
        print 'Generating ff plot'
        
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        n_replicas = len(self.cluster)
        self.best_cluster = 0

        nrows, ncols = 3, 3
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        # AX = [py.subplot(nrows, ncols, cnt) for cnt in range(1, 9)]
        AX = [py.subplot(nrows, ncols, cnt+1) for cnt in range(9)]

        ## FF
        for i in had:
          self.data = load('%s/data/ff%s-%d.dat' % (wdir, i, istep))

          if i=='pion':
            c='r'
            lab=r'$\pi$'
          elif i=='kaon': 
            c='b'
            lab=r'$K$'
          elif i=='hadron':
            c='k'
            lab=r'$h$'

          self.plot_lines(AX[0], 'u', c, n_replicas, lab)
          self.plot_lines(AX[1], 'd' , c, n_replicas)
          self.plot_lines(AX[2], 's', c, n_replicas)
          self.plot_lines(AX[3], 'ub', c, n_replicas, lab)
          self.plot_lines(AX[4], 'db' , c, n_replicas)
          self.plot_lines(AX[5], 'sb', c, n_replicas)
          self.plot_lines(AX[6], 'c', c, n_replicas)
          self.plot_lines(AX[7], 'b', c, n_replicas)
          self.plot_lines(AX[8], 'g', c, n_replicas)

        AX[0].text(0.82,0.8,r'$u$', color = 'k', size = 20)
        AX[0].set_ylabel(r'$zD(z)$', size = 20)
        AX[0].legend(loc='upper center')
        AX[1].text(0.82,0.8,r'$d$', color = 'k', size = 20)
        AX[2].text(0.82,0.8,r'$s$', color = 'k', size = 20)
        AX[3].text(0.82,0.8,r'$\bar{u}$', color = 'k', size = 20)
        AX[3].set_ylabel(r'$zD(z)$', size = 20)
        AX[4].text(0.82,0.8,r'$\bar{d}$', color = 'k', size = 20)
        AX[5].text(0.82,0.8,r'$\bar{s}$', color = 'k', size = 20)
        AX[6].text(0.82,0.8,r'$c,\bar{c}$', color = 'k', size = 20)
        AX[6].set_ylabel(r'$zD(z)$', size = 20)
        AX[6].set_xlabel(r'$z$', size = 20, x=0.9, va = 'baseline')
        AX[7].text(0.82,0.8,r'$b,\bar{b}$', color = 'k', size = 20)
        AX[7].set_xlabel(r'$z$', size = 20, x=0.9, va = 'baseline')
        AX[8].text(0.82,0.8,r'$g$', color = 'k', size = 20)
        AX[8].set_xlabel(r'$z$', size = 20, x=0.9, va = 'baseline')

        for ax in AX:
            ax.set_ylim(0.0,0.9)
            ax.set_xlim(0.15,0.9)
            ax.set_xticks([0.2,0.4,0.6,0.8])

        py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        py.savefig('%s/gallery/ff-line-plots2-%d.png' % (wdir, istep), dpi = 1000)

if __name__ == '__main__':
    pass
