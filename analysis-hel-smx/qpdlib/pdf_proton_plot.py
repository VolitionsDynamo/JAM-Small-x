#!/usr/bin/env python
import os, sys
import subprocess
import numpy as np
import scipy as sp
import pandas as pd
import copy
import re

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

class PDF_PLOT_CORE:

    def plot_lines(self, ax, flav, c, lab = '', all_cluster = False): ## plot all the PDF replicas as lines
        n_replicas = len(self.xf_data['XF'][flav])
        line_width = 10.0 / float(n_replicas)
        if line_width > 1.0: line_width = 1.0
        xs = self.xf_data['X']

        if all_cluster:
            clusters = list(set(self.cluster))
            colors = ['red', 'gold', 'dodgerblue', 'black']
            count_clusters = [0 for _ in clusters]
            for i in range(n_replicas):
                count_clusters[self.cluster[i]] += 1
                if count_clusters[self.cluster[i]] == 1:
                    ax.plot(xs, self.xf_data['XF'][flav][i], color = colors[self.cluster[i]], linewidth = line_width, \
                            label = r'$\overline{\chi _{%i}^2} = %.2f$' % (self.cluster[i], self.cluster_average_chi2[self.cluster[i]]))
                else:
                    ax.plot(xs, self.xf_data['XF'][flav][i], color = colors[self.cluster[i]], label = lab, linewidth = line_width)
        else:
            for i in range(n_replicas):
                if self.cluster[i] != self.best_cluster: continue
                ax.plot(xs, self.xf_data['XF'][flav][i], color = c, label = lab, linewidth = line_width)

    def plot_band_JAM(self, ax, flav, c, lab, scale = 1, distinction = {'alpha': 0.7, 'zorder': 8}):
        xs = self.xf_data['X']
        n_replicas = len(self.xf_data['XF'][flav])
        Y = []
        for i in range(n_replicas):
            if self.cluster[i] != self.best_cluster: continue
            Y.append(self.xf_data['XF'][flav][i])
        Y0 = np.mean(Y, axis = 0)
        dY = np.std(Y, axis = 0)
        alpha = distinction['alpha']
        zorder = distinction['zorder']
        return ax.fill_between(xs, (Y0 - dY) * scale, (Y0 + dY) * scale, color = c, label = lab, alpha = alpha, zorder = 1)

    def plot_band(self, ax, group, name, flav, alpha = 0.2, label = '', scale = 1):
        xs = group['X']
        D = group['XF'][name][flav]
        zorder = 2
        if name == 'NNPDF':
            alpha = 0.3
            zorder = 2
        y_down, y_up = scale * (D['center'] - D['difference']), scale * (D['center'] + D['difference'])
        return ax.fill_between(xs, y_down, y_up, color = group['XF'][name]['color'], alpha = alpha, label = label, zorder = zorder)

    def plot_ratios(self, ax, flav, c, lab = ''): ## plot all the PDF replicas as separate lines
        n_replicas = len(self.xf_data['XF'][flav])
        line_width = 10.0 / float(n_replicas)
        if line_width > 1.0: line_width = 1.0
        xs = self.xf_data['X']
        for i in range(n_replicas):
            if self.cluster[i] != self.best_cluster: continue
            ax.plot(xs, self.xf_data['XF'][flav][i] / self.base_xf_data['mean_xf'][flav], color = c, label = lab, linewidth = line_width)

    def plot_band_steps(self, ax, step, flav, color, alpha, label = None):
        xs = self.xf_data[step]['X']
        D = self.xf_data[step]['XF'][flav]
        Y = []
        n_replicas = len(self.xf_data[step]['XF'][flav])
        for i in range(n_replicas):
            if self.cluster[step][i] != self.best_cluster[step]: continue
            Y.append(self.xf_data[step]['XF'][flav][i])
        Y0 = np.mean(Y, axis = 0)
        dY = np.std(Y, axis = 0)
        return ax.fill_between(xs, (Y0 - dY), (Y0 + dY), color = color, alpha = alpha, zorder = step, label = label)

    def pdf_histogram(self, ax, flav, parameter, label): ## histogram for PDF normalization
        colors = ['orangered', 'limegreen', 'dodgerblue', 'fuchsia']
        shapes = self.pdf_parameter_data.keys()
        for shape in shapes:
            if flav in self.pdf_parameter_data[shape]:
                ys = self.pdf_parameter_data[shape][flav][parameter]
            else:
                continue
            if all([_ == ys[0] for _ in ys]):
                if shape == 1:
                    ax.text(0.1, 0.9, label.replace('$', r'$\mathrm{fixed}~', 1), fontsize = 17, horizontalalignment = 'left', verticalalignment = 'top', transform = ax.transAxes)
                continue
            ys_best = []
            for i in range(len(ys)):
                if self.cluster[i] != self.best_cluster:
                    continue
                else:
                    ys_best.append(ys[i])
            if (flav == 'g') or (flav == 'sea'):
                ax.hist(ys_best, len(ys_best), histtype = 'bar', color = colors[shape], label = label.rsplit('$', 1)[0] + '^{(%s)}$' % shape)
            else:
                ax.hist(ys_best, len(ys_best), histtype = 'bar', color = colors[shape], label = label)
            if ax.get_ylim()[1] > 100.0: ax.semilogy()
        return

class PLOTS(PDF_PLOT_CORE):

    def __init__(self, task, wdir, Q2 = None, dpi = 200):

        if  task == 1:
            self.line_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 2:
            self.groups_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 3:
            self.ratio_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains keys 'current' and 'base'
            ## 'self.ratio_figure' will plot ratio of 'current' 'xf' values over 'base' 'xf' values
            ## wdir = {'base': './results1/step04/', 'current': './results1/step06/'}
            ## input file of current step will be loaded

        if  task == 4:
            self.steps_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains integer keys
            ## the keys have to be the number of the steps
            ## 'xf' of larger key values will be plotted on top of that of smaller key values
            ## 'xf' of larger key values will also be plotted with less opacity
            ## wdir = {2: './results1/step02/', 3: './results1/step03/', 4: './results1/step04/'}
            ## input file of largest key value will be loaded

        if  task == 10:
            self.histogram_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

    def line_figure(self, wdir, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        self.cluster_average_chi2 = labels['cluster_average_chi2']
        self.best_cluster = np.argmin(self.cluster_average_chi2)

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        axs = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        # self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting PDF line figure from %s' % wdir
        self.plot_lines(axs[0], 'uv', 'r', all_cluster = True)
        self.plot_lines(axs[0], 'dv', 'b', all_cluster = True)
        self.plot_lines(axs[1], 'd/u', 'r', all_cluster = True)
        self.plot_lines(axs[2], 'db+ub', 'r', all_cluster = True)
        self.plot_lines(axs[3], 'db-ub', 'r', all_cluster = True)
        self.plot_lines(axs[4], 's+sb', 'r', all_cluster = True)
        self.plot_lines(axs[5], 'rs', 'r', all_cluster = True)
        self.plot_lines(axs[6], 'g', 'r', all_cluster = True)

        axs[0].set_ylim(0.0, 0.85)
        axs[0].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        axs[0].set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])
        axs[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = axs[0].transAxes, size = 28)
        axs[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'k', transform = axs[0].transAxes, size = 28)
        axs[1].set_ylim(0.0, 1.0)
        axs[1].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        axs[1].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])
        axs[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = axs[1].transAxes, size = 28)
        axs[2].set_ylim(0.0, 0.5)
        axs[2].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        axs[2].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])
        axs[2].text(0.07, 0.1, r'\boldmath$x(\overline{d}\!+\!\overline{u})$', color = 'k', transform = axs[2].transAxes, size = 28)
        axs[3].set_ylim(-0.05, 0.13)
        axs[3].set_yticks([-0.05, 0.0, 0.05, 0.1])
        axs[3].set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])
        axs[3].text(0.07, 0.81, r'\boldmath$x(\overline{d}\!-\!\overline{u})$', color = 'k', transform = axs[3].transAxes, size = 28)
        axs[4].set_ylim(0.0, 0.48)
        axs[4].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        axs[4].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])
        axs[4].text(0.6, 0.8, r'\boldmath$x(s\!+\!\overline{s})$', color = 'k', transform = axs[4].transAxes, size = 28)
        axs[5].set_ylim(0.0, 1.6)
        axs[5].text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = axs[5].transAxes, size = 30)
        axs[5].set_yticks([0.0, 0.5, 1.0, 1.5])
        axs[5].set_yticklabels([r'$0$', r'$0.5$', r'$1$', r'$1.5$'])
        axs[6].set_ylim(0.0, 2.5)
        axs[6].text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = axs[6].transAxes, size = 31)
        axs[6].set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])

        for ax in axs:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            # ax.legend(frameon = 0, loc = 'best')
        legends = axs[6].legend(frameon = 0, loc = 'lower left')
        for line in legends.get_lines():
            line.set_linewidth(1.0)

        for i in range(len(axs)):
            if (i % ncols) == 0: axs[i].set_ylabel(r'$xf(x)$', size = 20)
        axs[3].axhline(0.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
        for i in range(len(axs) - 1, len(axs) - 1 - ncols, -1):
            axs[i].set_xlabel(r'\boldmath$x$', size = 30)
            axs[i].set_xticklabels([r'$0.01$', r'$0.1$'])
            # axs[i].text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)

        py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-lines-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-lines-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def groups_comparison_figure(self, wdir, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        self.best_cluster = np.argmin(labels['cluster_average_chi2'])

        if Q2 == None: Q2 = conf['Q20']
        groups = load('%s/analysis/qpdlib/lhapdf-%f.dat' % (os.environ['FITPACK'], Q2))

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        axs = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting PDF band figure with following groups at Q2 = %.2f' % Q2
        print groups['groups']

        ax = axs[0]
        self.plot_band_JAM(ax, 'uv', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'uv', alpha = 0.2, label = r'$\rm %s$' % _)
        self.plot_band_JAM(ax, 'dv', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'dv', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.text(0.6, 0.1, r'\boldmath$xd_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.85)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])

        ax = axs[1]
        self.plot_band_JAM(ax, 'd/u', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'd/u', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
        ax.set_yticklabels([r'$0.0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$', r''])

        ax = axs[2]
        self.plot_band_JAM(ax, 'db+ub', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db+ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.13, r'\boldmath$x(\overline{d}\!+\!\overline{u})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.65)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        ax.set_yticklabels([r'$0.0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$'])

        ax = axs[3]
        self.plot_band_JAM(ax, 'db-ub', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db-ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.07, 0.8, r'\boldmath$x(\overline{d}\!-\!\overline{u})$', color = 'k', transform = ax.transAxes,size=28)
        ax.axhline(y = 0.0, c = 'k', ls = '-', lw = 1, alpha = 0.5)
        ax.set_ylim(-0.07, 0.11)
        ax.set_yticks([-0.05, 0.0, 0.05, 0.1])
        ax.set_yticklabels([r'-$0.05$', r'$0.0$', r'$0.05$', r'$0.1$'])

        ax = axs[4]
        self.plot_band_JAM(ax, 's+sb', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 's+sb', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.6, 0.8, r'\boldmath$x(s\!+\!\overline{s})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.48)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0.0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = axs[5]
        self.plot_band_JAM(ax, 'rs', 'r', r'$\rm JAM$')
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'rs', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 1.6)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5])
        ax.set_yticklabels([r'$0.0$', r'$0.5$', r'$1$', r'$1.5$'])

        ax = axs[6]
        self.plot_band_JAM(ax, 'g', 'r', r'$\rm JAM~Jet$')
        for _ in groups['XF']:
            if _ == 'CJ15nlo': label=r'$\rm CJ15$'
            elif _ == 'NNPDF31_nlo_as_0118': label = r'$\rm NNPDF3.1$'
            elif _ == 'MMHT2014nlo68cl': label = r'$\rm MMHT14$'
            elif _ == 'ABMP16_3_nlo': label = r'$\rm ABMP16$'
            elif _ == 'CSKK_nnlo_EIG': label = r'$\rm CSSK$'
            elif _ == 'JAM19PDF_proton_nlo': label = r'$\rm JAM19$'
            else: label = r'$\rm %s$' % _
            self.plot_band(ax, groups, _, 'g', alpha = 0.2, label = label)
        ax.text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.set_ylim(0.0, 3.9)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
        ax.set_yticklabels([r'$0.0$', r'', r'$1.0$', r'', r'$2.0$', r'', r'$3.0$', r'$3.5$'])
        ax.legend(frameon = 0, loc = 'best', bbox_to_anchor = [1.87, 0.9], fontsize = 15)

        for ax in axs:
            ax.semilogx()
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            if ax == axs[5] or ax == axs[6]:
                ax.set_xlabel(r'\boldmath$x$', size = 30)
                ax.set_xticklabels([r'$0.01$', r'$0.1$'])
                ax.text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)
            else:
                ax.set_xticklabels([r'', r''])

        for i in range(len(axs) - 1, len(axs) - 1 - ncols, -1):
            axs[i].set_xlabel(r'\boldmath$x$', size = 30)
            axs[i].set_xticklabels([r'$0.01$', r'$0.1$'])
            # axs[i].text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)

        py.tight_layout()
        py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-groups-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-groups-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def ratio_figure(self, wdirs, Q2, dpi):
        load_config('%s/input.py' % wdirs['current'])
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdirs['current'], istep))
        self.cluster = labels['cluster']
        self.best_cluster = np.argmin(labels['cluster_average_chi2'])
        base_step = int(re.findall(r'\d+', wdirs['base'].split('/')[2])[0])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        axs = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdirs['current'], istep))
            self.base_xf_data = load('%s/data/pdf-%d.dat' % (wdirs['base'], base_step))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['current'], istep, Q2))
            self.base_xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['base'], base_step, Q2))

        self.base_xf_data['mean_xf'] = {}
        for flavor in self.base_xf_data['XF']:
            self.base_xf_data['mean_xf'][flavor] = np.mean(self.base_xf_data['XF'][flavor], axis = 0)

        print '\nplotting PDF ratio figure of %s over %s at Q2 = %.2f' % (wdirs['current'], wdirs['base'], Q2)
        self.plot_ratios(axs[0], 'uv', 'r')
        self.plot_ratios(axs[0], 'dv', 'b')
        self.plot_ratios(axs[1], 'd/u', 'r')
        self.plot_ratios(axs[2], 'db+ub', 'r')
        self.plot_ratios(axs[3], 'db-ub', 'r')
        self.plot_ratios(axs[4], 's+sb', 'r')
        self.plot_ratios(axs[5], 'rs', 'r')
        self.plot_ratios(axs[6], 'g', 'r')

        axs[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'r', transform = axs[0].transAxes, size = 28)
        axs[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'b', transform = axs[0].transAxes, size = 28)
        axs[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = axs[1].transAxes, size = 28)
        axs[2].text(0.07, 0.1, r'\boldmath$x(\overline{d}\!+\!\overline{u})$', color = 'k', transform = axs[2].transAxes, size = 28)
        axs[3].text(0.07, 0.81, r'\boldmath$x(\overline{d}\!-\!\overline{u})$', color = 'k', transform = axs[3].transAxes, size = 28)
        axs[4].text(0.07, 0.8, r'\boldmath$x(s\!+\!\overline{s})$', color = 'k', transform = axs[4].transAxes, size = 28)
        axs[5].text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = axs[5].transAxes, size = 30)
        axs[6].text(0.07, 0.8, r'\boldmath$xg$', color = 'k', transform = axs[6].transAxes, size = 31)

        for ax in axs:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_ylim(0.8, 1.2)
            ax.set_yticks([0.9, 1.0, 1.1])
            ax.set_yticklabels([r'$0.9$', r'$1.0$', r'$1.1$'])
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.axhline(1.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
            ax.legend(frameon = 0, loc = 'best')

        for i in range(len(axs)):
            if (i % ncols) == 0: axs[i].set_ylabel(r'$xf(x)$', size = 20)
        for i in range(len(axs) - 1, len(axs) - 1 - ncols, -1):
            axs[i].set_xlabel(r'\boldmath$x$', size = 30)
            axs[i].set_xticklabels([r'$0.01$', r'$0.1$'])
            # axs[i].text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)


        py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-ratio-%d.png' % (wdirs['current'], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-ratio-%d-%f.png' % (wdirs['current'], istep, Q2), dpi = dpi)

    def steps_comparison_figure(self, wdirs, Q2, dpi):
        step_keys = sorted(wdirs.keys())
        load_config('%s/input.py' % wdirs[step_keys[-1]])
        istep = core.get_istep()

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        axs = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]
        colors = ['coral', 'grey', 'wheat', 'blueviolet', 'brown', 'olivedrab', 'deepskyblue']
        # colors = ['k', 'r', 'g', 'Yellow', 'b']

        ## PDF
        self.xf_data = {}
        self.cluster = {}
        self.best_cluster = {}
        labels = {}
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d.dat' % (wdirs[step_key], step_key))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = np.argmin(labels[step_key]['cluster_average_chi2'])
        else:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d-%f.dat' % (wdirs[step_key], step_key, Q2))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = np.argmin(labels[step_key]['cluster_average_chi2'])

        print '\nplotting PDF band figure with following steps at Q2 = %.2f' % Q2
        print [wdirs[_] for _ in wdirs]

        ax = axs[0]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'uv', color, alpha)
            self.plot_band_steps(ax, step, 'dv', color, alpha)
        ax.text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.text(0.6, 0.1, r'\boldmath$xd_v$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.85)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])

        ax = axs[1]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'd/u', color, alpha)
        ax.text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 0.9)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])

        ax = axs[2]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db+ub', color, alpha)
        ax.text(0.07, 0.13, r'\boldmath$x(\overline{d}\!+\!\overline{u})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.5)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = axs[3]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db-ub', color, alpha)
        ax.text(0.07, 0.8, r'\boldmath$x(\overline{d}\!-\!\overline{u})$', color = 'k', transform = ax.transAxes,size=28)
        ax.axhline(y = 0.0, c = 'k', linestyle = 'dashdot', linewidth = 0.5, alpha = 0.5)
        ax.set_ylim(-0.07, 0.13)
        ax.set_yticks([-0.05, 0.0, 0.05, 0.1])
        ax.set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])

        ax = axs[4]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 's+sb', color, alpha)
        ax.text(0.6, 0.8, r'\boldmath$x(s\!+\!\overline{s})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.set_ylim(0.0, 0.48)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        ax.set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])

        ax = axs[5]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'rs', color, alpha)
        ax.text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = ax.transAxes, size = 30)
        ax.set_ylim(0.0, 1.6)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5])
        ax.set_yticklabels([r'$0$', r'$0.5$', r'$1$', r'$1.5$'])

        ax = axs[6]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            label = r'$\mathrm{step~%02d}$' % step
            self.plot_band_steps(ax, step, 'g', color, alpha, label)
        ax.text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.set_ylim(0.0, 2.5)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.legend(frameon = 0, loc = 'best', bbox_to_anchor = [1.87, 0.9], fontsize = 20)

        for ax in axs:
            ax.semilogx()
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])

        for i in range(len(axs) - 1, len(axs) - 1 - ncols, -1):
            axs[i].set_xlabel(r'\boldmath$x$', size = 30)
            axs[i].set_xticklabels([r'$0.01$', r'$0.1$'])
            axs[i].text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)

        py.tight_layout()
        py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-steps-%d.png' % (wdirs[step_keys[-1]], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-steps-%d-%f.png' % (wdirs[step_keys[-1]], istep, Q2), dpi = dpi)

    def histogram_figure(self, wdir, Q2, dpi):
        ## plot histogram for parameters of all flavors
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        # jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        self.best_cluster = np.argmin(labels['cluster_average_chi2'])

        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.pdf_parameter_data = load('%s/data/pdf-parameters-%d.dat' % (wdir, istep))
        else:
            self.pdf_parameter_data = load('%s/data/pdf-parameters-%d-%f.dat' % (wdir, istep, Q2))

        shapes = self.pdf_parameter_data.keys()
        flavors = self.pdf_parameter_data[shapes[0]].keys()

        for _ in sorted(self.pdf_parameter_data[shapes[0]][flavors[0]]):
            print '\nplotting PDF parameter %s histograms from %s' % (_, wdir)
            if len(flavors) == 8:
                n_rows, n_columns = 4, 2
                figure = py.figure(figsize = (n_columns * 5.0, n_rows * 3.0))
                axs = [py.subplot(n_rows, n_columns, count + 1) for count in range(8)]
            elif len(flavors) == 7:
                n_rows, n_columns = 4, 2
                figure = py.figure(figsize = (n_columns * 5.0, n_rows * 3.0))
                axs = [py.subplot(n_rows, n_columns, count + 1) for count in range(8) if count != 7]
            else:
                sys.exit('please make a subplot configuration for your flavors')

            if sorted(flavors) == sorted(['g', 'sea', 'uv', 'ub', 'dv', 'db', 's', 'sb']):
                self.pdf_histogram(axs[0], 'g', _, r'\boldmath$g$')
                self.pdf_histogram(axs[1], 'sea', _, r'\boldmath$\mathrm{sea}$')
                self.pdf_histogram(axs[2], 'uv', _, r'\boldmath$u_v$')
                self.pdf_histogram(axs[3], 'ub', _, r'\boldmath$\overline{u}$')
                self.pdf_histogram(axs[4], 'dv', _, r'\boldmath$d_v$')
                self.pdf_histogram(axs[5], 'db', _, r'\boldmath$\overline{d}$')
                self.pdf_histogram(axs[6], 's', _, r'\boldmath$s$')
                self.pdf_histogram(axs[7], 'sb', _, r'\boldmath$\overline{s}$')
            elif sorted(flavors) == sorted(['g', 'uv', 'ub', 'dv', 'db', 's', 'sb']):
                self.pdf_histogram(axs[0], 'uv', _, r'\boldmath$u_v$')
                self.pdf_histogram(axs[1], 'ub', _, r'\boldmath$\overline{u}$')
                self.pdf_histogram(axs[2], 'dv', _, r'\boldmath$d_v$')
                self.pdf_histogram(axs[3], 'db', _, r'\boldmath$\overline{d}$')
                self.pdf_histogram(axs[4], 's', _, r'\boldmath$s$')
                self.pdf_histogram(axs[5], 'sb', _, r'\boldmath$\overline{s}$')
                self.pdf_histogram(axs[6], 'g', _, r'\boldmath$g$')

            for ax in axs:
                # ax.semilogx()
                # ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
                # ax.set_xticklabels([r'', r''])
                ax.legend(frameon = 0, loc = 'best', fontsize = 17)

            axs[0].title.set_text(r'$\mathrm{parameter}~%s~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV^2}$' % (_, Q2))
            axs[1].title.set_text(r'$\mathrm{parameter}~%s~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV^2}$' % (_, Q2))

            py.tight_layout()
            if Q2 == conf['Q20']:
                py.savefig('%s/gallery/pdf-histogram-%s-%d.png' % (wdir, _, istep), dpi = dpi)
            else:
                py.savefig('%s/gallery/pdf-histogram-%s-%d-%f.png' % (wdir, _, istep, Q2), dpi = dpi)

class JET(PDF_PLOT_CORE):

    def __init__(self, task, wdir, Q2 = None, dpi = 200):

        if  task == 1:
            self.line_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 2:
            self.groups_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 3:
            self.ratio_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains keys 'current' and 'base'
            ## 'self.ratio_figure' will plot ratio of 'current' 'xf' values over 'base' 'xf' values
            ## wdir = {'base': './results1/step04/', 'current': './results1/step06/'}
            ## input file of current step will be loaded

        if  task == 4:
            self.steps_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains integer keys
            ## the keys have to be the number of the steps
            ## 'xf' of larger key values will be plotted on top of that of smaller key values
            ## 'xf' of larger key values will also be plotted with less opacity
            ## wdir = {2: './results1/step02/', 3: './results1/step03/', 4: './results1/step04/'}
            ## input file of largest key value will be loaded

        if  task == 10:
            self.histogram_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

    def line_figure(self, wdir, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        self.cluster_average_chi2 = labels['cluster_average_chi2']
        self.best_cluster = np.argmin(self.cluster_average_chi2)

        nrows, ncols = 3, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        axs = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(6) if cnt != 5]

        ## PDF
        # self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting PDF line figure from %s' % wdir
        self.plot_lines(axs[0], 'uv', 'r', all_cluster = True)
        self.plot_lines(axs[0], 'dv', 'b', all_cluster = True)
        self.plot_lines(axs[1], 'd/u', 'r', all_cluster = True)
        self.plot_lines(axs[2], 'db+ub', 'r', all_cluster = True)
        self.plot_lines(axs[3], 'db-ub', 'r', all_cluster = True)
        self.plot_lines(axs[4], 'g', 'r', all_cluster = True)

        axs[0].set_ylim(0.0, None)
        # axs[0].set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        # axs[0].set_yticklabels([r'$0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])
        axs[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = axs[0].transAxes, size = 17)
        axs[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'k', transform = axs[0].transAxes, size = 17)
        axs[0].title.set_text(r'\boldmath$\mathrm{valance~quarks}~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV}^2$' % Q2)
        axs[1].set_ylim(0.0, 1.0)
        # axs[1].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        # axs[1].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])
        # axs[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = axs[1].transAxes, size = 28)
        axs[1].title.set_text(r'\boldmath$d/u~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV}^2$' % Q2)
        axs[2].set_ylim(-0.02, 0.5)
        # axs[2].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
        # axs[2].set_yticklabels([r'$0$', r'', r'$0.2$', r'', r'$0.4$'])
        # axs[2].text(0.07, 0.1, r'\boldmath$x(\overline{d}\!+\!\overline{u})$', color = 'k', transform = axs[2].transAxes, size = 28)
        axs[2].title.set_text(r'\boldmath$x\left(\overline{d}~+~\overline{u}\right)$')
        axs[3].set_ylim(-0.05, None)
        axs[3].set_yticks([-0.05, 0.0, 0.05, 0.1])
        axs[3].set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])
        # axs[3].text(0.07, 0.81, r'\boldmath$x(\overline{d}\!-\!\overline{u})$', color = 'k', transform = axs[3].transAxes, size = 28)
        axs[3].title.set_text(r'\boldmath$x\left(\overline{d}~-~\overline{u}\right)$')
        axs[4].set_ylim(0.0, 2.5)
        # axs[4].set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        # axs[4].text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = axs[4].transAxes, size = 31)
        axs[4].title.set_text(r'\boldmath$xg$')

        for ax in axs:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
        legends = axs[4].legend(frameon = 0, loc = 'lower left')
        for line in legends.get_lines():
            line.set_linewidth(1.0)

        axs[3].axhline(0.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
        axs[2].axhline(0.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
        for i in range(len(axs)):
            if (i % ncols) == 0: axs[i].set_ylabel(r'$xf(x)$', size = 20)
        for i in range(len(axs) - 1, len(axs) - 1 - ncols, -1):
            axs[i].set_xlabel(r'\boldmath$x$', size = 30)
            axs[i].set_xticklabels([r'$0.01$', r'$0.1$'])
            axs[i].text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)

        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-lines-jet-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-lines-jet-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def groups_comparison_figure(self, wdir, Q2, dpi):
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        self.best_cluster = np.argmin(labels['cluster_average_chi2'])

        if Q2 == None: Q2 = conf['Q20']
        groups = load('%s/analysis/qpdlib/lhapdf-%f.dat' % (os.environ['FITPACK'], Q2))

        nrows, ncols = 3, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        axs = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(6) if cnt != 5]

        ## PDF
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting PDF band figure with following groups at Q2 = %.2f' % Q2
        print groups['groups']

        jet_distinction = {'alpha': 0.4, 'zorder': 1}

        ax = axs[0]
        self.plot_band_JAM(ax, 'uv', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'uv', alpha = 0.2, label = r'$\rm %s$' % _)
        self.plot_band_JAM(ax, 'dv', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'dv', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.set_ylim(0.0, 0.8)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        ax.set_yticklabels([r'$0.0$', r'$0.2$', r'$0.4$', r'$0.6$', r'$0.8$'])
        ax.text(0.4, 0.81, r'\boldmath$xu_v$', color = 'k', transform = ax.transAxes, size = 17)
        ax.text(0.6, 0.1, r'\boldmath$xd_v$', color = 'k', transform = ax.transAxes, size = 17)
        ax.title.set_text(r'\boldmath$\mathrm{valance~quarks}~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV}^2$' % Q2)

        ax = axs[1]
        self.plot_band_JAM(ax, 'd/u', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'd/u', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        ax.set_yticklabels([r'$0.0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$', r'', r'$0.8$'])
        # ax.text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = ax.transAxes, size = 30)
        ax.title.set_text(r'\boldmath$d/u~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV}^2$' % Q2)

        ax = axs[2]
        self.plot_band_JAM(ax, 'db+ub', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db+ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.set_ylim(0.0, 0.65)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        ax.set_yticklabels([r'$0.0$', r'', r'$0.2$', r'', r'$0.4$', r'', r'$0.6$'])
        # ax.text(0.07, 0.13, r'\boldmath$x(\overline{d}\!+\!\overline{u})$', color = 'k', transform = ax.transAxes, size = 28)
        ax.title.set_text(r'\boldmath$x\left(\overline{d}~+~\overline{u}\right)$')

        ax = axs[3]
        self.plot_band_JAM(ax, 'db-ub', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            self.plot_band(ax, groups, _, 'db-ub', alpha = 0.2, label = r'$\rm %s$' % _)
        ax.axhline(y = 0.0, c = 'k', ls = '-', lw = 1, alpha = 0.5)
        ax.set_ylim(-0.07, 0.1)
        ax.set_yticks([-0.05, 0.0, 0.05, 0.1])
        ax.set_yticklabels([r'-$0.05$', r'$0$', r'$0.05$', r'$0.1$'])
        # ax.text(0.07, 0.8, r'\boldmath$x(\overline{d}\!-\!\overline{u})$', color = 'k', transform = ax.transAxes,size=28)
        ax.title.set_text(r'\boldmath$x\left(\overline{d}~-~\overline{u}\right)$')

        ax = axs[4]
        self.plot_band_JAM(ax, 'g', 'r', r'$\rm JAM~Jet$', distinction = jet_distinction)
        for _ in groups['XF']:
            if _ == 'CJ15nlo': label=r'$\rm CJ15$'
            elif _ == 'NNPDF31_nlo_as_0118': label = r'$\rm NNPDF3.1$'
            elif _ == 'MMHT2014nlo68cl': label = r'$\rm MMHT14$'
            elif _ == 'ABMP16_3_nlo': label = r'$\rm ABMP16$'
            elif _ == 'CSKK_nnlo_EIG': label = r'$\rm CSSK$'
            elif _ == 'JAM19PDF_proton_nlo': label = r'$\rm JAM19$'
            else: label = r'$\rm %s$' % _
            self.plot_band(ax, groups, _, 'g', alpha = 0.2, label = label)
        ax.set_ylim(0.0, 3.9)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
        ax.set_yticklabels([r'$0.0$', r'', r'$1.0$', r'', r'$2.0$', r'', r'$3.0$', r'$3.5$'])
        # ax.text(0.8, 0.8, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.title.set_text(r'\boldmath$xg$')
        # ax.legend(frameon = 0, loc = 'best', bbox_to_anchor = [1.87, 0.9], fontsize = 15)
        ax.legend(frameon = 0, loc = 'best', fontsize = 10)

        for ax in axs:
            ax.semilogx()
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1, 0.5])

        for i in range(len(axs) - 1, len(axs) - 1 - ncols, -1):
            axs[i].set_xlabel(r'\boldmath$x$', size = 30)
            axs[i].set_xticklabels([r'$0.01$', r'$0.1$'])
            # axs[i].text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)

        # fig.suptitle(r'$Q^2 ~=~ %.2f$' % Q2) ## this does not work with 'tight_layout'
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-groups-jet-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-groups-jet-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def ratio_figure(self, wdirs, Q2, dpi):
        load_config('%s/input.py' % wdirs['current'])
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdirs['current'], istep))
        self.cluster = labels['cluster']
        self.best_cluster = np.argmin(labels['cluster_average_chi2'])
        base_step = int(re.findall(r'\d+', wdirs['base'].split('/')[2])[0])

        nrows, ncols = 4, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        axs = [py.subplot(nrows, ncols, cnt+1) for cnt in range(8) if cnt != 7]

        ## PDF
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/pdf-%d.dat' % (wdirs['current'], istep))
            self.base_xf_data = load('%s/data/pdf-%d.dat' % (wdirs['base'], base_step))
        else:
            self.xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['current'], istep, Q2))
            self.base_xf_data = load('%s/data/pdf-%d-%f.dat' % (wdirs['base'], base_step, Q2))

        self.base_xf_data['mean_xf'] = {}
        for flavor in self.base_xf_data['XF']:
            self.base_xf_data['mean_xf'][flavor] = np.mean(self.base_xf_data['XF'][flavor], axis = 0)

        print '\nplotting PDF ratio figure of %s over %s at Q2 = %.2f' % (wdirs['current'], wdirs['base'], Q2)
        self.plot_ratios(axs[0], 'uv', 'r')
        self.plot_ratios(axs[0], 'dv', 'b')
        self.plot_ratios(axs[1], 'd/u', 'r')
        self.plot_ratios(axs[2], 'db+ub', 'r')
        self.plot_ratios(axs[3], 'db-ub', 'r')
        self.plot_ratios(axs[4], 's+sb', 'r')
        self.plot_ratios(axs[5], 'rs', 'r')
        self.plot_ratios(axs[6], 'g', 'r')

        axs[0].title.set_text(r'\boldmath$\mathrm{valance~quarks}~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV}^2$' % Q2)
        axs[1].title.set_text(r'\boldmath$d/u$')
        axs[2].title.set_text(r'\boldmath$x\left(\overline{d}~+~\overline{u}\right)$')
        axs[3].title.set_text(r'\boldmath$x\left(\overline{d}~-~\overline{u}\right)$')
        axs[4].title.set_text(r'\boldmath$x\left(s~+~\overline{s}\right)$')
        axs[5].title.set_text(r'\boldmath$R_s$')
        axs[6].title.set_text(r'\boldmath$xg$')

        axs[0].text(0.4, 0.81, r'\boldmath$xu_v$', color = 'r', transform = axs[0].transAxes, size = 17)
        axs[0].text(0.4, 0.07, r'\boldmath$xd_v$', color = 'b', transform = axs[0].transAxes, size = 17)
        # axs[1].text(0.07, 0.13, r'\boldmath$d/u$', color = 'k', transform = axs[1].transAxes, size = 28)
        # axs[2].text(0.07, 0.1, r'\boldmath$x(\overline{d}\!+\!\overline{u})$', color = 'k', transform = axs[2].transAxes, size = 28)
        # axs[3].text(0.07, 0.81, r'\boldmath$x(\overline{d}\!-\!\overline{u})$', color = 'k', transform = axs[3].transAxes, size = 28)
        # axs[4].text(0.07, 0.8, r'\boldmath$x(s\!+\!\overline{s})$', color = 'k', transform = axs[4].transAxes, size = 28)
        # axs[5].text(0.05, 0.78, r'\boldmath$R_s$', color = 'k', transform = axs[5].transAxes, size = 30)
        # axs[6].text(0.07, 0.8, r'\boldmath$xg$', color = 'k', transform = axs[6].transAxes, size = 31)

        for ax in axs:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_ylim(0.8, 1.2)
            ax.set_yticks([0.9, 1.0, 1.1])
            ax.set_yticklabels([r'$0.9$', r'$1.0$', r'$1.1$'])
            ax.set_xticks([1e-2, 1e-1])
            ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.axhline(1.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
            ax.legend(frameon = 0, loc = 'best')

        for i in range(len(axs)):
            if (i % ncols) == 0: axs[i].set_ylabel(r'$xf(x)$', size = 20)
        for i in range(len(axs) - 1, len(axs) - 1 - ncols, -1):
            axs[i].set_xlabel(r'\boldmath$x$', size = 30)
            axs[i].set_xticklabels([r'$0.01$', r'$0.1$'])
            # axs[i].text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)

        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-ratio-jet-%d.png' % (wdirs['current'], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-ratio-jet-%d-%f.png' % (wdirs['current'], istep, Q2), dpi = dpi)

    def steps_comparison_figure(self, wdirs, Q2, dpi):
        step_keys = sorted(wdirs.keys())
        load_config('%s/input.py' % wdirs[step_keys[-1]])
        istep = core.get_istep()

        nrows, ncols = 1, 1
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        # axs = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(6) if cnt != 5]
        ax = py.subplot(nrows, ncols, 1)
        colors = ['coral', 'grey', 'wheat', 'blueviolet', 'brown', 'olivedrab', 'deepskyblue']
        # colors = ['k', 'r', 'g', 'Yellow', 'b']

        ## PDF
        self.xf_data = {}
        self.cluster = {}
        self.best_cluster = {}
        labels = {}
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d.dat' % (wdirs[step_key], step_key))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = np.argmin(labels[step_key]['cluster_average_chi2'])
        else:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/pdf-%d-%f.dat' % (wdirs[step_key], step_key, Q2))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = np.argmin(labels[step_key]['cluster_average_chi2'])

        print '\nplotting PDF band figure with following steps at Q2 = %.2f' % Q2
        print [wdirs[_] for _ in wdirs]

        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            label = r'$\mathrm{step~%02d}$' % step
            self.plot_band_steps(ax, step, 'g', color, alpha, label)
        # ax.text(0.3, 0.3, r'\boldmath$xg$', color = 'k', transform = ax.transAxes, size = 31)
        ax.set_ylim(0.0, None)
        # ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.legend(frameon = 0, loc = 'upper right', fontsize = 10)
        ax.title.set_text(r'\boldmath$xg~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV}^2$' % Q2)

        ax.semilogx()
        ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
        ax.xaxis.set_label_coords(0.7, -0.05)
        ax.set_xlim(8e-3, 9e-1)
        ax.set_xticks([1e-2, 1e-1])
        ax.set_xlabel(r'\boldmath$x$', size = 30)
        ax.set_xticklabels([r'$0.01$', r'$0.1$'])

        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/pdf-steps-jet-%d.png' % (wdirs[step_keys[-1]], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/pdf-steps-jet-%d-%f.png' % (wdirs[step_keys[-1]], istep, Q2), dpi = dpi)

    def histogram_figure(self, wdir, Q2, dpi):
        ## plot histogram for parameters of all flavors
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        self.best_cluster = np.argmin(labels['cluster_average_chi2'])

        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.pdf_parameter_data = load('%s/data/pdf-parameters-%d.dat' % (wdir, istep))
        else:
            self.pdf_parameter_data = load('%s/data/pdf-parameters-%d-%f.dat' % (wdir, istep, Q2))

        shapes = self.pdf_parameter_data.keys()
        flavors = self.pdf_parameter_data[shapes[0]].keys()

        for _ in sorted(self.pdf_parameter_data[shapes[0]][flavors[0]]):
            print '\nplotting PDF parameter %s histograms from %s' % (_, wdir)
            if len(flavors) == 8:
                n_rows, n_columns = 4, 2
                figure = py.figure(figsize = (n_columns * 5.0, n_rows * 3.0))
                axs = [py.subplot(n_rows, n_columns, count + 1) for count in range(8)]
            elif len(flavors) == 7:
                n_rows, n_columns = 4, 2
                figure = py.figure(figsize = (n_columns * 5.0, n_rows * 3.0))
                axs = [py.subplot(n_rows, n_columns, count + 1) for count in range(8) if count != 7]
            else:
                sys.exit('please make a subplot configuration for your flavors')

            if sorted(flavors) == sorted(['g', 'sea', 'uv', 'ub', 'dv', 'db', 's', 'sb']):
                self.pdf_histogram(axs[0], 'g', _, r'\boldmath$g$')
                self.pdf_histogram(axs[1], 'sea', _, r'\boldmath$\mathrm{sea}$')
                self.pdf_histogram(axs[2], 'uv', _, r'\boldmath$u_v$')
                self.pdf_histogram(axs[3], 'ub', _, r'\boldmath$\overline{u}$')
                self.pdf_histogram(axs[4], 'dv', _, r'\boldmath$d_v$')
                self.pdf_histogram(axs[5], 'db', _, r'\boldmath$\overline{d}$')
                self.pdf_histogram(axs[6], 's', _, r'\boldmath$s$')
                self.pdf_histogram(axs[7], 'sb', _, r'\boldmath$\overline{s}$')
            elif sorted(flavors) == sorted(['g', 'uv', 'ub', 'dv', 'db', 's', 'sb']):
                self.pdf_histogram(axs[0], 'uv', _, r'\boldmath$u_v$')
                self.pdf_histogram(axs[1], 'ub', _, r'\boldmath$\overline{u}$')
                self.pdf_histogram(axs[2], 'dv', _, r'\boldmath$d_v$')
                self.pdf_histogram(axs[3], 'db', _, r'\boldmath$\overline{d}$')
                self.pdf_histogram(axs[4], 's', _, r'\boldmath$s$')
                self.pdf_histogram(axs[5], 'sb', _, r'\boldmath$\overline{s}$')
                self.pdf_histogram(axs[6], 'g', _, r'\boldmath$g$')

            for ax in axs:
                # ax.semilogx()
                # ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
                # ax.set_xticklabels([r'', r''])
                ax.legend(frameon = 0, loc = 'best', fontsize = 17)

            axs[0].title.set_text(r'$\mathrm{parameter}~%s~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV^2}$' % (_, Q2))
            axs[1].title.set_text(r'$\mathrm{parameter}~%s~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV^2}$' % (_, Q2))

            py.tight_layout()
            if Q2 == conf['Q20']:
                py.savefig('%s/gallery/pdf-histogram-jet-%s-%d.png' % (wdir, _, istep), dpi = dpi)
            else:
                py.savefig('%s/gallery/pdf-histogram-jet-%s-%d-%f.png' % (wdir, _, istep, Q2), dpi = dpi)

if __name__ == '__main__':
    pass
