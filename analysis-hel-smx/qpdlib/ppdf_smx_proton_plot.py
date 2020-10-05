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

class PPDF_PLOT_CORE:

    def plot_lines(self, ax, flav, c, lab = '', all_cluster = 1): ## plot all the PDF replicas as lines
        n_replicas = len(self.xf_data['XF'][flav])
        line_width = 10.0 / float(n_replicas)
        if line_width > 1.0: line_width = 1.0
        xs = self.xf_data['X']

        if all_cluster:
            clusters = list(set(self.cluster))
            colors = ['red', 'gold', 'dodgerblue', 'black']
            if all([round(_, 2) == round(self.cluster_average_chi2[0], 2) for _ in self.cluster_average_chi2]):
                colors = ['red' for _ in colors]
            count_clusters = [0 for _ in clusters]
            for i in range(n_replicas):
                if len(list(set(colors))) == 1:
                    ax.plot(xs, self.xf_data['XF'][flav][i], color = colors[self.cluster[i]], label = lab, linewidth = line_width)
                else:
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

    def ppdf_histogram(self, ax, flav, parameter, label): ## histogram for PDF normalization
        colors = ['orangered', 'limegreen', 'dodgerblue', 'fuchsia']
        shapes = self.ppdf_parameter_data.keys()
        for shape in shapes:
            if flav in self.ppdf_parameter_data[shape]:
                ys = self.ppdf_parameter_data[shape][flav][parameter]
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

class PLOTS(PPDF_PLOT_CORE):

    def __init__(self, task, wdir, Q2 = None, dpi = 200):

        if  task == 1:
            self.line_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 4:
            self.steps_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains integer keys
            ## the keys have to be the number of the steps
            ## 'xf' of larger key values will be plotted on top of that of smaller key values
            ## 'xf' of larger key values will also be plotted with less opacity
            ## wdir = {6: './results1/step06/', 7: './results1/step07/', 9: './results1/step09/'}
            ## input file of largest key value will be loaded

        if  task == 10:
            self.histogram_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

    def line_figure(self, wdir, Q2, dpi):
        ## this block may contain something not useful
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        # jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        self.cluster_average_chi2 = labels['cluster_average_chi2']
        self.best_cluster = np.argmin(self.cluster_average_chi2)

        nrows, ncols = 2, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        axs = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(4)]

        ## polarized PDF
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/ppdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/ppdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting polarized PDF line figure from %s at Q2 = %.2f' % (wdir, Q2)
        self.plot_lines(axs[0], 'up', 'r', all_cluster = 2)
        self.plot_lines(axs[1], 'dp', 'r', all_cluster = 2)
        self.plot_lines(axs[2], 'sp', 'r', all_cluster = 2)
        self.plot_lines(axs[3], 'g', 'r', all_cluster = 2)

        axs[0].set_ylim(0.0, 0.5)
        axs[0].set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        axs[0].set_yticklabels([r'$0.1$', r'$0.2$', r'$0.3$', r'$0.4$', r'$0.5$'])
        # axs[0].text(0.6, 0.3, r'\boldmath$x\Delta u^+$', color = 'k', transform = axs[0].transAxes, size = 28)
        axs[0].title.set_text(r'\boldmath$x\Delta u^+~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV}^2$' % Q2)
        axs[1].set_ylim(-0.2, 0.0)
        axs[1].set_yticks([-0.05, -0.1, -0.15])
        axs[1].set_yticklabels([r'$-0.05$', r'$-0.10$', r'$-0.15$'])
        # axs[1].text(0.05, 0.15, r'\boldmath$x\Delta d^+$', color = 'k', transform = axs[1].transAxes, size = 28)
        axs[1].title.set_text(r'\boldmath$x\Delta d^+~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV}^2$' % Q2)
        axs[2].set_ylim(-0.1, 0.05)
        axs[2].set_yticks([-0.08, -0.04, 0.0, 0.04])
        axs[2].set_yticklabels([r'$-0.08$', r'$-0.04$', r'$0.00$', r'$0.04$'])
        # axs[2].text(0.1, 0.7, r'\boldmath$x\Delta s^+$', color = 'k', transform = axs[2].transAxes, size = 28)
        axs[2].title.set_text(r'\boldmath$x\Delta s^+$')
        axs[3].set_ylim(-0.11, 0.3)
        axs[3].set_yticks([-0.1, 0.0, 0.1, 0.2])
        axs[3].set_yticklabels([r'$-0.1$', r'$0.0$', r'$0.1$', r'$0.2$'])
        # axs[3].text(0.7, 0.8, r'\boldmath$x\Delta g$', color = 'k', transform = axs[3].transAxes, size = 28)
        axs[3].title.set_text(r'\boldmath$x\Delta g$')

        for ax in axs:
            ax.semilogx()
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            # ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])

        axs[2].axhline(0.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
        axs[3].axhline(0.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
        # for i in range(len(axs)):
        #     if (i % ncols) == 0: axs[i].set_ylabel(r'$xf(x)$', size = 20)
        for i in range(len(axs) - 1, len(axs) - 1 - ncols, -1):
            axs[i].set_xlabel(r'\boldmath$x$', size = 30)
            # axs[i].text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)
            axs[i].set_xticks([0.01, 0.1, 0.5, 0.8])
            axs[i].set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$', r'$0.8$'])

        legends = axs[0].legend(frameon = 0, loc = 'best')
        for line in legends.get_lines():
            line.set_linewidth(1.0)

        py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/ppdf-lines-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-lines-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def steps_comparison_figure(self, wdirs, Q2, dpi):
        step_keys = sorted(wdirs.keys())
        load_config('%s/input.py' % wdirs[step_keys[-1]])
        istep = core.get_istep()

        nrows, ncols = 2, 2
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        axs = [py.subplot(nrows, ncols, cnt + 1) for cnt in range(4)]
        colors = ['k', 'r', 'g', 'Yellow', 'b']

        ## PPDF
        self.xf_data = {}
        self.cluster = {}
        self.best_cluster = {}
        labels = {}
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/ppdf-%d.dat' % (wdirs[step_key], step_key))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = np.argmin(labels[step_key]['cluster_average_chi2'])
        else:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/ppdf-%d-%f.dat' % (wdirs[step_key], step_key, Q2))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = np.argmin(labels[step_key]['cluster_average_chi2'])

        print '\nplotting PPDF band figure with following steps at Q2 = %.2f' % Q2
        print [wdirs[_] for _ in wdirs]

        ax = axs[0]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'uv', color, alpha)
            self.plot_band_steps(ax, step, 'dv', color, alpha)
        ax.set_ylim(0.0, 0.5)
        ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        ax.set_yticklabels([r'$0.1$', r'$0.2$', r'$0.3$', r'$0.4$', r'$0.5$'])
        ax.text(0.6, 0.3, r'\boldmath$x\Delta u^+$', color = 'k', transform = axs[0].transAxes, size = 28)

        ax = axs[1]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'd/u', color, alpha)
        ax.set_ylim(-0.2, 0.0)
        ax.set_yticks([-0.05, -0.1, -0.15])
        ax.set_yticklabels([r'$-0.05$', r'$-0.10$', r'$-0.15$'])
        ax.text(0.05, 0.15, r'\boldmath$x\Delta d^+$', color = 'k', transform = axs[1].transAxes, size = 28)

        ax = axs[2]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db+ub', color, alpha)
        ax.set_ylim(-0.1, 0.05)
        ax.set_yticks([-0.08, -0.04, 0.0, 0.04])
        ax.set_yticklabels([r'$-0.08$', r'$-0.04$', r'$0.00$', r'$0.04$'])
        ax.text(0.1, 0.7, r'\boldmath$x\Delta s^+$', color = 'k', transform = axs[2].transAxes, size = 28)

        ax = axs[3]
        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            self.plot_band_steps(ax, step, 'db-ub', color, alpha)
        ax.set_ylim(-0.11, 0.3)
        ax.set_yticks([-0.1, 0.0, 0.1, 0.2])
        ax.set_yticklabels([r'$-0.1$', r'$0.0$', r'$0.1$', r'$0.2$'])
        ax.text(0.7, 0.8, r'\boldmath$x\Delta g$', color = 'k', transform = axs[3].transAxes, size = 28)

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
            ax.set_xlim(8e-3, 9e-1)
            ax.set_xticks([1e-2, 1e-1])
            # ax.xaxis.set_label_coords(0.95, -0.05)
            ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
            ax.set_xticklabels([r'', r''])
            ax.xaxis.set_label_coords(0.7, -0.05)
            ax.legend(frameon = 0, loc = 'best')

        axs[2].axhline(0.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
        axs[3].axhline(0.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
        for i in range(len(axs) - 1, len(axs) - 1 - ncols, -1):
            axs[i].set_xlabel(r'\boldmath$x$', size = 30)
            # axs[i].text(0.84, -0.095, r'$0.5$', color = 'k', transform = ax.transAxes, size = 20)
            axs[i].set_xticks([0.01, 0.1, 0.5, 0.8])
            axs[i].set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$', r'$0.8$'])

        py.tight_layout()
        py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/ppdf-steps-%d.png' % (wdirs[step_keys[-1]], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-steps-%d-%f.png' % (wdirs[step_keys[-1]], istep, Q2), dpi = dpi)

class PJET(PPDF_PLOT_CORE):

    def __init__(self, task, wdir, Q2 = None, dpi = 200):

        if  task == 1:
            self.line_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

        if  task == 4:
            self.steps_comparison_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)
            ## 'wdir' for 'self.ratio_figure' is a dictionary that contains integer keys
            ## the keys have to be the number of the steps
            ## 'xf' of larger key values will be plotted on top of that of smaller key values
            ## 'xf' of larger key values will also be plotted with less opacity
            ## wdir = {6: './results1/step06/', 7: './results1/step07/', 9: './results1/step09/'}
            ## input file of largest key value will be loaded

        if  task == 10:
            self.histogram_figure(wdir, Q2, dpi)
            # self.loop_over_steps(wdir,'simple',last)

    def line_figure(self, wdir, Q2, dpi):
        ## this block may contain something not useful
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        # jar = load('%s/data/jar-%d.dat' % (wdir, istep))
        labels = load('%s/data/labels-%d.dat' % (wdir, istep))
        self.cluster = labels['cluster']
        self.cluster_average_chi2 = labels['cluster_average_chi2']
        self.best_cluster = np.argmin(self.cluster_average_chi2)

        nrows, ncols = 1, 1
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        ax = py.subplot(nrows, ncols, 1)

        ## polarized PDF
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            self.xf_data = load('%s/data/ppdf-%d.dat' % (wdir, istep))
        else:
            self.xf_data = load('%s/data/ppdf-%d-%f.dat' % (wdir, istep, Q2))

        print '\nplotting polarized PDF line figure from %s at Q2 = %.2f' % (wdir, Q2)
        self.plot_lines(ax, 'g', 'r', all_cluster = 2)

        ax.set_ylim(-0.11, 0.3)
        ax.set_yticks([-0.1, 0.0, 0.1, 0.2])
        ax.set_yticklabels([r'$-0.1$', r'$0.0$', r'$0.1$', r'$0.2$'])
        # ax.text(0.7, 0.8, r'\boldmath$x\Delta g$', color = 'k', transform = axs[3].transAxes, size = 28)
        ax.title.set_text(r'\boldmath$x\Delta g~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV}^2$' % Q2)

        legends = ax.legend(frameon = 0, loc = 'best')
        for line in legends.get_lines():
            line.set_linewidth(1.0)

        ax.semilogx()
        ax.set_xlim(8e-3, 9e-1)
        ax.set_xticks([1e-2, 1e-1])
        ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
        ax.set_xticklabels([r'', r''])

        ax.axhline(0.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
        ax.set_xticks([0.01, 0.1, 0.5, 0.8])
        ax.set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$', r'$0.8$'])

        # py.subplots_adjust(left = 0.12, bottom = 0.08, right = 0.99, top = 0.97, wspace = None, hspace = 0.2)
        py.tight_layout()
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/ppdf-lines-pjet-%d.png' % (wdir, istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-lines-pjet-%d-%f.png' % (wdir, istep, Q2), dpi = dpi)

    def steps_comparison_figure(self, wdirs, Q2, dpi):
        step_keys = sorted(wdirs.keys())
        load_config('%s/input.py' % wdirs[step_keys[-1]])
        istep = core.get_istep()

        nrows, ncols = 1, 1
        fig = py.figure(figsize = (ncols * 5.0, nrows * 3.0))
        ax = py.subplot(nrows, ncols, 1)
        colors = ['y', 'r', 'b', 'g', 'm']

        ## PPDF
        self.xf_data = {}
        self.cluster = {}
        self.best_cluster = {}
        labels = {}
        if Q2 == None: Q2 = conf['Q20']
        if Q2 == conf['Q20']:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/ppdf-%d.dat' % (wdirs[step_key], step_key))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = np.argmin(labels[step_key]['cluster_average_chi2'])
        else:
            for step_key in step_keys:
                self.xf_data[step_key] = load('%s/data/ppdf-%d-%f.dat' % (wdirs[step_key], step_key, Q2))
                labels[step_key] = load('%s/data/labels-%d.dat' % (wdirs[step_key], step_key))
                self.cluster[step_key] = labels[step_key]['cluster']
                self.best_cluster[step_key] = np.argmin(labels[step_key]['cluster_average_chi2'])

        print '\nplotting PPDF band figure with following steps at Q2 = %.2f' % Q2
        print [wdirs[_] for _ in wdirs]

        for i in range(len(step_keys)):
            step = step_keys[i]
            # alpha = (1.0 / len(step_keys)) * (len(step_keys) - i)
            alpha = 0.5
            color = colors[i]
            label = r'$\mathrm{step~%02d}$' % step
            self.plot_band_steps(ax, step, 'g', color, alpha, label)
        ax.set_ylim(-0.11, 0.3)
        ax.set_yticks([-0.1, 0.0, 0.1, 0.2])
        ax.set_yticklabels([r'$-0.1$', r'$0.0$', r'$0.1$', r'$0.2$'])
        # ax.text(0.1, 0.1, r'\boldmath$x\Delta g$', color = 'k', transform = ax.transAxes, size = 28)
        ax.legend(frameon = 0, loc = 'lower left', fontsize = 7)

        ax.semilogx()
        ax.set_xlim(8e-3, 9e-1)
        ax.set_xticks([1e-2, 1e-1])
        ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
        ax.set_xticklabels([r'', r''])

        ax.axhline(0.0, color = 'k', linestyle = 'dashdot', linewidth = 0.5)
        ax.set_xticks([0.01, 0.1, 0.5])
        ax.set_xticklabels([r'$0.01$', r'$0.1$', r'$0.5$'])
        ax.set_ylabel(r'$x \Delta g$', size = 15)
        ax.set_xlabel(r'\boldmath$x$', size = 15)
        ax.xaxis.set_label_coords(1.05, 0.0)

        py.tight_layout()
        # py.subplots_adjust(left=0.05, bottom=0.04, right=0.99, top=0.99, wspace=0.17, hspace=0.04)
        if Q2 == conf['Q20']:
            py.savefig('%s/gallery/ppdf-steps-pjet-%d.png' % (wdirs[step_keys[-1]], istep), dpi = dpi)
        else:
            py.savefig('%s/gallery/ppdf-steps-pjet-%d-%f.png' % (wdirs[step_keys[-1]], istep, Q2), dpi = dpi)

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
            self.ppdf_parameter_data = load('%s/data/ppdf-parameters-%d.dat' % (wdir, istep))
        else:
            self.ppdf_parameter_data = load('%s/data/ppdf-parameters-%d-%f.dat' % (wdir, istep, Q2))

        shapes = self.ppdf_parameter_data.keys()
        flavors = self.ppdf_parameter_data[shapes[0]].keys()

        for _ in sorted(self.ppdf_parameter_data[shapes[0]][flavors[0]]):
            print '\nplotting PPDF parameter %s histograms from %s' % (_, wdir)
            if len(flavors) == 8:
                n_rows, n_columns = 4, 2
                figure = py.figure(figsize = (n_columns * 5.0, n_rows * 3.0))
                axs = [py.subplot(n_rows, n_columns, count + 1) for count in range(8)]
                sys.exit('please make a subplot configuration for your flavors')
            elif len(flavors) == 7:
                n_rows, n_columns = 4, 2
                figure = py.figure(figsize = (n_columns * 5.0, n_rows * 3.0))
                axs = [py.subplot(n_rows, n_columns, count + 1) for count in range(8) if count != 7]
            else:
                sys.exit('please make a subplot configuration for your flavors')

            if sorted(flavors) == sorted(['g', 'up', 'ub', 'dp', 'db', 'sp', 'sb']):
                self.ppdf_histogram(axs[0], 'up', _, r'\boldmath$u_+$')
                self.ppdf_histogram(axs[1], 'ub', _, r'\boldmath$\overline{u}$')
                self.ppdf_histogram(axs[2], 'dp', _, r'\boldmath$d_+$')
                self.ppdf_histogram(axs[3], 'db', _, r'\boldmath$\overline{d}$')
                self.ppdf_histogram(axs[4], 'sp', _, r'\boldmath$s_+$')
                self.ppdf_histogram(axs[5], 'sb', _, r'\boldmath$\overline{s}$')
                self.ppdf_histogram(axs[6], 'g', _, r'\boldmath$g$')

            for ax in axs:
                # ax.semilogx()
                # ax.tick_params(axis = 'both', which = 'both', right = True, top = True, direction = 'in', labelsize = 20)
                # ax.set_xticklabels([r'', r''])
                ax.legend(frameon = 0, loc = 'best', fontsize = 17)

            axs[0].title.set_text(r'$\mathrm{parameter}~%s~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV^2}$' % (_, Q2))
            axs[1].title.set_text(r'$\mathrm{parameter}~%s~\mathrm{at}~Q^2 = %.2f~\mathrm{GeV^2}$' % (_, Q2))

            py.tight_layout()
            if Q2 == conf['Q20']:
                py.savefig('%s/gallery/ppdf-histogram-pjet-%s-%d.png' % (wdir, _, istep), dpi = dpi)
            else:
                py.savefig('%s/gallery/ppdf-histogram-pjet-%s-%d-%f.png' % (wdir, _, istep, Q2), dpi = dpi)

if __name__ == '__main__':
    pass
