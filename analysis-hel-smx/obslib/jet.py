#!/usr/bin/env python
import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
import scipy as sp
# import time
# from glob import glob

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

## from fitpack obslib
from obslib.jets.jet_tools import find_eta_bins

def plot_data(wdir, data, istep, dpi, plot_with_factor):
    if len(data) == 4:
        n_rows, n_columns = 2, 2
    else:
        print('number of dataset is different from 4')
        print('you need to specify number of rows and columns for subplots in "plot_data"')
        return
    figure, ax = py.subplots(n_rows, n_columns)

    plot_styles = ['-']
    plot_colors = ['b', 'g', 'm', 'y', 'c', 'r']
    # plot_styles = ['-', '-.', '--', ':', '-', '--']
    plot_with_factor_input = plot_with_factor

    count = 0
    for idx in data:
        i_row, i_column = count / n_columns, count % n_columns
        data_idx = data[idx]

        plot_with_factor = plot_with_factor_input
        if not plot_with_factor:
            data_idx['plot-factor'] = np.ones(len(data_idx['value']))
        if all(_ == 1.0 for _ in data_idx['plot-factor']):
            ## this counts for the case where there is no plot factor from the experimental group
            ## value of plot factor is all set to 1.0 in such case
            plot_with_factor = False

        if plot_with_factor:
            ## find the base of the plot factors
            plot_factor = data_idx['plot-factor'][0]
            for i in [2, 3, 5, 6, 7, 10]:
                plot_factor_power = np.log(plot_factor) / np.log(i)
                if abs(round(plot_factor_power) - plot_factor_power) < 1e-10:
                    plot_factor_base = i
                    break

        _, absolute_eta, eta_bins = find_eta_bins(data_idx)

        for i in xrange(len(eta_bins)):
            eta_range = ''
            if plot_with_factor:
                if absolute_eta:
                    eta_min, eta_max = str(data_idx['eta-abs-min'][eta_bins[i][0] - 1]), str(data_idx['eta-abs-max'][eta_bins[i][0] - 1])
                    plot_factor_power = round(np.log(data_idx['plot-factor'][eta_bins[i][0] - 1]) / np.log(plot_factor_base))
                    eta_range += r'$|\eta|~\in~[%s,~%s]~(\times %i^{%i})$' % (eta_min, eta_max, plot_factor_base, plot_factor_power)
                else:
                    eta_min, eta_max = str(data_idx['eta-min'][eta_bins[i][0] - 1]), str(data_idx['eta-max'][eta_bins[i][0] - 1])
                    plot_factor_power = round(np.log(data_idx['plot-factor'][eta_bins[i][0] - 1]) / np.log(plot_factor_base))
                    eta_range += r'$\eta~\in~[%s,~%s]~(\times %i^{%i})$' % (eta_min, eta_max, plot_factor_base, plot_factor_power)
            else:
                if absolute_eta:
                    eta_min, eta_max = str(data_idx['eta-abs-min'][eta_bins[i][0] - 1]), str(data_idx['eta-abs-max'][eta_bins[i][0] - 1])
                    eta_range += r'$|\eta|~\in~[%s,~%s]$' % (eta_min, eta_max)
                else:
                    eta_min, eta_max = str(data_idx['eta-min'][eta_bins[i][0] - 1]), str(data_idx['eta-max'][eta_bins[i][0] - 1])
                    eta_range += r'$\eta~\in~[%s,~%s]$' % (eta_min, eta_max)

            pt = np.array(data_idx['pT'][eta_bins[i][0] - 1 : eta_bins[i][1]])
            value = np.array(data_idx['plot-factor'][eta_bins[i][0] - 1 : eta_bins[i][1]]) * \
                np.array(data_idx['value'][eta_bins[i][0] - 1 : eta_bins[i][1]])
            alpha = np.array(data_idx['plot-factor'][eta_bins[i][0] - 1 : eta_bins[i][1]]) * \
                np.array(data_idx['alpha'][eta_bins[i][0] - 1 : eta_bins[i][1]])
            ax[i_row, i_column].errorbar(pt, value, alpha, fmt = '%s.' % plot_colors[i])
            theory = np.array(data_idx['plot-factor'][eta_bins[i][0] - 1 : eta_bins[i][1]]) * \
                np.array(data_idx['thy'][eta_bins[i][0] - 1 : eta_bins[i][1]])
            ax[i_row, i_column].plot(pt, theory, '%s%s' % (plot_colors[i], plot_styles[0]), label = eta_range)
            ax[i_row, i_column].axhline(0.0, color = 'black', ls = ':')

        ax[i_row, i_column].set_yscale('log')
        ax[i_row, i_column].relim()
        ax[i_row, i_column].autoscale_view()
        if (idx != 10001) and (idx != 10002):
            ax[i_row, i_column].legend(fontsize = 10, loc = 'best')
        if idx == 10001:
            year = 'run II'
        elif idx == 10002:
            year = 'run II'
        elif idx == 10003:
            year = 'MB'
        elif idx == 10004:
            year = 'HT'
        ax[i_row, i_column].title.set_text(r'$\mathrm{%s~%s}$' % (data_idx['col'][0].upper(), year.replace(' ', '~')))
        count += 1

    ax[1, 0].set_xlabel(r'$p_T~(\mathrm{GeV})$', size = 15)
    ax[1, 1].set_xlabel(r'$p_T~(\mathrm{GeV})$', size = 15)
    # plot_temp.xaxis.set_label_coords(1.05, 0.2)
    ax[0, 0].set_ylabel(r'$\frac{\mathrm{d}^2\sigma}{\mathrm{d}y \mathrm{d} p_T}~{\scriptstyle [\mathrm{nb} / (\mathrm{GeV} / c)]}$', size = 17)
    ## the units of D0 and CDF can be different
    ax[1, 0].set_ylabel(r'$\frac{\mathrm{d}^2\sigma}{2 \pi \mathrm{d}y \mathrm{d} p_T}~{\scriptstyle [\mathrm{nb} / (\mathrm{GeV} / c)]}$', size = 17)
    ## y labeling has to be more automatic

    py.tight_layout()
    figure.savefig('%s/gallery/jet-%d.png' % (wdir, istep), dpi = dpi)

def plot_data_separate(wdir, data_idx, istep, idx, dpi, plot_with_factor):
    if not plot_with_factor:
        data_idx['plot-factor'] = np.ones(len(data_idx['value']))
    if all(_ == 1.0 for _ in data_idx['plot-factor']):
        ## this counts for the case where there is no plot factor from the experimental group
        ## value of plot factor is all set to 1.0 in such case
        plot_with_factor = False

    if plot_with_factor:
        ## find the base of the plot factors
        plot_factor = data_idx['plot-factor'][0]
        for i in [2, 3, 5, 6, 7, 10]:
            plot_factor_power = np.log(plot_factor) / np.log(i)
            if abs(round(plot_factor_power) - plot_factor_power) < 1e-10:
                plot_factor_base = i
                break

    _, absolute_eta, eta_bins = find_eta_bins(data_idx)

    figure = py.figure(figsize = [10.0, 5.0], dpi = 500)
    plot_temp = figure.add_subplot(1, 1, 1)

    plot_styles = ['-']
    plot_colors = ['b', 'g', 'm', 'y', 'c', 'r']
    # plot_styles = ['-', '-.', '--', ':', '-', '--']
    for i in xrange(len(eta_bins)):
        eta_range = ''
        if plot_with_factor:
            if absolute_eta:
                eta_min, eta_max = str(data_idx['eta-abs-min'][eta_bins[i][0] - 1]), str(data_idx['eta-abs-max'][eta_bins[i][0] - 1])
                plot_factor_power = round(np.log(data_idx['plot-factor'][eta_bins[i][0] - 1]) / np.log(plot_factor_base))
                eta_range += r'$|\eta|~\in~[%s,~%s]~(\times %i^{%i})$' % (eta_min, eta_max, plot_factor_base, plot_factor_power)
            else:
                eta_min, eta_max = str(data_idx['eta-min'][eta_bins[i][0] - 1]), str(data_idx['eta-max'][eta_bins[i][0] - 1])
                plot_factor_power = round(np.log(data_idx['plot-factor'][eta_bins[i][0] - 1]) / np.log(plot_factor_base))
                eta_range += r'$\eta~\in~[%s,~%s]~(\times %i^{%i})$' % (eta_min, eta_max, plot_factor_base, plot_factor_power)
        else:
            if absolute_eta:
                eta_min, eta_max = str(data_idx['eta-abs-min'][eta_bins[i][0] - 1]), str(data_idx['eta-abs-max'][eta_bins[i][0] - 1])
                eta_range += r'$|\eta|~\in~[%s,~%s]$' % (eta_min, eta_max)
            else:
                eta_min, eta_max = str(data_idx['eta-min'][eta_bins[i][0] - 1]), str(data_idx['eta-max'][eta_bins[i][0] - 1])
                eta_range += r'$\eta~\in~[%s,~%s]$' % (eta_min, eta_max)

        pt = np.array(data_idx['pT'][eta_bins[i][0] - 1 : eta_bins[i][1]])
        value = np.array(data_idx['plot-factor'][eta_bins[i][0] - 1 : eta_bins[i][1]]) * \
            np.array(data_idx['value'][eta_bins[i][0] - 1 : eta_bins[i][1]])
        alpha = np.array(data_idx['plot-factor'][eta_bins[i][0] - 1 : eta_bins[i][1]]) * \
            np.array(data_idx['alpha'][eta_bins[i][0] - 1 : eta_bins[i][1]])
        plot_temp.errorbar(pt, value, alpha, fmt = '%s.' % plot_colors[i])
        theory = np.array(data_idx['plot-factor'][eta_bins[i][0] - 1 : eta_bins[i][1]]) * \
            np.array(data_idx['thy'][eta_bins[i][0] - 1 : eta_bins[i][1]])
        plot_temp.plot(pt, theory, '%s%s' % (plot_colors[i], plot_styles[0]), label = eta_range)

    plot_temp.legend(fontsize = 10, loc = 'best')

    plot_temp.set_xlabel(r'$p_T~(\mathrm{GeV})$', size = 15)
    # plot_temp.xaxis.set_label_coords(1.05, 0.2)
    if data_idx['obs'][0].replace('<', '').replace('>', '') == 'd2_sigma_over_d_y_d_pt':
        if data_idx['units'][0] == 'nb':
            plot_temp.set_ylabel(r'$\frac{\mathrm{d}^2\sigma}{\mathrm{d}y \mathrm{d} p_T}~{\scriptstyle [\mathrm{nb} / (\mathrm{GeV} / c)]}$', size = 17)
        elif data_idx['units'][0] == 'pb':
            plot_temp.set_ylabel(r'$\frac{\mathrm{d}^2\sigma}{\mathrm{d}y \mathrm{d} p_T}~{\scriptstyle [\mathrm{pb} / (\mathrm{GeV} / c)]}$', size = 17)
    elif data_idx['obs'][0].replace('<', '').replace('>', '') == 'd2_sigma_over_2_pi_d_y_d_pt':
        if data_idx['units'][0] == 'nb':
            plot_temp.set_ylabel(r'$\frac{\mathrm{d}^2\sigma}{2 \pi \mathrm{d}y \mathrm{d} p_T}~{\scriptstyle [\mathrm{nb} / (\mathrm{GeV} / c)]}$', size = 17)
        elif data_idx['units'][0] == 'pb':
            plot_temp.set_ylabel(r'$\frac{\mathrm{d}^2\sigma}{2 \pi \mathrm{d}y \mathrm{d} p_T}~{\scriptstyle [\mathrm{pb} / (\mathrm{GeV} / c)]}$', size = 17)
    plot_temp.set_yscale('log')

    if data_idx['col'][0] == 'star':
        if idx == 10003:
            plot_temp.set_title(r'$\mathrm{STAR~MB}$')
            py.tight_layout()
            figure.savefig('%s/gallery/jet-%s-MB-%d.png' % (wdir, data_idx['col'][0], istep), dpi = dpi)
        elif idx == 10004:
            plot_temp.set_title(r'$\mathrm{STAR~HT}$')
            py.tight_layout()
            figure.savefig('%s/gallery/jet-%s-HT-%d.png' % (wdir, data_idx['col'][0], istep), dpi = dpi)
    else:
        plot_temp.set_title(r'$\mathrm{%s}$' % data_idx['col'][0].upper())
        py.tight_layout()
        figure.savefig('%s/gallery/jet-%s-%d.png' % (wdir, data_idx['col'][0], istep), dpi = dpi)

def plot_data_on_theory(wdir, data, istep, dpi):
    if (10003 in data) and (10004 in data):
        nrows, ncols = 5, 3
    else:
        nrows, ncols = 4, 3
    fig = py.figure(figsize = (ncols * 4.5, nrows * 3.0))
    count = 0
    ax = {}
    eta_range = {}

    for key in sorted(data):
        _, absolute_eta, eta_bins = find_eta_bins(data[key])
        for i in range(len(eta_bins)):
            count += 1
            if absolute_eta:
                eta_min, eta_max = data[key]['eta-abs-min'][eta_bins[i][0] - 1], data[key]['eta-abs-max'][eta_bins[i][0] - 1]
            else:
                eta_min, eta_max = data[key]['eta-min'][eta_bins[i][0] - 1], data[key]['eta-max'][eta_bins[i][0] - 1]

            ax[count] = py.subplot(nrows, ncols, count)
            eta_range[count] = ''
            theory = data[key]['thy'][eta_bins[i][0] - 1 : eta_bins[i][1]]
            pt = data[key]['pT'][eta_bins[i][0] - 1 : eta_bins[i][1]]
            value = data[key]['value'][eta_bins[i][0] - 1 : eta_bins[i][1]]
            alpha = data[key]['alpha'][eta_bins[i][0] - 1 : eta_bins[i][1]]

            if absolute_eta:
                eta_range[count] += r'$|\eta|~\in~[%s,~%s]$' % (str(eta_min), str(eta_max))
            else:
                eta_range[count] += r'$\eta~\in~[%s,~%s]$' % (str(eta_min), str(eta_max))

            if key == 10001:
                ax[count].errorbar(pt, value / theory, alpha / theory, color = 'firebrick', marker = '.', linestyle = 'none', label = r'$\mathrm{%s}$' % data[key]['col'][0].upper())
            elif key == 10002:
                ax[count].errorbar(pt, value / theory, alpha / theory, color = 'darkorange', marker = 'v', linestyle = 'none', label = r'$\mathrm{%s}$' % data[key]['col'][0].upper())
            elif key == 10003:
                ax[count].errorbar(pt, value / theory, alpha / theory, color = 'olivedrab', marker = 's', linestyle = 'none', label = r'$\mathrm{%s}$' % (data[key]['col'][0].upper() + '~MB'))
            elif key == 10004:
                ax[count].errorbar(pt, value / theory, alpha / theory, color = 'darkcyan', marker = '*', linestyle = 'none', label = r'$\mathrm{%s}$' % (data[key]['col'][0].upper() + '~HT'))
        if key == 10002: count += 1 ## make a white space to fill up a row for CDF

    ax[1].legend(fontsize = 20, loc = 'upper left')
    ax[4].legend(fontsize = 20, loc = 'upper center')
    ax[7].legend(fontsize = 20, loc = 'upper left')
    ax[10].legend(fontsize = 20, loc = 'upper center')
    if (10003 in data) and (10004 in data):
        ax[13].legend(fontsize = 20, loc = 'upper center')
        ax[14].legend(fontsize = 20, loc = 'upper center')

    minorLocator = MultipleLocator(0.05)
    majorLocator = MultipleLocator(0.2)
    for i in ax:
        ax[i].axhline(1.0, color = 'black', ls = ':')
        ax[i].set_ylim(0.8, 1.2)
        ax[i].yaxis.set_minor_locator(minorLocator)
        ax[i].yaxis.set_major_locator(majorLocator)

    ax[1].text(100.0, 0.85, eta_range[1], fontsize = 22)
    ax[2].text(50.0, 0.85, eta_range[2], fontsize = 22)
    ax[3].text(100.0, 0.85, eta_range[3], fontsize = 22)
    ax[4].text(100.0, 0.85, eta_range[4], fontsize = 22)
    ax[5].text(75.0, 0.85, eta_range[5], fontsize = 22)
    ax[6].text(75.0, 0.85, eta_range[6], fontsize = 22)

    ax[7].text(200.0, 0.85, eta_range[7], fontsize = 22)
    ax[8].text(100.0, 0.85, eta_range[8], fontsize = 22)
    ax[9].text(100.0, 0.85, eta_range[9], fontsize = 22)
    ax[10].text(100.0, 0.85, eta_range[10], fontsize = 22)
    ax[11].text(100.0, 0.85, eta_range[11], fontsize = 22)

    if (10003 in data) and (10004 in data):
        ax[13].text(11.0, 0.9, eta_range[13], fontsize = 22)
        ax[14].text(20.0, 0.85, eta_range[14], fontsize = 22)

    ax[1].set_ylabel(r'\boldmath$\mathrm{data/theory}$', size = 24)
    ax[1].yaxis.set_label_coords(-0.17, 0.5)
    ax[4].set_ylabel(r'\boldmath$\mathrm{data/theory}$', size = 24)
    ax[4].yaxis.set_label_coords(-0.17, 0.5)
    ax[7].set_ylabel(r'\boldmath$\mathrm{data/theory}$', size = 24)
    ax[7].yaxis.set_label_coords(-0.17, 0.5)
    ax[10].set_ylabel(r'\boldmath$\mathrm{data/theory}$', size = 24)
    ax[10].yaxis.set_label_coords(-0.17, 0.5)
    if (10003 in data) and (10004 in data):
        ax[13].set_ylabel(r'\boldmath$\mathrm{data/theory}$', size = 24)
        ax[13].yaxis.set_label_coords(-0.17, 0.5)

    if (10003 in data) and (10004 in data):
        ax[9].set_xlabel(r'\boldmath$p_{T}~(\mathrm{GeV})$', size = 24)
        ax[9].xaxis.set_label_coords(0.90, -0.12)
        ax[13].set_xlabel(r'\boldmath$p_{T}~(\mathrm{GeV})$', size = 24)
        ax[13].xaxis.set_label_coords(0.90, -0.12)
        ax[14].set_xlabel(r'\boldmath$p_{T}~(\mathrm{GeV})$', size = 24)
        ax[14].xaxis.set_label_coords(0.90, -0.12)
    else:
        ax[9].set_xlabel(r'\boldmath$p_{T}~(\mathrm{GeV})$', size = 24)
        ax[9].xaxis.set_label_coords(0.90, -0.12)
        ax[10].set_xlabel(r'\boldmath$p_{T}~(\mathrm{GeV})$', size = 24)
        ax[10].xaxis.set_label_coords(0.90, -0.12)
        ax[11].set_xlabel(r'\boldmath$p_{T}~(\mathrm{GeV})$', size = 24)
        ax[11].xaxis.set_label_coords(0.90, -0.12)

    py.tight_layout()
    # py.subplots_adjust(left = 0.08, bottom = 0.08, right = 0.99, top = 0.97, wspace = 0.2, hspace = 0.1)
    py.savefig('%s/gallery/jet-ratio-%d.png' % (wdir, istep), dpi = dpi)
    py.close()

def make_figure(wdir, task, plot_with_factor = True, only_best_cluster = True, dpi = 200):

    if task == 1:
        print('\nplotting JET data from %s' % (wdir))
    elif task == 2:
        print('\nplotting JET data over theory from %s' % (wdir))

    load_config('%s/input.py' % wdir)
    istep = core.get_istep()
    replicas = core.get_replicas(wdir)
    core.mod_conf(istep, replicas[0]) #--set conf as specified in istep

    predictions = load('%s/data/predictions-%d.dat' % (wdir, istep))
    if 'jet' not in predictions['reactions']:
        print('JET is not in data file')
        return
    labels  = load('%s/data/labels-%d.dat' % (wdir, istep))
    cluster = labels['cluster']
    cluster_average_chi2 = labels['cluster_average_chi2']
    best_cluster = np.argmin(cluster_average_chi2)

    data = predictions['reactions']['jet']

    for idx in data:
        predictions = copy.copy(data[idx]['prediction-rep'])
        del data[idx]['prediction-rep']
        del data[idx]['residuals-rep']
        if only_best_cluster:
            best_predictions = [predictions[i] for i in range(len(predictions)) if cluster[i] == best_cluster]
            data[idx]['thy'] = np.mean(best_predictions, axis = 0)
            data[idx]['dthy'] = np.std(best_predictions, axis = 0) ** 0.5
        else:
            all_predictions = [predictions[i] for i in range(len(predictions))]
            data[idx]['thy'] = np.mean(all_predictions, axis = 0)
            data[idx]['dthy'] = np.std(all_predictions, axis = 0) ** 0.5

    if task == 1:
        plot_data(wdir, data, istep, dpi, plot_with_factor)
        # for idx in data:
        #     plot_data_separate(wdir, data[idx], istep, idx, dpi, plot_with_factor)
    elif task == 2:
        plot_data_on_theory(wdir, data, istep, dpi)

    return

if __name__ == "__main__":
    pass
