#!/usr/bin/env python
import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde
from scipy.integrate import quad
from mpmath import fp
import tempfile

## matplotlib
import matplotlib
matplotlib.use('Agg')
import pylab as py
import matplotlib.gridspec as gridspec
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Times-Roman']})
rc('text',usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
from matplotlib.ticker import FixedLocator

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint, tex
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
from analysis.corelib import classifier

def generate_maps():
    maps = {}

    maps['experiment'] = {}
    maps['experiment']['EMC']=tex('EMC')
    maps['experiment']['SMC']=tex('SMC')
    maps['experiment']['SLACE80E130']=tex('E80/E130')
    maps['experiment']['SLAC(E142)']=tex('E142')
    maps['experiment']['SLAC(E143)']=tex('E143')
    maps['experiment']['SLAC(E154)']=tex('E154')
    maps['experiment']['SLAC(E155)']=tex('E155')
    maps['experiment']['SLAC(E155x)']=tex('E155x')
    maps['experiment']['HERMES']=tex('HERMES')
    maps['experiment']['COMPASS']=tex('COMPASS')
    maps['experiment']['JLabHA(E06014)']=tex('E06014')
    maps['experiment']['JLabHA(E99117)']=tex('E99117')
    maps['experiment']['JLabHB(EG1DVCS)']=tex('eg1')+'-'+tex('DVCS')
    maps['experiment']['JLabHB(EG1b)']=tex('eg1b')

    maps['target'] = {}
    maps['target']['proton']='p'
    maps['target']['neutron']='n'
    maps['target']['helium']='He'
    maps['target']['deuteron']='d'

    maps['observable'] = {}
    maps['observable']['A1']  = r'$A_1^{\mathrm{%s}}$'
    maps['observable']['A2']  = r'$A_2^{\mathrm{%s}}$'
    maps['observable']['Apa'] = r'$A_\parallel^{\mathrm{%s}}$'
    maps['observable']['Ape'] = r'$A_\perp^{\mathrm{%s}}$'
    maps['observable']['Atpe']= r'$\widetilde{A}_\perp^{\mathrm{%s}}$'

    maps['color'] = {}
    maps['color']['p']='r'
    maps['color']['He']='b'
    maps['color']['d']='g'

    return maps

def pick_subsets(data):
    set_headers = {}
    for dataset in data:
        set_heater_temp = [k for k in data[dataset].keys() if 'set' in k]
        if any(set_heater_temp): set_headers[dataset] = set_heater_temp[0]

    subset_bins = {}
    for dataset in data:
        if dataset in set_headers:
            subset_bins[dataset] = {}
            set_header = set_headers[dataset]
            subsets = [int(i) for i in data[dataset][set_header]] ## some sets can be absent due to 'Q2' cut and 'W2' cut
            subsets = sorted(set(subsets))
            for subset in subsets:
                subset_temp = []
                for i in range(len(data[dataset][set_header])):
                    if data[dataset][set_header][i] == subset: subset_temp.append(i)
                subset_bins[dataset][subset] = subset_temp
        else:
            subset_bins[dataset] = {}
            subset_bins[dataset][0] = range(len(data[dataset]['X']))

    return subset_bins

def data_Apa_proton(wdir, data, indices, subset_bins, maps, istep, dpi):
    n_columns, n_rows = 4, 5
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.7))
    axs = {}
    count = 0
    for index in indices:
        correlation = [_ for _ in data[index] if ('_c' in _) and ('norm' not in _)]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if (data[index]['col'] == 'SLACE80E130') and (i_set == 5): continue
            if len(idx) <= 2: continue
            count += 1
            ax = py.subplot(n_rows, n_columns, count)
            axs[count] = ax

            xs_temp = data[index]['X'][idx]
            idx_order = np.argsort(xs_temp)
            xs = data[index]['X'][idx][idx_order]
            Q2s = data[index]['Q2'][idx][idx_order]
            experimental_values = data[index]['value'][idx][idx_order]

            alphas = data[index]['alpha'][idx][idx_order]
            shifts = data[index]['shift'][idx][idx_order]

            theory_mean = data[index]['thy'][idx][idx_order]
            theory_deviation = data[index]['dthy'][idx][idx_order]

            ax.errorbar(xs, experimental_values, yerr = alphas, fmt = 'k.', markersize = 5, capsize = 0, alpha = 0.8)
            ax.plot(xs, theory_mean, 'r-')
            ax.fill_between(xs, theory_mean - theory_deviation, theory_mean + theory_deviation, color = 'r', alpha = 0.5)

            if len(xs) > 6:
                ax.semilogx()

            if any(correlation):
                pass
                experimental_bottom = min(experimental_values - alphas)
                experimental_top = max(experimental_values - alphas)
                experimental_bottom = experimental_bottom - 0.2 * (experimental_top - experimental_bottom)
                correlation_center = np.full(len(alphas), experimental_bottom)
                correlation_positive, correlation_negative = np.zeros(len(alphas)), np.zeros(len(alphas))
                for i in range(len(shifts)):
                    if shifts[i] > 0.0: correlation_positive[i] = shifts[i]
                    else: correlation_negative[i] = shifts[i]
                ax.fill_between(xs, correlation_center, correlation_center + correlation_positive, color = 'g', alpha = 0.5)
                ax.fill_between(xs, correlation_center, correlation_center + correlation_negative, color = 'b', alpha = 0.5)

            # matplotlib.pyplot.locator_params(axis = 'x', nbins = 2)
            # ax.locator_params(axis = 'x', tight = True)
            # ax.set_xticks(ax.get_xticks()[::2])
            # matplotlib.ticker.LogLocator(base = 10.0, subs = 1.0, numdecs = 4, numticks = 3)
            observable = maps['observable'][np.array(data[index]['obs'])[idx][0]] % np.array(data[index]['target'])[idx][0]
            title = maps['experiment'][np.array(data[index]['col'])[idx][0]]
            ax.title.set_text(title + r'$,~$' + observable)
            ax.text(0.07, 0.77, r'$Q^2 \in \left[%0.1f,~%0.1f \right]$' % (min(Q2s), max(Q2s)), transform = ax.transAxes)

    x_label_counts = range(count, count - n_columns, -1)
    for x_label_count in x_label_counts:
        axs[x_label_count].set_xlabel(r'$x$', size = 15)

    py.tight_layout()
    py.savefig('%s/gallery/pidis-apa-proton-%d.png' % (wdir, istep), dpi = dpi)
    return

def data_Ape_proton(wdir, data, indices, subset_bins, maps, istep, dpi):
    n_columns, n_rows = 4, 3
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.7))
    axs = {}
    count = 0
    for index in indices:
        correlation = [_ for _ in data[index] if ('_c' in _) and ('norm' not in _)]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if i_set == 5: continue
            if len(idx) <= 2: continue
            count += 1
            ax = py.subplot(n_rows, n_columns, count)
            axs[count] = ax

            xs_temp = data[index]['X'][idx]
            idx_order = np.argsort(xs_temp)
            xs = data[index]['X'][idx][idx_order]
            Q2s = data[index]['Q2'][idx][idx_order]
            experimental_values = data[index]['value'][idx][idx_order]

            alphas = data[index]['alpha'][idx][idx_order]
            shifts = data[index]['shift'][idx][idx_order]

            theory_mean = data[index]['thy'][idx][idx_order]
            theory_deviation = data[index]['dthy'][idx][idx_order]

            ax.errorbar(xs, experimental_values, yerr = alphas, fmt = 'k.', markersize = 5, capsize = 0, alpha = 0.8)
            ax.plot(xs, theory_mean, 'r-')
            ax.fill_between(xs, theory_mean - theory_deviation, theory_mean + theory_deviation, color = 'r', alpha = 0.5)

            if len(xs) > 6:
                ax.semilogx()

            if any(correlation):
                pass
                experimental_bottom = min(experimental_values - alphas)
                experimental_top = max(experimental_values - alphas)
                experimental_bottom = experimental_bottom - 0.2 * (experimental_top - experimental_bottom)
                correlation_center = np.full(len(alphas), experimental_bottom)
                correlation_positive, correlation_negative = np.zeros(len(alphas)), np.zeros(len(alphas))
                for i in range(len(shifts)):
                    if shifts[i] > 0.0: correlation_positive[i] = shifts[i]
                    else: correlation_negative[i] = shifts[i]
                ax.fill_between(xs, correlation_center, correlation_center + correlation_positive, color = 'g', alpha = 0.5)
                ax.fill_between(xs, correlation_center, correlation_center + correlation_negative, color = 'b', alpha = 0.5)

            observable = maps['observable'][np.array(data[index]['obs'])[idx][0]] % np.array(data[index]['target'])[idx][0]
            title = maps['experiment'][np.array(data[index]['col'])[idx][0]]
            ax.title.set_text(title + r'$,~$' + observable)
            ax.text(0.07, 0.77, r'$Q^2 \in \left[%0.1f,~%0.1f \right]$' % (min(Q2s), max(Q2s)), transform = ax.transAxes)

    x_label_counts = range(count, count - n_columns, -1)
    for x_label_count in x_label_counts:
        axs[x_label_count].set_xlabel(r'$x$', size = 15)

    py.tight_layout()
    py.savefig('%s/gallery/pidis-ape-proton-%d.png' % (wdir, istep), dpi = dpi)
    return

def data_Apa_proton_DVCS(wdir, data, indices, subset_bins, maps, istep, dpi):
    n_columns, n_rows = 4, 3
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.7))
    axs = {}
    count = 0
    for index in indices:
        correlation = [_ for _ in data[index] if ('_c' in _) and ('norm' not in _)]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if len(idx) <= 2: continue
            count += 1
            ax = py.subplot(n_rows, n_columns, count)
            axs[count] = ax

            xs_temp = data[index]['X'][idx]
            idx_order = np.argsort(xs_temp)
            xs = data[index]['X'][idx][idx_order]
            Q2s = data[index]['Q2'][idx][idx_order]
            energys = data[index]['E'][idx][idx_order]
            experimental_values = data[index]['value'][idx][idx_order]

            alphas = data[index]['alpha'][idx][idx_order]
            shifts = data[index]['shift'][idx][idx_order]

            theory_mean = data[index]['thy'][idx][idx_order]
            theory_deviation = data[index]['dthy'][idx][idx_order]

            ax.errorbar(xs, experimental_values, yerr = alphas, fmt = 'k.', markersize = 5, capsize = 0, alpha = 0.8)
            ax.plot(xs, theory_mean, 'r-')
            ax.fill_between(xs, theory_mean - theory_deviation, theory_mean + theory_deviation, color = 'r', alpha = 0.5)

            if len(xs) > 6:
                ax.semilogx()

            if any(correlation):
                pass
                experimental_bottom = min(experimental_values - alphas)
                experimental_top = max(experimental_values - alphas)
                experimental_bottom = experimental_bottom - 0.2 * (experimental_top - experimental_bottom)
                correlation_center = np.full(len(alphas), experimental_bottom)
                correlation_positive, correlation_negative = np.zeros(len(alphas)), np.zeros(len(alphas))
                for i in range(len(shifts)):
                    if shifts[i] > 0.0: correlation_positive[i] = shifts[i]
                    else: correlation_negative[i] = shifts[i]
                ax.fill_between(xs, correlation_center, correlation_center + correlation_positive, color = 'g', alpha = 0.5)
                ax.fill_between(xs, correlation_center, correlation_center + correlation_negative, color = 'b', alpha = 0.5)

            observable = maps['observable'][np.array(data[index]['obs'])[idx][0]] % np.array(data[index]['target'])[idx][0]
            title = maps['experiment'][np.array(data[index]['col'])[idx][0]]
            ax.title.set_text(title + r'$,~$' + observable)
            ax.text(0.07, 0.77, r'$E = %.1f,~Q^2 \in \left[%0.1f,~%0.1f \right]$' % (energys[0], min(Q2s), max(Q2s)), transform = ax.transAxes)

    x_label_counts = range(count, count - n_columns, -1)
    for x_label_count in x_label_counts:
        axs[x_label_count].set_xlabel(r'$x$', size = 15)

    py.tight_layout()
    py.savefig('%s/gallery/pidis-apa-proton-DVCS-%d.png' % (wdir, istep), dpi = dpi)
    return

def data_Apa_proton_EG1b(wdir, data, indices, subset_bins, maps, istep, dpi):
    n_columns, n_rows = 4, 3
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.7))
    axs = {}
    count = 0
    for index in indices:
        correlation = [_ for _ in data[index] if ('_c' in _) and ('norm' not in _)]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if len(idx) <= 2: continue
            count += 1
            ax = py.subplot(n_rows, n_columns, count)
            axs[count] = ax

            xs_temp = data[index]['X'][idx]
            idx_order = np.argsort(xs_temp)
            xs = data[index]['X'][idx][idx_order]
            Q2s = data[index]['Q2'][idx][idx_order]
            energys = data[index]['E'][idx][idx_order]
            experimental_values = data[index]['value'][idx][idx_order]

            alphas = data[index]['alpha'][idx][idx_order]
            shifts = data[index]['shift'][idx][idx_order]

            theory_mean = data[index]['thy'][idx][idx_order]
            theory_deviation = data[index]['dthy'][idx][idx_order]

            ax.errorbar(xs, experimental_values, yerr = alphas, fmt = 'k.', markersize = 5, capsize = 0, alpha = 0.8)
            ax.plot(xs, theory_mean, 'r-')
            ax.fill_between(xs, theory_mean - theory_deviation, theory_mean + theory_deviation, color = 'r', alpha = 0.5)

            if len(xs) > 6:
                ax.semilogx()

            if any(correlation):
                pass
                experimental_bottom = min(experimental_values - alphas)
                experimental_top = max(experimental_values - alphas)
                experimental_bottom = experimental_bottom - 0.2 * (experimental_top - experimental_bottom)
                correlation_center = np.full(len(alphas), experimental_bottom)
                correlation_positive, correlation_negative = np.zeros(len(alphas)), np.zeros(len(alphas))
                for i in range(len(shifts)):
                    if shifts[i] > 0.0: correlation_positive[i] = shifts[i]
                    else: correlation_negative[i] = shifts[i]
                ax.fill_between(xs, correlation_center, correlation_center + correlation_positive, color = 'g', alpha = 0.5)
                ax.fill_between(xs, correlation_center, correlation_center + correlation_negative, color = 'b', alpha = 0.5)

            observable = maps['observable'][np.array(data[index]['obs'])[idx][0]] % np.array(data[index]['target'])[idx][0]
            title = maps['experiment'][np.array(data[index]['col'])[idx][0]]
            ax.title.set_text(title + r'$,~$' + observable)
            ax.text(0.07, 0.77, r'$E = %.1f,~Q^2 \in \left[%0.1f,~%0.1f \right]$' % (energys[0], min(Q2s), max(Q2s)), transform = ax.transAxes)

    x_label_counts = range(count, count - n_columns, -1)
    for x_label_count in x_label_counts:
        axs[x_label_count].set_xlabel(r'$x$', size = 15)

    py.tight_layout()
    py.savefig('%s/gallery/pidis-apa-proton-EG1b-%d.png' % (wdir, istep), dpi = dpi)
    return

def data_Apa_deuteron(wdir, data, indices, subset_bins, maps, istep, dpi):
    n_columns, n_rows = 4, 3
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.7))
    axs = {}
    count = 0
    for index in indices:
        correlation = [_ for _ in data[index] if ('_c' in _) and ('norm' not in _)]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if len(idx) <= 2: continue
            count += 1
            ax = py.subplot(n_rows, n_columns, count)
            axs[count] = ax

            xs_temp = data[index]['X'][idx]
            idx_order = np.argsort(xs_temp)
            xs = data[index]['X'][idx][idx_order]
            Q2s = data[index]['Q2'][idx][idx_order]
            experimental_values = data[index]['value'][idx][idx_order]

            alphas = data[index]['alpha'][idx][idx_order]
            shifts = data[index]['shift'][idx][idx_order]

            theory_mean = data[index]['thy'][idx][idx_order]
            theory_deviation = data[index]['dthy'][idx][idx_order]

            ax.errorbar(xs, experimental_values, yerr = alphas, fmt = 'k.', markersize = 5, capsize = 0, alpha = 0.8)
            ax.plot(xs, theory_mean, 'r-')
            ax.fill_between(xs, theory_mean - theory_deviation, theory_mean + theory_deviation, color = 'r', alpha = 0.5)

            if len(xs) > 6:
                ax.semilogx()

            if any(correlation):
                pass
                experimental_bottom = min(experimental_values - alphas)
                experimental_top = max(experimental_values - alphas)
                experimental_bottom = experimental_bottom - 0.2 * (experimental_top - experimental_bottom)
                correlation_center = np.full(len(alphas), experimental_bottom)
                correlation_positive, correlation_negative = np.zeros(len(alphas)), np.zeros(len(alphas))
                for i in range(len(shifts)):
                    if shifts[i] > 0.0: correlation_positive[i] = shifts[i]
                    else: correlation_negative[i] = shifts[i]
                ax.fill_between(xs, correlation_center, correlation_center + correlation_positive, color = 'g', alpha = 0.5)
                ax.fill_between(xs, correlation_center, correlation_center + correlation_negative, color = 'b', alpha = 0.5)

            observable = maps['observable'][np.array(data[index]['obs'])[idx][0]] % np.array(data[index]['target'])[idx][0]
            title = maps['experiment'][np.array(data[index]['col'])[idx][0]]
            ax.title.set_text(title + r'$,~$' + observable)
            ax.text(0.07, 0.77, r'$Q^2 \in \left[%0.1f,~%0.1f \right]$' % (min(Q2s), max(Q2s)), transform = ax.transAxes)

    x_label_counts = range(count, count - n_columns, -1)
    for x_label_count in x_label_counts:
        axs[x_label_count].set_xlabel(r'$x$', size = 15)

    py.tight_layout()
    py.savefig('%s/gallery/pidis-apa-deuteron-%d.png' % (wdir, istep), dpi = dpi)
    return

def data_Ape_deuteron(wdir, data, indices, subset_bins, maps, istep, dpi):
    n_columns, n_rows = 4, 3
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.7))
    axs = {}
    count = 0
    for index in indices:
        correlation = [_ for _ in data[index] if ('_c' in _) and ('norm' not in _)]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if len(idx) <= 2: continue
            count += 1
            ax = py.subplot(n_rows, n_columns, count)
            axs[count] = ax

            xs_temp = data[index]['X'][idx]
            idx_order = np.argsort(xs_temp)
            xs = data[index]['X'][idx][idx_order]
            Q2s = data[index]['Q2'][idx][idx_order]
            # energys = data[index]['E'][idx][idx_order]
            experimental_values = data[index]['value'][idx][idx_order]

            alphas = data[index]['alpha'][idx][idx_order]
            shifts = data[index]['shift'][idx][idx_order]

            theory_mean = data[index]['thy'][idx][idx_order]
            theory_deviation = data[index]['dthy'][idx][idx_order]

            ax.errorbar(xs, experimental_values, yerr = alphas, fmt = 'k.', markersize = 5, capsize = 0, alpha = 0.8)
            ax.plot(xs, theory_mean, 'r-')
            ax.fill_between(xs, theory_mean - theory_deviation, theory_mean + theory_deviation, color = 'r', alpha = 0.5)

            if len(xs) > 6:
                ax.semilogx()

            if any(correlation):
                pass
                experimental_bottom = min(experimental_values - alphas)
                experimental_top = max(experimental_values - alphas)
                experimental_bottom = experimental_bottom - 0.2 * (experimental_top - experimental_bottom)
                correlation_center = np.full(len(alphas), experimental_bottom)
                correlation_positive, correlation_negative = np.zeros(len(alphas)), np.zeros(len(alphas))
                for i in range(len(shifts)):
                    if shifts[i] > 0.0: correlation_positive[i] = shifts[i]
                    else: correlation_negative[i] = shifts[i]
                ax.fill_between(xs, correlation_center, correlation_center + correlation_positive, color = 'g', alpha = 0.5)
                ax.fill_between(xs, correlation_center, correlation_center + correlation_negative, color = 'b', alpha = 0.5)

            observable = maps['observable'][np.array(data[index]['obs'])[idx][0]] % np.array(data[index]['target'])[idx][0]
            title = maps['experiment'][np.array(data[index]['col'])[idx][0]]
            ax.title.set_text(title + r'$,~$' + observable)
            # ax.text(0.07, 0.77, r'$E = %.1f,~Q^2 \in \left[%0.1f,~%0.1f \right]$' % (energys[0], min(Q2s), max(Q2s)), transform = ax.transAxes)
            ax.text(0.07, 0.77, r'$Q^2 \in \left[%0.1f,~%0.1f \right]$' % (min(Q2s), max(Q2s)), transform = ax.transAxes)

    x_label_counts = range(count, count - n_columns, -1)
    for x_label_count in x_label_counts:
        axs[x_label_count].set_xlabel(r'$x$', size = 15)

    py.tight_layout()
    py.savefig('%s/gallery/pidis-ape-deuteron-%d.png' % (wdir, istep), dpi = dpi)
    return

def data_Apa_deuteron_DVCS(wdir, data, indices, subset_bins, maps, istep, dpi):
    n_columns, n_rows = 4, 3
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.7))
    axs = {}
    count = 0
    for index in indices:
        correlation = [_ for _ in data[index] if ('_c' in _) and ('norm' not in _)]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if len(idx) <= 2: continue
            count += 1
            ax = py.subplot(n_rows, n_columns, count)
            axs[count] = ax

            xs_temp = data[index]['X'][idx]
            idx_order = np.argsort(xs_temp)
            xs = data[index]['X'][idx][idx_order]
            Q2s = data[index]['Q2'][idx][idx_order]
            experimental_values = data[index]['value'][idx][idx_order]

            alphas = data[index]['alpha'][idx][idx_order]
            shifts = data[index]['shift'][idx][idx_order]

            theory_mean = data[index]['thy'][idx][idx_order]
            theory_deviation = data[index]['dthy'][idx][idx_order]

            ax.errorbar(xs, experimental_values, yerr = alphas, fmt = 'k.', markersize = 5, capsize = 0, alpha = 0.8)
            ax.plot(xs, theory_mean, 'r-')
            ax.fill_between(xs, theory_mean - theory_deviation, theory_mean + theory_deviation, color = 'r', alpha = 0.5)

            if len(xs) > 6:
                ax.semilogx()

            if any(correlation):
                pass
                experimental_bottom = min(experimental_values - alphas)
                experimental_top = max(experimental_values - alphas)
                experimental_bottom = experimental_bottom - 0.2 * (experimental_top - experimental_bottom)
                correlation_center = np.full(len(alphas), experimental_bottom)
                correlation_positive, correlation_negative = np.zeros(len(alphas)), np.zeros(len(alphas))
                for i in range(len(shifts)):
                    if shifts[i] > 0.0: correlation_positive[i] = shifts[i]
                    else: correlation_negative[i] = shifts[i]
                ax.fill_between(xs, correlation_center, correlation_center + correlation_positive, color = 'g', alpha = 0.5)
                ax.fill_between(xs, correlation_center, correlation_center + correlation_negative, color = 'b', alpha = 0.5)

            observable = maps['observable'][np.array(data[index]['obs'])[idx][0]] % np.array(data[index]['target'])[idx][0]
            title = maps['experiment'][np.array(data[index]['col'])[idx][0]]
            ax.title.set_text(title + r'$,~$' + observable)
            ax.text(0.07, 0.77, r'$Q^2 \in \left[%0.1f,~%0.1f \right]$' % (min(Q2s), max(Q2s)), transform = ax.transAxes)

    x_label_counts = range(count, count - n_columns, -1)
    for x_label_count in x_label_counts:
        axs[x_label_count].set_xlabel(r'$x$', size = 15)

    py.tight_layout()
    py.savefig('%s/gallery/pidis-apa-deuteron-DVCS-%d.png' % (wdir, istep), dpi = dpi)
    return

def data_Apa_deuteron_EG1b(wdir, data, indices, subset_bins, maps, istep, dpi):
    n_columns, n_rows = 4, 3
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.7))
    axs = {}
    count = 0
    for index in indices:
        correlation = [_ for _ in data[index] if ('_c' in _) and ('norm' not in _)]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if len(idx) <= 2: continue
            count += 1
            ax = py.subplot(n_rows, n_columns, count)
            axs[count] = ax

            xs_temp = data[index]['X'][idx]
            idx_order = np.argsort(xs_temp)
            xs = data[index]['X'][idx][idx_order]
            Q2s = data[index]['Q2'][idx][idx_order]
            energys = data[index]['E'][idx][idx_order]
            experimental_values = data[index]['value'][idx][idx_order]

            alphas = data[index]['alpha'][idx][idx_order]
            shifts = data[index]['shift'][idx][idx_order]

            theory_mean = data[index]['thy'][idx][idx_order]
            theory_deviation = data[index]['dthy'][idx][idx_order]

            ax.errorbar(xs, experimental_values, yerr = alphas, fmt = 'k.', markersize = 5, capsize = 0, alpha = 0.8)
            ax.plot(xs, theory_mean, 'r-')
            ax.fill_between(xs, theory_mean - theory_deviation, theory_mean + theory_deviation, color = 'r', alpha = 0.5)

            if len(xs) > 6:
                ax.semilogx()

            if any(correlation):
                pass
                experimental_bottom = min(experimental_values - alphas)
                experimental_top = max(experimental_values - alphas)
                experimental_bottom = experimental_bottom - 0.2 * (experimental_top - experimental_bottom)
                correlation_center = np.full(len(alphas), experimental_bottom)
                correlation_positive, correlation_negative = np.zeros(len(alphas)), np.zeros(len(alphas))
                for i in range(len(shifts)):
                    if shifts[i] > 0.0: correlation_positive[i] = shifts[i]
                    else: correlation_negative[i] = shifts[i]
                ax.fill_between(xs, correlation_center, correlation_center + correlation_positive, color = 'g', alpha = 0.5)
                ax.fill_between(xs, correlation_center, correlation_center + correlation_negative, color = 'b', alpha = 0.5)

            observable = maps['observable'][np.array(data[index]['obs'])[idx][0]] % np.array(data[index]['target'])[idx][0]
            title = maps['experiment'][np.array(data[index]['col'])[idx][0]]
            ax.title.set_text(title + r'$,~$' + observable)
            ax.text(0.07, 0.77, r'$E = %.1f,~Q^2 \in \left[%0.1f,~%0.1f \right]$' % (energys[0], min(Q2s), max(Q2s)), transform = ax.transAxes)

    x_label_counts = range(count, count - n_columns, -1)
    for x_label_count in x_label_counts:
        axs[x_label_count].set_xlabel(r'$x$', size = 15)

    py.tight_layout()
    py.savefig('%s/gallery/pidis-apa-deuteron-EG1b-%d.png' % (wdir, istep), dpi = dpi)
    return

def data_helium(wdir, data, indices, subset_bins, maps, istep, dpi):
    n_columns, n_rows = 4, 3
    py.figure(figsize = (n_columns * 2.5, n_rows * 1.7))
    axs = {}
    count = 0
    for index in indices:
        correlation = [_ for _ in data[index] if ('_c' in _) and ('norm' not in _)]
        i_sets = np.sort(subset_bins[index].keys())
        for i_set in i_sets:
            idx = subset_bins[index][i_set]
            if len(idx) <= 2: continue
            count += 1
            ax = py.subplot(n_rows, n_columns, count)
            axs[count] = ax

            xs_temp = data[index]['X'][idx]
            idx_order = np.argsort(xs_temp)
            xs = data[index]['X'][idx][idx_order]
            Q2s = data[index]['Q2'][idx][idx_order]
            experimental_values = data[index]['value'][idx][idx_order]

            alphas = data[index]['alpha'][idx][idx_order]
            shifts = data[index]['shift'][idx][idx_order]

            theory_mean = data[index]['thy'][idx][idx_order]
            theory_deviation = data[index]['dthy'][idx][idx_order]

            ax.errorbar(xs, experimental_values, yerr = alphas, fmt = 'k.', markersize = 5, capsize = 0, alpha = 0.8)
            ax.plot(xs, theory_mean, 'r-')
            ax.fill_between(xs, theory_mean - theory_deviation, theory_mean + theory_deviation, color = 'r', alpha = 0.5)

            if len(xs) > 6:
                ax.semilogx()

            if any(correlation):
                pass
                experimental_bottom = min(experimental_values - alphas)
                experimental_top = max(experimental_values - alphas)
                experimental_bottom = experimental_bottom - 0.2 * (experimental_top - experimental_bottom)
                correlation_center = np.full(len(alphas), experimental_bottom)
                correlation_positive, correlation_negative = np.zeros(len(alphas)), np.zeros(len(alphas))
                for i in range(len(shifts)):
                    if shifts[i] > 0.0: correlation_positive[i] = shifts[i]
                    else: correlation_negative[i] = shifts[i]
                ax.fill_between(xs, correlation_center, correlation_center + correlation_positive, color = 'g', alpha = 0.5)
                ax.fill_between(xs, correlation_center, correlation_center + correlation_negative, color = 'b', alpha = 0.5)

            observable = maps['observable'][np.array(data[index]['obs'])[idx][0]] % np.array(data[index]['target'])[idx][0]
            title = maps['experiment'][np.array(data[index]['col'])[idx][0]]
            ax.title.set_text(title + r'$,~$' + observable)
            ax.text(0.07, 0.77, r'$Q^2 \in \left[%0.1f,~%0.1f \right]$' % (min(Q2s), max(Q2s)), transform = ax.transAxes)

    x_label_counts = range(count, count - n_columns, -1)
    for x_label_count in x_label_counts:
        axs[x_label_count].set_xlabel(r'$x$', size = 15)

    py.tight_layout()
    py.savefig('%s/gallery/pidis-helium-%d.png' % (wdir, istep), dpi = dpi)
    return

def make_figure(wdir, task, only_best_cluster = True, dpi = 200):

    if task == 1:
        print('\nplotting PIDIS data from %s' % (wdir))
    elif task == 2:
        print('\nplotting PIDIS data over theory from %s' % (wdir))

    load_config('%s/input.py' % wdir)
    istep = core.get_istep()
    replicas = core.get_replicas(wdir)
    core.mod_conf(istep, replicas[0]) #--set conf as specified in istep

    predictions = load('%s/data/predictions-%d.dat' % (wdir, istep))
    if 'pidis' not in predictions['reactions']:
        print('PIDIS is not in data file')
        return
    labels  = load('%s/data/labels-%d.dat' % (wdir, istep))
    cluster = labels['cluster']
    cluster_average_chi2 = labels['cluster_average_chi2']
    best_cluster = np.argmin(cluster_average_chi2)

    data = predictions['reactions']['pidis']

    for idx in data:
        predictions = copy.copy(data[idx]['prediction-rep'])
        del data[idx]['prediction-rep']
        del data[idx]['residuals-rep']
        if only_best_cluster:
            best_predictions = [predictions[i] for i in range(len(predictions)) if cluster[i] == best_cluster]
            data[idx]['thy'] = np.mean(best_predictions, axis = 0)
            data[idx]['dthy'] = np.std(best_predictions, axis = 0)
        else:
            all_predictions = [predictions[i] for i in range(len(predictions))]
            data[idx]['thy'] = np.mean(all_predictions, axis = 0)
            data[idx]['dthy'] = np.std(all_predictions, axis = 0)

    maps = generate_maps()
    subset_bins = pick_subsets(data)
    # log_list = []

    if task == 1:
        indices = {}
        for idx in data:
            observable = data[idx]['obs'][0]
            target = data[idx]['target'][0]
            collaboration = data[idx]['col'][0]
            if (observable == 'A1') or (observable == 'Apa'):
                if target == 'p':
                    if collaboration.startswith('JLab') != True:
                        if 'Apa_proton' not in indices:
                            indices['Apa_proton'] = []
                        indices['Apa_proton'].append(idx)
                    elif collaboration.startswith('JLabHB(EG1DVCS)') == True:
                        if 'Apa_proton_DVCS' not in indices:
                            indices['Apa_proton_DVCS'] = []
                        indices['Apa_proton_DVCS'].append(idx)
                    elif collaboration.startswith('JLabHB(EG1b)') == True:
                        if 'Apa_proton_EG1b' not in indices:
                            indices['Apa_proton_EG1b'] = []
                        indices['Apa_proton_EG1b'].append(idx)
                elif target == 'd':
                    if collaboration.startswith('JLab') != True:
                        if 'Apa_deuteron' not in indices:
                            indices['Apa_deuteron'] = []
                        indices['Apa_deuteron'].append(idx)
                    elif collaboration.startswith('JLabHB(EG1DVCS)') == True:
                        if 'Apa_deuteron_DVCS' not in indices:
                            indices['Apa_deuteron_DVCS'] = []
                        indices['Apa_deuteron_DVCS'].append(idx)
                    elif collaboration.startswith('JLabHB(EG1b)') == True:
                        if 'Apa_deuteron_EG1b' not in indices:
                            indices['Apa_deuteron_EG1b'] = []
                        indices['Apa_deuteron_EG1b'].append(idx)
            elif (observable == 'A2') or (observable == 'Ape') or (observable == 'Atpe'):
                if target == 'p':
                    if collaboration.startswith('JLab') != True:
                        if 'Ape_proton' not in indices:
                            indices['Ape_proton'] = []
                        indices['Ape_proton'].append(idx)
                elif target == 'd':
                    if collaboration.startswith('JLab') != True:
                        if 'Ape_deuteron' not in indices:
                            indices['Ape_deuteron'] = []
                        indices['Ape_deuteron'].append(idx)

            if target == 'h':
                if 'helium' not in indices:
                    indices['helium'] = []
                indices['helium'].append(idx)

        for key in indices:
            print(key)
            if key == 'Apa_proton':
                data_Apa_proton(wdir, data, indices[key], subset_bins, maps, istep, dpi)
            elif key == 'Apa_proton_DVCS':
                data_Apa_proton_DVCS(wdir, data, indices[key], subset_bins, maps, istep, dpi)
            elif key == 'Apa_proton_EG1b':
                data_Apa_proton_EG1b(wdir, data, indices[key], subset_bins, maps, istep, dpi)
            elif key == 'Apa_deuteron':
                data_Apa_deuteron(wdir, data, indices[key], subset_bins, maps, istep, dpi)
            elif key == 'Apa_deuteron_DVCS':
                data_Apa_deuteron_DVCS(wdir, data, indices[key], subset_bins, maps, istep, dpi)
            elif key == 'Apa_deuteron_EG1b':
                data_Apa_deuteron_EG1b(wdir, data, indices[key], subset_bins, maps, istep, dpi)
            elif key == 'Ape_proton':
                data_Ape_proton(wdir, data, indices[key], subset_bins, maps, istep, dpi)
            elif key == 'Ape_deuteron':
                data_Ape_deuteron(wdir, data, indices[key], subset_bins, maps, istep, dpi)
            elif key == 'helium':
                data_helium(wdir, data, indices[key], subset_bins, maps, istep, dpi)

    elif task == 2:
        plot_data_on_theory(wdir, data, istep, dpi)

    return

if __name__ == "__main__":
    pass
