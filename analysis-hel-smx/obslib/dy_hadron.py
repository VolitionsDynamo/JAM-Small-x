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

def plot_ratio(wdir, data, istep, Q2_bins, dpi = 200):
    ratio_1 = 2.5
    ratio_2 = -3.0
    n_rows, n_columns = 4, 4
    fig = py.figure(figsize = (n_columns * 4.5, n_rows * 3.0))
    count = 0
    ax = {}
    label_1 = r'$\mathrm{E866}~pp$'
    label_2 = r'$\mathrm{E866}~pd$'

    for _ in data:
        for i in range(len(Q2_bins)):
            Q2_min, Q2_max = Q2_bins[i]
            thyk = data[_]['thy']
            ratio= data[_]['value'] / thyk
            data[_].update({'ratio' : ratio})
            df = pd.DataFrame(data[_])
            df = df.query('Q2>%f and Q2<%f' % (Q2_min, Q2_max))
            count += 1
            ax[count] = py.subplot(n_rows, n_columns, count)
            df = df.query('ratio<%f' % (ratio_1))
            df = df.query('ratio>%f' % (ratio_2))
            thy = df['thy']
            xF = df.xF
            value = df.value
            alpha = df.alpha
            if _ == 10001:
                ax[count].errorbar(xF, value/thy, alpha/thy, color = 'firebrick', marker = '.', linestyle = 'none', label = label_1)
            if _ == 10002:
                ax[count].errorbar(xF, value/thy, alpha/thy, color = 'darkgreen', marker = '^', linestyle = 'none', label = label_2)
    count_max = count

    ax[1].legend(fontsize = 20, loc = 9)
    ax[9].legend(fontsize = 20, loc = 'upper left')

    ## proton
    ax[1].set_ylim(0.5, 1.6)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[1].yaxis.set_minor_locator(minorLocator)
    ax[1].yaxis.set_major_locator(majorLocator)
    ax[1].text(0.1, 0.6, r'$Q^2~\in~[37,42]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[2].set_ylim(0.5, 1.3)
    minorLocator = MultipleLocator(0.05)
    majorLocator = MultipleLocator(0.25)
    ax[2].yaxis.set_minor_locator(minorLocator)
    ax[2].yaxis.set_major_locator(majorLocator)
    ax[2].text(0.07, 0.57, r'$Q^2~\in~[47,49]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[3].set_ylim(0.5, 2.1)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[3].yaxis.set_minor_locator(minorLocator)
    ax[3].yaxis.set_major_locator(majorLocator)
    ax[3].text(0.1, 1.8, r'$Q^2~\in~[54,58]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[4].set_ylim(0.7, 2.1)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[4].yaxis.set_minor_locator(minorLocator)
    ax[4].yaxis.set_major_locator(majorLocator)
    ax[4].text(0.1, 1.8, r'$Q^2~\in~[63,66]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[5].set_ylim(0.3, 2.5)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[5].yaxis.set_minor_locator(minorLocator)
    ax[5].yaxis.set_major_locator(majorLocator)
    ax[5].text(0.1, 2.1, r'$Q^2~\in~[70,73]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[6].set_ylim(0.5, 2.5)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[6].yaxis.set_minor_locator(minorLocator)
    ax[6].yaxis.set_major_locator(majorLocator)
    ax[6].text(0.15, 2.1, r'$Q^2~\in~[118,131]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[7].set_ylim(0.3, 2.5)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[7].yaxis.set_minor_locator(minorLocator)
    ax[7].yaxis.set_major_locator(majorLocator)
    ax[7].text(0.1, 2.1, r'$Q^2~\in~[148,156]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[8].set_ylim(0.4,  2.7)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[8].yaxis.set_minor_locator(minorLocator)
    ax[8].yaxis.set_major_locator(majorLocator)
    ax[8].text(0.32, 2,r' $Q^2~\in~[173,222]\,\mathrm{GeV^2}$', fontsize = 22)

    ## deuteron
    ax[9].set_ylim(0.5, 1.7)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[9].yaxis.set_minor_locator(minorLocator)
    ax[9].yaxis.set_major_locator(majorLocator)
    ax[9].text(0.1, 0.6, r'$Q^2~\in~[37,42]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[10].set_ylim(0.5, 1.6)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[10].yaxis.set_minor_locator(minorLocator)
    ax[10].yaxis.set_major_locator(majorLocator)
    ax[10].text(0.2, 1.4, r'$Q^2~\in~[47,49]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[11].set_ylim(0.4, 1.72)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[11].yaxis.set_minor_locator(minorLocator)
    ax[11].yaxis.set_major_locator(majorLocator)
    ax[11].text(0.2, 1.5, r'$Q^2~\in~[52,57]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[12].set_ylim(0.6, 1.42)
    minorLocator = MultipleLocator(0.05)
    majorLocator = MultipleLocator(0.25)
    ax[12].yaxis.set_minor_locator(minorLocator)
    ax[12].yaxis.set_major_locator(majorLocator)
    ax[12].text(0.1, 1.27, r'$Q^2~\in~[62,65]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[13].set_ylim(0.45, 1.35)
    minorLocator = MultipleLocator(0.05)
    majorLocator = MultipleLocator(0.25)
    ax[13].yaxis.set_minor_locator(minorLocator)
    ax[13].yaxis.set_major_locator(majorLocator)
    ax[13].text(0.2, 1.2, r'$Q^2~\in~[70,72]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[14].set_ylim(0.45, 2.5)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[14].yaxis.set_minor_locator(minorLocator)
    ax[14].yaxis.set_major_locator(majorLocator)
    ax[14].text(0.2, 2.2, r'$Q^2~\in~[124,129]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[15].set_ylim(0.45, 1.75)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[15].yaxis.set_minor_locator(minorLocator)
    ax[15].yaxis.set_major_locator(majorLocator)
    ax[15].text(0.05, 1.55, r'$Q^2~\in~[145,154]\,\mathrm{GeV^2}$', fontsize = 22)

    ax[16].set_ylim(0.0, 2.0)
    minorLocator = MultipleLocator(0.1)
    majorLocator = MultipleLocator(0.5)
    ax[16].yaxis.set_minor_locator(minorLocator)
    ax[16].yaxis.set_major_locator(majorLocator)
    ax[16].text(0.05, 0.2, r'$Q^2~\in~[73,280]\,\mathrm{GeV^2}$', fontsize = 22)

    for count in range(1, count_max + 1):
        if count % 4 == 1:
            ax[count].set_ylabel(r'\boldmath$\mathrm{data/theory}$', size = 24)
            ax[count].yaxis.set_label_coords(-0.17, 0.5)

    for count in range(count_max, count_max - n_columns, -1):
        ax[count].set_xlabel(r'\boldmath$x_{\rm F}$', size = 24)

    for i in range(1, 13):
        ax[i].set_xticklabels([])

    for _ in ax:
        ax[_].axhline(1, color = 'black', ls = ':')
        #AX[_].set_ylim(0.8,1.2)
        ax[_].set_xlim(0.0, 1.0)
        ax[_].xaxis.set_label_coords(0.90, -0.12)
        ax[_].tick_params(axis = 'x', which = 'major', labelsize = 18)
        ax[_].tick_params(axis = 'y', which = 'major', labelsize = 18)
        ax[_].xaxis.set_tick_params(which = 'major', length = 6)
        ax[_].xaxis.set_tick_params(which = 'minor', length = 3)
        ax[_].yaxis.set_tick_params(which = 'major', length = 6)
        ax[_].yaxis.set_tick_params(which = 'minor', length = 3)
        minorLocator = MultipleLocator(0.05)
        majorLocator = MultipleLocator(0.2)
        ax[_].xaxis.set_minor_locator(minorLocator)
        ax[_].xaxis.set_major_locator(majorLocator)

    py.tight_layout()
    py.subplots_adjust(left = 0.08, bottom = 0.08, right = 0.99, top = 0.97, wspace = 0.2, hspace = 0.1)
    py.savefig('%s/gallery/dy-%d.png' % (wdir, istep), dpi = dpi)
    py.close()

def plot_obs(wdir, only_best_cluster = True, dpi = 200):

    print('\nplotting dy data from %s' % (wdir))

    load_config('%s/input.py' % wdir)
    istep = core.get_istep()
    replicas = core.get_replicas(wdir)
    core.mod_conf(istep, replicas[0]) #--set conf as specified in istep

    predictions = load('%s/data/predictions-%d.dat' % (wdir, istep))
    if 'dy' not in predictions['reactions']:
        print('DY is not in data file')
        return
    labels  = load('%s/data/labels-%d.dat' % (wdir, istep))
    cluster = labels['cluster']
    cluster_average_chi2 = labels['cluster_average_chi2']
    best_cluster = np.argmin(cluster_average_chi2)

    data = predictions['reactions']['dy']

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

    Q2_bins = []
    # Q2_bins.append([30, 40])
    Q2_bins.append([35, 45])    ## label pp:  37 < Q2 < 42   pd: idem
    Q2_bins.append([45, 52])    ## label pp:  47 < Q2 < 49   pd: idem
    Q2_bins.append([52, 60])    ## label pp:  54 < Q2 < 58   pd:  55 < Q2 < 57
    Q2_bins.append([60, 68])    ## label pp:  63 < Q2 < 66   pd:  62 < Q2 < 65
    Q2_bins.append([68, 75])    ## label pp:  70 < Q2 < 73   pd:  70 < Q2 < 72
    Q2_bins.append([100, 140])  ## label pp: 118 < Q2 < 131  pd: 124 < Q2 < 129
    Q2_bins.append([140, 160])  ## label pp: 148 < Q2 < 156  pd: 145 < Q2 < 154
    Q2_bins.append([160, 280])  ## label pp: 173 < Q2 < 222  pd:  73 < Q2 < 280

    plot_ratio(wdir, data, istep, Q2_bins, dpi)

    return
