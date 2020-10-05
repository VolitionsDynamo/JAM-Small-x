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

def plot_proton(wdir, data, istep, dpi):
    datanmc = data[10020]  ## NMC p
    dfnmc = pd.DataFrame(datanmc)

    databcdms = data[10016]  ## BCDMS p
    dfbcdms = pd.DataFrame(databcdms)

    dataslac = data[10010]  ## SLAC p
    dfslac = pd.DataFrame(dataslac)

    ## NMC proton
    df_nmc_query = {}
    df_nmc_query[24] = dfnmc.query('X>4.8e-3 and X < 5.7e-3') ## x = 0.005 2^24
    df_nmc_query[23] = dfnmc.query('X>7.2e-3 and X < 9.3e-3') ## x = 0.008 2^23
    df_nmc_query[22] = dfnmc.query('X>1.15e-2 and X < 1.4e-2') ## x = 0.013 2^22
    df_nmc_query[21] = dfnmc.query('X>1.65e-2 and X < 1.9e-2') ## x = 0.018 2^21
    df_nmc_query[20] = dfnmc.query('X>2.3e-2 and X < 2.9e-2') ## x = 0.025 2^20
    df_nmc_query[19] = dfnmc.query('X>2.9e-2 and X < 3.8e-2') ## x = 0.035 2^19
    df_nmc_query[18] = dfnmc.query('X>4.6e-2 and X < 5.4e-2') ## x = 0.05 2^18
    df_nmc_query[17] = dfnmc.query('X>6.5e-2 and X < 7.3e-2') ## x = 0.07 2^17
    df_nmc_query[15] = dfnmc.query('X>8.5e-2 and X < 9.2e-2') ## x = 0.09 2^15
    df_nmc_query[14] = dfnmc.query('X>9.8e-2 and X < 10.3e-2') ## x = 0.10 2^14
    df_nmc_query[13] = dfnmc.query('X>10.7e-2 and X < 11.3e-2') ## x = 0.11 2^13
    df_nmc_query[12] = dfnmc.query('X>13.6e-2 and X < 14.6e-2') ## x = 0.14 2^12
    df_nmc_query[11] = dfnmc.query('X>17.1e-2 and X < 18.7e-2') ## x = 0.18 2^11
    df_nmc_query[10] = dfnmc.query('X>21.7e-2 and X < 23.7e-2') ## x = 0.23 2^10
    df_nmc_query[9] = dfnmc.query('X>26.0e-2 and X < 29.0e-2') ## x = 0.28 2^9
    df_nmc_query[8] = dfnmc.query('X>33.0e-2 and X < 36.0e-2') ## x = 0.35 2^8
    df_nmc_query[7] = dfnmc.query('X>42.0e-2 and X < 47.6e-2') ## x = 0.45 2^7
    df_nmc_query[6] = dfnmc.query('X>47.6e-2 and X < 52.5e-2') ## x = 0.50 2^6

    df_nmc_plot = {}
    df_nmc_plot['theory'] = {}
    df_nmc_plot['Q2'] = {}
    df_nmc_plot['value'] = {}
    df_nmc_plot['alpha'] = {}
    for key in df_nmc_query:
        df_nmc_plot['theory'][key] = df_nmc_query[key]['thy']
        df_nmc_plot['Q2'][key] = df_nmc_query[key]['Q2']
        df_nmc_plot['value'][key] = df_nmc_query[key]['value']
        df_nmc_plot['alpha'][key] = df_nmc_query[key]['alpha']

    ## SLAC proton
    df_slac_query = {}
    df_slac_query[16] = dfslac.query('X>7.6e-2 and X < 8.4e-2') ## x = 0.08 2^16
    df_slac_query[14] = dfslac.query('X>9.6e-2 and X < 10.7e-2') ## x = 0.10 2^14
    df_slac_query[12] = dfslac.query('X>13.5e-2 and X < 14.6e-2') ## x = 0.14 2^12
    df_slac_query[11] = dfslac.query('X>17.1e-2 and X < 18.7e-2') ## x = 0.18 2^11
    df_slac_query[10] = dfslac.query('X>21.5e-2 and X < 23.8e-2') ## x = 0.23 2^10
    df_slac_query[9] = dfslac.query('X>26.0e-2 and X < 29.0e-2') ## x = 0.28 2^9
    df_slac_query[8] = dfslac.query('X>33.0e-2 and X < 36.0e-2') ## x = 0.35 2^8
    df_slac_query[7] = dfslac.query('X>42.0e-2 and X < 47.6e-2') ## x = 0.45 2^7
    df_slac_query[6] = dfslac.query('X>47.6e-2 and X < 52.5e-2') ## x = 0.50 2^6
    df_slac_query[5] = dfslac.query('X>53.0e-2 and X < 57.7e-2') ## x = 0.55 2^5
    df_slac_query[4] = dfslac.query('X>57.7e-2 and X < 61.5e-2') ## x = 0.58 2^4
    df_slac_query[3] = dfslac.query('X>61.5e-2 and X < 62.5e-2') ## x = 0.61 2^3
    df_slac_query[2] = dfslac.query('X>62.5e-2 and X < 67.0e-2') ## x = 0.65 2^2
    df_slac_query[1] = dfslac.query('X>68.0e-2 and X < 73.0e-2') ## x = 0.70 2

    df_slac_plot = {}
    df_slac_plot['theory'] = {}
    df_slac_plot['Q2'] = {}
    df_slac_plot['value'] = {}
    df_slac_plot['alpha'] = {}
    for key in df_slac_query:
        df_slac_plot['theory'][key] = df_slac_query[key]['thy']
        df_slac_plot['Q2'][key] = df_slac_query[key]['Q2']
        df_slac_plot['value'][key] = df_slac_query[key]['value']
        df_slac_plot['alpha'][key] = df_slac_query[key]['alpha']

    ## BCDMS proton
    df_bcdms_query = {}
    df_bcdms_query[17] = dfbcdms.query('X>6.8e-2 and X < 7.2e-2') ## x = 0.07 2^17
    df_bcdms_query[14] = dfbcdms.query('X>9.6e-2 and X < 10.7e-2') ## x = 0.10 2^14
    df_bcdms_query[12] = dfbcdms.query('X>13.5e-2 and X < 14.6e-2') ## x = 0.14 2^12
    df_bcdms_query[11] = dfbcdms.query('X>17.1e-2 and X < 18.7e-2') ## x = 0.18 2^11
    df_bcdms_query[10] = dfbcdms.query('X>21.7e-2 and X < 23.7e-2') ## x = 0.23 2^10
    df_bcdms_query[9] = dfbcdms.query('X>26.0e-2 and X < 29.0e-2') ## x = 0.28 2^9
    df_bcdms_query[8] = dfbcdms.query('X>33.0e-2 and X < 36.0e-2') ## x = 0.35 2^8
    df_bcdms_query[7] = dfbcdms.query('X>42.0e-2 and X < 48.0e-2') ## x = 0.45 2^7
    df_bcdms_query[5] = dfbcdms.query('X>53.0e-2 and X < 58.0e-2') ## x = 0.55 2^5
    df_bcdms_query[2] = dfbcdms.query('X>62.0e-2 and X < 66.0e-2') ## x = 0.65 2^2
    df_bcdms_query[1] = dfbcdms.query('X>73.0e-2 and X < 76.0e-2') ## x = 0.75 2^0

    df_bcdms_plot = {}
    df_bcdms_plot['theory'] = {}
    df_bcdms_plot['Q2'] = {}
    df_bcdms_plot['value'] = {}
    df_bcdms_plot['alpha'] = {}
    for key in df_bcdms_query:
        df_bcdms_plot['theory'][key] = df_bcdms_query[key]['thy']
        df_bcdms_plot['Q2'][key] = df_bcdms_query[key]['Q2']
        df_bcdms_plot['value'][key] = df_bcdms_query[key]['value']
        df_bcdms_plot['alpha'][key] = df_bcdms_query[key]['alpha']

    nrows, ncols = 1, 1
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax = py.subplot(nrows, ncols, 1)

    for key in df_nmc_plot['value']:
        NMC = ax.errorbar(df_nmc_plot['Q2'][key], 2.0 ** float(key) * df_nmc_plot['value'][key], df_nmc_plot['alpha'][key], \
                          marker = '*', color = 'firebrick', markersize = 6, linestyle = 'none')
        ax.plot(df_nmc_plot['Q2'][key], 2.0 ** float(key) * df_nmc_plot['theory'][key], linestyle = 'solid', color = 'firebrick')

    for key in df_slac_plot['value']:
        SLAC = ax.errorbar(df_slac_plot['Q2'][key], 2.0 ** float(key) * df_slac_plot['value'][key], df_slac_plot['alpha'][key], \
                           marker = '^', color = 'darkgreen', markersize = 6, linestyle = 'none')
        ax.plot(df_slac_plot['Q2'][key], 2.0 ** float(key) * df_slac_plot['theory'][key], linestyle = 'solid', color = 'darkgreen')

    for key in df_bcdms_plot['value']:
        BCDMS = ax.errorbar(df_bcdms_plot['Q2'][key], 2.0 ** float(key) * df_bcdms_plot['value'][key], df_bcdms_plot['alpha'][key], \
                            marker = 'o', color = 'royalblue', markersize = 6, linestyle = 'none')
        ax.plot(df_bcdms_plot['Q2'][key], 2.0 ** float(key) * df_bcdms_plot['theory'][key], linestyle = 'solid', color = 'royalblue')

    ax.text(2.4, 6e6, r'$x=0.005\, (i=24)$', fontsize = 18)
    ax.text(3.7, 3.3e6, r'$x=0.008$', fontsize = 18)
    ax.text(5.8, 1.8e6, r'$x=0.0013$', fontsize = 18)
    ax.text(7.5, 9e5, r'$x=0.018$', fontsize = 18)
    ax.text(12.0, 4.5e5, r'$x=0.035$', fontsize = 18)
    ax.text(15.0, 2.1e5, r'$x=0.035$', fontsize = 18)
    ax.text(21.0, 1.05e5, r'$x=0.05$', fontsize = 18)
    ax.text(28.0, 49.5e3, r'$x=0.07$', fontsize = 18)
    ax.text(2.7, 22e3, r'$x=0.08$', fontsize = 18)
    ax.text(28.0, 1.2e4, r'$x=0.09$', fontsize = 18)
    ax.text(41.0, 6200.0, r'$x=0.10$', fontsize = 18)
    ax.text(38.0, 3000.0, r'$x=0.11$', fontsize = 18)
    ax.text(62.0, 1400.0, r'$x=0.14$', fontsize = 18)
    ax.text(70.0, 640.0, r'$x=0.18$', fontsize = 18)
    ax.text(93.0, 290.0, r'$x=0.23$', fontsize = 18)
    ax.text(125.0, 119.0, r'$x=0.28$', fontsize = 18)
    ax.text(150.0, 44.0, r'$x=0.35$', fontsize = 18)
    ax.text(188.0, 13.0, r'$x=0.45$', fontsize = 18)
    ax.text(70.0, 6.0, r'$x=0.50$', fontsize = 18)
    ax.text(250.0, 2.0, r'$x=0.55$', fontsize = 18)
    ax.text(23.0, 1.0, r'$x=0.58$', fontsize = 18)
    ax.text(23.0, 0.4, r'$x=0.61$', fontsize = 18)
    ax.text(250.0, 0.1, r'$x=0.65$', fontsize = 18)
    ax.text(28.0, 0.038, r'$x=0.70$', fontsize = 18)
    ax.text(250.0, 0.008, r'$x=0.75$', fontsize = 18)
    ax.text(17.0, 0.012, r'$(i=0)$', fontsize = 18)

    ax.semilogy()
    ax.semilogx()
    ax.set_xlim(1.5, 5e2)
    ax.set_ylim(0.004, 1.5e7)

    ax.tick_params(axis = 'both', labelsize = 18)

    ax.yaxis.set_tick_params(which = 'major', length = 6)
    ax.set_yticks([0.01, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
    ax.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$', r'$10^7$'])
    locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.yaxis.set_tick_params(which = 'minor', length = 3)

    ax.xaxis.set_tick_params(which = 'major', length = 6)
    ax.xaxis.set_tick_params(which = 'minor', length = 3)
    ax.set_xticks([2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2])
    xtick_labels = [2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2]
    ax.set_xticklabels(['$%0.0f$' % x for x in xtick_labels])

    ax.set_ylabel(r'$F_2^p$'+r'$\, \times\, 2^{\, i}$', size = 24)
    ax.set_xlabel(r'$Q^2 \: \rm{(GeV^2)}$', size = 24)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = 'on', labelright = 'off')
    ax.tick_params(axis = 'x', which = 'both', labeltop = 'off', labelbottom = 'on')

    ax.legend([NMC, SLAC, BCDMS], [r'$\mathrm{NMC}$', r'$\mathrm{SLAC}$', r'$\mathrm{BCDMS}$'], fontsize = 24, frameon = 0)
    py.tight_layout()
    py.savefig('%s/gallery/dis-proton-%d.png' % (wdir, istep), dpi = dpi)
    py.close()

def plot_deuteron(wdir, data, istep, dpi):
    datanmc = data[10021] ## NMC d/p
    dfnmc = pd.DataFrame(datanmc)

    databcdms = data[10017] ## BCDMS d
    dfbcdms = pd.DataFrame(databcdms)

    dataslac = data[10011] ## SLAC d
    dfslac = pd.DataFrame(dataslac)

    ## NMC proton
    df_nmc_query = {}
    df_nmc_query[17] = dfnmc.query('X>4.8e-3 and X < 5.7e-3') ## x = 0.005 2^17
    df_nmc_query[16] = dfnmc.query('X>7.2e-3 and X < 9.3e-3') ## x = 0.008 2^16
    df_nmc_query[15] = dfnmc.query('X>1.15e-2 and X < 1.4e-2') ## x = 0.013 2^15
    df_nmc_query[14] = dfnmc.query('X>1.65e-2 and X < 1.9e-2') ## x = 0.018 2^14
    df_nmc_query[13] = dfnmc.query('X>2.3e-2 and X < 2.9e-2') ## x = 0.025 2^13
    df_nmc_query[12] = dfnmc.query('X>2.9e-2 and X < 3.8e-2') ## x   =0.035 2^12
    df_nmc_query[11] = dfnmc.query('X>4.6e-2 and X < 5.4e-2') ## x = 0.05 2^11
    df_nmc_query[10] = dfnmc.query('X>6.5e-2 and X < 7.3e-2') ## x = 0.07 2^10
    df_nmc_query[9] = dfnmc.query('X>8.5e-2 and X < 9.2e-2') ## x = 0.09 2^9
    df_nmc_query[8] = dfnmc.query('X>10.7e-2 and X < 11.3e-2') ## x = 0.11 2^8
    df_nmc_query[7] = dfnmc.query('X>13.6e-2 and X < 14.6e-2') ## x = 0.14 2^7
    df_nmc_query[6] = dfnmc.query('X>17.1e-2 and X < 18.7e-2') ## x = 0.18 2^6
    df_nmc_query[5] = dfnmc.query('X>21.7e-2 and X < 23.7e-2') ## x = 0.23 2^5
    df_nmc_query[4] = dfnmc.query('X>26.0e-2 and X < 29.0e-2') ## x = 0.28 2^4
    df_nmc_query[3] = dfnmc.query('X>33.0e-2 and X < 36.0e-2') ## x = 0.35 2^3
    df_nmc_query[2] = dfnmc.query('X>42.0e-2 and X < 47.6e-2') ## x = 0.45 2^2
    df_nmc_query[1] = dfnmc.query('X>47.6e-2 and X < 57.5e-2') ## x = 0.55 2^1
    df_nmc_query[0] = dfnmc.query('X>57.5e-2 and X < 70.5e-2') ## x = 0.68 2^0

    df_nmc_plot = {}
    df_nmc_plot['theory'] = {}
    df_nmc_plot['Q2'] = {}
    df_nmc_plot['value'] = {}
    df_nmc_plot['alpha'] = {}
    for key in df_nmc_query:
        df_nmc_plot['theory'][key] = df_nmc_query[key]['thy']
        df_nmc_plot['Q2'][key] = df_nmc_query[key]['Q2']
        df_nmc_plot['value'][key] = df_nmc_query[key]['value']
        df_nmc_plot['alpha'][key] = df_nmc_query[key]['alpha']

    ## SLAC proton
    df_slac_query = {}
    df_slac_query[14] = dfslac.query('X>7.6e-2 and X < 8.4e-2') ## x = 0.08 2^14
    df_slac_query[13] = dfslac.query('X>9.6e-2 and X < 10.7e-2') ## x = 0.10 2^13
    df_slac_query[12] = dfslac.query('X>13.5e-2 and X < 14.6e-2') ## x = 0.14 2^12
    df_slac_query[11] = dfslac.query('X>17.1e-2 and X < 18.7e-2') ## x = 0.18 2^11
    df_slac_query[10] = dfslac.query('X>21.5e-2 and X < 23.8e-2') ## x = 0.23 2^10
    df_slac_query[9] = dfslac.query('X>26.0e-2 and X < 29.0e-2') ## x = 0.28 2^9
    df_slac_query[8] = dfslac.query('X>33.0e-2 and X < 36.0e-2') ## x = 0.35 2^8
    df_slac_query[7] = dfslac.query('X>42.0e-2 and X < 47.6e-2') ## x = 0.45 2^7
    df_slac_query[6] = dfslac.query('X>47.6e-2 and X < 52.5e-2') ## x = 0.50 2^6
    df_slac_query[5] = dfslac.query('X>53.0e-2 and X < 57.7e-2') ## x = 0.55 2^5
    df_slac_query[4] = dfslac.query('X>57.7e-2 and X < 61.5e-2') ## x = 0.58 2^4
    df_slac_query[3] = dfslac.query('X>61.5e-2 and X < 62.5e-2') ## x = 0.61 2^3
    df_slac_query[2] = dfslac.query('X>62.5e-2 and X < 67.0e-2') ## x = 0.65 2^2
    df_slac_query[1] = dfslac.query('X>68.0e-2 and X < 73.0e-2') ## x = 0.70 2

    df_slac_plot = {}
    df_slac_plot['theory'] = {}
    df_slac_plot['Q2'] = {}
    df_slac_plot['value'] = {}
    df_slac_plot['alpha'] = {}
    for key in df_slac_query:
        df_slac_plot['theory'][key] = df_slac_query[key]['thy']
        df_slac_plot['Q2'][key] = df_slac_query[key]['Q2']
        df_slac_plot['value'][key] = df_slac_query[key]['value']
        df_slac_plot['alpha'][key] = df_slac_query[key]['alpha']

    ## BCDMS proton
    df_bcdms_query = {}
    df_bcdms_query[15] = dfbcdms.query('X>6.8e-2 and X < 7.2e-2') ## x = 0.07 2^15
    df_bcdms_query[13] = dfbcdms.query('X>9.6e-2 and X < 10.7e-2') ## x = 0.10 2^13
    df_bcdms_query[12] = dfbcdms.query('X>13.5e-2 and X < 14.6e-2') ## x = 0.14 2^12
    df_bcdms_query[11] = dfbcdms.query('X>17.1e-2 and X < 18.7e-2') ## x = 0.18 2^11
    df_bcdms_query[10] = dfbcdms.query('X>21.7e-2 and X < 23.7e-2') ## x = 0.23 2^10
    df_bcdms_query[9] = dfbcdms.query('X>26.0e-2 and X < 29.0e-2') ## x = 0.28 2^9
    df_bcdms_query[8] = dfbcdms.query('X>33.0e-2 and X < 36.0e-2') ## x = 0.35 2^8
    df_bcdms_query[7] = dfbcdms.query('X>42.0e-2 and X < 48.0e-2') ## x = 0.45 2^7
    df_bcdms_query[5] = dfbcdms.query('X>53.0e-2 and X < 58.0e-2') ## x = 0.55 2^5
    df_bcdms_query[2] = dfbcdms.query('X>62.0e-2 and X < 66.0e-2') ## x = 0.65 2^2
    df_bcdms_query[0] = dfbcdms.query('X>73.0e-2 and X < 76.0e-2') ## x = 0.75 2^0

    df_bcdms_plot = {}
    df_bcdms_plot['theory'] = {}
    df_bcdms_plot['Q2'] = {}
    df_bcdms_plot['value'] = {}
    df_bcdms_plot['alpha'] = {}
    for key in df_bcdms_query:
        df_bcdms_plot['theory'][key] = df_bcdms_query[key]['thy']
        df_bcdms_plot['Q2'][key] = df_bcdms_query[key]['Q2']
        df_bcdms_plot['value'][key] = df_bcdms_query[key]['value']
        df_bcdms_plot['alpha'][key] = df_bcdms_query[key]['alpha']

    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 8.5, nrows * 10.0))
    ax = {}
    ax[1] = py.subplot(nrows, ncols, 1)
    ax[2] = py.subplot(nrows, ncols, 2)

    for key in df_nmc_plot['value']:
        NMC = ax[2].errorbar(df_nmc_plot['Q2'][key], 2.0 ** float(key) * df_nmc_plot['value'][key], df_nmc_plot['alpha'][key], \
                          marker = '*', color = 'firebrick', markersize = 6, linestyle = 'none')
        ax[2].plot(df_nmc_plot['Q2'][key], 2.0 ** float(key) * df_nmc_plot['theory'][key], linestyle = 'solid', color = 'firebrick')

    for key in df_slac_plot['value']:
        SLAC = ax[1].errorbar(df_slac_plot['Q2'][key], 2.0 ** float(key) * df_slac_plot['value'][key], df_slac_plot['alpha'][key], \
                           marker = '^', color = 'darkgreen', markersize = 6, linestyle = 'none')
        ax[1].plot(df_slac_plot['Q2'][key], 2.0 ** float(key) * df_slac_plot['theory'][key], linestyle = 'solid', color = 'darkgreen')

    for key in df_bcdms_plot['value']:
        BCDMS = ax[1].errorbar(df_bcdms_plot['Q2'][key], 2.0 ** float(key) * df_bcdms_plot['value'][key], df_bcdms_plot['alpha'][key], \
                            marker = 'o', color = 'royalblue', markersize = 6, linestyle = 'none')
        ax[1].plot(df_bcdms_plot['Q2'][key], 2.0 ** float(key) * df_bcdms_plot['theory'][key], linestyle = 'solid', color = 'royalblue')

    ax[1].text(1.8, 1.2e4, r'$x=0.07 \,(i=15)$', fontsize = 18)
    ax[1].text(3.0, 0.5e4, r'$x=0.08$', fontsize = 18)
    ax[1].text(45.0, 2780.0, r'$x=0.10$', fontsize = 18)
    ax[1].text(70.0, 1250.0, r'$x=0.14$', fontsize = 18)
    ax[1].text(80.0, 580.0, r'$x=0.18$', fontsize = 18)
    ax[1].text(100.0, 250.0, r'$x=0.23$', fontsize = 18)
    ax[1].text(140.0, 105.0, r'$x=0.28$', fontsize = 18)
    ax[1].text(170.0, 37.0, r'$x=0.35$', fontsize = 18)
    ax[1].text(2.9, 17.0, r'$x=0.45$', fontsize = 18)
    ax[1].text(3.5, 6.5, r'$x=0.50$', fontsize = 18)
    ax[1].text(4.5, 2.2, r'$x=0.55$', fontsize = 18)
    ax[1].text(5.1, 0.8, r'$x=0.58$', fontsize = 18)
    ax[1].text(5.8, 0.33, r'$x=0.61$', fontsize = 18)
    ax[1].text(7.0, 0.12, r'$x=0.65$', fontsize = 18)
    ax[1].text(8.9, 0.029, r'$x=0.70$', fontsize = 18)
    ax[1].text(8.0, 0.01, r'$x=0.75 \,(i=0)$', fontsize = 18)

    ax[2].text(2.0, 1.2e5, r'$x=0.005 \, (i=17)$', fontsize = 18)
    ax[2].text(4.0, 0.63e5, r'$x=0.008$', fontsize = 18)
    ax[2].text(6.3, 3.1e4, r'$x=0.013$', fontsize = 18)
    ax[2].text(8.2, 1.5e4, r'$x=0.018$', fontsize = 18)
    ax[2].text(10.5, 7.8e3, r'$x=0.025$', fontsize = 18)
    ax[2].text(16.5, 3.6e3, r'$x=0.035$', fontsize = 18)
    ax[2].text(22.0, 1.9e3, r'$x=0.05$', fontsize = 18)
    ax[2].text(30.0, 0.9e3, r'$x=0.07$', fontsize = 18)
    ax[2].text(40, 4.7e2, r'$x=0.09$', fontsize = 18)
    ax[2].text(52.0, 2.2e2, r'$x=0.11$', fontsize = 18)
    ax[2].text(56.0, 1.1e2, r'$x=0.14$', fontsize = 18)
    ax[2].text(74.0, 5.4e1, r'$x=0.18$', fontsize = 18)
    ax[2].text(75.0, 2.5e1, r'$x=0.23$', fontsize = 18)
    ax[2].text(1.7, 12.2, r'$x=0.28$', fontsize = 18)
    ax[2].text(2.8, 6.0, r'$x=0.35$', fontsize = 18)
    ax[2].text(4.5, 2.9, r'$x=0.45$', fontsize = 18)
    ax[2].text(5.8, 1.35, r'$x=0.55$', fontsize = 18)
    ax[2].text(6.2, 0.7, r'$x=0.68 \,(i=0)$', fontsize = 18)

    ax[1].semilogy()
    ax[1].semilogx()
    ax[2].semilogy()
    ax[2].semilogx()

    ax[1].set_xlim(1.5, 5e2)
    ax[1].set_ylim(0.003, 4e4)
    ax[2].set_xlim(1.5, 150.0)
    ax[2].set_ylim(0.5, 3e5)

    ax[1].tick_params(axis = 'both', labelsize = 18)
    ax[2].tick_params(axis = 'both', labelsize = 18)

    ax[1].yaxis.set_tick_params(which = 'major', length = 6)
    ax[1].yaxis.set_tick_params(which = 'minir', length = 3)
    ax[1].set_yticks([0.01, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4])
    ax[1].set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
    ax[1].yaxis.set_minor_locator(locmin)
    ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[1].yaxis.set_tick_params(which = 'minor', length = 3)

    ax[2].yaxis.set_tick_params(which = 'major', length = 6)
    ax[2].yaxis.set_tick_params(which = 'minir', length = 3)
    ax[2].set_yticks([1.0, 10.0, 1e2, 1e3, 1e4, 1e5])
    ax[2].set_yticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])
    locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
    ax[2].yaxis.set_minor_locator(locmin)
    ax[2].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[2].yaxis.set_tick_params(which = 'minor', length = 3)

    ax[1].xaxis.set_tick_params(which = 'major', length = 6)
    ax[1].xaxis.set_tick_params(which = 'minor', length = 3)
    ax[1].set_xticks([2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2])
    xticklabels=[2e0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2]
    ax[1].set_xticklabels(['$%0.0f$'%x for x in xticklabels])

    ax[2].xaxis.set_tick_params(which = 'major', length = 6)
    ax[2].xaxis.set_tick_params(which = 'minor', length = 3)
    ax[2].set_xticks([2e0, 5e0, 1e1, 2e1, 5e1, 1e2])
    xticklabels=[2e0, 5e0, 1e1, 2e1, 5e1, 1e2]
    ax[2].set_xticklabels(['$%0.0f$'%x for x in xticklabels])

    ax[1].set_ylabel(r'$F_2^d$'+r'$\, \times\, 2^{\, i}$', size = 24)
    ax[1].set_xlabel(r'$Q^2 \: \rm{(GeV^2)}$', size = 24)
    ax[2].set_ylabel(r'$F_2^d/F_2^p$'+r'$\, \times\, 2^{\, i}$', size = 24)
    ax[2].set_xlabel(r'$Q^2 \: \rm{(GeV^2)}$', size = 24)

    ax[1].yaxis.set_ticks_position('both')
    ax[1].xaxis.set_ticks_position('both')
    ax[1].tick_params(axis = 'y', which = 'both', labelleft = 'on', labelright = 'off')
    ax[1].tick_params(axis = 'x', which = 'both', labeltop = 'off', labelbottom = 'on')

    ax[2].yaxis.set_ticks_position('both')
    ax[2].xaxis.set_ticks_position('both')
    ax[2].tick_params(axis = 'y', which = 'both', labelleft = 'on', labelright = 'off')
    ax[2].tick_params(axis = 'x', which = 'both', labeltop = 'off', labelbottom = 'on')

    ax[1].legend([SLAC, BCDMS], [r'$\mathrm{SLAC}$', r'$\mathrm{BCDMS}$'], fontsize = 24, frameon = 0)
    py.tight_layout()

    ax[2].legend([NMC], [r'$\mathrm{NMC}$'], fontsize = 24, frameon = 0)
    py.savefig('%s/gallery/dis-deuteron-%d.png' % (wdir, istep), dpi = dpi)
    py.close()

def plot_hera_318(wdir, data, istep, dpi):
    data_hera = {}
    try:
        data_hera[10026] = data[10026] ## HERA e+p \sqrt(s)=318 NC
        data_hera[10030] = data[10030] ## HERA e-p \sqrt(s)=318 NC
        data_hera[10031] = data[10031] ## HERA e+p \sqrt(s)=318 CC
        data_hera[10032] = data[10032] ## HERA e-p \sqrt(s)=318 CC
    except:
        sys.exit('HERA keys 10026, 10030, 10031 and 10032 do not all exist in data["idis"]')

    df_hera = {}
    for key in data_hera:
        df_hera[key] = pd.DataFrame(data_hera[key])

    df_hera_query = {}
    df_hera_plot = {}
    for key in df_hera:
        df_hera_query[key] = {}
        df_hera_plot[key] = {}

    # HERA e+p NC \sqrt{s} = 318 GeV
    df_hera_query[10026] = {}
    df_hera_query[10026][24] = df_hera[10026].query('X < 3.5e-5') ## 2.47e-05,2.928e-05,3.088e-05 2^24
    df_hera_query[10026][23] = df_hera[10026].query('X>3.5e-5 and X < 5.5e-5') ## x = 3.66e-05,4.06e-05,4.09e-05,4.323e-05,4.60e-05,5e-05,5.124e-05,5.22e-05,5.31e-05 2^23
    df_hera_query[10026][22] = df_hera[10026].query('X>5.5e-5 and X < 9.5e-5') ## x = 5.92e-05,6.176e-05,6.83e-05,7.32e-05,7.54e-05,8.e-05,8.029e-05,8.18e-05,8.55e-05,
    df_hera_query[10026][21] = df_hera[10026].query('X>9.5e-5 and X < 1.5e-4') ## x = 9.515e-05,9.86e-05,1.0499e-04,1.118e-04,1.2443e-04,1.29e-04,1.30e-04,1.39e-04,1.392e-04
    df_hera_query[10026][20] = df_hera[10026].query('X>1.5e-4 and X < 2.5e-4') ## x = 1.61e-04,1.741e-04,1.821e-04,2.0e-04,2.09e-04,2.276e-04,2.37e-04 2^20
    df_hera_query[10026][19] = df_hera[10026].query('X>2.5e-4 and X < 3.5e-4') ## x = 2.68e-04,2.9e-04,3.14e-04,3.2e-04,3.28e-04,3.45e-04
    df_hera_query[10026][18] = df_hera[10026].query('X>3.5e-4 and X < 5.5e-4') ## x = 3.55e-04,3.88e-04,4.1e-04,4.603e-04,5e-04,5.31e-04 2^18
    df_hera_query[10026][17] = df_hera[10026].query('X>5.5e-4 and X < 6.6e-4') ## x = 5.90e-04,5.92e-04,6.16e-04,6.341e-04,6.57e-04
    df_hera_query[10026][16] = df_hera[10026].query('X>7.5e-4 and X < 8.5e-4') ## x = 7.0e-04,8.0e-04 2^16
    df_hera_query[10026][15] = df_hera[10026].query('X>8.5e-4 and X < 0.98e-3') ## x = 8.6e-04,9.2e-04,9.22e-04,9.4e-04 2^15
    df_hera_query[10026][14] = df_hera[10026].query('X>0.98e-3 and X < 1.9e-3') ## x = 0.001,0.0011,0.00124,0.0013,0.0015, 0.0016, 0.00172,0.00188
    df_hera_query[10026][13] = df_hera[10026].query('X>0.0019 and X < 0.0030') ## x = 0.002,0.00212,0.0025,0.0026,0.0027 2^13
    df_hera_query[10026][12] = df_hera[10026].query('X>0.0030 and X < 0.0040') ## x = 00.0032,0.033,0.0039 2^12
    df_hera_query[10026][11] = df_hera[10026].query('X>4.8e-3 and X < 5.6e-3') ## x = 0.0050 x=0.053, 0.0066 2^11
    df_hera_query[10026][10] = df_hera[10026].query('X>7.9e-3 and X < 8.6e-3') ## x = 0.0080 x=0.0085 !check if it is possible to see both
    df_hera_query[10026][9] = df_hera[10026].query('X>1.0e-2 and X < 1.5e-2') ## x = 1.05 x=0.013 x=0.014 2^9
    df_hera_query[10026][8] = df_hera[10026].query('X>1.9e-2 and X < 2.2e-2') ## x = 0.020 x=0.0219  !check if it is possible to see both
    df_hera_query[10026][7] = df_hera[10026].query('X>3.0e-2 and X < 3.4e-2') ## x = 0.032 2^7
    df_hera_query[10026][6] = df_hera[10026].query('X>4.8e-2 and X < 5.6e-2') ## x = 0.050, x = 0.047 2^6
    df_hera_query[10026][5] = df_hera[10026].query('X>7.9e-2 and X < 8.8e-2') ## x = 0.080 x=0.0875 !check if it is possible to see both
    df_hera_query[10026][4] = df_hera[10026].query('X>12e-2 and X < 14e-2') ## x = 0.13 2^4
    df_hera_query[10026][3] = df_hera[10026].query('X>17e-2 and X < 19e-2') ## x = 0.18 2^3
    df_hera_query[10026][2] = df_hera[10026].query('X>23e-2 and X < 26e-2') ## x = 0.25 2^2
    df_hera_query[10026][1] = df_hera[10026].query('X>38.0e-2 and X < 42.0e-2') ## x = 0.40 2^1
    df_hera_query[10026][0] = df_hera[10026].query('X>62.0e-2 and X < 66e-2') ## x = 0.65 2^0

    df_hera_plot[10026] = {}
    df_hera_plot[10026]['theory'] = {}
    df_hera_plot[10026]['Q2'] = {}
    df_hera_plot[10026]['value'] = {}
    df_hera_plot[10026]['alpha'] = {}
    for key in df_hera_query[10026]:
        df_hera_plot[10026]['theory'][key] = df_hera_query[10026][key]['thy']
        df_hera_plot[10026]['Q2'][key] = df_hera_query[10026][key]['Q2']
        df_hera_plot[10026]['value'][key] = df_hera_query[10026][key]['value']
        df_hera_plot[10026]['alpha'][key] = df_hera_query[10026][key]['alpha']

    ## HERA 5 e-p NC \sqrt{s} = 318 GeV
    df_hera_query[10030] = {}
    df_hera_query[10030][16] = df_hera[10030].query('X>7.5e-4 and X < 8.5e-4') #x=0.0008 2^16
    df_hera_query[10030][14] = df_hera[10030].query('X>0.00098 and X < 0.002') #x=0.0013,0.0015, 0.0016 2^14
    df_hera_query[10030][13] = df_hera[10030].query('X>0.0019 and X < 0.0030') #x=0.002,0.0026, 2^12
    df_hera_query[10030][12] = df_hera[10030].query('X>0.0030 and X < 0.0040') #x=0.0032,0.033,0.0039 2^13
    df_hera_query[10030][11] = df_hera[10030].query('X>0.0048 and X < 0.0068') #x=0.005,0.0053, 0.0066 2^11
    df_hera_query[10030][10] = df_hera[10030].query('X>7.9e-3 and X < 8.6e-3') #x=0.0080 x=0.0085 2^10
    df_hera_query[10030][9] = df_hera[10030].query('X>1.0e-2 and X < 1.5e-2') #x=1.05 x=0.013  x=0.014 2^9
    df_hera_query[10030][8] = df_hera[10030].query('X>1.9e-2 and X < 2.2e-2') #x=0.020 x=0.0219 2^8
    df_hera_query[10030][7] = df_hera[10030].query('X>3.0e-2 and X < 3.4e-2') #x=0.032 2^7
    df_hera_query[10030][6] = df_hera[10030].query('X>4.8e-2 and X < 5.6e-2') #x=0.050 #x=0.0547  2^6
    df_hera_query[10030][5] = df_hera[10030].query('X>7.9e-2 and X < 8.8e-2') #x=0.080 x=0.0875 2^5
    df_hera_query[10030][4] = df_hera[10030].query('X>12e-2 and X < 14e-2') #x=0.13 2^4
    df_hera_query[10030][3] = df_hera[10030].query('X>17e-2 and X < 19e-2') #x=0.18 2^3
    df_hera_query[10030][2] = df_hera[10030].query('X>23e-2 and X < 26e-2') #x=0.25   2^2
    df_hera_query[10030][1] = df_hera[10030].query('X>38.0e-2 and X < 42.0e-2') #x=0.40 2^1
    df_hera_query[10030][0] = df_hera[10030].query('X>62.0e-2 and X < 66e-2') #x=0.65 2^0

    df_hera_plot[10030] = {}
    df_hera_plot[10030]['theory'] = {}
    df_hera_plot[10030]['Q2'] = {}
    df_hera_plot[10030]['value'] = {}
    df_hera_plot[10030]['alpha'] = {}
    for key in df_hera_query[10030]:
        df_hera_plot[10030]['theory'][key] = df_hera_query[10030][key]['thy']
        df_hera_plot[10030]['Q2'][key] = df_hera_query[10030][key]['Q2']
        df_hera_plot[10030]['value'][key] = df_hera_query[10030][key]['value']
        df_hera_plot[10030]['alpha'][key] = df_hera_query[10030][key]['alpha']

    ## HERA 6 e+p \sqrt{s} = 318 GeV CC
    df_hera_query[10031] = {}
    df_hera_query[10031][7] = df_hera[10031].query('X>7.9e-3 and X < 8.1e-3') ## x = 0.008
    df_hera_query[10031][6] = df_hera[10031].query('X>1.0e-2 and X < 1.5e-2') ## x = 0.013
    df_hera_query[10031][5] = df_hera[10031].query('X>3.0e-2 and X < 3.4e-2') ## x = 0.032
    df_hera_query[10031][4] = df_hera[10031].query('X>7.9e-2 and X < 8.6e-2') ## x = 0.080
    df_hera_query[10031][3] = df_hera[10031].query('X>12e-2 and X < 14e-2') ## x = 0.13
    df_hera_query[10031][2] = df_hera[10031].query('X>23e-2 and X < 26e-2') ## x = 0.25
    df_hera_query[10031][1] = df_hera[10031].query('X>38.0e-2 and X < 42.0e-2') ## x = 0.40

    df_hera_plot[10031] = {}
    df_hera_plot[10031]['theory'] = {}
    df_hera_plot[10031]['Q2'] = {}
    df_hera_plot[10031]['value'] = {}
    df_hera_plot[10031]['alpha'] = {}
    for key in df_hera_query[10031]:
        df_hera_plot[10031]['theory'][key] = df_hera_query[10031][key]['thy']
        df_hera_plot[10031]['Q2'][key] = df_hera_query[10031][key]['Q2']
        df_hera_plot[10031]['value'][key] = df_hera_query[10031][key]['value']
        df_hera_plot[10031]['alpha'][key] = df_hera_query[10031][key]['alpha']

    #HERA7 e-p CC \sqrt{s}=318 GeV
    df_hera_query[10032] = {}
    df_hera_query[10032][7] = df_hera[10032].query('X>0.0079 and X < 0.0081') ## x = 0.008 2^7
    df_hera_query[10032][6] = df_hera[10032].query('X>0.010 and X < 0.015') ## x = 0.013 2^6
    df_hera_query[10032][5] = df_hera[10032].query('X>0.030 and X < 0.034') ## x = 0.032 2^5
    df_hera_query[10032][4] = df_hera[10032].query('X>0.079 and X < 0.081') ## x = 0.080 2^4
    df_hera_query[10032][3] = df_hera[10032].query('X>0.12 and X < 0.14') ## x = 0.13  2^3
    df_hera_query[10032][2] = df_hera[10032].query('X>0.23 and X < 0.27') ## x = 0.25 2^2
    df_hera_query[10032][1] = df_hera[10032].query('X>0.38 and X < 0.42') ## x = 0.40 2^1
    df_hera_query[10032][0] = df_hera[10032].query('X>0.63 and X < 0.66') ## x = 0.65 2^0

    df_hera_plot[10032] = {}
    df_hera_plot[10032]['theory'] = {}
    df_hera_plot[10032]['Q2'] = {}
    df_hera_plot[10032]['value'] = {}
    df_hera_plot[10032]['alpha'] = {}
    for key in df_hera_query[10032]:
        df_hera_plot[10032]['theory'][key] = df_hera_query[10032][key]['thy']
        df_hera_plot[10032]['Q2'][key] = df_hera_query[10032][key]['Q2']
        df_hera_plot[10032]['value'][key] = df_hera_query[10032][key]['value']
        df_hera_plot[10032]['alpha'][key] = df_hera_query[10032][key]['alpha']

    nrows, ncols = 1, 2
    fig = py.figure(figsize = (ncols * 9.0, nrows * 11.0))
    ax = {}
    ax[1] = fig.add_subplot(nrows, ncols, 1)
    ax[2] = fig.add_subplot(nrows, ncols, 2)

    ## left subplot plots neutral current
    ## HERA e+p NC 318
    for key in df_hera_plot[10026]['value']:
        HERA_10026 = ax[1].errorbar(df_hera_plot[10026]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10026]['value'][key], df_hera_plot[10026]['alpha'][key], \
                                    marker = '.', color = 'black', markersize = 7, linestyle = 'none')
        ax[1].plot(df_hera_plot[10026]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10026]['theory'][key], linestyle = 'solid', color = 'black')

    ## HERA5 e-p Nc 318
    for key in df_hera_plot[10030]['value']:
        HERA_10030 = ax[1].errorbar(df_hera_plot[10030]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10030]['value'][key], df_hera_plot[10030]['alpha'][key], \
                                    marker = 'D', color = 'firebrick', markersize = 8, linestyle = 'none')
        ax[1].plot(df_hera_plot[10030]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10030]['theory'][key], linestyle = 'solid', color = 'firebrick')

    ## right subplot plots charged current
    ## HERA 6 e+p CC 318
    for key in df_hera_plot[10031]['value']:
        HERA_10031 = ax[2].errorbar(df_hera_plot[10031]['Q2'][key], 5.0 ** float(key) * df_hera_plot[10031]['value'][key], df_hera_plot[10031]['alpha'][key], \
                                    marker = '^', color = 'royalblue', fillstyle = 'none', markersize = 8, linestyle = 'none')
        ax[2].plot(df_hera_plot[10031]['Q2'][key], 5.0 ** float(key) * df_hera_plot[10031]['theory'][key], linestyle = 'solid', color = 'royalblue')

    ## HERA 7 e-p CC 318
    for key in df_hera_plot[10032]['value']:
        HERA_10032 = ax[2].errorbar(df_hera_plot[10032]['Q2'][key], 5.0 ** float(key) * df_hera_plot[10032]['value'][key], df_hera_plot[10032]['alpha'][key], \
                                    marker = 's', color = 'darkgreen', fillstyle = 'none', markersize = 8, linestyle = 'none')
        ax[2].plot(df_hera_plot[10032]['Q2'][key], 5.0 ** float(key) * df_hera_plot[10032]['theory'][key], linestyle = 'solid', color = 'darkgreen')

    ax[1].text(3.4,1.65e7, r'$x=2.8\cdot 10^{-5}\, (i=24)$', fontsize = 18)
    ax[1].text(5.5,8.65e6, r'$x=4.6\cdot 10^{-5}$', fontsize = 18)
    ax[1].text(8.0,4.6e6, r'$x=7.3\cdot 10^{-5}$', fontsize = 18)
    ax[1].text(15.0,2.7e6, r'$x=1.2\cdot 10^{-4}$', fontsize = 18)
    ax[1].text(22.0,1.4e6, r'$x=1.9\cdot 10^{-4}$', fontsize = 18)
    ax[1].text(35.0,7.5e5, r'$x=3.1\cdot 10^{-4}$', fontsize = 18)
    ax[1].text(45.0,3.9e5, r'$x=4.4\cdot 10^{-4}$', fontsize = 18)
    ax[1].text(60.0,1.8e5, r'$x=6.2\cdot 10^{-4}$', fontsize = 18)
    ax[1].text(80.0,8.5e4, r'$x=7.7\cdot 10^{-4}$', fontsize = 18)
    ax[1].text(86.0,4e4, r'$x=9.1 \cdot 10^{-4}$', fontsize = 18)
    ax[1].text(160.0,2e4, r'$x=0.0014$', fontsize = 18)
    ax[1].text(260.0,1e4, r'$x=0.0024$', fontsize = 18)
    ax[1].text(440.0,4500, r'$x=0.0035$', fontsize = 18)
    ax[1].text(740.0,2000, r'$x=0.0056$', fontsize = 18)
    ax[1].text(1.5e3,970, r'$x=0.0083$', fontsize = 18)
    ax[1].text(1.7e3,410, r'$x=0.013$', fontsize = 18)
    ax[1].text(2.7e3,190, r'$x=0.021$', fontsize = 18)
    ax[1].text(3.9e3,77, r'$x=0.032$', fontsize = 18)
    ax[1].text(6.5e3,35, r'$x=0.052$', fontsize = 18)
    ax[1].text(1.15e4,19, r'$x=0.084$', fontsize = 18)
    ax[1].text(12.0,5.5, r'$x=0.13$', fontsize = 18)
    ax[1].text(18.0,2.25, r'$x=0.18$', fontsize = 18)
    ax[1].text(36.0,0.97, r'$x=0.25$', fontsize = 18)
    ax[1].text(43.0,0.25, r'$x=0.40$', fontsize = 18)
    ax[1].text(30.0,0.02, r'$x=0.65 \, (i=0)$', fontsize = 18)

    ax[2].text(370.0, 1.4e5, r'$x=0.008 \, (i=7)$', fontsize = 18, color = 'darkgreen')
    ax[2].text(1200.0, 1.3e4, r'$x=0.013 \, (i=6)$', fontsize = 18, color = 'darkgreen')
    ax[2].text(3700.0, 2e3, r'$x=0.032 \, (i=5)$', fontsize = 18, color = 'darkgreen')
    ax[2].text(6000.0, 400, r'$x=0.08 \, (i=4)$', fontsize = 18, color = 'darkgreen')
    ax[2].text(9800.0, 70, r'$x=0.13 \, (i=3)$', fontsize = 18, color = 'darkgreen')
    ax[2].text(280.0, 11.5, r'$x=0.25 \, (i=2)$', fontsize = 18, color = 'darkgreen')
    ax[2].text(400.0, 1.5, r'$x=0.40 \, (i=1)$', fontsize = 18, color = 'darkgreen')
    ax[2].text(1300.0, 0.035, r'$x=0.65 \, (i=0) $', fontsize = 18, color = 'darkgreen')

    ax[2].text(550.0, 0.6e5, r'$x=0.008 \, (i=7)$', fontsize = 18, color = 'royalblue')
    ax[2].text(1200.0, 0.8e4, r'$x=0.013 \, (i=6)$', fontsize = 18, color = 'royalblue')
    ax[2].text(3700.0, 1000, r'$x=0.032 \, (i=5)$', fontsize = 18, color = 'royalblue')
    ax[2].text(4200.0, 180, r'$x=0.08 \, (i=4)$', fontsize = 18, color = 'royalblue')
    ax[2].text(0.7e4,23, r'$x=0.13 \, (i=3)$', fontsize = 18, color = 'royalblue')
    ax[2].text(270.0, 5.2, r'$x=0.25 \, (i=2)$', fontsize = 18, color = 'royalblue')
    ax[2].text(430.0, 0.4, r'$x=0.40 \, (i=1)$', fontsize = 18, color = 'royalblue')

    ax[1].semilogy()
    ax[1].semilogx()
    ax[2].semilogy()
    ax[2].semilogx()

    ax[1].set_xlim(1.0, 8e4)
    ax[1].set_ylim(0.004, 4e7)
    # ax[2].set_xlim(1.5, 150.0)
    # ax[2].set_ylim(0.5, 3e5)

    ax[1].set_xticks([1.0, 10.0, 1e2, 1e3, 1e4, 8e4])
    ax[1].set_xticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$8\cdot10^4$'])
    locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
    ax[1].yaxis.set_minor_locator(locmin)
    ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[1].yaxis.set_tick_params(which = 'minor', length = 3)
    ax[1].set_yticks([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
    ax[1].set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$', r'$10^7$'])

    ax[2].set_xticks([2e2, 1e3, 1e4, 4e4])
    ax[2].set_xticklabels([r'$200$', r'$10^3$', r'$10^4$', r'$4\cdot10^4$'])
    ax[2].set_yticks([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5])
    ax[2].set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])

    ax[1].tick_params(axis = 'both', labelsize = 22)
    ax[2].tick_params(axis = 'both', labelsize = 22)

    ax[1].set_ylabel(r'$\sigma_r^{p,NC}$'+r'$\, \times\, 2^{\, i}$', size = 28)
    ax[1].set_xlabel(r'$Q^2 \: \rm{(GeV^2)}$', size = 28)
    ax[2].set_ylabel(r'$\sigma_r^{p,CC}$'+r'$\, \times\, 5^{\, i}$', size = 28)
    ax[2].set_xlabel(r'$Q^2 \: \rm{(GeV^2)}$', size = 28)

    ax[1].xaxis.set_tick_params(which = 'major', length = 6)
    ax[1].xaxis.set_tick_params(which = 'minor', length = 3)
    ax[1].yaxis.set_tick_params(which = 'major', length = 6)
    ax[1].yaxis.set_tick_params(which = 'minor', length = 3)
    ax[2].xaxis.set_tick_params(which = 'major', length = 6)
    ax[2].xaxis.set_tick_params(which = 'minor', length = 3)
    ax[2].yaxis.set_tick_params(which = 'major', length = 6)
    ax[2].yaxis.set_tick_params(which = 'minor', length = 3)

    ax[1].yaxis.set_ticks_position('both')
    ax[1].xaxis.set_ticks_position('both')
    ax[1].tick_params(axis = 'y', which = 'both', labelleft = 'on', labelright = 'off')
    ax[1].tick_params(axis = 'x', which = 'both', labeltop = 'off', labelbottom = 'on')

    ax[2].yaxis.set_ticks_position('both')
    ax[2].xaxis.set_ticks_position('both')
    ax[2].tick_params(axis = 'y', which = 'both', labelleft = 'on', labelright = 'off')
    ax[2].tick_params(axis = 'x', which = 'both', labeltop = 'off', labelbottom = 'on')

    ax[1].text(2000.0, 0.9e6, r'$\sqrt{s}=318 \,\rm{GeV}$', fontsize = 24)
    ax[2].text(0.72e4, 2e4, r'$\sqrt{s}=318 \,\rm{GeV}$', fontsize = 24)
    ax[1].legend([HERA_10026, HERA_10030], [r'$\mathrm{HERA~NC}~e^+p$', r'$\mathrm{HERA~NC}~e^-p$'], fontsize = 24, frameon = 0)
    ax[2].legend([HERA_10031, HERA_10032], [r'$\mathrm{HERA~CC}~e^+p$', r'$\mathrm{HERA~CC}~e^-p$'], fontsize = 24, frameon = 0)
    py.tight_layout()

    py.savefig('%s/gallery/dis-hera-318-%d.png' % (wdir, istep), dpi = dpi)
    py.close()

def plot_hera_other_than_318(wdir, data, istep, dpi):
    data_hera = {}
    try:
        data_hera[10027] = data[10027] ## HERA e+p \sqrt(s)=300 NC
        data_hera[10028] = data[10028] ## HERA e+p \sqrt(s)=251 NC
        data_hera[10029] = data[10029] ## HERA e+p \sqrt(s)=225 NC
    except:
        sys.exit('HERA keys 10027, 10028 and 10029 do not all exist in data["idis"]')

    df_hera = {}
    for key in data_hera:
        df_hera[key] = pd.DataFrame(data_hera[key])

    df_hera_query = {}
    df_hera_plot = {}
    for key in df_hera:
        df_hera_query[key] = {}
        df_hera_plot[key] = {}

    ## HERA e+p NC \sqrt{s} = 300 GeV
    df_hera_query[10027] = {}
    df_hera_query[10027][29] = df_hera[10027].query('X<3.5e-5 ') ## x = 3.27e-05
    df_hera_query[10027][28] = df_hera[10027].query('X>3.5e-5 and X < 9.5e-5') ## x = 4.09e-05,5e-05,5.73e-05,8e-05,8.18e-05
    df_hera_query[10027][27] = df_hera[10027].query('X>9.5e-5 and X < 2.5e-4') ## x = 9.86e-05,1.3e-04,1.39e-04,1.61e-04,2e-04,2.46e-04
    df_hera_query[10027][26] = df_hera[10027].query('X>2.5e-4 and X < 3.5e-4') ## x = 2.68e-04,3.2e-04,3.28e-04,3.35e-04
    df_hera_query[10027][25] = df_hera[10027].query('X>3.5e-4 and X < 5.5e-4') ## x = 4.1e-04,5e-04
    df_hera_query[10027][24] = df_hera[10027].query('X>5.5e-4 and X < 6.6e-4') ## x = 5.74e-04
    df_hera_query[10027][23] = df_hera[10027].query('X>6.9e-4 and X < 8.5e-4') ## x = 8e-04
    df_hera_query[10027][21] = df_hera[10027].query('X>0.98e-3 and X < 0.0014') ## x = 1.3e-03
    df_hera_query[10027][19] = df_hera[10027].query('X>0.0019 and X < 0.0029') ## x = 2e-03,2.12e-03
    df_hera_query[10027][18] = df_hera[10027].query('X>0.0029 and X < 0.0040') ## x = 3.2e-03
    df_hera_query[10027][17] = df_hera[10027].query('X>4.0e-3 and X < 5.55e-3') ## x = 5e-03
    df_hera_query[10027][15] = df_hera[10027].query('X>7.95e-3 and X < 8.9e-3') ## x = 8e-03,8.5e-03
    df_hera_query[10027][12] = df_hera[10027].query('X>1.25e-2 and X < 1.5e-2') ## x = 1.3e-02,1.4e-02
    df_hera_query[10027][9] = df_hera[10027].query('X>1.99e-2 and X < 2.2e-2') ## x = 2e-02
    df_hera_query[10027][7] = df_hera[10027].query('X>3.0e-2 and X < 3.4e-2') ## x = 3.20e-02
    df_hera_query[10027][6] = df_hera[10027].query('X>4.8e-2 and X < 5.6e-2') ## x = 0.050
    df_hera_query[10027][5] = df_hera[10027].query('X>7.9e-2 and X < 8.8e-2') ## x = 0.080
    df_hera_query[10027][4] = df_hera[10027].query('X>12e-2 and X < 14e-2') ## x = 0.13
    df_hera_query[10027][3] = df_hera[10027].query('X>17e-2 and X < 19e-2') ## x = 0.18
    df_hera_query[10027][2] = df_hera[10027].query('X>23e-2 and X < 26e-2') ## x = 0.25
    df_hera_query[10027][1] = df_hera[10027].query('X>30e-2 and X < 45e-2') ## x = 0.4

    df_hera_plot[10027] = {}
    df_hera_plot[10027]['theory'] = {}
    df_hera_plot[10027]['Q2'] = {}
    df_hera_plot[10027]['value'] = {}
    df_hera_plot[10027]['alpha'] = {}
    for key in df_hera_query[10027]:
        df_hera_plot[10027]['theory'][key] = df_hera_query[10027][key]['thy']
        df_hera_plot[10027]['Q2'][key] = df_hera_query[10027][key]['Q2']
        df_hera_plot[10027]['value'][key] = df_hera_query[10027][key]['value']
        df_hera_plot[10027]['alpha'][key] = df_hera_query[10027][key]['alpha']

    ## HERA e+p NC \sqrt{s} = 251 GeV
    df_hera_query[10028] = {}
    df_hera_query[10028][28] = df_hera[10028].query('X>3.5e-5 and X < 9.5e-5') ## x = 3.72e-05,4.15e-05,4.65e-05,5.19e-05,,5.8e-05,6.51e-05,7.27e-05,8.12e-05,9.21e-05,9.31e-05
    df_hera_query[10028][27] = df_hera[10028].query('X>9.5e-5 and X < 2.5e-4') ## x = 1.038e-04,1.062e-04,1.16e-04,1.21e-04,1.315e-04,1.35e-04, x=1.509e-04,1.517e-04,1.582e-04,1.63e-04,1.71e-04,1.765e-04,1.83e-04,1.973e-04,2.013e-04,2.233e-04,2.236e-04,2.28e-04,2.4e-04,2.492e-04
    df_hera_query[10028][26] = df_hera[10028].query('X>2.5e-4 and X < 3.5e-4') ## x = 2.58e-04,2.617e-04,2.785e-04,2.792e-04,2.954e-04,2.99e-04,3.115e-04,3.156e-04,3.422e-04,3.43e-04,3.481e-04
    df_hera_query[10028][25] = df_hera[10028].query('X>3.5e-4 and X < 5.5e-4') ## x = 3.642e-04,3.722e-04,3.945e-04,4.153e-04,4.22e-04,4.23e-04,4.33e-04,4.552e-04,4.642e-04,4.831e-04,4.932e-04,5.261e-04,5.29e-04,5.409e-04
    df_hera_query[10028][24] = df_hera[10028].query('X>5.5e-4 and X < 6.6e-4') ## x = 5.70e-04,5.92e-04,6.039e-04,6.07e-04,6.165e-04,6.514e-04
    df_hera_query[10028][23] = df_hera[10028].query('X>6.9e-4 and X < 8.5e-4') ## x = 7.e-04,7.268e-04,7.587e-04,7.63e-04,7.636e-04,8.052e-04,8.123e-04,8.375e-04,8.384e-04
    df_hera_query[10028][22] = df_hera[10028].query('X>8.5e-4 and X < 0.98e-3') ## x =  8.8e-04,9.1e-04,9.206e-04,9.344e-04,  9.54500000e-04
    df_hera_query[10028][21] = df_hera[10028].query('X>0.98e-3 and X < 0.0013') ## x = 1.e-03,1.007e-03,1.030e-03,1.044e-03,1.062e-03,1.117e-03,1.14e-03,1.184e-03,1.23e-03,1.246e-03,1.27e-03,1.273e-03
    df_hera_query[10028][20] = df_hera[10028].query('X>1.3e-3 and X < 1.9e-03') ## 1 . 3660e-03,1.392e-03,1.397e-03,1.409-03,1.46e-03,1.479e-03,1.578e-03,1.585e-03,1.591e-03
    df_hera_query[10028][19] = df_hera[10028].query('X>0.0019 and X < 0.0029') ## x =  1.97300000e-03,2.08900000e-03,2.11000000e-03,2.20000000e-03,2.22700000e-03, 2.29000000e-03, 2.36700000e-03, 2.41600000e-03, 2.46600000e-03, 2.50000000e-03, 2.70000000e-03,   2.80000000e-03, 2.86400000e-03
    df_hera_query[10028][18] = df_hera[10028].query('X>0.0029 and X < 0.0040') ## x = 3e-03,3.05e-03,3.1e-03,3.2e-03,3.288e-03,3.3e-03, 3.4e-03, 3.452e-03, 3.5e-03, 3.60000000e-03, 3.62300000e-03, 3.70000000e-03, 3.81800000e-03, 3.90000000e-03
    df_hera_query[10028][17] = df_hera[10028].query('X>4.0e-3 and X < 5.55e-3') ## x = 4.1e-03, 4.11e-03, 4.2e-03, 4.3e-03, 4.439e-03, 4.5e-03, 4.6e-03, 4.8e-03, 5.2e-03, 5.5e-03
    df_hera_query[10028][16] = df_hera[10028].query('X>5.55e-3 and X < 7.95e-3') ## x = 5.60000000e-03, 5.72700000e-03, 5.75400000e-03, 5.8e-03, 5.9e-03, 5.918e-03, 6e-03, 6.1e-03,6.2e-03,6.4e-03,6.6e-03,6.9e-03,7.30000000e-03,7.39800000e-03,7.40000000e-03,7.60000000e-03,7.90000000e-03
    df_hera_query[10028][15] = df_hera[10028].query('X>7.95e-3 and X < 8.9e-3') ## x = 0.0080 x=0.0083 ,0.008877
    df_hera_query[10028][14] = df_hera[10028].query('X>9e-3 and X < 9.9e-3') ## x = 9.1e-03,9.3e-03,9.864e-03
    df_hera_query[10028][13] = df_hera[10028].query('X>9.9e-3 and X < 1.25e-2') ## x = 1e-02,1.04e-02,1.05e-02,1.09e-02,1.16e-02,1.21e-02
    df_hera_query[10028][12] = df_hera[10028].query('X>1.25e-2 and X < 1.5e-2') ## x = 1.3e-02,1.31e-02,1.350e-02,1.48e-02,1.49e-02
    df_hera_query[10028][11] = df_hera[10028].query('X>1.5e-2 and X < 1.7e-2') ## x = 1.51e-02,1.52e-02,1.61e-02,1.660e-02
    df_hera_query[10028][10] = df_hera[10028].query('X>1.7e-2 and X < 1.99e-2') ## x = 1.71e-02,1.85e-02,1.97e-02
    df_hera_query[10028][9] = df_hera[10028].query('X>1.99e-2 and X < 2.2e-2') ## x = 2e-02,2.01e-02,2.1e-02
    df_hera_query[10028][8] = df_hera[10028].query('X>2.3e-2 and X < 2.7e-2') ## x = 2.420e-02,2.61e-02
    df_hera_query[10028][7] = df_hera[10028].query('X>3.0e-2 and X < 3.4e-2') ## x = 3.20e-02,3.22e-02
    df_hera_query[10028][6] = df_hera[10028].query('X>4.8e-2 and X < 5.6e-2') ## x = 0.050
    df_hera_query[10028][5] = df_hera[10028].query('X>7.9e-2 and X < 8.8e-2') ## x = 0.080
    df_hera_query[10028][4] = df_hera[10028].query('X>12e-2 and X < 14e-2') ## x = 0.13
    df_hera_query[10028][3] = df_hera[10028].query('X>17e-2 and X < 19e-2') ## x = 0.18
    df_hera_query[10028][2] = df_hera[10028].query('X>23e-2 and X < 26e-2') ## x = 0.25

    df_hera_plot[10028] = {}
    df_hera_plot[10028]['theory'] = {}
    df_hera_plot[10028]['Q2'] = {}
    df_hera_plot[10028]['value'] = {}
    df_hera_plot[10028]['alpha'] = {}
    for key in df_hera_query[10028]:
        df_hera_plot[10028]['theory'][key] = df_hera_query[10028][key]['thy']
        df_hera_plot[10028]['Q2'][key] = df_hera_query[10028][key]['Q2']
        df_hera_plot[10028]['value'][key] = df_hera_query[10028][key]['value']
        df_hera_plot[10028]['alpha'][key] = df_hera_query[10028][key]['alpha']

    ## HERA e+p NC \sqrt{s} = 225 GeV
    df_hera_query[10029] = {}
    df_hera_query[10029][28] = df_hera[10029].query('X>3.5e-5 and X < 9.5e-5') ## x = 4.64e-05,5.26e-05, 5.8e-05,6.58e-05,7.59e-05,8.12e-05,9.21e-05
    df_hera_query[10029][27] = df_hera[10029].query('X>9.5e-5 and X < 2.5e-4') ## x = 1.062e-04,1.16e-04,1.315e-04,1.409e-04,1.43e-04,1.509e-04,1.517e-04,1.710e-04,1.830e-04,1.973e-04,2.0130e-04,2.236e-04,2.28e-04,2.40e-04
    df_hera_query[10029][26] = df_hera[10029].query('X>2.5e-4 and X < 3.5e-4') ## 2.58e-04,2.617e-04,2.785e-04,2.99e-04,3.156e-04,3.3e-04,3.422e-04,3.43e-04,3.481e-04
    df_hera_query[10029][25] = df_hera[10029].query('X>3.5e-4 and X < 5.5e-4') ## 3.642e-04,3.945e-04,4.136e-04,4.22e-04,4.23e-04,4.33e-04,4.552e-04,4.642e-04, 4.831e-04, 5.261e-04, 5.29e-04, 5.409e-04
    df_hera_query[10029][24] = df_hera[10029].query('X>5.5e-4 and X < 6.6e-4') ## 5.92e-04, 6.039e-04, 6.07e-04, 6.165e-04, 6.570e-04
    df_hera_query[10029][23] = df_hera[10029].query('X>6.9e-4 and X < 8.5e-4') ## 7e-04, 7.587e-04, 7.63e-04, 7.636e-04, 8.052e-04, 8.123e-04, 8.384e-04
    df_hera_query[10029][22] = df_hera[10029].query('X>8.5e-4 and X < 0.98e-3') ##  8.80000000e-04, 9.20600000e-04, 9.54500000e-04
    df_hera_query[10029][21] = df_hera[10029].query('X>0.98e-3 and X < 0.0013') ## 1e-03,1.007e-03, 1.044e-03,1.062e-03,1.184e-03,1.23e-03,1.270e-03,1.273e-03
    df_hera_query[10029][20] = df_hera[10029].query('X>1.3e-3 and X < 1.9e-03') ## 1.366e-03,1.392e-03,1.397e-03,1.4090e-03,1.479e-03,1.578e-03,1.585e-03,1.591e-03,1.71e-03,1.8e-03,1.812e-03,1.821e-03
    df_hera_query[10029][19] = df_hera[10029].query('X>0.0019 and X < 0.0029') ## 1.973e-03,2.089e-03,2.11e-03,2.2270e-03,2.29e-03,2.3670e-03,2.416e-03,2.466e-03,2.5e-03,2.731e-03,2.8e-03,2.864e-03
    df_hera_query[10029][18] = df_hera[10029].query('X>0.0029 and X < 0.0040') ## 3.05e-03,3.1e-03,3.17e-03,3.2e-03,3.288e-03,3.4e-03,3.452e-03,3.5e-03,3.6e-03,3.623e-03,3.818e-03,3.9e-03
    df_hera_query[10029][17] = df_hera[10029].query('X>4.0e-3 and X < 5.55e-3') ## 4.11e-03,4.20e-03,4.439e-03,4.5e-03,4.6e-03,4.8e-03,5e-03,5.2e-03
    df_hera_query[10029][16] = df_hera[10029].query('X>5.55e-3 and X < 7.95e-3') ## 5.727e-03,5.754e-03,5.8e-03,5.9e-03, 5.918e-03,6e-03,6.1e-03,6.6e-03,6.9e-03,7.398e-03,7.6e-03,7.9e-03
    df_hera_query[10029][14] = df_hera[10029].query('X>9e-3 and X < 9.9e-3') ## x = 9.1e-03, 9.3e-03, 9.864e-03
    df_hera_query[10029][13] = df_hera[10029].query('X>9.9e-3 and X < 1.25e-2') ## x = 1e-02, 1.05e-02, 1.16e-02, 1.21e-02
    df_hera_query[10029][12] = df_hera[10029].query('X>1.25e-2 and X < 1.5e-2') ## x = 1.3e-02, 1.31e-02, 1.36e-02, 1.48e-02
    df_hera_query[10029][11] = df_hera[10029].query('X>1.5e-2 and X < 1.7e-2') ## x = 1.51e-02, 1.52e-02, 1.61e-02
    df_hera_query[10029][10] = df_hera[10029].query('X>1.7e-2 and X < 1.99e-2') ## x = 1.71e-02,1.80e-02,1.85e-02,1.97e-02
    df_hera_query[10029][9] = df_hera[10029].query('X>1.99e-2 and X < 2.2e-2') ## x = 2e-02,2.01e-02,2.1e-02
    df_hera_query[10029][8] = df_hera[10029].query('X>2.3e-2 and X < 2.7e-2') ## x = 2.420e-02,2.61e-02
    df_hera_query[10029][7] = df_hera[10029].query('X>3.0e-2 and X < 3.4e-2') ## x = 3.20e-02,3.22e-02
    df_hera_query[10029][6] = df_hera[10029].query('X>4.8e-2 and X < 5.6e-2') ## x = 0.050
    df_hera_query[10029][5] = df_hera[10029].query('X>7.9e-2 and X < 8.8e-2') ## x = 0.080
    df_hera_query[10028][4] = df_hera[10028].query('X>12e-2 and X < 14e-2') ## x = 0.13
    df_hera_query[10028][3] = df_hera[10028].query('X>17e-2 and X < 19e-2') ## x = 0.18
    df_hera_query[10028][2] = df_hera[10028].query('X>23e-2 and X < 26e-2') ## x = 0.25
    df_hera_query[10028][1] = df_hera[10028].query('X>30e-2 and X < 45e-2') ## x = 0.4
    df_hera_query[10028][0] = df_hera[10028].query('X>60e-2 and X < 70e-2') ## x = 0.65

    df_hera_plot[10029] = {}
    df_hera_plot[10029]['theory'] = {}
    df_hera_plot[10029]['Q2'] = {}
    df_hera_plot[10029]['value'] = {}
    df_hera_plot[10029]['alpha'] = {}
    for key in df_hera_query[10029]:
        df_hera_plot[10029]['theory'][key] = df_hera_query[10029][key]['thy']
        df_hera_plot[10029]['Q2'][key] = df_hera_query[10029][key]['Q2']
        df_hera_plot[10029]['value'][key] = df_hera_query[10029][key]['value']
        df_hera_plot[10029]['alpha'][key] = df_hera_query[10029][key]['alpha']

    for key in df_hera_query[10028]:
        df_hera_plot[10028]['theory'][key] = df_hera_query[10028][key]['thy']
        df_hera_plot[10028]['Q2'][key] = df_hera_query[10028][key]['Q2']
        df_hera_plot[10028]['value'][key] = df_hera_query[10028][key]['value']
        df_hera_plot[10028]['alpha'][key] = df_hera_query[10028][key]['alpha']

    nrows, ncols = 1, 1
    fig = py.figure(figsize=(ncols * 12.0, nrows * 14.0))
    ax = fig.add_subplot(nrows, ncols, 1)

    ## HERA e+p NC \sqrt{s} = 300 GeV
    for key in df_hera_plot[10027]['value']:
        if key == 1: continue
        HERA_10027 = ax.errorbar(df_hera_plot[10027]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10027]['value'][key], df_hera_plot[10027]['alpha'][key], \
                                 marker = 's', color = 'darkgreen', fillstyle = 'none', markersize = 8, linestyle = 'none')
        ax.plot(df_hera_plot[10027]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10027]['theory'][key], linestyle = 'solid', color = 'darkgreen')

    ## HERA5 e+p NC \sqrt{s} = 251 GeV
    for key in df_hera_plot[10028]['value']:
        if (key == 22) or (key == 20) or (key <= 4): continue
        color = 'firebrick'
        if key == 24: color = 'black'
        HERA_10028 = ax.errorbar(df_hera_plot[10028]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10028]['value'][key], df_hera_plot[10028]['alpha'][key], \
                                 marker = '.', color = color, markersize = 8, linestyle = 'none')
        ax.plot(df_hera_plot[10028]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10028]['theory'][key], linestyle = 'solid', color = color)

    ## HERA e+p NC \sqrt{s} = 225 GeV
    for key in df_hera_plot[10029]['value']:
        HERA_10029 = ax.errorbar(df_hera_plot[10029]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10029]['value'][key], df_hera_plot[10029]['alpha'][key], \
                                 marker = '^', color = 'royalblue', fillstyle = 'none', markersize = 8, linestyle = 'none')
        ax.plot(df_hera_plot[10029]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10029]['theory'][key], linestyle = 'solid', color = 'royalblue')

    for key in range(5):
        ax.errorbar(df_hera_plot[10028]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10028]['value'][key], df_hera_plot[10028]['alpha'][key], \
                    marker = '^', color = 'royalblue', fillstyle = 'none', markersize = 8, linestyle = 'none')
        ax.plot(df_hera_plot[10028]['Q2'][key], 2.0 ** float(key) * df_hera_plot[10028]['theory'][key], linestyle = 'solid', color = 'royalblue')

    ax.text(2.5, 4.5e8, r'$x=3.3 \cdot 10^{-5}\, (i=29)$', fontsize = 18)
    ax.text(5.7, 2.5e8, r'$x=5.9 \cdot 10^{-5}$', fontsize = 18)
    ax.text(19.0, 1.7e8, r'$x=1.6 \cdot 10^{-4}$', fontsize = 18)
    ax.text(35.0, 8e7, r'$x=3.1 \cdot 10^{-4}$', fontsize = 18)
    ax.text(35.0, 4.4e7, r'$x=3.1 \cdot 10^{-4}$', fontsize = 18)
    ax.text(45.0, 2.1e7, r'$x=4.5 \cdot 10^{-4}$', fontsize = 18)
    ax.text(57.0, 1.1e7, r'$x=6.2 \cdot 10^{-4}$', fontsize = 18)
    ax.text(0.5e2, 5e6, r'$x=7.8 \cdot 10^{-4}$', fontsize = 18)
    ax.text(0.9e2, 2.7e6, r'$x=9.1 \cdot 10^{-4}$', fontsize = 18)
    ax.text(1e2, 1.15e6, r'$x=0.0011$', fontsize = 18)
    ax.text(1.9e2, 6e5, r'$x=0.0015$', fontsize = 18)
    ax.text(2.5e2, 2.9e5, r'$x=0.0022$', fontsize = 18)
    ax.text(3.7e2, 1.2e5, r'$x=0.0034$', fontsize = 18)
    ax.text(0.45e3, 5.6e4, r'$x=0.0046$', fontsize = 18)
    ax.text(0.8e3, 2.8e4, r'$x=0.0064$', fontsize = 18)
    ax.text(0.6e3, 1.2e4, r'$x=0.0094$', fontsize = 18)
    ax.text(0.8e3, 6200.0, r'$x=0.011$', fontsize = 18)
    ax.text(1.5e3, 3200.0, r'$x=0.0138$', fontsize = 18)
    ax.text(1e3, 1400.0, r'$x=0.0156$', fontsize = 18)
    ax.text(1e3, 700.0, r'$x=0.0183$', fontsize = 18)
    ax.text(1.9e3, 320.0, r'$x=0.020$', fontsize = 18)
    ax.text(1.1e3, 150.0, r'$x=0.025$', fontsize = 18)
    ax.text(2.5e3, 67.0, r'$x=0.032$', fontsize = 18)
    ax.text(60.0, 28.0, r'$x=0.05$', fontsize = 18)
    ax.text(60.0, 12.0, r'$x=0.08$', fontsize = 18)
    ax.text(60.0, 5.3, r'$x=0.13$', fontsize = 18)
    ax.text(60.0, 2.15, r'$x=0.18$', fontsize = 18)
    ax.text(160.0, 0.85, r'$x=0.25$', fontsize = 18)
    ax.text(60.0, 0.3, r'$x=0.40$', fontsize = 18)
    ax.text(70.0, 0.022, r'$x=0.65 \, (i=0)$', fontsize = 18)

    ax.text(4e3, 28.0, r'$x=0.05$', fontsize = 18, color = 'darkgreen')
    ax.text(6e3, 11.0, r'$x=0.08$', fontsize = 18, color = 'darkgreen')
    ax.text(1.6e3, 4.6, r'$x=0.13$', fontsize = 18, color = 'darkgreen')
    ax.text(2.3e3, 2.0, r'$x=0.18$', fontsize = 18, color = 'darkgreen')
    ax.text(1.1e3, 0.8, r'$x=0.25\, (i=2)$', fontsize = 18, color = 'darkgreen')

    ax.semilogy()
    ax.semilogx()
    ax.set_xlim(1.0, 4e4)
    ax.set_ylim(1e-2, 2e9)

    ax.set_ylabel(r'$\sigma_r^{p,NC}$'+r'$\, \times\, 2^{\, i}$', size = 24)
    ax.set_xlabel(r'$Q^2 \: \rm{(GeV^2)}$', size = 24)

    ax.set_xticks([1.0, 10.0, 1e2, 1e3, 1e4])
    ax.set_xticklabels([r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$'])
    locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks = 12)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.yaxis.set_tick_params(which = 'minor', length = 3)
    ax.set_yticks([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9])
    ax.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$', r'$10^7$', r'$10^8$', r'$10^9$'])

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = 'on', labelright = 'off')
    ax.tick_params(axis = 'x', which = 'both', labeltop = 'off', labelbottom = 'on')

    ax.xaxis.set_tick_params(which = 'major', length = 6)
    ax.xaxis.set_tick_params(which = 'minor', length = 3)
    ax.yaxis.set_tick_params(which = 'major', length = 6)
    ax.yaxis.set_tick_params(which = 'minor', length = 3)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = 'on', labelright = 'off')
    ax.tick_params(axis = 'x', which = 'both', labeltop = 'off', labelbottom = 'on')

    ax.tick_params(axis='both',labelsize=22)
    ax.legend([HERA_10027, HERA_10028, HERA_10029], \
              [r'$\mathrm{HERA}~\sqrt{s}=300\, \mathrm{GeV}$', \
               r'$\mathrm{HERA}~\sqrt{s}=251\, \mathrm{GeV}$', \
               r'$\mathrm{HERA}~\sqrt{s}=225\, \mathrm{GeV}$'], fontsize = 24, frameon = 0)

    py.tight_layout()

    py.savefig('%s/gallery/dis-hera-other-%d.png' % (wdir, istep), dpi = dpi)
    py.close()

def plot_obs(wdir, plot_HERA = True, only_best_cluster = True, dpi = 200):

    print('\nplotting dis data from %s' % (wdir))

    load_config('%s/input.py' % wdir)
    istep = core.get_istep()
    replicas = core.get_replicas(wdir)
    core.mod_conf(istep, replicas[0]) #--set conf as specified in istep

    predictions = load('%s/data/predictions-%d.dat' % (wdir, istep))
    if 'idis' not in predictions['reactions']:
        print('inclusive DIS is not in data file')
        return
    labels  = load('%s/data/labels-%d.dat' % (wdir, istep))
    cluster = labels['cluster']
    cluster_average_chi2 = labels['cluster_average_chi2']
    best_cluster = np.argmin(cluster_average_chi2)

    data = predictions['reactions']['idis']
    if data.keys() == []: return

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

        if 'X' in data[idx]: data[idx]['x'] = data[idx]['X']
        data[idx]['rQ2'] = np.around(data[idx]['Q2'], decimals = 0)
        data[idx]['rx'] = np.around(data[idx]['x'], decimals = 2)

    plot_proton(wdir, data, istep, dpi)
    plot_deuteron(wdir, data, istep, dpi)
    if plot_HERA:
        plot_hera_318(wdir, data, istep, dpi)
        plot_hera_other_than_318(wdir, data, istep, dpi)

    return
