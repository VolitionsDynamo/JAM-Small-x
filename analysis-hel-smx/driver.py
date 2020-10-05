!/usr/bin/env python
import os, sys
import argparse
import re

import kmeanconf as kc
import all_sim

from analysis_smx.corelib import core
from analysis_smx.corelib import inspect
from analysis_smx.corelib import predict
from analysis_smx.corelib import classifier
from analysis_smx.corelib import jar
from analysis_smx.corelib import mlsamples
from analysis_smx.corelib import summary

from analysis_smx.qpdlib  import tolhapdf
from analysis_smx.qpdlib  import benchmark
from analysis_smx.qpdlib  import ppdf_proton
from analysis_smx.qpdlib  import ppdf_proton_plot
from analysis_smx.qpdlib  import ppdf_smx_proton
from analysis_smx.qpdlib  import ppdf_smx_proton_plot

from analysis_smx.parlib import params

from analysis_smx.obslib  import idis, pidis, pidis_smx

## use as
## ./driver.py results1/step01

working_directory = sys.argv[1]
dpi = 200

print working_directory

#####################
## Simulation for ALL
#####################
#--target
tar = 'p'
#--estimation of errors
#--'opt': optimistic
#--'mod': moderate
#--'pes': pessimistic
est = 'opt'
#--luminosity (keep as 100 fb-1)
lum = '100:fb-1'
#--force predictions to be regenerated
force = True

FILT=[]
#FILT.append(('g1 N', 0.0, 'less'))

## initial process
inspect.get_msr_inspected(working_directory) #,limit=2.0,FILT=FILT)
#all_sim.ALL(working_directory,tar=tar,est=est,lum=lum,force=force)

predict.get_predictions(working_directory, cores = 2)
classifier.gen_labels(working_directory, kc)
jar.gen_jar_file(working_directory, kc)
summary.print_summary(working_directory, kc)

## PDFs
## plot PDF distributions with lines
#flavors = ['g', 'uv', 'dv', 'd/u', 'db+ub', 'db-ub', 's+sb', 'rs']
#pdf_proton.gen_xf(working_directory, flavors)
#pdf_proton_plot.PLOTS(1, working_directory, dpi = dpi)
## compare PDF distributions with other groups
#groups = ['JAM19PDF_proton_nlo', 'NNPDF31_nlo_as_0118', 'MMHT2014nlo68cl', 'CSKK_nnlo_EIG']
#pdf_proton.gen_xf(working_directory, flavors, 4.0)
#lhapdf_data.generate_xf(working_directory, groups, flavors, 4.0) ## after you upate 'groups', this only has to be run once
#pdf_proton_plot.PLOTS(2, working_directory, Q2 = 4.0, dpi = dpi)

## polarized PDFs
flavors = ['g','up', 'dp', 'sp']
#ppdf_proton.gen_xf(working_directory, flavors)
#ppdf_proton_plot.PLOTS(5, working_directory, dpi = dpi)
#ppdf_proton_plot.PLOTS(1, working_directory, dpi = dpi)
ppdf_smx_proton.gen_xf(working_directory, flavors,Q2=10)
ppdf_smx_proton_plot.PLOTS(5, working_directory,Q2=10,dpi = dpi)  #5 for ppdf with error band
ppdf_smx_proton_plot.PLOTS(1, working_directory,Q2=10,dpi = dpi)  #1 for ppdf all replicas plotted
#ppdf_smx_proton_plot.PLOTS(6, working_directory,Q2=10,dpi = dpi)  #6 for g1 with error band

## others related to PDFs
#pdf_and_ppdf_proton.get_pdf_spin_breakdown(working_directory, ['g', 'u', 'd', 's'])
#pdf_and_ppdf_proton_plot.RATIO(1, working_directory)
#pdf_and_ppdf_proton_plot.SPIN(1, working_directory)
#pdf_and_ppdf_proton_plot.SPIN(2, working_directory)
#params.plot_params(working_directory, 'pdf', kc)
#params.plot_params(working_directory, 'ppdf', kc)
params.plot_params(working_directory, 'ppdf_smx', kc)

## observables DIS, DY, JET
#idis.plot_obs(working_directory, dpi)
#dy_hadron.plot_obs(working_directory, dpi)
#jet.make_figure(working_directory, 2, plot_with_factor = True, replica_lines = True, dpi = dpi)
#pidis.make_figure(working_directory, 1, only_best_cluster = False, dpi = dpi)
pidis_smx.make_figure(working_directory, 1, only_best_cluster = False, dpi = dpi)
#pjet.make_figure(working_directory, 1, plot_with_factor = False, replica_lines = True, dpi = dpi)
