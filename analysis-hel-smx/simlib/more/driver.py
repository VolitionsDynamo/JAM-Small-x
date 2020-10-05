#!/usr/bin/env python
import os,sys

from analysis.corelib import core
from analysis.corelib import inspect
from analysis.corelib import predict
from analysis.corelib import classifier
from analysis.corelib import optpriors
from analysis.corelib import jar
from analysis.corelib import mlsamples
from analysis.corelib import summary
from analysis.qpdlib  import tolhapdf
from analysis.qpdlib  import benchmark
from analysis.obslib  import idis

from analysis.simlib import sim


wdir=sys.argv[1]

#--gen pseudo data
herakeys = [10026,10027,10028,10029,10030,10031,10032]
os.system('mv {:s}/input.py {:s}/orig_input.py'.format(wdir,wdir))
os.system('cp ./simul/input.py {:s}'.format(wdir))

EICkeys  = [30006,30007,30008,30009]                  #eic.run(reaction,tableid)
EICkeys2  = [230006,230007,230008,230009]                  #eic.run(reaction,tableid)
sim.gen_pseudo_template_data(herakeys)

predict.get_predictions(wdir)

#sim.gen_pseudo_from_predict(wdir,EICkeys)      ##--observables #dy_pion.plot_obs(wdir,kc)
#sim.gen_pseudo_template_EIC(EICkeys,'p','idis',20,20)  


#ln.plot_obs(wdir,kc)
#dy_pT_pion.plot_obs(wdir,kc)

#--initial analysis
#
#inspect.get_msr_inspected(wdir)
#classifier.gen_labels(wdir,kc)
#jar.gen_jar_file(wdir,kc)
#summary.print_summary(wdir,kc)

