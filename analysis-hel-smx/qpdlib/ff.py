#!/usr/bin/env python
import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT

#--matplotlib
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#matplotlib.rc('text',usetex=True)
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import pylab as py
import matplotlib.gridspec as gridspec

#--from scipy stack 
from scipy.integrate import quad

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core
from analysis.corelib import classifier

def gen_xf(wdir, had, flavors = ['g', 'u', 'ub', 'd', 'db', 's', 'sb'], Q2 = None):
    
    fflabel='ff%s'%had
    print('\ngenerating ff-%s from %s' % (had,wdir))

    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) ## set conf as specified in istep

    if fflabel not in conf['steps'][istep]['active distributions']:
        print('ff-%s not in active distribution' % had)
        return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman

    jar = load('%s/data/jar-%d.dat' % (wdir, istep))
    ## 'jar' contains parameters and their orders
    replicas = jar['replicas']
    ## 'jar['replicas']' is a list, in which each element is a list of parameters obeying 'jar['order']'

    ff = conf[fflabel]

    ## setup kinematics
    X = np.linspace(0.01, 0.99, 100)
    if Q2 == None: Q2 = conf['Q20']

    ## compute XF for all replicas
    XF = {}
    cnt = 0
    for par in replicas:
        cnt += 1
        lprint('%d/%d' % (cnt, len(replicas)))

        ## filter
        #flag=False
        #params=replica['params'][istep]
        #order=replica['order'][istep]
        #for i in range(len(order)):
        #    if order[i][0]!=1:continue
        #    if order[i][1]!='pdf':continue
        #    #if order[i][2]=='s1 a':
        #    #   if params[i]<-0.9: flag=True
        #if flag: continue

        parman.set_new_params(par, initial = True)

        #print
        #print conf['ffpion'].get_xF(0.5,10.0,'u')
        #print conf['ffkaon'].get_xF(0.5,10.0,'u')
        #print ff.get_xF(0.5,10.0,'u')

        for flavor in flavors:
            if flavor not in XF: XF[flavor] = []
            if   flavor=='c' or flavor=='cb' or flavor=='c+cb': _Q2=conf['aux'].mc2+1
            elif flavor=='b' or flavor=='bb' or flavor=='b+bb': _Q2=conf['aux'].mb2+1
            else:           _Q2=Q2
            if  flavor=='u+ub':
                func=lambda x: ff.get_xF(x,_Q2,'u')+ff.get_xF(x,_Q2,'ub')
            elif flavor=='d+db':
                func=lambda x: ff.get_xF(x,_Q2,'d')+ff.get_xF(x,_Q2,'db')
            elif flavor=='s+sb':
                func=lambda x: ff.get_xF(x,_Q2,'s')+ff.get_xF(x,_Q2,'cb')
            elif flavor=='c+cb':
                func=lambda x: ff.get_xF(x,_Q2,'c')+ff.get_xF(x,_Q2,'cb')
            elif flavor=='b+bb':
                func=lambda x: ff.get_xF(x,_Q2,'b')+ff.get_xF(x,_Q2,'bb')
            else:
                func=lambda x: ff.get_xF(x,_Q2,flavor)

            XF[flavor].append([func(x) for x in X])
    #print func(0.5)
    print
    checkdir('%s/data' % wdir)
    if Q2 == conf['Q20']:
        save({'X': X, 'Q2': Q2, 'XF': XF}, '%s/data/ff%s-%d.dat' % (wdir, had, istep))
    else:
        save({'X': X, 'Q2': Q2, 'XF': XF}, '%s/data/ff%s-%d-%d.dat' % (wdir, had, istep, int(Q2)))


