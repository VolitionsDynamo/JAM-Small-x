import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import pylab as py


#--from tools
from tools.tools     import checkdir,save,load
import tools.config
from tools.config    import load_config, conf, options
from tools.inputmod  import INPUTMOD

#--from local
from analysis.corelib import core
from analysis.corelib import classifier


def plot_params(wdir,dist,kc):
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas = core.get_replicas(wdir) 
    core.mod_conf(istep,replicas[0])

    clusters,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc)

    _order = replicas[0]['order'][istep]

    #--get correct order from dist
    order  = []
    idx    = []
    for i in range(len(_order)):
        if _order[i][1] != dist: continue
        order.append(_order[i][2])
        idx.append(i)

    #--get correct params from dist
    params = np.zeros((len(order),len(replicas)))
    for i in range(len(order)):
        for j in range(len(replicas)):
            params[i][j] = replicas[j]['params'][istep][idx[i]]

    #--create plot with enough space for # of parameters
    nrows, ncols = np.ceil(len(order)/5.0), 5
    fig = py.figure(figsize=(ncols*7,nrows*4))
    X = np.linspace(1,len(replicas),len(replicas))

    #--create plot
    for i in range(len(order)):
        ax = py.subplot(nrows,ncols, i+1)
        ax.set_title('%s'%(order[i]), size=20)
        #params = [replicas[j]['params'][istep][i] for j in range(len(replicas))]
        color  = [colors[clusters[j]] for j in range(len(replicas))]
        ax.scatter(X,params[i],color=color) 
        ax.plot(X,np.ones(len(X))*np.average(params[i]),'k--',alpha=0.5)

    filename='%s/gallery/%s-params.png'%(wdir,dist)
    checkdir('%s/gallery'%wdir)
    py.savefig(filename)
    print 'Saving figure to %s'%filename
     



