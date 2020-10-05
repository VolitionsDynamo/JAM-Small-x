#!/usr/bin/env python
import os,sys
import subprocess

#--matplotlib
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('text',usetex=True)
import pylab as py


#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

def get_msr_inspected(wdir):
    print('\nget msr inspected (filtered msr files) using %s\n'%wdir)

    replicas=sorted(os.listdir('%s/msr'%wdir))
    
    X=[]
    for i in range(len(replicas)):
        lprint('progress: %d/%d'%(i+1,len(replicas)))
        replica=load('%s/msr/%s'%(wdir,replicas[i]))
        istep=sorted(replica['params'].keys())[-1] #--pick last step
        params=replica['params'][istep]
        if params is None: continue
        data=replica['chi2']#[istep]
        chi2,npts=0,0
        for reaction in data:
            #print 
            for idx in  data[reaction]:
                chi2+=data[reaction][idx]['chi2']
                npts+=data[reaction][idx]['npts']
        if chi2/npts-1<3:
            X.append(chi2/npts-1)
            checkdir('%s/msr-inspected'%wdir)
            cmd=['cp']
            cmd.append('%s/msr/%s'%(wdir,replicas[i]))
            cmd.append('%s/msr-inspected/%s'%(wdir,replicas[i]))
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    print()
    print('original  num. samples :%d'%len(replicas))
    print('inspected num. samples :%d'%len(X))

    nrows=1
    ncols=1

    #--plot labeled residuals
    ax=py.subplot(nrows,ncols,1)
    ax.hist(X,bins=10)

    #py.tight_layout()
    checkdir('%s/gallery'%wdir)
    py.savefig('%s/gallery/chi2-dist-inspect.pdf'%(wdir))
    py.close()


