#!/usr/bin/env python
import os,sys
import numpy as np
import scipy as sp
import pandas as pd
import copy

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config
from tools.inputmod  import INPUTMOD
from tools.randomstr import id_generator


def set_passive_distributions(conf,istep):

    #--needs to be revises NS(09/02/19)

    print 'setting passive for istep ',istep
    step=conf['steps'][istep]
    for dist in step['passive distributions']:
        for par in conf['params'][dist]:
            #print dist, par
            conf['params'][dist][par]['fixed']=True
    return conf

def fix_parameters(conf,istep):

    #--needs to be revises NS(09/02/19)

    step=conf['steps'][istep]
    for dist in step['fix parameters']:
        for par in step['fix parameters'][dist]:
            conf['params'][dist][par]['fixed']=True
    return conf

def set_passive_params(istep,prior):

    #--needs to be revises NS(09/02/19)

    step=conf['steps'][istep]
    if 'passive distributions' not in step: return 
    for dist in step['passive distributions']:
        for par in conf['params'][dist]:
            #conf['params'][dist][par]['fixed']=True
            for idep in step['dep']:
                for i in range(len(prior['order'][idep])):
                    _,_dist,_par = prior['order'][idep][i]
                    if  dist==_dist and par==_par:
                        conf['params'][dist][par]['value']=prior['params'][idep][i]

        if dist=='pdf'    : parman.set_pdf_params()
        if dist=='ppdf'   : parman.set_ppdf_params()
        if dist=='ffpion' : parman.set_ffpion_params()
        if dist=='ffkaon' : parman.set_ffkaon_params()

def mod_conf(istep,replica=None):

    step=conf['steps'][istep]

    #--remove pdf/ff that is not in the step
    distributions=conf['params'].keys()  #--pdf,ppdf,ffpion,ffkaon,...
    for dist in distributions:
        if  dist in step['active distributions']:  
            continue
        else:
            del conf['params'][dist] 

    if replica!=None:
        #--set fixed==True for passive distributions
        if 'passive distributions' in step:
            for dist in step['passive distributions']:
                for par in conf['params'][dist]:

                    if conf['params'][dist][par]['fixed']==False:
                        conf['params'][dist][par]['fixed']=True

                    #--set prior parameters values for passive distributions
                    for istep in step['dep']:
                        prior_order=replica['order'][istep]
                        prior_params=replica['params'][istep]
                        for i in range(len(prior_order)):
                            _,_dist,_par = prior_order[i]
                            if  dist==_dist and par==_par:
                                conf['params'][dist][par]['value']=prior_params[i]

        #--another version for fixed parameters 
        if 'fix parameters' in step:
            for dist in step['fix parameters']:
                for par in step['fix parameters'][dist]:
                    conf['params'][dist][par]['fixed']=True
                    #--set prior parameters values for passive distributions
                    for istep in step['dep']:
                        prior_order=replica['order'][istep]
                        prior_params=replica['params'][istep]
                        for i in range(len(prior_order)):
                            _,_dist,_par = prior_order[i]
                            if  dist==_dist and par==_par:
                                conf['params'][dist][par]['value']=prior_params[i]
                        
    #--remove datasets not in the step
    datasets=conf['datasets'].keys() #--idis,dy,....
    for dataset in datasets:
        if  dataset in step['datasets']:  

            #--remove entry from xlsx
            xlsx=conf['datasets'][dataset]['xlsx'].keys()
            for idx in xlsx:
                if  idx in step['datasets'][dataset]:
                    continue
                else:
                    del conf['datasets'][dataset]['xlsx'][idx]

            #--remove entry from norm
            norm=conf['datasets'][dataset]['norm'].keys()
            for idx in norm:
                if  idx in step['datasets'][dataset]:
                    continue
                else:
                    del conf['datasets'][dataset]['norm'][idx]
        else:
            del conf['datasets'][dataset]                

def get_replicas(wdir):
    """
    load the msr files
    """
    replicas=sorted(os.listdir('%s/msr-inspected'%wdir))
    replicas=[load('%s/msr-inspected/%s'%(wdir,_)) for _ in replicas]
    return replicas

def get_replicas_names(wdir):
    replicas=sorted(os.listdir('%s/msr-inspected'%wdir))
    return replicas



def get_istep():
    #--pick last step
    return sorted(conf['steps'])[-1] 



