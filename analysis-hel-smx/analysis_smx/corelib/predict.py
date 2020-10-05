#!/usr/bin/env python
import os, sys
import subprocess
import numpy as np
import scipy as sp
import pandas as pd
import copy

#--from tools
from tools           import config
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config
from tools.inputmod  import INPUTMOD
from tools.randomstr import id_generator

#--from fitlib
from fitlib_smx.resman import RESMAN

#--from local
import core

def ask_replica_mismatch():
    if sys.version_info[0] == 2:
        while True:
            flag = raw_input('o (overwrite)/e (exit): ')
            if flag == '':
                print('please type "o" or "e"')
            elif flag.lower()[0] == 'o':
                part_flag = False
                exit_flag = False
                break
            elif flag.lower()[0] == 'e':
                part_flag = False
                exit_flag = True
                break
            else:
                print('please type "o" or "e"')
    elif sys.version_info[0] == 3:
        while True:
            flag = input('y/n: ')
            if flag == '':
                print('please type "o" or "e"')
            elif flag.lower()[0] == 'o':
                part_flag = False
                exit_flag = False
                break
            elif flag.lower()[0] == 'e':
                part_flag = False
                exit_flag = True
                break
            else:
                print('please type "o" or "e"')

    return part_flag, exit_flag

def get_predictions(wdir, cores = 2, mod_conf = None, pre_flags = {}):
    ## 'pre_flags' can be defined for example like the following
    ## pre_flags = {'regenerate_same': False, 'regenerate_part': False, 'regenerate_mismatch': True}
    if mod_conf == None:
        load_config('%s/input.py' % wdir)
    else:
        config.conf = conf
    conf['bootstrap'] = False
    istep = core.get_istep()

    current_replica_names = sorted(os.listdir('%s/msr-inspected' % wdir))
    previous_replica_names_path = '%s/data/predicted-replicas-%d.dat' % (wdir, istep)
    previous_predictions_path = '%s/data/predictions-%d.dat' % (wdir, istep)
    part_flag = False
    if os.path.exists(previous_predictions_path) and os.path.exists(previous_replica_names_path):
        previous_replica_names = load(previous_replica_names_path)
        if current_replica_names == previous_replica_names:
            print('The replicas names are the same since predictions were generated last time,')
            print('do you wish to regenrate and overwrite previous predictions?')
            if 'regenerate_same' in pre_flags:
                if pre_flags['regenerate_same']:
                    print('regenerating everything based on predifined flag...')
                    part_flag = False
                else:
                    print('exit generating predictions based on predifined flag.')
                    return
            else:
                if sys.version_info[0] == 2:
                    while True:
                        same_flag = raw_input('y (regenerate)/n (skip generation of predictions): ')
                        if same_flag == '':
                            print('please type "y" or "n"')
                        elif same_flag.lower()[0] == 'n':
                            return
                        elif same_flag.lower()[0] == 'y':
                            part_flag = False
                            break
                        else:
                            print('please type "y" or "n"')
                elif sys.version_info[0] == 3:
                    while True:
                        same_flag = input('y (regenerate)/n (skip generation of predictions): ')
                        if same_flag == '':
                            print('please type "y" or "n"')
                        elif same_flag.lower()[0] == 'n':
                            return
                        elif same_flag.lower()[0] == 'y':
                            part_flag = False
                            break
                        else:
                            print('please type "y" or "n"')
        elif all(_ in current_replica_names for _ in previous_replica_names):
            print('Based on the replica names,')
            print('part of the predictions you are trying to generate has already been generated last time,')
            print('do you wish to regenrate with all replicas or only generate with new replicas?')
            if 'regenerate_part' in pre_flags:
                if pre_flags['regenerate_part']:
                    print('regenerating with all replicas based on predifined flag...')
                    part_flag = False
                else:
                    print('generating only on new replicas based on predifined flag...')
                    part_flag = True
            else:
                if sys.version_info[0] == 2:
                    while True:
                        part_flag = raw_input('w (regenerate with all replicas)/p (generate only with new replicas): ')
                        if part_flag == '':
                            print('please type "w" or "p"')
                        elif part_flag.lower()[0] == 'w':
                            part_flag = False
                            break
                        elif part_flag.lower()[0] == 'p':
                            part_flag = True
                            break
                        else:
                            print('please type "y" or "n"')
                elif sys.version_info[0] == 3:
                    while True:
                        part_flag = input('w (regenerate with all replicas)/p (generate only with new replicas): ')
                        if part_flag == '':
                            print('please type "w" or "p"')
                        elif part_flag.lower()[0] == 'w':
                            part_flag = False
                            break
                        elif part_flag.lower()[0] == 'p':
                            part_flag = True
                            break
                        else:
                            print('please type "y" or "n"')
        else:
            part_flag = False
    else:
        pass

    if part_flag:
        replicas_to_load = np.setdiff1d(current_replica_names, previous_replica_names)
        replicas = [load('%s/msr-inspected/%s' % (wdir, _)) for _ in replicas_to_load]
    else:
        replicas = core.get_replicas(wdir)

    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep

    resman=RESMAN(nworkers=cores,parallel=True,datasets=True)
    parman=resman.parman
    order=replicas[0]['order'][istep]
    ## order and parameters are updated later for each replicas individually
    # parman.order = order
    # order = parman.order

    if part_flag:
        previous_data = load(previous_predictions_path)
        if order != previous_data['order']:
            print('Previous replicas have different parameter orders from new replicas,')
            print('previously: ', previous_data['order'])
            print('currently: ', order)
            print('do you wish to overwrite the previous predictions or exit to check?')
            if 'regenerate_mismatch' in pre_flags:
                if pre_flags['regenerate_mismatch']:
                    print('regenerating everything based on predifined flag...')
                    part_flag = False
                    exit_flag = False
                else:
                    print('exit generating predictions based on predifined flag')
                    exit_flag = True
            else:
                part_flag, exit_flag = ask_replica_mismatch()
            if exit_flag: sys.exit()

    obsres={}
    if 'idis'     in conf['datasets'] : obsres['idis']     = resman.idis_res
    if 'pidis'    in conf['datasets'] : obsres['pidis']    = resman.pidis_res
    if 'sidis'    in conf['datasets'] : obsres['sidis']    = resman.sidis_res
    if 'psidis'   in conf['datasets'] : obsres['psidis']   = resman.psidis_res
    if 'dy'       in conf['datasets'] : obsres['dy']       = resman.dy_res
    if 'sia'      in conf['datasets'] : obsres['sia']      = resman.sia_res
    if 'qpdf'     in conf['datasets'] : obsres['qpdf']     = resman.qpdf_res
    if 'dy-pion'  in conf['datasets'] : obsres['dy-pion']  = resman.dy_pion_res
    if 'pion_qT'  in conf['datasets'] : obsres['pion_qT']  = resman.pion_qTres
    if 'ln'       in conf['datasets'] : obsres['ln']       = resman.ln_res
    if 'jet'      in conf['datasets'] : obsres['jet']      = resman.jet_res
    if 'pjet'     in conf['datasets'] : obsres['pjet']     = resman.pjet_res

    if part_flag:
        # previous_data = load(previous_predictions_path)
        if obsres.keys() != previous_data['reactions'].keys():
            print('Previous replicas have different observables from new replicas,')
            print('previously: ', previous_data['reactions'].keys())
            print('currently: ', obsres.keys())
            print('do you wish to overwrite the previous predictions or exit to check?')
            if 'regenerate_mismatch' in pre_flags:
                if pre_flags['regenerate_mismatch']:
                    print('regenerating everything based on predifined flag...')
                    part_flag = False
                    exit_flag = False
                else:
                    print('exit generating predictions based on predifined flag')
                    exit_flag = True
            else:
                part_flag, exit_flag = ask_replica_mismatch()
            if exit_flag: sys.exit()
        if part_flag:
            for _ in obsres.keys():
                if obsres[_].tabs.keys() != previous_data['reactions'][_].keys():
                    print('Previous replicas have different datasets from new replicas in observable %s,' % _)
                    print('previously: ', previous_data['reactions'][_].keys())
                    print('currently: ', obsres[_].tabs.keys())
                    print('do you wish to overwrite the previous predictions or exit to check?')
                    if 'regenerate_mismatch' in pre_flags:
                        if pre_flags['regenerate_mismatch']:
                            print('regenerating everything based on predifined flag...')
                            part_flag = False
                            exit_flag = False
                        else:
                            print('exit generating predictions based on predifined flag')
                            exit_flag = True
                    else:
                        part_flag, exit_flag = ask_replica_mismatch()
                    if exit_flag: sys.exit()

    #--setup big table to store all we want
    data={}
    data['order']=order
    data['params']=[]
    data['reactions']={}
    data['res']=[]
    data['rres']=[]
    data['nres']=[]

    for _ in obsres:
        tabs=copy.copy(obsres[_].tabs)
        #--create a space to store all the predictions from replicas
        for idx in tabs:
            tabs[idx]['prediction-rep']=[]
            tabs[idx]['residuals-rep']=[]
        data['reactions'][_]=tabs

    print('\ngen predictions using %s\n'%wdir)

    cnt=0
    if part_flag:
        n_all_replicas = len(current_replica_names)
        n_available_replicas = len(previous_replica_names)
        cnt += n_available_replicas
    else:
        replicas = core.get_replicas(wdir)

    for replica in replicas:
        cnt+=1
        if part_flag:
            lprint('progress: %d/%d'%(cnt, n_all_replicas))
        else:
            lprint('progress: %d/%d'%(cnt, len(replicas)))

        core.mod_conf(istep,replica) #--will update passive dist

        parman.par=copy.copy(replica['params'][istep])
        parman.order=copy.copy(replica['order'][istep])
        data['params']=np.append(data['params'],parman.par)

        #for i in range(len(parman.par)):
        #    print parman.order[i],parman.par[i],replica['order'][istep][i]

        #--compute residuals (==theory)
        res,rres,nres=resman.get_residuals(parman.par, initial = True)
        data['res'].append(res)
        data['rres'].append(rres)
        data['nres'].append(nres)

        #--save predictions of the current step and current replica at data
        for _ in obsres:
            for idx in data['reactions'][_]:
                prediction=copy.copy(obsres[_].tabs[idx]['prediction'])
                residuals=copy.copy(obsres[_].tabs[idx]['residuals'])
                data['reactions'][_][idx]['prediction-rep'].append(prediction)
                data['reactions'][_][idx]['residuals-rep'].append(residuals)
    print

    #--close resman
    resman.shutdown()

    ## append newly generated parts if 'part_flag' is True
    if part_flag:
        data['params'] = np.append(previous_data['params'], data['params'])
        data['res'] = np.concatenate((previous_data['res'], data['res']))
        data['rres'] = np.concatenate((previous_data['rres'], data['rres']))
        data['nres'] = np.concatenate((previous_data['nres'], data['nres']))

        for _ in data['reactions']:
            for dataset in data['reactions'][_]:
                data['reactions'][_][dataset]['prediction-rep'] = previous_data['reactions'][_][dataset]['prediction-rep'] + data['reactions'][_][dataset]['prediction-rep']
                data['reactions'][_][dataset]['residuals-rep'] = previous_data['reactions'][_][dataset]['prediction-rep'] + data['reactions'][_][dataset]['residuals-rep']

    #--convert tables to numpy array before saving
    for _ in ['res','rres','nres']:
        data[_]=np.array(data[_])

    checkdir('%s/data'%wdir)
    if mod_conf==None:
        save(data,'%s/data/predictions-%d.dat'%(wdir,istep))
        save(current_replica_names, '%s/data/predicted-replicas-%d.dat' % (wdir, istep))
    else:
        save(data,'%s/data/predictions-%d-sim.dat'%(wdir,istep))

def get_summary(self,wdir,istep):

    #--needs revision NS (09/02/19)

    data=load('%s/data/predictions-%d.dat'%(wdir,istep))

    summary=[]
    dic={}
    global_chi2=0
    global_npts=0
    for reaction in data['reactions']:
        dic[reaction]={}
        for idx in data['reactions'][reaction]:
            dic[reaction][idx]={}
            tab=data['reactions'][reaction][idx]
            value=tab['value']
            alpha=tab['alpha']
            CHI2=[]
            for i in range(len(tab['predictions'])):
                chi2tot=np.sum(data['res'][i]**2)/data['res'][i].size
                #if chi2tot>2: continue
                prediction=tab['predictions'][i]
                CHI2.append(np.sum(((value-prediction)/alpha)**2))
            npts=len(value)
            col=tab['col'][0]
            chi2min=np.amin(CHI2)/npts
            chi2max=np.amax(CHI2)/npts
            chi2ave=np.mean(CHI2)/npts

            global_chi2+=np.amin(CHI2)
            global_npts+=npts

            msg ='reaction: %8s'
            msg+=' idx: %7d'
            msg+=' col: %10s'
            msg+=' chi2/npts (min,ave,max): %6.2f %6.2f %6.2f'
            msg+=' npts: %5d'
            msg=msg%(reaction,idx,col[:10],chi2min,chi2ave,chi2max,npts)
            dic[reaction][idx]['col']=col
            dic[reaction][idx]['chi2']=np.amin(CHI2)
            dic[reaction][idx]['chi2/npts']=np.amin(CHI2)/npts
            dic[reaction][idx]['npts']=npts
            if  reaction=='dy':
                msg+=' rea: %s'%tab['reaction'][0]
            if  reaction=='sia':
                msg+=' had: %s'%tab['hadron'][0]
            if  reaction=='idis':
                msg+=' tar: %s'%tab['target'][0]
                msg+=' obs: %s'%tab['obs'][0]
            if  reaction=='sidis':
                msg+=' tar: %1s'%tab['target'][0][0]
                msg+=' had: %5s'%tab['hadron'][0]
                msg+=' obs: %s'%tab['obs'][0]
            if  reaction=='pidis':
                msg+=' tar: %s'%tab['target'][0]
                msg+=' obs: %s'%tab['obs'][0]
            if  reaction=='psidis':
                msg+=' tar: %1s'%tab['target'][0][0]
                msg+=' had: %5s'%tab['hadron'][0]
                msg+=' obs: %s'%tab['obs'][0]
            if  reaction == 'jet':
                msg += ' rea: %s' % tab['reaction'][0]
                msg += ' col: %s' % tab['col'][0]
            if  reaction == 'pjet':
                msg += ' rea: %s' % tab['reaction'][0]
                msg += ' col: %s' % tab['col'][0]
            #print msg
            summary.append(msg)
        #print tab.keys()

    print 'global summary'
    print 'chi2/npts =',global_chi2/global_npts
    print 'npts      =',global_npts

    summary=[_+'\n' for _ in summary]
    F=open('%s/data/summary-%d.txt'%(wdir,istep),'w')
    F.writelines(summary)
    F.close()
    save(dic,'%s/data/summary-%d.dat'%(wdir,istep))

def gen_exp_dict(self,wdir,istep):

    #--needs revision NS (09/02/19)

    data=load('%s/data/predictions-%d.dat'%(wdir,istep))
    summary=[]
    for reaction in data['reactions']:
        for idx in data['reactions'][reaction]:
            tab=data['reactions'][reaction][idx]
            col=tab['col'][0]
            value=tab['value']
            npts=len(value)

            msg ='reaction: %8s,'
            msg+=' idx: %7d,'
            msg+=' col: %10s,'
            msg+=' npts: %5d,'
            msg=msg%(reaction,idx,col[:10],npts)
            if  reaction=='dy':
                msg+=' rea: %s,'%tab['reaction'][0]
            if  reaction=='sia':
                msg+=' had: %s,'%tab['hadron'][0]
            if  reaction=='idis':
                msg+=' tar: %s,'%tab['target'][0]
                msg+=' obs: %s,'%tab['obs'][0]
            if  reaction=='sidis':
                msg+=' tar: %1s,'%tab['target'][0][0]
                msg+=' had: %5s,'%tab['hadron'][0]
                msg+=' obs: %s,'%tab['obs'][0]
            if  reaction=='pidis':
                msg+=' tar: %s,'%tab['target'][0]
                msg+=' obs: %s,'%tab['obs'][0]
            if  reaction=='psidis':
                msg+=' tar: %1s,'%tab['target'][0][0]
                msg+=' had: %5s,'%tab['hadron'][0]
                msg+=' obs: %s,'%tab['obs'][0]
            if  reaction == 'jet':
                msg += ' rea: %s' % tab['reaction'][0]
                msg += ' col: %s' % tab['col'][0]
            if  reaction == 'pjet':
                msg += ' rea: %s' % tab['reaction'][0]
                msg += ' col: %s' % tab['col'][0]
            summary.append(msg)

    summary=[_+'\n' for _ in summary]
    F=open('%s/data/exp-dict-%d.txt'%(wdir,istep),'w')
    F.writelines(summary)
    F.close()






