#!/usr/bin/env python
import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core

def gen_cj_grid():

    Q2=[1.30000E+00,1.50159E+00,1.75516E+00,2.07810E+00\
            ,2.49495E+00,3.04086E+00,3.76715E+00,4.50000E+00\
            ,4.75000E+00,6.23113E+00,8.37423E+00,1.15549E+01\
            ,1.64076E+01,2.40380E+01,3.64361E+01,5.73145E+01\
            ,9.38707E+01,1.60654E+02,2.88438E+02,5.45587E+02\
            ,1.09231E+03,2.32646E+03,5.30043E+03,1.29956E+04\
            ,3.45140E+04,1.00000E+05]

    X=[1.00000E-06,1.28121E-06,1.64152E-06,2.10317E-06\
           ,2.69463E-06,3.45242E-06,4.42329E-06,5.66715E-06\
           ,7.26076E-06,9.30241E-06,1.19180E-05,1.52689E-05\
           ,1.95617E-05,2.50609E-05,3.21053E-05,4.11287E-05\
           ,5.26863E-05,6.74889E-05,8.64459E-05,1.10720E-04\
           ,1.41800E-04,1.81585E-04,2.32503E-04,2.97652E-04\
           ,3.80981E-04,4.87518E-04,6.26039E-04,8.00452E-04\
           ,1.02297E-03,1.30657E-03,1.66759E-03,2.12729E-03\
           ,2.71054E-03,3.44865E-03,4.37927E-03,5.54908E-03\
           ,7.01192E-03,8.83064E-03,1.10763E-02,1.38266E-02\
           ,1.71641E-02,2.11717E-02,2.59364E-02,3.15062E-02\
           ,3.79623E-02,4.53425E-02,5.36750E-02,6.29705E-02\
           ,7.32221E-02,8.44039E-02,9.64793E-02,1.09332E-01\
           ,1.23067E-01,1.37507E-01,1.52639E-01,1.68416E-01\
           ,1.84794E-01,2.01731E-01,2.19016E-01,2.36948E-01\
           ,2.55242E-01,2.73927E-01,2.92954E-01,3.12340E-01\
           ,3.32036E-01,3.52019E-01,3.72282E-01,3.92772E-01\
           ,4.13533E-01,4.34326E-01,4.55495E-01,4.76836E-01\
           ,4.98342E-01,5.20006E-01,5.41818E-01,5.63773E-01\
           ,5.85861E-01,6.08077E-01,6.30459E-01,6.52800E-01\
           ,6.75387E-01,6.98063E-01,7.20830E-01,7.43683E-01\
           ,7.66623E-01,7.89636E-01,8.12791E-01,8.35940E-01\
           ,8.59175E-01,8.82485E-01,9.05866E-01,9.29311E-01\
           ,9.52817E-01,9.76387E-01,1.00000E+00]

    return X,Q2

def gen_dss_grid():

    Q2=[1.3e0, 1.5e0, 2.5e0, 
             4.0e0, 6.4e0, 1.0e1, 1.5e1, 2.5e1, 4.0e1, 6.4e1,
             1.0e2, 1.8e2, 3.2e2, 5.8e2, 1.0e3, 1.8e3,
             3.2e3, 5.8e3, 1.0e4, 1.8e4, 3.2e4, 5.8e4, 1.0e5]
    #X=[0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
    #        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
    #        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5, 
    #        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
    #        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 
    #        0.925, 0.95, 0.975, 1.0] 
    X=np.linspace(0.01,1,200)

    return X,Q2

def gen_grid(dist):
    if    dist.startswith('pdf') : return gen_cj_grid()
    elif  dist.startswith('ppdf'): return gen_cj_grid()
    #elif  dist.startswith('ff')  : return gen_dss_grid()
    elif  dist.startswith('ff')  : return gen_cj_grid()

def _gen_table(dist):

    X,Q2=gen_grid(dist)
    qpd=conf[dist]

    nx=len(X)
    nQ2=len(Q2)

    #--fill table
    table={iflav:[]  for iflav in [-5,-4,-3,-2,-1,1,2,3,4,5,21]}  
    npts=nQ2*nx
    for iQ2 in range(nQ2):
        for ix in range(nx):
            table[21].append(qpd.get_xF(X[ix],Q2[iQ2],'g'))
            table[ 1].append(qpd.get_xF(X[ix],Q2[iQ2],'d'))
            table[ 2].append(qpd.get_xF(X[ix],Q2[iQ2],'u'))
            table[ 3].append(qpd.get_xF(X[ix],Q2[iQ2],'s'))
            table[ 4].append(qpd.get_xF(X[ix],Q2[iQ2],'c'))
            table[ 5].append(qpd.get_xF(X[ix],Q2[iQ2],'b'))
            table[-1].append(qpd.get_xF(X[ix],Q2[iQ2],'db'))
            table[-2].append(qpd.get_xF(X[ix],Q2[iQ2],'ub'))
            table[-3].append(qpd.get_xF(X[ix],Q2[iQ2],'sb'))
            table[-4].append(qpd.get_xF(X[ix],Q2[iQ2],'c'))
            table[-5].append(qpd.get_xF(X[ix],Q2[iQ2],'b'))

    #--remap tables to match with lhapdf format
    for iflav in [-5,-4,-3,-2,-1,1,2,3,4,5,21]: 
        new_list=[]
        for ix in range(nx):
            for inQ2 in range(nQ2):
                idx=ix+inQ2*nx
                new_list.append(table[iflav][idx])
        table[iflav]=new_list

    return X,Q2,table

def gen_lhapdf_dat_file(X,Q2,table, wdir, file_name,setlabel):
    lines=[]  
    lines.append('PdfType: replica')
    lines.append('Format: lhagrid1')
    lines.append('---')
    line=''
    for _ in X: line+=('%10.5e '%_).upper()
    lines.append(line)
    line=''
    for _ in Q2: line+=('%10.5e '%_**0.5).upper()
    lines.append(line)
    lines.append('-5 -4 -3 -2 -1 1 2 3 4 5 21')

    nx=len(X)
    nQ2=len(Q2)

    for i in range(nx*nQ2):
        line=''
        for iflav in [-5,-4,-3,-2,-1,1,2,3,4,5,21]:
            line+=('%10.5e '%table[iflav][i]).upper()
        lines.append(line)
    lines.append('---')
    lines=[l+'\n' for l in lines]
    idx=str(setlabel).zfill(4)
    tab=open('%s/data/%s/%s_%s.dat'%(wdir, file_name, file_name,idx),'w')
    tab.writelines(lines)
    tab.close()

def gen_lhapdf_info_file(X,Q2,nrep, wdir, file_name,info):

    aS=[conf['alphaS'].get_alphaS(_) for _ in Q2]
    mZ=conf['aux'].mZ
    mb=conf['aux'].mb
    mc=conf['aux'].mc
    alphaSMZ=conf['aux'].alphaSMZ
    xmin=X[0]
    xmax=X[-1]
    Qmin=Q2[0]**0.5
    Qmax=Q2[-1]**0.5


    lines=[]
    lines.append('SetDesc:         "<description>"')
    lines.append('SetIndex:        <index>')
    lines.append('Authors:         <authors>')
    lines.append('Reference:       <reference>')
    lines.append('Format:          lhagrid1')
    lines.append('DataVersion:     1')
    lines.append('NumMembers:      %d'%nrep)
    lines.append('Particle:        <particle>')
    lines.append('Flavors:         [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]')
    lines.append('OrderQCD:        1')
    lines.append('FlavorScheme:    variable')
    lines.append('NumFlavors:      5')
    lines.append('ErrorType:       replicas')
    lines.append('XMin:            %0.2e'%xmin)
    lines.append('XMax:            %0.2e'%xmax)
    lines.append('QMin:            %0.2e'%Qmin)
    lines.append('QMax:            %0.2e'%Qmax)
    lines.append('MZ:              %f'%mZ)
    lines.append('MUp:             0.0')
    lines.append('MDown:           0.0')
    lines.append('MStrange:        0.0')
    lines.append('MCharm:          %f'%mc)
    lines.append('MBottom:         %f'%mb)
    lines.append('MTop:            180.0')
    lines.append('AlphaS_MZ:       %f'%alphaSMZ)
    lines.append('AlphaS_OrderQCD: 1')
    lines.append('AlphaS_Type:     ipol')
    line='AlphaS_Qs: ['
    for _ in Q2: line+=('%10.5e, '%_**0.5).upper()
    line=line.rstrip(',')+']'
    lines.append(line)
    line='AlphaS_Vals: ['
    for _ in aS: line+=('%10.5e, '%_).upper()
    line=line.rstrip(',')+']'
    lines.append(line)
    lines.append('AlphaS_Lambda4: 0')
    lines.append('AlphaS_Lambda5: 0')

    for i in range(len(lines)):
        for _ in info:
            lines[i]=lines[i].replace(_,info[_])


    lines=[l+'\n' for l in lines]
    tab=open('%s/data/%s/%s.info'%(wdir, file_name, file_name),'w')
    tab.writelines(lines)
    tab.close()

def gen_tables(wdir,dist, file_name,info,info_only=False):

    print('\ngenerating LHAPDF tables for %s using %s'%(dist,wdir))

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    core.mod_conf(istep) #--set conf as specified in istep   

    resman=RESMAN(nworkers=1,parallel=False,datasets=False)
    parman=resman.parman

    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas=jar['replicas']

    #--check order consistency
    order=jar['order']
    parman.order = order

    #--create output dir
    checkdir(wdir + '/data/')
    checkdir(wdir + '/data/%s/' % file_name)

    #--gen lhapdf_data_files
    if info_only==False:
        cnt=0
        for par in replicas:
            lprint('progress: %d/%d'%(cnt+1,len(replicas)))
            parman.set_new_params(par)
            X,Q2,table = _gen_table(dist)
            gen_lhapdf_dat_file(X,Q2,table, wdir, file_name,cnt)
            cnt+=1
        print

    #--gen_lhapdf_info_file
    X,Q2=gen_grid(dist)
    nrep=len(replicas)
    gen_lhapdf_info_file(X,Q2,nrep, wdir, file_name,info)

def rename_tables(dirname,newname):

    print('\nrenaming LHAPDF tables for %s to %s'%(dirname,newname))

    checkdir(newname)
    
    F=os.listdir(dirname)
    cnt=0
    for f in F:
        cnt+=1
        lprint('progress %d/%d'%(cnt,len(F)))
        cmd=['cp','%s/%s'%(dirname,f),'%s/%s'%(newname,f.replace(dirname,newname))]
        p=Popen(cmd,  stdout=PIPE, stderr=STDOUT)
        output = p.stdout.read()
    print
    








