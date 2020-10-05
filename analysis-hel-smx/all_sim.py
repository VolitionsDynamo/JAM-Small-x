import os,sys
#--matplotlib
import matplotlib
matplotlib.use('Agg')
import pylab  as py
import pandas as pd
import numpy as np
from subprocess import Popen, PIPE
from scipy.integrate import quad

#import lhapdf

#--from analysis
from analysis_smx.corelib import core
from analysis_smx.corelib import predict

#--from tools
from tools           import config
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config
from tools.inputmod  import INPUTMOD
from tools.randomstr import id_generator
from tools.config    import load_config,conf

#--from fitlib
from fitlib_smx.resman import RESMAN
from fitlib_smx.parman import PARMAN

#--from qcdlib
from qcdlib.aux import AUX
from qcdlib.alphaS import ALPHAS
from qcdlib.eweak import EWEAK

def ALL(wdir,tar='p',est='opt',lum='100:fb-1',force=True):

    #--generate initial data file
    gen_all_xlsx(wdir,tar,est)

    ##--modify conf with new data file
    conf = gen_conf(wdir,tar,est)

    ##--get predictions on new data file if not already done
    print('Generating predictions...')
    name = 'all-%s-%s'%(tar,est)
    #predict.get_predictions(wdir,force=force,mod_conf=conf,name=name)
    predict.get_predictions(wdir)#,force=force,mod_conf=conf,name=name)    

    ##--update tables
    update_tabs(wdir,tar,est,lum)

    #--plot asymmetry and errors
    plot(wdir,tar,est,lum)

    #--generate lhapdf info and data files
    gen_lhapdf_info_file(wdir,tar,est)
    gen_lhapdf_dat_file (wdir,tar,est)

#--generate pseudo-data
def gen_all_xlsx(wdir,tar,est):

    checkdir('%s/sim'%wdir)

    #-- the kinem. var.
    data={_:[] for _ in ['col','target','X','Xdo','Xup','Q2','Q2do','Q2up','obs','value','stat_u','syst_u','pole','polh','pole_u','polh_u','lum_u','RS']}

    #--get specific points from data file at fitpack/database/pvdis/expdata/1000.xlsx
    fdir = os.environ['FITPACK']
    grid = pd.read_excel(fdir + '/database/EIC/expdata/1000.xlsx')
    grid = grid.to_dict(orient='list')
    data['X']    = grid['X']
    data['Q2']   = grid['Q2']
    data['Xup']  = grid['Xup']
    data['Xdo']  = grid['Xdo']
    data['Q2up'] = grid['Q2up']
    data['Q2do'] = grid['Q2do']
    data['RS']   = grid['RS']

    obs = 'Apa'

    pole = 0.7
    polh = 0.7

    for i in range(len(data['X'])):
        data['col']   .append('JAM4EIC')
        data['target'].append(tar)
        data['obs']   .append(obs)
        data['value'] .append(0.0)
        data['stat_u'].append(1e-10)
        data['syst_u'].append(0.0)
        data['pole']  .append(pole)
        data['polh']  .append(polh)
        data['pole_u'].append(0.0)
        data['polh_u'].append(0.0)
        data['lum_u'] .append(0.0)

    df=pd.DataFrame(data)
    filename = '%s/sim/all-%s-%s.xlsx'%(wdir,tar,est)
    df.to_excel(filename, index=False)
    print('Generating xlsx file and saving to %s'%filename)

def gen_conf(wdir,tar,est):

    print('Modifying config with new experimental data file...')

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    conf['steps'][istep]['datasets'] = {}
    conf['steps'][istep]['datasets']['pidis']=[]
    conf['datasets']['pidis']['filters']=[]

    #--placeholder index
    idx = 80000
    conf['datasets']['pidis']['xlsx'][idx]='./%s/sim/all-%s-%s.xlsx'%(wdir,tar,est)
    conf['steps'][istep]['datasets']['pidis'].append(idx)

    fdir = os.environ['FITPACK']
    fn   = [fdir + '/database/EIC/expdata/1000.xlsx']

    conf['pidis grid']  = {}
    if tar =='p': conf['pidis grid']['overwrite'] = True
    if tar =='d': conf['pidis grid']['overwrite'] = False
    conf['pidis grid']['xlsx']  = fn
    conf['idis grid']   = {}
    if tar == 'p': conf['idis grid']['overwrite'] = True
    if tar == 'd': conf['idis grid']['overwrite'] = False
    conf['idis grid']['xlsx']   = fn

    return conf

def update_tabs(wdir,tar,est,lum):

    istep=core.get_istep()
    data=load('%s/data/predictions-%d-all-%s-%s.dat'%(wdir,istep,tar,est))

    blist=[]
    blist.append('thy')
    blist.append('shift')
    blist.append('residuals')
    blist.append('prediction')
    blist.append('N')
    blist.append('Shift')
    blist.append('W2')
    blist.append('alpha')
    blist.append('residuals-rep')
    blist.append('r-residuals')
    blist.append('L')
    blist.append('H')

    #--placeholder index
    idx = 80000
    tab=data['reactions']['pidis'][idx]

    #--delete unnecessary data
    for k in blist: 
        try:    del tab[k]
        except: continue

    #--save mean value
    tab['value']=np.mean(tab['prediction-rep'],axis=0)

    #--save individual values
    for i in range(len(tab['prediction-rep'])):
        tab['value%s'%(i+1)] = tab['prediction-rep'][i]

    del tab['prediction-rep']

    tab['stat_u'],tab['pole_u'],tab['polh_u'],tab['lum_u'],tab['syst_u'] = all_errors(wdir,tar,est,tab['value'],lum)

    df=pd.DataFrame(tab)
    filename = '%s/sim/all-%s-%s.xlsx'%(wdir,tar,est)
    df.to_excel(filename, index=False)
    print('Updating xlsx file and saving to %s'%filename)

def all_errors(wdir,tar,est,value,lum):

    conf['aux'] = AUX()
    data = pd.read_excel('%s/sim/all-%s-%s.xlsx'%(wdir,tar,est), index=False)
    data = data.to_dict(orient='list')
    target=data['target'][0]
    l    = len(data['value'])
    X    = np.array(data['X'])
    Q2   = np.array(data['Q2'])
    Xup  = np.array(data['Xup'])
    Xdo  = np.array(data['Xdo'])
    Q2up = np.array(data['Q2up'])
    Q2do = np.array(data['Q2do'])
    dx  = Xup  - Xdo
    dQ2 = Q2up - Q2do
    bins = dx*dQ2

    pole = data['pole'][0]
    polh = data['polh'][0]

    RS = data['RS'][0]
    S  = RS**2

    M2 = conf['aux'].M2

    if tar=='d': M2 = 4*M2

    #--luminosity
    lum = convert_lum(lum)

    GF = conf['aux'].GF

    #--get structure functions
    resman=RESMAN(parallel=False,datasets=False)
    parman = resman.parman
    resman.setup_idis()
    resman.setup_pidis()
    pidis = resman.pidis_thy
    idis  = resman.idis_thy
    pidis.data[tar]['g1'] = np.zeros(pidis.X.size)
    pidis.data[tar]['g2'] = np.zeros(pidis.X.size)
    idis.data [tar]['F2']  = np.zeros(idis.X.size)
    idis.data [tar]['FL']  = np.zeros(idis.X.size)
    idis   = resman.idis_thy
    pidis  = resman.pidis_thy

    idis._update()
    pidis._update()

    g1 = lambda x,q2: pidis.get_stf(x,q2,stf='g1',tar=tar) 
    g2 = lambda x,q2: pidis.get_stf(x,q2,stf='g2',tar=tar) 
    F2 = lambda x,q2: idis.get_stf(x,q2, stf='F2' ,tar=tar) 
    FL = lambda x,q2: idis.get_stf(x,q2, stf='FL' ,tar=tar) 

    rho2 = lambda x,q2: 1 + 4*x**2*M2/q2
    gamma2 = lambda x,q2: 4*x**2*M2/q2
    F1 = lambda x,q2: (gamma2(x,q2)*F2(x,q2) - FL(x,q2))/(2*x) 

    y=lambda x,q2: (q2/2/x)/((S-M2)/2)

    alpha = lambda q2: conf['eweak'].get_alpha(q2)

    C = lambda x,q2: 4*np.pi*alpha(q2)**2*S/q2

    T1 = lambda x,q2: x*F1(x,q2)
    T2 = lambda x,q2: (1-y(x,q2)-y(x,q2)**2/4*gamma2(x,q2))*F2(x,q2)/y(x,q2)

    T3 = lambda x,q2: x*(2-y(x,q2)-y(x,q2)**2/2*gamma2(x,q2))*g1(x,q2)
    T4 = lambda x,q2: x*y(x,q2)*gamma2(x,q2)*g2(x,q2)
    
    _ddsigR = lambda x,q2: C(x,q2)*(T1(x,q2) + T2(x,q2) - T3(x,q2) + T4(x,q2))
    _ddsigL = lambda x,q2: C(x,q2)*(T1(x,q2) + T2(x,q2) + T3(x,q2) - T4(x,q2))

    #--integrate over bin
    z1,w1 = np.polynomial.legendre.leggauss(3)
    z2,w2 = np.polynomial.legendre.leggauss(3)

    ddsigR = np.zeros((len(X),len(z1),len(z2)))
    ddsigL = np.zeros((len(X),len(z1),len(z2)))
    for i in range(len(X)):
        _x   = 0.5*((Xup[i] -Xdo[i])*z1  + Xup[i]  + Xdo[i])
        _q2  = 0.5*((Q2up[i]-Q2do[i])*z2 + Q2up[i] + Q2do[i])
        xjac  = 0.5*(Xup[i] -Xdo[i])
        q2jac = 0.5*(Q2up[i]-Q2do[i])
        for j in range(len(_x)):
            for k in range(len(_q2)):
                ddsigR[i][j][k] = _ddsigR(_x[j],_q2[k])*xjac*q2jac
                ddsigL[i][j][k] = _ddsigL(_x[j],_q2[k])*xjac*q2jac
   
    #--integrate over Q2
    dsigR = np.sum(w2*ddsigR,axis=2) 
    dsigL = np.sum(w2*ddsigL,axis=2) 

    #--integrate over X
    sigR = np.sum(w1*dsigR,axis=1) 
    sigL = np.sum(w1*dsigL,axis=1) 

    NR = lum*sigR
    NL = lum*sigL

    #--theoretical asymmetry
    A = np.array(value)
 
    #--absolute uncertainty (not divided by sigma) 
    stat2 = (1 + A**2)/(NR + NL)

    stat = np.sqrt(stat2)

    #--statistical uncertainties
    data['stat_u'] = stat

    Y = np.array(data['Q2'])/2/np.array(data['X'])/((S-M2)/2)

    #--add systemic uncertainties beyond polarization
    #--optimistic: 1% on lum, 1% on Pe, 2% on Ph, no further errors
    if est == 'opt':
        data['lum_u']  = A*0.01
        data['pole_u'] = A*0.01
        data['polh_u'] = A*0.02
        data['syst_u'] = A*0.0

    #--moderate: do no have yet
    elif est == 'mod':
        return
        #data['syst_u'] = np.zeros(l)
        #for i in range(l):
        #    if Y[i] <= 0.01: data['syst_u'][i] = Am[i]*0.020
        #    if Y[i] >  0.01: data['syst_u'][i] = Am[i]*0.025
    #--pessimistic: 2% on lum, 2% on Pe, 4% on Ph, no further errors
    elif est == 'pes':
        data['lum_u']  = A*0.02
        data['pole_u'] = A*0.02
        data['polh_u'] = A*0.04
        data['syst_u'] = A*0.0

    else:
        print('est must be opt, mod, or pes')
        return

    return data['stat_u'],data['pole_u'],data['polh_u'],data['lum_u'],data['syst_u']

def convert_lum(lum):
    one=0.3893793721  #--GeV2 mbarn from PDG
    lum,units=lum.split(':')
    lum=float(lum)
    units=units.strip()
    if units=='fb-1':   return lum*one*1e12
    else:               sys.exit('units not convertible!')

def plot(wdir,tar,est,lum):

    nrows,ncols=1,2
    fig = py.figure(figsize=(ncols*8,nrows*5))
    ax11=py.subplot(nrows,ncols,1)
    ax12=py.subplot(nrows,ncols,2)

    tab   = pd.read_excel('%s/sim/all-%s-%s.xlsx'%(wdir,tar,est))
    tab   = tab.to_dict(orient='list')
    X     = np.array(tab['X'])
    value = np.array(tab['value'])
    stat  = np.abs(np.array(tab['stat_u'])/value)
    pole  = np.abs(np.array(tab['pole_u'])/value)
    polh  = np.abs(np.array(tab['polh_u'])/value)
    lum   = np.abs(np.array(tab['lum_u']) /value)
    syst  = np.abs(np.array(tab['syst_u'])/value)

    #--check for zero systematic uncertainty
    zero = True
    for i in range(len(syst)):
        if syst[i] != 0: zero = False

    alpha = np.sqrt((stat**2 + pole**2 + polh**2 + lum**2 + syst**2))

    hand = {}
    hand['alpha']                 = ax11.scatter(X,alpha,color='red'    ,s=20,marker='s')
    hand['stat']                  = ax11.scatter(X,stat ,color='green'  ,s=10,marker='o')
    hand['pole']                  = ax11.scatter(X,pole ,color='blue'   ,s=10,marker='*')
    hand['lum']                   = ax11.scatter(X,lum  ,color='black'  ,s=10,marker='^')
    hand['polh']                  = ax11.scatter(X,polh ,color='orange' ,s=10,marker='*')
    if zero != True: hand['syst'] = ax11.scatter(X,syst ,color='magenta',s=10,marker='v')

    #--plot asymmetry
    x = X
    q2 = np.array(tab['Q2'])
    alpha  = alpha*value
    l = len(value)
    #--get fixed values for Q2
    Q2 = []
    for j in range(len(q2)):
        if q2[j] in Q2: continue
        Q2.append(q2[j])
        #--get corresponding x, asymmetry values
    for j in range(len(Q2)):
        Xp, valp, alpp = [],[],[]
        Xn, valn, alpn = [],[],[]
        for i in range(l):
            if q2[i] != Q2[j]: continue
            if value[i] > 0: 
                valp .append(value[i])
                Xp   .append(x[i])
                alpp .append(alpha[i])
            else:
                valn .append(value[i])
                Xn   .append(x[i])
                alpn .append(alpha[i])


        hand['pos'] = ax12.errorbar(Xp,valp        ,yerr=alpp,color='firebrick',fmt='o',ms=3.0)
        hand['neg'] = ax12.errorbar(Xn,np.abs(valn),yerr=alpn,color='darkgreen',fmt='^',ms=3.0)
        ax12.plot(Xp,valp        ,color='firebrick',ls='-',alpha=0.5)
        ax12.plot(Xn,np.abs(valn),color='darkgreen',ls='-',alpha=0.5)


    ax11.set_xlim(5e-5,1)
    ax11.semilogx()
    ax11.semilogy()

    ax11.set_ylim(8e-3,3e-1)
    ax11.set_ylabel(r'$|\sigma_{A_{LL}}/A_{LL}|$',size=30)
    if est == 'opt': ax11.text(0.6,0.6,r'\textrm{Optimistic}' ,transform = ax11.transAxes,size=30)
    if est == 'mod': ax11.text(0.6,0.6,r'\textrm{Moderate}'   ,transform = ax11.transAxes,size=30)
    if est == 'pes': ax11.text(0.6,0.6,r'\textrm{Pessimistic}',transform = ax11.transAxes,size=30)


    ax11.set_xlabel(r'$x$',size=30)

    ax11.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)


    if tar == 'p':   ax11.text(0.7,0.85,r'\textrm{Proton}'     ,transform = ax11.transAxes,size=30)
    if tar == 'd':   ax11.text(0.7,0.85,r'\textrm{Deuteron}'   ,transform = ax11.transAxes,size=30)

    handles = [hand['alpha'],hand['stat'],hand['pole'],hand['polh'],hand['lum']]
    label1 = r'\textbf{\textrm{Total}}'
    label2 = r'\textbf{\textrm{Stat}}'
    label3 = r'\textbf{\textrm{Pole}}'
    label4 = r'\textbf{\textrm{Polh}}'
    label5 = r'\textbf{\textrm{Lum}}'
    labels = [label1,label2,label3,label4,label5]

    if zero != True:
        handles.append(hand['syst'])
        labels.append(r'\textbf{\textrm{Syst}}')

    ax11.legend(handles,labels,loc='lower left', fontsize = 20, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    ax12.semilogx()
    ax12.semilogy()
    ax12.set_xlabel(r'$x$',size=30)
    ax12.set_xlim(3e-4,4.0)
    ax12.set_ylim(2e-5,5e-1)
    ax12.set_ylabel(r'$A_{LL}$',size=30)

    ax12.text(0.60,0.20,r'$Q^2=%s$'%np.round(Q2[0],1) + ' ' + r'\textrm{GeV}' + r'$^2$',transform=ax12.transAxes,size=14)
    ax12.text(0.66,0.28,r'$Q^2=%s$'%np.round(Q2[1],1)                                  ,transform=ax12.transAxes,size=14)
    ax12.text(0.73,0.38,r'$Q^2=%s$'%np.round(Q2[2],1)                                  ,transform=ax12.transAxes,size=14)
    ax12.text(0.80,0.45,r'$Q^2=%s$'%np.round(Q2[3],1)                                  ,transform=ax12.transAxes,size=14)
    ax12.text(0.80,0.55,r'$Q^2=%s$'%np.round(Q2[4],1)                                  ,transform=ax12.transAxes,size=14)
    ax12.text(0.80,0.65,r'$Q^2=%s$'%np.round(Q2[5],1)                                  ,transform=ax12.transAxes,size=14)
    ax12.text(0.80,0.75,r'$Q^2=%s$'%np.round(Q2[6],1)                                  ,transform=ax12.transAxes,size=14)
    ax12.text(0.80,0.85,r'$Q^2=%s$'%np.round(Q2[7],1)                                  ,transform=ax12.transAxes,size=14)
    ax12.text(0.80,0.95,r'$Q^2=%s$'%np.round(Q2[8],1)                                  ,transform=ax12.transAxes,size=14)

    if est == 'opt': ax12.text(0.4,0.9,r'\textrm{Optimistic}' ,transform=ax12.transAxes,size=20)
    if est == 'mod': ax12.text(0.4,0.9,r'\textrm{Moderate}'   ,transform=ax12.transAxes,size=20)
    if est == 'pes': ax12.text(0.4,0.9,r'\textrm{Pessimistic}',transform=ax12.transAxes,size=20)

    ax12.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)

    handles = [hand['pos'],hand['neg']]
    label1  = r'\textbf{\textrm{JAM4EIC(+)}}'
    label2  = r'\textbf{\textrm{JAM4EIC(-)}}'
    labels  = [label1,label2]
    ax12.legend(handles,labels,loc='upper left', fontsize = 20, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.tight_layout()
    checkdir('%s/gallery'%wdir)
    filename = '%s/gallery/all-sim-%s-%s'%(wdir,tar,est)
    py.savefig(filename)
    print('Saving A_LL plot to %s'%filename)
    py.clf()

#--generate lhapdf info and data files
def gen_lhapdf_info_file(wdir,tar,est):

    info={}
    info['<description>'] = 'Double longitudinal asymmetry'
    info['<index>']       = '0'
    info['<authors>']     = 'JAM Collaboration'
    info['<reference>']   = ''
    info['<particle>']    = '%s'%tar

    #--get tables
    X,Q2,table,replicas=get_tables(wdir,tar,est)

    #--kinematic limits
    xmin=X[0]
    xmax=X[-1]
    Qmin=Q2[0]**0.5
    Qmax=Q2[-1]**0.5

    #--qcd params
    load_config('%s/input.py'%wdir)
    RESMAN(nworkers=1,parallel=False,datasets=False)
    aS=[conf['alphaS'].get_alphaS(_) for _ in Q2]
    mZ=conf['aux'].mZ
    mb=conf['aux'].mb
    mc=conf['aux'].mc
    alphaSMZ=conf['aux'].alphaSMZ

    #--begin lhapdf info file
    lines=[]
    lines.append('SetDesc:         "<description>"')
    lines.append('SetIndex:        <index>')
    lines.append('Authors:         <authors>')
    lines.append('Reference:       <reference>')
    lines.append('Format:          lhagrid1')
    lines.append('DataVersion:     1')
    lines.append('NumMembers:      1')
    lines.append('Particle:        <particle>')
    lines.append('Flavors:         [80000]')
    lines.append('OrderQCD:        1')
    lines.append('FlavorScheme:    <flav scheme>')
    lines.append('NumFlavors:      1')
    lines.append('ErrorType:       no error')
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
    dirname = 'all-%s'%(tar)
    checkdir('%s/lhapdf/%s'%(wdir,dirname))
    tab=open('%s/lhapdf/%s/%s.info'%(wdir,dirname,dirname),'w')
    tab.writelines(lines)
    tab.close()
    print('Saving lhapdf info file to %s/lhapdf/%s/%s.info'%(wdir,dirname,dirname))

def gen_lhapdf_dat_file(wdir,tar,est):

    #--get tables
    X,Q2,central,replicas=get_tables(wdir,tar,est)
    nx=len(X)
    nQ2=len(Q2)
    nrep = len(replicas)
    L = nx


    #--central value
    #--start lhapdf file
    lines=[]
    lines.append('PdfType: central')
    lines.append('Format: lhagrid1')
    lines.append('---')
    line=''
    for _ in X: line+=('%10.5e '%_).upper()
    lines.append(line)
    line=''
    for _ in Q2: line+=('%10.5e '%_**0.5).upper()
    lines.append(line)
    flavs=''
    flavs+='90000 '
    lines.append(flavs)

    for i in range(L):
        line=''
        line+=('%10.5e '%central[i]).upper()
        lines.append(line)
    lines.append('---')
    lines=[l+'\n' for l in lines]
    if central: idx=str(0).zfill(4)
    else:       idx=str(0).zfill(4)

    dirname = 'all-%s'%(tar)
    checkdir('%s/lhapdf/%s'%(wdir,dirname))
    tab=open('%s/lhapdf/%s/%s_%s.dat'%(wdir,dirname,dirname,idx),'w')
    tab.writelines(lines)
    tab.close()

    #--individual replicas
    for j in range(nrep): 
        #--start lhapdf file
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
        flavs=''
        flavs+='90000 '
        lines.append(flavs)

        for i in range(L):
            line=''
            line+=('%10.5e '%replicas[j][i]).upper()
            lines.append(line)
        lines.append('---')
        lines=[l+'\n' for l in lines]
        l = len(str(j+1))
        idx=str(0).zfill(4-l)+str((j+1))

        tab=open('%s/lhapdf/%s/%s_%s.dat'%(wdir,dirname,dirname,idx),'w')
        tab.writelines(lines)
        tab.close()

    print('Saving lhapdf data files inside %s/lhapdf/%s'%(wdir,dirname))

def get_tables(wdir,tar,est):
    replicas = []
    tab=pd.read_excel('%s/sim/all-%s-%s.xlsx'%(wdir,tar,est))
    tab=tab.to_dict(orient='list')
    central = tab['value']
    _replicas = {}
    for key in tab:
        if key[:5] != 'value': continue
        if key     == 'value': continue
        _replicas[key] = tab[key]
    for i in range(len(_replicas)):
        replicas.append(_replicas['value%s'%(i+1)])

    #--get specific points from data file at fitpack/database/pvdis/expdata/1000.xlsx
    fdir = os.environ['FITPACK']
    grid = pd.read_excel(fdir + '/database/EIC/expdata/1000.xlsx')
    grid = grid.to_dict(orient='list')
    X    = grid['X']
    Q2   = grid['Q2']

    return X,Q2,central,replicas







 
