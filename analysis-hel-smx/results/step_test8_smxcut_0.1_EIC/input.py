import os
import numpy as np
conf={}

## block comment can not appear as this file will be executed by 'exec'
## setup posterior sampling

conf['bootstrap'] = True
# conf['flat par'] = True ## this has to be set to 'True' for random parameter generation
conf['flat par'] = True
conf['ftol'] = 1e-8

## setup qcd evolution

conf['dglap mode'] = 'truncated'
conf['alphaSmode'] = 'backward'
conf['order'] = 'LO'
conf['Q20']   = 1.27 ** 2.0

## setup for idis

conf['tmc']   = False
conf['ht']  = False
conf['nuc'] = False
conf['sidis nuc smearing'] = False
conf['hq'] = False

## grids

conf['path2idistab']   = '%s/grids/grids-idis/notmc4/' % os.environ['FITPACK']
conf['path2pidistab']  = '%s/grids/grids-pidis/notmc4/' % os.environ['FITPACK']
conf['path2dytab']     = '%s/grids/grids-dy/' % os.environ['FITPACK']
conf['path2jettab']    = '%s/grids/grids-jets' % os.environ['FITPACK']
conf['path2pjettab']    = '%s/grids/grids-pjets' % os.environ['FITPACK']

##small-x helicity
conf['ww']=False #False = set g2 to zero

conf['smx_hel']= True #True = include small-x helicity evo

conf['smx']={}
conf['smx']['alphas_fixed']=0.3
conf['smx']['smxcut']=0.1
conf['smx']['x0']=1.0
conf['smx']['deta']=0.02 #0.03
conf['smx']['eta_max']=6.5 #11.0
conf['smx']['ds']=0.02 #0.03
conf['smx']['s_max']=6.5 #11.0
conf['smx']['var']=['s10','eta','1']

smxcut=conf['smx']['smxcut']

## datasets

conf['datasets'] = {}

## lepton-hadron reactions

Q2cut = 1.3 ** 2.0
W2cut = 4.0 #10.0
#jet_pt_cut = 10.0 ## pt cut for unpolarized JET dataset
#pjet_pt_cut = 10.0 ## pt cut for polarized JET dataset

# IDIS
#conf['datasets']['idis'] = {}
#conf['datasets']['idis']['filters'] = []
#conf['datasets']['idis']['filters'].append("Q2>%f" % Q2cut)
#conf['datasets']['idis']['filters'].append("W2>%f" % W2cut)
#conf['datasets']['idis']['xlsx'] = {}
#conf['datasets']['idis']['xlsx'][10010] = 'idis/expdata/10010.xlsx' # proton   | F2            | SLAC
#conf['datasets']['idis']['xlsx'][10011] = 'idis/expdata/10011.xlsx' # deuteron | F2            | SLAC
#conf['datasets']['idis']['xlsx'][10016] = 'idis/expdata/10016.xlsx' # proton   | F2            | BCDMS
#conf['datasets']['idis']['xlsx'][10017] = 'idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS
#conf['datasets']['idis']['xlsx'][10020] = 'idis/expdata/10020.xlsx' # proton   | F2            | NMC
#conf['datasets']['idis']['xlsx'][10021] = 'idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC
##conf['datasets']['idis']['xlsx'][10026] = 'idis/expdata/10026.xlsx' # proton   | sigma red     | HERA II NC e+ (1)
##conf['datasets']['idis']['xlsx'][10027] = 'idis/expdata/10027.xlsx' # proton   | sigma red     | HERA II NC e+ (2)
##conf['datasets']['idis']['xlsx'][10028] = 'idis/expdata/10028.xlsx' # proton   | sigma red     | HERA II NC e+ (3)
##conf['datasets']['idis']['xlsx'][10029] = 'idis/expdata/10029.xlsx' # proton   | sigma red     | HERA II NC e+ (4)
##conf['datasets']['idis']['xlsx'][10030] = 'idis/expdata/10030.xlsx' # proton   | sigma red     | HERA II NC e-
##conf['datasets']['idis']['xlsx'][10031] = 'idis/expdata/10031.xlsx' # proton   | sigma red     | HERA II CC e+
##conf['datasets']['idis']['xlsx'][10032] = 'idis/expdata/10032.xlsx' # proton   | sigma red     | HERA II NC e-
#conf['datasets']['idis']['norm'] = {}
#conf['datasets']['idis']['norm'][10010] = {'value': 1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10011] = {'value': 1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10016] = {'value': 1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10017] = {'value': 1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10020] = {'value': 1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}

# PIDIS
conf['datasets']['pidis'] = {}
conf['datasets']['pidis']['filters'] = []
conf['datasets']['pidis']['filters'].append("Q2>%f" %Q2cut)
conf['datasets']['pidis']['filters'].append("W2>%f" % W2cut)
conf['datasets']['pidis']['filters'].append("X<=%f" % smxcut)
conf['datasets']['pidis']['xlsx'] = {}
## --------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['xlsx'][10002] = 'pidis/expdata/10002.xlsx' # 10002 | proton   | A1   | COMPASS         |          |
conf['datasets']['pidis']['xlsx'][10003] = 'pidis/expdata/10003.xlsx' # 10003 | proton   | A1   | COMPASS         |          |
conf['datasets']['pidis']['xlsx'][10004] = 'pidis/expdata/10004.xlsx' # 10004 | proton   | A1   | EMC             |          |
conf['datasets']['pidis']['xlsx'][10007] = 'pidis/expdata/10007.xlsx' # 10007 | proton   | Apa  | HERMES          |          |
#conf['datasets']['pidis']['xlsx'][10008] = 'pidis/expdata/10008.xlsx' # 10008 | proton   | A2   | HERMES          |          |
conf['datasets']['pidis']['xlsx'][10017] = 'pidis/expdata/10017.xlsx' # 10017 | proton   | Apa  | JLabHB(EG1DVCS) |          |
conf['datasets']['pidis']['xlsx'][10022] = 'pidis/expdata/10022.xlsx' # 10022 | proton   | Apa  | SLAC(E143)      |          |
#conf['datasets']['pidis']['xlsx'][10023] = 'pidis/expdata/10023.xlsx' # 10023 | proton   | Ape  | SLAC(E143)      |          |
#conf['datasets']['pidis']['xlsx'][10028] = 'pidis/expdata/10028.xlsx' # 10028 | proton   | Ape  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10029] = 'pidis/expdata/10029.xlsx' # 10029 | proton   | Apa  | SLAC(E155)      |          |
#conf['datasets']['pidis']['xlsx'][10031] = 'pidis/expdata/10031.xlsx' # 10031 | proton   | Atpe | SLAC(E155x)     |          |
conf['datasets']['pidis']['xlsx'][10032] = 'pidis/expdata/10032.xlsx' # 10032 | proton   | Apa  | SLACE80E130     |          |
conf['datasets']['pidis']['xlsx'][10035] = 'pidis/expdata/10035.xlsx' # 10035 | proton   | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10036] = 'pidis/expdata/10036.xlsx' # 10036 | proton   | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10041] = 'pidis/expdata/10041.xlsx' # 10041 | proton   | Apa  | JLabHB(EG1b)    | E =1 GeV |
conf['datasets']['pidis']['xlsx'][10042] = 'pidis/expdata/10042.xlsx' # 10042 | proton   | Apa  | JLabHB(EG1b)    | E =2 GeV |
conf['datasets']['pidis']['xlsx'][10043] = 'pidis/expdata/10043.xlsx' # 10043 | proton   | Apa  | JLabHB(EG1b)    | E =4 GeV |
conf['datasets']['pidis']['xlsx'][10044] = 'pidis/expdata/10044.xlsx' # 10044 | proton   | Apa  | JLabHB(EG1b)    | E =5 GeV |
conf['datasets']['pidis']['xlsx'][10005] = 'pidis/expdata/10005.xlsx' # 10005 | neutron  | A1   | HERMES          |          |
## --------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['xlsx'][10001] = 'pidis/expdata/10001.xlsx' # 10001 | deuteron | A1   | COMPASS         |          |
conf['datasets']['pidis']['xlsx'][10006] = 'pidis/expdata/10006.xlsx' # 10006 | deuteron | Apa  | HERMES          |          |
conf['datasets']['pidis']['xlsx'][10016] = 'pidis/expdata/10016.xlsx' # 10016 | deuteron | Apa  | JLabHB(EG1DVCS) |          |
#conf['datasets']['pidis']['xlsx'][10020] = 'pidis/expdata/10020.xlsx' # 10020 | deuteron | Ape  | SLAC(E143)      |          |
conf['datasets']['pidis']['xlsx'][10021] = 'pidis/expdata/10021.xlsx' # 10021 | deuteron | Apa  | SLAC(E143)      |          |
#conf['datasets']['pidis']['xlsx'][10026] = 'pidis/expdata/10026.xlsx' # 10026 | deuteron | Ape  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10027] = 'pidis/expdata/10027.xlsx' # 10027 | deuteron | Apa  | SLAC(E155)      |          |
#conf['datasets']['pidis']['xlsx'][10030] = 'pidis/expdata/10030.xlsx' # 10030 | deuteron | Atpe | SLAC(E155x)     |          |
conf['datasets']['pidis']['xlsx'][10033] = 'pidis/expdata/10033.xlsx' # 10033 | deuteron | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10034] = 'pidis/expdata/10034.xlsx' # 10034 | deuteron | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10037] = 'pidis/expdata/10037.xlsx' # 10037 | deuteron | Apa  | JLabHB(EG1b)    | E =1 GeV |
conf['datasets']['pidis']['xlsx'][10038] = 'pidis/expdata/10038.xlsx' # 10038 | deuteron | Apa  | JLabHB(EG1b)    | E =2 GeV |
conf['datasets']['pidis']['xlsx'][10039] = 'pidis/expdata/10039.xlsx' # 10039 | deuteron | Apa  | JLabHB(EG1b)    | E =4 GeV |
conf['datasets']['pidis']['xlsx'][10040] = 'pidis/expdata/10040.xlsx' # 10040 | deuteron | Apa  | JLabHB(EG1b)    | E =5 GeV |

conf['datasets']['pidis']['xlsx'][70000] = 'pidis/expdata/70000.xlsx' # 70000 | proton | Apa  | EIC pseudo-data    |

## --------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['norm'] = {}
conf['datasets']['pidis']['norm'][10002] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10003] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10004] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10022] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['pidis']['norm'][10023] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10029] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['pidis']['norm'][10031] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10041] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
## ---------------------------------------------------------------------------------------------------------------------------
#conf['datasets']['pidis']['norm'][10020] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10021] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10001] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10027] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['pidis']['norm'][10030] = {'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}

## parameters
conf['params'] = {}

## pdf parameters

conf['pdf parametrization'] = 4
conf['params']['pdf'] = {}

# first shape of PDF
conf['params']['pdf']['g1 N']    ={'value':    4.57778e-01, 'min':  None, 'max':  None, 'fixed': True}
conf['params']['pdf']['g1 a']    ={'value':   -6.93460e-01, 'min':  -1.9, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['g1 b']    ={'value':    4.86243e+00, 'min':     0, 'max':    10, 'fixed': True} #False}

conf['params']['pdf']['uv1 N']   ={'value':    2.40182e-01, 'min':  None, 'max':  None, 'fixed': True}
conf['params']['pdf']['uv1 a']   ={'value':   -5.00000e-01, 'min':  -0.5, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['uv1 b']   ={'value':    2.66351e+00, 'min':     0, 'max':    10, 'fixed': True} #False}

conf['params']['pdf']['dv1 N']   ={'value':    5.72486e-02, 'min':  None, 'max':  None, 'fixed': True}
conf['params']['pdf']['dv1 a']   ={'value':   -4.63416e-01, 'min':  -0.5, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['dv1 b']   ={'value':    7.83629e+00, 'min':     0, 'max':    10, 'fixed': True} #False}

conf['params']['pdf']['db1 N']   ={'value':    4.89659e-02, 'min':     0, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['db1 a']   ={'value':    3.30067e-01, 'min':    -1, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['db1 b']   ={'value':    4.30919e+00, 'min':     0, 'max':    10, 'fixed': True} #False}

conf['params']['pdf']['ub1 N']   ={'value':    3.79371e-02, 'min':     0, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['ub1 a']   ={'value':    1.00000e+00, 'min':    -1, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['ub1 b']   ={'value':    4.03210e+00, 'min':     0, 'max':    10, 'fixed': True} #False}

conf['params']['pdf']['s1 N']    ={'value':    2.11838e-20, 'min':     0, 'max':     1, 'fixed': True}
conf['params']['pdf']['s1 a']    ={'value':    4.68215e-01, 'min':    -1, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['s1 b']    ={'value':    9.99871e+00, 'min':     0, 'max':    10, 'fixed': True} #False}

conf['params']['pdf']['sb1 N']   ={'value':    6.16314e-20, 'min':     0, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['sb1 a']   ={'value':    2.11638e-01, 'min':    -1, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['sb1 b']   ={'value':    1.32463e+00, 'min':     0, 'max':    10, 'fixed': True} #False}

conf['params']['pdf']['sea1 N']  ={'value':    1.18309e-02, 'min':     0, 'max':     1, 'fixed': True} #False}
conf['params']['pdf']['sea1 a']  ={'value':   -1.21623e+00, 'min':  -1.9, 'max':    -1, 'fixed': True} #False}
conf['params']['pdf']['sea1 b']  ={'value':    1.00000e+01, 'min':     0, 'max':    10, 'fixed': True} #False}

conf['params']['pdf']['sea2 N']  ={'value':    1.18309e-02, 'min':     0, 'max':     1, 'fixed': 'sea1 N'}
conf['params']['pdf']['sea2 a']  ={'value':   -1.21623e+00, 'min':  -1.9, 'max':    -1, 'fixed': 'sea1 a'}
conf['params']['pdf']['sea2 b']  ={'value':    1.00000e+01, 'min':     0, 'max':    10, 'fixed': 'sea1 b'}

# second shape of PDF
#conf['params']['pdf']['g2 N']    = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': True}
#conf['params']['pdf']['g2 a']    = {'value': -0.5, 'min':  -1.9, 'max':     1, 'fixed': True}
#conf['params']['pdf']['g2 b']    = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['uv2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': True}
#conf['params']['pdf']['uv2 a']   = {'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': True}
#conf['params']['pdf']['uv2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['dv2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': True}
#conf['params']['pdf']['dv2 a']   = {'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': True}
#conf['params']['pdf']['dv2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['db2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': True}
#conf['params']['pdf']['db2 a']   = {'value': -0.5, 'min':    -1, 'max':     1, 'fixed': True}
#conf['params']['pdf']['db2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['ub2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': True}
#conf['params']['pdf']['ub2 a']   = {'value': -0.5, 'min':    -1, 'max':     1, 'fixed': True}
#conf['params']['pdf']['ub2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['s2 N']    = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': True}
#conf['params']['pdf']['s2 a']    = {'value':    0, 'min':    -1, 'max':     1, 'fixed': True}
#conf['params']['pdf']['s2 b']    = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['sb2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': True}
#conf['params']['pdf']['sb2 a']   = {'value':    0, 'min':    -1, 'max':     1, 'fixed': True}
#conf['params']['pdf']['sb2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}

#conf['params']['pdf']['g2 N']    = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': False, 'zero': True}
#conf['params']['pdf']['g2 a']    = {'value': -0.5, 'min':  -1.9, 'max':     1, 'fixed': False, 'prior': 'g1 a'}
#conf['params']['pdf']['g2 b']    = {'value':    6, 'min':     0, 'max':    10, 'fixed': False, 'prior': 'g1 b'}
#
#conf['params']['pdf']['uv2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': False, 'zero': True}
#conf['params']['pdf']['uv2 a']   = {'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': False, 'prior': 'uv1 a'}
#conf['params']['pdf']['uv2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': False, 'prior': 'uv1 b'}
#
#conf['params']['pdf']['dv2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': False, 'zero': True}
#conf['params']['pdf']['dv2 a']   = {'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': False, 'prior': 'dv1 a'}
#conf['params']['pdf']['dv2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': False, 'prior': 'dv1 b'}
#
#conf['params']['pdf']['db2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': False, 'zero': True}
#conf['params']['pdf']['db2 a']   = {'value': -0.5, 'min':    -1, 'max':     1, 'fixed': False, 'prior': 'db1 a'}
#conf['params']['pdf']['db2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': False, 'prior': 'db1 b'}
#
#conf['params']['pdf']['ub2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': False, 'zero': True}
#conf['params']['pdf']['ub2 a']   = {'value': -0.5, 'min':    -1, 'max':     1, 'fixed': False, 'prior': 'ub1 a'}
#conf['params']['pdf']['ub2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': False, 'prior': 'ub1 b'}
#
#conf['params']['pdf']['s2 N']    = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': False, 'zero': True}
#conf['params']['pdf']['s2 a']    = {'value':    0, 'min':    -1, 'max':     1, 'fixed': False, 'prior': 's1 a'}
#conf['params']['pdf']['s2 b']    = {'value':    6, 'min':     0, 'max':    10, 'fixed': False, 'prior': 's1 b'}
#
#conf['params']['pdf']['sb2 N']   = {'value':  0.0, 'min':  -0.5, 'max':   0.5, 'fixed': False, 'zero': True}
#conf['params']['pdf']['sb2 a']   = {'value':    0, 'min':    -1, 'max':     1, 'fixed': False, 'prior': 'sb1 a'}
#conf['params']['pdf']['sb2 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': False, 'prior': 'sb1 b'}

# third shape of PDF
#conf['params']['pdf']['g3 N']    = {'value':  0.0, 'min':  None, 'max':  None, 'fixed': True}
#conf['params']['pdf']['g3 a']    = {'value': -0.5, 'min':  -1.9, 'max':     1, 'fixed': True}
#conf['params']['pdf']['g3 b']    = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['uv3 N']   = {'value':  0.0, 'min':  None, 'max':  None, 'fixed': True}
#conf['params']['pdf']['uv3 a']   = {'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': True}
#conf['params']['pdf']['uv3 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['dv3 N']   = {'value':  0.0, 'min':  None, 'max':  None, 'fixed': True}
#conf['params']['pdf']['dv3 a']   = {'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': True}
#conf['params']['pdf']['dv3 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['db3 N']   = {'value':  0.0, 'min':     0, 'max':     1, 'fixed': True}
#conf['params']['pdf']['db3 a']   = {'value': -0.5, 'min':    -1, 'max':     1, 'fixed': True}
#conf['params']['pdf']['db3 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['ub3 N']   = {'value':  0.0, 'min':     0, 'max':     1, 'fixed': True}
#conf['params']['pdf']['ub3 a']   = {'value': -0.5, 'min':    -1, 'max':     1, 'fixed': True}
#conf['params']['pdf']['ub3 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['s3 N']    = {'value':  0.0, 'min':     0, 'max':     1, 'fixed': True}
#conf['params']['pdf']['s3 a']    = {'value':    0, 'min':    -1, 'max':     1, 'fixed': True}
#conf['params']['pdf']['s3 b']    = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#conf['params']['pdf']['sb3 N']   = {'value':  0.0, 'min':     0, 'max':     1, 'fixed': True}
#conf['params']['pdf']['sb3 a']   = {'value':    0, 'min':    -1, 'max':     1, 'fixed': True}
#conf['params']['pdf']['sb3 b']   = {'value':    6, 'min':     0, 'max':    10, 'fixed': True}
#
#ppdf parameters
conf['ppdf parametrization'] = 0
conf['su2+su3'] = True
conf['params']['ppdf'] = {}

conf['params']['ppdf']['g1 N']    ={'value':    0.70000e+00, 'min':   0.5, 'max':    2, 'fixed': True}
conf['params']['ppdf']['g1 a']    ={'value':    7.00000e+00, 'min':   5, 'max':     10, 'fixed': True}
conf['params']['ppdf']['g1 b']    ={'value':    3.00000e+00, 'min':     2, 'max':    5, 'fixed': True}
conf['params']['ppdf']['g1 c']    ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed':True}
conf['params']['ppdf']['g1 d']    ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed':True}

conf['params']['ppdf']['up1 N']   ={'value':    0e+00, 'min':    None, 'max': None, 'fixed':True}
conf['params']['ppdf']['up1 a']   ={'value':   -2.82670e-01, 'min':    -0.3, 'max':     0.2, 'fixed':True}
conf['params']['ppdf']['up1 b']   ={'value':    1.41541e+00, 'min':     1, 'max':    2.5, 'fixed':True}
conf['params']['ppdf']['up1 c']    ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed':True}
conf['params']['ppdf']['up1 d']   ={'value':    0.0000e+00, 'min':   -10, 'max':    10, 'fixed': True}

conf['params']['ppdf']['dp1 N']   ={'value':    0, 'min':None, 'max':None, 'fixed':True}
conf['params']['ppdf']['dp1 a']   ={'value':   -4.52161e-01, 'min':    -0.5, 'max':     0.2, 'fixed':True}
conf['params']['ppdf']['dp1 b']   ={'value':    3.99563e+00, 'min':     1, 'max':    5, 'fixed':True}
conf['params']['ppdf']['dp1 c']   ={'value':   0.0, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf']['dp1 d']   ={'value':   0.0, 'min':   -10, 'max':    10, 'fixed': True}

conf['params']['ppdf']['sp1 N']   ={'value':    -1.52130e-01, 'min':   -0.2, 'max':    -0.1, 'fixed':True}
conf['params']['ppdf']['sp1 a']   ={'value':    1.27210e+00, 'min':    0.5, 'max': 2, 'fixed':True}
conf['params']['ppdf']['sp1 b']   ={'value':    1.00000, 'min':     0, 'max': 2, 'fixed':True}
conf['params']['ppdf']['sp1 c']   ={'value':    0.000e+00, 'min':   -10, 'max':    10, 'fixed':True}
conf['params']['ppdf']['sp1 d']   ={'value':    0.000e+00, 'min':   -10, 'max':    10, 'fixed':True}

conf['params']['ppdf']['um1 N']   ={'value':    0.00000e+00, 'min':   -1, 'max':    1, 'fixed': True}
conf['params']['ppdf']['um1 a']   ={'value':    0.00000e+00, 'min':    -1, 'max':     2, 'fixed': True}
conf['params']['ppdf']['um1 b']   ={'value':    0.00000e+00, 'min':     0, 'max':    10, 'fixed': True}
conf['params']['ppdf']['um1 c']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf']['um1 d']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}

conf['params']['ppdf']['dm1 N']   ={'value':    0.00000e+00, 'min':   -1, 'max':    1, 'fixed': True}
conf['params']['ppdf']['dm1 a']   ={'value':    0.00000e+00, 'min':    -1, 'max':     2, 'fixed': True}
conf['params']['ppdf']['dm1 b']   ={'value':    0.00000e+00, 'min':     0, 'max':    10, 'fixed': True}
conf['params']['ppdf']['dm1 c']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf']['dm1 d']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}

conf['params']['ppdf']['sm1 N']   ={'value':    0.00000e+00, 'min':   -1, 'max':    1, 'fixed': True}
conf['params']['ppdf']['sm1 a']   ={'value':    0.00000e+00, 'min':    -1, 'max':     2, 'fixed': True}
conf['params']['ppdf']['sm1 b']   ={'value':    0.00000e+00, 'min':     0, 'max':    10, 'fixed': True}
conf['params']['ppdf']['sm1 c']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf']['sm1 d']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}

#ppdf_smx parameters
conf['params']['ppdf_smx'] = {}

conf['params']['ppdf_smx']['g1 a']   ={'value':    0.00000e+00, 'min':    -10, 'max':     10, 'fixed': True}
conf['params']['ppdf_smx']['g1 b']   ={'value':    0.00000e+00, 'min':     -10, 'max':    10, 'fixed': True}
conf['params']['ppdf_smx']['g1 c']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf_smx']['g1 lam2smx']   ={'value':    0.04, 'min':   0.01, 'max':    1., 'fixed': True}

conf['params']['ppdf_smx']['up1 a']   ={'value':    0, 'min':    -25, 'max':  7, 'fixed':False}
conf['params']['ppdf_smx']['up1 b']   ={'value':    0., 'min':     -10, 'max':    25, 'fixed':False}
conf['params']['ppdf_smx']['up1 c']   ={'value':    0., 'min':   -15, 'max':    10, 'fixed':False}
conf['params']['ppdf_smx']['up1 lam2smx']   ={'value':    1.49500e-02, 'min':   0.01, 'max':    1.0, 'fixed':False}

conf['params']['ppdf_smx']['dp1 a']   ={'value':    0, 'min':    -10, 'max':     25, 'fixed':False}
conf['params']['ppdf_smx']['dp1 b']   ={'value':    0, 'min':     -10, 'max':    10, 'fixed':False}
conf['params']['ppdf_smx']['dp1 c']   ={'value':    0., 'min':   -10, 'max':    15.0, 'fixed':False}
conf['params']['ppdf_smx']['dp1 lam2smx']   ={'value':    1.49500e-02, 'min':   0.01, 'max':    1.0, 'fixed': 'up1 lam2smx'}

conf['params']['ppdf_smx']['sp1 a']   ={'value':    0, 'min':    -10, 'max':    2, 'fixed':False}
conf['params']['ppdf_smx']['sp1 b']   ={'value':    0., 'min':     -100, 'max':    25, 'fixed':False}
conf['params']['ppdf_smx']['sp1 c']   ={'value':    0., 'min':   -10, 'max':    50, 'fixed':False}
conf['params']['ppdf_smx']['sp1 lam2smx']   ={'value':    1.49500e-02, 'min':   0.01, 'max':    1.0, 'fixed': 'up1 lam2smx'}

conf['params']['ppdf_smx']['um1 a']   ={'value':    0.00000e+00, 'min':    -10, 'max':     10, 'fixed':True}
conf['params']['ppdf_smx']['um1 b']   ={'value':    0.00000e+00, 'min':     -10, 'max':    10, 'fixed':True}
conf['params']['ppdf_smx']['um1 c']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf_smx']['um1 lam2smx']   ={'value':    1.49500e-02, 'min':   0.01, 'max':    1., 'fixed': 'up1 lam2smx'}

conf['params']['ppdf_smx']['dm1 a']   ={'value':    0.00000e+00, 'min':    -10, 'max':     10, 'fixed': True}
conf['params']['ppdf_smx']['dm1 b']   ={'value':    0.00000e+00, 'min':     -10, 'max':    10, 'fixed': True}
conf['params']['ppdf_smx']['dm1 c']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf_smx']['dm1 lam2smx']   ={'value':    4.00000e-02, 'min':   0.01, 'max':    1., 'fixed': 'dm1 lam2smx'}

conf['params']['ppdf_smx']['sm1 a']   ={'value':    0.00000e+00, 'min':    -10, 'max':     10, 'fixed': True}
conf['params']['ppdf_smx']['sm1 b']   ={'value':    0.00000e+00, 'min':     -10, 'max':    10, 'fixed': True}
conf['params']['ppdf_smx']['sm1 c']   ={'value':    0.00000e+00, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf_smx']['sm1 lam2smx']   ={'value':    1.49500e-02, 'min':   0.01, 'max':    1., 'fixed': 'up1 lam2smx'}



# steps
conf['steps'] = {}

# DIS without HERA
# fit only first shape of PDF
#conf['steps'][1] = {}
#conf['steps'][1]['dep'] = []
#conf['steps'][1]['active distributions'] = ['pdf']
#conf['steps'][1]['passive distributions'] = []
#conf['steps'][1]['datasets'] = {}
#conf['steps'][1]['datasets']['idis'] = []
#conf['steps'][1]['datasets']['idis'].append(10010) ## proton   | F2            | SLAC
#conf['steps'][1]['datasets']['idis'].append(10011) ## deuteron | F2            | SLAC
#conf['steps'][1]['datasets']['idis'].append(10016) ## proton   | F2            | BCDMS
#conf['steps'][1]['datasets']['idis'].append(10017) ## deuteron | F2            | BCDMS
#conf['steps'][1]['datasets']['idis'].append(10020) ## proton   | F2            | NMC
#conf['steps'][1]['datasets']['idis'].append(10021) ## d/p      | F2d/F2p       | NMC

# PIDIS
# fit only PPDF
conf['steps'][1] = {}
conf['steps'][1]['dep'] = []
conf['steps'][1]['active distributions'] = ['pdf','ppdf','ppdf_smx']
conf['steps'][1]['passive distributions'] = []
conf['steps'][1]['datasets'] = {}

#conf['steps'][1]['datasets']['idis'] = []
#conf['steps'][1]['datasets']['idis'].append(10010) ## proton   | F2            | SLAC
#conf['steps'][1]['datasets']['idis'].append(10011) ## deuteron | F2            | SLAC
#conf['steps'][1]['datasets']['idis'].append(10016) ## proton   | F2            | BCDMS
#conf['steps'][1]['datasets']['idis'].append(10017) ## deuteron | F2            | BCDMS
#conf['steps'][1]['datasets']['idis'].append(10020) ## proton   | F2            | NMC
#conf['steps'][1]['datasets']['idis'].append(10021) ## d/p      | F2d/F2p       | NMC

conf['steps'][1]['datasets']['pidis'] = []
conf['steps'][1]['datasets']['pidis'].append(10002) # 10002 | proton   | A1   | COMPASS         |          |
conf['steps'][1]['datasets']['pidis'].append(10003) # 10003 | proton   | A1   | COMPASS         |          |
conf['steps'][1]['datasets']['pidis'].append(10004) # 10004 | proton   | A1   | EMC             |          |
conf['steps'][1]['datasets']['pidis'].append(10007) # 10007 | proton   | Apa  | HERMES          |          |
#conf['steps'][1]['datasets']['pidis'].append(10008) # 10008 | proton   | A2   | HERMES          |          |
conf['steps'][1]['datasets']['pidis'].append(10017) # 10017 | proton   | Apa  | JLabHB(EG1DVCS) |          |
conf['steps'][1]['datasets']['pidis'].append(10022) # 10022 | proton   | Apa  | SLAC(E143)      |          |
#conf['steps'][1]['datasets']['pidis'].append(10023) # 10023 | proton   | Ape  | SLAC(E143)      |          |
#conf['steps'][1]['datasets']['pidis'].append(10028) # 10028 | proton   | Ape  | SLAC(E155)      |          |
conf['steps'][1]['datasets']['pidis'].append(10029) # 10029 | proton   | Apa  | SLAC(E155)      |          |
#conf['steps'][1]['datasets']['pidis'].append(10031) # 10031 | proton   | Atpe | SLAC(E155x)     |          |
conf['steps'][1]['datasets']['pidis'].append(10032) # 10032 | proton   | Apa  | SLACE80E130     |          |
conf['steps'][1]['datasets']['pidis'].append(10035) # 10035 | proton   | A1   | SMC             |          |
conf['steps'][1]['datasets']['pidis'].append(10036) # 10036 | proton   | A1   | SMC             |          |
conf['steps'][1]['datasets']['pidis'].append(10041) # 10041 | proton   | Apa  | JLabHB(EG1b)    | E =1 GeV |
conf['steps'][1]['datasets']['pidis'].append(10042) # 10042 | proton   | Apa  | JLabHB(EG1b)    | E =2 GeV |
conf['steps'][1]['datasets']['pidis'].append(10043) # 10043 | proton   | Apa  | JLabHB(EG1b)    | E =4 GeV |
conf['steps'][1]['datasets']['pidis'].append(10044) # 10044 | proton   | Apa  | JLabHB(EG1b)    | E =5 GeV |
conf['steps'][1]['datasets']['pidis'].append(10005) # 10005 | neutron  | A1   | HERMES          |          |
conf['steps'][1]['datasets']['pidis'].append(10001) # 10001 | deuteron | A1   | COMPASS         |          |
conf['steps'][1]['datasets']['pidis'].append(10006) # 10006 | deuteron | Apa  | HERMES          |          |
conf['steps'][1]['datasets']['pidis'].append(10016) # 10016 | deuteron | Apa  | JLabHB(EG1DVCS) |          |
#conf['steps'][1]['datasets']['pidis'].append(10020) # 10020 | deuteron | Ape  | SLAC(E143)      |          |
conf['steps'][1]['datasets']['pidis'].append(10021) # 10021 | deuteron | Apa  | SLAC(E143)      |          |
#conf['steps'][1]['datasets']['pidis'].append(10026) # 10026 | deuteron | Ape  | SLAC(E155)      |          |
conf['steps'][1]['datasets']['pidis'].append(10027) # 10027 | deuteron | Apa  | SLAC(E155)      |          |
#conf['steps'][1]['datasets']['pidis'].append(10030) # 10030 | deuteron | Atpe | SLAC(E155x)     |          |
conf['steps'][1]['datasets']['pidis'].append(10033) # 10033 | deuteron | A1   | SMC             |          |
conf['steps'][1]['datasets']['pidis'].append(10034) # 10034 | deuteron | A1   | SMC             |          |
conf['steps'][1]['datasets']['pidis'].append(10037) # 10037 | deuteron | Apa  | JLabHB(EG1b)    | E =1 GeV |
conf['steps'][1]['datasets']['pidis'].append(10038) # 10038 | deuteron | Apa  | JLabHB(EG1b)    | E =2 GeV |
conf['steps'][1]['datasets']['pidis'].append(10039) # 10039 | deuteron | Apa  | JLabHB(EG1b)    | E =4 GeV |
conf['steps'][1]['datasets']['pidis'].append(10040) # 10040 | deuteron | Apa  | JLabHB(EG1b)    | E =5 GeV |

conf['steps'][1]['datasets']['pidis'].append(70000) # 70000 | proton | Apa  | EIC pseudo-data    |
