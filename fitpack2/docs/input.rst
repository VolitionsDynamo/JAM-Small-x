Input file
==========

In this page we provide a template to build an ``input.py`` 
file to be read by `resman.py`



Headers
-------

.. code-block:: python

   import os
   conf={}
   
   
   #--setup qcd evolution
   
   conf['dglap mode']='truncated'
   conf['alphaSmode']='backward'
   conf['order'] = 'NLO'
   conf['Q20']   = 1.27**2


   #--setup posterior sampling
   
   conf['bootstrap']=False
   conf['flat par']=False
   conf['ftol']=1e-8
   
   #--setup for idis
   
   conf['tmc']   = False
   conf['ht']    = False
   conf['nuc smearing']=True
   conf['sidis nuc smearing']=False
   conf['hq']=False
   
   #--grids
   
   conf['path2idistab'] = '%s/grids/grids-idis/distab/'%os.environ['FITPACK']
   conf['path2dytab']   = '%s/grids/grids-dy/'%os.environ['FITPACK']
   conf['path2jettab']   = '%s/grids/grids-jets/'%os.environ['FITPACK']
   conf['path2pidistab']   = '%s/grids/grids-pidis/'%os.environ['FITPACK']



Datasets
--------

First we create a key for all the datasets 

.. code-block:: python
   
   conf['datasets']={}
 

Unpolarized DIS
...............

.. code-block:: python

   #--lepton-hadron reactions
   
   Q2cut=1.3**2
   W2cut=4.0
   pt_cut = 4.0
   
   ##--IDIS
   conf['datasets']['idis']={}
   conf['datasets']['idis']['filters']=[]
   conf['datasets']['idis']['filters'].append("Q2>%f"%Q2cut)
   conf['datasets']['idis']['filters'].append("W2>%f"%W2cut)
   conf['datasets']['idis']['xlsx']={}
   conf['datasets']['idis']['xlsx'][10010]='idis/expdata/10010.xlsx' # proton   | F2            | SLAC
   conf['datasets']['idis']['xlsx'][10011]='idis/expdata/10011.xlsx' # deuteron | F2            | SLAC
   conf['datasets']['idis']['xlsx'][10016]='idis/expdata/10016.xlsx' # proton   | F2            | BCDMS
   conf['datasets']['idis']['xlsx'][10017]='idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS
   conf['datasets']['idis']['xlsx'][10020]='idis/expdata/10020.xlsx' # proton   | F2            | NMC
   conf['datasets']['idis']['xlsx'][10021]='idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC
   conf['datasets']['idis']['xlsx'][10026]='idis/expdata/10026.xlsx' # proton   | sigma red     | HERA II NC e+ (1)
   conf['datasets']['idis']['xlsx'][10027]='idis/expdata/10027.xlsx' # proton   | sigma red     | HERA II NC e+ (2)
   conf['datasets']['idis']['xlsx'][10028]='idis/expdata/10028.xlsx' # proton   | sigma red     | HERA II NC e+ (3)
   conf['datasets']['idis']['xlsx'][10029]='idis/expdata/10029.xlsx' # proton   | sigma red     | HERA II NC e+ (4)
   conf['datasets']['idis']['xlsx'][10030]='idis/expdata/10030.xlsx' # proton   | sigma red     | HERA II NC e-
   conf['datasets']['idis']['xlsx'][10031]='idis/expdata/10031.xlsx' # proton   | sigma red     | HERA II CC e+
   conf['datasets']['idis']['xlsx'][10032]='idis/expdata/10032.xlsx' # proton   | sigma red     | HERA II NC e-
   conf['datasets']['idis']['norm']={}
   conf['datasets']['idis']['norm'][10010]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10011]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10016]={'value':    9, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10017]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10020]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}

Polarized DIS
.............

.. code-block:: python
   
   conf['datasets']['pidis']={}
   conf['datasets']['pidis']['filters']=[]
   conf['datasets']['pidis']['filters'].append("Q2>%f"%Q2cut)
   conf['datasets']['pidis']['filters'].append("W2>%f"%W2cut)
   conf['datasets']['pidis']['xlsx']={}
   conf['datasets']['pidis']['xlsx'][10001]='pidis/expdata/10001.xlsx' # 10001 | deuteron | A1   | COMPASS         |          
   conf['datasets']['pidis']['xlsx'][10006]='pidis/expdata/10006.xlsx' # 10006 | deuteron | Apa  | HERMES          |          
   conf['datasets']['pidis']['xlsx'][10020]='pidis/expdata/10020.xlsx' # 10020 | deuteron | Ape  | SLAC(E143)      |          
   conf['datasets']['pidis']['xlsx'][10021]='pidis/expdata/10021.xlsx' # 10021 | deuteron | Apa  | SLAC(E143)      |          
   conf['datasets']['pidis']['xlsx'][10026]='pidis/expdata/10026.xlsx' # 10026 | deuteron | Ape  | SLAC(E155)      |          
   conf['datasets']['pidis']['xlsx'][10027]='pidis/expdata/10027.xlsx' # 10027 | deuteron | Apa  | SLAC(E155)      |          
   #conf['datasets']['pidis']['xlsx'][10030]='pidis/expdata/10030.xlsx' # 10030 | deuteron | Atpe | SLAC(E155x)     |          
   conf['datasets']['pidis']['xlsx'][10033]='pidis/expdata/10033.xlsx' # 10033 | deuteron | A1   | SMC             |          
   conf['datasets']['pidis']['xlsx'][10034]='pidis/expdata/10034.xlsx' # 10034 | deuteron | A1   | SMC             |          
   conf['datasets']['pidis']['xlsx'][10002]='pidis/expdata/10002.xlsx' # 10002 | proton   | A1   | COMPASS         |          
   conf['datasets']['pidis']['xlsx'][10003]='pidis/expdata/10003.xlsx' # 10003 | proton   | A1   | COMPASS         |          
   conf['datasets']['pidis']['xlsx'][10004]='pidis/expdata/10004.xlsx' # 10004 | proton   | A1   | EMC             |          
   conf['datasets']['pidis']['xlsx'][10007]='pidis/expdata/10007.xlsx' # 10007 | proton   | Apa  | HERMES          |          
   conf['datasets']['pidis']['xlsx'][10008]='pidis/expdata/10008.xlsx' # 10008 | proton   | A2   | HERMES          |          
   conf['datasets']['pidis']['xlsx'][10022]='pidis/expdata/10022.xlsx' # 10022 | proton   | Apa  | SLAC(E143)      |          
   conf['datasets']['pidis']['xlsx'][10023]='pidis/expdata/10023.xlsx' # 10023 | proton   | Ape  | SLAC(E143)      |          
   conf['datasets']['pidis']['xlsx'][10028]='pidis/expdata/10028.xlsx' # 10028 | proton   | Ape  | SLAC(E155)      |          
   conf['datasets']['pidis']['xlsx'][10029]='pidis/expdata/10029.xlsx' # 10029 | proton   | Apa  | SLAC(E155)      |          
   #conf['datasets']['pidis']['xlsx'][10031]='pidis/expdata/10031.xlsx' # 10031 | proton   | Atpe | SLAC(E155x)     |          
   conf['datasets']['pidis']['xlsx'][10032]='pidis/expdata/10032.xlsx' # 10032 | proton   | Apa  | SLACE80E130     |          
   conf['datasets']['pidis']['xlsx'][10035]='pidis/expdata/10035.xlsx' # 10035 | proton   | A1   | SMC             |          
   conf['datasets']['pidis']['xlsx'][10036]='pidis/expdata/10036.xlsx' # 10036 | proton   | A1   | SMC             |          
   conf['datasets']['pidis']['xlsx'][10005]='pidis/expdata/10005.xlsx' # 10005 | neutron  | A1   | HERMES          |          
   #conf['datasets']['pidis']['xlsx'][10018]='pidis/expdata/10018.xlsx' # 10018 | helium   | A1   | SLAC(E142)      |          
   #conf['datasets']['pidis']['xlsx'][10019]='pidis/expdata/10019.xlsx' # 10019 | helium   | A2   | SLAC(E142)      |          
   #conf['datasets']['pidis']['xlsx'][10024]='pidis/expdata/10024.xlsx' # 10024 | helium   | Ape  | SLAC(E154)      |          
   #conf['datasets']['pidis']['xlsx'][10025]='pidis/expdata/10025.xlsx' # 10025 | helium   | Apa  | SLAC(E154)      |          
   conf['datasets']['pidis']['norm']={}
   conf['datasets']['pidis']['norm'][10001]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10006]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10020]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10021]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10026]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10027]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   #conf['datasets']['pidis']['norm'][10030]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10033]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10034]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10002]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10003]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10004]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10007]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10008]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10022]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10023]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10028]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10029]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   #conf['datasets']['pidis']['norm'][10031]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10032]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10035]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10036]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['pidis']['norm'][10005]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   #conf['datasets']['pidis']['norm'][10018]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   #conf['datasets']['pidis']['norm'][10019]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   #conf['datasets']['pidis']['norm'][10024]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   #conf['datasets']['pidis']['norm'][10025]={'value':    1.00000, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   

Unpolarized SIDIS
.................

.. code-block:: python
   
   conf['datasets']['sidis']={}
   conf['datasets']['sidis']['filters']=[]
   conf['datasets']['sidis']['filters'].append("Q2>%f"%Q2cut) 
   conf['datasets']['sidis']['filters'].append("W2>%f"%W2cut) 
   conf['datasets']['sidis']['filters'].append('Z>0.2 and Z<0.8')
   conf['datasets']['sidis']['xlsx']={}
   conf['datasets']['sidis']['xlsx'][1005]='sidis/expdata/1005.xlsx' # deuteron , mult , pi+ , COMPASS
   conf['datasets']['sidis']['xlsx'][1006]='sidis/expdata/1006.xlsx' # deuteron , mult , pi- , COMPASS
   conf['datasets']['sidis']['xlsx'][2005]='sidis/expdata/2005.xlsx' # deuteron , mult , K+  , COMPASS
   conf['datasets']['sidis']['xlsx'][2006]='sidis/expdata/2006.xlsx' # deuteron , mult , K-  , COMPASS
   conf['datasets']['sidis']['xlsx'][1007]='sidis/expdata/1007.xlsx' # proton   , mult , pi+ , HERMES, Q2-z
   conf['datasets']['sidis']['xlsx'][1008]='sidis/expdata/1008.xlsx' # proton   , mult , pi- , HERMES, Q2-z
   conf['datasets']['sidis']['xlsx'][1009]='sidis/expdata/1009.xlsx' # deuteron , mult , pi+ , HERMES, Q2-z
   conf['datasets']['sidis']['xlsx'][1010]='sidis/expdata/1010.xlsx' # deuteron , mult , pi- , HERMES, Q2-z
   conf['datasets']['sidis']['xlsx'][2007]='sidis/expdata/2007.xlsx' # proton   , mult , K+  , HERMES, Q2-z
   conf['datasets']['sidis']['xlsx'][2008]='sidis/expdata/2008.xlsx' # proton   , mult , K-  , HERMES, Q2-z
   conf['datasets']['sidis']['xlsx'][2009]='sidis/expdata/2009.xlsx' # deuteron , mult , K+  , HERMES, Q2-z
   conf['datasets']['sidis']['xlsx'][2010]='sidis/expdata/2010.xlsx' # deuteron , mult , K-  , HERMES, Q2-z
   conf['datasets']['sidis']['norm']={}
   conf['datasets']['sidis']['norm'][1005]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][1006]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][1007]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][1008]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][1009]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][1010]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][2005]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][2006]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][2007]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][2008]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][2009]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sidis']['norm'][2010]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}

Polarized SIDIS
.................

.. code-block:: python
   
   conf['datasets']['psidis']={}
   conf['datasets']['psidis']['filters']=[]
   conf['datasets']['psidis']['filters'].append("Q2>%f"%Q2cut) 
   conf['datasets']['psidis']['filters'].append("W2>%f"%W2cut) 
   conf['datasets']['psidis']['filters'].append('Z>0.2 and Z<0.8')
   conf['datasets']['psidis']['xlsx']={}
   conf['datasets']['psidis']['xlsx'][20004]='psidis/expdata/20004.xlsx' # 20004 | proton   | A1pi+  | HERMES  
   conf['datasets']['psidis']['xlsx'][20005]='psidis/expdata/20005.xlsx' # 20005 | proton   | A1pi-  | HERMES  
   conf['datasets']['psidis']['xlsx'][20008]='psidis/expdata/20008.xlsx' # 20008 | deuteron | A1pi+  | HERMES  
   conf['datasets']['psidis']['xlsx'][20009]='psidis/expdata/20009.xlsx' # 20009 | deuteron | A1pi-  | HERMES  
   conf['datasets']['psidis']['xlsx'][20012]='psidis/expdata/20012.xlsx' # 20012 | deuteron | A1K+   | HERMES  
   conf['datasets']['psidis']['xlsx'][20013]='psidis/expdata/20013.xlsx' # 20013 | deuteron | A1K-   | HERMES  
   conf['datasets']['psidis']['xlsx'][20014]='psidis/expdata/20014.xlsx' # 20014 | deuteron | A1Ksum | HERMES  
   conf['datasets']['psidis']['xlsx'][20017]='psidis/expdata/20017.xlsx' # 20017 | proton   | A1pi+  | COMPASS 
   conf['datasets']['psidis']['xlsx'][20018]='psidis/expdata/20018.xlsx' # 20018 | proton   | A1pi-  | COMPASS 
   conf['datasets']['psidis']['xlsx'][20019]='psidis/expdata/20019.xlsx' # 20019 | proton   | A1K+   | COMPASS 
   conf['datasets']['psidis']['xlsx'][20020]='psidis/expdata/20020.xlsx' # 20020 | proton   | A1K-   | COMPASS 
   conf['datasets']['psidis']['xlsx'][20021]='psidis/expdata/20021.xlsx' # 20021 | deuteron | A1pi+  | COMPASS 
   conf['datasets']['psidis']['xlsx'][20022]='psidis/expdata/20022.xlsx' # 20022 | deuteron | A1pi-  | COMPASS 
   conf['datasets']['psidis']['xlsx'][20025]='psidis/expdata/20025.xlsx' # 20025 | deuteron | A1K+   | COMPASS 
   conf['datasets']['psidis']['xlsx'][20026]='psidis/expdata/20026.xlsx' # 20026 | deuteron | A1K-   | COMPASS 
   conf['datasets']['psidis']['norm']={}
   conf['datasets']['psidis']['norm'][20004]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20005]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20008]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20009]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20012]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20013]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20014]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20017]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20018]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20019]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20020]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20021]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20022]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20025]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['psidis']['norm'][20026]={'value':    1.000000,'fixed':False,'min':0.5,'max':1.5}

Drell-Yan
.........

.. code-block:: python
   
   conf['datasets']['dy']={}
   conf['datasets']['dy']['filters']=[]
   conf['datasets']['dy']['filters'].append("Q2>36") 
   conf['datasets']['dy']['xlsx']={}
   conf['datasets']['dy']['xlsx'][10001]='dy/expdata/10001.xlsx'
   conf['datasets']['dy']['xlsx'][10002]='dy/expdata/10002.xlsx'
   conf['datasets']['dy']['norm']={}
   conf['datasets']['dy']['norm'][10001]={'value':    1,'fixed':False,'min':   0.8,'max':    1.2}
   conf['datasets']['dy']['norm'][10002]={'value':    1,'fixed':False,'min':   0.8,'max':    1.2}
   
SIA pion
........

.. code-block:: python
   
   conf['datasets']['sia']={}
   conf['datasets']['sia']['filters']=[]
   conf['datasets']['sia']['filters'].append('z>0.2 and z<0.9') 
   conf['datasets']['sia']['xlsx']={}
   conf['datasets']['sia']['xlsx'][1001]='sia/expdata/1001.xlsx'  # hadron: pion exp: TASSO
   conf['datasets']['sia']['xlsx'][1002]='sia/expdata/1002.xlsx'  # hadron: pion exp: TASSO
   conf['datasets']['sia']['xlsx'][1003]='sia/expdata/1003.xlsx'  # hadron: pion exp: TASSO
   conf['datasets']['sia']['xlsx'][1004]='sia/expdata/1004.xlsx'  # hadron: pion exp: TASSO
   conf['datasets']['sia']['xlsx'][1005]='sia/expdata/1005.xlsx'  # hadron: pion exp: TASSO
   conf['datasets']['sia']['xlsx'][1006]='sia/expdata/1006.xlsx'  # hadron: pion exp: TASSO
   conf['datasets']['sia']['xlsx'][1007]='sia/expdata/1007.xlsx'  # hadron: pion exp: TPC
   conf['datasets']['sia']['xlsx'][1008]='sia/expdata/1008.xlsx'  # hadron: pion exp: TPC
   conf['datasets']['sia']['xlsx'][1012]='sia/expdata/1012.xlsx'  # hadron: pion exp: HRS
   conf['datasets']['sia']['xlsx'][1013]='sia/expdata/1013.xlsx'  # hadron: pion exp: TOPAZ
   conf['datasets']['sia']['xlsx'][1014]='sia/expdata/1014.xlsx'  # hadron: pion exp: SLD
   conf['datasets']['sia']['xlsx'][1018]='sia/expdata/1018.xlsx'  # hadron: pion exp: ALEPH
   conf['datasets']['sia']['xlsx'][1019]='sia/expdata/1019.xlsx'  # hadron: pion exp: OPAL
   conf['datasets']['sia']['xlsx'][1025]='sia/expdata/1025.xlsx'  # hadron: pion exp: DELPHI
   conf['datasets']['sia']['xlsx'][1028]='sia/expdata/1028.xlsx'  # hadron: pion exp: BABAR
   conf['datasets']['sia']['xlsx'][1029]='sia/expdata/1029.xlsx'  # hadron: pion exp: BELL
   conf['datasets']['sia']['xlsx'][1030]='sia/expdata/1030.xlsx'  # hadron: pion exp: ARGUS
   conf['datasets']['sia']['norm']={}
   conf['datasets']['sia']['norm'][1001]={'value':    1.10478e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1002]={'value':    9.82581e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1003]={'value':    1.03054e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1005]={'value':    1.03419e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1006]={'value':    9.79162e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1014]={'value':    9.97770e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1018]={'value':    1.02378e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1028]={'value':    9.76001e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1029]={'value':    8.68358e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1030]={'value':    1.01509e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1007]={'value':    1.06806e+00,'fixed':False,'min':0.5,'max':1.5,'dN':0.1}
   conf['datasets']['sia']['norm'][1008]={'value':    1.06603e+00,'fixed':False,'min':0.5,'max':1.5,'dN':0.1}
   conf['datasets']['sia']['norm'][1019]={'value':    1.03100e+00,'fixed':False,'min':0.5,'max':1.5,'dN':0.1}
   
   ##--SIA pion HQ
   conf['datasets']['sia']['xlsx'][1010]='sia/expdata/1010.xlsx'  # hadron: pion exp: TPC(c)
   conf['datasets']['sia']['xlsx'][1011]='sia/expdata/1011.xlsx'  # hadron: pion exp: TPC(b)
   conf['datasets']['sia']['xlsx'][1016]='sia/expdata/1016.xlsx'  # hadron: pion exp: SLD(c)
   conf['datasets']['sia']['xlsx'][1017]='sia/expdata/1017.xlsx'  # hadron: pion exp: SLD(b)
   conf['datasets']['sia']['xlsx'][1023]='sia/expdata/1023.xlsx'  # hadron: pion exp: OPAL(c)
   conf['datasets']['sia']['xlsx'][1024]='sia/expdata/1024.xlsx'  # hadron: pion exp: OPAL(b)
   conf['datasets']['sia']['xlsx'][1027]='sia/expdata/1027.xlsx'  # hadron: pion exp: DELPHI(b)
   conf['datasets']['sia']['norm'][1016]={'value':    1.18920e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1017]={'value':    1.00345e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1023]={'value':    1.33434e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][1024]={'value':    1.19151e+00,'fixed':False,'min':0.5,'max':1.5}

SIA Kaon
........

.. code-block:: python
   
   conf['datasets']['sia']['xlsx'][2030]='sia/expdata/2030.xlsx'  # hadron: kaon exp: ARGUS
   conf['datasets']['sia']['xlsx'][2028]='sia/expdata/2028.xlsx'  # hadron: kaon exp: BABAR
   conf['datasets']['sia']['xlsx'][2029]='sia/expdata/2029.xlsx'  # hadron: kaon exp: BELL
   conf['datasets']['sia']['xlsx'][2001]='sia/expdata/2001.xlsx'  # hadron: kaon exp: TASSO
   conf['datasets']['sia']['xlsx'][2002]='sia/expdata/2002.xlsx'  # hadron: kaon exp: TASSO
   conf['datasets']['sia']['xlsx'][2003]='sia/expdata/2003.xlsx'  # hadron: kaon exp: TASSO
   conf['datasets']['sia']['xlsx'][2004]='sia/expdata/2004.xlsx'  # hadron: kaon exp: TASSO
   conf['datasets']['sia']['xlsx'][2005]='sia/expdata/2005.xlsx'  # hadron: kaon exp: TASSO
   conf['datasets']['sia']['xlsx'][2006]='sia/expdata/2006.xlsx'  # hadron: kaon exp: TASSO
   conf['datasets']['sia']['xlsx'][2007]='sia/expdata/2007.xlsx'  # hadron: kaon exp: TPC
   conf['datasets']['sia']['xlsx'][2008]='sia/expdata/2008.xlsx'  # hadron: kaon exp: TPC
   conf['datasets']['sia']['xlsx'][2012]='sia/expdata/2012.xlsx'  # hadron: kaon exp: HRS
   conf['datasets']['sia']['xlsx'][2013]='sia/expdata/2013.xlsx'  # hadron: kaon exp: TOPAZ
   conf['datasets']['sia']['xlsx'][2014]='sia/expdata/2014.xlsx'  # hadron: kaon exp: SLD
   conf['datasets']['sia']['xlsx'][2018]='sia/expdata/2018.xlsx'  # hadron: kaon exp: ALEPH
   conf['datasets']['sia']['xlsx'][2019]='sia/expdata/2019.xlsx'  # hadron: kaon exp: OPAL
   conf['datasets']['sia']['xlsx'][2025]='sia/expdata/2025.xlsx'  # hadron: kaon exp: DELPHI
   conf['datasets']['sia']['xlsx'][2031]='sia/expdata/2031.xlsx'  # hadron: kaon exp: DELPHI
   conf['datasets']['sia']['norm'][2030]={'value':    1.00482e+00,'fixed':False,'min':0.5,'max':1.5,'dN':0.1}
   conf['datasets']['sia']['norm'][2028]={'value':    9.97435e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2029]={'value':    1.00000e+00,'fixed':True,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2002]={'value':    9.83394e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2003]={'value':    9.94421e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2005]={'value':    9.92876e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2014]={'value':    9.12186e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2018]={'value':    8.88453e-01,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2019]={'value':    8.40305e-01,'fixed':False,'min':0.5,'max':1.5,'dN':0.1}
   conf['datasets']['sia']['norm'][2007]={'value':    1.20488e+00,'fixed':False,'min':0.5,'max':1.5,'dN':0.1}
   conf['datasets']['sia']['norm'][2008]={'value':    1.10314e+00,'fixed':False,'min':0.5,'max':1.5,'dN':0.1}
   
   ##--SIA kaon HQ
   conf['datasets']['sia']['xlsx'][2016]='sia/expdata/2016.xlsx'  # hadron: kaon exp: SLD(c)
   conf['datasets']['sia']['xlsx'][2017]='sia/expdata/2017.xlsx'  # hadron: kaon exp: SLD(b)
   conf['datasets']['sia']['xlsx'][2023]='sia/expdata/2023.xlsx'  # hadron: kaon exp: OPAL(c)
   conf['datasets']['sia']['xlsx'][2024]='sia/expdata/2024.xlsx'  # hadron: kaon exp: OPAL(b)
   conf['datasets']['sia']['xlsx'][2027]='sia/expdata/2027.xlsx'  # hadron: kaon exp: DELPHI(b)
   conf['datasets']['sia']['norm'][2016]={'value':    1.03996e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2017]={'value':    1.00015e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2023]={'value':    1.35922e+00,'fixed':False,'min':0.5,'max':1.5}
   conf['datasets']['sia']['norm'][2024]={'value':    1.38163e+00,'fixed':False,'min':0.5,'max':1.5}

Unpolarized JETs
................

.. code-block:: python
   
   conf['datasets']['jet'] = {}
   conf['datasets']['jet']['xlsx'] = {}
   conf['datasets']['jet']['xlsx'][10001] = 'jets/expdata/10001.xlsx' ## D0 dataset
   conf['datasets']['jet']['xlsx'][10002] = 'jets/expdata/10002.xlsx' ## CDF dataset
   conf['datasets']['jet']['xlsx'][10003] = 'jets/expdata/10003.xlsx' ## STAR MB dataset
   conf['datasets']['jet']['xlsx'][10004] = 'jets/expdata/10004.xlsx' ## STAR HT dataset
   conf['datasets']['jet']['norm'] = {}
   conf['datasets']['jet']['norm'][10001] = {'value': 1, 'fixed': True, 'min': 0.8, 'max': 1.2}
   conf['datasets']['jet']['norm'][10002] = {'value': 1, 'fixed': False, 'min': 0.8, 'max': 1.2}
   conf['datasets']['jet']['norm'][10003] = {'value': 1, 'fixed': True, 'min': 0.8, 'max': 1.2}
   conf['datasets']['jet']['norm'][10004] = {'value': 1, 'fixed': True, 'min': 0.8, 'max': 1.2}
   conf['channels'] = {'q_qp': 0, 'q_qbp': 1, 'q_q': 2, 'q_qb': 3, 'q_g': 4, 'g_g': 5}
   conf['qr'] = {}
   conf['qr']['parameter'] = {'epsilon': 1e-11, 'block_size': 10, 'power_scheme': 1}
   conf['fit'] = {'method': 'fixed', 'f_scale': 0.25, 'r_scale': 0.25}

Polarized JETs
................

.. code-block:: python

   conf['datasets']['pjet'] = {}
   conf['datasets']['pjet']['filters'] = []
   conf['datasets']['pjet']['filters'].append("pT>%f" % pt_cut)
   conf['datasets']['pjet']['xlsx'] = {}
   conf['datasets']['pjet']['xlsx'][20001] = 'pjets/expdata/20001.xlsx' ## STAR 2006 paper on 2003 and 2004 data
   conf['datasets']['pjet']['xlsx'][20002] = 'pjets/expdata/20002.xlsx' ## STAR 2012 paper on 2005 data
   conf['datasets']['pjet']['xlsx'][20003] = 'pjets/expdata/20003.xlsx' ## STAR 2012 paper on 2006 data
   conf['datasets']['pjet']['xlsx'][20004] = 'pjets/expdata/20004.xlsx' ## STAR 2015 paper on 2009 data
   conf['datasets']['pjet']['xlsx'][20005] = 'pjets/expdata/20005.xlsx' ## PHENIX 2011 paper on 2005 data
   conf['datasets']['pjet']['xlsx'][99999] = 'pjets/expdata/99999.xlsx' ## STAR 2019 paper on 2012 data
   conf['datasets']['pjet']['norm'] = {}
   conf['datasets']['pjet']['norm'][20001] = {'value': 1, 'fixed': True, 'min': 0.8, 'max': 1.2}
   conf['datasets']['pjet']['norm'][20002] = {'value': 1, 'fixed': False, 'min': 0.8, 'max': 1.2}
   conf['datasets']['pjet']['norm'][20003] = {'value': 1, 'fixed': False, 'min': 0.8, 'max': 1.2}
   conf['datasets']['pjet']['norm'][20004] = {'value': 1, 'fixed': False, 'min': 0.8, 'max': 1.2}
   # conf['datasets']['pjet']['norm'][20005] = {'value': 1, 'fixed': False, 'min': 0.8, 'max': 1.2}
   conf['datasets']['pjet']['norm'][20006] = {'value': 1, 'fixed': False, 'min': 0.8, 'max': 1.2}
   conf['channels'] = {'q_qp': 0, 'q_qbp': 1, 'q_q': 2, 'q_qb': 3, 'q_g': 4, 'g_g': 5}
   conf['qr'] = {}
   conf['qr']['parameter'] = {'epsilon': 1e-11, 'block_size': 10, 'power_scheme': 1}
   conf['fit'] = {'method': 'fixed', 'f_scale': 0.25, 'r_scale': 0.25}


Parameters
----------

.. code-block:: python

   #--parameters
   conf['params'] = {}

pdfs
....

.. code-block:: python
   
   conf['pdf parametrization'] = 2
   conf['params']['pdf'] = {}
   
   conf['params']['pdf']['g1 N']    ={'value':    3.87592e-01, 'min':  None, 'max':  None, 'fixed': True }
   conf['params']['pdf']['g1 a']    ={'value':   -6.60217e-01, 'min':  -1.9, 'max':     1, 'fixed': False}
   conf['params']['pdf']['g1 b']    ={'value':    8.12091e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['uv1 N']   ={'value':    3.47549e-01, 'min':  None, 'max':  None, 'fixed': True }
   conf['params']['pdf']['uv1 a']   ={'value':   -1.14055e-01, 'min':  -0.5, 'max':     1, 'fixed': False}
   conf['params']['pdf']['uv1 b']   ={'value':    3.21230e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['dv1 N']   ={'value':    1.52089e-01, 'min':  None, 'max':  None, 'fixed': True }
   conf['params']['pdf']['dv1 a']   ={'value':   -3.80099e-02, 'min':  -0.5, 'max':     1, 'fixed': False}
   conf['params']['pdf']['dv1 b']   ={'value':    4.36319e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['db1 N']   ={'value':    2.08640e-02, 'min':     0, 'max':     1, 'fixed': False}
   conf['params']['pdf']['db1 a']   ={'value':   -4.05372e-01, 'min':    -1, 'max':     1, 'fixed': False}
   conf['params']['pdf']['db1 b']   ={'value':    1.00000e+01, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['ub1 N']   ={'value':    1.64626e-02, 'min':     0, 'max':     1, 'fixed': False}
   conf['params']['pdf']['ub1 a']   ={'value':    0.000000000, 'min':    -1, 'max':     1, 'fixed': False}
   conf['params']['pdf']['ub1 b']   ={'value':    1.00000e+01, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['s1 N']    ={'value':    0.00000e+00, 'min':     0, 'max':     1, 'fixed': True }
   conf['params']['pdf']['s1 a']    ={'value':   -0.19319e+00, 'min':  -0.5, 'max':     1, 'fixed': False}
   conf['params']['pdf']['s1 b']    ={'value':    4.99999e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['sb1 N']   ={'value':    0.00000e+00, 'min':     0, 'max':     1, 'fixed': False}
   conf['params']['pdf']['sb1 a']   ={'value':   -0.19317e+00, 'min':  -0.5, 'max':     1, 'fixed': False}
   conf['params']['pdf']['sb1 b']   ={'value':    4.99988e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['sea1 N']  ={'value':    2.08640e-02, 'min':     0, 'max':     1, 'fixed': False}
   conf['params']['pdf']['sea1 a']  ={'value':   -1.500000000, 'min':  -1.9, 'max':    -1, 'fixed': False}
   conf['params']['pdf']['sea1 b']  ={'value':    1.00000e+01, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['sea2 N']  ={'value':    2.08640e-02, 'min':     0, 'max':     1, 'fixed': 'sea1 N'}
   conf['params']['pdf']['sea2 a']  ={'value':   -1.500000000, 'min':  -1.9, 'max':    -1, 'fixed': 'sea1 a'}
   conf['params']['pdf']['sea2 b']  ={'value':    1.00000e+01, 'min':     0, 'max':    10, 'fixed': 'sea1 b'}
   
spin pdfs
.........

.. code-block:: python
   
   conf['ppdf parametrization']=0
   conf['su2+su3']=False
   conf['params']['ppdf']={}
   
   conf['params']['ppdf']['g1 N']    ={'value':    3.87592e-01, 'min':   -10, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['g1 a']    ={'value':   -6.60217e-01, 'min':    -1, 'max':     2, 'fixed': False}
   conf['params']['ppdf']['g1 b']    ={'value':    8.12091e+00, 'min':     0, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['g1 d']    ={'value':    0.21230e+00, 'min':   -10, 'max':    10, 'fixed': False}
   
   conf['params']['ppdf']['up1 N']   ={'value':    3.47549e-01, 'min':   -10, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['up1 a']   ={'value':   -1.14055e-01, 'min':    -1, 'max':     2, 'fixed': False}
   conf['params']['ppdf']['up1 b']   ={'value':    3.21230e+00, 'min':     0, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['up1 d']   ={'value':    0.21230e+00, 'min':   -10, 'max':    10, 'fixed': False}
   
   conf['params']['ppdf']['dp1 N']   ={'value':    1.52089e-01, 'min':   -10, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['dp1 a']   ={'value':   -3.80099e-02, 'min':    -1, 'max':     2, 'fixed': False}
   conf['params']['ppdf']['dp1 b']   ={'value':    4.36319e+00, 'min':     0, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['dp1 d']   ={'value':    0.36319e+00, 'min':   -10, 'max':    10, 'fixed': False}
   
   conf['params']['ppdf']['sp1 N']   ={'value':    1.52089e-01, 'min':   -10, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['sp1 a']   ={'value':   -3.80099e-02, 'min':    -1, 'max':     2, 'fixed': False}
   conf['params']['ppdf']['sp1 b']   ={'value':    4.36319e+00, 'min':     0, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['sp1 d']   ={'value':    0.36319e+00, 'min':   -10, 'max':    10, 'fixed': False}
   
   conf['params']['ppdf']['um1 N']   ={'value':    3.47549e-01, 'min':   -10, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['um1 a']   ={'value':   -1.14055e-01, 'min':    -1, 'max':     2, 'fixed': False}
   conf['params']['ppdf']['um1 b']   ={'value':    3.21230e+00, 'min':     0, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['um1 d']   ={'value':    0.21230e+00, 'min':   -10, 'max':    10, 'fixed': False}
   
   conf['params']['ppdf']['dm1 N']   ={'value':    1.52089e-01, 'min':   -10, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['dm1 a']   ={'value':   -3.80099e-02, 'min':    -1, 'max':     2, 'fixed': False}
   conf['params']['ppdf']['dm1 b']   ={'value':    4.36319e+00, 'min':     0, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['dm1 d']   ={'value':    0.36319e+00, 'min':   -10, 'max':    10, 'fixed': False}
   
   conf['params']['ppdf']['sm1 N']   ={'value':    1.52089e-01, 'min':   -10, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['sm1 a']   ={'value':   -3.80099e-02, 'min':    -1, 'max':     2, 'fixed': False}
   conf['params']['ppdf']['sm1 b']   ={'value':    4.36319e+00, 'min':     0, 'max':    10, 'fixed': False}
   conf['params']['ppdf']['sm1 d']   ={'value':    0.36319e+00, 'min':   -10, 'max':    10, 'fixed': False}
   
pion fragamentation functions
..............................

.. code-block:: python
   
   conf['ffpion parametrization']=0
   conf['params']['ffpion']={}
   
   conf['params']['ffpion']['g1 N']  ={'value':    2.95437e-01, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffpion']['g1 a']  ={'value':    1.00469e+00, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffpion']['g1 b']  ={'value':    6.85766e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                   
   conf['params']['ffpion']['u1 N']  ={'value':    2.67821e-02, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffpion']['u1 a']  ={'value':    1.76877e-01, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffpion']['u1 b']  ={'value':    4.81521e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                      
   conf['params']['ffpion']['d1 N']  ={'value':    2.99974e-01, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffpion']['d1 a']  ={'value':   -6.89477e-01, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffpion']['d1 b']  ={'value':    4.79992e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                      
   conf['params']['ffpion']['s1 N']  ={'value':    1.54863e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
   conf['params']['ffpion']['s1 a']  ={'value':    3.00305e-01, 'min': -1.8 , 'max':     2,'fixed':'d1 a'}
   conf['params']['ffpion']['s1 b']  ={'value':    1.83178e+00, 'min':  0   , 'max':    10,'fixed':'d1 b'}
                                                                                                      
   conf['params']['ffpion']['c1 N']  ={'value':    1.84550e-01, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffpion']['c1 a']  ={'value':   -5.05798e-02, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffpion']['c1 b']  ={'value':    3.19952e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                      
   conf['params']['ffpion']['b1 N']  ={'value':    3.74125e-01, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffpion']['b1 a']  ={'value':   -1.59541e+00, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffpion']['b1 b']  ={'value':    4.50102e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                      
   conf['params']['ffpion']['ub1 N'] ={'value':    2.99974e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
   conf['params']['ffpion']['ub1 a'] ={'value':   -6.89477e-01, 'min': -1.8 , 'max':     2,'fixed':'d1 a'}
   conf['params']['ffpion']['ub1 b'] ={'value':    4.79992e+00, 'min':  0   , 'max':    10,'fixed':'d1 b'}
                                                                                                     
   conf['params']['ffpion']['db1 N'] ={'value':    2.67821e-02, 'min':  0.0 , 'max':     1,'fixed':'u1 N'}
   conf['params']['ffpion']['db1 a'] ={'value':    1.76877e-01, 'min': -1.8 , 'max':     2,'fixed':'u1 a'}
   conf['params']['ffpion']['db1 b'] ={'value':    4.81521e+00, 'min':  0   , 'max':    10,'fixed':'u1 b'}
                                                                                                     
   conf['params']['ffpion']['sb1 N'] ={'value':    1.54863e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
   conf['params']['ffpion']['sb1 a'] ={'value':    3.00305e-01, 'min': -1.8 , 'max':     2,'fixed':'d1 a'}
   conf['params']['ffpion']['sb1 b'] ={'value':    1.83178e+00, 'min':  0   , 'max':    10,'fixed':'d1 b'}
                                                                                                     
   conf['params']['ffpion']['cb1 N'] ={'value':    1.84550e-01, 'min':  0.0 , 'max':     1,'fixed':'c1 N'}
   conf['params']['ffpion']['cb1 a'] ={'value':   -5.05798e-02, 'min': -1.8 , 'max':     2,'fixed':'c1 a'}
   conf['params']['ffpion']['cb1 b'] ={'value':    3.19952e+00, 'min':  0   , 'max':    10,'fixed':'c1 b'}
   
   conf['params']['ffpion']['bb1 N'] ={'value':    3.74125e-01, 'min':  0.0 , 'max':     1,'fixed':'b1 N'}
   conf['params']['ffpion']['bb1 a'] ={'value':   -1.59541e+00, 'min': -1.8 , 'max':     2,'fixed':'b1 a'}
   conf['params']['ffpion']['bb1 b'] ={'value':    4.50102e+00, 'min':  0   , 'max':    10,'fixed':'b1 b'}
   

kaon fragamentation functions
..............................

.. code-block:: python
   
   conf['params']['ffkaon']={}
   conf['ffkaon parametrization']=0
   
   conf['params']['ffkaon']['g1 N']  ={'value':   2.33320e-01, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffkaon']['g1 a']  ={'value':   1.48737e+00, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffkaon']['g1 b']  ={'value':   9.62755e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                   
   conf['params']['ffkaon']['u1 N']  ={'value':   4.03672e-02, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffkaon']['u1 a']  ={'value':   1.26356e+00, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffkaon']['u1 b']  ={'value':   1.62596e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                      
   conf['params']['ffkaon']['d1 N']  ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffkaon']['d1 a']  ={'value':  -1.43444e+00, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffkaon']['d1 b']  ={'value':   1.65143e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                      
   conf['params']['ffkaon']['s1 N']  ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
   conf['params']['ffkaon']['s1 a']  ={'value':  -1.43444e+00, 'min': -1.8 , 'max':     2,'fixed':'d1 a'}
   conf['params']['ffkaon']['s1 b']  ={'value':   1.65143e+00, 'min':  0   , 'max':    10,'fixed':'d1 b'}
                                                                                                      
   conf['params']['ffkaon']['c1 N']  ={'value':   7.76923e-01, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffkaon']['c1 a']  ={'value':  -1.80000e+00, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffkaon']['c1 b']  ={'value':   2.50452e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                      
   conf['params']['ffkaon']['b1 N']  ={'value':   5.66971e-01, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffkaon']['b1 a']  ={'value':  -1.80000e+00, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffkaon']['b1 b']  ={'value':   3.41727e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                      
   conf['params']['ffkaon']['ub1 N'] ={'value':   4.03672e-02, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
   conf['params']['ffkaon']['ub1 a'] ={'value':   2.26356e+00, 'min': -1.8 , 'max':     2,'fixed':'d1 a'}
   conf['params']['ffkaon']['ub1 b'] ={'value':   1.62596e+00, 'min':  0   , 'max':    10,'fixed':'d1 b'}
                                                                                                  
   conf['params']['ffkaon']['db1 N'] ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
   conf['params']['ffkaon']['db1 a'] ={'value':  -1.43444e+00, 'min': -1.8 , 'max':     2,'fixed':'d1 a'}
   conf['params']['ffkaon']['db1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    10,'fixed':'d1 b'}
                                                                                                  
   conf['params']['ffkaon']['sb1 N'] ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':False}
   conf['params']['ffkaon']['sb1 a'] ={'value':  -1.43444e+00, 'min': -1.8 , 'max':     2,'fixed':False}
   conf['params']['ffkaon']['sb1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    10,'fixed':False}
                                                                                                     
   conf['params']['ffkaon']['cb1 N'] ={'value':   7.76923e-01, 'min':  0.0 , 'max':     1,'fixed':'c1 N'}
   conf['params']['ffkaon']['cb1 a'] ={'value':  -1.80000e+00, 'min': -1.8 , 'max':     2,'fixed':'c1 a'}
   conf['params']['ffkaon']['cb1 b'] ={'value':   2.50452e+00, 'min':  0   , 'max':    10,'fixed':'c1 b'}
   
   conf['params']['ffkaon']['bb1 N'] ={'value':   5.66971e-01, 'min':  0.0 , 'max':     1,'fixed':'b1 N'}
   conf['params']['ffkaon']['bb1 a'] ={'value':  -1.80000e+00, 'min': -1.8 , 'max':     2,'fixed':'b1 a'}
   conf['params']['ffkaon']['bb1 b'] ={'value':   3.41727e+00, 'min':  0   , 'max':    10,'fixed':'b1 b'}



Steps
-----

We provide here examples to define a step in the multistep strategy.
The steps are currentlty interpretable by the script ``maxlike.py``. 
Some explanation are in order:

- **dep**: a  list of the steps that are used as priors. 
- **active distributions**: a list of distributions like pdfs or ffs that are needed
- **passive distributions**: remove from the active set the distributions whose parameters
  are fixed according from the priors specified at ``dep`` 


.. code-block:: python

   #--steps
   conf['steps'] = {}
   
   #--idis and no hera
   conf['steps'][1]={}
   conf['steps'][1]['dep']=[] 
   conf['steps'][1]['active distributions']=['pdf']
   conf['steps'][1]['datasets']={}
   conf['steps'][1]['datasets']['idis']=[]
   conf['steps'][1]['datasets']['idis'].append(10010) # proton   | F2            | SLAC
   conf['steps'][1]['datasets']['idis'].append(10011) # deuteron | F2            | SLAC
   conf['steps'][1]['datasets']['idis'].append(10016) # proton   | F2            | BCDMS
   conf['steps'][1]['datasets']['idis'].append(10017) # deuteron | F2            | BCDMS
   conf['steps'][1]['datasets']['idis'].append(10020) # proton   | F2            | NMC
   conf['steps'][1]['datasets']['idis'].append(10021) # d/p      | F2d/F2p       | NMC
   
   #--idis + hera
   conf['steps'][2]={}
   conf['steps'][2]['dep']=[1]
   conf['steps'][2]['active distributions']=['pdf']
   conf['steps'][2]['datasets']={}
   conf['steps'][2]['datasets']['idis']=[]
   conf['steps'][2]['datasets']['idis'].append(10010) # proton   | F2            | SLAC
   conf['steps'][2]['datasets']['idis'].append(10011) # deuteron | F2            | SLAC
   conf['steps'][2]['datasets']['idis'].append(10016) # proton   | F2            | BCDMS
   conf['steps'][2]['datasets']['idis'].append(10017) # deuteron | F2            | BCDMS
   conf['steps'][2]['datasets']['idis'].append(10020) # proton   | F2            | NMC
   conf['steps'][2]['datasets']['idis'].append(10021) # d/p      | F2d/F2p       | NMC
   conf['steps'][2]['datasets']['idis'].append(10026) # proton   | sigma red     | HERA II NC e+ (1)
   conf['steps'][2]['datasets']['idis'].append(10027) # proton   | sigma red     | HERA II NC e+ (2)
   conf['steps'][2]['datasets']['idis'].append(10028) # proton   | sigma red     | HERA II NC e+ (3)
   conf['steps'][2]['datasets']['idis'].append(10029) # proton   | sigma red     | HERA II NC e+ (4)
   conf['steps'][2]['datasets']['idis'].append(10030) # proton   | sigma red     | HERA II NC e-
   conf['steps'][2]['datasets']['idis'].append(10031) # proton   | sigma red     | HERA II CC e+
   conf['steps'][2]['datasets']['idis'].append(10032) # proton   | sigma red     | HERA II NC e-
   
   #--idis + dy
   conf['steps'][3]={}
   conf['steps'][3]['dep']=[2]
   conf['steps'][3]['active distributions']=['pdf']
   conf['steps'][3]['datasets']={}
   conf['steps'][3]['datasets']['idis']=[]
   conf['steps'][3]['datasets']['idis'].append(10010) # proton   | F2            | SLAC
   conf['steps'][3]['datasets']['idis'].append(10011) # deuteron | F2            | SLAC
   conf['steps'][3]['datasets']['idis'].append(10016) # proton   | F2            | BCDMS
   conf['steps'][3]['datasets']['idis'].append(10017) # deuteron | F2            | BCDMS
   conf['steps'][3]['datasets']['idis'].append(10020) # proton   | F2            | NMC
   conf['steps'][3]['datasets']['idis'].append(10021) # d/p      | F2d/F2p       | NMC
   conf['steps'][3]['datasets']['idis'].append(10026) # proton   | sigma red     | HERA II NC e+ (1)
   conf['steps'][3]['datasets']['idis'].append(10027) # proton   | sigma red     | HERA II NC e+ (2)
   conf['steps'][3]['datasets']['idis'].append(10028) # proton   | sigma red     | HERA II NC e+ (3)
   conf['steps'][3]['datasets']['idis'].append(10029) # proton   | sigma red     | HERA II NC e+ (4)
   conf['steps'][3]['datasets']['idis'].append(10030) # proton   | sigma red     | HERA II NC e-
   conf['steps'][3]['datasets']['idis'].append(10031) # proton   | sigma red     | HERA II CC e+
   conf['steps'][3]['datasets']['idis'].append(10032) # proton   | sigma red     | HERA II NC e-
   conf['steps'][3]['datasets']['dy']=[]
   conf['steps'][3]['datasets']['dy'].append(10001) 
   conf['steps'][3]['datasets']['dy'].append(10002) 







