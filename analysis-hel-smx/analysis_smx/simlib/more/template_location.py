import os
conf={}

##--IDIS
conf['datasets']['idis']={}
conf['datasets']['idis']['filters']=[]
conf['datasets']['idis']['xlsx']={}
conf['datasets']['idis']['xlsx'][10010]= 'template/10010.xlsx' # proton   | F2            | SLAC
conf['datasets']['idis']['xlsx'][10011]= 'template/10011.xlsx' # deuteron | F2            | SLAC
conf['datasets']['idis']['xlsx'][10016]= 'template/10016.xlsx' # proton   | F2            | BCDMS
conf['datasets']['idis']['xlsx'][10017]= 'template/10017.xlsx' # deuteron | F2            | BCDMS
conf['datasets']['idis']['xlsx'][10020]= 'template/10020.xlsx' # proton   | F2            | NMC
conf['datasets']['idis']['xlsx'][10021]= 'template/10021.xlsx' # d/p      | F2d/F2p       | NMC
conf['datasets']['idis']['xlsx'][10026]= 'template/10026.xlsx' # proton   | sigma red     | HERA II NC e+ (1)
conf['datasets']['idis']['xlsx'][10027]= 'template/10027.xlsx' # proton   | sigma red     | HERA II NC e+ (2)
conf['datasets']['idis']['xlsx'][10028]= 'template/10028.xlsx' # proton   | sigma red     | HERA II NC e+ (3)
conf['datasets']['idis']['xlsx'][10029]= 'template/10029.xlsx' # proton   | sigma red     | HERA II NC e+ (4)
conf['datasets']['idis']['xlsx'][10030]= 'template/10030.xlsx' # proton   | sigma red     | HERA II NC e-
conf['datasets']['idis']['xlsx'][10031]= 'template/10031.xlsx' # proton   | sigma red     | HERA II CC e+
conf['datasets']['idis']['xlsx'][10032]= 'template/10032.xlsx' # proton   | sigma red     | HERA II NC e-

##--PIDIS
conf['datasets']['pidis']={}
conf['datasets']['pidis']['filters']=[]
conf['datasets']['pidis']['xlsx']={}
#-------------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['xlsx'][10002]='template/10002.xlsx' # 10002 | proton   | A1   | COMPASS         |          |
conf['datasets']['pidis']['xlsx'][10003]='template/10003.xlsx' # 10003 | proton   | A1   | COMPASS         |          |
conf['datasets']['pidis']['xlsx'][10004]='template/10004.xlsx' # 10004 | proton   | A1   | EMC             |          |
conf['datasets']['pidis']['xlsx'][10007]='template/10007.xlsx' # 10007 | proton   | Apa  | HERMES          |          |
conf['datasets']['pidis']['xlsx'][10008]='template/10008.xlsx' # 10008 | proton   | A2   | HERMES          |          |
conf['datasets']['pidis']['xlsx'][10017]='template/10017.xlsx' # 10017 | proton   | Apa  | JLabHB(EG1DVCS) |          |
conf['datasets']['pidis']['xlsx'][10022]='template/10022.xlsx' # 10022 | proton   | Apa  | SLAC(E143)      |          |
conf['datasets']['pidis']['xlsx'][10023]='template/10023.xlsx' # 10023 | proton   | Ape  | SLAC(E143)      |          |
conf['datasets']['pidis']['xlsx'][10028]='template/10028.xlsx' # 10028 | proton   | Ape  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10029]='template/10029.xlsx' # 10029 | proton   | Apa  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10031]='template/10031.xlsx' # 10031 | proton   | Atpe | SLAC(E155x)     |          |
conf['datasets']['pidis']['xlsx'][10032]='template/10032.xlsx' # 10032 | proton   | Apa  | SLACE80E130     |          |
conf['datasets']['pidis']['xlsx'][10035]='template/10035.xlsx' # 10035 | proton   | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10036]='template/10036.xlsx' # 10036 | proton   | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10041]='template/10041.xlsx' # 10041 | proton   | Apa  | JLabHB(EG1b)    | E =1 GeV |
conf['datasets']['pidis']['xlsx'][10042]='template/10042.xlsx' # 10042 | proton   | Apa  | JLabHB(EG1b)    | E =2 GeV |
conf['datasets']['pidis']['xlsx'][10043]='template/10043.xlsx' # 10043 | proton   | Apa  | JLabHB(EG1b)    | E =4 GeV |
conf['datasets']['pidis']['xlsx'][10044]='template/10044.xlsx' # 10044 | proton   | Apa  | JLabHB(EG1b)    | E =5 GeV |
conf['datasets']['pidis']['xlsx'][10005]='template/10005.xlsx' # 10005 | neutron  | A1   | HERMES          |          |
#-------------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['xlsx'][10001]='template/10001.xlsx' # 10001 | deuteron | A1   | COMPASS         |          |
conf['datasets']['pidis']['xlsx'][10006]='template/10006.xlsx' # 10006 | deuteron | Apa  | HERMES          |          |
conf['datasets']['pidis']['xlsx'][10016]='template/10016.xlsx' # 10016 | deuteron | Apa  | JLabHB(EG1DVCS) |          |
conf['datasets']['pidis']['xlsx'][10020]='template/10020.xlsx' # 10020 | deuteron | Ape  | SLAC(E143)      |          |
conf['datasets']['pidis']['xlsx'][10021]='template/10021.xlsx' # 10021 | deuteron | Apa  | SLAC(E143)      |          |
conf['datasets']['pidis']['xlsx'][10026]='template/10026.xlsx' # 10026 | deuteron | Ape  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10027]='template/10027.xlsx' # 10027 | deuteron | Apa  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10030]='template/10030.xlsx' # 10030 | deuteron | Atpe | SLAC(E155x)     |          |
conf['datasets']['pidis']['xlsx'][10033]='template/10033.xlsx' # 10033 | deuteron | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10034]='template/10034.xlsx' # 10034 | deuteron | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10037]='template/10037.xlsx' # 10037 | deuteron | Apa  | JLabHB(EG1b)    | E =1 GeV |
conf['datasets']['pidis']['xlsx'][10038]='template/10038.xlsx' # 10038 | deuteron | Apa  | JLabHB(EG1b)    | E =2 GeV |
conf['datasets']['pidis']['xlsx'][10039]='template/10039.xlsx' # 10039 | deuteron | Apa  | JLabHB(EG1b)    | E =4 GeV |
conf['datasets']['pidis']['xlsx'][10040]='template/10040.xlsx' # 10040 | deuteron | Apa  | JLabHB(EG1b)    | E =5 GeV |
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['xlsx'][230010]='template/30010.xlsx' # 30005 | neutron  | Apa   | EIC simul
conf['datasets']['pidis']['xlsx'][230001]='template/30001.xlsx' # 30001 | neutron  | Apa   | EIC simul
conf['datasets']['pidis']['xlsx'][230002]='template/30002.xlsx' # 30002 | neutron  | Apa   | EIC simul
conf['datasets']['pidis']['xlsx'][230003]='template/30003.xlsx' # 30003 | neutron  | Apa   | EIC simul
conf['datasets']['pidis']['xlsx'][230004]='template/30004.xlsx' # 30004 | neutron  | Apa   | EIC simul
conf['datasets']['pidis']['xlsx'][230005]='template/30005.xlsx' # 30005 | neutron  | Apa   | EIC simul
conf['datasets']['pidis']['xlsx'][230006]='template/30006.xlsx' # 30006 | neutron  | Apa   | EIC simul
conf['datasets']['pidis']['xlsx'][230007]='template/30007.xlsx' # 30007 | neutron  | Apa   | EIC simul
conf['datasets']['pidis']['xlsx'][230008]='template/30008.xlsx' # 30008 | neutron  | Apa   | EIC simul
conf['datasets']['pidis']['xlsx'][230009]='template/30009.xlsx' # 30009 | neutron  | Apa   | EIC simul
#---------------------------------------------------------------------------------------------------------------------------
