import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
from scipy.integrate import quad

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
from analysis.corelib import classifier

def gen_xf(wdir, flavors, Q2 = None):
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    ## 'conf' will be modified for each replica individually later in the loop over 'replicas'
    ## the reason for doing this is that 'fix parameters' has to be set correctly for each replica

    if 'ppdf' not in conf['steps'][istep]['active distributions']:
        print('ppdf-proton not in active distribution')
        return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep]
    ## make sure 'parman' uses the same order for active distributions as all the replicas do
    # print parman.order

    ppdf = conf['ppdf']

    ## setup kinematics
    xs = 10.0 ** np.linspace(-3, -1, 100)
    xs = np.append(xs, np.linspace(0.1, 0.99, 100))
    if Q2 == None: Q2 = conf['Q20']
    print('\ngenerating polarized pdf-proton from %s at Q2 = %f' % (wdir, Q2))

    ## compute xf for all replicas
    xfs = {}
    n_replicas = len(replicas)
    for i in range(n_replicas):
        lprint('%d/%d' % (i + 1, n_replicas))

        ## filter
        #flag=False
        #params=replica['params'][istep]
        #order=replica['order'][istep]
        #for i in range(len(order)):
        #    if order[i][0]!=1:continue
        #    if order[i][1]!='ppdf':continue
        #    #if order[i][2]=='s1 a':
        #    #   if params[i]<-0.9: flag=True
        #if flag: continue

        parman.order = copy.copy(replicas[i]['order'][istep])
        core.mod_conf(istep, replicas[i])
        parman.set_new_params(replicas[i]['params'][istep], initial = True)

        for flavor in flavors:
            if flavor not in xfs: xfs[flavor] = []
            if flavor == 'up':
                func = lambda x: ppdf.get_xF(x, Q2, 'u') + ppdf.get_xF(x, Q2, 'ub')
            elif flavor == 'dp':
                func = lambda x: ppdf.get_xF(x, Q2, 'd') + ppdf.get_xF(x, Q2, 'sb')
            elif flavor == 'sp':
                func = lambda x: ppdf.get_xF(x, Q2, 's') + ppdf.get_xF(x, Q2, 'sb')
            else:
                func = lambda x: ppdf.get_xF(x, Q2, flavor)

            xfs[flavor].append([func(x) for x in xs])
    print
    checkdir('%s/data' % wdir)
    if Q2 == conf['Q20']:
        save({'X': xs, 'Q2': Q2, 'XF': xfs}, '%s/data/ppdf-%d.dat' % (wdir, istep))
    else:
        save({'X': xs, 'Q2': Q2, 'XF': xfs}, '%s/data/ppdf-%d-%f.dat' % (wdir, istep, Q2))

def get_parameters(wdir, Q2 = None):
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    ## 'conf' will be modified for each replica individually later in the loop over 'replicas'
    ## the reason for doing this is that 'fix parameters' has to be set correctly for each replica

    if 'ppdf' not in conf['steps'][istep]['active distributions']:
        print('ppdf-proton not in active distribution')
        return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep] ## make sure 'parman' uses the same order as all the replicas do
    # print parman.order

    ppdf = conf['ppdf']

    ## setup kinematics
    if Q2 == None: Q2 = conf['Q20']
    ## get parameters for all flavors of PDF
    print('\ngetting ppdf-parameters from %s at Q2 = %.2f' % (wdir, Q2))

    ## check if any of the shapes are fixed
    shape_1 = []
    shape_2 = []
    shape_3 = []
    for parameter in conf['params']['ppdf']:
        if '1' in parameter:
            shape_1.append(conf['params']['ppdf'][parameter]['fixed'])
        elif '2' in parameter:
            shape_2.append(conf['params']['ppdf'][parameter]['fixed'])
        elif '3' in parameter:
            shape_3.append(conf['params']['ppdf'][parameter]['fixed'])
        else:
            print('there seems to be more than three shapes')
            print(parameter)
            sys.exit('please update script %s' % __file__)

    fixed_shape_1 = False
    fixed_shape_2 = False
    fixed_shape_3 = False
    if (len(shape_1) != 0) and all([_ == True for _ in shape_1]):
        fixed_shape_1 = True
    elif len(shape_1) == 0:
        fixed_shape_1 = True
    if (len(shape_2) != 0) and all([_ == True for _ in shape_2]):
        fixed_shape_2 = True
    elif len(shape_2) == 0:
        fixed_shape_2 = True
    if (len(shape_3) != 0) and all([_ == True for _ in shape_3]):
        fixed_shape_3 = True
    elif len(shape_3) == 0:
        fixed_shape_3 = True

    parameters = {}
    if not fixed_shape_1: parameters[1] = {}
    if not fixed_shape_2: parameters[2] = {}
    if not fixed_shape_3: parameters[3] = {}

    n_replicas = len(replicas)
    for i in range(n_replicas):
        lprint('%d/%d' % (i + 1, n_replicas))

        parman.order = copy.copy(replicas[i]['order'][istep])
        core.mod_conf(istep, replicas[i])
        parman.set_new_params(replicas[i]['params'][istep], initial = True)

        for flavor, value in ppdf.params.iteritems():
            for j in parameters:
                if str(j) in flavor:
                    flavor_key = flavor.replace(str(j), '')
                    if flavor_key not in parameters[j]:
                        parameters[j][flavor_key] = {'n': [], 'a': [], 'b': [], 'c': [], 'd': []}

                    parameters[j][flavor_key]['n'].append(value[0])
                    parameters[j][flavor_key]['a'].append(value[1])
                    parameters[j][flavor_key]['b'].append(value[2])
                    parameters[j][flavor_key]['c'].append(value[3])
                    parameters[j][flavor_key]['d'].append(value[4])
                else:
                    pass

    ## remove c or d parameters if fixed in all flavors and shapes
    ## this is to avoid producing empty figures
    shapes_fixed_c = []
    shapes_fixed_d = []
    for shape in parameters:
        flavors_fixed_c = []
        flavors_fixed_d = []
        for flavor in parameters[shape]:
            cs = parameters[shape][flavor]['c']
            ds = parameters[shape][flavor]['d']
            if all([_ == cs[0] for _ in cs]):
                flavors_fixed_c.append(True)
            else:
                flavors_fixed_c.append(False)
            if all([_ == ds[0] for _ in ds]):
                flavors_fixed_d.append(True)
            else:
                flavors_fixed_d.append(False)

        if all(flavors_fixed_c):
            shapes_fixed_c.append(True)
        else:
            shapes_fixed_c.append(False)
        if all(flavors_fixed_d):
            shapes_fixed_d.append(True)
        else:
            shapes_fixed_d.append(False)

    if all(shapes_fixed_c):
        for shape in parameters:
            for flavor in parameters[shape]:
                del parameters[shape][flavor]['c']
        print('parameter c is fixed for shape %s' % shape)
    if all(shapes_fixed_d):
        for shape in parameters:
            for flavor in parameters[shape]:
                del parameters[shape][flavor]['d']
        print('parameter d is fixed for shape %s' % shape)

    print
    checkdir('%s/data' % wdir)
    if Q2 == conf['Q20']:
        save(parameters, '%s/data/ppdf-parameters-%d.dat' % (wdir, istep))
    else:
        save(parameters, '%s/data/ppdf-parameters-%d-%f.dat' % (wdir, istep, Q2))

def get_moments(wdir, flavors, x_1, x_2, Q2 = None):
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    ## 'conf' will be modified for each replica individually later in the loop over 'replicas'
    ## the reason for doing this is that 'fix parameters' has to be set correctly for each replica

    if 'ppdf' not in conf['steps'][istep]['active distributions']:
        print('ppdf-proton not in active distribution')
        return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep]
    ## make sure 'parman' uses the same order for active distributions as all the replicas do
    # print parman.order

    ppdf = conf['ppdf']

    ## setup kinematics
    # X = 10.0 ** np.linspace(-3, -1, 100)
    # X = np.append(X, np.linspace(0.1, 0.99, 100))
    if Q2 == None: Q2 = conf['Q20']
    print('\ngenerating momentus for polarized pdf-proton from %s at Q2 = %f' % (wdir, Q2))

    ## compute moments for all replicas
    moments = {}
    statistics = {}
    n_replicas = len(replicas)
    for i in range(n_replicas):
        lprint('%d/%d' % (i + 1, n_replicas))

        parman.order = copy.copy(replicas[i]['order'][istep])
        core.mod_conf(istep, replicas[i])
        parman.set_new_params(replicas[i]['params'][istep], initial = True)

        for flavor in flavors:
            if flavor not in moments: moments[flavor] = []
            if flavor not in statistics: statistics[flavor] = {}
            if flavor == 'up':
                func = lambda x: (ppdf.get_xF(x, Q2, 'u') + ppdf.get_xF(x, Q2, 'ub')) / x
            elif flavor == 'dp':
                func = lambda x: (ppdf.get_xF(x, Q2, 'd') + ppdf.get_xF(x, Q2, 'sb')) / x
            elif flavor == 'sp':
                func = lambda x: (ppdf.get_xF(x, Q2, 's') + ppdf.get_xF(x, Q2, 'sb')) / x
            else:
                func = lambda x: ppdf.get_xF(x, Q2, flavor) / x

            moments[flavor].append(quad(func, x_1, x_2)[0])

    for flavor in flavors:
        statistics[flavor]['mean'] = np.mean(moments[flavor])
        statistics[flavor]['standard_deviation'] = np.std(moments[flavor])

    print
    for flavor in flavors:
        print flavor, ': ', statistics[flavor]

    checkdir('%s/data' % wdir)
    if Q2 == conf['Q20']:
        save({'Q2': Q2, 'moments': moments, 'statistics': statistics}, '%s/data/ppdf-moments-%s-to-%s-%d.dat' % (wdir, x_1, x_2, istep))
    else:
        save({'Q2': Q2, 'moments': moments, 'statistics': statistics}, '%s/data/ppdf-moments-%s-to-%s-%d-%f.dat' % (wdir, x_1, x_2, istep, Q2))

if __name__ == '__main__':
    pass
