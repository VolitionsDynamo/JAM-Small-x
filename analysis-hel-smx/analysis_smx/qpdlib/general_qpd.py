from tools.tools import load
from tools.config import conf, load_config
from analysis.corelib import core
from fitlib.resman import RESMAN

def print_grouped_parameters(working_directory, distribution_name, i_replica = 1):
    ## print parameters from replicas and PDF class to compare
    load_config('%s/input.py' % working_directory)
    istep = core.get_istep()

    replicas = core.get_replicas(working_directory)
    core.mod_conf(istep, replicas[0]) ## set conf as specified in istep

    if distribution_name not in conf['steps'][istep]['active distributions']:
        print('%s-proton is not an active distribution' % distribution_name)
        return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep] ## make sure 'parman' uses the same order as all the replicas do
    # print parman.order

    distribution = conf[distribution_name]

    print('%s parameters in class' % distribution_name)
    core.mod_conf(istep, replicas[i_replica - 1])
    parman.set_new_params(replicas[i_replica - 1]['params'][istep], initial = True)
    for name, value in distribution.params.iteritems():
        if value[0] != 0.0:
            print '%7s: %.5e, %.5e, %.5e, %.5e, %.5e' % (name, value[0], value[1], value[2], value[3], value[4])

    print('%s parameters in replicas' % distribution_name)
    orders = []
    unique_orders = []
    values = []
    for i in range(len(replicas[i_replica - 1]['order'][istep])):
        order = replicas[i_replica - 1]['order'][istep][i]
        value = replicas[i_replica - 1]['params'][istep][i]
        if (order[0] == 1) and (order[1] == distribution_name):
            orders.append(order[2])
            unique_orders.append(order[2].split(' ')[0])
            values.append(value)

    unique_orders = list(set(unique_orders))
    for unique_order in unique_orders:
        parameters = []
        all_parameters = {'N': None, 'a': 0.0, 'b': 0.0, 'c': 0.0, 'd': 0.0}
        for i in range(len(orders)):
            if orders[i].split(' ')[0] == unique_order:
                all_parameters[orders[i].split(' ')[1]] = values[i]
        if all_parameters == {'N': None, 'a': 0.0, 'b': 0.0, 'c': 0.0, 'd': 0.0}:
            continue
        if all_parameters['N'] == None:
            print '%7s: None, %.5e, %.5e, %.5e, %.5e' % (unique_order, all_parameters['a'], all_parameters['b'], all_parameters['c'], all_parameters['d'])
        else:
            print '%7s: %.5e, %.5e, %.5e, %.5e, %.5e' % (unique_order, all_parameters['N'], all_parameters['a'], all_parameters['b'], all_parameters['c'], all_parameters['d'])
        # print order[2].ljust(13), ':', value

def print_individual_parameters(working_directory, distribution, norm, index = 0):
    load_config('%s/input.py' % working_directory)
    istep = core.get_istep()
    jar = load('%s/data/jar-%d.dat' % (working_directory, istep))

    if norm:
        for i in range(len(jar['order'])):
            if jar['order'][i][0] == 2:
                print jar['order'][i][1].ljust(8), '%s' % str(jar['order'][i][2]).ljust(10), ':', jar['replicas'][index][i]
    else:
        for i in range(len(jar['order'])):
            if (jar['order'][i][0] == 1) and (jar['order'][i][1] == distribution):
                print jar['order'][i][2].ljust(13), ':', jar['replicas'][index][i]

if __name__ == '__main__':
    pass
