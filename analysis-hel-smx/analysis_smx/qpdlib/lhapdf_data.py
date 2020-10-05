import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import lhapdf

## from scipy stack
from scipy.integrate import quad

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

def generate_xf(wdir, groups, flavors, Q2 = None):
    load_config('%s/input.py' % wdir)
    ## setup kinematics
    xs = np.logspace(-2.1, -0.04, 100)
    if Q2 == None: Q2 = conf['Q20']

    print('\ngenerating pdf-proton from LHAPDF at Q2 = %f' % Q2)

    ## compute xf for all replicas
    xfs = {}
    for group in groups:
        xfs[group] = lhapdf_xfs(group, xs, flavors, Q2)
        if group == 'JAM19PDF_proton_nlo':
            xfs[group]['color'] = 'g'
        elif group == 'NNPDF31_nlo_as_0118':
            xfs[group]['color'] = 'm'
        elif group == 'CJ15nlo':
            xfs[group]['color'] = 'k'
        elif group == 'MMHT2014nlo68cl':
            xfs[group]['color'] = 'Yellow'
        elif group == 'CSKK_nnlo_EIG':
            xfs[group]['color'] = 'chocolate'
        elif group == 'ABMP16_3_nlo':
            xfs[group]['color'] = 'b'
    print
    save({'X': xs, 'Q2': Q2, 'XF': xfs, 'groups': groups}, '%s/analysis/qpdlib/lhapdf-%f.dat' % (os.environ['FITPACK'], Q2))

def lhapdf_xfs(pdf_name, xs, flavors, Q2):
    pdf_information = lhapdf.getPDFSet(pdf_name)
    n_xs = len(xs)
    Q = np.sqrt(Q2)

    if (pdf_information.errorType == 'hessian') or (pdf_information.errorType == 'symmhessian'):
        ## 'lhapdf.getPDFSet("ABMP16_3_nlo").errorType' is 'symmhessian'
        pdf_0 = lhapdf.mkPDF(pdf_name, 0)
        xfs = {}
        for flavor in flavors:
            xfs[flavor] = {'center': np.zeros(n_xs), 'difference': np.zeros(n_xs)}
        for i in range(1, ((pdf_information.size - 1) / 2) + 1):
            pdf_plus = lhapdf.mkPDF(pdf_name, 2 * i)
            pdf_minus = lhapdf.mkPDF(pdf_name, 2 * i - 1)
            for j in range(n_xs):
                for flavor in flavors:
                    if flavor == 'uv':
                        xfs[flavor]['center'][j] = pdf_0.xfxQ(2, xs[j], Q) - pdf_0.xfxQ(-2, xs[j], Q)
                        xfs[flavor]['difference'][j] += (pdf_plus.xfxQ(2, xs[j], Q) - pdf_plus.xfxQ(-2, xs[j], Q) - (pdf_minus.xfxQ(2, xs[j], Q) - pdf_minus.xfxQ(-2, xs[j], Q))) ** 2.0
                    elif flavor == 'dv':
                        xfs[flavor]['center'][j] = pdf_0.xfxQ(1, xs[j], Q) - pdf_0.xfxQ(-1, xs[j], Q)
                        xfs[flavor]['difference'][j] += (pdf_plus.xfxQ(1, xs[j], Q) - pdf_plus.xfxQ(-1, xs[j], Q) - (pdf_minus.xfxQ(1, xs[j], Q) - pdf_minus.xfxQ(-1, xs[j], Q))) ** 2.0
                    elif flavor == 'u':
                        xfs[flavor]['center'][j] = pdf_0.xfxQ(2, xs[j], Q)
                        xfs[flavor]['difference'][j] += (pdf_plus.xfxQ(2, xs[j], Q) - pdf_minus.xfxQ(2, xs[j], Q)) ** 2.0
                    elif flavor == 'd':
                        xfs[flavor]['center'][j] = pdf_0.xfxQ(1, xs[j], Q)
                        xfs[flavor]['difference'][j] += (pdf_plus.xfxQ(1, xs[j], Q) - pdf_minus.xfxQ(1, xs[j], Q)) ** 2.0
                    elif flavor == 'd/u':
                        xfs[flavor]['center'][j] = pdf_0.xfxQ(1, xs[j], Q) / pdf_0.xfxQ(2, xs[j], Q)
                        xfs[flavor]['difference'][j] += ((pdf_plus.xfxQ(1, xs[j], Q) - pdf_minus.xfxQ(1, xs[j], Q)) / pdf_0.xfxQ(2, xs[j], Q)) ** 2.0
                    elif flavor == 'db+ub':
                        xfs[flavor]['center'][j] = pdf_0.xfxQ(-1, xs[j], Q) + pdf_0.xfxQ(-2, xs[j], Q)
                        xfs[flavor]['difference'][j] += (pdf_plus.xfxQ(-1, xs[j], Q) + pdf_plus.xfxQ(-2, xs[j], Q) - (pdf_minus.xfxQ(-1, xs[j], Q) + pdf_minus.xfxQ(-2, xs[j], Q))) ** 2.0
                    elif flavor == 'db-ub':
                        xfs[flavor]['center'][j] = pdf_0.xfxQ(-1, xs[j], Q) - pdf_0.xfxQ(-2, xs[j], Q)
                        xfs[flavor]['difference'][j] += (pdf_plus.xfxQ(-1, xs[j], Q) - pdf_plus.xfxQ(-2, xs[j], Q) - (pdf_minus.xfxQ(-1, xs[j], Q) - pdf_minus.xfxQ(-2, xs[j], Q))) ** 2.0
                    elif flavor == 's+sb':
                        xfs[flavor]['center'][j] = pdf_0.xfxQ(3, xs[j], Q) + pdf_0.xfxQ(-3, xs[j], Q)
                        xfs[flavor]['difference'][j] += (pdf_plus.xfxQ(3, xs[j], Q) + pdf_plus.xfxQ(-3, xs[j], Q) - (pdf_minus.xfxQ(3, xs[j], Q) + pdf_minus.xfxQ(-3, xs[j], Q))) ** 2.0
                    elif flavor == 'g':
                        xfs[flavor]['center'][j] = pdf_0.xfxQ(21, xs[j], Q)
                        xfs[flavor]['difference'][j] += (pdf_plus.xfxQ(21, xs[j], Q) - pdf_minus.xfxQ(21, xs[j], Q)) ** 2.0
                    elif flavor == 'rs':
                        if pdf_name == 'CJ15nlo':
                            ## R_s is set to be 0.4 in CJ15
                            xfs[flavor]['center'][j] = 0.4
                            xfs[flavor]['difference'][j] += 0.0
                        else:
                            xfs[flavor]['center'][j] = (pdf_0.xfxQ(3, xs[j], Q) + pdf_0.xfxQ(-3, xs[j], Q)) / (pdf_0.xfxQ(-1, xs[j], Q) + pdf_0.xfxQ(-2, xs[j], Q))
                            xfs[flavor]['difference'][j] += ((pdf_plus.xfxQ(3, xs[j], Q) + pdf_plus.xfxQ(-3, xs[j], Q) - (pdf_minus.xfxQ(3, xs[j], Q) + pdf_minus.xfxQ(-3, xs[j], Q))) / \
                                (pdf_0.xfxQ(-1, xs[j], Q) + pdf_0.xfxQ(-2, xs[j], Q))) ** 2.0

        for flavor in flavors:
            xfs[flavor]['difference'] = 0.5 * np.sqrt(xfs[flavor]['difference'])

    elif pdf_information.errorType == 'replicas':
        pdf_all = lhapdf.mkPDFs(pdf_name)

        xfs = {}
        xfs_temp = {}
        for flavor in flavors:
            xfs[flavor] = {'center': [], 'difference': []}
        for i in range(n_xs):
            for flavor in flavors:
                xfs_temp[flavor] = []
            for pdf_individual in pdf_all:
                for flavor in flavors:
                    if flavor == 'uv':
                        xfs_temp[flavor].append(pdf_individual.xfxQ(2, xs[i], Q) - pdf_individual.xfxQ(-2, xs[i], Q))
                    elif flavor == 'dv':
                        xfs_temp[flavor].append(pdf_individual.xfxQ(1, xs[i], Q) - pdf_individual.xfxQ(-1, xs[i], Q))
                    elif flavor == 'u':
                        xfs_temp[flavor].append(pdf_individual.xfxQ(2, xs[i], Q))
                    elif flavor == 'd':
                        xfs_temp[flavor].append(pdf_individual.xfxQ(1, xs[i], Q))
                    elif flavor == 'd/u':
                        xfs_temp[flavor].append(pdf_individual.xfxQ(1, xs[i], Q) / pdf_individual.xfxQ(2, xs[i], Q))
                    elif flavor == 'db+ub':
                        xfs_temp[flavor].append(pdf_individual.xfxQ(-1, xs[i], Q) + pdf_individual.xfxQ(-2, xs[i], Q))
                    elif flavor == 'db-ub':
                        xfs_temp[flavor].append(pdf_individual.xfxQ(-1, xs[i], Q) - pdf_individual.xfxQ(-2, xs[i], Q))
                    elif flavor == 's+sb':
                        xfs_temp[flavor].append(pdf_individual.xfxQ(3, xs[i], Q) + pdf_individual.xfxQ(-3, xs[i], Q))
                    elif flavor == 'g':
                        xfs_temp[flavor].append(pdf_individual.xfxQ(21, xs[i], Q))
                    elif flavor == 'rs':
                        xfs_temp[flavor].append((pdf_individual.xfxQ(3, xs[i], Q) + pdf_individual.xfxQ(-3, xs[i], Q)) / (pdf_individual.xfxQ(-1, xs[i], Q) + pdf_individual.xfxQ(-2, xs[i], Q)))

            for flavor in flavors:
                xfs[flavor]['center'].append(np.mean(xfs_temp[flavor], axis = 0))
                xfs[flavor]['difference'].append(np.std(xfs_temp[flavor], axis = 0))

        for flavor in flavors:
            xfs[flavor]['center'] = np.array(xfs[flavor]['center'])
            xfs[flavor]['difference'] = np.array(xfs[flavor]['difference'])

    return xfs

if __name__ == '__main__':
    pass
