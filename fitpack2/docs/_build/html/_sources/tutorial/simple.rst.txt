Simple 
======

- Create a works space folder  anywhere in your system orinside fitpack2.
  We will call such folder as ``workspace`` in this tutorial.

- We need to create two files in the workspace

  **input.py**: this file will define the scope of what is going to be fitted
  **driver.py**: this file wil serve as our main script  


input.py
::::::::

First we add some general QCD flags

.. code-block:: python

   #--fitting setups
   conf['bootstrap'] = False
   conf['flat par']  = False
   conf['ftol']      = 1e-8
   
   #--setups for DGLAP
   conf['alphaSmode'] = 'backward'
   #conf['alphaSmode'] = 'forward'
   conf['order']      = 'NLO'
   conf['Q20']        = 1.27**2
   

We now add addtional flags specific to DIS.  
In python we call DIS as IDIS to avoid clash 
with python's internal library that is also called
dis. 

.. code-block:: python

   conf['tmc']      = False
   conf['ht']       = False
   conf['nuc']      = False
   conf['offshell'] = False
   conf['hq']       = False

Next we create a space for the datasets

.. code-block:: python

   #--datasets
   conf['datasets']={}

and fill the entries with world  DIS datasets


.. code-block:: python
   
   
   #--IDIS
   conf['datasets'] = {}
   conf['datasets']['idis'] = {}
   conf['datasets']['idis']['xlsx'] = {}
   conf['datasets']['idis']['xlsx'][10010] = 'idis/expdata/10010.xlsx' # proton   | F2            | SLAC
   conf['datasets']['idis']['xlsx'][10016] = 'idis/expdata/10016.xlsx' # proton   | F2            | BCDMS
   conf['datasets']['idis']['xlsx'][10020] = 'idis/expdata/10020.xlsx' # proton   | F2            | NMC
   conf['datasets']['idis']['xlsx'][10011] = 'idis/expdata/10011.xlsx' # deuteron | F2            | SLAC
   conf['datasets']['idis']['xlsx'][10017] = 'idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS
   conf['datasets']['idis']['xlsx'][10021] = 'idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC

The datasets are located at `fitpack/database` and can be viewed with 
a spreadsheed reader. Next we specify DIS cuts 

.. code-block:: python

   Q2cut=1.3**2
   W2cut=10.0

   conf['datasets']['idis']['filters']=[]
   conf['datasets']['idis']['filters'].append("Q2>%f"%Q2cut)
   conf['datasets']['idis']['filters'].append("W2>%f"%W2cut)

In addition we add normalization parameters to be fitted. 


.. code-block:: python

   conf['datasets']['idis']['norm'] = {}
   conf['datasets']['idis']['norm'][10010] ={'value':    1.03591e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10016] ={'value':    9.88788e-01, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10020] ={'value':    1.02603e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10011] ={'value':    1.03121e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10017] ={'value':    1.01588e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   
The `min` and `max` specify the allowed ranges for these parameters while
the actual normalization uncertaninty is specified within each `xlsx` table. 
We proceed next to specify the paramters for the `pdfs`. 

   
.. code-block:: python

   #--parameters
   conf['params'] = {}
   
   #--pdf parameters
   conf['params']['pdf'] = {}
   
   conf['params']['pdf']['g1 N']    ={'value':    3.09994e-01, 'min':  None, 'max':  None, 'fixed': True }
   conf['params']['pdf']['g1 a']    ={'value':   -5.20900e-01, 'min':  -1.9, 'max':     1, 'fixed': False}
   conf['params']['pdf']['g1 b']    ={'value':    4.29360e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['uv1 N']   ={'value':    3.25322e-01, 'min':  None, 'max':  None, 'fixed': True }
   conf['params']['pdf']['uv1 a']   ={'value':   -2.14402e-01, 'min':  -0.6, 'max':     1, 'fixed': False}
   conf['params']['pdf']['uv1 b']   ={'value':    3.04406e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['dv1 N']   ={'value':    1.06672e-01, 'min':  None, 'max':  None, 'fixed': True }
   conf['params']['pdf']['dv1 a']   ={'value':   -3.45404e-01, 'min':  -0.6, 'max':     1, 'fixed': False}
   conf['params']['pdf']['dv1 b']   ={'value':    4.48193e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['db1 N']   ={'value':    3.65346e-02, 'min':     0, 'max':     1, 'fixed': False}
   conf['params']['pdf']['db1 a']   ={'value':   -9.35028e-01, 'min':    -1, 'max':     1, 'fixed': False}
   conf['params']['pdf']['db1 b']   ={'value':    4.48545e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['ub1 N']   ={'value':    1.70043e-02, 'min':     0, 'max':     1, 'fixed': False}
   conf['params']['pdf']['ub1 a']   ={'value':   -1.00000e+00, 'min':    -1, 'max':     1, 'fixed': False}
   conf['params']['pdf']['ub1 b']   ={'value':    1.00000e+01, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['s1 N']    ={'value':    9.91077e-02, 'min':     0, 'max':     1, 'fixed': True}
   conf['params']['pdf']['s1 a']    ={'value':    1.00000e+00, 'min':  -0.6, 'max':     1, 'fixed': False}
   conf['params']['pdf']['s1 b']    ={'value':    4.43290e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['sb1 N']   ={'value':    2.96987e-02, 'min':     0, 'max':     1, 'fixed': False}
   conf['params']['pdf']['sb1 a']   ={'value':   -6.00000e-01, 'min':  -0.6, 'max':     1, 'fixed': False}
   conf['params']['pdf']['sb1 b']   ={'value':    3.56087e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['sea1 N']  ={'value':    3.68792e-03, 'min':     0, 'max':     1, 'fixed': False}
   conf['params']['pdf']['sea1 a']  ={'value':   -1.87906e+00, 'min':  -1.9, 'max':    -1, 'fixed': False}
   conf['params']['pdf']['sea1 b']  ={'value':    8.07746e+00, 'min':     0, 'max':    10, 'fixed': False}
   
   conf['params']['pdf']['sea2 N']  ={'value':    3.68792e-03, 'min':     0, 'max':     1, 'fixed': 'sea1 N'}
   conf['params']['pdf']['sea2 a']  ={'value':   -1.87906e+00, 'min':  -1.9, 'max':    -1, 'fixed': 'sea1 a'}
   conf['params']['pdf']['sea2 b']  ={'value':    8.07746e+00, 'min':     0, 'max':    10, 'fixed': 'sea1 b'}
   

The entries for ``fixed`` are set as follows:

- **False**: the parameter is free to vary 
- **True**: the parameter is fixed to the numerical value in entry  `value`
- **other**: the parameter is fixed to be the same as another parmeter that has that entry name



driver.py
:::::::::

First make this a python executable by adding the following line 
at the begining

.. code-block:: python

   #!/usr/bin/env python
   import sys,os
   import numpy as np
   import pylab as py
   import numpy as np
   
   import matplotlib
   matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
   matplotlib.rc('text',usetex=True)
   import pylab  as py
   
   #--local
   from tools.config   import load_config,conf
   from fitlib.resman  import RESMAN
   from fitlib         import simple

Make the ``driver.py`` an excecutable 

.. code-block:: shell

   chmod +x driver.py

Then you should be able execute it as 

.. code-block:: shell

   ./driver.py

If the imports does not work, it measn you have not soruces the ``setup.sh`` script.
Lets proceed to add some main routines

.. code-block:: python

   def main00():

       np.random.seed(seed=100)
       load_config('input.py')
       nworkers=2
       resman=RESMAN(nworkers,parallel=True)
       resman.test(10)
       resman.shutdown()

   if __name__== "__main__":
       
       main00()

This code will load the ``input.py`` and run the ``test`` method inside ``fitlib/resman.py`` 
``RESMAN`` is the master class controlling everything. The test evalute 10 times the residuals.
This gives you an idea for the speed at which the chi2 is evalauted.

Lets proceed to carryout some fits. For this add the following main in the code

.. code-block:: python

   def main01():
       simple.MAXLIKE('input.py').run()

and replace the main part with 

.. code-block:: python

   if __name__== "__main__":
       
       main01()

The code will run showing some metrics for the chi2s, the values of the parameters etc. 
Once completed it will generate an ``output.py`` which is essentially the same as the ``input.py``
but with the parameters in the dictionaries updated. 
For MC procedures, we will run the code slighly differently. 

The simple example shows is very instructive when testing new parametrization, addition of new datasets or new observables. 
Using the ``output.py`` it is simple to craft scripts to produce plots for the PDFs, observables etc.

For instance to examine and/or make plots for the observables that have been fitted we can craft the following

.. code-block:: python

   def main02():
       load_config('output.py')
       nworkers=2
       resman=RESMAN(nworkers,parallel=True)
       par=resman.parman.par
       res,rres,nres=resman.get_residuals(par)
       
       for idx in resman.idis_tabs:
           tab=resman.idis_tabs[idx]
           tab=pd.DataFrame(tab)
           print(tab)
       resman.shutdown()

Notice that we load the ``output.py`` instead of the ``input.py`` as we want to use the most update parameters from a fit.
The ``resman`` object will have the member ``idis_tabs`` provided IDIS datasets are present in the data entries of the 
input file. This member is a dictionary containing each data set as a dictionary labeled by the ``idx`` variable. 
Typically the ``tab`` will contain the following columns 

- **obs**: the lable for the experimental observables e.g. :math:`F_2, \sigma_{\rm red}` etc 
- **value**: the quoted experimetal value of the observable
- **alpha**: uncorrelated uncertainties added in quadrature
- **predictions**: the theory prediction using the fitted pdfs

For each observable there are additional columns specifiying the kinematics of the observable. 
For insrance, DIS tables will have the columns ``X`` and ``Q2``. 


In order to evaluate the PDFs that we have fitted we proceed as follow 

.. code-block:: python
   
   def main03():
       load_config('output.py')
       resman=RESMAN(datasets=False)
       pdf=conf['pdf']
       xi=0.5
       mu2=10.0
       flav='i'
       print(pdf.get_xF(xi,mu2,flav))

The PDF class (see ``qcdlib/pdf.py``) is loaded from ``RESMAN`` into the global dictionary ``conf``.
The parameters of for the PDFs are internally updated by the values in the input file. 

The dictionarry ``conf`` is defined at ``tools/config.py``. This conf dictioary is used as a global 
object similar to common blocks in Fortran.








