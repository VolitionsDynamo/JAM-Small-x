Simple 
======

Objectives 
----------

1. Fit unpolarized PDFs from world inclusive DIS data (*single fit*). 
2. Use the fitted parameters to compute PDFs and DIS observables 
3. Next steps  
 
1. Fit unpolarized PDFs
-----------------------

Create a empty directory **workspace** (or with any other name)
which does not need to be inside the fitpack.
e.g.

.. code-block:: shell

   jam
   jam/fitpack
   jam/workspace/
      
In what follow we will create python scripts under *jam/workspace*.
Specifically we will create two scripts 

- **input.py**: this specifies what to fit (datasets, pdfs, etc)
- **driver.py**: this will controls the JAM machinery


input.py
::::::::

First we add some general QCD flags

.. code-block:: python

   import os
   conf={}

   #--setup posterior sampling 
   conf['bootstrap']=False
   conf['flat par']=False
   conf['ftol']=1e-8

   #--setup qcd evolution
   conf['dglap mode']='truncated'
   conf['alphaSmode']='backward'
   conf['order'] = 'NLO'
   conf['Q20']   = 1.27**2

We now add addtional flags specific to DIS.  
In python we call DIS as IDIS to avoid clash 
with python's internal library that is also called
dis. 

.. code-block:: python

   #--setup for idis
   conf['tmc']   = False
   conf['ht']    = False
   conf['nuc smearing']=True
   conf['sidis nuc smearing']=False
   conf['hq']=False

   #--grids
   conf['path2idistab']   = '%s/grids/grids-idis/distab/'%os.environ['FITPACK']

Notice that the path is relative to an enviroment variable `FITPACK` which 
need to be set in our commadline rc file.   

Next we create a space for the datasets

.. code-block:: python

   #--datasets
   conf['datasets']={}

and fill the entries with world  DIS datasets

.. code-block:: python

   #--lepton-hadron reactions
   ##--IDIS
   conf['datasets']['idis']={}
   conf['datasets']['idis']['xlsx']={}
   conf['datasets']['idis']['xlsx'][10010]='idis/expdata/10010.xlsx' # proton   | F2            | SLAC
   conf['datasets']['idis']['xlsx'][10011]='idis/expdata/10011.xlsx' # deuteron | F2            | SLAC
   conf['datasets']['idis']['xlsx'][10016]='idis/expdata/10016.xlsx' # proton   | F2            | BCDMS
   conf['datasets']['idis']['xlsx'][10017]='idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS
   conf['datasets']['idis']['xlsx'][10020]='idis/expdata/10020.xlsx' # proton   | F2            | NMC
   conf['datasets']['idis']['xlsx'][10021]='idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC

The datasets are located at `fitpack/database` and can be viewed with 
a spreadsheed reader. Next we specify DIS cuts 

.. code-block:: python

   Q2cut=1.3**2
   W2cut=10.0

   conf['datasets']['idis']['filters']=[]
   conf['datasets']['idis']['filters'].append("Q2>%f"%Q2cut)
   conf['datasets']['idis']['filters'].append("W2>%f"%W2cut)

In addition we have to add normalization parameters to be fitted. 

.. code-block:: python

   conf['datasets']['idis']['norm']={}
   conf['datasets']['idis']['norm'][10010]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10011]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10016]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10017]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
   conf['datasets']['idis']['norm'][10020]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}

The `min` and `max` specify the allowed ranges for these parameters while
the actual normalization uncertaninty is specified within each `xlsx` table. 

We proceed now to specify the paramters for the `pdfs`. Firts we need
to specify the parametrization. This is done as 

.. code-block:: python

    conf['pdf parametrization'] = 2

This indicates that the fitter will load the parametrization inside 
`fitpack/qcdlib/pdf2.py`.  You can create your own version e.g. 
fitpack/qcdlib/pdfX.py` and adjust 
`fitpack/fitlib/parman.py`
`fitpack/fitlib/resman.py`
accordingly

We next create an entry called `params` followed by entries for the pdf parameters

.. code-block:: python

    #--parameters
    conf['params'] = {}

    conf['params']['pdf'] = {}
    
    conf['params']['pdf']['g1 N']    ={'value':    1, 'min':  None, 'max':  None, 'fixed': True }
    conf['params']['pdf']['g1 a']    ={'value': -0.5, 'min':  -1.9, 'max':     1, 'fixed': False}
    conf['params']['pdf']['g1 b']    ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}
    
    conf['params']['pdf']['uv1 N']   ={'value':    1, 'min':  None, 'max':  None, 'fixed': True }
    conf['params']['pdf']['uv1 a']   ={'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': False}
    conf['params']['pdf']['uv1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}
    
    conf['params']['pdf']['dv1 N']   ={'value':    1, 'min':  None, 'max':  None, 'fixed': True }
    conf['params']['pdf']['dv1 a']   ={'value': -0.5, 'min':  -0.5, 'max':     1, 'fixed': False}
    conf['params']['pdf']['dv1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}
    
    conf['params']['pdf']['db1 N']   ={'value':    1, 'min':     0, 'max':     1, 'fixed': False}
    conf['params']['pdf']['db1 a']   ={'value': -0.5, 'min':    -1, 'max':     1, 'fixed': False}
    conf['params']['pdf']['db1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}
    
    conf['params']['pdf']['ub1 N']   ={'value':    1, 'min':     0, 'max':     1, 'fixed': False}
    conf['params']['pdf']['ub1 a']   ={'value': -0.5, 'min':    -1, 'max':     1, 'fixed': False}
    conf['params']['pdf']['ub1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}
    
    conf['params']['pdf']['s1 N']    ={'value':    1, 'min':     0, 'max':     1, 'fixed': False}
    conf['params']['pdf']['s1 a']    ={'value':    0, 'min':    -1, 'max':     1, 'fixed': False}
    conf['params']['pdf']['s1 b']    ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}
    
    conf['params']['pdf']['sb1 N']   ={'value':    1, 'min':     0, 'max':     1, 'fixed': 's1 N'}
    conf['params']['pdf']['sb1 a']   ={'value':    0, 'min':    -1, 'max':     1, 'fixed': 's1 a'}
    conf['params']['pdf']['sb1 b']   ={'value':    6, 'min':     0, 'max':    10, 'fixed': 's1 b'}
    
    conf['params']['pdf']['sea1 N']  ={'value':  0.5, 'min':     0, 'max':     1, 'fixed': False}
    conf['params']['pdf']['sea1 a']  ={'value': -1.5, 'min':  -1.9, 'max':    -1, 'fixed': False}
    conf['params']['pdf']['sea1 b']  ={'value':    6, 'min':     0, 'max':    10, 'fixed': False}
    
    conf['params']['pdf']['sea2 N']  ={'value':    1, 'min':     0, 'max':     1, 'fixed': 'sea1 N'}
    conf['params']['pdf']['sea2 a']  ={'value': -1.5, 'min':  -1.9, 'max':    -1, 'fixed': 'sea1 a'}
    conf['params']['pdf']['sea2 b']  ={'value':    6, 'min':     0, 'max':    10, 'fixed': 'sea1 b'}


`min` and `max` determines the allowed ranges for the parameters to vary.
The entry  `fixed` works as follow:

- **False**: the parameter is free to vary 
- **True**: the parameter is fixed to the numerical value in entry  `value`
- **other**: the parameter is fixed to be the same as another parmeter that has that entry


driver.py
:::::::::

First make this a python executable by adding the following line 
at the begining

.. code-block:: python

   #!/usr/bin/env python

Then in the shell  type

.. code-block:: shell

   chmod +x driver.py

Then you should be able execute it as 

.. code-block:: shell

   ./driver.py

Now we start crafting the driver. First add the folloing 
imports

.. code-block:: python

   import sys,os
   import numpy as np
   import pylab as py

   #--from fitlib
   import tools.config
   from  fitlib.simple import MAXLIKE

In this example, we use the fitter code from `fitpack/fitlib/simple.py`. 
The class `MAXLIKE` orchestrates all the additionals setups
that is needed to make the fit with the specs from the `input.py`

We now proceed to create a main to call the run 
the `MAXLIKE` class.

.. code-block:: python

   def main1():

       Input='./input.py' #--path to inputfile
       ncores=2   #--number of cores to be used         
       MAXLIKE(Input,ncores).run()

Comments:

- **Input**: this is the path to the `input.py`. 
  You can name the input file as with other name as needed.

- **ncores**: physical cpu cores where the code will run in parallel.

Add the execution block at the end 

.. code-block:: python

   if __name__=="__main__":

       main1()

Run the script and the code should start the fitting process. 
Upon completion the code generates two files 

- **output.py**: This code is the same as the `input.py`, 
  but with the parameters of fit. You can run the fitting 
  process again using the new input file and see that the chi2 
  are already small.  

- **summary**: The chi2 summary of the fit.


2. Compute PDFs and DIS observables 
-----------------------------------

We proceed to use the results from the `output.py` to compute 
PDF and DIS observables. There are many ways you can do this, 
and it depends on your knowdlege of the JAM ecosystem. 
Here we will show a simple way to do this, and craft 
the code lines in the existing `driver.py`. 

Firts with need to add some additional imports

.. code-block:: python

   #--from tools
   from tools.tools     import checkdir,save,load
   from tools.config    import load_config, conf, options
   from tools.inputmod  import INPUTMOD
   from tools.randomstr import id_generator
   
   #--from fitlib
   from fitlib.resman import RESMAN

The strcture of the code should look like this

.. code-block:: python

   def main2():

       #--load input
       load_config('./output.py')

       #--initialize resman
       nworkers=1
       parallel=False
       datasets=False
       resman=RESMAN(nworkers,parallel,datasets)
       parman=resman.parman

       #..DO SOME CALCS..

       #--close resman
       resman.shutdown()



Here we first load the `output.py` into the dictionary `conf`.
In JAM `conf` is like analogous to a global common block in 
fortran. Unlike fortran where we can only store numbers, 
strings, in `conf` is a dictionary and in JAM we store 
instantons of classes such as `pdf` from which we can 
compute values for the pdfs. 


Next we launch the `RESMAN` class where we pass the following 
parameters: `nworkers`, `parallel` and `datasets`. 
The first two specify the parallelization in currently the codes
only parallelizes the calculation of the chi2 (which involves)
`datasets`. So if `datasets=False`, there is not point in 
using paralleization resources. In that case we just use
`nworkers=1`, `parallel=False`. 
On the other hand if you want to be able to compute chi2, 
theory predictions for the datasets then you can set those
parameters as  
`nworkers=4`, `parallel=True` and `datasets=True`

It is important that you reade the 

- `fitpack/fitlib/resman.py`
- `fitpack/fitlib/parman.py`. 

where Abehind the scenes all the setups occur when loading 
the class `RESMAN` 
 
To complete the task for this section lets first solve DGLAP
and evalute PDFs for some :math:`x,Q^2`. To do that we simple 
modify the `main2` in the `driver.py` as follow

.. code-block:: python

   def main2():

       #--load input
       load_config('./output.py')

       #--initialize resman
       nworkers=1
       parallel=False
       datasets=False
       resman=RESMAN(nworkers,parallel,datasets)
       parman=resman.parman

       x=0.5
       Q2=10.0
       conf['pdf'].evolve(Q2)
       print(conf['pdf'].get_xF(x,Q2,'u'))
       print(conf['pdf'].get_xF(x,Q2,'d'))
       print(conf['pdf'].get_xF(x,Q2,'g'))

       #--close resman
       resman.shutdown()

The entry `pdf` contains an instanton of the class PDF. The class 
is located inside `fitpack/qcdlib/pdfX.py`, where `X` is the specific
parametrization indicated in your `input.py`
The method `evolve` evolves the pdfs from the inputscale to the specified 
value of :math:`Q^2` 

Finally, lets compute DIS observables. This lets create `main3` as follows


.. code-block:: python

   def main3():
   
       #--load input
       load_config('./output.py')
   
       #--initialize resman
       nworkers=2
       parallel=True
       datasets=True
       resman=RESMAN(nworkers,parallel,datasets)
       parman=resman.parman
   
       par=parman.par
       res,rres,nres = resman.get_residuals(par)
   
       print('npts=',res.size)
       print('chi2=',np.sum(res**2))
       
       tabs=resman.idisres.tabs
      
       for _ in tabs:
           print()
           print(_) 
           print('X=',tabs[_]['X'][:5])
           print('Q2=',tabs[_]['Q2'][:5])
           print('value=',tabs[_]['value'][:5])
           print('alpha=',tabs[_]['alpha'][:5])
           print('theory=',tabs[_]['thy'][:5])
   
   
       #--close resman
       resman.shutdown()


3. Next steps 
-------------

Before engauging with more complicated setups, here we provide 
some homework.

- Cross sections, structure function level. 

  + Plot :math:`F_2^p` as a function of :math:`Q^2` for fixed values 
    :math:`x`. Do you see scaling? Hint: try to tap into 
    the method `get_FXN` that is defined at `fitpack/obslib/idis/theory.py`. 
    This can be accessed by an entry in the `conf` after loading `RESMAN`

  + Make data vs theory for one of the datasets e.g. SLACp. See how good 
    is the agreement. 

  + Decrease the :math:`W^2_{\rm cut}` in the input file and see how good 
    or bad is the agreement in predicting kinematic regions where the data
    was not included. Do you understand what is going on physically ? 

- Plots at the pdf level

  + Make plots for pdfs. Try to show plot different flavors at different scales.
    See if the PDF looks like those found at https://inspirehep.net/record/1734309

  + Compute the momentum sum rule and show that is scale independent.

  + Compute the quark sum rules and show that they are also scale independet.   

- Changing the parametrizarion

  + Create a new pdfX.py script using the pdf2.py as a template and try to change
    the parametrization and make new fits. 
    Hints: You will have to addapt the `input.py`, `paraman.py` 
    and `resman.py` accordingly. 




 











