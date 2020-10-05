Inclusive DIS 
=============

in progress


Theory
------

TMC
...


.. math::

   \frac{ \sum_{t=0}^{N}f(t,k) }{N}


:math:`\frac{ \sum_{t=0}^{N}f(t,k) }{N}`




The AOT is implement... at ...

.. code-block:: python

    def get_FXN(self,x,Q2,stf='F2',twist=2,nucleon='proton')
  



Nuclear Corrections
...................

 
1. Fit unpolarized PDFs
-----------------------

Create a empty directory **workspace** (or with any other name)
which does not need to be inside the fitpack.
e.g.



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

