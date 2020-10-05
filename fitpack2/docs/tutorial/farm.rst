using the farm 
=======================

Stept 1: prepare the input.py 
-----------------------------

We introduce a new element called steps in the `input.py`. 
It looks like this

.. code-block:: python

   #--steps
   conf['steps']={}
   
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

There are some advantages of using this which we will explain in detail.
But for testing purposes, use the `input.py` template from 
`fitpack/workspace/`


Stept 2: test locally if the program runs under the input.py   
------------------------------------------------------------

Open the `test` script and look for the lines 

.. code-block:: python

   #--ATTENTION: user dependent paths
   user    ='username'  
   fitpack ='path-to-fitpack'    
   wdir    ='path-to-working-directory'
   python  ='path-to-your-python-binary'

Setup the entries according to your username. i.e.
   
.. code-block:: python

   user    ='nsato'
   path    ='/w/general-scifs17exp/JAM/'
   fitpack ='%s/nsato/jam/fitpack'%path
   wdir    ='%s/nsato/jam/pions-qt'%path
   python  ='%s/apps/anaconda2/bin/python'%path
 
Run the script 

.. code-block:: shell

   ./test

If it works, the program should start the fitting process. 

Stept 3: send jobs to farm 
--------------------------

Here you have two options depending on your username permissions


option 1: send jobs to farm via farm
....................................

- As in the test case,  adujust the paths and usernames. 

- example:  ...


option 2: send jobs to farm via sfarm (special permissions are needed) 
......................................................................


1) Check that the setup works  
    
    `./sfarm  0  -d results/stepXX   -p results/stepYY`

    where

   - `results/stepXX`: new directory where the current run will store data

   - `results/stepYY`: priors from which the current run will get guess parameters
     If the code seems to run, then it works 

2) Send all the jobs to the farm  

    ./sfarm  1  -d results/stepXX   -p results/stepYY
    
3) To monitor the status you can just check the number of msr files you have 
   collected

   `watch -n 3 "ls results/stepXX/msr | wc -l"`
              
   you can use the following:

   `squeue -u <username>`

   To cancel all the jobs

   `scancel -u <username>`









