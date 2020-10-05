Gettin Started
==============

dependencies
------------

- Only Linux and OSX is supported
- We recommend to install anaconda (python3) 

installation
------------

1. Create and activate a new python enviroment

.. code-block:: shell

   conda create --name  jam  python=3.6 
   conda activate jam

2. Clone the fitpack repository

.. code-block:: shell

   git clone  git@github.com:QCDHUB/fitpack2
   cd fitpack2

3. Install dependencies on your enviroment

.. code-block:: shell

   pip install -r deps

4. Install core and default libraries

.. code-block:: shell

   ./pacman.py install

Additional repos for observables can be installed by 
adjusting the file ``repos.tex``. Also to updates on all 
repos can be done via

.. code-block:: shell

   ./pacman.py update


setups
------

Some environmental variables  need to be set. For csh 

.. code-block:: shell

   setenv FITPACK $PWD
   setenv PATH ${FITPACK}/bin:${PATH}
   setenv PYTHONPATH ${FITPACK}

For bash

.. code-block:: shell

   path=`pwd`
   export FITPACK=$path
   export PATH=$FITPACK/bin:$PATH
   export PYTHONPATH=$FITPACK:$PYTHONPATH

Alternatively source the setup file ``fitpack2/setup.sh`` or 
``fitpack2/setup.csh``

.. code-block:: shell

   source setup.sh

To test the installation do 

.. code-block:: shell

   cd fitlib
   ./driver.py 


Alternative installation
------------------------

The following repos can only be modified 
via pull request. 

- fitlib
- tools

You can forked them into your account 
and before running ``./pacman.py install``
replace ``QCDHUB`` by your github username
in order to clone from you account.

To update your fork with recent changes at the fitpack upstream 
you need to do the following within fitlib and tools  

.. code-block:: shell

   git remote add upstream git@github.com:JeffersonLab/fitpack.git
   git fetch upstream

This is done only once. After that you can sync the fork using 

.. code-block:: shell

   git pull upstream master


Next steps
----------

Checkout the tutorials














