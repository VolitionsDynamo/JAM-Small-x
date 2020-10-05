.. JAM documentation master file, created by
   sphinx-quickstart on Wed Dec  4 16:10:03 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: ./logos/jam.jpg

JAM ECOSYSTEM
=============

Program that provides a framework to study nucleon structure and hadronization
in QCD. It provides the tools to extract QFT objects suchs as parton
distribution functons or fragmentation functions from experimental cross
sections

.. toctree::
   :hidden:

   self

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   start
   tutorial/index
   input
   observables/index
   farm
   simulation/index

.. Indices and tables
   ==================
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

Contribute to the documentation
===============================

The documentation are located at `fitpack/docs`. 
If you want to help improving the documentation
you can do it simple as modifying/adding rst files 
inside the `docs` and the do a pull request. 
You need to have sphinx and the rtd them 
installed in your system 

.. code-block:: shell
   
   pip install -U sphinx
   pip install sphinx-rtd-theme

Once you modify the `docs` type

.. code-block:: shell
   
   cd docs 
   make html
   firefox index.html

You will need to reload the page once you run the `make`
again. Finally add your author ship at the end of rst files 
so that people can blame you :).  

Authors
==================

- Carlota Andres (JLab) 
- Nobuo Sato (JLab)
- Jake Ethier (Nikef)




