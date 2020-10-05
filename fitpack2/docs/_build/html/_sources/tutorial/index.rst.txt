Tutorial
========


Short summary
-------------
The repositories inside fitpack are organized as follow

- **qcdlib**: collection of script to solve QCD RGE equations such as running of strong coupling, DGLAP and TMD evoltion. 

- **database**: collection of excel spreed sheets for experimental data from differentfacilities and different reactions.

- **obslib**: collection of subrepos containing dedicated scripts to compute observables. Typically these repos contain:
  
  - **theory.py**: a script to compute theoretical observables
  - **reader.py**: a script to load the excel files containing expermental observables from the database
  - **residuals.py**: a script to compute the residuals  

- **fitlib**: collection of scripts to carry out chi-square miminization as well as MC sampling 

- **tools**: dedicated scripts including base codes for reader, residuals as well as other utilities suchs as parallelization, code loading. The script ``config`` place an important role (see below).

- **analysis**: dedicated scripts to analyze the fits including a collection of plotting routines as well as simulation routines

- **nuclib**: a repository with nuclear smearing tables needed for nuclear reactions

- **grids**: a collection of subrepos for tabulated mellin grids


.. toctree::
   :maxdepth: 1
   :caption: Examples:

   simple
   analysis

.. Indices and tables
   ==================
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`





