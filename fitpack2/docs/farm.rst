JLab computing for JAM
======================

- Get a JLab account

- Contact N. Sato to give you access to the unix JAM group 
  as well as to additional computing resources at the lab.

Getting started
:::::::::::::::

From your local machine ssh into Jlab systems 
 
.. code-block:: shell

    ssh -4 -XY  <username>@login.jlab.org
    ssh -XY     jlabl1
    ssh -XY     thypc21

After this you have landed inside ``thypc21``. As usual 
you need to customize your shell enviroiment. 

.. code-block:: shell

    cd  ~
    cp ../.cshrc  ./ 
    cp ../.tmux-completion.tcsh  ./
    cp ../.tmux.conf ./
    cp ../.vimrc ./
    cp -r ../.vim  ./

Fireup a tmux session 

.. code-block:: shell

    tmux new-session -s test 

Next, access to JAM's working directory

.. code-block:: shell

    cd /work/JAM

Create a folder under your name if is not present and 
use this space for all of your work. 

**ATTENTION**: DO NOT execute any heavy computing on ``thypc21``. Use this 
machine to use tmux and other simple things as vim.

Running heavy programs
::::::::::::::::::::::

As mentioned before, use ``thypc21`` as your portal to 
JLab computing systems. You should have access to 

1. **ifarm**: ``ssh ifarm``

2. **vertex**: ``ssh vertex``

The optimal workflow is to use ssh into any of these systems 
from a tmux window launched from ``thypc21``

When sshing into ``ifarm``, you land on arbitrary machine. However
all the machines will use the same shell setups. In contrast 
when sshing into ``vertex``, you land on a different home directory 
hence you need to makesure to adjust the shell setups accordingly. 


Getting a node from the farm
::::::::::::::::::::::::::::

The ``ifarm`` is shared by many users who might be running lots 
of programs. Also the cpu usage is limited. For instance the system 
will kill your jam script python execution if running it with 20 slaves.

To access to more cpu resorces you need to use salloc  

.. code-block:: shell

    salloc --partition theory -n 1
    srun --pty csh

For more option, consult the page https://slurm.schedmd.com/salloc.html.

In addition you can access to gpu resources via

.. code-block:: shell

    salloc --partition gpu  --gres gpu:1
    srun --pty csh


Again, make sure your ``tmux`` session/windows/panels are all launched from 
``thypc21``.  







