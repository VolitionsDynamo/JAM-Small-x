Analysis (OUTDATED)
========


At the workspace:

0) Add a step entry at the kemeanconf.py  (choose number of cluster = 1)

1) In order to analyze the data e.g. compute quatities,  make plots etc, 
   use the driver.py (modify as it is needed). Run the script as
   
   `./driver.py results/stepXX`
     
   the results should all be stored inside `results/stepXX`

2) The driver calls modules located at fitpack/analysis
   
   `analysis/corelib`:  core scripts and it is unlikely that you will touch  these

   `analysis/obslib`: here you will place plots associate with data and theory, 
                      kinematics, etc.

   `analysis/qpdlib`: here we make the plots associated with pdfs, moments, ffs,  etc.


3) You should create new scripts under your forked `fitpack/analysis/`. 
   Once your script is ready to shared among other, you can do a pull request.










