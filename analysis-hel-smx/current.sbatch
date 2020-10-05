#!/bin/csh
#SBATCH --account=theory
#SBATCH --nodes 1
#SBATCH --partition=production
#SBATCH --cpus-per-task 5
#SBATCH --mem=2G
#SBATCH --time=0-96:00:00
#SBATCH --constraint=general
#SBATCH --job-name="jam-01-8TXUD3DO6O"
#SBATCH --output=/w/jam-sciwork18/dadamiak/analysis-hel-smx/results/test1/out/jam-01-8TXUD3DO6O.out
#SBATCH --error=/w/jam-sciwork18/dadamiak/analysis-hel-smx/results/test1/out/jam-01-8TXUD3DO6O.err
setenv FITPACK /w/jam-sciwork18/dadamiak/fitpack2
setenv PYTHONPATH /w/jam-sciwork18/dadamiak/fitpack2:/work/JAM/apps/lhapdf2/lib/python2.7/site-packages/
/w/jam-sciwork18/apps/anaconda2/bin/python  /w/jam-sciwork18/dadamiak/fitpack2/fitlib_smx/maxlike.py -n 2 -v 5 input.py -msrh "mv <<fname>> /w/jam-sciwork18/dadamiak/analysis-hel-smx/results/test1/msr/" 
