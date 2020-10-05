### predict.py

The function `get_predictions` generates predictions on observables with the `.msr` files.

It also saves the names of the replicas used to generate predictions last time, it will check the names of the replicas each time it runs, for the following cases

> 1. if the replica names are the same, it will ask you if you would like to regenerate or exit generation of predictions
> 2. if the replica names contains the replica names from last time, it will ask you if you would like to regenerate with all replicas or with only the new replicas and append to previous predictions
> 3. none of above, you will not be asked anything and the program will regenerate with all replicas

Once the program is generating predictions with only new replicas, it also checks the following

> 1. does the order of the parameters from one of the previous and new replicas match
> 2. does the observables from one of the previous and new replicas match
> 3. does the datasets of each observable from one of the previous and new replicas match

if the answer is no for any of the above question, you will be asked if you would like to regenerate with all replicas or exit the program and check.

You can also predefine flags to manage the above situations with input `pre_flags`.

An example would be `{'regenerate_same': False, 'regenerate_part': False, 'regenerate_mismatch': True}`, which would do the following
> 1. setting `regenerate_same` to `False` would select to exit generation of predictions when replica names are the same, i.e. not to regenerate
> 2. setting `regenerate_part` to `False` would select to generate only with new replicas when replica names contains that of last time, i.e. not to regenerate with all replicas
> 3. setting `regenerate_mismatch` to `False` would select to regenerate with all replicas when there is a mismatch
