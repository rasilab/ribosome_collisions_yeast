#!/usr/bin/env python

"""Choose parameters for simulation
This program will set the parameter combinations that are considered in this
simulation run. See README.org for background on why these parameters were
chosen.
"""

import numpy as np
import pandas as pd
import itertools as it

# include each of these params and one of its values in each simulation
and_params = {
    'k_elong_stall': [0.1],
    'l_mrna': list(range(100,700,100)),
    'x_stall': [80],
}

# create a list of all parameter combinations from above dict
# there will be as many list elements as the product of all list lengths above
and_params = [list(x) for x in it.product(*[[(p, v)
                                             for v in and_params[p]]
                                            for p in and_params])]

# combine the 'and' parameters and parameter combinations from above
simcount = 0
temp = dict()
for params in it.product(and_params):
    temp[simcount] = dict(it.chain.from_iterable(params))
    simcount += 1

# convert to pandas dataframe
input_params = pd.DataFrame.from_dict(temp, orient='index')

# sort  the paramters by these parameter combinations
input_params = input_params.sort_values(
    by=['l_mrna']).reset_index(
    drop=True)

input_params.to_csv('sim.params.tsv', sep='\t')  # write to tab-delimited file
input_params.info()  # display the table of input parameters
