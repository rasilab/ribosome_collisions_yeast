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
    # canonical translation params
    'k_init': list(2.0**np.arange(-8, 1, 1)),
    # no stalling
    'k_elong_stall_1': [10],

    # set all quality control params to be zero to simulate canonical model
    'k_preterm_no_hit_intact': [0],
    'k_preterm_5_hit_intact': [0],
    'k_preterm_3_hit_intact': [0],
    'k_preterm_both_hit_intact': [0],
    'k_preterm_no_hit_endocleaved': [0],
    'k_preterm_5_hit_endocleaved': [0],
    'k_preterm_3_hit_endocleaved': [0],
    'k_preterm_both_hit_endocleaved': [0],
    'k_cleave_no_hit': [0],
    'k_cleave_5_hit': [0],
    'k_cleave_3_hit': [0],
    'k_cleave_both_hit': [0],
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
input_params = input_params.reset_index(drop=True)

input_params.to_csv('sim.params.tsv', sep='\t')  # write to tab-delimited file
input_params.info()  # display the table of input parameters
