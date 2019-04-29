#!/usr/bin/env python

"""Choose parameters for simulation
"""


import numpy as np
import pandas as pd
import itertools as it

# include each of these params and one of its values in each simulation,
# they are changed from the default value
k_elong_stall_value = 0.1
and_params = {
    # translation params
    'k_init': list(2.0**np.arange(-8, 1, 1)),
    'k_preterm_no_hit_intact': [0],
    'k_preterm_5_hit_intact': [0],
    'k_preterm_3_hit_intact': [0],
    'k_preterm_both_hit_intact': [0],
    'preterm_intact_model': ['trafficjam'],
    'preterm_intact_rate': [0],
    'n_stall': [6],
    'x_stall': [','.join([str(x) for x in range(400, 406)])],
    'k_elong_stall': [','.join([str(np.round(k_elong_stall * 6,3))] * 6) 
                      for k_elong_stall in [k_elong_stall_value]],
}

# create a list of all parameter combinations from above dict
# there will be as many list elements as the product of all list lengths above
and_params = [list(x) for x in it.product(*[[(p, v)
                                             for v in and_params[p]]
                                            for p in and_params])]

preterm_intact_rate_list = [0, 0.02, 0.1, 0.5, 1, 4] 
cleave_rate = 0.001
preterm_intact_params = list()
for preterm_intact_rate in preterm_intact_rate_list:
    scale_factor = (k_elong_stall_value + preterm_intact_rate) / k_elong_stall_value
    preterm_intact_params.append({
        # this is the CSAT model
        'k_preterm_no_hit_intact': 0,
        'k_preterm_5_hit_intact': preterm_intact_rate,
        'k_preterm_3_hit_intact': 0,
        'k_preterm_both_hit_intact': preterm_intact_rate,
        'preterm_intact_model': 'hit5',
        'preterm_intact_rate': preterm_intact_rate,
        # co-translational cleavage params
        # endocleavage occurs these many codons behind a-site
        # this is the simple cleavage model
        'l_cleave': 5,
        # scale cleavage rate to maintain the same cleavage efficacy despite
        # the abortive termination relieving the stall
        'k_cleave_no_hit': 0,
        'k_cleave_5_hit': cleave_rate * scale_factor,
        'k_cleave_3_hit': 0,
        'k_cleave_both_hit': cleave_rate * scale_factor,
        'cleave_model': 'hit5',
        'cleave_rate': cleave_rate * scale_factor,
    })

# convert each param combination from dict to list of tuples
preterm_intact_params = [list(x.items()) for x in preterm_intact_params]

# combine the 'and' parameters and parameter combinations from above
simcount = 0
temp = dict()
for params in it.product(and_params, preterm_intact_params):
    temp[simcount] = dict(it.chain.from_iterable(params))
    simcount += 1

# convert to pandas dataframe
input_params = pd.DataFrame.from_dict(temp, orient='index')

# sort  the paramters by these parameter combinations
input_params = input_params.sort_values(
    by=['n_stall', 'k_init']).reset_index(drop=True)

input_params.to_csv('sim.params.tsv', sep='\t')  # write to tab-delimited file
input_params.info()  # display the table of input parameters
