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

    # canonical mrna metabolism params
    'k_transcription': [0.000],  # set to 0 for single mrna sim
    'k_deadenylation': [0.00],  # set to 0 for single mrna sim

    # initial conditions
    'n_mrna_0': [1],

    # no cleaving of mRNA
    'k_cleave_no_hit': [0],
    'k_cleave_5_hit': [0],
    'k_cleave_3_hit': [0],
    'k_cleave_both_hit': [0],
    'cleave_model': ['simple'],
    'cleave_rate': [0],
    'n_stall': [6],
    'x_stall': [','.join([str(x) for x in range(400, 406)])],
    'k_elong_stall': [','.join([str(np.round(k_elong_stall * 6,3))] * 6) 
                      for k_elong_stall in [0.02, 0.1, 0.5]],
}

# create a list of all parameter combinations from above dict
# there will be as many list elements as the product of all list lengths above
and_params = [list(x) for x in it.product(*[[(p, v)
                                             for v in and_params[p]]
                                            for p in and_params])]

# rates of premature termination when there is no
# mrna cleavage at the a-site
preterm_intact_rate_list = [0, 0.01, 0.02, 0.05, 0.1] 
preterm_intact_params = list()
for preterm_intact_rate in preterm_intact_rate_list:
    preterm_intact_params.append({
        # this is the SAT model
        'k_preterm_no_hit_intact': preterm_intact_rate,
        'k_preterm_5_hit_intact': preterm_intact_rate,
        'k_preterm_3_hit_intact': 0,
        'k_preterm_both_hit_intact': 0,
        'preterm_intact_model': 'simple',
        'preterm_intact_rate': preterm_intact_rate,
    })
    if preterm_intact_rate == 0:
        preterm_intact_params[-1]['preterm_intact_model'] = 'trafficjam'

preterm_intact_rate_list = [0.02, 0.1, 0.5, 1, 4] 
for preterm_intact_rate in preterm_intact_rate_list:
    preterm_intact_params.append({
        # this is the CSAT model
        'k_preterm_no_hit_intact': 0,
        'k_preterm_5_hit_intact': preterm_intact_rate,
        'k_preterm_3_hit_intact': 0,
        'k_preterm_both_hit_intact': preterm_intact_rate,
        'preterm_intact_model': 'hit5',
        'preterm_intact_rate': preterm_intact_rate,
    })
    preterm_intact_params.append({
        # this is the CAT model
        'k_preterm_no_hit_intact': 0,
        'k_preterm_5_hit_intact': 0,
        'k_preterm_3_hit_intact': preterm_intact_rate,
        'k_preterm_both_hit_intact': preterm_intact_rate,
        'preterm_intact_model': 'hit3',
        'preterm_intact_rate': preterm_intact_rate,
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
    by=['k_init', 'k_elong_stall', 'preterm_intact_rate']).reset_index(
    drop=True)

input_params.to_csv('sim.params.tsv', sep='\t')  # write to tab-delimited file
input_params.info()  # display the table of input parameters
