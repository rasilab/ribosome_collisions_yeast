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
    'k_init': [2**(-3)],
    'k_preterm_no_hit_intact': [0],
    'k_preterm_5_hit_intact': [0],
    'k_preterm_3_hit_intact': [0],
    'k_preterm_both_hit_intact': [0],
    'preterm_intact_model': ['trafficjam'],
    'preterm_intact_rate': [0],
    'n_stall': [6],
}

# create a list of all parameter combinations from above dict
# there will be as many list elements as the product of all list lengths above
and_params = [list(x) for x in it.product(*[[(p, v)
                                             for v in and_params[p]]
                                            for p in and_params])]

# include one of these param combinations in each of the simulations
n_list = range(6, 7, 1)
k_elong_stall_list = [0.1, 10]
combo_params = list()
for (n, k_elong_stall) in it.product(n_list, k_elong_stall_list):
    combo_params.append({
        # location of stall
        'x_stall': ','.join([str(x)
                             for x in range(400, 400 + n)]),
        # elongation rate at stall
        # do not allow it exceed 10s-1 which is the normal elongation rate
        'k_elong_stall': ','.join(
            [str(np.round(min(k_elong_stall * (n), 10), 1))] * (n)),
        # number of stalls
        'n_stall': n,
    })
# convert each param combination from dict to list of tuples
combo_params = [list(x.items()) for x in combo_params]

cleave_rate_list = list(np.round(0.001 * 2**(np.arange(-4.5,3,0.5)), 6))
mrna_cleave_params = list()
for cleave_rate in cleave_rate_list:
    mrna_cleave_params.append({
        # co-translational cleavage params
        # endocleavage occurs these many codons behind a-site
        # this is the simple cleavage model
        'l_cleave': 5,
        'k_cleave_no_hit': cleave_rate,
        'k_cleave_5_hit': cleave_rate,
        'k_cleave_3_hit': 0,
        'k_cleave_both_hit': 0,
        'cleave_model': 'simple',
        'cleave_rate': cleave_rate,
    })

    mrna_cleave_params.append({
        # co-translational cleavage params
        # endocleavage occurs these many codons behind a-site
        # this is the collision-stimulated cleavage model
        'l_cleave': 5,
        'k_cleave_no_hit': 0,
        'k_cleave_5_hit': cleave_rate,
        'k_cleave_3_hit': 0,
        'k_cleave_both_hit': cleave_rate,
        'cleave_model': 'hit5',
        'cleave_rate': cleave_rate,
    })

# convert each param combination from dict to list of tuples
mrna_cleave_params = [list(x.items()) for x in mrna_cleave_params]

# combine the 'and' parameters and parameter combinations from above
simcount = 0
temp = dict()
for params in it.product(and_params, combo_params, mrna_cleave_params):
    temp[simcount] = dict(it.chain.from_iterable(params))
    simcount += 1

# convert to pandas dataframe
input_params = pd.DataFrame.from_dict(temp, orient='index')

# sort  the paramters by these parameter combinations
input_params = input_params.sort_values(
    by=['n_stall', 'k_init', 'k_elong_stall']).reset_index(
    drop=True)

input_params.to_csv('sim.params.tsv', sep='\t')  # write to tab-delimited file
input_params.info()  # display the table of input parameters
