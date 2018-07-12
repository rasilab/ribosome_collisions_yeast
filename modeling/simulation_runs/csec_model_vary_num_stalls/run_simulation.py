#!/usr/bin/env python
# SBATCH --mem=8000

"""Program to run simulation with a specific parameter combination

This program

- reads the table of kinetic parameter combinations that we want to simulate
  in this run
- sets simulation time, recording intervals
- gets the specific simulation to run as a number from the command line
- creates a BioNetGen file for that simulation by exporting the PySB model
- runs NFSim simulation using the BioNetGen as input

"""


import sys
# insert path to the model definition file
sys.path.insert(0, '../../')

import numpy as np
import pandas as pd
from tasep import Tasep
from pysb.export import export
import subprocess as sp
import os

# extract the parameter combination for this simulation
# simindex = int(sys.argv[1])

input_params = pd.read_table('sim.params.tsv', index_col=0)

for simindex in range(90):
    params_for_this_job = input_params.loc[simindex].to_dict()
    # instantiate model with the correct parameters
    model = Tasep(**params_for_this_job)

#     # do not track reactions that occur too often and are not
#     # critical for our interpretation
#     for rule in model.rules:
#         if rule.name.startswith('deadenylation'):
#             rule.tag = False
#         elif rule.name.startswith('elongation'):
#             rule.tag = False
#         elif rule.name.startswith('exonucleolysis'):
#             rule.tag = False
#         else:
#             rule.tag = True
#         if rule.name == ('exonucleolysis_5end_0'):
#             rule.tag = True

    outdir = './output'
    bnglfile = f'{outdir}/tasep_{simindex}.bngl'
    xmlfile = bnglfile.replace('.bngl', '.xml')
    gdatfile = bnglfile.replace('.bngl', '.gdat')
    rxnfile = bnglfile.replace('.bngl', '.rxns.tsv')
    paramsfile =  bnglfile.replace('.bngl', '.params.tsv')

    # write all model parameters to a separate file
    with open(paramsfile, 'w') as file:
        file.write('parameter\tvalue\n')
        for param in model.parameters:
            file.write(f'{param.name}\t{param.value}\n')
    # compress the params file
    sp.run(['gzip', '-f', paramsfile.split('/')[-1]], cwd=outdir)

# # write BNGL file
# with open(bnglfile, 'w') as file:
#     file.write(export(model, 'bngl'))
# 
# # convert BNGL file to XML for NFSim input
# sp.run(['BNG2.pl', '--xml', '--outdir', outdir, bnglfile])
# 
# # parameters for NFsim
# equilibrium_time = 0  # seconds
# tstop = str(1000000)  # seconds
# maxcputime = str(100 * 60)  # seconds
# osteps = str(100)  # number of samples
# seed = str(111)  # random number initial seed
# gml = str(1000000)  # max num of mol allowed in simulation
# utl = '3'  # max number of bonds to traverse during simulation
# network = '-connect'  # whether to infer reaction network connectivity
# 
# # print NFSim command
# nfsim_command = [
#     'NFsim', '-xml', xmlfile, '-sim', tstop, '-oSteps', osteps,
#     '-seed', seed, '-o', gdatfile, '-rxnlog', rxnfile,
#     '-utl', utl,
#     '-gml', gml, '-maxcputime', maxcputime,
#     network
# ]
# print(' '.join(nfsim_command))
# 
# # simlate with NFSim
# sp.run(nfsim_command)
# 
# # calculate statistics of various reaction firings
# sp.run(['Rscript', '../../get_mrna_lifetime_and_psr.R', str(simindex), outdir])
# 
# # remove the files that are too big or are not needed for analysis
# # comment these lines out for troubleshooting
# os.remove(gdatfile)
# os.remove(bnglfile)
# os.remove(xmlfile)
# os.remove(rxnfile)
