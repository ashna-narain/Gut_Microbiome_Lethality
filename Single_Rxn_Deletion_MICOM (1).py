import os
import warnings
import re

import math
import sympy
import scipy
import cobra
import csv
import numpy as np
import pandas as pd
from micom import Community
from cobra.io import read_sbml_model
from micom.media import minimal_medium
from concurrent.futures import ThreadPoolExecutor, as_completed

file = open('/data/ashna/Gut_Microbiome/Single_Lethals_Communities_List_Write.csv', 'w')

def process_pair(pair):

    orgs = pair.split(',')
    file_1 = f'/data/ashna/Gut_Microbiome/AGORA2_Models/{orgs[0]}.xml'
    file_2 = f'/data/ashna/Gut_Microbiome/AGORA2_Models/{orgs[1]}.xml'

    model_1 = read_sbml_model(file_1)
    model_2 = read_sbml_model(file_2)

    table = {'id': ['model_1', 'model_2'], 'file': [file_1, file_2]}
    data = pd.DataFrame(table)
    taxonomy = data

    com = Community(taxonomy, progress=False)
    com.medium = minimal_medium(com, 1, min_growth=0.1)
    sol_WT = com.cooperative_tradeoff(fraction=0.7, fluxes=False)
    sol_df = sol_WT.members

    # Collect basic info
    row = [
        orgs[0], orgs[1],
        len(com.reactions), len(com.metabolites),
        sol_WT.growth_rate,
        sol_df._get_value('model_1', 'growth_rate'),
        sol_df._get_value('model_2', 'growth_rate'),
        sol_df._get_value('medium', 'reactions')
    ]

    # Perform single reaction deletion
    single_lethal_reactions_com = []
    for rxn in com.reactions:
        rxn_id = rxn.id
        if rxn_id not in list(com.medium.keys()):
            com_del = com
            with com_del:
                com_del.reactions.get_by_id(rxn_id).knock_out()
                try:
                    sol_del = com_del.cooperative_tradeoff(fraction = 0.7, fluxes=False)
                    sol_del_df = sol_del.members

                    if sol_del.growth_rate <= 0.1 * sol_WT.growth_rate or sol_del_df._get_value('model_1', 'growth_rate') <= sol_del.growth_rate * 0.1 or  sol_del_df._get_value('model_2', 'growth_rate') <= sol_del.growth_rate * 0.1:
                        single_lethal_reactions_com.append(rxn.id)

                except ValueError:
                    print('Error', rxn_id)
                    pass
    print(pair)
    print(single_lethal_reactions_com)
    file.write(pair, single_lethal_reactions_com)
    return row, pair, single_lethal_reactions_com


org_list_file = open('/data/ashna/Gut_Microbiome/Species.csv')
org_names = [(i.split())[0] for i in org_list_file.readlines()]

pairwise_orgs = []
for i in range(len(org_names)):
    for j in range(i + 1, len(org_names)):
        pairwise_orgs.append(org_names[i] + ',' + org_names[j])

        
# Main code block to run in parallel
rows = []
lethals_pd = pd.DataFrame(columns = pairwise_orgs, index = range(8000))

completed_pd = pd.read_csv("/data/ashna/Gut_Microbiome/Single_Reaction_Lethals_Communities_with_Min_Growths.csv")
completed_pd = completed_pd.drop('Unnamed: 0', axis=1)

incomplete = []
for i in completed_pd.columns:
    if (completed_pd[i].dropna()).empty:
        print(i)
        incomplete.append(i)

with ThreadPoolExecutor() as executor:
    futures = {executor.submit(process_pair, pair): pair for pair in incomplete}
    for future in as_completed(futures):
        row, pair, single_lethal_reactions = future.result()
        rows.append(row)
        lethals_pd[pair] = pd.Series(single_lethal_reactions)

        
lethals_pd.to_csv('/data/ashna/Gut_Microbiome/Single_Lethals_Communities_2.csv')
print(rows)