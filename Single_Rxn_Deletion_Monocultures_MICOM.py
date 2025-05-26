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
from cobra.medium import minimal_medium
from cobra.flux_analysis import (single_gene_deletion, single_reaction_deletion, 
                                 double_gene_deletion, double_reaction_deletion)
from concurrent.futures import ThreadPoolExecutor, as_completed

def process(org):

    model = read_sbml_model(f'/data/ashna/Gut_Microbiome/AGORA2_Models/{org}.xml')

    model.medium = minimal_medium(model, 1)
    
    sol_WT = model.slim_optimize()

    single_lethal_reactions_model = []
    
    # Perform single reaction deletion
    for rxn in model.reactions:
        if rxn not in model.exchanges:
            rxn_id = rxn.id
            if rxn_id not in list(model.medium.keys()):
                model_del = model
                with model_del:
                    model_del.reactions.get_by_id(rxn_id).knock_out()
                    try:
                        sol_del = model_del.slim_optimize()
                        if sol_del <= 0.1 * sol_WT:
                            single_lethal_reactions_model.append(rxn.id)
                    except ValueError:
                        print('Error', rxn_id)
                        pass
    
    print(single_lethal_reactions_model)
    return org, single_lethal_reactions_model

org_list_file = open('/data/ashna/Gut_Microbiome/Species.csv')
org_names = [(i.split())[0] for i in org_list_file.readlines()]
        
# Main code block to run in parallel
lethals_pd = pd.DataFrame(columns = org_names, index = range(1000))

with ThreadPoolExecutor() as executor:
    futures = {executor.submit(process, org): org for org in org_names}
    for future in as_completed(futures):
        org, single_lethal_reactions = future.result()
        lethals_pd[org] = pd.Series(single_lethal_reactions)
        print(len(single_lethal_reactions))
        
lethals_pd.to_csv('/data/ashna/Gut_Microbiome/Single_Lethals_Monocultures.csv')