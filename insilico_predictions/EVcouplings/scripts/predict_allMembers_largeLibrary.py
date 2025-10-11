"""
Script to predict fitness of all variants present in large Rubisco library
"""

import numpy as np
import numpy.random as npr
import pandas as pd
from evcouplings.couplings import CouplingsModel

model_path = "input/I_ab368e89ebe9421898840b5c014db5b3_TARGET_b0.3.model"
variants_of_interest_file_combi = "input/combinatorial_mutations.tsv"
variants_of_interest_file_satu = "input/saturational_mutations.tsv"

# Load a model for EVmutation
model = CouplingsModel(model_path)
# which positions are included in model?
positions = model.index_list

variants = {}
with open(variants_of_interest_file_combi) as f:
	for line in f:
		if line.strip("\n") == "WT":
			continue
		else: 
			variants[line.strip("\n")] = 0

with open(variants_of_interest_file_satu) as f:
	for line in f:
		if line.strip("\n") == "WT":
			continue
		else: 
			variants[line.strip("\n")] = 0

def scoring_func(mutation_list: list):
    variants_to_score = []
    for element in mutation_list:
        if int(element[1:-1]) in positions:
            variants_to_score.append((int(element[1:-1]), element[0], element[-1]))
        else: 
            return "NA"

    # Make the prediction, model loaded outside of function
    if len(variants_to_score) > 0:
        delta_E, _, _ = model.delta_hamiltonian(variants_to_score)
    else:
        delta_E = 0

    return delta_E

scored_variants = {}
for mut_variant in variants.keys():
	mut_list = mut_variant.split(",")
	scored_variants[mut_variant] = [scoring_func(mut_list)]

dataframe = pd.DataFrame(scored_variants).transpose()
dataframe.to_csv("variants_scored_EVcouplings.tsv", '\t', index=True)
