"""
written by Ute H., adapted to Python 2.7
Script to predict fitness of all variants present in large Rubisco library using DeepSequence
"""

import numpy as np
import numpy.random as npr
import argparse
import pandas as pd

import os
import sys
#os.environ["THEANO_FLAGS"] = "floatX=float32,device=cuda"

WORKING_DIR="../../../DeepSequence" # Put in the DeepSequence directory here
N_ELBO_SAMPLES=400
module_path = os.path.abspath(WORKING_DIR)
if module_path not in sys.path:
        sys.path.append(module_path)
from DeepSequence import helper
from DeepSequence.model import VariationalAutoencoder

from collections import defaultdict

# default parameters for VariationalAutencoder
model_params = {
    "bs"                :   100,
    "encode_dim_zero"   :   1500,
    "encode_dim_one"    :   1500,
    "decode_dim_zero"   :   100,
    "decode_dim_one"    :   500,
    "n_latent"          :   30,
    "logit_p"           :   0.001,
    "sparsity"          :   "logit",
    "final_decode_nonlin":  "sigmoid",
    "final_pwm_scale"   :   True,
    "n_pat"             :   4,
    "r_seed"            :   12345,
    "conv_pat"          :   True,
    "d_c_size"          :   40
}

DeepSequenceModelPath = "GallionellaRubisco_I"
alignmentFile = "input/TARGET_b0.3_I_.a2m"
output_dir = "results/" # where do I find params?
variants_of_interest_file_combi = "input/combinatorial_mutations.tsv"
variants_of_interest_file_satu = "input/saturational_mutations.tsv"

variants = {}
#with open(variants_of_interest_file_combi) as f:
#	for line in f:
#		if line.strip("\n") == "WT":
#			continue
#		else: 
#			variants[line.strip("\n")] = 0

with open(variants_of_interest_file_satu) as f:
	for line in f:
		if line.strip("\n") == "WT":
			continue
		else: 
			variants[line.strip("\n")] = 0

# create DeepSequence data_helper
data_helper = helper.DataHelper(working_dir=output_dir, alignment_file=alignmentFile, calc_weights=False,)

# create VAE + load parameters
vae_model   = VariationalAutoencoder(data_helper,
    batch_size                     =   model_params["bs"],
    encoder_architecture           =   [model_params["encode_dim_zero"],
                                            model_params["encode_dim_one"]],
    decoder_architecture           =   [model_params["decode_dim_zero"],
                                            model_params["decode_dim_one"]],
    n_latent                       =   model_params["n_latent"],
    logit_p                        =   model_params["logit_p"],
    sparsity                       =   model_params["sparsity"],
    encode_nonlinearity_type       =   "relu",
    decode_nonlinearity_type       =   "relu",
    final_decode_nonlinearity      =   model_params["final_decode_nonlin"],
    final_pwm_scale                =   model_params["final_pwm_scale"],
    conv_decoder_size              =   model_params["d_c_size"],
    convolve_patterns              =   model_params["conv_pat"],
    n_patterns                     =   model_params["n_pat"],
    random_seed                    =   model_params["r_seed"],
    working_dir                    =   output_dir,
    )
vae_model.load_parameters(file_prefix=DeepSequenceModelPath)
focuscols = data_helper.uniprot_focus_cols_list # needed later for scoring function

def scoring_func(mutation_list):
    mut_seq = ""
    variants_to_score = []
    for element in mutation_list:
        if int(element[1:-1]) in focuscols:
           variants_to_score.append((int(element[1:-1]), element[0], element[-1]))
        else:
            return "NA"
            
    # Make the prediction, model loaded outside of function
    if len(variants_to_score) > 0:
    	delta_elbos = data_helper.delta_elbo(vae_model, variants_to_score, N_pred_iterations=N_ELBO_SAMPLES) 
    	return delta_elbos
    else:
    	return "NA"

scored_variants = {}
n = 1
print(len(variants.keys()))
for mut_variant in variants.keys():
	print(n)
	print(mut_variant)
	n += 1
	if n == 40:
		break
	mut_list = mut_variant.split(",")
	scored_variants[mut_variant] = [scoring_func(mut_list)]
	print(scored_variants[mut_variant])

dataframe = pd.DataFrame(scored_variants).transpose()
dataframe.to_csv("variants_scored_EVcouplings.tsv", '\t', index=True)
