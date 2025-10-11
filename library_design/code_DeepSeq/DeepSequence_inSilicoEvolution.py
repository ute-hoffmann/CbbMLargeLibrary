"""
written by Ute H., adapted to Python 2.7

Loosely based on https://github.com/Asplund-Samuelsson/furee/blob/master/source/in_silico_evolution.py
EVmutation code from https://github.com/fhalab/MLDE/blob/39e9edccc346119a834c62677d16e39dd49dfbdb/code/zero_shot/zero_shot_predictor.py#L29
Adjusted in 10/2022 to use DeepSequence predictions instead of EVmutation.

Perform Metropolis-Hastings MCMC with DeepSequence model as scoring function to create zero-shot predictions
i.e.: combine approach of Biswas et al., 2021 (Metropolis algorithm to perform in silico directed evolution) and Wittmann et al., 2021 (DeepSequence, zero-shot predictions) using DeepSequence (Riesselman et al., 2018)
"""

import numpy as np
import numpy.random as npr
import argparse
import pandas as pd

import os
import sys
os.environ["THEANO_FLAGS"] = "floatX=float32,device=cuda"

# for local machine WORKING_DIR="../../DeepSequence/" # Put in the DeepSequence directory here
WORKING_DIR="../MLDE/deep_sequence/DeepSequence"
N_ELBO_SAMPLES=400
module_path = os.path.abspath(WORKING_DIR)
if module_path not in sys.path:
        sys.path.append(module_path)
from DeepSequence import helper
from DeepSequence.model import VariationalAutoencoder

from collections import defaultdict
#from tqdm.autonotebook import tqdm

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


# Read arguments from the commandline
parser = argparse.ArgumentParser()

# Required input: Infile and project output (positional)
parser.add_argument(
    'infile', type=str,
    help='Text file with one sequence on one line. Original fasta to start from.' # equivalent to wt.fasta for vae_inference ? 
)

parser.add_argument(
    'outfile', type=str,
    help='Tab-delimited file with sampled sequences and scores.'
)

parser.add_argument(
    '-s', '--steps', type=int, default=10,
    help='Number of MCMC steps [10].'
)
parser.add_argument(
    '-t', '--trust', type=int, default=7,
    help='Trust radius [7].'
)
parser.add_argument(
    '-T', '--temperature', type=float, default=0.1,
    help='MCMC temperature; lower is slower [0.1].'
)
parser.add_argument(
    '-R', '--ratio', action='store_true', default=False,
    help='Use ratio instead of difference for sequence proposal rejection.'
)
parser.add_argument( # adjust to DeepSequence model
    '-DS', "--DeepSequenceModelPath", type=str,
    help='filePrefix of model in working_dir, created with code from Hsu et al., 2022, train_vae.sh'
)
parser.add_argument( # a2m file
    '-a2m', "--a2mFile", type=str,
    help='Path to .a2m file, as created by EVcouplings webserver'
)

parser.add_argument( # where to save output
    '-o', "--output_dir", type=str,
    help="directory in which DeepSequence VAE assumes model parameters and will save output, if there is any"
)

# Parse arguments
args = parser.parse_args()

sequence_file = args.infile
outfile = args.outfile
steps = args.steps
trust = args.trust
temperature = args.temperature
ratio = args.ratio
DeepSequenceModelPath = args.DeepSequenceModelPath
alignmentFile = args.a2mFile
output_dir = args.output_dir


# Load target sequence
with open(sequence_file) as s:
    starting_sequence = s.read().strip()

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

# prepare function to propose new sequences
positions = data_helper.uniprot_focus_cols_list # only use positions which have good coverage in alignment, only these can be used for predictions. Plus: use Uniprot numbering - starts with position 1, not with index 0 !
aa_dict = data_helper.aa_dict
letters_sorted = []
for i in aa_dict.keys():
    letters_sorted.append(i)

def propose(sequence):
    """
    Given a string, return a proposed mutated string.
    The proposed mutant is generated as follows:
    - Pick a position given the pick probabilities for each position (uniform distribution)
    - Given that position,
    pick a letter given the pick probabilities for each letter at that position (uniform distribution)
    that is not the existing letter.
    :param sequence: The sequence to propose a new mutation on.
    :returns: A string.
    """
    npr.seed(None) # seed is set in DeepSequence helper.py function, with set seed: no real permutations anymore!
    if len(sequence) == 0:
        raise ValueError(
            "sequence passed into `propose` must be at least length 1."
        )

    # define uniform probability for positions and letters
    pos_prob = np.array([1.0/len(positions)] * len(positions))
    pwm = np.tile(np.array([[1.0/len(aa_dict)] * len(aa_dict)]), (len(sequence), 1))

    position = positions[np.argmax(npr.multinomial(1, pos_prob))]
    new_sequence = ""
    for i, letter in enumerate(sequence):
        if (position-1) != i: # position - i, because position does not refer to index, but position
            new_sequence += letter
        else:
            letter_idx = aa_dict[letter]
            pwm[i, letter_idx] = 0 # set prob. to return to original amino acid to 0
            new_letter_idx = np.argmax(npr.multinomial(1, pwm[i, :]))
            new_letter = letters_sorted[new_letter_idx]
            new_sequence += new_letter
    return new_sequence

# from Hsu et al., 2022
def seq2mutation_fromwt(seq, wt, ignore_gaps=False, sep=':', offset=1,
        focus_only=True): # offset of 1 guarantees to get positions, not indices, which makes function compatible with uniprot_focus_cols numbering
    mutations = []
    for i in range(offset, offset+len(seq)):
        if ignore_gaps and (seq[i-offset] == '-'):
            continue
        if wt[i-offset].islower() and focus_only:
            continue
        if seq[i-offset].upper() != wt[i-offset].upper():
            mutations.append((i, wt[i-offset].upper(), seq[i-offset].upper()))
    return mutations

# Define sequence scoring function # adapted from vae_inference.py from Hsu et al., 2022
def scoring_func(mut_sequence):
    mut_tups = seq2mutation_fromwt(seq=mut_sequence, wt=starting_sequence)
    mut_tups = [t for t in mut_tups if t[0] in focuscols]
    delta_elbos = data_helper.delta_elbo(vae_model, mut_tups, N_pred_iterations=N_ELBO_SAMPLES) 
    return delta_elbos

# Patch the is_accepted function to work with multiple scores
def ia_patch_diff(best, candidate, temperature):
    # Compare candidate based on difference in scores
    # basically exactly like ju.sampler.is_accepted, only minor difference (always calculation of np.random.uniform(0,1)
    # best: float
    # candidate: float
    # temperature: float
    # returns bool
    c = np.exp((candidate - best) / temperature)
    p = np.random.uniform(0, 1)
    return np.all(c >= p)

def ia_patch_ratio(best, candidate, temperature):
     # Compare candidate based on ratio of scores
    # best: float
    # candidate: float
    # temperature: float
    # returns bool
     c = np.exp(np.log(candidate / best) / temperature)
     p = np.random.uniform(0, 1)
     return np.all(c >= p)

# Use selected patch (difference or ratio of scores)
if not ratio:
    is_accepted = ia_patch_diff
else:
    is_accepted = ia_patch_ratio


# from https://github.com/ElArkk/jax-unirep/blob/master/jax_unirep/sampler.py
def hamming_distance(s1, s2):
    """
    s1: string
    s2: string
    Return hamming distance between two strings of the same length."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def sample_one_chain(starter_sequence, n_steps, scoring_func, is_accepted_kwargs = {}, trust_radius = 7, propose_kwargs = {},):
    """
    Return one chain of MCMC samples of new sequences.
    Given a `starter_sequence` (str),
    this function will sample one chain of protein sequences,
    scored using a user-provided `scoring_func`.
    Design choices made here include the following.
    Firstly, we record all sequences that were sampled,
    and not just the accepted ones.
    This behaviour differs from other MCMC samplers
    that record only the accepted values.
    We do this just in case sequences that are still "good"
    (but not better than current) are rejected.
    The effect here is that we get a cluster of sequences
    that are one-apart from newly accepted sequences.
    Secondly, we check the Hamming distance
    between the newly proposed sequences and the original.
    This corresponds to the "trust radius"
    specified in the [jax-unirep paper](https://doi.org/10.1101/2020.01.23.917682).
    If the hamming distance > trust radius,
    we reject the sequence outright.
    A dictionary containing the following key-value pairs are returned:
    - "sequences": All proposed sequences.
    - "scores": All scores from the scoring function.
    - "accept": Whether the sequence was accepted as the new 'current sequence'
        on which new sequences are proposed.
    This can be turned into a pandas DataFrame.
    ### Parameters
    - `starter_sequence`: The starting sequence.
    - `n_steps`: Number of steps for the MC chain to walk.
    - `scoring_func`: Scoring function for a new sequence.
        It should only accept a string `sequence`.
    - `is_accepted_kwargs`: Dictionary of kwargs to pass into
        `is_accepted` function.
        See `is_accepted` docstring for more details.
    - `trust_radius`: Maximum allowed number of mutations away from
        starter sequence.
    - `propose_kwargs`: Dictionary of kwargs to pass into
        `propose` function.
        See `propose` docstring for more details.
    - `verbose`: Whether or not to print iteration number
        and associated sequence + score. Defaults to False
    ### Returns
    A dictionary with `sequences`, `accept` and `score` as keys.
    """
    current_sequence = starter_sequence
    current_score = scoring_func(mut_sequence=starter_sequence) 

    chain_data = defaultdict(list)
    chain_data["sequences"].append(current_sequence)
    chain_data["scores"].append(current_score)
    chain_data["accept"].append(True)

    for i in range(n_steps):
        new_sequence = propose(current_sequence)
        new_score = scoring_func(mut_sequence=new_sequence)

        default_is_accepted_kwargs = {"temperature": 0.1}
        default_is_accepted_kwargs.update(is_accepted_kwargs)
        accept = is_accepted(
            best=current_score,
            candidate=new_score,
            **default_is_accepted_kwargs)

        # Check hamming distance
        if hamming_distance(starter_sequence, new_sequence) > trust_radius: ## where do we get this from?
            accept = False

        # Determine acceptance
        if accept:
            current_sequence = new_sequence
            current_score = new_score

        # Record data.
        chain_data["sequences"].append(new_sequence)
        chain_data["scores"].append(new_score)
        chain_data["accept"].append(accept)
    chain_data["scores"] = np.hstack(chain_data["scores"])
    return chain_data


# Perform sampling
sampled_sequences = sample_one_chain(
    starting_sequence, n_steps=steps, scoring_func=scoring_func,
    trust_radius=trust, is_accepted_kwargs={'temperature': temperature}
)

# Extract the scores
scores = sampled_sequences.pop('scores')

# Create data frame
sampled_seqs_df = pd.DataFrame(sampled_sequences)

# Add score columns
sampled_seqs_df = pd.concat([
    sampled_seqs_df,
    pd.DataFrame({"Scores": scores})],
    axis=1
)

# Re-order columns
cols = sampled_seqs_df.columns.values
cols = list(cols[cols != 'accept']) + ['accept']

sampled_seqs_df = sampled_seqs_df[cols]

# Add step numbers
sampled_seqs_df['step'] = list(range(0, steps + 1))

# Save sampled sequences
sampled_seqs_df.to_csv(outfile, '\t', index=False)
