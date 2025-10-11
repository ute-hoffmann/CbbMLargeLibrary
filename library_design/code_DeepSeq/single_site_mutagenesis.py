'''
Infers average ELBO values from DeepSequence VAE models.
Based on open source code from DeepSequence repo.
'''
import argparse
import numpy as np
import os
import sys

WORKING_DIR="../MLDE/deep_sequence/DeepSequence" # Put in the DeepSequence directory here
N_ELBO_SAMPLES=400

module_path = os.path.abspath(WORKING_DIR)
if module_path not in sys.path:
        sys.path.append(module_path)

from DeepSequence.model import VariationalAutoencoder
from DeepSequence import helper

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a2m", "--a2mFile", type=str, help="Path to .a2m file, as created by EVcouplings webserver")
    parser.add_argument("-f", "--filename_singleSite", type=str, help="output file name")
    parser.add_argument('-o', "--output_dir", type=str, help="directory in which DeepSequence VAE assumes model parameters and will save output, if there is any")
    parser.add_argument("-p", "--model_prefix", type=str, help="filePrefix of model in output_dir")

    args = parser.parse_args()

    # create DeepSequence data_helper
    data_helper = helper.DataHelper(working_dir=args.output_dir, alignment_file=args.a2mFile, calc_weights=False,)

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
        working_dir                    =   args.output_dir,
        )
    vae_model.load_parameters(file_prefix=args.model_prefix)
    data_helper.single_mutant_matrix(model=vae_model, filename_prefix=args.filename_singleSite)


if __name__ == "__main__":
    main()
