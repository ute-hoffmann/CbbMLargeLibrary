#!/bin/bash
#
# calculates ELBO values for single mutants, using DeepSequence VAE, Gallionella Rubisco
#
mkdir -p results/2023-04-20_inSilicoEvolution_I

LOGFILE="results/2023-04-20_inSilicoEvolution_I/2023-04-20_logfile_I_single-site.txt"
date > $LOGFILE 2>&1 # Log start time

THEANO_FLAGS='floatX=float32,device=cuda' python source/single_site_mutagenesis.py \
    -a2m input/TARGET_b0.3_I_.a2m \
    -f 2023-04-20_single_mutant_matrix_Gallionella_Rubisco_I\
    -o results/2023-04-20_inSilicoEvolution_I/ \
    -p GallionellaRubisco_I >> $LOGFILE 2>&1

date >> $LOGFILE 2>&1
