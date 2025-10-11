#!/bin/bash
#
# calculates ELBO values for variants present in large library, using DeepSequence VAE, CbbM(I) variant
#

LOGFILE="results/2023-04-20_inSilicoEvolution_I/2024-07-30_logfile_I_all_combinations.txt"
date > $LOGFILE 2>&1 # Log start time

THEANO_FLAGS='floatX=float32,device=cuda' python source/DeepSequence_predict_allMembers_largeLibrary_GallionellaI.py >> $LOGFILE 2>&1

date >> $LOGFILE 2>&1
