#!/usr/bin/env python3
"""
Read in file with best trajectories and return file with point mutations, number of mutations, fitness score, sequence
"""

import os
import pandas as pd

path_file = 'results/2023-04-20_inSilicoEvolution_I/2023-04-21_bestSequencesTrajectory_I.csv'

# Load target sequence
sequence_file = "input/Gallionella_Rubisco_I.txt"
with open(sequence_file) as s:
    starting_sequence = s.read().strip()
data_holder = {"Fitness":[], "NumberMutations":[], "PointMutations":[], "Sequence":[], "Temperature":[]}

with open(path_file) as f:
	linecounter = 1
	for line in f:
		if linecounter == 1:
			linecounter += 1
			continue
		# structure of file: Fitness,Step,Replicate,Sequence\n
		line_split = line.strip("\n").split(",")
		score = float(line_split[0])
		sequence=line_split[-2]
		temperature=line_split[-1]
		mutation_list = ""
		mutation_counter = 0
		for i in range(len(starting_sequence)):
        		if starting_sequence[i].capitalize() != sequence[i].capitalize():
              			mutation_list += " " + starting_sequence[i] + str(i+1) + sequence[i]
              			mutation_counter += 1
		data_holder["Fitness"].append(score)
		data_holder["Sequence"].append(sequence)
		data_holder["NumberMutations"].append(mutation_counter)
		data_holder["PointMutations"].append(mutation_list)	
		data_holder["Temperature"].append(temperature)
	
df = pd.DataFrame(data_holder)
df.to_csv('results/2023-04-20_inSilicoEvolution_I/2023-04-21_I_Sequences_Fitness_Point-Mutations.csv', index=False)
