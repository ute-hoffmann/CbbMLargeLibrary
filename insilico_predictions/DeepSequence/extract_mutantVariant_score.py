import pandas as pd

mutant_score_dict = {}
header = True
with open("2023-04-24_I_Sequences_Fitness_Point-Mutations.csv") as f:
	for line in f:
		if header:
			header = False
			continue
		fitness = line.split(",")[0]
		mutant = line.split(",")[2]
		mutant_name = ""
		for part in mutant.split(" "):
			if part == "":
				continue
			mutant_name += part + ","
		mutant_name = mutant_name[:-1]
		mutant_score_dict[mutant_name] = [fitness]

pd.DataFrame(mutant_score_dict).transpose().to_csv("quadruple_variants_scored_DeepSeq.tsv", '\t', index=True)
