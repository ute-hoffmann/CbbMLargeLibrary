import pandas as pd

fitness_data = pd.read_csv("../../WeightingStrategies_fitnessCalc/result_lfcSE_Pearson_Weighted_sameWeight.tsv", sep="\t")
satu_combi = pd.read_csv("input/combi_satur_mutatios.tsv", names = ["sgRNA_target", "category"], sep="\t")
fitness_data = fitness_data.merge(satu_combi, on="sgRNA_target", how="left")

new_set_of_variants = pd.read_csv("input/2025-02-09_predict_new_variants.csv", sep=",")
new_set_of_variants = new_set_of_variants[new_set_of_variants["train_test_split"]==1][["mutant"]].drop_duplicates()
print(len(new_set_of_variants[["mutant"]]))

conditions = set(fitness_data["condition"])
print(conditions)
for cond in conditions:
	subset_fitness_data = fitness_data[fitness_data["condition"]==cond]
	
	single_amino_exchanges = subset_fitness_data[subset_fitness_data["category"]=="saturational"][["sgRNA_target", "norm"]].drop_duplicates()
	single_amino_exchanges["norm"] = single_amino_exchanges["norm"]-1.0
	
	additive_scores = []
	for variant in new_set_of_variants["mutant"]:
		additive_fitness=1.0
		exchanges = variant.split(":")
		if "M154A" in exchanges:
			additive_scores.append("NaN")
		else:
			for aa_ex in exchanges:
				additive_fitness += float(single_amino_exchanges[single_amino_exchanges["sgRNA_target"]==aa_ex]["norm"])
			additive_scores.append(additive_fitness)
	new_set_of_variants["additive_score_" + cond] = additive_scores
		
new_set_of_variants.to_csv("results/new_variants_additiveScores.csv", sep="\t")
