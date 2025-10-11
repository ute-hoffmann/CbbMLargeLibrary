import pandas as pd

fitness_data = pd.read_csv("../../WeightingStrategies_fitnessCalc/result_lfcSE_Pearson_Weighted_sameWeight.tsv", sep="\t")
satu_combi = pd.read_csv("../../input/combi_satur_mutatios.tsv", names = ["sgRNA_target", "category"], sep="\t")
fitness_data = fitness_data.merge(satu_combi, on="sgRNA_target", how="left")

conditions = set(fitness_data["condition"])
for cond in conditions:
	subset_fitness_data = fitness_data[fitness_data["condition"]==cond]
	
	single_amino_exchanges = subset_fitness_data[subset_fitness_data["category"]=="combiANDsatur"][["sgRNA_target", "norm"]].drop_duplicates()
	single_amino_exchanges["norm"] = single_amino_exchanges["norm"]-1.0
	
	combinatorial_variants = subset_fitness_data[subset_fitness_data["category"]=="combinatorial"][["sgRNA_target", "norm"]].drop_duplicates()
	
	additive_scores = []
	for variant in combinatorial_variants["sgRNA_target"]:
		additive_fitness=1.0
		exchanges = variant.split(",")
		if "M154A" in exchanges:
			additive_scores.append("NaN")
		else:
			for aa_ex in exchanges:
				additive_fitness += float(single_amino_exchanges[single_amino_exchanges["sgRNA_target"]==aa_ex]["norm"])
			additive_scores.append(additive_fitness)
	combinatorial_variants["additive_score"] = additive_scores
		
	combinatorial_variants.to_csv("results/combinatorial_variants_additiveScores_" + cond + ".csv", sep="\t")
