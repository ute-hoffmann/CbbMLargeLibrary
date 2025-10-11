#with open("2025-01-13_proteinNPT_predictions_CLN2_allCombinatorial_basedOnSaturational_CL_N2.csv") as predict:
with open("2025-01-23_predictions_allCombinatorial_basedOnCopleteSaturational_CL_N2.csv") as predict:
#	with open("2025-01-13_proteinNPT_CLN2_predictions.tsv", "w") as new_f:
	with open("2025-01-23_proteinNPT_CLN2_predictions.tsv", "w") as new_f:
		header = True
		for line in predict:
			if header:
				header = False
				new_f.write("sgRNA_target\tproteinNPT_predict\n") # in .Rmd, mutants are called sgRNA_target
				continue
			split_line = line.strip("\n").split(",")
			mut_colon = split_line[0]
			prot_npt = split_line[1]
			mutant_name = ""
			for part in mut_colon.split(":"):
				if mutant_name:
					mutant_name += ","
					mutant_name += part
				else:
					mutant_name = part
			new_f.write(mutant_name + "\t" + prot_npt + "\n")
