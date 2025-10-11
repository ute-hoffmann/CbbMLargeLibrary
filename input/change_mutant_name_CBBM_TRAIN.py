with open("CBBM_TRAIN.csv") as MSA_predict:
	with open("CBBM_TRAIN_only_predictAndMutant.tsv", "w") as new_f:
		header = True
		for line in MSA_predict:
			if header:
				header = False
				new_f.write("sgRNA_target\tMSA_Transform\n") # in .Rmd, mutants are called sgRNA_target
				continue
			split_line = line.strip("\n").split(",")
			mut_colon = split_line[1]
			MSA_trans = split_line[-1]
			mutant_name = ""
			for part in mut_colon.split(":"):
				if mutant_name:
					mutant_name += ","
					mutant_name += part
				else:
					mutant_name = part
			new_f.write(mutant_name + "\t" + MSA_trans + "\n")
