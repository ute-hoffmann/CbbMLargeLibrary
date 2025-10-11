with open("2025-02-03_ProteinNPT_CBBM_TRAIN_fitness1_fitness2_fitness3_zero_shot_fitness_predictions_train_test_split_embed_MSA_Transformer_head_CNN_aug_auxiliary_froz_True_drop_0.0_val_False_2025-02-03_trainMultipleObjectives_fold-1.csv") as ProteinNPT_predict:
	with open("ProteinNPT_mutant_fitnessPredictions.tsv", "w") as new_f:
		header = True
		for line in ProteinNPT_predict:
			if header:
				header = False
				new_f.write("sgRNA_target\tCL_N2\tCL_O2\tLD\n") # in .Rmd, mutants are called sgRNA_target
				continue
			split_line = line.strip("\n").split(",")
			mut_colon = split_line[0]
			CL_N2 = split_line[1]
			CL_O2 = split_line[3]
			LD = split_line[5]
			mutant_name = ""
			for part in mut_colon.split(":"):
				if mutant_name:
					mutant_name += ","
					mutant_name += part
				else:
					mutant_name = part
			new_f.write(mutant_name + "\t" + CL_N2 + "\t" + CL_O2 + "\t" + LD + "\n")
