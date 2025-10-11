header = True
with open("all_variants_predict_new.csv") as f:
	with open("all_variants_predict_new_CL_N2.csv", "w") as CLN2:
		with open("all_variants_predict_new_CL_O2.csv", "w") as CLO2:
			with open("all_variants_predict_new_LD.csv", "w") as LD_f:
				CLN2.write("mutated_sequence,mutant,DMS_score,fold_random_5\n")
				CLO2.write("mutated_sequence,mutant,DMS_score,fold_random_5\n")
				LD_f.write("mutated_sequence,mutant,DMS_score,fold_random_5\n")
				for line in f:
					if header:
						header = False
						continue
					split_line = line.strip("\n").split(",")
					mutated_sequence = split_line[0]
					mutant = split_line[1]
					CL_N2 = split_line[2]
					CL_O2 = split_line[3]
					LD = split_line[4]
					fold = split_line[5]
					CLN2.write(mutated_sequence + "," + mutant + "," + CL_N2 + "," + fold + "\n")
					CLO2.write(mutated_sequence + "," + mutant + "," + CL_O2 + "," + fold + "\n")
					LD_f.write(mutated_sequence + "," + mutant + "," + LD + "," + fold + "\n")

header = True
with open("all_variants_predict_combinatorial.csv") as f:
	with open("all_variants_predict_combinatorial_CL_N2.csv", "w") as CLN2:
		with open("all_variants_predict_combinatorial_CL_O2.csv", "w") as CLO2:
			with open("all_variants_predict_combinatorial_LD.csv", "w") as LD_f:
				CLN2.write("mutated_sequence,mutant,DMS_score,train_test_split\n")
				CLO2.write("mutated_sequence,mutant,DMS_score,train_test_split\n")
				LD_f.write("mutated_sequence,mutant,DMS_score,train_test_split\n")
				for line in f:
					if header:
						header = False
						continue
					split_line = line.strip("\n").split(",")
					mutated_sequence = split_line[0]
					mutant = split_line[1]
					CL_N2 = split_line[2]
					CL_O2 = split_line[3]
					LD = split_line[4]
					fold = split_line[5]
					CLN2.write(mutated_sequence + "," + mutant + "," + CL_N2 + "," + fold + "\n")
					CLO2.write(mutated_sequence + "," + mutant + "," + CL_O2 + "," + fold + "\n")
					LD_f.write(mutated_sequence + "," + mutant + "," + LD + "," + fold + "\n")
header = True
with open("saturational_proteinNPT_with_foldChange.csv") as f:
	with open("saturational_proteinNPT_with_foldChange_CL_N2.csv", "w") as CLN2:
		with open("saturational_proteinNPT_with_foldChange_CL_O2.csv", "w") as CLO2:
			with open("saturational_proteinNPT_with_foldChange_LD.csv", "w") as LD_f:
				CLN2.write("mutated_sequence,mutant,DMS_score,fold_random_5\n")
				CLO2.write("mutated_sequence,mutant,DMS_score,fold_random_5\n")
				LD_f.write("mutated_sequence,mutant,DMS_score,fold_random_5\n")
				for line in f:
					if header:
						header = False
						continue
					split_line = line.strip("\n").split(",")
					mutated_sequence = split_line[0]
					mutant = split_line[1]
					CL_N2 = split_line[2]
					CL_O2 = split_line[3]
					LD = split_line[4]
					fold = split_line[5]
					CLN2.write(mutated_sequence + "," + mutant + "," + CL_N2 + "," + fold + "\n")
					CLO2.write(mutated_sequence + "," + mutant + "," + CL_O2 + "," + fold + "\n")
					LD_f.write(mutated_sequence + "," + mutant + "," + LD + "," + fold + "\n")
