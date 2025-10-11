# create file to train for predicting new variants:
header = True
with open("all_variants_predict_new.csv", "w") as f:
	with open("saturational_proteinNPT_with_foldChange.csv") as f2:
		f.write("mutated_sequence,mutant,CL_N2,CL_O2,LD,fold_random_5\n")
		for line in f2:
			if header:
				header = False
				continue
			f.write(line)
	header = True
	with open("combinatorial_proteinNPT_with_random-foldChange.csv") as f2:
		for line in f2:
			if header:
				header = False
				continue
			f.write(line)
	
header = True
# create file to predict combinatorial:
with open("all_variants_predict_combinatorial.csv", "w") as f:
	with open("saturational_proteinNPT_with_foldChange.csv") as f2:
		f.write("mutated_sequence,mutant,CL_N2,CL_O2,LD,train_test_split\n")
		for line in f2:
			if header:
				header = False
				continue
			split_line = line.split(",")
			new_line = ""
			for i in range(len(split_line)-1):
				new_line += split_line[i] + "," 
			new_line += "0\n"
			f.write(new_line)
	header = True
	with open("combinatorial_proteinNPT_with_random-foldChange.csv") as f2:
		for line in f2:
			if header:
				header = False
				continue
			split_line = line.split(",")
			new_line = ""
			for i in range(len(split_line)-1):
				new_line += split_line[i] + "," 
			new_line += "1\n"
			f.write(new_line)

