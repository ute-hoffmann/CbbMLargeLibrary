import random

# load original sequence
cbbM_seq = ""
with open("Gallionella_Rubisco_I.txt") as f:
	for line in f:
		cbbM_seq = line

list_mutant_names = []
with open("list_new_sequences.txt") as seq_f:
	# header is "mutated_sequence,mutant,CL_N2,CL_O2,LD,train_test_split\n"
	with open("data_for_new_variants.csv", "w") as new_f:
		for new_seq in seq_f:
			mutant_name = ""
			for pos in range(1,len(cbbM_seq)+1):
				if cbbM_seq[pos-1] != new_seq[pos-1]:
					WT_aa = cbbM_seq[pos-1]
					mut_aa = new_seq[pos-1]
					mutant_name += WT_aa + str(pos) + mut_aa + ":"
			if mutant_name == "":
				print("something weird with sequence "+ new_seq)
				continue
			new_f.write(new_seq.strip("\n") + "," + mutant_name[:-1] + "," + str(random.random()) + "," + str(random.random()) + "," + str(random.random()) + ",1\n")
			list_mutant_names.append(mutant_name)
			
print("Number mutant names: " + str(len(list_mutant_names)))
print("Unique mutant names: " + str(len(set(list_mutant_names))))
			
mutants_already_in_file = []

with open("2025-02-09_predict_new_variants.csv", "w") as new_f:
	with open("all_variants_predict_new.csv") as old_f:
		for line in old_f:
			new_f.write(line)
			mutant = line.split(",")[1]
			mutants_already_in_file.append(mutant)
	with open("data_for_new_variants.csv") as old_f:
		for line in old_f:
			mutant = line.split(",")[1]
			if mutant in mutants_already_in_file:
				print(mutant + " already in file")
				continue
			mutants_already_in_file.append(mutant)	
			new_f.write(line)
