import random

# load original sequence
cbbM_seq = ""
with open("Gallionella_Rubisco_I.txt") as f:
	for line in f:
		cbbM_seq = line

fold_random_dict = {0:0, 1:0, 2:0, 3:0, 4:0}
with open("combinatorial_for_proteinNPT.tsv") as f:
	with open("combinatorial_proteinNPT_with_random-foldChange.csv", "w") as f2:
		f2.write("mutated_sequence,mutant,CL_N2,CL_O2,LD,fold_random_5\n")
		header = True
		for line in f:
			if header:
				header = False
				continue
			if "STOP" in line:
				continue
			split_line = line.strip("\n").split("\t")
			aa_combination = split_line[0]
			CL_N2 = split_line[1]
			CL_O2 = split_line[2]
			LD = split_line[3]
			
			original_aas = []
			positions = []
			new_aas = []
			
			mutant = ""
			for aa_exchange in aa_combination.strip('"').split(","):
				original_aas.append(aa_exchange[0])
				new_aas.append(aa_exchange[-1])
				positions.append(int(aa_exchange[1:-1]))
				if mutant:
					mutant += ":" + aa_exchange
				else:
					mutant += aa_exchange
			
			new_seq = ""
			for pos in range(1,len(cbbM_seq)+1): # STOP!
				if pos in positions:
					index = positions.index(pos)
					if original_aas[index] != cbbM_seq[pos-1]:
						print(original_aas[index])
						print(cbbM_seq[pos-1])
						print(pos)
						break
					new_seq += new_aas[index]
				else:
					new_seq += cbbM_seq[pos-1]
			if len(new_seq) != len(cbbM_seq):
				print("something is weird for " + aa_combination + " and new seq " + new_seq)
			else:
				fold_random = random.randint(0,4)
				fold_random_dict[fold_random] += 1
				f2.write(new_seq.strip("\n") + "," + mutant + "," + CL_N2 + "," + CL_O2 + "," + LD + "," + str(fold_random) + "\n")
				
print(fold_random_dict)
