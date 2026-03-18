# load original sequence
cbbM_seq = ""
with open("Gallionella_Rubisco_I.txt") as f:
	for line in f:
		cbbM_seq = line

printed=False
with open("variant_stdev.csv") as f:
	with open("mutant_code.tsv", "w") as f_mut:
		header = True
		for line in f:
			if header:
				header = False
				f_mut.write("aa_seq\tname\n")
				continue
			if "STOP" in line:
				continue
			split_line = line.strip("\n").split("\t")
			condition = split_line[0]
			aa_combination = split_line[1]
			
			original_aas = []
			positions = []
			new_aas = []
			
			mutant = ""
			if not aa_combination:
				continue
			WT = False
			if aa_combination == "WT":
				WT = True
				if condition=="CL_N2":
					f_mut.write(cbbM_seq.strip("\n") + "*" + "\tWT\n")
				
			combinatorial = False
			if len(aa_combination.strip('"').split(",")) > 1:
				combinatorial = True
			
			if not WT:
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
					continue


			if condition == "CL_N2" and not WT:
				f_mut.write(new_seq.strip("\n") + "*" + "\t" + aa_combination + "\n")
