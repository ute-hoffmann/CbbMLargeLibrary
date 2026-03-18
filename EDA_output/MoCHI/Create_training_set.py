# load original sequence
cbbM_seq = ""
with open("Gallionella_Rubisco_I.txt") as f:
	for line in f:
		cbbM_seq = line

header_for_mochi = "aa_seq\tNham_aa\tWT\tfitness\tsigma\n"

with open("CL_N2_variant_stdev.tsv", "w") as f2:
	f2.write(header_for_mochi)
	
with open("predict_CL_N2_variant_stdev.tsv", "w") as f2:
	f2.write(header_for_mochi)
	
with open("all_CL_N2_variant_stdev.tsv", "w") as f2:
	f2.write(header_for_mochi)
	
with open("CL_O2_variant_stdev.tsv", "w") as f2:
	f2.write(header_for_mochi)

with open("predict_CL_O2_variant_stdev.tsv", "w") as f2:
	f2.write(header_for_mochi)
	
with open("all_CL_O2_variant_stdev.tsv", "w") as f2:
	f2.write(header_for_mochi)
	
with open("LD_variant_stdev.tsv", "w") as f2:
	f2.write(header_for_mochi)

with open("predict_LD_variant_stdev.tsv", "w") as f2:
	f2.write(header_for_mochi)

with open("all_LD_variant_stdev.tsv", "w") as f2:
	f2.write(header_for_mochi)
	
with open("predict_only.tsv", "w") as f2:
	f2.write(header_for_mochi)

printed=False
with open("variant_stdev.csv") as f:
		header = True
		for line in f:
			if header:
				header = False
				continue
			if "STOP" in line:
				continue
			split_line = line.strip("\n").split("\t")
			condition = split_line[0]
			aa_combination = split_line[1]
			norm = split_line[2]
			sd_fitness = split_line[3]
			
			original_aas = []
			positions = []
			new_aas = []
			
			mutant = ""
			if not aa_combination:
				continue
			WT = False
			if aa_combination == "WT":
				WT = True
				
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

			Nham_aa = len(aa_combination.strip('"').split(","))
			if not WT:
				file_n = "predict_" + condition + "_variant_stdev.tsv"
				
				if condition == "CL_N2":
					with open("predict_only.tsv", "a") as f2:
						f2.write(new_seq.strip("\n") + "*" + "\t" + str(Nham_aa) + "\t\t0\t0\n")
			elif WT and condition == "CL_N2":
				with open("predict_only.tsv", "a") as f2:
					f2.write(cbbM_seq.strip("\n") + "*\t0\tTrue\t" + norm + "\t" + sd_fitness + "\n")
			if not combinatorial:
				file_n = condition + "_variant_stdev.tsv"
				with open(file_n, "a") as f2:
					if WT:
						f2.write(cbbM_seq.strip("\n") + "*" + "\t0\tTrue\t" + norm + "\t" + sd_fitness + "\n")
					else:
						f2.write(new_seq.strip("\n") + "*" + "\t" + str(Nham_aa) + "\t\t" + norm + "\t" + sd_fitness + "\n")
			file_n_2 = "all_" + condition + "_variant_stdev.tsv"
			with open(file_n_2, "a") as f2:
				if WT:
					f2.write(cbbM_seq.strip("\n") + "*" + "\t0\tTrue\t" + norm + "\t" + sd_fitness + "\n")
				else:
					f2.write(new_seq.strip("\n") + "*" + "\t" + str(Nham_aa) + "\t\t" + norm + "\t" + sd_fitness + "\n")
