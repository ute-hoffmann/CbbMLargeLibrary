# got saturational_for_proteinNPT.csv
## assign fold (name: 

# load original sequence
cbbM_seq = ""
with open("Gallionella_Rubisco_I.txt") as f:
	for line in f:
		cbbM_seq = line

fold_random = 4
baseAA = "A0"
with open("saturational_for_proteinNPT.csv") as f:
	with open("saturational_proteinNPT_with_foldChange.csv", "w") as f2:
		f2.write("mutated_sequence,mutant,CL_N2,CL_O2,LD,fold_random_5\n")
		header = True
		for line in f:
			if header:
				header = False
				continue
			split_line = line.strip("\n").split(",")
			aa_exchange = split_line[0]
			if baseAA != split_line[1]:
				fold_random += 1
				if fold_random == 5:
					fold_random = 0 
			baseAA = split_line[1]
			CL_N2 = split_line[2]
			CL_O2 = split_line[3]
			LD = split_line[4]
			original_aa = aa_exchange[0]
			new_aa = aa_exchange[-1]
			position = int(aa_exchange[1:-1])
			new_seq = ""
			for pos in range(1,len(cbbM_seq)+1):
				if pos == position:
					if original_aa != cbbM_seq[pos-1]:
						print(original_aa)
						print(cbbM_seq[pos-1])
						print(position)
						print(pos)
						break
					new_seq += new_aa
				else:
					new_seq += cbbM_seq[pos-1]
			if len(new_seq) != len(cbbM_seq):
				print("something is weird for " + aa_exchange + " and new seq " + new_seq)
			else:
				f2.write(new_seq.strip("\n") + "," + aa_exchange + "," + CL_N2 + "," + CL_O2 + "," + LD + "," + str(fold_random) + "\n")
