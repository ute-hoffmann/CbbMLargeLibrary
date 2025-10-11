# load original sequence
cbbM_seq = ""
with open("Gallionella_Rubisco_I.txt") as f:
	for line in f:
		cbbM_seq = line

list_of_seqs = []

with open("list_combinations.txt") as table_variants:
	number_mutants = 0
	for line in table_variants:
		mutant = line.strip("\n")
		if mutant == "":
			continue
		number_mutants += 1
		
		original_aas = []
		positions = []
		new_aas = []
		
		for aa_exchange in mutant.strip('"').split(":"):
			original_aas.append(aa_exchange[0])
			new_aas.append(aa_exchange[-1])
			positions.append(int(aa_exchange[1:-1]))
		
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
			print("something is weird for " + mutant + " and new seq " + new_seq)
		else:
			list_of_seqs.append(new_seq)

print("Number variants " + str(number_mutants))
list_of_seqs = set(list_of_seqs)
print("Number unique variants " + str(len(list_of_seqs)))


with open("list_new_sequences.txt", "w") as new_f:
	for seq in list_of_seqs:
		new_f.write(seq)


