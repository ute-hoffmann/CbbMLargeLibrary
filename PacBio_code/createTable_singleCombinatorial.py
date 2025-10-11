expected_variants = {}
combi_variants = {}

with open("interm_file/saturational_mutations.tsv") as f:
	for line in f:
		expected_variants[line.strip("\n")] = 0
		
with open("interm_file/combinatorial_mutations.tsv") as f:
	for line in f:
		expected_variants[line.strip("\n")] = 0
		combi_variants[line.strip("\n")] = 0

all_mutations = {}
all_barcodes = {}
total = 0
unique = 0
with open("output/selected_barcodes_Nextflow.fa") as mut_file:
	for line in mut_file:
		if line[0] == ">":
			total += 1
			barcode = line
			mutant = line.strip("\n").strip(">").split("|")[0]
			if barcode not in all_barcodes.keys():
				all_barcodes[barcode] = 1
				unique += 1
			if mutant not in all_mutations.keys():
				all_mutations[mutant] = 1
			else:
				all_mutations[mutant] += 1
print(total)
print(unique)
with open("../input/combi_satur_mutatios.tsv", "w") as f:
	for mut in all_mutations.keys():
		if mut == "":
			continue
		mutant = str(mut)
		if mutant == "WT":
			f.write(mutant + "\tWT\n")
		elif mutant not in expected_variants.keys():
			f.write(mutant + "\tnotExpected\n")
		elif len(mutant.split(",")) > 1:
			f.write(mutant + "\tcombinatorial\n")
		elif mutant in combi_variants.keys():
			f.write(mutant + "\tcombiANDsatur\n")
		else: 
			f.write(mutant + "\tsaturational\n")
			
with open("../input/number_mutations_variants.tsv", "w") as f:
	f.write("sgRNA_target\tnumber_muts\tnum_barcodes\n")
	for mut in all_mutations.keys():
		if mut == "":
			continue
		mutant = str(mut)
		if mutant == "WT":
			f.write(mutant + "\t0\t" + str(all_mutations[mut]) + "\n")
		else:
			f.write(mutant + "\t" + str(len(mutant.split(","))) + "\t" + str(all_mutations[mut]) + "\n")
