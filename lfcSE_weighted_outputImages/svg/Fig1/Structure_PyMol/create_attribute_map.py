import statistics

single_muts_dict = {}

with open("../heatMap/single_muts.csv") as f:
	header = True
	for line in f:
		if header:
			header = False
			continue
		split_line = line.strip("\n").split(",")
		AA_pos = split_line[0][:-1]
		if AA_pos not in single_muts_dict.keys():
			single_muts_dict[AA_pos] = [float(split_line[1])]
		else:
			single_muts_dict[AA_pos].append(float(split_line[1]))
			
mean_single_muts = {}
median_single_muts = {}
max_single_muts = {}

for aa_pos in single_muts_dict.keys():
	aa_number = int(str(aa_pos)[1:])-14
	mean_single_muts[aa_number] = statistics.mean(single_muts_dict[aa_pos])
	median_single_muts[aa_number] = statistics.median(single_muts_dict[aa_pos])
	max_single_muts[aa_number] = max(single_muts_dict[aa_pos])
	
with open("mean_single_muts.defattr", "w") as f:
	f.write("attribute: mean_norm_fitness\n")
	for aa_pos in mean_single_muts.keys():
		f.write("\t/A:" + str(aa_pos) + "\t" + str(mean_single_muts[aa_pos]) + "\n")
		
with open("median_single_muts.defattr", "w") as f:
	f.write("attribute: median_norm_fitness\n")
	for aa_pos in median_single_muts.keys():
		f.write("\t/A:" + str(aa_pos) + "\t" + str(median_single_muts[aa_pos]) + "\n")
		
with open("max_single_muts.defattr", "w") as f:
	f.write("attribute: max_norm_fitness\n")
	for aa_pos in max_single_muts.keys():
		f.write("\t/A:" + str(aa_pos) + "\t" + str(max_single_muts[aa_pos]) + "\n")
